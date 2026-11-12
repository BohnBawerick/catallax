//! Binary-noise cleanup: a recorded repair for files written from binary
//! floating point, such as a seam vertex at `y = -7.3e-16` or a rotation entry
//! of `2.4e-16`, where the writer meant zero.
//!
//! While a cleaned reading is active, every decimal written with at least ten
//! significant digits is replaced by the simplest rational within a relative
//! `1e-12` of it, when that rational has a denominator of at most 10,000.
//! Shorter decimals are taken as written. The reading reports the largest
//! change and the largest magnitude it read, from which the repair proves a
//! whole-file bound on how far any face moved. The exact reading always comes
//! first; a cleaned one is tried only after it refuses for incidence, an inexact
//! rotation or a self-intersection, and only its own exact checks admit it.

use std::cell::RefCell;

use catallax_numbers::{Integer, Rational};

use crate::{Curve, ImportRefusal, SourceFormat, StepFile, Surface, Vec3};

/// What a cleaned reading changed, in file units.
#[derive(Clone, Debug, Default)]
pub(crate) struct Cleaned {
    /// Largest change to any number.
    pub max_change: Rational,
    /// Largest magnitude of any number read.
    pub max_magnitude: Rational,
    /// Numbers changed.
    pub changed: usize,
}

thread_local! {
    static ACTIVE: RefCell<Option<Cleaned>> = const { RefCell::new(None) };
}

/// Significant digits of a decimal's mantissa.
fn significant(text: &str) -> usize {
    let mantissa = text.split(['E', 'e', 'D', 'd']).next().unwrap_or("");
    mantissa
        .chars()
        .filter(char::is_ascii_digit)
        .skip_while(|c| *c == '0')
        .count()
}

/// The simplest rational in `[lo, hi]`, for `0 < lo <= hi`, by continued
/// fractions: the one with the smallest denominator.
fn simplest_positive(lo: &Rational, hi: &Rational) -> Rational {
    let floor = Rational::from_integer(lo.floor());
    if floor == *lo {
        return lo.clone();
    }
    if &(&floor + &Rational::one()) <= hi {
        return &floor + &Rational::one();
    }
    // Both ends share the integer part; recurse on the reciprocal fractions.
    let inner = simplest_positive(
        &(hi - &floor).recip().expect("nonzero fraction"),
        &(lo - &floor).recip().expect("nonzero fraction"),
    );
    &floor + &inner.recip().expect("nonzero continued fraction")
}

/// The simplest rational within `tolerance` of `x`.
fn simplest(x: &Rational, tolerance: &Rational) -> Rational {
    let (lo, hi) = (x - tolerance, x + tolerance);
    if lo.sign() != catallax_numbers::Sign::Positive
        && hi.sign() != catallax_numbers::Sign::Negative
    {
        return Rational::zero();
    }
    if lo.sign() == catallax_numbers::Sign::Positive {
        simplest_positive(&lo, &hi)
    } else {
        -simplest_positive(&-hi, &-lo)
    }
}

/// Applies the active cleanup, if any, to one decimal read from `text`.
pub(crate) fn adjust(text: &str, value: Rational) -> Rational {
    ACTIVE.with(|active| {
        let mut active = active.borrow_mut();
        let Some(stats) = active.as_mut() else {
            return value;
        };
        stats.max_magnitude = stats.max_magnitude.clone().max(value.abs());
        if significant(text) < 10 {
            return value;
        }
        let tolerance = &value.abs().max(Rational::one())
            * &Rational::from_ratio(1, 1_000_000_000_000).expect("nonzero");
        let candidate = simplest(&value, &tolerance);
        if candidate == value || candidate.denominator() > &Integer::from_i64(10_000) {
            return value;
        }
        stats.max_change = stats.max_change.clone().max((&candidate - &value).abs());
        stats.changed += 1;
        candidate
    })
}

/// Reads `bytes` with cleanup active. `None` when nothing changed.
pub(crate) fn read(
    format: SourceFormat,
    bytes: &[u8],
) -> Result<Option<(StepFile, Cleaned)>, ImportRefusal> {
    ACTIVE.with(|active| *active.borrow_mut() = Some(Cleaned::default()));
    let file = format.read(bytes);
    let cleaned = ACTIVE
        .with(|active| active.borrow_mut().take())
        .expect("cleanup state set above");
    let file = file?;
    Ok((cleaned.changed > 0).then_some((file, cleaned)))
}

/// Whether a refusal of the exact reading may be binary noise.
pub(crate) fn retries(error: &ImportRefusal) -> bool {
    match error {
        ImportRefusal::NotSolid { detail } => detail.contains("exactly"),
        ImportRefusal::Unsupported { type_name, .. } => {
            type_name.contains("not exactly orthonormal")
        }
        ImportRefusal::SelfIntersecting { .. }
        | ImportRefusal::SelfIntersectionNotVerified { .. } => true,
        _ => false,
    }
}

fn norm_squared(v: &Vec3) -> Rational {
    v.iter().fold(Rational::zero(), |sum, c| sum + c * c)
}

/// A proven bound, in millimetres, on how far any point of any face moved,
/// and a conservative squared nearest-feature size of the cleaned file.
/// `None` when the file carries spline data, which this bound does not cover.
pub(crate) fn bound(file: &StepFile, cleaned: &Cleaned) -> Option<(Rational, Rational)> {
    let mut directions = Vec::new();
    let mut feature: Option<Rational> = None;
    let mut smaller = |q: Rational| {
        if q.sign() == catallax_numbers::Sign::Positive {
            feature = Some(feature.take().map_or(q.clone(), |f| f.min(q)));
        }
    };
    for face in file.solids.iter().flat_map(|s| &s.faces) {
        match &face.surface {
            Surface::Plane(p) => directions.extend([&p.axis, &p.ref_dir]),
            Surface::Cylinder(p, r)
            | Surface::Cone45(p, r)
            | Surface::Cone(p, r, _)
            | Surface::Sphere(p, r) => {
                directions.extend([&p.axis, &p.ref_dir]);
                smaller(r * r);
            }
            Surface::Torus(p, major, minor) => {
                directions.extend([&p.axis, &p.ref_dir]);
                smaller(major * major);
                smaller(minor * minor);
            }
            Surface::BSpline(_) => return None,
        }
        for bound in &face.loops {
            for e in &bound.edges {
                match &e.curve {
                    Curve::Line(_, d) => directions.push(d),
                    Curve::Circle(p, r) => {
                        directions.extend([&p.axis, &p.ref_dir]);
                        smaller(r * r);
                    }
                    Curve::Ellipse(p, a, b) => {
                        directions.extend([&p.axis, &p.ref_dir]);
                        smaller(a * a);
                        smaller(b * b);
                    }
                    Curve::BSpline(_) => return None,
                }
                let chord: Vec3 = std::array::from_fn(|i| &e.end.point[i] - &e.start.point[i]);
                smaller(norm_squared(&chord));
            }
        }
    }
    // Directions are ratios: the smallest squared norm, less what cleanup could
    // have taken from it, bounds how much a change turns them.
    let shortest = directions
        .iter()
        .map(|d| norm_squared(d))
        .min()
        .unwrap_or_else(Rational::one);
    let quarter = Rational::from_ratio(1, 4).expect("nonzero");
    let floor = shortest.min(Rational::one()) * quarter;
    let change = &cleaned.max_change;
    if &(change * change) * &Rational::from_i64(64) > floor {
        return None;
    }
    // A changed length moves a face by at most the change times the unit; a
    // changed ratio turns a direction of squared norm at least `floor` by at
    // most 4 change / sqrt(floor) <= 8 change / floor (floor <= 1/4), over
    // distances of at most 4 times the largest magnitude read.
    let scale = file.mm_per_unit.clone().max(Rational::one());
    let reach = &Rational::from_i64(4) * &(&cleaned.max_magnitude * &scale);
    let turned = &(&Rational::from_i64(8) * change) / &floor;
    let moved = &(change * &scale) + &(&turned * &reach);
    Some((moved, feature?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn noise_snaps_and_written_values_stay() {
        ACTIVE.with(|a| *a.borrow_mut() = Some(Cleaned::default()));
        assert_eq!(
            adjust("-7.347880795E-16", q("-7.347880795E-16")),
            Rational::zero()
        );
        assert_eq!(
            adjust("2.9999999999999996", q("2.9999999999999996")),
            q("3")
        );
        assert_eq!(
            adjust("0.33333333333333331", q("0.33333333333333331")),
            q("1/3")
        );
        // Ten digits that are not noise near a simple value stay.
        assert_eq!(
            adjust("0.7071067811865476", q("0.7071067811865476")),
            q("0.7071067811865476")
        );
        // Short decimals are what the writer meant.
        assert_eq!(adjust("1.0E-12", q("1.0E-12")), q("1.0E-12"));
        let stats = ACTIVE.with(|a| a.borrow_mut().take()).unwrap();
        assert_eq!(stats.changed, 3);
        assert!(stats.max_change < q("1e-15"));
        // Inactive: nothing changes.
        assert_eq!(
            adjust("2.9999999999999996", q("2.9999999999999996")),
            q("2.9999999999999996")
        );
    }

    #[test]
    fn simplest_rationals() {
        assert_eq!(simplest(&q("0.3333334"), &q("0.000001")), q("1/3"));
        assert_eq!(simplest(&q("-2.5"), &q("0.6")), q("-2"));
        assert_eq!(simplest(&q("1e-20"), &q("1e-12")), Rational::zero());
    }
}
