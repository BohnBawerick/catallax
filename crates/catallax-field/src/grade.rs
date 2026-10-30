//! The supported position functions of RECIPE.md section 4.16, shared by the
//! lattice radius and the gyroid level.
//!
//! A grade is a constant or `v0 + clamp(t, 0, 1) (v1 - v0)` with
//! `t = ((p - f) · (g - f)) / |g - f|²`. Equal endpoints are the constant they
//! equal, so a constant written either way compiles to the same definition. The
//! anchors are part coordinates: nothing here reads a body's bounds, so editing
//! the body never moves them.
//!
//! A terminal field certifies the subset whose anchors differ along exactly one
//! part axis, with family limits on the endpoints and on the slope
//! `K = |v1 - v0| / |g - f|`. Other valid grades stay valid definitions.

use crate::geometry;
use crate::metric::{PRECISION, Scalar};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_recipe::Grade;
use catallax_shapes::Vector3;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn f(r: &Rational) -> f64 {
    r.to_f64_nearest()
}
fn rational_dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    &(&a[0] * &b[0]) + &(&(&a[1] * &b[1]) + &(&a[2] * &b[2]))
}
fn imin(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().min(b.lo()), a.hi().min(b.hi())).unwrap()
}
fn imax(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().max(b.lo()), a.hi().max(b.hi())).unwrap()
}
fn clamp01(t: Interval) -> Interval {
    imin(
        imax(t, Interval::point(0.0).unwrap()),
        Interval::point(1.0).unwrap(),
    )
}

/// A clamped linear grade with distinct anchors and distinct endpoint values.
#[derive(Clone, Debug)]
pub(crate) struct Linear {
    pub(crate) from: [Rational; 3],
    pub(crate) direction: [Rational; 3],
    pub(crate) length_sq: Rational,
    pub(crate) start: Rational,
    pub(crate) end: Rational,
}

/// A compiled grade.
#[derive(Clone, Debug)]
pub(crate) enum GradeDef {
    Constant(Rational),
    Linear(Box<Linear>),
}

/// The piece of a clamped linear grade that holds over a whole box.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Piece {
    /// The constant start value, at or before the start anchor plane.
    Start,
    /// The linear ramp, between the anchor planes inclusive.
    Ramp,
    /// The constant end value, at or beyond the end anchor plane.
    End,
}

impl GradeDef {
    /// Compiles a constant or linear grade. A linear grade with coincident anchors
    /// refuses with `impossible`; equal endpoint values give the constant. A
    /// distance grade is not a position function here and returns `None`.
    pub(crate) fn compile(
        grade: &Grade,
        impossible: &'static str,
    ) -> Result<Option<Self>, Refusal> {
        Ok(Some(match grade {
            Grade::Constant(c) => Self::Constant(c.value().clone()),
            Grade::Linear {
                from,
                to,
                start,
                end,
            } => {
                if from == to {
                    return Err(Refusal::new(
                        impossible,
                        "a linear grade needs two distinct points",
                    ));
                }
                if start == end {
                    return Ok(Some(Self::Constant(start.value().clone())));
                }
                let from = from.0.each_ref().map(|n| n.value().clone());
                let to = to.0.each_ref().map(|n| n.value().clone());
                let direction: [Rational; 3] = std::array::from_fn(|i| &to[i] - &from[i]);
                Self::Linear(Box::new(Linear {
                    length_sq: rational_dot(&direction, &direction),
                    from,
                    direction,
                    start: start.value().clone(),
                    end: end.value().clone(),
                }))
            }
            Grade::Distance { .. } => return Ok(None),
        }))
    }

    /// The constant value, or `None` for a ramp.
    pub(crate) fn constant(&self) -> Option<&Rational> {
        match self {
            Self::Constant(c) => Some(c),
            Self::Linear(_) => None,
        }
    }

    /// The part axis the anchors differ along, when they differ along exactly one.
    pub(crate) fn axis(&self) -> Option<usize> {
        let Self::Linear(linear) = self else {
            return None;
        };
        let mut axes = (0..3).filter(|&i| !linear.direction[i].is_zero());
        let axis = axes.next()?;
        axes.next().is_none().then_some(axis)
    }

    /// A rational upper bound on the slope `|v1 - v0| / |g - f|`, exact when the
    /// anchors differ along one axis. Zero for a constant.
    pub(crate) fn slope_bound(&self) -> Rational {
        let Self::Linear(linear) = self else {
            return q(0);
        };
        let rise = (&linear.end - &linear.start).abs();
        let length = match self.axis() {
            Some(axis) => linear.direction[axis].abs(),
            None => Ball::from_rational(&linear.length_sq, PRECISION)
                .sqrt()
                .expect("positive squared length")
                .lo(),
        };
        rise.checked_div(&length)
            .expect("a positive lower bound on a nonzero length")
    }

    fn parameter(linear: &Linear, p: &[Rational; 3]) -> Rational {
        let rel: [Rational; 3] = std::array::from_fn(|i| &p[i] - &linear.from[i]);
        rational_dot(&rel, &linear.direction)
            .checked_div(&linear.length_sq)
            .expect("distinct grade points")
    }

    /// The exact value at a rational point.
    pub(crate) fn value(&self, p: &[Rational; 3]) -> Rational {
        match self {
            Self::Constant(c) => c.clone(),
            Self::Linear(linear) => {
                let t = Self::parameter(linear, p).max(q(0)).min(q(1));
                &linear.start + &(&t * &(&linear.end - &linear.start))
            }
        }
    }

    /// The tentative value in binary64.
    pub(crate) fn fast(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Constant(c) => f(c),
            Self::Linear(linear) => {
                let Linear {
                    from,
                    direction,
                    length_sq,
                    start,
                    end,
                } = linear.as_ref();
                let rel: [f64; 3] = std::array::from_fn(|i| p[i] - f(&from[i]));
                let direction = direction.each_ref().map(f);
                let t = rel[0]
                    .mul(&direction[0])
                    .add(&rel[1].mul(&direction[1]))
                    .add(&rel[2].mul(&direction[2]))
                    / f(length_sq);
                let t = t.clamp(0.0, 1.0);
                f(start) + t * (f(end) - f(start))
            }
        }
    }

    fn interval_parameter(linear: &Linear, p: &[Interval; 3]) -> Interval {
        let rel: [Interval; 3] =
            std::array::from_fn(|i| p[i] - Interval::from_rational(&linear.from[i]));
        (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
            s + rel[i] * Interval::from_rational(&linear.direction[i])
        }) / Interval::from_rational(&linear.length_sq)
    }

    /// An enclosure of every value over a box.
    pub(crate) fn interval(&self, p: &[Interval; 3]) -> Interval {
        match self {
            Self::Constant(c) => Interval::from_rational(c),
            Self::Linear(linear) => {
                let t = clamp01(Self::interval_parameter(linear, p));
                Interval::from_rational(&linear.start)
                    + t * (Interval::from_rational(&linear.end)
                        - Interval::from_rational(&linear.start))
            }
        }
    }

    /// The piece that holds over the whole box, or `None` when the box meets an
    /// anchor plane's two sides. A constant is one ramp-free piece everywhere.
    pub(crate) fn piece(&self, p: &[Interval; 3]) -> Option<Piece> {
        let Self::Linear(linear) = self else {
            return Some(Piece::Start);
        };
        let t = Self::interval_parameter(linear, p);
        if t.hi() <= 0.0 {
            Some(Piece::Start)
        } else if t.lo() >= 1.0 {
            Some(Piece::End)
        } else if t.lo() >= 0.0 && t.hi() <= 1.0 {
            Some(Piece::Ramp)
        } else {
            None
        }
    }

    /// The unclamped affine extension of the ramp over a box: `v0 + t (v1 - v0)`.
    /// It equals the grade between the anchor planes. A constant returns itself.
    pub(crate) fn affine(&self, p: &[Interval; 3]) -> Interval {
        match self {
            Self::Constant(c) => Interval::from_rational(c),
            Self::Linear(linear) => {
                Interval::from_rational(&linear.start)
                    + Self::interval_parameter(linear, p)
                        * (Interval::from_rational(&linear.end)
                            - Interval::from_rational(&linear.start))
            }
        }
    }

    /// The gradient of the affine extension, enclosed. Zero for a constant.
    pub(crate) fn gradient(&self) -> [Interval; 3] {
        match self {
            Self::Constant(_) => [Interval::point(0.0).unwrap(); 3],
            Self::Linear(linear) => {
                let scale = (&linear.end - &linear.start)
                    .checked_div(&linear.length_sq)
                    .expect("distinct grade points");
                std::array::from_fn(|i| Interval::from_rational(&(&linear.direction[i] * &scale)))
            }
        }
    }

    /// The gradient of one piece, enclosed: the ramp's, or zero on a clamp.
    pub(crate) fn piece_gradient(&self, piece: Piece) -> [Interval; 3] {
        match piece {
            Piece::Ramp => self.gradient(),
            Piece::Start | Piece::End => [Interval::point(0.0).unwrap(); 3],
        }
    }

    /// The exact range of `t` over a closed rational box.
    fn parameter_range(linear: &Linear, lo: &[Rational; 3], hi: &[Rational; 3]) -> [Rational; 2] {
        // `t` is affine, so its extremes are at the corners selected by the signs of
        // the direction's components.
        let corner = |upper: bool| -> [Rational; 3] {
            std::array::from_fn(|i| {
                if (linear.direction[i].sign() == Sign::Positive) == upper {
                    hi[i].clone()
                } else {
                    lo[i].clone()
                }
            })
        };
        [corner(false), corner(true)].map(|p| Self::parameter(linear, &p))
    }

    /// The piece that holds over a whole closed rational box, decided exactly, or
    /// `None` when the box meets both sides of an anchor plane.
    pub(crate) fn exact_piece(&self, lo: &[Rational; 3], hi: &[Rational; 3]) -> Option<Piece> {
        let Self::Linear(linear) = self else {
            return Some(Piece::Start);
        };
        let [first, last] = Self::parameter_range(linear, lo, hi);
        if last <= q(0) {
            Some(Piece::Start)
        } else if first >= q(1) {
            Some(Piece::End)
        } else if first >= q(0) && last <= q(1) {
            Some(Piece::Ramp)
        } else {
            None
        }
    }

    /// Whether the grade is one affine function over a closed rational box: between
    /// its anchor planes inclusive, or wholly on one side of them.
    pub(crate) fn affine_over(&self, lo: &[Rational; 3], hi: &[Rational; 3]) -> bool {
        self.exact_piece(lo, hi).is_some()
    }

    /// The axis and positions of the two anchor planes of a grade whose anchors
    /// differ along one part axis.
    pub(crate) fn planes(&self) -> Option<(usize, [Rational; 2])> {
        let axis = self.axis()?;
        let Self::Linear(linear) = self else {
            return None;
        };
        let from = linear.from[axis].clone();
        let to = &from + &linear.direction[axis];
        Some((axis, [from, to]))
    }

    /// The exact range of the grade over a closed rational box.
    pub(crate) fn range(&self, lo: &[Rational; 3], hi: &[Rational; 3]) -> [Rational; 2] {
        match self {
            Self::Constant(c) => [c.clone(), c.clone()],
            Self::Linear(linear) => {
                let [a, b] = Self::parameter_range(linear, lo, hi).map(|t| {
                    let t = t.max(q(0)).min(q(1));
                    &linear.start + &(&t * &(&linear.end - &linear.start))
                });
                [a.clone().min(b.clone()), a.max(b)]
            }
        }
    }

    /// The rate of change along `d` in the positive-infinitesimal direction from `p`,
    /// or `None` where the grade is constant that way.
    pub(crate) fn rate(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let Self::Linear(linear) = self else {
            return Ok(None);
        };
        let Linear {
            direction,
            length_sq,
            start,
            end,
            ..
        } = linear.as_ref();
        let t = Self::parameter(linear, p);
        let along = &(&Real::from(direction[0].clone()) * &d.x)
            + &(&(&Real::from(direction[1].clone()) * &d.y)
                + &(&Real::from(direction[2].clone()) * &d.z));
        let along_sign = geometry::sign(&along, dec)?;
        let ramp = if t > q(0) && t < q(1) {
            true
        } else if t == q(0) {
            along_sign == Sign::Positive
        } else if t == q(1) {
            along_sign == Sign::Negative
        } else {
            false
        };
        Ok(ramp.then(|| {
            &Real::from(
                (end - start)
                    .checked_div(length_sq)
                    .expect("distinct grade points"),
            ) * &along
        }))
    }

    /// Composes a rational rigid motion `x -> M x + t` into the anchors.
    pub(crate) fn transform(&mut self, m: &[[Rational; 3]; 3], t: &[Rational; 3]) {
        if let Self::Linear(linear) = self {
            let apply = |v: &[Rational; 3]| -> [Rational; 3] {
                std::array::from_fn(|i| rational_dot(&m[i], v))
            };
            let moved = apply(&linear.from);
            linear.from = std::array::from_fn(|i| &moved[i] + &t[i]);
            linear.direction = apply(&linear.direction);
        }
    }
}

/// Checks the terminal certification subset of one field family. A constant, or
/// a linear grade whose endpoints are equal, is always in it. Otherwise the
/// anchors must differ along exactly one part axis (`{family}.grade-function-unsupported`),
/// both endpoints must lie in `[low, high]` and the slope must be at most
/// `max_slope` (`{family}.grade-range-unsupported`). Nothing is clamped.
pub(crate) fn check_certified(
    grade: &GradeDef,
    [function, range]: [&'static str; 2],
    [low, high]: [&Rational; 2],
    max_slope: &Rational,
) -> Result<(), Refusal> {
    let GradeDef::Linear(linear) = grade else {
        return Ok(());
    };
    if grade.axis().is_none() {
        return Err(Refusal::new(
            function,
            "a terminal field certifies linear grades whose anchors differ along one part axis",
        ));
    }
    for value in [&linear.start, &linear.end] {
        if value < low || value > high {
            return Err(Refusal::new(
                range,
                format!("grade endpoint {value} is outside the certified range [{low}, {high}]"),
            ));
        }
    }
    let slope = grade.slope_bound();
    if &slope > max_slope {
        return Err(Refusal::new(
            range,
            format!("grade slope {slope} per mm exceeds the certified {max_slope} per mm"),
        ));
    }
    Ok(())
}

/// The exact range of a constant or linear grade over a closed rational box, as a
/// terminal field reads it. A box across an anchor plane holds the clamped value,
/// never an unclamped endpoint or the value at the box centre.
///
/// Refusals: `field.grade-box` when a lower corner exceeds its upper corner,
/// `field.grade-unsupported` for a distance grade, the family-neutral
/// `field.grade-impossible` for coincident anchors, and the decider's work and
/// deadline outcomes.
pub fn grade_range(
    grade: &Grade,
    lo: &[Rational; 3],
    hi: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    dec.check_deadline()?;
    if (0..3).any(|i| lo[i] > hi[i]) {
        return Err(Refusal::new(
            "field.grade-box",
            "each lower box corner must not exceed its upper corner",
        ));
    }
    let Some(compiled) = GradeDef::compile(grade, "field.grade-impossible")? else {
        return Err(Refusal::new(
            "field.grade-unsupported",
            "distance grades need trimmed-face distance; unfinished coverage",
        ));
    };
    Ok(compiled.range(lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grade(text: &str) -> Grade {
        serde_json::from_str(text).unwrap()
    }
    fn r(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn exact_range_clamps_across_an_anchor_plane() {
        let level = grade(
            r#"{"linear":{"from":["0","0","-1"],"to":["0","0","1"],"start":"-1/10","end":"1/10"}}"#,
        );
        let range = grade_range(
            &level,
            &[r("0"), r("0"), r("999/1000")],
            &[r("0"), r("0"), r("1001/1000")],
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(range, [r("999/10000"), r("1/10")]);
    }

    #[test]
    fn equal_endpoints_compile_to_the_constant() {
        let equal = grade(
            r#"{"linear":{"from":["0","0","20"],"to":["0","0","30"],"start":"1","end":"1"}}"#,
        );
        let compiled = GradeDef::compile(&equal, "lattice.impossible")
            .unwrap()
            .unwrap();
        assert_eq!(compiled.constant(), Some(&r("1")));
    }
}
