//! Closed cubic normal-offset candidates, before topology and field acceptance.

use catallax_numbers::{Decider, EvidenceBound, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{BSplineSurface, EllipticOffsetDefinition, OffsetSide, Surface};

fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("features.offset-fit", detail)
}
fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(true),
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(error) => Err(error),
    }
}

/// One complete interpolation span's analytic correspondence bound.
#[derive(Clone, Debug)]
pub struct OffsetSpanBound {
    /// Exact normal-offset quarter chart, from zero through three.
    pub quarter: usize,
    /// Closed local chart interval, within zero through one.
    pub parameter: [Rational; 2],
    /// Uniform position-error upper bound, derived from fourth derivatives.
    pub error: EvidenceBound,
}

/// An unpublished fitted carrier with its retained exact construction target.
/// Field correspondence, convexity, common trims and body publication are
/// separate checks; this object cannot be used as a certified body.
#[derive(Clone, Debug)]
pub struct EllipticOffsetFit {
    /// Cubic cross-section swept exactly with degree one along the original axis.
    pub surface: BSplineSurface,
    /// Retained exact original support, offset distance and normal side.
    pub definition: EllipticOffsetDefinition,
    /// Exact axial coordinate interval relative to the original support frame.
    pub axial: [Real; 2],
    /// Complete ordered quarter-chart span cover and its interpolation bounds.
    pub spans: Vec<OffsetSpanBound>,
    /// Maximum of the complete parameter-correspondence bounds, not a body metric.
    pub position_error: EvidenceBound,
}

/// Fits a closed C1 cubic cross-section from exact positions and derivatives,
/// then sweeps it exactly along the original axis. Refinement is deterministic.
/// Inward distances at or beyond the normal neighbourhood remain a named gap,
/// not a claim that the exact erosion is impossible or empty.
pub fn fit_elliptic_offset(
    definition: &EllipticOffsetDefinition,
    axial: [Real; 2],
    target_mm: &Rational,
    dec: &mut Decider,
) -> Result<EllipticOffsetFit, Refusal> {
    dec.check_deadline()?;
    if target_mm <= &Rational::zero() || !positive(&(&axial[1] - &axial[0]), dec)? {
        return Err(refusal(
            "positive fit target and ordered axial interval required",
        ));
    }
    if definition.side() == OffsetSide::Inward {
        let Surface::EllipticCylinder(original) = definition.original() else {
            unreachable!("validated elliptic-offset definition invariant");
        };
        let reach = &original.minor_sq / &original.major_sq.sqrt();
        if !positive(&(&reach - definition.distance()), dec)? {
            return Err(refusal(
                "inward offset needs a different branch at the ellipse's normal-neighbourhood limit",
            ));
        }
    }
    let chart = definition.chart();
    let mut subdivisions = 1_usize;
    let spans = loop {
        let n = i64::try_from(subdivisions)
            .map_err(|_| refusal("offset subdivision count overflow"))?;
        let mut bounds = Vec::new();
        let mut fits = true;
        'cover: for quarter in 0..4 {
            for i in 0..n {
                dec.check_deadline()?;
                let lo = Rational::from_ratio(i, n).unwrap();
                let hi = Rational::from_ratio(i + 1, n).unwrap();
                let error = chart.hermite_error_bound(quarter, [&lo, &hi], dec)?;
                if &error.upper > target_mm {
                    fits = false;
                    break 'cover;
                }
                bounds.push(OffsetSpanBound {
                    quarter,
                    parameter: [lo, hi],
                    error,
                });
            }
        }
        if fits {
            break bounds;
        }
        subdivisions = subdivisions
            .checked_mul(2)
            .ok_or_else(|| refusal("offset subdivision count overflow"))?;
    };
    let count = spans.len();
    let denominator =
        i64::try_from(subdivisions).map_err(|_| refusal("offset subdivision count overflow"))?;
    let third_step =
        Real::from_rational(Rational::from_ratio(1, denominator).unwrap() / Rational::from_i64(3));
    let sweep = chart.axial_direction().scaled(&(&axial[1] - &axial[0]));
    let mut nodes = Vec::with_capacity(count);
    for i in 0..count {
        dec.check_deadline()?;
        let quarter = i / subdivisions;
        let parameter = Rational::from_ratio((i % subdivisions) as i64, denominator).unwrap();
        nodes.push(chart.point_and_tangent(quarter, &parameter, &axial[0], dec)?);
    }
    let mut controls = Vec::with_capacity((count * 3 + 1) * 2);
    let mut append = |point: catallax_shapes::Point3| {
        controls.push(point.clone());
        controls.push(point.plus(&sweep));
    };
    append(nodes[0].0.clone());
    for i in 0..count {
        dec.check_deadline()?;
        let (point, tangent) = &nodes[i];
        let (next, next_tangent) = &nodes[(i + 1) % count];
        append(point.plus(&tangent.scaled(&third_step)));
        append(next.plus(&next_tangent.scaled(&(-&third_step))));
        append(next.clone());
    }
    let mut knots = vec![Rational::zero(); 4];
    for i in 1..count {
        let parameter = Rational::from_ratio(
            i64::try_from(i).map_err(|_| refusal("offset knot count overflow"))?,
            denominator,
        )
        .unwrap();
        knots.extend(std::iter::repeat_n(parameter, 3));
    }
    knots.extend(std::iter::repeat_n(Rational::from_i64(4), 4));
    let control_count = controls.len();
    let surface = BSplineSurface::new(
        [3, 1],
        [
            knots,
            vec![
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::one(),
            ],
        ],
        [control_count / 2, 2],
        controls,
        vec![Real::one(); control_count],
        dec,
    )?;
    let upper = spans
        .iter()
        .map(|span| &span.error.upper)
        .max()
        .unwrap()
        .clone();
    Ok(EllipticOffsetFit {
        surface,
        definition: definition.clone(),
        axial,
        position_error: EvidenceBound {
            upper,
            unit: catallax_numbers::EvidenceUnit::Millimetres,
        },
        spans,
    })
}
