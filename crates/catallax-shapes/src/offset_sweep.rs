//! Checked whole-domain correspondence for finite elliptic-offset sweeps.

use crate::{BSplineSurface, EllipticOffsetDefinition, OffsetSide, Point3, Surface, Vector3};
use catallax_numbers::{
    Decider, EvidenceBound, EvidenceUnit, Outcome, Rational, Real, Refusal, StrictSign,
};

/// Verified regular convex sweep and its complete signed-definition correspondence.
/// Private fields bind this proof to immutable inputs. It is not a body proof:
/// cap domains, shared incidence and outward material orientation remain required.
#[derive(Debug)]
pub struct OffsetSweepProof<'a> {
    surface: &'a BSplineSurface,
    definition: &'a EllipticOffsetDefinition,
    axial: [&'a Real; 2],
    bound: EvidenceBound,
    spans: usize,
}

impl<'a> OffsetSweepProof<'a> {
    /// Immutable carrier whose full knot cover was checked.
    pub fn surface(&self) -> &'a BSplineSurface {
        self.surface
    }
    /// Exact signed normal-offset definition whose interpolation was checked.
    pub fn definition(&self) -> &'a EllipticOffsetDefinition {
        self.definition
    }
    /// Exact lower and upper axial coordinates of the checked sweep.
    pub fn axial(&self) -> [&'a Real; 2] {
        self.axial
    }
    /// Uniform two-way position-correspondence bound, in millimetres.
    pub fn bound(&self) -> &EvidenceBound {
        &self.bound
    }
    /// Complete number of nonzero cross-section spans.
    pub fn span_count(&self) -> usize {
        self.spans
    }
    /// Consumes the proof to retain its metric bound without copying it.
    pub fn into_bound(self) -> EvidenceBound {
        self.bound
    }
}
fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.offset-sweep", detail)
}
fn same(a: &Real, b: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    dec.check_deadline()?;
    if a.same_definition(b) {
        return Ok(());
    }
    match dec.compare(a, b) {
        Outcome::ProvenCoincident(_) => Ok(()),
        Outcome::Proven(_) => Err(refusal(
            "carrier does not interpolate its retained signed offset definition",
        )),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(error) => Err(error),
    }
}
fn same_point(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<(), Refusal> {
    for (a, b) in [&a.x, &a.y, &a.z].into_iter().zip([&b.x, &b.y, &b.z]) {
        same(a, b, dec)?;
    }
    Ok(())
}
fn positive(value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(()),
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => Err(refusal(
            "regular convex offset-sweep condition is not satisfied",
        )),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(error) => Err(error),
    }
}

fn convex_span(
    controls: [&Point3; 4],
    reference: &Vector3,
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let differences: [_; 3] = std::array::from_fn(|i| controls[i + 1].minus(controls[i]));
    // The entire derivative lies in one open half-plane. Combined with strict
    // curvature, this rules out an extra tangent winding between endpoint data.
    for difference in &differences {
        positive(&difference.dot(reference), dec)?;
    }
    let second = [
        differences[1].minus(&differences[0]),
        differences[2].minus(&differences[1]),
    ];
    let twice = Real::from_i64(2);
    // Bernstein coefficients of C' cross C'', with positive common factors
    // omitted (and the two interior coefficients multiplied by three).
    let curvature = [
        differences[0].cross(&second[0]),
        differences[0]
            .cross(&second[1])
            .plus(&differences[1].cross(&second[0]).scaled(&twice)),
        differences[1]
            .cross(&second[1])
            .scaled(&twice)
            .plus(&differences[2].cross(&second[0])),
        differences[2].cross(&second[1]),
    ];
    for coefficient in curvature {
        positive(&coefficient.dot(axis), dec)?;
    }
    Ok(())
}

/// Verifies a closed cubic Hermite cross-section swept exactly along the
/// retained elliptic-cylinder axis. It recomputes interpolation, regularity,
/// convexity and full-span fourth-derivative bounds; producer bounds are not input.
/// The common [0,4] x [0,1] parameter domain supplies both directions of coverage.
/// Strict positive curvature, the endpoint tangent data and the open derivative
/// half-planes give a single full tangent winding, hence an embedded convex loop.
///
/// This first checker covers polynomial degree-(3,1) composite Bezier sweeps,
/// unit weights and the regular normal branch. Other carriers are coverage gaps,
/// not declarations that their geometry is impossible. It never accepts a
/// cropped domain, detached control column, wrong normal side or forged seam.
pub fn certify_elliptic_offset_sweep<'a>(
    surface: &'a BSplineSurface,
    definition: &'a EllipticOffsetDefinition,
    axial: [&'a Real; 2],
    dec: &mut Decider,
) -> Result<OffsetSweepProof<'a>, Refusal> {
    dec.check_deadline()?;
    let domain = surface.domain();
    let zero = Rational::zero();
    let one = Rational::one();
    let four = Rational::from_i64(4);
    if surface.degrees() != [3, 1]
        || surface.control_count()[1] != 2
        || domain != [[&zero, &four], [&zero, &one]]
    {
        return Err(refusal(
            "whole polynomial composite-Bezier sweep domain required",
        ));
    }
    let real_one = Real::one();
    for weight in surface.weights() {
        same(weight, &real_one, dec)?;
    }
    positive(&(axial[1] - axial[0]), dec)?;
    let Surface::EllipticCylinder(original) = definition.original() else {
        unreachable!("validated elliptic-offset definition invariant");
    };
    if definition.side() == OffsetSide::Inward {
        let reach = &original.minor_sq / &original.major_sq.sqrt();
        positive(&(&reach - definition.distance()), dec)?;
    }
    let chart = definition.chart();
    let axis = chart.axial_direction();
    let sweep = axis.scaled(&(axial[1] - axial[0]));
    let mut upper = Rational::zero();
    let mut count = 0;
    let mut previous_end: Option<(Point3, Vector3)> = None;
    let mut quarter = 0;
    let mut previous_parameter = Rational::zero();
    for (index, span) in surface.spans(0).enumerate() {
        dec.check_deadline()?;
        let controls = surface
            .controls()
            .get(index * 6..index * 6 + 8)
            .ok_or_else(|| refusal("every nonzero span must have four cubic controls"))?;
        if span[0] != &previous_parameter {
            return Err(refusal(
                "missing or out-of-order span in full correspondence cover",
            ));
        }
        if span[0] == &Rational::from_i64((quarter + 1) as i64) {
            quarter += 1;
        }
        if quarter >= 4 || span[1] > &Rational::from_i64((quarter + 1) as i64) {
            return Err(refusal(
                "a correspondence span crosses an exact quarter-chart boundary",
            ));
        }
        let origin = Rational::from_i64(quarter as i64);
        let interval = [span[0] - &origin, span[1] - &origin];
        // The exact quarter-chart contract includes its shared endpoint and
        // tangent across quarter cuts and closure; reuse that endpoint once.
        let start = match previous_end.take() {
            Some(endpoint) => endpoint,
            None => chart.point_and_tangent(quarter, &interval[0], axial[0], dec)?,
        };
        let end = chart.point_and_tangent(quarter, &interval[1], axial[0], dec)?;
        let step = Real::from_rational((span[1] - span[0]) / Rational::from_i64(3));
        let expected = [
            start.0.clone(),
            start.0.plus(&start.1.scaled(&step)),
            end.0.plus(&end.1.scaled(&(-&step))),
            end.0.clone(),
        ];
        for (i, point) in expected.iter().enumerate() {
            same_point(&controls[2 * i], point, dec)?;
            same_point(&controls[2 * i + 1], &controls[2 * i].plus(&sweep), dec)?;
        }
        let reference = start.1.plus(&end.1);
        convex_span(
            [&controls[0], &controls[2], &controls[4], &controls[6]],
            &reference,
            axis,
            dec,
        )?;
        let bound = chart.hermite_error_bound(quarter, [&interval[0], &interval[1]], dec)?;
        upper = upper.max(bound.upper);
        previous_parameter = span[1].clone();
        previous_end = Some(end);
        count += 1;
    }
    if previous_parameter != four || surface.control_count()[0] != count * 3 + 1 {
        return Err(refusal(
            "incomplete composite-Bezier control or parameter cover",
        ));
    }
    Ok(OffsetSweepProof {
        surface,
        definition,
        axial,
        bound: EvidenceBound {
            upper,
            unit: EvidenceUnit::Millimetres,
        },
        spans: count,
    })
}
