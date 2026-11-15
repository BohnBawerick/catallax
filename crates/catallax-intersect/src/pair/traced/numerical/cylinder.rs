//! One canonical finite cylinder source, independent of a cone query.

use super::*;

/// A retained numerical cylinder source about the original analytic origin.
/// Source controls and checks remain numerical, with no exact incidence authority.
#[derive(Clone, Debug)]
pub struct NumericalNativeCylinderSource {
    original: Cylinder,
    source: BSplineSurface,
    axis_parameters_mm: [Rational; 2],
    evidence: NumericalEvidence,
}
impl NumericalNativeCylinderSource {
    /// Original carrier, including its axis and seam direction.
    pub fn original(&self) -> &Cylinder {
        &self.original
    }
    /// Complete retained definition over native u=[0,4], v=[0,1].
    pub fn source(&self) -> &BSplineSurface {
        &self.source
    }
    /// Signed millimetres at native v=0 and v=1 from the original origin.
    pub fn axis_parameters_mm(&self) -> &[Rational; 2] {
        &self.axis_parameters_mm
    }
    /// Preset tolerance and actual source-control/analytic-carrier checks.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
}

/// Constructs the same native cylinder definition used by finite cone queries.
/// The interval is signed millimetres along the unit axis from the original
/// cylinder origin. No cone-dependent frame or parameter remapping is introduced.
pub fn numerical_native_cylinder_source(
    cylinder: &Cylinder,
    axis_parameters_mm: [Rational; 2],
    tolerance_mm: Rational,
    dec: &mut Decider,
) -> Result<NumericalNativeCylinderSource, Refusal> {
    let run = |dec: &mut Decider| build(cylinder, axis_parameters_mm, tolerance_mm, dec);
    if operation_budget_active() {
        return run(dec);
    }
    let (result, work) = with_operation_budget(20_000_000, || run(dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "numerical cylinder source exhausted its operation allowance",
        ));
    }
    result
}
fn build(
    cylinder: &Cylinder,
    parameters: [Rational; 2],
    tolerance: Rational,
    dec: &mut Decider,
) -> Result<NumericalNativeCylinderSource, Refusal> {
    dec.charge_work(1)?;
    if tolerance.sign() != Sign::Positive || parameters[0] >= parameters[1] {
        return Err(bad(
            "canonical cylinder source requires increasing finite coordinates and a positive tolerance",
        ));
    }
    if dec.policy().max_precision_bits < BITS {
        return Err(bad("native numerical sources require 128-bit enclosures"));
    }
    let radius = cylinder.radius();
    if enclosure(&radius, dec)?.lo().sign() != Sign::Positive {
        return Err(bad("canonical cylinder source requires a positive radius"));
    }
    let frame = Frame::new(
        cylinder.origin.clone(),
        &cylinder.axis,
        &cylinder.ref_dir,
        dec,
    )?;
    let retained = source(
        &frame,
        &parameters,
        [&radius, &radius],
        "cylinder",
        &tolerance,
        dec,
    )?;
    let mut checks = vec![retained.check];
    carrier_checks(
        &retained.surface,
        &frame,
        Some(&radius),
        None,
        "cylinder",
        &mut checks,
        dec,
    )?;
    checks.push(NumericalCheck::Passed {
        name: "native.cylinder-original-origin-parameter-frame".into(),
        scope: "canonical-cylinder".into(),
    });
    let evidence = NumericalEvidence::new(tolerance, checks, vec![])?;
    dec.check_deadline()?;
    Ok(NumericalNativeCylinderSource {
        original: cylinder.clone(),
        source: retained.surface,
        axis_parameters_mm: parameters,
        evidence,
    })
}
