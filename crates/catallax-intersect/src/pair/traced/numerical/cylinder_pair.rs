//! Coupled numerical curves on two retained finite native cylinders.

use super::*;

/// All reported fitted pieces, with both original cylinder parameter frames.
/// Samples carry no exact incidence or complete branch-count authority.
#[derive(Clone, Debug)]
pub struct NumericalCylinderIntersection {
    sources: [NumericalNativeCylinderSource; 2],
    fits: checked::NativeFits,
}

impl NumericalCylinderIntersection {
    /// Original cylinders, native sources and signed axis intervals, in operand order.
    pub fn sources(&self) -> &[NumericalNativeCylinderSource; 2] {
        &self.sources
    }
    /// Retained 3D curves and both coupled UV definitions.
    pub fn pieces(&self) -> &[NumericalIntersectionPiece] {
        &self.fits.pieces
    }
    /// Source, finite-domain and sampled coupled-fit checks at 1e-7 mm.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.fits.evidence
    }
    /// Actual request work by stage, including source construction.
    pub fn work_stages(&self) -> NumericalNativeWork {
        self.fits.work_stages
    }
}

/// Searches all sixteen native quarter pairs at 1e-7 mm, retaining every
/// reported 3D/UV fit. Axis intervals are signed millimetres from each original
/// cylinder origin. This produces no Boolean topology or exact geometry facts.
pub fn numerical_native_cylinder_intersection(
    cylinders: [&Cylinder; 2],
    axis_parameters_mm: [[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<NumericalCylinderIntersection, Refusal> {
    let run = |dec: &mut Decider| build(cylinders, axis_parameters_mm, dec);
    if operation_budget_active() {
        return run(dec);
    }
    let (result, work) = with_operation_budget(20_000_000, || run(dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "finite cylinder intersection exhausted its operation allowance",
        ));
    }
    result
}

fn build(
    cylinders: [&Cylinder; 2],
    axis_parameters_mm: [[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<NumericalCylinderIntersection, Refusal> {
    let started = catallax_numbers::current_operation_work().used;
    let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
    let [first, second] = axis_parameters_mm;
    let sources = [
        numerical_native_cylinder_source(cylinders[0], first, tolerance.clone(), dec)?,
        numerical_native_cylinder_source(cylinders[1], second, tolerance.clone(), dec)?,
    ];
    let mut checks = Vec::new();
    for (side, source) in sources.iter().enumerate() {
        for check in source.evidence().checks() {
            dec.charge_work(1)?;
            let mut check = check.clone();
            match &mut check {
                NumericalCheck::Passed { scope, .. } | NumericalCheck::Samples { scope, .. } => {
                    *scope = format!("cylinder-{side}/{scope}");
                }
            }
            checks.push(check);
        }
    }
    let evidence = NumericalEvidence::new(tolerance, checks, vec![])?;
    let after_sources = catallax_numbers::current_operation_work().used;
    let fits = checked::search_sources(
        &[sources[0].source().clone(), sources[1].source().clone()],
        &cylinders.map(|c| Surface::Cylinder(c.clone())),
        &evidence,
        true,
        started,
        after_sources,
        dec,
    )?;
    Ok(NumericalCylinderIntersection { sources, fits })
}
