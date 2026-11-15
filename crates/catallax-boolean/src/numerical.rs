//! Numerical sources for open face consumers, separate from exact incidence.

use catallax_numbers::{Decider, NumericalCheck, NumericalEvidence, Rational, Refusal};
use catallax_shapes::{BSplineSurface, Cylinder, RationalPatch, SplineSupport, Surface};
use catallax_topology::numerical::NumericalOptions;

/// A canonical cylinder chart, independent of a cone face or assembly key.
#[derive(Clone, Debug)]
pub struct NumericalCylinderQuarter {
    /// Native quadrant, 0 through 3.
    pub id: usize,
    /// Deterministic descriptor-local key. The caller supplies assembly identities.
    pub key: String,
    /// Retained embedded patch over native u=[id,id+1], v=[0,1].
    pub surface: Surface,
    /// Signed millimetres at native v=0 and v=1 from the original tool origin.
    pub axis_parameters_mm: [Rational; 2],
}

/// A retained cylinder source and all four native quarter patches.
/// Constructing it supplies numerical source/control checks and patch embedding
/// checks. It does not construct a Boolean, admitted body or STEP result.
#[derive(Clone, Debug)]
pub struct NumericalCylinderSource {
    /// Original analytic carrier, with the caller's origin, axis and seam.
    pub original: Cylinder,
    /// Complete finite numerical definition, without native knot changes.
    pub source: BSplineSurface,
    /// Signed millimetres at native v=0 and v=1 from the original origin.
    pub axis_parameters_mm: [Rational; 2],
    /// Canonical quarters in increasing native-u order.
    pub quarters: Vec<NumericalCylinderQuarter>,
    /// Preset tolerance and actual source/control/chart checks.
    pub evidence: NumericalEvidence,
    /// Actual work for source construction and its checks.
    pub source_work: u64,
    /// Actual work for admitting every retained quarter patch.
    pub chart_work: u64,
}
impl NumericalCylinderSource {
    /// Creates the same cylinder source used by finite cone intersections.
    /// Coordinates are signed millimetres along the unit cylinder axis from its
    /// original origin. Reuse the same carrier and interval across all face calls.
    pub fn new(
        cylinder: &Cylinder,
        axis_parameters_mm: [Rational; 2],
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let run = |dec: &mut Decider| build(cylinder, axis_parameters_mm, options, dec);
        if catallax_numbers::operation_budget_active() {
            return run(dec);
        }
        let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || run(dec));
        if work.exhausted {
            return Err(Refusal::new(
                "numbers.operation-budget",
                "numerical cylinder quarter construction exhausted its operation allowance",
            ));
        }
        result
    }
    /// Maps a native v in [0,1] to signed millimetres about the original origin.
    pub fn signed_mm_at_v(&self, v: &Rational) -> Result<Rational, Refusal> {
        if v < &Rational::zero() || v > &Rational::one() {
            return Err(Refusal::new(
                "boolean.numerical-cylinder-source",
                "native cylinder v leaves [0,1]",
            ));
        }
        let parameters = &self.axis_parameters_mm;
        Ok(&parameters[0] + &((&parameters[1] - &parameters[0]) * v.clone()))
    }
}
fn build(
    cylinder: &Cylinder,
    parameters: [Rational; 2],
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalCylinderSource, Refusal> {
    let before = catallax_numbers::current_operation_work().used;
    let retained = catallax_intersect::numerical_native_cylinder_source(
        cylinder,
        parameters.clone(),
        options.tolerance_mm,
        dec,
    )?;
    let after_source = catallax_numbers::current_operation_work().used;
    let source = retained.source();
    let knots = source.knots();
    let mut quarters = Vec::new();
    let mut checks = retained.evidence().checks().to_vec();
    for (id, span) in source.span_indices().enumerate() {
        dec.charge_work(1)?;
        let domain =
            std::array::from_fn(|axis| [&knots[axis][span[axis]], &knots[axis][span[axis] + 1]]);
        let carrier = source.bezier_patch(span, domain, dec)?;
        let patch = RationalPatch::new(carrier.clone(), dec)?;
        if !patch.carrier().same_definition(&carrier) {
            return Err(Refusal::new(
                "boolean.numerical-cylinder-source",
                "quarter admission changed its retained native carrier",
            ));
        }
        let key = format!("cylinder/quarter-{id}");
        checks.push(NumericalCheck::Passed {
            name: "native.cylinder-retained-quarter-embedding".into(),
            scope: key.clone(),
        });
        quarters.push(NumericalCylinderQuarter {
            id,
            key,
            surface: Surface::BSpline(SplineSupport::RationalPatch(patch)),
            axis_parameters_mm: parameters.clone(),
        });
    }
    if quarters.len() != 4 {
        return Err(Refusal::new(
            "boolean.numerical-cylinder-source",
            "retained cylinder requires exactly four native quarters",
        ));
    }
    dec.charge_work(checks.len() as u64)?;
    let evidence =
        NumericalEvidence::new(retained.evidence().tolerance_mm().clone(), checks, vec![])?;
    dec.check_deadline()?;
    let finished = catallax_numbers::current_operation_work().used;
    Ok(NumericalCylinderSource {
        original: retained.original().clone(),
        source: source.clone(),
        axis_parameters_mm: parameters,
        quarters,
        evidence,
        source_work: after_source - before,
        chart_work: finished - after_source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Real, with_operation_budget};
    use catallax_shapes::{Cone, Point3, Vector3};
    fn q(i: i64) -> Rational {
        Rational::from_i64(i)
    }
    fn options() -> NumericalOptions {
        NumericalOptions {
            tolerance_mm: Rational::from_ratio(1, 10_000_000).unwrap(),
        }
    }
    fn cylinder(skew: bool) -> Cylinder {
        Cylinder {
            origin: Point3::new(
                Real::from(Rational::from_ratio(21, 2).unwrap()),
                12.into(),
                0.into(),
            ),
            axis: if skew {
                Vector3::new(1.into(), (-1).into(), 1.into())
            } else {
                Vector3::new(0.into(), (-1).into(), 0.into())
            },
            radius_sq: 1.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        }
    }
    #[test]
    fn standalone_cylinder_sources_match_finite_queries_in_native_parameters() {
        for skew in [false, true] {
            let cylinder = cylinder(skew);
            let parameters = [
                -Rational::from_ratio(99, 2).unwrap(),
                Rational::from_ratio(99, 2).unwrap(),
            ];
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let source = NumericalCylinderSource::new(
                    &cylinder,
                    parameters.clone(),
                    options(),
                    &mut dec,
                )
                .unwrap();
                assert_eq!(
                    source.source_work + source.chart_work,
                    catallax_numbers::current_operation_work().used
                );
                assert_eq!(source.signed_mm_at_v(&q(0)).unwrap(), parameters[0]);
                assert_eq!(source.signed_mm_at_v(&q(1)).unwrap(), parameters[1]);
                assert_eq!(
                    source
                        .signed_mm_at_v(&Rational::from_ratio(1, 2).unwrap())
                        .unwrap(),
                    q(0)
                );
                for upper in [false, true] {
                    let cone = Surface::Cone(Cone {
                        apex: Point3::new(
                            0.into(),
                            if upper { 23.into() } else { (-11).into() },
                            0.into(),
                        ),
                        axis: Vector3::new(
                            0.into(),
                            if upper { (-1).into() } else { 1.into() },
                            0.into(),
                        ),
                        slope: 1.into(),
                        ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
                    });
                    let tool = Surface::Cylinder(cylinder.clone());
                    let query = catallax_intersect::numerical_native_finite_intersection(
                        [&cone, &tool],
                        [q(11), q(12)],
                        parameters.clone(),
                        &mut dec,
                    )
                    .unwrap()
                    .unwrap();
                    assert!(source.source.same_definition(&query.windows().sources()[1]));
                }
                for (id, quarter) in source.quarters.iter().enumerate() {
                    assert_eq!(quarter.id, id);
                    assert_eq!(quarter.key, format!("cylinder/quarter-{id}"));
                    assert_eq!(quarter.axis_parameters_mm, parameters);
                    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = &quarter.surface
                    else {
                        unreachable!()
                    };
                    assert_eq!(
                        patch.carrier().domain(),
                        [[&q(id as i64), &q(id as i64 + 1)], [&q(0), &q(1)]]
                    );
                    assert!(patch.carrier().has_rational_definition());
                }
                eprintln!(
                    "canonical cylinder skew={skew}: {} work",
                    source.source_work + source.chart_work
                );
            });
            assert!(!work.exhausted, "{work:?}");
        }
    }
    #[test]
    fn canonical_source_rejects_invalid_domains_and_honors_the_caller_limits() {
        let cylinder = cylinder(false);
        for parameters in [[q(1), q(1)], [q(2), q(1)]] {
            assert!(
                NumericalCylinderSource::new(
                    &cylinder,
                    parameters,
                    options(),
                    &mut Decider::default()
                )
                .is_err()
            );
        }
        let (result, work) = with_operation_budget(0, || {
            NumericalCylinderSource::new(
                &cylinder,
                [q(-49), q(49)],
                options(),
                &mut Decider::default(),
            )
        });
        assert!(result.is_err());
        assert!(work.exhausted);
        let mut options = options();
        options.tolerance_mm = q(0);
        assert!(
            NumericalCylinderSource::new(
                &cylinder,
                [q(-49), q(49)],
                options,
                &mut Decider::default()
            )
            .is_err()
        );
    }
}
