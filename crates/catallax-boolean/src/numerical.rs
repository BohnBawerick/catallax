//! Numerical sources for open face consumers, separate from exact incidence.

use catallax_numbers::{Decider, NumericalCheck, NumericalEvidence, Rational, Refusal};
use catallax_shapes::{BSplineSurface, Cylinder, RationalPatch, SplineSupport, Surface};
use catallax_topology::numerical::NumericalOptions;

/// A canonical cylinder chart, independent of a cone face or assembly key.
#[derive(Clone, Debug)]
pub struct NumericalCylinderQuarter {
    /// Native quadrant, 0 through 3.
    id: usize,
    /// Deterministic descriptor-local key. The caller supplies assembly identities.
    key: String,
    /// Retained embedded patch over native u=[id,id+1], v=[0,1].
    surface: Surface,
    /// Signed millimetres at native v=0 and v=1 from the original tool origin.
    axis_parameters_mm: [Rational; 2],
}

impl NumericalCylinderQuarter {
    /// Native quadrant, 0 through 3.
    pub fn id(&self) -> usize {
        self.id
    }
    /// Descriptor-local key; callers provide their own assembly identity.
    pub fn key(&self) -> &str {
        &self.key
    }
    /// Retained patch over native u=[id,id+1], v=[0,1].
    pub fn surface(&self) -> &Surface {
        &self.surface
    }
    /// Signed millimetres at native v=0 and v=1 about the original origin.
    pub fn axis_parameters_mm(&self) -> &[Rational; 2] {
        &self.axis_parameters_mm
    }
}

/// A retained cylinder source and all four native quarter patches.
/// Constructing it supplies numerical source/control checks and patch embedding
/// checks. It does not construct a Boolean, admitted body or STEP result.
/// Construction keeps the retained definition, interval, quarters and evidence
/// together. Read-only accessors preserve that association.
///
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// use catallax_numbers::Rational;
/// fn change(mut source: NumericalCylinderSource) {
///     source.axis_parameters_mm[1] = Rational::one();
/// }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// fn change(mut source: NumericalCylinderSource) { source.quarters.clear(); }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// fn change(mut source: NumericalCylinderSource) {
///     source.evidence = source.evidence().clone();
/// }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderQuarter;
/// fn change(mut quarter: NumericalCylinderQuarter) { quarter.key.clear(); }
/// ```
#[derive(Clone, Debug)]
pub struct NumericalCylinderSource {
    /// Original analytic carrier, with the caller's origin, axis and seam.
    original: Cylinder,
    /// Complete finite numerical definition, without native knot changes.
    source: BSplineSurface,
    /// Signed millimetres at native v=0 and v=1 from the original origin.
    axis_parameters_mm: [Rational; 2],
    /// Canonical quarters in increasing native-u order.
    quarters: Vec<NumericalCylinderQuarter>,
    /// Preset tolerance and actual source/control/chart checks.
    evidence: NumericalEvidence,
    /// Actual work for source construction and its checks.
    source_work: u64,
    /// Actual work for admitting every retained quarter patch.
    chart_work: u64,
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
    /// Original analytic cylinder, including its authored origin and seam.
    pub fn original(&self) -> &Cylinder {
        &self.original
    }
    /// Complete retained finite source, without changing native knots.
    pub fn source(&self) -> &BSplineSurface {
        &self.source
    }
    /// Signed millimetres at native v=0 and v=1 about the original origin.
    pub fn axis_parameters_mm(&self) -> &[Rational; 2] {
        &self.axis_parameters_mm
    }
    /// All four retained quarters in increasing native-u order.
    pub fn quarters(&self) -> &[NumericalCylinderQuarter] {
        &self.quarters
    }
    /// Actual source/control/chart checks at the preset tolerance.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Measured source construction and check work.
    pub fn source_work(&self) -> u64 {
        self.source_work
    }
    /// Measured retained quarter admission and check work.
    pub fn chart_work(&self) -> u64 {
        self.chart_work
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
    fn cylinder(offset: &str, skew: bool) -> Cylinder {
        let axis = if skew {
            Vector3::new(1.into(), (-1).into(), 1.into())
                .scaled(&(&Real::one() / &Real::from(3).sqrt()))
        } else {
            Vector3::new(0.into(), (-1).into(), 0.into())
        };
        Cylinder {
            origin: Point3::new(Rational::parse(offset).unwrap().into(), 12.into(), 0.into()),
            axis,
            radius_sq: 1.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        }
    }
    #[test]
    fn standalone_cylinder_sources_match_finite_queries_in_native_parameters() {
        for (offset, skew, extent) in [
            ("10", false, "49"),
            ("10.5", false, "99/2"),
            ("10.5", true, "99/2"),
        ] {
            let cylinder = cylinder(offset, skew);
            let extent = Rational::parse(extent).unwrap();
            let parameters = [-&extent, extent];
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
                assert!(source.signed_mm_at_v(&q(-1)).is_err());
                assert!(source.signed_mm_at_v(&q(2)).is_err());
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
                    "canonical cylinder offset={offset} skew={skew}: {} work",
                    source.source_work + source.chart_work
                );
            });
            assert!(!work.exhausted, "{work:?}");
        }
    }
    #[test]
    fn existing_outer_wall_uses_its_original_origin_radius_and_axis_span() {
        let outer = Cylinder {
            origin: Point3::new(0.into(), 0.into(), 0.into()),
            axis: Vector3::new(0.into(), 1.into(), 0.into()),
            radius_sq: 144.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        };
        let mut dec = Decider::default();
        let descriptor =
            NumericalCylinderSource::new(&outer, [q(0), q(12)], options(), &mut dec).unwrap();
        for (v, y) in [(q(0), q(0)), (q(1), q(12))] {
            let point = descriptor
                .source
                .evaluate_box([&q(0), &q(0)], [&v, &v], 128, &dec)
                .unwrap();
            for (actual, expected) in point.iter().zip([q(12), y, q(0)]) {
                assert!(actual.lo() <= expected && actual.hi() >= expected)
            }
        }
        assert_eq!(descriptor.quarters.len(), 4);
        assert_eq!(
            descriptor
                .signed_mm_at_v(&Rational::from_ratio(1, 2).unwrap())
                .unwrap(),
            q(6)
        );
    }

    #[test]
    fn standalone_descriptor_is_identical_at_one_two_and_four_workers() {
        let jobs = [
            ("10", false, "49"),
            ("10.5", false, "99/2"),
            ("10.5", true, "99/2"),
        ];
        let build = |(offset, skew, extent): (&str, bool, &str)| {
            let extent = Rational::parse(extent).unwrap();
            NumericalCylinderSource::new(
                &cylinder(offset, skew),
                [-&extent, extent],
                options(),
                &mut Decider::default(),
            )
            .unwrap()
        };
        let baseline: Vec<_> = jobs.iter().copied().map(build).collect();
        for workers in [1, 2, 4] {
            let rows = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|rank| {
                        let jobs = &jobs;
                        let build = &build;
                        scope.spawn(move || {
                            jobs.iter()
                                .copied()
                                .enumerate()
                                .filter(|(i, _)| i % workers == rank)
                                .map(|(i, job)| (i, build(job)))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|handle| handle.join().unwrap())
                    .collect::<Vec<_>>()
            });
            for (i, actual) in rows {
                assert!(actual.source.same_definition(&baseline[i].source));
                assert_eq!(actual.evidence, baseline[i].evidence);
                assert_eq!(
                    (actual.source_work, actual.chart_work),
                    (baseline[i].source_work, baseline[i].chart_work)
                );
                for (a, b) in actual.quarters.iter().zip(&baseline[i].quarters) {
                    assert_eq!(a.key, b.key);
                    assert_eq!(a.axis_parameters_mm, b.axis_parameters_mm);
                    let Surface::BSpline(SplineSupport::RationalPatch(a)) = &a.surface else {
                        unreachable!()
                    };
                    let Surface::BSpline(SplineSupport::RationalPatch(b)) = &b.surface else {
                        unreachable!()
                    };
                    assert!(a.carrier().same_definition(b.carrier()));
                }
            }
        }
    }
    #[test]
    fn canonical_source_honors_an_expired_caller_deadline() {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            dec.set_deadline(Some(std::time::Instant::now()));
            NumericalCylinderSource::new(
                &cylinder("10", false),
                [q(-49), q(49)],
                options(),
                &mut dec,
            )
        });
        assert_eq!(result.unwrap_err().code, "numbers.deadline");
        assert!(!work.exhausted);
    }
    #[test]
    fn canonical_source_rejects_invalid_domains_and_honors_the_caller_limits() {
        let cylinder = cylinder("10", false);
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
