//! Retained numerical fits and finite sampled incidence checks.

use super::*;
use catallax_shapes::BSplineCurve;
use catallax_trace::{
    global_solver::SolverRefusal,
    spline::{BSplineCurve as FloatCurve, PointForm},
    surface_curves::{CurveKind, IntersectionCurve},
    surface_intersector::Intersector,
    surface_work::WorkReport,
};
use std::collections::BTreeSet;

/// One fitted piece, with its native parameter and both retained UV definitions.
/// Its kind is the numerical tracer's classification, not an exact model fact.
#[derive(Clone, Debug)]
pub struct NumericalIntersectionPiece {
    curve: BSplineCurve,
    uv: [BSplineCurve; 2],
    kind: CurveKind,
}

impl NumericalIntersectionPiece {
    /// Retained 3D fit over its original native domain.
    pub fn curve(&self) -> &BSplineCurve {
        &self.curve
    }
    /// Retained UV fits in the source operand order and the same native domain.
    pub fn uv(&self) -> &[BSplineCurve; 2] {
        &self.uv
    }
    /// The tracer's numerical meeting classification.
    pub fn kind(&self) -> CurveKind {
        self.kind
    }
}

/// Every piece returned by one finite numerical search, with actual checks.
/// Sampling establishes no complete branch-count or whole-curve error proof.
#[derive(Clone, Debug)]
pub struct NumericalNativeIntersection {
    windows: NumericalNativeWindows,
    pieces: Vec<NumericalIntersectionPiece>,
    evidence: NumericalEvidence,
    trace_work: WorkReport,
    work_stages: NumericalNativeWork,
}

/// Actual operation-work deltas, including callback deadline checks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NumericalNativeWork {
    /// Analytic frames, finite windows and rounded source checks.
    pub windows: u64,
    /// Float conversion and intersector setup.
    pub setup: u64,
    /// Floating search, stepping and output, including nested helpers.
    pub tracing: u64,
    /// Retained definitions and independent sample/evidence checks.
    pub fit_checks: u64,
}

impl NumericalNativeIntersection {
    /// Original analytic carriers and rounded native source definitions.
    pub fn windows(&self) -> &NumericalNativeWindows {
        &self.windows
    }
    /// All reported pieces, without candidate-seed selection or deduplication.
    pub fn pieces(&self) -> &[NumericalIntersectionPiece] {
        &self.pieces
    }
    /// Source checks, finite sampled UV ownership and both incidence residuals.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Floating-trace checkpoint visits. Callback overhead and exact validation
    /// contribute additional work to the same operation allowance.
    pub fn trace_work(&self) -> WorkReport {
        self.trace_work
    }
    /// Measured work for the complete query, excluding any caller prefix.
    pub fn work_stages(&self) -> NumericalNativeWork {
        self.work_stages
    }
}

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.numerical-native-fit", detail)
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value)
        .ok_or_else(|| invalid("a retained fit contains a nonfinite value"))
}

fn retained(
    curve: &FloatCurve,
    dimension: usize,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    if dimension == 3 && curve.dimension() != 3 {
        return Err(invalid(
            "a retained 3D fit has the wrong coordinate dimension",
        ));
    }
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for i in 0..curve.control_points().len() {
        dec.charge_work(1)?;
        let (point, weight) = curve
            .control_point(i, PointForm::Euclidean)
            .map_err(|fault| invalid(format!("invalid fitted control: {fault:?}")))?;
        if dimension == 2 && point.z != 0.0 {
            return Err(invalid("a retained UV fit has a nonzero third coordinate"));
        }
        controls.push(Point3::new(
            exact(point.x)?.into(),
            exact(point.y)?.into(),
            if dimension == 2 {
                Real::zero()
            } else {
                exact(point.z)?.into()
            },
        ));
        weights.push(exact(weight)?.into());
    }
    BSplineCurve::new(
        curve.degree(),
        curve
            .knots()
            .iter()
            .copied()
            .map(exact)
            .collect::<Result<Vec<_>, _>>()?,
        controls,
        weights,
        dec,
    )
}

struct Carrier {
    origin: [Ball; 3],
    axis: [Ball; 3],
    radius: Option<Ball>,
    slope: Option<Ball>,
}

impl Carrier {
    fn new(surface: &Surface, dec: &mut Decider) -> Result<Self, Refusal> {
        let (origin, axis, radius, slope) = match surface {
            Surface::Cone(cone) => (&cone.apex, &cone.axis, None, Some(cone.slope.clone())),
            Surface::Cylinder(cylinder) => (
                &cylinder.origin,
                &cylinder.axis,
                Some(cylinder.radius()),
                None,
            ),
            _ => return Err(invalid("unexpected native carrier family")),
        };
        let axis = unit(axis, dec)?;
        let mut enclose = |values: [&Real; 3]| -> Result<[Ball; 3], Refusal> {
            Ok([
                enclosure(values[0], dec)?,
                enclosure(values[1], dec)?,
                enclosure(values[2], dec)?,
            ])
        };
        let origin = enclose([&origin.x, &origin.y, &origin.z])?;
        let axis = enclose([&axis.x, &axis.y, &axis.z])?;
        Ok(Self {
            origin,
            axis,
            radius: radius.map(|v| enclosure(&v, dec)).transpose()?,
            slope: slope.map(|v| enclosure(&v, dec)).transpose()?,
        })
    }

    fn residual(&self, point: &[Ball; 3]) -> Result<Rational, Refusal> {
        let offset: [Ball; 3] = std::array::from_fn(|k| &point[k] - &self.origin[k]);
        let height = (&offset[0] * &self.axis[0])
            + (&offset[1] * &self.axis[1])
            + (&offset[2] * &self.axis[2]);
        let radial: [Ball; 3] = std::array::from_fn(|k| &offset[k] - &(&height * &self.axis[k]));
        let distance = (radial[0].square() + radial[1].square() + radial[2].square())
            .sqrt()
            .ok_or_else(|| invalid("a sampled radial distance was unresolved"))?;
        let residual = if let Some(radius) = &self.radius {
            &distance - radius
        } else {
            if height.lo().sign() != Sign::Positive {
                return Err(invalid(
                    "a sampled fit reaches the cone apex or opposite nappe",
                ));
            }
            &distance - &(&height * self.slope.as_ref().expect("cone slope"))
        };
        Ok(absolute_upper(&residual))
    }
}

fn check_piece(
    piece: &NumericalIntersectionPiece,
    sources: &[BSplineSurface; 2],
    carriers: &[Carrier; 2],
    index: usize,
    dec: &mut Decider,
) -> Result<Vec<NumericalCheck>, Refusal> {
    for uv in &piece.uv {
        if uv.domain() != piece.curve.domain() {
            return Err(invalid("a UV fit changed the native 3D curve interval"));
        }
    }
    let mut parameters = BTreeSet::new();
    for (_, span) in piece.curve.spans() {
        for i in 0..=4 {
            parameters
                .insert(span[0] + &((span[1] - span[0]) * Rational::from_ratio(i, 4).unwrap()));
        }
    }
    let mut coupled = [Vec::new(), Vec::new()];
    let mut analytic = [Vec::new(), Vec::new()];
    for t in parameters {
        dec.charge_work(1)?;
        let point = piece.curve.evaluate_box([&t, &t], BITS, dec)?;
        for k in 0..2 {
            // These rational UV definitions can be evaluated exactly at the
            // rational sample. No parameter clamp or boundary projection occurs.
            let at = piece.uv[k].point_at(&Real::from(t.clone()), dec)?;
            let u =
                at.x.rational_ref()
                    .ok_or_else(|| invalid("a retained UV coefficient lost rational identity"))?;
            let v =
                at.y.rational_ref()
                    .ok_or_else(|| invalid("a retained UV coefficient lost rational identity"))?;
            for (parameter, domain) in [u, v].into_iter().zip(sources[k].domain()) {
                if parameter < domain[0] || parameter > domain[1] {
                    return Err(invalid("a sampled UV fit leaves its finite source domain"));
                }
            }
            let on = sources[k].evaluate_box([u, u], [v, v], BITS, dec)?;
            let delta: [Ball; 3] = std::array::from_fn(|j| &point[j] - &on[j]);
            let residual = (delta[0].square() + delta[1].square() + delta[2].square())
                .sqrt()
                .ok_or_else(|| invalid("a coupled sampled distance was unresolved"))?
                .hi();
            let location = vec![t.clone(), u.clone(), v.clone()];
            coupled[k].push(NumericalSample {
                parameters: location.clone(),
                distance_upper_mm: residual,
            });
            analytic[k].push(NumericalSample {
                parameters: location,
                distance_upper_mm: carriers[k].residual(&point)?,
            });
        }
    }
    let mut checks = Vec::new();
    for (k, (coupled, analytic)) in coupled.into_iter().zip(analytic).enumerate() {
        let scope = format!("piece-{index}/source-{k}");
        checks.push(NumericalCheck::Samples {
            name: "native.fit-sampled-coupled-residual".into(),
            scope: scope.clone(),
            samples: coupled,
        });
        checks.push(NumericalCheck::Samples {
            name: "native.fit-sampled-analytic-carrier-residual".into(),
            scope,
            samples: analytic,
        });
    }
    Ok(checks)
}

/// Traces both full native windows at a preset 1e-7 mm tolerance. Every returned
/// 3D/UV piece is retained and checked at endpoints and quarter-span locations.
/// No exact incidence, complete cover proof or Boolean topology is published.
pub fn numerical_native_intersection(
    surfaces: [&Surface; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeIntersection>, Refusal> {
    if operation_budget_active() {
        return build_intersection(surfaces, dec);
    }
    let (result, work) = with_operation_budget(20_000_000, || build_intersection(surfaces, dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "numerical native intersection exhausted its operation allowance",
        ));
    }
    result
}

pub(super) fn pay_trace_work(units: u64, deadline: Option<std::time::Instant>) -> bool {
    let mut meter = Decider::default();
    meter.set_deadline(deadline);
    meter.charge_work(units).is_ok()
}

fn trace_with_caller(
    intersector: &Intersector<'_>,
    dec: &Decider,
) -> (Result<Vec<IntersectionCurve>, SolverRefusal>, WorkReport) {
    intersector.do_intersection_with_work(dec.deadline(), pay_trace_work)
}

fn build_intersection(
    surfaces: [&Surface; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeIntersection>, Refusal> {
    let started = catallax_numbers::current_operation_work().used;
    let Some(windows) = numerical_native_windows(surfaces, dec)? else {
        return Ok(None);
    };
    let after_windows = catallax_numbers::current_operation_work().used;
    let a =
        crate::pair::trace_candidates::surface(&windows.sources[0], &Isometry::identity(), dec)?;
    let b =
        crate::pair::trace_candidates::surface(&windows.sources[1], &Isometry::identity(), dec)?;
    dec.charge_work(1)?;
    let mut intersector = Intersector::new(&a, a.natural_domain(), &b, b.natural_domain());
    intersector.approx_tolerance = windows.evidence.tolerance_mm().to_f64_nearest() / 4.;
    intersector.check_output_spans = true;
    let after_setup = catallax_numbers::current_operation_work().used;
    let (fits, trace_work) = trace_with_caller(&intersector, dec);
    let after_trace = catallax_numbers::current_operation_work().used;
    if trace_work.cancelled {
        dec.charge_work(0)?;
        return Err(invalid(
            "numerical native tracing cancelled without a retained result",
        ));
    }
    let fits = fits.map_err(|fault| invalid(format!("native numerical trace: {fault:?}")))?;
    let carriers = [
        Carrier::new(&windows.originals[0], dec)?,
        Carrier::new(&windows.originals[1], dec)?,
    ];
    let mut pieces = Vec::new();
    let mut checks = windows.evidence.checks().to_vec();
    for (index, fit) in fits.into_iter().enumerate() {
        let [Some(a), Some(b)] = fit.uv_curves else {
            return Err(invalid(
                "a reported piece needs both retained UV definitions",
            ));
        };
        let piece = NumericalIntersectionPiece {
            curve: retained(&fit.curve, 3, dec)?,
            uv: [retained(&a, 2, dec)?, retained(&b, 2, dec)?],
            kind: fit.kind,
        };
        checks.extend(check_piece(
            &piece,
            &windows.sources,
            &carriers,
            index,
            dec,
        )?);
        pieces.push(piece);
    }
    checks.push(NumericalCheck::Passed {
        name: "native.retained-every-reported-piece".into(),
        scope: "finite-window-search".into(),
    });
    let samples: usize = checks
        .iter()
        .map(|check| match check {
            NumericalCheck::Samples { samples, .. } => samples.len(),
            _ => 1,
        })
        .sum();
    dec.charge_work(samples as u64)?;
    let evidence = NumericalEvidence::new(windows.evidence.tolerance_mm().clone(), checks, vec![])?;
    let finished = catallax_numbers::current_operation_work().used;
    let work_stages = NumericalNativeWork {
        windows: after_windows - started,
        setup: after_setup - after_windows,
        tracing: after_trace - after_setup,
        fit_checks: finished - after_trace,
    };
    Ok(Some(NumericalNativeIntersection {
        windows,
        pieces,
        evidence,
        trace_work,
        work_stages,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axial() -> [Surface; 2] {
        [
            Surface::Cone(Cone {
                apex: Point3::new(s(0), s(-11), s(0)),
                axis: Vector3::new(s(0), s(1), s(0)),
                slope: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
            Surface::Cylinder(Cylinder {
                origin: Point3::new(s(10), s(12), s(0)),
                axis: Vector3::new(s(0), s(-1), s(0)),
                radius_sq: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
        ]
    }

    #[test]
    fn numerical_native_query_keeps_all_fits_and_checks_both_carriers() {
        for reverse in [false, true] {
            let surfaces = axial();
            let ordered = if reverse {
                [&surfaces[1], &surfaces[0]]
            } else {
                [&surfaces[0], &surfaces[1]]
            };
            let (result, work) = with_operation_budget(20_000_000, || {
                numerical_native_intersection(ordered, &mut Decider::default())
            });
            let result = result.unwrap().unwrap();
            assert!(!work.exhausted, "{work:?}");
            assert_eq!(result.pieces().len(), 2);
            assert!(!result.trace_work().cancelled);
            assert!(result.trace_work().accepted > 1000);
            assert_eq!(
                result.evidence().tolerance_mm(),
                &Rational::from_ratio(1, 10_000_000).unwrap()
            );
            assert_eq!(
                matches!(result.windows().originals()[0], Surface::Cylinder(_)),
                reverse
            );
            for piece in result.pieces() {
                assert_eq!(piece.curve().domain(), piece.uv()[0].domain());
                assert_eq!(piece.curve().domain(), piece.uv()[1].domain());
            }
            assert_eq!(result.evidence().checks().len(), 5 + 4 * 2 + 1);
            eprintln!(
                "numerical native query reverse={reverse} total_work={} trace={:?}",
                work.used,
                result.trace_work()
            );
            let (negative, work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let mut piece = result.pieces()[0].clone();
                let delta = Vector3::new(
                    Real::from(Rational::from_ratio(1, 1000).unwrap()),
                    s(0),
                    s(0),
                );
                piece.curve = BSplineCurve::new(
                    piece.curve.degree(),
                    piece.curve.knots().to_vec(),
                    piece.curve.controls().map(|p| p.plus(&delta)).collect(),
                    piece.curve.weights().cloned().collect(),
                    &mut dec,
                )
                .unwrap();
                let carriers = [
                    Carrier::new(&surfaces[usize::from(reverse)], &mut dec).unwrap(),
                    Carrier::new(&surfaces[usize::from(!reverse)], &mut dec).unwrap(),
                ];
                let checks =
                    check_piece(&piece, result.windows.sources(), &carriers, 0, &mut dec).unwrap();
                NumericalEvidence::new(result.evidence().tolerance_mm().clone(), checks, vec![])
            });
            assert_eq!(
                negative.unwrap_err().code,
                "numbers.numerical-tolerance-exceeded"
            );
            assert!(!work.exhausted, "{work:?}");
        }
    }

    #[test]
    fn a_caller_deadline_expires_during_native_tracing_without_spending_the_work_limit() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let surfaces = axial();
            let mut dec = Decider::default();
            let windows = numerical_native_windows([&surfaces[0], &surfaces[1]], &mut dec)
                .unwrap()
                .unwrap();
            let a = crate::pair::trace_candidates::surface(
                &windows.sources[0],
                &Isometry::identity(),
                &mut dec,
            )
            .unwrap();
            let b = crate::pair::trace_candidates::surface(
                &windows.sources[1],
                &Isometry::identity(),
                &mut dec,
            )
            .unwrap();
            let mut intersector = Intersector::new(&a, a.natural_domain(), &b, b.natural_domain());
            intersector.approx_tolerance = 2.5e-8;
            intersector.check_output_spans = true;
            let (baseline, complete) = trace_with_caller(&intersector, &dec);
            let pieces = baseline.unwrap().len();
            let callbacks = complete.accepted - complete.source_items + 2;
            let clock = catallax_numbers::TestDeadline::after(
                catallax_numbers::TestDeadlinePoint::Enclosure,
                usize::try_from(callbacks / 2).unwrap(),
            );
            dec.set_deadline(Some(clock.instant()));
            let (expired, cancelled) = trace_with_caller(&intersector, &dec);
            eprintln!("mid-trace deadline report: {cancelled:?}");
            assert!(
                expired.is_err(),
                "trace retained output after caller expiry: {cancelled:?}"
            );
            assert!(cancelled.cancelled);
            assert!(cancelled.search_branches > 0 && cancelled.trace_visits > 0);
            assert!(cancelled.accepted < complete.accepted);
            assert!(!catallax_numbers::current_operation_work().exhausted);
            clock.assert_expired();
            assert_eq!(dec.charge_work(0).unwrap_err().code, "numbers.deadline");
            dec.set_deadline(None);
            let (recovered, report) = trace_with_caller(&intersector, &dec);
            assert_eq!(recovered.unwrap().len(), pieces);
            assert_eq!(report, complete);
        });
        assert!(!work.exhausted, "{work:?}");
    }

    #[test]
    fn moved_rotated_and_scaled_native_queries_keep_the_preset_tolerance() {
        for scale in [Rational::from_ratio(1, 2).unwrap(), r(3)] {
            let surfaces = [
                Surface::Cone(Cone {
                    apex: Point3::new(Real::from(r(3) + &r(11) * &scale), s(-7), s(2)),
                    axis: Vector3::new(s(-1), s(0), s(0)),
                    slope: s(1),
                    ref_dir: Vector3::new(s(0), s(1), s(0)),
                }),
                Surface::Cylinder(Cylinder {
                    origin: Point3::new(
                        Real::from(r(3) - &r(12) * &scale),
                        Real::from(r(-7) + &r(10) * &scale),
                        s(2),
                    ),
                    axis: Vector3::new(s(1), s(0), s(0)),
                    radius_sq: Real::from(&scale * &scale),
                    ref_dir: Vector3::new(s(0), s(1), s(0)),
                }),
            ];
            let (result, work) = with_operation_budget(20_000_000, || {
                numerical_native_intersection([&surfaces[0], &surfaces[1]], &mut Decider::default())
            });
            let result = result.unwrap().unwrap();
            assert!(!work.exhausted, "scale={scale}: {work:?}");
            assert_eq!(result.pieces().len(), 2);
            assert_eq!(
                result.evidence().tolerance_mm(),
                &Rational::from_ratio(1, 10_000_000).unwrap()
            );
            let stages = result.work_stages();
            assert_eq!(
                stages.windows + stages.setup + stages.tracing + stages.fit_checks,
                work.used
            );
            eprintln!(
                "moved/rotated native query scale={scale} total_work={} stages={stages:?}",
                work.used
            );
        }
    }
}
