#![allow(missing_docs)]
//! Public radius-one rational cylinder; no donor execution or saved expectations.
use catallax_trace::{
    extent::Extent2,
    global_solver::SolverRefusal,
    spline::{BSplineCurve, BSplineSurface, ControlPoint},
    surface_curves::{CurveKind, IntersectionCurve},
    surface_eval::evaluate_point,
    surface_intersector::Intersector,
    vector::{Vector2, Vector3},
};

fn cylinder(shift: f64, closed_v: bool, reverse: bool) -> BSplineSurface {
    let r = 0.5_f64.sqrt();
    let ring = [
        (1., 0., 1.),
        (1., 1., r),
        (0., 1., 1.),
        (-1., 1., r),
        (-1., 0., 1.),
        (-1., -1., r),
        (0., -1., 1.),
        (1., -1., r),
        (1., 0., 1.),
    ];
    let zs = if reverse { [2., 0.] } else { [0., 2.] };
    let point = |i: usize, j: usize| {
        let (x, y, w) = ring[i];
        ControlPoint {
            coordinates: [x * w, y * w, zs[j] * w, w],
        }
    };
    let points = if closed_v {
        (0..2)
            .flat_map(|j| (0..9).map(move |i| point(i, j)))
            .collect()
    } else {
        (0..9)
            .flat_map(|i| (0..2).map(move |j| point(i, j)))
            .collect()
    };
    let periodic = vec![0., 0., 0., 1., 1., 2., 2., 3., 3., 4., 4., 4.]
        .into_iter()
        .map(|k| k + shift)
        .collect();
    let axial = vec![0., 0., 2., 2.];
    let (degrees, knots, count) = if closed_v {
        ([1, 2], [axial, periodic], [2, 9])
    } else {
        ([2, 1], [periodic, axial], [9, 2])
    };
    BSplineSurface::new(degrees, knots, count, points).unwrap()
}

#[derive(Clone, Copy, Debug)]
struct Case {
    closed_v: bool,
    reverse: bool,
    swap: bool,
    trim: usize,
    shift: f64,
}

fn domain(case: Case) -> Extent2 {
    let (lo, hi) = [(0., 4.), (0., 2.), (2., 4.), (0.25, 3.75)][case.trim];
    if case.closed_v {
        Extent2::new(
            Vector2::new(0.25, lo + case.shift),
            Vector2::new(1.75, hi + case.shift),
        )
    } else {
        Extent2::new(
            Vector2::new(lo + case.shift, 0.25),
            Vector2::new(hi + case.shift, 1.75),
        )
    }
}

fn curve_bits(curve: &BSplineCurve, out: &mut Vec<u64>) {
    out.extend([
        curve.degree() as u64,
        curve.dimension() as u64,
        curve.knots().len() as u64,
        curve.control_points().len() as u64,
    ]);
    out.extend(curve.knots().iter().map(|x| x.to_bits()));
    out.extend(
        curve
            .control_points()
            .iter()
            .flat_map(|p| p.coordinates.map(f64::to_bits)),
    );
}

fn complete_bits(curves: &[IntersectionCurve]) -> Vec<u64> {
    let mut out = vec![curves.len() as u64];
    for curve in curves {
        out.extend([curve.kind as u64, curve.deviation.to_bits()]);
        curve_bits(&curve.curve, &mut out);
        for uv in &curve.uv_curves {
            out.push(u64::from(uv.is_some()));
            if let Some(uv) = uv {
                curve_bits(uv, &mut out);
            }
        }
    }
    out
}

#[derive(Debug, PartialEq)]
struct Outcome {
    raw: Vec<u64>,
    normalized: Vec<u64>,
    refusals: Vec<String>,
}

fn run(case: Case) -> Outcome {
    let c = cylinder(case.shift, case.closed_v, case.reverse);
    let p = BSplineSurface::canonical_plane(
        Vector3::new(1., 0., 0.),
        Vector3::new(0., 1., 0.),
        Vector3::new(0., 0., 1.),
    )
    .unwrap();
    // Both operands clip the seam's axial range; the plane is more restrictive.
    let cd = domain(case);
    let pd = Extent2::new(Vector2::new(-1., 0.5), Vector2::new(1., 1.5));
    let intersector = if case.swap {
        Intersector::new(&p, pd, &c, cd)
    } else {
        Intersector::new(&c, cd, &p, pd)
    };
    assert_eq!(
        intersector.domains,
        if case.swap { [pd, cd] } else { [cd, pd] }
    );
    let seam = c.iso_curve(usize::from(case.closed_v), case.shift).unwrap();
    let refusals: Vec<_> = [false, true]
        .into_iter()
        .map(|on_surface| {
            let result = catallax_trace::surface_drop_curve::drop_and_trim_curve(
                &c,
                c.natural_domain(),
                &seam,
                seam.natural_interval(),
                1e-3,
                true,
                on_surface,
            );
            assert!(matches!(result, Err(SolverRefusal::NotPorted(_))));
            format!("{:?}", result.unwrap_err())
        })
        .collect();
    let curves = intersector
        .do_intersection()
        .unwrap_or_else(|e| panic!("{case:?}: {e:?}"));
    if case.trim == 3 {
        assert!(curves.is_empty(), "excluded seam: {case:?}");
        return Outcome {
            raw: complete_bits(&curves),
            normalized: complete_bits(&curves),
            refusals,
        };
    }
    assert_eq!(curves.len(), 1, "duplicate/lost seam: {case:?}");
    assert_eq!(curves[0].kind, CurveKind::Tangent);
    let curve = &curves[0];
    let iv = curve.curve.natural_interval();
    let endpoints = [iv.min, iv.max].map(|t| curve.curve.evaluate_point(t).unwrap().z);
    assert!((endpoints[0].min(endpoints[1]) - 0.5).abs() <= 1e-12);
    assert!((endpoints[0].max(endpoints[1]) - 1.5).abs() <= 1e-12);
    assert_eq!(
        endpoints[1] > endpoints[0],
        !case.reverse,
        "trace direction: {case:?}"
    );
    for s in 0..2 {
        for point in curve.uv_curves[s].as_ref().unwrap().control_points() {
            let [u, v, _, _] = point.coordinates;
            assert!(intersector.domains[s].contains_point(Vector2::new(u, v), 0.0));
        }
    }
    for step in 0..=16 {
        let t = iv.evaluate(f64::from(step) / 16.);
        let point = curve.curve.evaluate_point(t).unwrap();
        assert!((point.x - 1.).abs() <= 1e-12 && point.y.abs() <= 1e-12);
        assert!((0.5 - 1e-12..=1.5 + 1e-12).contains(&point.z));
        for s in 0..2 {
            let uv = curve.uv_curves[s]
                .as_ref()
                .unwrap()
                .evaluate_point(t)
                .unwrap();
            let uv = Vector2::new(uv.x, uv.y);
            assert!(
                intersector.domains[s].contains_point(uv, 1e-12),
                "{case:?} {uv:?}"
            );
            let image = evaluate_point(intersector.surfaces[s], uv).unwrap();
            assert!(point.distance(image) <= 1e-12);
        }
    }
    // Normalize only the periodic coordinate of the owner's UV control points.
    // Everything else, including ordering, knots, weights and signed zeros, is compared by bits.
    let mut normalized = curves.clone();
    let owner = usize::from(case.swap);
    let uv = normalized[0].uv_curves[owner].as_mut().unwrap();
    let mut points = uv.control_points().to_vec();
    for point in &mut points {
        point.coordinates[usize::from(case.closed_v)] -= case.shift;
    }
    *uv = BSplineCurve::new(uv.degree(), uv.knots().to_vec(), points).unwrap();
    Outcome {
        raw: complete_bits(&curves),
        normalized: complete_bits(&normalized),
        refusals,
    }
}

fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    for closed_v in [false, true] {
        for reverse in [false, true] {
            for swap in [false, true] {
                for trim in 0..4 {
                    for shift in [0., 4., -4., 8.] {
                        out.push(Case {
                            closed_v,
                            reverse,
                            swap,
                            trim,
                            shift,
                        });
                    }
                }
            }
        }
    }
    out
}

#[test]
fn seam_ranges_are_complete_clipped_and_period_translated() {
    let cases = cases();
    for family in cases.chunks_exact(4) {
        let expected = run(family[0]).normalized;
        for &case in &family[1..] {
            assert_eq!(run(case).normalized, expected, "{case:?}");
        }
    }
}

#[test]
fn complete_outcomes_match_at_actual_one_two_and_four_workers() {
    let cases = cases();
    let mut expected = None;
    for workers in [1, 2, 4] {
        let barrier = std::sync::Barrier::new(workers);
        let mut rows = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|worker| {
                    let cases = &cases;
                    let barrier = &barrier;
                    scope.spawn(move || {
                        barrier.wait();
                        cases
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| i % workers == worker)
                            .map(|(i, &case)| (i, run(case)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|handle| handle.join().unwrap())
                .collect::<Vec<_>>()
        });
        rows.sort_by_key(|(i, _)| *i);
        assert_eq!(rows.len(), cases.len());
        if let Some(expected) = &expected {
            assert_eq!(&rows, expected);
        } else {
            expected = Some(rows);
        }
        println!(
            "B02 workers={workers} cases={} partitions={:?} complete_outcomes=identical; application_work=unmetered",
            cases.len(),
            (0..workers)
                .map(|w| cases
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| i % workers == w)
                    .count())
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn arbitrary_seam_drop_keeps_its_typed_refusal() {
    let c = cylinder(0., false, false);
    let side = c.iso_curve(0, 0.).unwrap();
    for on_surface in [false, true] {
        let answer = catallax_trace::surface_drop_curve::drop_and_trim_curve(
            &c,
            c.natural_domain(),
            &side,
            side.natural_interval(),
            1e-3,
            true,
            on_surface,
        );
        assert!(matches!(answer, Err(SolverRefusal::NotPorted(_))));
    }
}
