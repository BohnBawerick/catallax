//! Differential checks against the pinned strict x86 SMLib.
mod common;

use catallax_trace::{
    extent::Extent2,
    spline::{BSplineCurve, BSplineSurface, ControlPoint, PointForm},
};
use common::{Input, Output};

fn point(input: &mut Input) -> ControlPoint {
    ControlPoint {
        coordinates: [input.d(), input.d(), input.d(), input.d()],
    }
}

fn knots(values: &[f64], out: &mut Output) {
    let unique = catallax_trace::spline::unique_knots(values, None);
    out.i(unique.len() as i64);
    for (value, multiplicity) in unique {
        out.d(value);
        out.i(multiplicity as i64);
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "curve" => {
            let p = input.i() as usize;
            let n = input.i() as usize;
            let _dim = input.i();
            let u = (0..n + p + 1).map(|_| input.d()).collect();
            let points = (0..n).map(|_| point(input)).collect();
            let curve = BSplineCurve::new(p, u, points).unwrap();
            out.i(curve.degree() as i64);
            out.i(n as i64);
            out.i(curve.knots().len() as i64);
            out.flag(curve.is_rational());
            out.i(curve.dimension() as i64);
            let interval = curve.natural_interval();
            out.d(interval.min);
            out.d(interval.max);
            for (point, weight) in curve.ends() {
                out.v3(point);
                out.d(weight);
            }
            for i in 0..=n {
                out.d(curve.greville(i));
            }
            knots(curve.knots(), out);
            for form in [
                PointForm::NonRational,
                PointForm::Homogeneous,
                PointForm::Euclidean,
            ] {
                for i in 0..n {
                    if let Some((point, weight)) = out.status(curve.control_point(i, form)) {
                        out.v3(point);
                        out.d(weight);
                    }
                }
            }
            let tol = input.d();
            out.flag(curve.is_degenerate(tol));
            out.flag(curve.has_excess_knot_multiplicity());
            out.flag(curve.is_piecewise_c0());
            out.flag(curve.has_repeated_end_control_points(0.0));
            if let Some(bbox) = out.status(curve.bounding_box(interval)) {
                out.v3(bbox.min);
                out.v3(bbox.max);
            }
            if let Some(bbox) =
                out.status(curve.bounding_box(catallax_trace::extent::Extent1::new(
                    interval.evaluate(0.2),
                    interval.evaluate(0.8),
                )))
            {
                out.v3(bbox.min);
                out.v3(bbox.max);
            }
            let t = input.d();
            let der = input.i() as usize;
            let left = input.i() != 0;
            let nonzero = input.i() != 0;
            out.i(curve.find_span(t, left).unwrap() as i64);
            if let Some(values) = out.status(curve.evaluate(t, der, left, nonzero)) {
                for point in values {
                    out.v3(point);
                }
            }
            let split = curve.split(interval.mid());
            out.i(i64::from(split.is_err()));
            if let Ok(curves) = split {
                for c in curves {
                    out.i(c.degree() as i64);
                    out.i(c.control_points().len() as i64);
                    out.i(c.knots().len() as i64);
                    for &u in c.knots() {
                        out.d(u);
                    }
                    for cp in c.control_points() {
                        for &x in &cp.coordinates {
                            out.d(x);
                        }
                    }
                }
            }
        }
        "surface" => {
            let p = input.i() as usize;
            let q = input.i() as usize;
            let n = input.i() as usize;
            let m = input.i() as usize;
            let u = (0..n + p + 1).map(|_| input.d()).collect();
            let v = (0..m + q + 1).map(|_| input.d()).collect();
            let points = (0..n * m).map(|_| point(input)).collect();
            let surface = BSplineSurface::new([p, q], [u, v], [n, m], points).unwrap();
            for axis in 0..2 {
                out.i(surface.degrees()[axis] as i64);
                out.i(surface.counts()[axis] as i64);
                out.i(surface.knots()[axis].len() as i64);
            }
            out.flag(surface.is_rational());
            out.flag(surface.is_bilinear());
            out.flag(surface.has_excess_knot_multiplicity());
            out.flag(surface.has_repeated_end_control_points(0.0));
            let domain = surface.natural_domain();
            out.v2(domain.min);
            out.v2(domain.max);
            for axis in 0..2 {
                for i in 0..=surface.counts()[axis] {
                    out.d(surface.greville(axis, i).unwrap());
                }
            }
            for values in surface.knots() {
                knots(values, out);
            }
            for form in [
                PointForm::NonRational,
                PointForm::Homogeneous,
                PointForm::Euclidean,
            ] {
                for i in 0..n {
                    for j in 0..m {
                        if let Some((point, weight)) =
                            out.status(surface.control_point([i, j], form))
                        {
                            out.v3(point);
                            out.d(weight);
                        }
                    }
                }
            }
            let net = surface.euclidean_net([0, 0], [n - 1, m - 1]);
            out.i(i64::from(net.is_err()));
            if let Ok(points) = net {
                for point in points {
                    out.v3(point);
                }
            }
            let tol = input.d();
            out.flag(surface.is_degenerate(tol));
            let target = input.v3();
            let count = input.i() as usize;
            if let Some(extrema) = out.status(surface.net_extrema(target, count)) {
                out.i(extrema.nearest.len() as i64);
                out.i(extrema.farthest.len() as i64);
                for (distance, index) in extrema.nearest.iter().chain(&extrema.farthest) {
                    out.d(*distance);
                    out.i(index[0] as i64);
                    out.i(index[1] as i64);
                }
            }
            let low = input.v2();
            let high = input.v2();
            for index in surface
                .find_spans(low, [false, false])
                .unwrap()
                .into_iter()
                .chain(surface.find_spans(high, [true, true]).unwrap())
            {
                out.i(index as i64);
            }
            if let Some(partial) = out.status(surface.partial_net(Extent2::new(low, high))) {
                for row in partial.corners {
                    for point in row {
                        out.v3(point);
                    }
                }
                out.v3(partial.bounding_box.min);
                out.v3(partial.bounding_box.max);
                for value in partial.chord_heights.into_iter().chain(partial.angles) {
                    out.d(value);
                }
            }
            for axis in 0..2 {
                let t = if axis == 0 {
                    domain.u_interval().mid()
                } else {
                    domain.v_interval().mid()
                };
                let iv = if axis == 0 {
                    domain.u_interval()
                } else {
                    domain.v_interval()
                };
                for t in [
                    iv.min - 5e-7,
                    iv.min,
                    t,
                    iv.max,
                    iv.max + 5e-7,
                    iv.max + 2e-6,
                ] {
                    if let Some(iso) = out.status(surface.iso_curve(axis, t)) {
                        out.i(iso.degree() as i64);
                        out.i(iso.control_points().len() as i64);
                        for cp in iso.control_points() {
                            let (point, weight) = cp.get(PointForm::Homogeneous);
                            out.v3(point);
                            out.d(weight);
                        }
                    }
                }
            }
            for axis in 0..2 {
                let t = if axis == 0 {
                    domain.u_interval().mid()
                } else {
                    domain.v_interval().mid()
                };
                let split = surface.split(axis, t);
                out.i(i64::from(split.is_err()));
                if let Ok(surfaces) = split {
                    for s in surfaces {
                        for value in s
                            .degrees()
                            .into_iter()
                            .chain(s.counts())
                            .chain(s.knots().map(|v| v.len()))
                        {
                            out.i(value as i64);
                        }
                        for knots in s.knots() {
                            for &value in knots {
                                out.d(value);
                            }
                        }
                        for cp in s.control_points() {
                            for &value in &cp.coordinates {
                                out.d(value);
                            }
                        }
                    }
                }
            }
        }
        _ => panic!("unknown spline fixture operation {op}"),
    }
}

#[test]
fn spline_objects_match_strict_smlib() {
    common::run_gate(include_str!("fixtures/spline/cases.txt"), |_| 4, evaluate).assert_passed(&[
        "regular",
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

#[test]
fn spline_candidates_are_bit_identical_with_one_two_and_four_workers() {
    let fixture = include_str!("fixtures/spline/cases.txt");
    let budget = |_: &str| 4;
    let baseline = common::run_gate(fixture, budget, evaluate);
    assert!(baseline.failures.is_empty());
    let lines: Vec<_> = fixture.lines().collect();
    for workers in [2, 4] {
        let runs = std::thread::scope(|scope| {
            let jobs: Vec<_> = lines
                .chunks(lines.len().div_ceil(workers))
                .map(|chunk| {
                    let fixture = chunk.join("\n");
                    scope.spawn(move || common::run_gate(&fixture, budget, evaluate))
                })
                .collect();
            jobs.into_iter()
                .map(|job| job.join().unwrap())
                .collect::<Vec<_>>()
        });
        let mut bits = Vec::new();
        let mut cases = 0;
        for report in runs {
            assert!(report.failures.is_empty());
            cases += report.cases;
            bits.extend(report.candidate_bits);
        }
        assert_eq!(cases, baseline.cases);
        assert_eq!(
            bits, baseline.candidate_bits,
            "{workers} workers changed candidate bits"
        );
    }
}
