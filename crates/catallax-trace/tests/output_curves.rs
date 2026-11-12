//! Differential checks of the output-curve unit against the pinned strict x86
//! SMLib: `SmLine`, `SmHermiteCurve`, `SmBSplineCurve::CreateCanonical`,
//! `EditParameterization` and `Trim`, the IGES 126 curve helpers they use,
//! and `SmExtent2d::ClipLine2d`.
mod common;

use catallax_trace::{
    Fault,
    curve_edit::TrimRefusal,
    extent::{Extent1, Extent2},
    hermite::Hermite,
    iges_curve::{MAX_DEGREE, iges_126_curve},
    line::Line,
    pseudo_box::PseudoBox,
    spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z},
    vector::{Vector2, Vector3},
};
use common::{Input, Output, SUCCESS, status_code};

/// The zone tolerance SMLib's system default gives a free curve.
const ZONE_TOLERANCE: f64 = 1.0e-5;

fn point(input: &mut Input) -> ControlPoint {
    ControlPoint {
        coordinates: [input.d(), input.d(), input.d(), input.d()],
    }
}

/// `degree count knots... (x y z w)...` with the points as stored.
fn stored_curve(input: &mut Input) -> BSplineCurve {
    let degree = input.i() as usize;
    let count = input.i() as usize;
    let knots = (0..count + degree + 1).map(|_| input.d()).collect();
    let points = (0..count).map(|_| point(input)).collect();
    BSplineCurve::new(degree, knots, points).expect("valid test curve")
}

fn dump_curve(out: &mut Output, degree: usize, knots: &[f64], points: &[ControlPoint]) {
    out.i(degree as i64);
    out.i(points.len() as i64);
    out.i(knots.len() as i64);
    for &knot in knots {
        out.d(knot);
    }
    for point in points {
        for value in point.coordinates {
            out.d(value);
        }
    }
}

fn dump_bspline(out: &mut Output, curve: &BSplineCurve) {
    dump_curve(out, curve.degree(), curve.knots(), curve.control_points());
}

fn pseudo(out: &mut Output, b: &PseudoBox) {
    for basis in b.basis {
        out.v3(basis);
    }
    for interval in b.intervals {
        out.d(interval.min);
        out.d(interval.max);
    }
}

fn values(out: &mut Output, vectors: &[Vector3]) {
    for &v in vectors {
        out.v3(v);
    }
}

fn line(input: &mut Input, out: &mut Output) {
    let mode = input.i();
    let mut line = if mode == 0 {
        let (point, direction) = (input.v3(), input.v3());
        let (a, b, scale) = (input.d(), input.d(), input.d());
        let (dimension, inside_out) = (input.i(), input.i());
        out.i(SUCCESS);
        Line::new(
            point,
            direction,
            Extent1::new(a, b),
            scale,
            dimension == 2,
            inside_out != 0,
        )
    } else {
        let (start, end) = (input.v3(), input.v3());
        let (dimension, has_range) = (input.i(), input.i());
        let (a, b) = (input.d(), input.d());
        let range = (has_range != 0).then(|| Extent1::new(a, b));
        match out.status(Line::segment(start, end, dimension == 2, range)) {
            Some(line) => line,
            None => return,
        }
    };
    let tolerance = input.d();
    let mut t = input.d();
    let derivatives = input.i() as usize;
    let mut step_t = input.d();
    let mut part = Extent1::new(input.d(), input.d());
    let reparameterize = input.i();
    let new_range = Extent1::new(input.d(), input.d());

    let natural = line.natural_interval();
    out.d(natural.min);
    out.d(natural.max);
    out.flag(line.is_bounded());
    out.flag(line.is_degenerate(tolerance, None));
    out.flag(line.is_degenerate(tolerance, Some(part)));
    for pass in 0..2 {
        if pass == 1 {
            if reparameterize == 0 {
                break;
            }
            line.edit_parameterization(new_range);
            out.i(SUCCESS);
            let now = line.natural_interval();
            out.d(now.min);
            out.d(now.max);
            t = now.evaluate(0.3);
            step_t = now.evaluate(0.7);
            part = Extent1::new(now.evaluate(0.2), now.evaluate(0.9));
        }
        dump_curve(out, 1, &line.nurb_knots(), &line.nurb_points());
        out.i(SUCCESS);
        out.d(line.nurbs_to_step(t));
        if let Some(vectors) = out.status(line.evaluate(t, derivatives)) {
            values(out, &vectors);
        }
        if let Some(vectors) = out.status(line.evaluate_step(step_t, derivatives)) {
            values(out, &vectors);
        }
        out.i(SUCCESS);
        out.v3(line.evaluate_point(t));
        let (bbox, pseudo_box) = line.bounding_boxes(part);
        out.i(SUCCESS);
        out.v3(bbox.min);
        out.v3(bbox.max);
        pseudo(out, &pseudo_box);
    }
}

fn hermite(input: &mut Input, out: &mut Output) {
    let planar = input.i() != 0;
    let (p1, d1, p2, d2) = (input.v3(), input.v3(), input.v3(), input.v3());
    let (t, derivatives, nonzero) = (input.d(), input.i() as usize, input.i() != 0);
    let flat = |v: Vector3| Vector2::new(v.x, v.y);
    let curve = if planar {
        Hermite::planar(flat(p1), flat(d1), flat(p2), flat(d2))
    } else {
        Hermite::new(p1, d1, p2, d2)
    };
    let natural = curve.natural_interval();
    out.d(natural.min);
    out.d(natural.max);
    if let Some(vectors) = out.status(curve.evaluate(t, derivatives, nonzero)) {
        values(out, &vectors);
    }
    out.i(SUCCESS);
    out.v3(curve.evaluate_point(t));
    let (bbox, pseudo_box) = curve.bounding_boxes();
    out.i(SUCCESS);
    out.v3(bbox.min);
    out.v3(bbox.max);
    pseudo(out, &pseudo_box);
}

fn canon(input: &mut Input, out: &mut Output) {
    let dimension = input.i() as usize;
    let degree = input.i() as usize;
    let count = input.i() as usize;
    let rational = input.i() != 0;
    let range = (input.i() != 0).then(|| Extent1::new(input.d(), input.d()));
    let unique = input.i() as usize;
    let points: Vec<Vector3> = (0..count).map(|_| input.v3()).collect();
    let weights: Vec<f64> = if rational {
        (0..count).map(|_| input.d()).collect()
    } else {
        Vec::new()
    };
    let (mut knots, mut multiplicities) = (Vec::new(), Vec::new());
    for _ in 0..unique {
        knots.push(input.d());
        multiplicities.push(input.i() as usize);
    }
    let result = BSplineCurve::create_canonical(
        dimension,
        degree,
        &points,
        &multiplicities,
        &knots,
        rational.then_some(weights.as_slice()),
        range,
    );
    if let Some(curve) = out.status(result) {
        dump_bspline(out, &curve);
    }
}

fn trim(op: &str, input: &mut Input, out: &mut Output) {
    let curve = stored_curve(input);
    let requested = Extent1::new(input.d(), input.d());
    let unsure = input.i() != 0;
    if unsure {
        if let Some(interval) = out.status(curve.snapped_trim_interval(requested, ZONE_TOLERANCE)) {
            out.d(interval.min);
            out.d(interval.max);
            out.i(0);
        }
        return;
    }
    match curve.trim(requested, ZONE_TOLERANCE) {
        // `trim-short` marks the cases where SMLib removed short segments,
        // `trim-whole` those where they made up the whole cut curve and it
        // removed its interior knots instead; the port must enter the same
        // branch (the flag) and give the same curve.
        Ok(trimmed) if op == "trim" || op == "trim-short" || op == "trim-whole" => {
            out.i(SUCCESS);
            out.d(trimmed.interval.min);
            out.d(trimmed.interval.max);
            out.i(i64::from(trimmed.short_spans));
            dump_bspline(out, &trimmed.curve);
        }
        Err(TrimRefusal::Failed(fault)) => out.i(status_code(fault)),
        // A disagreement over the short-span branch fails the comparison.
        _ => out.i(-1),
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "line" => line(input, out),
        "hermite" => hermite(input, out),
        "canon" => canon(input, out),
        "reparam" => {
            let curve = stored_curve(input);
            let range = Extent1::new(input.d(), input.d());
            if let Some(edited) = out.status(curve.edit_parameterization(range)) {
                dump_bspline(out, &edited);
            }
        }
        "trim" | "trim-short" | "trim-whole" => trim(op, input, out),
        "clip" => {
            let (lo, hi) = (input.v2(), input.v2());
            let (point, vector, param) = (input.v2(), input.v2(), input.d());
            out.d(Extent2::new(lo, hi).clip_line(point, vector, param));
        }
        _ => panic!("unknown operation {op}"),
    }
}

#[test]
fn output_curves_match_smlib() {
    let fixture = include_str!("fixtures/output-curves/cases.txt");
    let report = common::run_gate(fixture, |_| 4, evaluate);
    report.assert_passed(&[
        "regular",
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

#[test]
fn output_curve_candidates_are_bit_identical_with_one_two_and_four_workers() {
    let fixture = include_str!("fixtures/output-curves/cases.txt");
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

// Refusals the C++ gate cannot compare: SMLib reads past its arrays, builds
// a curve with collapsed knots, or leaves a value unset.

fn quadratic() -> (Vec<Vector3>, Vec<f64>, Vec<usize>) {
    (
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(3.0, 0.0, 1.0),
        ],
        vec![0.0, 1.0],
        vec![3, 3],
    )
}

#[test]
fn canonical_refuses_storage_whose_counts_disagree() {
    let (points, knots, multiplicities) = quadratic();
    let make = |points: &[Vector3], multiplicities: &[usize], weights: Option<&[f64]>| {
        BSplineCurve::create_canonical(3, 2, points, multiplicities, &knots, weights, None)
    };
    assert!(make(&points, &multiplicities, None).is_ok());
    assert_eq!(
        make(&points[..2], &multiplicities, None),
        Err(Fault::InvalidInput)
    );
    assert_eq!(make(&points, &[3, 2], None), Err(Fault::InvalidInput));
    assert_eq!(make(&points, &[3], None), Err(Fault::InvalidInput));
    assert_eq!(
        make(&points, &multiplicities, Some(&[1.0, 2.0])),
        Err(Fault::InvalidInput)
    );
    assert_eq!(make(&[], &multiplicities, None), Err(Fault::InvalidInput));
    let nan = [Vector3::new(f64::NAN, 0.0, 0.0), points[1], points[2]];
    assert_eq!(make(&nan, &multiplicities, None), Err(Fault::InvalidInput));
    assert_eq!(
        BSplineCurve::create_canonical(
            3,
            2,
            &points,
            &multiplicities,
            &knots,
            None,
            Some(Extent1::new(0.0, f64::NAN))
        ),
        Err(Fault::InvalidInput)
    );
}

#[test]
fn canonical_keeps_planar_points_planar_in_two_dimensions() {
    let planar: Vec<Vector3> = (0..3)
        .map(|i| Vector3::new(f64::from(i), 1.0, NO_Z))
        .collect();
    let flat =
        BSplineCurve::create_canonical(2, 2, &planar, &[3, 3], &[0.0, 1.0], None, None).unwrap();
    assert_eq!(flat.dimension(), 2);
    let lifted =
        BSplineCurve::create_canonical(3, 2, &planar, &[3, 3], &[0.0, 1.0], None, None).unwrap();
    assert_eq!(lifted.dimension(), 3);
    assert!(
        lifted
            .control_points()
            .iter()
            .all(|p| p.coordinates[2] == 0.0)
    );
    assert!(
        lifted
            .control_points()
            .iter()
            .all(|p| p.coordinates[3] == NO_WEIGHT)
    );
}

#[test]
fn iges_126_refuses_what_smlib_reads_past() {
    let knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let points = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(2.0, 0.0, 0.0),
    ];
    let ok = |degree, knots: &[f64], weights: &[f64], bounds| {
        iges_126_curve(degree, knots, weights, &points, bounds)
    };
    assert!(ok(2, &knots, &[1.0; 3], [0.0, 1.0]).is_ok());
    assert_eq!(
        ok(2, &knots[..5], &[1.0; 3], [0.0, 1.0]),
        Err(Fault::InvalidInput)
    );
    assert_eq!(
        ok(2, &knots, &[1.0; 2], [0.0, 1.0]),
        Err(Fault::InvalidInput)
    );
    assert_eq!(
        ok(2, &knots, &[1.0; 3], [0.0, f64::INFINITY]),
        Err(Fault::InvalidInput)
    );
    assert_eq!(
        ok(0, &[0.0, 0.0, 1.0, 1.0], &[1.0; 3], [0.0, 1.0]),
        Err(Fault::Degenerate)
    );
    assert_eq!(ok(2, &knots, &[1.0; 3], [0.5, 0.5]), Err(Fault::Degenerate));
    assert_eq!(
        ok(2, &knots, &[1.0; 3], [-0.5, 1.0]),
        Err(Fault::Degenerate)
    );
    assert_eq!(
        ok(2, &knots, &[1.0, 1000.0, 1.0], [0.0, 1.0]),
        Err(Fault::Degenerate)
    );
    let reversed = [0.0, 0.0, 1.0, 0.5, 1.0, 1.0];
    assert_eq!(
        ok(2, &reversed, &[1.0; 3], [0.0, 1.0]),
        Err(Fault::Degenerate)
    );
    assert_eq!(MAX_DEGREE, 17);
}

#[test]
fn reparameterization_refuses_a_range_of_no_length() {
    let curve = BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![
            ControlPoint {
                coordinates: [0.0, 0.0, 0.0, NO_WEIGHT],
            },
            ControlPoint {
                coordinates: [1.0, 0.0, 0.0, NO_WEIGHT],
            },
        ],
    )
    .unwrap();
    assert_eq!(
        curve.edit_parameterization(Extent1::new(2.0, 2.0)),
        Err(Fault::InvalidInput)
    );
    assert_eq!(
        curve.edit_parameterization(Extent1::new(0.0, f64::NAN)),
        Err(Fault::InvalidInput)
    );
    let moved = curve.edit_parameterization(Extent1::new(2.0, 4.0)).unwrap();
    assert_eq!(moved.knots(), [2.0, 2.0, 4.0, 4.0]);
    assert_eq!(
        curve
            .edit_parameterization(Extent1::new(0.0, 1.0 + 1e-13))
            .unwrap(),
        curve
    );
}

#[test]
fn a_line_of_no_length_has_no_valid_curve_but_keeps_its_interval() {
    let degenerate = Line::new(
        Vector3::new(1.0, 2.0, 3.0),
        Vector3::new(0.0, 0.0, 0.0),
        Extent1::new(0.0, 1.0),
        1.0,
        false,
        false,
    );
    assert_eq!(degenerate.make_nurb(), Err(Fault::InvalidInput));
    assert_eq!(degenerate.natural_interval(), Extent1::new(0.0, 0.0));
    assert!(degenerate.is_degenerate(0.0, None));
    // With no parameter range the NURBS parameter is the STEP parameter.
    assert_eq!(degenerate.evaluate_point(0.0), Vector3::new(1.0, 2.0, 3.0));
    assert_eq!(degenerate.evaluate_point(0.3), Vector3::new(1.3, 2.0, 3.0));
}

#[test]
fn a_line_and_its_curve_agree() {
    let line = Line::segment(
        Vector3::new(0.0, 1.0, 2.0),
        Vector3::new(4.0, 3.0, -2.0),
        false,
        Some(Extent1::new(5.0, 9.0)),
    )
    .unwrap();
    let curve = line.make_nurb().unwrap();
    assert_eq!(curve.natural_interval(), line.natural_interval());
    for t in [5.0, 6.0, 7.5, 9.0] {
        let (a, b) = (line.evaluate_point(t), curve.evaluate_point(t).unwrap());
        assert!((a - b).length() < 1e-14, "{a:?} {b:?}");
    }
}

#[test]
fn hermite_refuses_derivatives_it_leaves_unset() {
    let curve = Hermite::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    );
    assert!(curve.evaluate(0.5, 3, true).is_ok());
    assert_eq!(curve.evaluate(0.5, 4, true), Err(Fault::InvalidInput));
    let end = curve.evaluate_point(1.0);
    assert!((end - Vector3::new(1.0, 1.0, 0.0)).length() < 1e-15);
}

#[test]
fn trim_cuts_the_same_curve_and_removes_a_short_span() {
    let curve = BSplineCurve::new(
        3,
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0],
        [(0.0, 0.0), (1.0, 2.0), (3.0, 2.0), (4.0, 0.0), (5.0, 1.0)]
            .map(|(x, y)| ControlPoint {
                coordinates: [x, y, 0.0, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap();
    let cut = curve.trim(Extent1::new(0.3, 1.7), ZONE_TOLERANCE).unwrap();
    assert!(cut.cut);
    assert_eq!(cut.interval, Extent1::new(0.3, 1.7));
    for u in [0.3, 0.8, 1.0, 1.2, 1.7] {
        let (a, b) = (
            curve.evaluate_point(u).unwrap(),
            cut.curve.evaluate_point(u).unwrap(),
        );
        assert!((a - b).length() < 1e-12, "u {u}: {a:?} {b:?}");
    }
    let whole = curve.trim(Extent1::new(-1.0, 3.0), ZONE_TOLERANCE).unwrap();
    assert!(!whole.cut);
    assert_eq!(whole.curve, curve);
    // The cut leaves a sliver span past the knot at 1; SMLib's short-span
    // removal drops it, keeping the parameter range and the end point.
    let sliver = curve
        .trim(Extent1::new(0.0, 1.0000001), ZONE_TOLERANCE)
        .unwrap();
    assert!(sliver.short_spans);
    assert_eq!(sliver.curve.degree(), 3);
    assert_eq!(
        sliver.curve.knots(),
        [
            0.0, 0.0, 0.0, 0.0, 1.0000001, 1.0000001, 1.0000001, 1.0000001
        ]
    );
    let end = curve.evaluate_point(1.0000001).unwrap();
    let kept = sliver.curve.evaluate_point(1.0000001).unwrap();
    assert!((end - kept).length() < 1e-12, "{end:?} {kept:?}");
    assert_eq!(
        curve.trim(Extent1::new(1.0, 1.0), ZONE_TOLERANCE),
        Err(TrimRefusal::Failed(Fault::InvalidInput))
    );
    assert_eq!(
        curve.trim(Extent1::new(-1.0, 1.0), ZONE_TOLERANCE),
        Err(TrimRefusal::Failed(Fault::Degenerate))
    );
}

// Review witnesses for #429.

// Review witness for the #487 short-span route: NLib's basis evaluator holds
// degrees through NL_MAXDEG, but the curve constructor accepts more.
#[test]
fn a_short_span_trim_refuses_a_degree_nlib_cannot_evaluate() {
    let mut knots = vec![0.0; 19];
    knots.push(1.0);
    knots.extend([2.0; 19]);
    let curve = BSplineCurve::new(
        18,
        knots,
        vec![
            ControlPoint {
                coordinates: [0.0, 0.0, 0.0, NO_WEIGHT],
            };
            20
        ],
    )
    .unwrap();
    assert_eq!(
        curve.trim(Extent1::new(0.0, 1.0000001), ZONE_TOLERANCE),
        Err(TrimRefusal::Failed(Fault::InvalidInput))
    );
}

#[test]
fn reparameterization_refuses_unclamped_knots_it_would_reshape() {
    let curve = BSplineCurve::new(
        2,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [0.0, 1.0, 0.0]
            .map(|x| ControlPoint {
                coordinates: [x, 0.0, 0.0, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap();
    assert_eq!(
        curve.edit_parameterization(Extent1::new(0.0, 1.0)),
        Err(Fault::InvalidInput)
    );
    // The same range changes nothing, so it stays accepted.
    assert_eq!(
        curve.edit_parameterization(Extent1::new(0.0, 5.0)).unwrap(),
        curve
    );
    // Trimming assumes clamped knots as well.
    assert_eq!(
        curve.trim(Extent1::new(2.2, 2.8), ZONE_TOLERANCE),
        Err(TrimRefusal::Failed(Fault::InvalidInput))
    );
}

#[test]
fn a_clamped_reparameterization_keeps_the_curve() {
    let curve = BSplineCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0],
        [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)]
            .map(|(x, y)| ControlPoint {
                coordinates: [x, y, 0.0, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap();
    let moved = curve
        .edit_parameterization(Extent1::new(10.0, 14.0))
        .unwrap();
    for i in 0..=8 {
        let fraction = f64::from(i) / 8.0;
        let a = curve.evaluate_point(2.0 * fraction).unwrap();
        let b = moved.evaluate_point(10.0 + 4.0 * fraction).unwrap();
        assert!((a - b).length() < 1e-12, "{a:?} {b:?}");
    }
}

#[test]
fn line_evaluation_refuses_derivative_counts_it_cannot_return() {
    let line = Line::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Extent1::new(0.0, 1.0),
        1.0,
        false,
        false,
    );
    for count in [usize::MAX, usize::MAX - 1, 5, 1 << 40] {
        assert_eq!(line.evaluate(0.5, count), Err(Fault::InvalidInput));
        assert_eq!(line.evaluate_step(0.5, count), Err(Fault::InvalidInput));
    }
    assert_eq!(line.evaluate(0.5, 4).unwrap().len(), 5);
    assert_eq!(line.evaluate_step(0.5, 0).unwrap().len(), 1);
}

#[test]
fn canonical_multiplicities_are_summed_before_they_are_expanded() {
    let points = [Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
    let make = |multiplicities: &[usize]| {
        BSplineCurve::create_canonical(3, 1, &points, multiplicities, &[0.0, 1.0], None, None)
    };
    assert!(make(&[2, 2]).is_ok());
    assert_eq!(make(&[usize::MAX, 2]), Err(Fault::InvalidInput));
    assert_eq!(make(&[usize::MAX, usize::MAX]), Err(Fault::InvalidInput));
    assert_eq!(make(&[1 << 40, 1]), Err(Fault::InvalidInput));
    assert_eq!(make(&[3, 2]), Err(Fault::InvalidInput));
}

#[test]
fn iges_126_refuses_nonfinite_weights_and_unrepresentable_degrees() {
    let points = [Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
    let knots = [0.0, 0.0, 1.0, 1.0];
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for weights in [[bad, 1.0], [1.0, bad]] {
            let result = iges_126_curve(1, &knots, &weights, &points, [0.0, 1.0]);
            assert!(result.is_err(), "{weights:?} accepted");
        }
    }
    assert_eq!(
        iges_126_curve(1, &knots, &[f64::NAN, 1.0], &points, [0.0, 1.0]),
        Err(Fault::InvalidInput)
    );
    for degree in [
        usize::MAX,
        usize::MAX - 1,
        usize::MAX - 2,
        1 << 62,
        MAX_DEGREE + 1,
    ] {
        assert_eq!(
            iges_126_curve(degree, &knots, &[1.0, 1.0], &points, [0.0, 1.0]),
            Err(Fault::Degenerate),
            "degree {degree}"
        );
    }
}
