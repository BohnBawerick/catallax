//! Finite domains, knot-side tangents, periodic seams and ambiguous preimages.
mod common;
use common::*;

fn polynomial(degree: usize, knots: &[&str], controls: &[[&str; 3]]) -> BSplineCurve {
    BSplineCurve::new(
        degree,
        knots.iter().map(|q| r(q)).collect(),
        controls.iter().map(|q| p(*q).to_real()).collect(),
        vec![Real::one(); controls.len()],
        &mut Decider::default(),
    )
    .unwrap()
}

fn contains(actual: [Ball; 3], expected: [&str; 3]) {
    for (actual, expected) in actual.iter().zip(expected.map(r)) {
        assert!(actual.contains(&expected), "{actual:?} excludes {expected:?}");
    }
}

#[test]
fn finite_corner_keeps_both_tangents_and_clips_trace_spans() {
    let curve = polynomial(
        1,
        &["0", "0", "1", "2", "2"],
        &[["0", "0", "0"], ["1", "0", "0"], ["1", "2", "0"]],
    );
    let mut geom = GeomStore::new();
    let curve = geom.add_curve(Curve::BSpline(curve));
    let start = explicit(&mut geom, p(["0", "0", "0"]));
    let corner = explicit(&mut geom, p(["1", "0", "0"]));
    let end = explicit(&mut geom, p(["1", "2", "0"]));
    let outside = explicit(&mut geom, p(["1", "3", "0"]));
    let mut dec = Decider::default();
    coincident(on_curve(&geom, curve, corner, &mut dec));
    assert!(matches!(on_curve(&geom, curve, outside, &mut dec), Outcome::Proven(Apart)));
    assert!(matches!(order_on_curve(&geom, curve, start, end, &mut dec), Outcome::Proven(StrictOrder::Less)));
    assert!(matches!(order_on_curve(&geom, curve, end, start, &mut dec), Outcome::Proven(StrictOrder::Greater)));
    let [arriving, leaving] = tangent_at(&geom, curve, corner, &mut dec).unwrap();
    vector_eq(&arriving, &vi([1, 0, 0]).to_real());
    vector_eq(&leaving, &vi([0, 2, 0]).to_real());
    let traced = exact(trace(&geom, curve, Some(start), Some(end), &mut dec));
    assert_eq!(traced.pieces().len(), 2);
    let half = Ball::from_rational(&r("1/2"), 128);
    contains(traced.pieces()[0].evaluate(&half, &dec).unwrap(), ["1/2", "0", "0"]);
    contains(traced.pieces()[1].evaluate(&half, &dec).unwrap(), ["1", "1", "0"]);
    contains(traced.pieces()[0].derivative(&half, &dec).unwrap(), ["1", "0", "0"]);
    contains(traced.pieces()[1].derivative(&half, &dec).unwrap(), ["0", "2", "0"]);
    let a = explicit(&mut geom, p(["1/2", "0", "0"]));
    let b = explicit(&mut geom, p(["1", "1", "0"]));
    let partial = exact(trace(&geom, curve, Some(a), Some(b), &mut dec));
    contains(partial.pieces()[0].evaluate(&half, &dec).unwrap(), ["3/4", "0", "0"]);
    contains(partial.pieces()[1].evaluate(&half, &dec).unwrap(), ["1", "1/2", "0"]);
}

#[test]
fn a_self_crossing_point_is_incident_but_has_no_arbitrarily_selected_tangent() {
    // C(t)=((t-1/4)(t-3/4), t(t-1/4)(t-3/4), 0).
    let curve = polynomial(
        3,
        &["0", "0", "0", "0", "1", "1", "1", "1"],
        &[
            ["3/16", "0", "0"],
            ["-7/48", "1/16", "0"],
            ["-7/48", "-5/24", "0"],
            ["3/16", "3/16", "0"],
        ],
    );
    let mut dec = Decider::default();
    let pose = Pose::new(
        Rotation::from_quaternion(n(3), n(0), n(1), n(0)).unwrap(),
        vi([17, -23, 31]).to_real(),
        &mut dec,
    ).unwrap();
    for (curve, point) in [
        (curve.clone(), Point3::origin()),
        (curve.transformed(&pose), pose.apply_point(&Point3::origin())),
    ] {
        let mut geom = GeomStore::new();
        let curve = geom.add_curve(Curve::BSpline(curve));
        let point = geom.add_point(PointDef::Explicit(point));
        coincident(on_curve(&geom, curve, point, &mut dec));
        assert!(matches!(tangent_at(&geom, curve, point, &mut dec), Err(Refusal { code: "shapes.spline-ambiguous-parameter", .. })));
    }
}

#[test]
fn closed_seam_has_one_parameter_two_tangents_and_wrapping_arcs() {
    // C(t)=(3t(1-t)^2, 3t^2(1-t), 0); the only repeated point is its seam.
    let raw = polynomial(
        3,
        &["0", "0", "0", "0", "1", "1", "1", "1"],
        &[["0", "0", "0"], ["1", "0", "0"], ["0", "1", "0"], ["0", "0", "0"]],
    );
    let mut dec = Decider::default();
    let seam = raw.unique_parameter_of(&Point3::origin(), &mut dec).unwrap().unwrap();
    real_eq(seam.value(), &Real::zero());
    let mut geom = GeomStore::new();
    let curve = geom.add_curve(Curve::BSpline(raw));
    let seam = explicit(&mut geom, pi([0, 0, 0]));
    let quarter = explicit(&mut geom, p(["27/64", "9/64", "0"]));
    let half = explicit(&mut geom, p(["3/8", "3/8", "0"]));
    let last_quarter = explicit(&mut geom, p(["9/64", "27/64", "0"]));
    let [arriving, leaving] = tangent_at(&geom, curve, seam, &mut dec).unwrap();
    vector_eq(&arriving, &vi([0, -3, 0]).to_real());
    vector_eq(&leaving, &vi([3, 0, 0]).to_real());
    assert!(matches!(cyclic_order(&geom, curve, [seam, quarter, half], &mut dec), Outcome::Proven(Cyclic::Forward)));
    assert!(matches!(cyclic_order(&geom, curve, [seam, half, quarter], &mut dec), Outcome::Proven(Cyclic::Backward)));
    let whole = exact(trace(&geom, curve, None, None, &mut dec));
    contains(whole.pieces()[0].evaluate(&Ball::from_rational(&r("1/4"), 128), &dec).unwrap(), ["27/64", "9/64", "0"]);
    let wrapped = exact(trace(&geom, curve, Some(last_quarter), Some(quarter), &mut dec));
    assert_eq!(wrapped.pieces().len(), 2);
    contains(wrapped.pieces()[0].evaluate(&Ball::from_rational(&n(1), 128), &dec).unwrap(), ["0", "0", "0"]);
    contains(wrapped.pieces()[1].evaluate(&Ball::from_rational(&n(0), 128), &dec).unwrap(), ["0", "0", "0"]);
}

#[test]
fn distinct_parameterisations_are_not_declared_different_point_sets() {
    let forward = polynomial(1, &["0", "0", "1", "1"], &[["0", "0", "0"], ["1", "0", "0"]]);
    let reverse = polynomial(1, &["0", "0", "1", "1"], &[["1", "0", "0"], ["0", "0", "0"]]);
    let mut geom = GeomStore::new();
    let a = geom.add_curve(Curve::BSpline(forward));
    let b = geom.add_curve(Curve::BSpline(reverse));
    assert!(!matches!(same_curve(&geom, a, b, &mut Decider::default()), Outcome::Proven(Apart)));
}
