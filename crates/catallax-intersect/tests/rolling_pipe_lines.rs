//! Lines against the rolling pipe, through the same call point classification uses.

use catallax_intersect::{CurveSurface, intersect_curve_surface};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_shapes::{
    Curve, Cylinder, GeomStore, Line, Plane, Point3, RollingPipe, SplineSupport, Surface, Vector3,
};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn real(value: &str) -> Real {
    Real::from(Rational::parse(value).unwrap())
}
fn vector(v: [&str; 3]) -> Vector3 {
    let [x, y, z] = v.map(real);
    Vector3::new(x, y, z)
}
fn point(v: [&str; 3]) -> Point3 {
    let [x, y, z] = v.map(real);
    Point3::new(x, y, z)
}
// A frame direction square to `v`, for whichever of its first two parts is nonzero.
fn across(v: &Vector3) -> Vector3 {
    Vector3::new(v.y.clone(), -v.x.clone(), Real::zero())
}
fn pipe(origin: [&str; 3], normal: [&str; 3], base: [&str; 3], axis: [&str; 3]) -> RollingPipe {
    let normal = vector(normal);
    let axis = vector(axis);
    let plane = Plane {
        origin: point(origin),
        ref_dir: across(&normal),
        normal,
    };
    let cylinder = Cylinder {
        origin: point(base),
        ref_dir: across(&axis),
        axis,
        radius_sq: real("16"),
    };
    RollingPipe::new(
        &plane,
        &cylinder,
        [true, true],
        real("1"),
        true,
        &mut decider(),
    )
    .unwrap()
}
fn meet(pipe: RollingPipe, start: Point3, direction: Vector3) -> Outcome<usize> {
    let mut geom = GeomStore::new();
    let mut dec = decider();
    let line = Line::new(start, direction, &mut dec).unwrap();
    let curve = geom.add_curve(Curve::Line(line));
    let surface = geom.add_surface(Surface::BSpline(SplineSupport::RollingPipe(pipe)));
    let answer = intersect_curve_surface(&mut geom, curve, surface, &mut dec);
    match answer {
        Ok(answer) => match answer.value {
            CurveSurface::Meets(hits) => Outcome::Proven(hits.len()),
            other => panic!("the line should cross the band: {other:?}"),
        },
        Err(refusal) => Outcome::Refused(refusal),
    }
}

#[test]
fn a_rational_law_answers_a_line_through_its_band() {
    let pipe = pipe(
        ["0", "0", "6"],
        ["3/5", "0", "4/5"],
        ["0", "0", "0"],
        ["0", "0", "1"],
    );
    let start = pipe.point(&Rational::from_i64(1), &half());
    let direction = start.minus(&pipe.centre(&Rational::from_i64(1)));
    assert!(matches!(
        meet(pipe, start, direction),
        Outcome::Proven(count) if count >= 1
    ));
}

#[test]
fn an_irrational_law_refuses_a_line_it_cannot_rule_out() {
    let pipe = pipe(
        ["1", "2", "3"],
        ["1", "2", "3"],
        ["1", "1", "1"],
        ["0", "1", "1"],
    );
    let start = pipe.point(&Rational::from_i64(1), &half());
    let direction = vector(["1/3", "1/5", "1/7"]);
    let Outcome::Refused(refusal) = meet(pipe, start, direction) else {
        panic!("an irrational law has no line solver yet");
    };
    assert_eq!(refusal.code, "intersect.parametric-surface");
}

fn half() -> Rational {
    Rational::parse("1/2").unwrap()
}
