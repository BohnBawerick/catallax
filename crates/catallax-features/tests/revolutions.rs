//! Full and partial revolution mass properties, names and impossible axes.

use catallax_features::{Profile, revolve};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{
    Axis, ElementId, Loop, Point2, Point3, Pose, Segment, StepId, Turn, Vector3,
};
use catallax_topology::{BodyId, Model};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64) -> Point2 {
    Point2([x.into(), y.into()])
}
fn polygon(points: &[(i64, i64)]) -> Loop {
    Loop(
        points
            .iter()
            .enumerate()
            .map(|(i, (x, y))| Segment::Line {
                id: ElementId::new(format!("s{i}")).unwrap(),
                to: p(*x, *y),
            })
            .collect(),
    )
}
fn axis() -> Axis {
    Axis {
        point: Point3([0.into(), 0.into(), 0.into()]),
        direction: Vector3([0.into(), 1.into(), 0.into()]),
    }
}
fn build(loops: &[Loop], angle: i64) -> (Model, BodyId) {
    let mut dec = Decider::default();
    let profile = Profile::new(&Pose::identity(), loops, &mut dec).unwrap();
    let mut model = Model::new();
    let body = revolve(
        &mut model,
        &profile,
        &StepId::new("turn").unwrap(),
        &axis(),
        &r(angle),
        &mut dec,
    )
    .unwrap();
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    let origins: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .flat_map(|f| model.face_origins(f))
        .collect();
    assert!(origins.iter().all(|o| o.step.as_str() == "turn"));
    assert_eq!(
        origins
            .iter()
            .filter(|o| matches!(
                o.role,
                catallax_recipe::Role::Start | catallax_recipe::Role::End
            ))
            .count(),
        if angle == 360 { 0 } else { 2 }
    );
    (model, body)
}
fn mass(model: &Model, body: BodyId) -> catallax_measure::MassProperties {
    let answer = catallax_measure::mass_properties(
        model,
        body,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    answer.value
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&r(n), 192)
}
fn contains(e: &catallax_measure::Enclosure, value: Ball) {
    assert!(
        e.lo <= value.lo() && e.hi >= value.hi(),
        "{e:?} versus {value:?}"
    );
    assert!(e.width() <= Rational::parse("1e-9").unwrap());
}

#[test]
fn full_and_partial_annular_revolutions() {
    let ring = polygon(&[(3, 0), (3, 2), (1, 2), (1, 0)]);
    for angle in [90, 180, 270, 360] {
        let (model, body) = build(std::slice::from_ref(&ring), angle);
        let mass = mass(&model, body);
        let fraction = Ball::from_rational(&Rational::from_ratio(angle, 360).unwrap(), 192);
        contains(&mass.volume, b(16) * Ball::pi(192) * fraction.clone());
        contains(
            &mass.area,
            b(32) * Ball::pi(192) * fraction + b(if angle == 360 { 0 } else { 8 }),
        );
        let centroid = mass.centroid.unwrap();
        assert!(centroid[1].contains(&r(1)));
        if angle == 360 {
            assert!(centroid[0].contains(&r(0)));
            assert!(centroid[2].contains(&r(0)));
        }
    }
}

#[test]
fn triangular_profile_revolves_to_a_cone_with_an_apex() {
    let (model, body) = build(&[polygon(&[(3, 0), (0, 4), (0, 0)])], 360);
    let mass = mass(&model, body);
    contains(&mass.volume, b(12) * Ball::pi(192));
    contains(&mass.area, b(24) * Ball::pi(192));
    for (value, expected) in mass.centroid.unwrap().iter().zip([0, 1, 0]) {
        assert!(value.contains(&r(expected)));
    }
}

#[test]
fn circle_revolves_to_a_torus_without_seams() {
    let circle = Loop(vec![Segment::Circle {
        id: ElementId::new("tube").unwrap(),
        centre: p(4, 0),
        radius: 1.into(),
    }]);
    let (model, body) = build(&[circle], 360);
    assert_eq!(model.counts(body).faces, 1);
    assert_eq!(model.counts(body).edges, 0);
    let mass = mass(&model, body);
    let pi = Ball::pi(192);
    contains(&mass.volume, b(8) * pi.square());
    contains(&mass.area, b(16) * pi.square());
}

#[test]
fn semicircle_revolves_to_a_sphere_without_pole_edges() {
    let profile = Loop(vec![
        Segment::Arc {
            id: ElementId::new("round").unwrap(),
            centre: p(0, 0),
            to: p(0, 2),
            turn: Turn::Ccw,
        },
        Segment::Line {
            id: ElementId::new("axis").unwrap(),
            to: p(0, -2),
        },
    ]);
    let (model, body) = build(&[profile], 360);
    let mass = mass(&model, body);
    contains(
        &mass.volume,
        Ball::from_rational(&Rational::from_ratio(32, 3).unwrap(), 192) * Ball::pi(192),
    );
    contains(&mass.area, b(16) * Ball::pi(192));
}

#[test]
fn impossible_revolutions_report_the_recipe_reason() {
    let profile = Profile::new(
        &Pose::identity(),
        &[polygon(&[(2, 0), (2, 2), (-2, 2), (-2, 0)])],
        &mut Decider::default(),
    )
    .unwrap();
    let step = StepId::new("turn").unwrap();
    let mut model = Model::new();
    for angle in [0, -1, 361] {
        assert_eq!(
            revolve(
                &mut model,
                &profile,
                &step,
                &axis(),
                &r(angle),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "revolve.angle"
        );
    }
    assert_eq!(
        revolve(
            &mut model,
            &profile,
            &step,
            &axis(),
            &r(360),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "revolve.axis-interior"
    );
    let mut invalid = axis();
    invalid.direction = Vector3([0.into(), 0.into(), 0.into()]);
    assert_eq!(
        revolve(
            &mut model,
            &profile,
            &step,
            &invalid,
            &r(360),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "revolve.zero-axis"
    );
    invalid.direction = Vector3([0.into(), 0.into(), 1.into()]);
    assert_eq!(
        revolve(
            &mut model,
            &profile,
            &step,
            &invalid,
            &r(360),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "revolve.axis-plane"
    );
}
