use super::*;
use catallax_boolean::{BoolOp, boolean};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Rational, with_operation_budget};
use catallax_shapes::{Pose, Rotation};
use catallax_topology::{block, cone, cylinder, sphere, torus};

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn contours(model: &Model, body: BodyId, d: RVector3) -> Contours {
    contour_generators(model, body, &d, &mut Decider::default())
        .unwrap()
        .value
}

fn is_zero(value: &Real) -> bool {
    sign(value, &mut Decider::default()).unwrap() == 0
}

fn curved_faces(model: &Model, body: BodyId) -> Vec<FaceId> {
    model
        .faces_of_body(body)
        .into_iter()
        .filter(|&face| {
            !matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Plane(_)
            )
        })
        .collect()
}

/// Ends of a line arc, as exact points.
fn ends(c: &Contours, arc: &ContourArc) -> [Point3; 2] {
    let CurveForm::Line { .. } =
        curve_form(c.model.geom(), arc.curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected a line generator");
    };
    [arc.start.unwrap(), arc.end.unwrap()].map(|p| c.model.geom().position(p).clone())
}

/// The surface normal at an exact point of a quadric face, unnormalised.
fn normal_at(surface: &Surface, p: &Point3) -> Vector3 {
    match surface {
        Surface::Cylinder(s) => {
            let r = p.minus(&s.origin);
            r.scaled(&s.axis.norm_squared())
                .minus(&s.axis.scaled(&r.dot(&s.axis)))
        }
        Surface::Cone(s) => {
            let r = p.minus(&s.apex);
            let steep = Real::one() + s.slope.square();
            r.scaled(&s.axis.norm_squared())
                .minus(&s.axis.scaled(&(&steep * &r.dot(&s.axis))))
        }
        Surface::Sphere(s) => p.minus(&s.center),
        other => panic!("no test normal for {:?}", other.kind()),
    }
}

/// Every end of every arc lies on its face's carrier with a normal
/// perpendicular to the view direction.
fn assert_generators(model: &Model, c: &Contours, d: &RVector3) {
    for arc in &c.arcs {
        let surface = model.geom().surface(model.face(arc.face).surface);
        for end in ends(c, arc) {
            assert!(is_zero(&normal_at(surface, &end).dot(&d.to_real())));
        }
    }
}

#[test]
fn cylinder_generators_are_two_axis_lines_and_none_along_the_axis() {
    let mut model = Model::new();
    let body = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let side = curved_faces(&model, body);
    assert_eq!(side.len(), 1);
    let front = contours(&model, body, vector(0, -1, 0));
    assert_eq!(front.arcs.len(), 2);
    let mut xs = Vec::new();
    for arc in &front.arcs {
        assert_eq!(vec![arc.face], side);
        assert_eq!(arc.on_faces, side);
        let [a, b] = ends(&front, arc);
        assert!(is_zero(&a.y) && is_zero(&b.y));
        assert!(is_zero(&(&a.x - &b.x)));
        let mut z = [a.z.as_rational().unwrap(), b.z.as_rational().unwrap()];
        z.sort();
        assert_eq!(z, [Rational::zero(), Rational::from(5)]);
        xs.push(a.x.as_rational().unwrap());
    }
    xs.sort();
    assert_eq!(xs, [Rational::from(-2), Rational::from(2)]);
    assert!(contours(&model, body, vector(0, 0, 1)).arcs.is_empty());
    assert!(contours(&model, body, vector(0, 0, -3)).arcs.is_empty());
    let d = vector(1, 2, 3);
    let oblique = contours(&model, body, d.clone());
    assert_eq!(oblique.arcs.len(), 2);
    assert_generators(&model, &oblique, &d);
}

#[test]
fn cone_generators_meet_the_apex_plane_and_vanish_inside_the_cone_angle() {
    let mut model = Model::new();
    // Radius 3 at z = 0 and 1 at z = 4: slope 1/2, apex at z = 6.
    let body = cone(&mut model, 3.into(), 1.into(), 4.into(), &Pose::identity()).unwrap();
    let front = contours(&model, body, vector(0, -1, 0));
    assert_eq!(front.arcs.len(), 2);
    for arc in &front.arcs {
        let [a, b] = ends(&front, arc);
        let (low, high) = if a.z.as_rational().unwrap().is_zero() {
            (a, b)
        } else {
            (b, a)
        };
        assert!(is_zero(&low.y) && is_zero(&high.y) && is_zero(&low.z));
        assert_eq!(high.z.as_rational().unwrap(), Rational::from(4));
        assert_eq!(low.x.as_rational().unwrap().abs(), Rational::from(3));
        assert_eq!(high.x.as_rational().unwrap().abs(), Rational::one());
    }
    // Along the axis, and inside the half angle (|d_perp| = 1 < 10 / 2), the
    // polar plane meets the nappe only at the apex, which is off the frustum.
    assert!(contours(&model, body, vector(0, 0, 1)).arcs.is_empty());
    assert!(contours(&model, body, vector(0, 1, 10)).arcs.is_empty());
    let d = vector(0, 1, 1);
    let oblique = contours(&model, body, d.clone());
    assert_eq!(oblique.arcs.len(), 2);
    assert_generators(&model, &oblique, &d);
}

#[test]
fn sphere_generator_is_the_whole_great_circle_facing_the_viewer() {
    let mut model = Model::new();
    let body = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let d = vector(1, 2, 2);
    let answer = contours(&model, body, d.clone());
    assert_eq!(answer.arcs.len(), 1);
    let arc = &answer.arcs[0];
    assert_eq!((arc.start, arc.end), (None, None));
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = curve_form(answer.model.geom(), arc.curve, &mut Decider::default()).unwrap()
    else {
        panic!("expected a circle");
    };
    assert!(is_zero(&center.minus(&Point3::origin()).norm_squared()));
    assert!(is_zero(&normal.cross(&d.to_real()).norm_squared()));
    assert!(is_zero(&(&radius_sq - &Real::from(25))));
}

#[test]
fn a_generator_on_a_coplanar_face_runs_along_its_edges() {
    // A half cylinder: its flat face lies in the front view's polar plane.
    let mut model = Model::new();
    let round = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let cutter = block(
        &mut model,
        [10.into(), 10.into(), 10.into()],
        &Pose::from_rational(Rotation::identity(), vector(-5, -10, -1)),
    )
    .unwrap();
    let half = boolean(
        &mut model,
        round,
        cutter,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let side = curved_faces(&model, half);
    assert_eq!(side.len(), 1);
    let front = contours(&model, half, vector(0, 1, 0));
    assert_eq!(front.arcs.len(), 2);
    for arc in &front.arcs {
        assert_eq!(vec![arc.face], side);
        assert_eq!(arc.on_faces.len(), 2);
        let [a, b] = ends(&front, arc);
        assert!(is_zero(&a.y) && is_zero(&b.y));
        assert_eq!(a.x.as_rational().unwrap().abs(), Rational::from(2));
    }
}

#[test]
fn planar_bodies_have_no_generators_and_spend_no_work() {
    let mut model = Model::new();
    let body = block(
        &mut model,
        [1.into(), 2.into(), 3.into()],
        &Pose::identity(),
    )
    .unwrap();
    let (answer, work) = with_operation_budget(1_000_000, || {
        contour_generators(&model, body, &vector(1, 1, 1), &mut Decider::default())
    });
    assert!(answer.unwrap().value.arcs.is_empty());
    assert!(!work.exhausted);
}

#[test]
fn tori_zero_directions_and_exhausted_allowances_refuse_with_counts() {
    let mut model = Model::new();
    let ring = torus(&mut model, 5.into(), 1.into(), &Pose::identity()).unwrap();
    let refusal =
        contour_generators(&model, ring, &vector(0, 1, 1), &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.refusal.code, "app.contour-surface");
    assert_eq!(refusal.operands[0].faces, 1);
    assert_eq!(refusal.operands[0].curved_faces, 1);
    let round = cylinder(&mut model, 2.into(), 5.into(), &Pose::identity()).unwrap();
    let zero =
        contour_generators(&model, round, &vector(0, 0, 0), &mut Decider::default()).unwrap_err();
    assert_eq!(zero.refusal.code, "app.contour-direction");
    let (starved, work) = with_operation_budget(10, || {
        contour_generators(&model, round, &vector(0, 1, 0), &mut Decider::default())
    });
    assert!(starved.is_err());
    assert!(work.exhausted);
}
