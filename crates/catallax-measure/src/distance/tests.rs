use super::*;
use catallax_shapes::{Pose, RPoint3, Surface};
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64, z: i64) -> Point3 {
    RPoint3::new(q(x), q(y), q(z)).into()
}
fn controls() -> Options {
    Options {
        width: Rational::from_ratio(1, 1000).unwrap(),
        max_cells: 4096,
    }
}
#[test]
fn cylinder_arc_trims_exclude_the_nearest_untrimmed_support_point() {
    let (mut model, body) = crate::fixtures::half_cylinder(&Pose::identity());
    let wall = model
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, wall, &mut dec).unwrap();
    let answer = prepared
        .point_to_face(&mut model, &p(0, -2, 1), &controls(), &mut dec)
        .unwrap();
    let want = Ball::from_rational(&q(8), 256).sqrt().unwrap();
    assert!(answer.lower <= want.lo(), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() >= &want.hi(), "{answer:?}");
    assert!(answer.lower > q(2), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() < &q(3), "{answer:?}");
    let [_, witness] = answer.witnesses.unwrap();
    assert!(catallax_boolean::face_contains(&mut model, wall, &witness, &mut dec).unwrap());
}
#[test]
fn cap_pair_minimum_distance_and_cylinder_height_boundary_are_certified() {
    let mut model = Model::new();
    let body = catallax_topology::cylinder(&mut model, q(1), q(2), &Pose::identity()).unwrap();
    let faces = model.faces_of_body(body);
    let wall = faces
        .iter()
        .copied()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    let caps: Vec<_> = faces
        .into_iter()
        .filter(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_)
            )
        })
        .collect();
    let mut dec = Decider::default();
    let wall = PreparedFace::new(&model, wall, &mut dec).unwrap();
    let answer = wall
        .point_to_face(&mut model, &p(1, 0, 3), &controls(), &mut dec)
        .unwrap();
    assert!(
        answer.lower <= q(1) && answer.lower >= Rational::parse("0.999").unwrap(),
        "{answer:?}"
    );
    assert!(
        answer.upper.as_ref().unwrap() >= &q(1)
            && answer.upper.as_ref().unwrap() <= &Rational::parse("1.001").unwrap(),
        "{answer:?}"
    );
    let a = PreparedFace::new(&model, caps[0], &mut dec).unwrap();
    let b = PreparedFace::new(&model, caps[1], &mut dec).unwrap();
    let answer = a
        .face_to_face(&b, &mut model, &controls(), &mut dec)
        .unwrap();
    assert_eq!(answer.lower, q(2));
    assert_eq!(answer.upper, Some(q(2)));
}
#[test]
fn spherical_cap_does_not_use_the_removed_south_pole() {
    let (mut model, body) = crate::fixtures::cap(5, 3);
    let face = model
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Sphere(_)
            )
        })
        .unwrap();
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    let answer = prepared
        .point_to_face(&mut model, &p(0, 0, -5), &controls(), &mut dec)
        .unwrap();
    // The boundary circle has radius four and z=3, so the minimum is sqrt(80).
    let want = Ball::from_rational(&q(80), 256).sqrt().unwrap();
    assert!(answer.lower <= want.lo(), "{answer:?}");
    assert!(answer.upper.as_ref().unwrap() >= &want.hi(), "{answer:?}");
    assert!(answer.lower > q(8), "{answer:?}");
}

#[test]
fn elliptic_cylinder_preserves_euclidean_distance_and_oblique_circle_trims() {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, EllipticCylinder, GeomKey, RVector3};
    use catallax_topology::{Coedge, Edge, Face, Loop, Sense};

    let mut model = Model::new();
    let surface = model
        .geom_mut()
        .add_surface(Surface::EllipticCylinder(EllipticCylinder {
            origin: p(0, 0, 0),
            axis: RVector3::new(q(3), q(0), q(4)).into(),
            major_dir: RVector3::new(q(0), q(1), q(0)).into(),
            major_sq: q(1).into(),
            minor_sq: Rational::from_ratio(16, 25).unwrap().into(),
        }));
    // The unit circles at (0,0,0) and (3,0,4), both parallel to xy,
    // lie on this carrier. They bound an oblique circular sweep.
    let mut loops = Vec::new();
    for (center, sense) in [(p(0, 0, 0), Sense::Same), (p(3, 0, 4), Sense::Reversed)] {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center,
            normal: RVector3::new(q(0), q(0), q(1)).into(),
            radius_sq: q(1).into(),
            ref_dir: RVector3::new(q(1), q(0), q(0)).into(),
        }));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Curve(curve),
            b: GeomKey::Surface(surface),
            at: None,
            origin: FactOrigin::Construction("measure.distance-elliptic-cylinder"),
        });
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let coedge = model.add_coedge(Coedge { edge, sense });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops,
    });
    let mut dec = Decider::default();
    let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
    for (xyz, expected) in [
        (["29/10", "0", "17/10"], "1/2"),
        (["23/5", "0", "24/5"], "1"),
    ] {
        let point = RPoint3::new(
            Rational::parse(xyz[0]).unwrap(),
            Rational::parse(xyz[1]).unwrap(),
            Rational::parse(xyz[2]).unwrap(),
        )
        .to_real();
        let answer = prepared
            .point_to_face(&mut model, &point, &controls(), &mut dec)
            .unwrap();
        let expected = Rational::parse(expected).unwrap();
        let upper = answer.upper.as_ref().unwrap();
        assert!(answer.lower <= expected && upper >= &expected, "{answer:?}");
        assert!(
            answer.converged && upper - &answer.lower <= controls().width,
            "{answer:?}"
        );
        let [_, witness] = answer.witnesses.unwrap();
        assert!(catallax_boolean::face_contains(&mut model, face, &witness, &mut dec).unwrap());
    }
}
