#![allow(missing_docs)]
use catallax_numbers::{BudgetPolicy, Decider, FactKind, Outcome, Rational, ZeroProof};
use catallax_shapes::*;

fn point(x: i64, y: i64, z: i64) -> Point3 {
    RPoint3::new(x.into(), y.into(), z.into()).to_real()
}
fn vector(x: i64, y: i64, z: i64) -> Vector3 {
    RVector3::new(x.into(), y.into(), z.into()).to_real()
}
fn plane(z: i64, normal: i64) -> Surface {
    Surface::Plane(Plane {
        origin: point(0, 0, z),
        normal: vector(0, 0, normal),
        ref_dir: vector(1, 0, 0),
    })
}
fn no_work() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 0,
        ..BudgetPolicy::default()
    })
}

#[test]
fn flush_construction_is_decided_without_arithmetic_in_either_insertion_order() {
    for point_first in [false, true] {
        let mut g = GeomStore::new();
        let (p, a) = if point_first {
            let p = g.add_point(PointDef::Explicit(point(3, 4, 2)));
            (p, g.add_surface(plane(2, 1)))
        } else {
            let a = g.add_surface(plane(2, 1));
            (g.add_point(PointDef::Explicit(point(3, 4, 2))), a)
        };
        let b = g.add_surface(plane(2, -1));
        let mut d = no_work();
        for result in [
            g.side(a, p, &mut d).map(|_| Apart),
            g.same_surface(a, b, &mut d),
        ] {
            assert!(
                matches!(result, Outcome::ProvenCoincident(c) if c.proof == ZeroProof::ConstructionFact)
            );
        }
        assert_eq!(d.stats().work, 0);
        assert_eq!(d.stats().coincident_fact, 2);
        let near = g.add_point(PointDef::Explicit(Point3::new(
            3.into(),
            4.into(),
            Rational::parse("2.0000000000000000000000000001")
                .unwrap()
                .into(),
        )));
        assert!(!g.facts().holds(
            FactKind::Incident,
            GeomKey::Point(near),
            GeomKey::Surface(a)
        ));
        assert!(matches!(
            g.side(a, near, &mut Decider::default()),
            Outcome::Proven(Side::Front)
        ));
    }
}

#[test]
fn shared_axes_and_centres_do_not_imply_equal_radii() {
    let mut g = GeomStore::new();
    let cylinder = |radius: i64| {
        Surface::Cylinder(Cylinder {
            origin: point(0, 0, 0),
            axis: vector(0, 0, 1),
            radius: radius.into(),
            ref_dir: vector(1, 0, 0),
        })
    };
    let sphere = |radius: i64| {
        Surface::Sphere(Sphere {
            center: point(0, 0, 0),
            axis: vector(0, 0, 1),
            radius: radius.into(),
            ref_dir: vector(1, 0, 0),
        })
    };
    let a = g.add_surface(cylinder(2));
    let b = g.add_surface(cylinder(3));
    let c = g.add_surface(sphere(2));
    let d = g.add_surface(sphere(3));
    assert!(
        g.facts()
            .holds(FactKind::Coaxial, GeomKey::Surface(a), GeomKey::Surface(b))
    );
    assert!(g.facts().holds(
        FactKind::Concentric,
        GeomKey::Surface(c),
        GeomKey::Surface(d)
    ));
    let mut dec = Decider::default();
    assert_eq!(g.same_surface(a, b, &mut dec), Outcome::Proven(Apart));
    assert_eq!(g.same_surface(c, d, &mut dec), Outcome::Proven(Apart));
}
