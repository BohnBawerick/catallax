use super::*;
use crate::geometry::canonical_point;
use catallax_intersect::{
    Component, CurveSurface, SurfacePair, intersect_curve_surface, intersect_surfaces,
};
use catallax_numbers::{Fact, FactOrigin, Real, StrictSign, with_operation_budget};
use catallax_shapes::{Cylinder, IntersectionCurve, Meet, Vector3};
use catallax_topology::Model;
use std::collections::BTreeMap;

fn r(value: i64) -> Real {
    Real::from_i64(value)
}
fn p(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(r(x), r(y), r(z))
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(r(x), r(y), r(z))
}
fn cylinder(origin: Point3, axis: Vector3, radius: i64) -> Surface {
    Surface::Cylinder(Cylinder {
        origin,
        axis,
        radius_sq: r(radius * radius),
        ref_dir: v(1, 0, 0),
    })
}
fn roots(
    geom: &mut GeomStore,
    pair: [SurfaceId; 2],
    target: SurfaceId,
    dec: &mut Decider,
) -> Vec<PointId> {
    let SurfacePair::Meets(parts) = intersect_surfaces(geom, pair[0], pair[1], dec)
        .unwrap()
        .value
    else {
        panic!("the cylinder supports intersect");
    };
    let mut points = Vec::new();
    for part in parts {
        if let Component::Curve { curve, .. } = part
            && let CurveSurface::Meets(hits) = intersect_curve_surface(geom, curve, target, dec)
                .unwrap()
                .value
        {
            points.extend(hits.into_iter().map(|hit| hit.point));
        }
    }
    points
}
fn quadrant(geom: &GeomStore, point: PointId, dec: &mut Decider) -> [i8; 2] {
    let point = geom.position(point);
    [&point.y - &r(25), &point.z - &r(18)].map(|value| match dec.sign(&value) {
        Outcome::Proven(StrictSign::Negative) => -1,
        Outcome::Proven(StrictSign::Positive) => 1,
        other => panic!("regular triple point quadrant: {other:?}"),
    })
}
type QuadrantPoints = BTreeMap<[i8; 2], PointId>;

fn fixture() -> (Model, [SurfaceId; 3], QuadrantPoints, QuadrantPoints) {
    let mut model = Model::new();
    let geom = model.geom_mut();
    let boss = geom.add_surface(cylinder(p(30, 25, 5), v(0, 0, 1), 8));
    let old = geom.add_surface(cylinder(p(23, 6, 18), v(0, -1, 0), 3));
    let new = geom.add_surface(cylinder(p(25, 6, 18), v(0, -1, 0), 3));
    let mut dec = Decider::default();
    let before = roots(geom, [boss, new], old, &mut dec);
    let after = roots(geom, [old, new], boss, &mut dec);
    let before = before
        .into_iter()
        .map(|point| (quadrant(geom, point, &mut dec), point))
        .collect();
    let after = after
        .into_iter()
        .map(|point| (quadrant(geom, point, &mut dec), point))
        .collect();
    (model, [boss, old, new], before, after)
}

#[test]
fn equivalent_rur_routes_keep_the_existing_four_witnesses_within_the_request() {
    let (mut model, _, before, after) = fixture();
    assert_eq!(
        before.keys().copied().collect::<Vec<_>>(),
        vec![[-1, -1], [-1, 1], [1, -1], [1, 1]]
    );
    assert_eq!(
        before.keys().collect::<Vec<_>>(),
        after.keys().collect::<Vec<_>>()
    );
    let expected: Vec<_> = before.values().copied().collect();
    let mut pool = expected.clone();
    let mut dec = Decider::default();
    let (answers, work) = with_operation_budget(90_000, || {
        after
            .iter()
            .map(|(quadrant, &point)| {
                let canonical = canonical_point(
                    &mut model,
                    &mut pool,
                    point,
                    &mut dec,
                    "boolean.imprint-ambiguous",
                )?;
                Ok::<_, Refusal>((*quadrant, point, canonical))
            })
            .collect::<Result<Vec<_>, _>>()
    });
    let answers = answers.unwrap();
    assert!(!work.exhausted);
    assert_eq!(pool, expected);
    for (quadrant, incoming, canonical) in answers {
        assert_eq!(canonical, before[&quadrant]);
        assert!(model.geom().facts().holds(
            FactKind::Coincident,
            GeomKey::Point(canonical),
            GeomKey::Point(incoming)
        ));
    }
}

#[test]
fn two_subnanometre_roots_in_the_same_certificate_hull_stay_distinct() {
    let mut model = Model::new();
    let geom = model.geom_mut();
    let gap = Rational::one().mul_pow2(-160);
    let distance = Real::from_rational(&Rational::from_i64(6) - &gap);
    let x = &distance / &r(2);
    let z = (&r(9) - &x.square()).sqrt();
    let y = (&r(4) - &(&x - &r(3)).square()).sqrt();
    let boss = geom.add_surface(cylinder(p(3, 0, 0), v(0, 0, 1), 2));
    let old = geom.add_surface(cylinder(p(0, 0, 18), v(0, -1, 0), 3));
    let new = geom.add_surface(cylinder(Point3::new(distance, r(0), r(18)), v(0, -1, 0), 3));
    let low = Point3::new(x.clone(), y.clone(), &r(18) - &z);
    let high = Point3::new(x, y, &r(18) + &z);
    let witness = geom.add_point(PointDef::Explicit(low.clone()));
    let curve = geom.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [boss, new],
        witness,
    }));
    let low = geom.add_point(PointDef::Meet(Meet {
        of: MeetOf::CurveSurface(curve, old),
        position: Position::Closed(low),
    }));
    let high = geom.add_point(PointDef::Meet(Meet {
        of: MeetOf::CurveSurface(curve, old),
        position: Position::Closed(high),
    }));
    let mut pool = Vec::new();
    pool.push(low);
    let mut dec = Decider::default();
    let (answer, work) = with_operation_budget(1_000_000, || {
        canonical_point(
            &mut model,
            &mut pool,
            high,
            &mut dec,
            "boolean.imprint-ambiguous",
        )
    });
    assert_eq!(answer.unwrap(), high);
    assert!(!work.exhausted);
    assert_eq!(pool, vec![low, high]);
    assert!(!model.geom().facts().holds(
        FactKind::Coincident,
        GeomKey::Point(low),
        GeomKey::Point(high)
    ));
}

#[test]
fn two_incidence_facts_do_not_discharge_the_third_cylinder() {
    let (mut model, [boss, old, new], before, _) = fixture();
    let incoming = before[&[-1, -1]];
    let near = Point3::new(
        r(24),
        &(&r(25) - &r(28).sqrt()) + &Real::from_rational(Rational::one().mul_pow2(-80)),
        &r(18) - &r(8).sqrt(),
    );
    let near = model.geom_mut().add_point(PointDef::Explicit(near));
    let mut dec = Decider::default();
    for surface in [old, new] {
        let Outcome::ProvenCoincident(proof) = model.geom().side(surface, near, &mut dec) else {
            panic!("the explicit point is exactly on both slot cylinders");
        };
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(near),
            b: GeomKey::Surface(surface),
            at: None,
            origin: FactOrigin::Proof(proof.proof),
        });
    }
    assert!(matches!(
        model.geom().side(boss, near, &mut dec),
        Outcome::Proven(_)
    ));
    assert!(!model.geom().facts().holds(
        FactKind::Incident,
        GeomKey::Point(near),
        GeomKey::Surface(boss)
    ));
    let mut pool = Vec::new();
    pool.push(near);
    let (answer, work) = with_operation_budget(1_000_000, || {
        canonical_point(
            &mut model,
            &mut pool,
            incoming,
            &mut dec,
            "boolean.imprint-ambiguous",
        )
    });
    assert_eq!(answer.unwrap(), incoming);
    assert!(!work.exhausted);
    assert_eq!(pool, vec![near, incoming]);
    assert!(!model.geom().facts().holds(
        FactKind::Coincident,
        GeomKey::Point(near),
        GeomKey::Point(incoming)
    ));
}
