//! Exact curve queries retain their finite section proof.

use catallax_intersect::{
    Component, Contact, CurveSurface, SurfacePair, intersect_curve_surface, intersect_surfaces,
    same_curve,
};
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Outcome, Rational, ZeroProof};
use catallax_shapes::{
    Apart, BezierExtrusion, Curve, CurveId, GeomKey, GeomStore, Line, Plane, RPoint3, RVector3,
    Surface, SurfaceId,
};

fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}
fn plane(
    geom: &mut GeomStore,
    origin: RPoint3,
    normal: RVector3,
    reference: RVector3,
    dec: &mut Decider,
) -> SurfaceId {
    geom.add_surface(Surface::Plane(
        Plane::new(origin.to_real(), normal.to_real(), reference.to_real(), dec).unwrap(),
    ))
}
fn wall(geom: &mut GeomStore, dec: &mut Decider) -> SurfaceId {
    geom.add_surface(Surface::BSpline(
        BezierExtrusion::new(
            vec![p(0, 0, 0), p(0, 4, 0), p(4, 4, 0), p(4, 0, 0)],
            v(0, 0, 2),
            dec,
        )
        .unwrap()
        .into(),
    ))
}
fn section(
    geom: &mut GeomStore,
    wall: SurfaceId,
    height: i64,
    flip: bool,
    dec: &mut Decider,
) -> CurveId {
    let plane = plane(
        geom,
        p(0, 0, height),
        v(0, 0, if flip { -1 } else { 1 }),
        v(1, 0, 0),
        dec,
    );
    let SurfacePair::Meets(parts) = intersect_surfaces(geom, plane, wall, dec).unwrap().value
    else {
        panic!("section expected")
    };
    let Component::Curve { curve, .. } = parts[0] else {
        panic!("curve expected")
    };
    curve
}

#[test]
fn repeated_and_reversed_sections_have_exact_point_set_equality() {
    let mut dec = Decider::default();
    let mut geom = GeomStore::new();
    let wall = wall(&mut geom, &mut dec);
    let first = section(&mut geom, wall, 1, false, &mut dec);
    let same = section(&mut geom, wall, 1, false, &mut dec);
    let reversed = section(&mut geom, wall, 1, true, &mut dec);
    let lower = section(&mut geom, wall, 0, false, &mut dec);
    assert!(same_curve(&geom, first, same, &mut dec).is_coincident());
    assert!(same_curve(&geom, first, reversed, &mut dec).is_coincident());
    assert!(matches!(
        same_curve(&geom, first, lower, &mut dec),
        Outcome::Proven(Apart)
    ));
    let line = geom.add_curve(Curve::Line(
        Line::new(p(0, 0, 1).to_real(), v(1, 0, 0).to_real(), &mut dec).unwrap(),
    ));
    assert!(matches!(
        same_curve(&geom, first, line, &mut dec),
        Outcome::Proven(Apart)
    ));
    assert!(matches!(
        same_curve(&geom, line, first, &mut dec),
        Outcome::Proven(Apart)
    ));
    assert!(matches!(
        intersect_curve_surface(&mut geom, first, wall, &mut dec)
            .unwrap()
            .value,
        CurveSurface::Contained
    ));
}

#[test]
fn plane_roots_keep_crossings_tangencies_and_finite_endpoints() {
    let mut dec = Decider::default();
    let mut geom = GeomStore::new();
    let wall = wall(&mut geom, &mut dec);
    let curve = section(&mut geom, wall, 1, false, &mut dec);
    for (origin, normal, reference, count, contact) in [
        (p(2, 0, 0), v(1, 0, 0), v(0, 1, 0), 1, Contact::Crossing),
        (p(0, 3, 0), v(0, 1, 0), v(1, 0, 0), 1, Contact::Tangent),
        (p(0, 2, 0), v(0, 1, 0), v(1, 0, 0), 2, Contact::Crossing),
        (p(0, 0, 0), v(1, 0, 0), v(0, 1, 0), 1, Contact::Tangent),
    ] {
        let plane = plane(&mut geom, origin, normal, reference, &mut dec);
        let answer = intersect_curve_surface(&mut geom, curve, plane, &mut dec).unwrap();
        assert!(answer.marks.is_empty());
        let CurveSurface::Meets(hits) = answer.value else {
            panic!("finite hits expected")
        };
        assert_eq!(hits.len(), count);
        for hit in hits {
            assert_eq!(hit.contact, contact);
            let point = geom.position(hit.point);
            assert!(dec.sign(&(&point.z - &1.into())).is_coincident());
        }
    }
    for (height, contained) in [(1, true), (0, false)] {
        let plane = plane(&mut geom, p(0, 0, height), v(0, 0, 1), v(1, 0, 0), &mut dec);
        assert_eq!(
            intersect_curve_surface(&mut geom, curve, plane, &mut dec)
                .unwrap()
                .value,
            if contained {
                CurveSurface::Contained
            } else {
                CurveSurface::Disjoint
            }
        );
    }
}

#[test]
fn incidence_without_matching_controls_retains_the_native_curve_definition() {
    let mut dec = Decider::default();
    let mut geom = GeomStore::new();
    let wall = wall(&mut geom, &mut dec);
    let false_curve = BezierExtrusion::new(
        vec![p(0, 0, 0), p(0, 5, 0), p(4, 5, 0), p(4, 0, 0)],
        v(0, 0, 2),
        &mut dec,
    )
    .unwrap()
    .lower_edge();
    let curve = geom.add_curve(Curve::BSpline(false_curve));
    geom.record_fact(Fact {
        kind: FactKind::Incident,
        a: GeomKey::Curve(curve),
        b: GeomKey::Surface(wall),
        at: None,
        origin: FactOrigin::Construction("test.false-incidence"),
    });
    let plane = plane(&mut geom, p(2, 0, 0), v(1, 0, 0), v(0, 1, 0), &mut dec);
    let answer = intersect_curve_surface(&mut geom, curve, plane, &mut dec).unwrap();
    assert!(answer.marks.is_empty());
    let CurveSurface::Meets(hits) = answer.value else {
        panic!("the native curve needs its complete finite crossing")
    };
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].contact, Contact::Crossing);
    let point = geom.position(hits[0].point);
    let expected = RPoint3::new(2.into(), Rational::from_ratio(15, 4).unwrap(), 0.into()).to_real();
    for (actual, expected) in [
        (&point.x, &expected.x),
        (&point.y, &expected.y),
        (&point.z, &expected.z),
    ] {
        assert!(matches!(
            dec.sign(&(actual - expected)),
            Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
        ));
    }
}

#[test]
fn finite_section_line_meetings_include_tangent_and_generator_hits() {
    let mut dec = Decider::default();
    let mut geom = GeomStore::new();
    let wall = wall(&mut geom, &mut dec);
    let curve = section(&mut geom, wall, 1, false, &mut dec);
    for (origin, direction, count) in [
        (p(0, 0, 1), v(1, 0, 0), 2),
        (p(0, 3, 1), v(1, 0, 0), 1),
        (p(0, 4, 1), v(1, 0, 0), 0),
        (p(2, 3, 0), v(0, 0, 1), 1),
        (p(2, 2, 0), v(0, 0, 1), 0),
        (p(0, 2, 0), v(1, 0, 0), 0),
        (p(0, 0, 0), v(0, 0, 1), 1),
    ] {
        let line = geom.add_curve(Curve::Line(
            Line::new(origin.to_real(), direction.to_real(), &mut dec).unwrap(),
        ));
        let points = catallax_intersect::bezier_section_line_points(&geom, curve, line, &mut dec)
            .unwrap()
            .unwrap();
        assert_eq!(points.len(), count);
        for point in points {
            assert!(dec.sign(&(&point.z - &1.into())).is_coincident());
        }
    }
}

#[test]
fn section_normal_acceleration_matches_independent_arch_curvatures() {
    use catallax_numbers::{Rational, Real};
    for reverse in [false, true] {
        let mut dec = Decider::default();
        let mut geom = GeomStore::new();
        let wall = wall(&mut geom, &mut dec);
        let curve = section(&mut geom, wall, 1, reverse, &mut dec);
        for (position, tangent, expected) in [
            (p(0, 0, 1), v(0, 1, 0), ["1/6", "0", "0"]),
            (p(2, 3, 1), v(1, 0, 0), ["0", "-2/3", "0"]),
            (p(4, 0, 1), v(0, -1, 0), ["-1/6", "0", "0"]),
        ] {
            for scale in [-3, 1, 2] {
                let tangent = tangent.to_real().scaled(&Real::from_i64(scale));
                let expected = expected.map(|x| Rational::parse(x).unwrap());
                let expected = RVector3::new(
                    expected[0].clone(),
                    expected[1].clone(),
                    expected[2].clone(),
                )
                .to_real()
                .scaled(&Real::from_i64(scale * scale));
                let actual = catallax_intersect::bezier_section_normal_acceleration(
                    &geom,
                    curve,
                    &position.to_real(),
                    &tangent,
                    &mut dec,
                )
                .unwrap()
                .unwrap();
                assert!(
                    dec.sign(&actual.minus(&expected).norm_squared())
                        .is_coincident()
                );
            }
        }
        assert_eq!(
            catallax_intersect::bezier_section_normal_acceleration(
                &geom,
                curve,
                &p(2, 3, 1).to_real(),
                &v(0, 1, 0).to_real(),
                &mut dec
            )
            .unwrap_err()
            .code,
            "intersect.bezier-tangent-mismatch"
        );
        assert_eq!(
            catallax_intersect::bezier_section_normal_acceleration(
                &geom,
                curve,
                &p(6, 3, 1).to_real(),
                &v(1, 0, 0).to_real(),
                &mut dec
            )
            .unwrap_err()
            .code,
            "intersect.bezier-point-off-curve"
        );
    }
}
