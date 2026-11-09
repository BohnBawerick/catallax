//! Finite spline spans share endpoints only after a complete separation proof.

use catallax_intersect::{BezierSectionPair, bezier_section_meetings};
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Rational};
use catallax_shapes::{
    BezierExtrusion, Curve, CurveId, GeomKey, GeomStore, Pose, RPoint3, RVector3, Rotation, Surface,
};

fn section(geom: &mut GeomStore, controls: &[[i64; 3]], dec: &mut Decider) -> CurveId {
    section_pose(geom, controls, None, dec)
}

fn section_pose(
    geom: &mut GeomStore,
    controls: &[[i64; 3]],
    pose: Option<&Pose>,
    dec: &mut Decider,
) -> CurveId {
    let support = BezierExtrusion::new(
        controls
            .iter()
            .map(|p| RPoint3::new(p[0].into(), p[1].into(), p[2].into()))
            .collect(),
        RVector3::new(0.into(), 0.into(), 2.into()),
        dec,
    )
    .unwrap();
    let support = pose.map_or_else(|| support.clone(), |pose| support.transformed(pose));
    let curve = geom.add_curve(Curve::BSpline(support.lower_edge()));
    let surface = geom.add_surface(Surface::BSpline(support.into()));
    geom.record_fact(Fact {
        kind: FactKind::Incident,
        a: GeomKey::Curve(curve),
        b: GeomKey::Surface(surface),
        at: None,
        origin: FactOrigin::Construction("test.bezier-section"),
    });
    curve
}

#[test]
fn adjacent_spans_and_layers_have_complete_endpoint_meetings() {
    let mut geom = GeomStore::new();
    let mut dec = Decider::default();
    let a = section(
        &mut geom,
        &[[0, 0, 0], [0, 2, 0], [1, 3, 0], [2, 3, 0]],
        &mut dec,
    );
    for (controls, count) in [
        ([[2, 3, 0], [3, 3, 0], [4, 2, 0], [4, 0, 0]], 1),
        ([[2, 3, 1], [3, 3, 1], [4, 2, 1], [4, 0, 1]], 0),
    ] {
        let b = section(&mut geom, &controls, &mut dec);
        for (first, second) in [(a, b), (b, a)] {
            let Some(BezierSectionPair::Points(points)) =
                bezier_section_meetings(&geom, first, second, &mut dec).unwrap()
            else {
                panic!("finite meetings expected");
            };
            assert_eq!(points.len(), count);
        }
    }
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn same_endpoint_arcs_are_distinct_and_exact_reversals_coincide() {
    let mut geom = GeomStore::new();
    let mut dec = Decider::default();
    let upper = section(
        &mut geom,
        &[[0, 0, 0], [0, 4, 0], [4, 4, 0], [4, 0, 0]],
        &mut dec,
    );
    let lower = section(
        &mut geom,
        &[[4, 0, 0], [4, -4, 0], [0, -4, 0], [0, 0, 0]],
        &mut dec,
    );
    let reversed = section(
        &mut geom,
        &[[4, 0, 0], [4, 4, 0], [0, 4, 0], [0, 0, 0]],
        &mut dec,
    );
    let Some(BezierSectionPair::Points(points)) =
        bezier_section_meetings(&geom, upper, lower, &mut dec).unwrap()
    else {
        panic!("two finite endpoints expected");
    };
    assert_eq!(points.len(), 2);
    assert!(matches!(
        bezier_section_meetings(&geom, upper, reversed, &mut dec).unwrap(),
        Some(BezierSectionPair::Coincident)
    ));
}

#[test]
fn interior_contact_is_refused_instead_of_omitted() {
    let mut geom = GeomStore::new();
    let mut dec = Decider::default();
    let arch = section(
        &mut geom,
        &[[0, 0, 0], [0, 4, 0], [4, 4, 0], [4, 0, 0]],
        &mut dec,
    );
    let tangent = section(
        &mut geom,
        &[[0, 3, 0], [1, 3, 0], [3, 3, 0], [4, 3, 0]],
        &mut dec,
    );
    let error = bezier_section_meetings(&geom, arch, tangent, &mut dec).unwrap_err();
    assert_eq!(error.code, "intersect.bezier-curve-pair");
    assert!(error.detail.contains("away from shared original endpoints"));
}

#[test]
fn subdivision_and_rigid_transport_preserve_complete_separation() {
    let pose = Pose::from_rational(
        Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap(),
        RVector3::new(7.into(), (-11).into(), 13.into()),
    );
    for pose in [None, Some(&pose)] {
        let mut geom = GeomStore::new();
        let mut dec = Decider::default();
        let arch = section_pose(
            &mut geom,
            &[[0, 0, 0], [0, 8, 0], [8, 8, 0], [8, 0, 0]],
            pose,
            &mut dec,
        );
        // The original hulls overlap in projection, but the arch reaches only y=6.
        let above = section_pose(
            &mut geom,
            &[[0, 7, 0], [1, 7, 0], [7, 7, 0], [8, 7, 0]],
            pose,
            &mut dec,
        );
        let Some(BezierSectionPair::Points(points)) =
            bezier_section_meetings(&geom, arch, above, &mut dec).unwrap()
        else {
            panic!("separated curves expected");
        };
        assert!(points.is_empty());
        assert_eq!(dec.stats().coincident_separation, 0);
        let left = section_pose(
            &mut geom,
            &[[0, 0, 0], [0, 2, 0], [1, 3, 0], [2, 3, 0]],
            pose,
            &mut dec,
        );
        let right = section_pose(
            &mut geom,
            &[[2, 3, 0], [3, 3, 0], [4, 2, 0], [4, 0, 0]],
            pose,
            &mut dec,
        );
        let Some(BezierSectionPair::Points(points)) =
            bezier_section_meetings(&geom, left, right, &mut dec).unwrap()
        else {
            panic!("shared endpoint expected");
        };
        assert_eq!(points.len(), 1);
        let expected = RPoint3::new(
            Rational::from_i64(2),
            Rational::from_i64(3),
            Rational::zero(),
        )
        .to_real();
        let expected = pose.map_or_else(|| expected.clone(), |pose| pose.apply_point(&expected));
        assert!(
            dec.sign(&points[0].minus(&expected).norm_squared())
                .is_coincident()
        );
    }
}

#[test]
fn section_queries_keep_neighbouring_wall_endpoints_and_finite_height() {
    use catallax_intersect::{Contact, CurveSurface, intersect_curve_surface};
    let mut geom = GeomStore::new();
    let mut dec = Decider::default();
    let left = section(
        &mut geom,
        &[[0, 0, 0], [0, 2, 0], [1, 3, 0], [2, 3, 0]],
        &mut dec,
    );
    for height in [0, 1] {
        let support = BezierExtrusion::new(
            vec![
                RPoint3::new(2.into(), 3.into(), height.into()),
                RPoint3::new(3.into(), 3.into(), height.into()),
                RPoint3::new(4.into(), 2.into(), height.into()),
                RPoint3::new(4.into(), 0.into(), height.into()),
            ],
            RVector3::new(0.into(), 0.into(), 2.into()),
            &mut dec,
        )
        .unwrap();
        let wall = geom.add_surface(Surface::BSpline(support.into()));
        let answer = intersect_curve_surface(&mut geom, left, wall, &mut dec).unwrap();
        assert!(answer.marks.is_empty());
        if height == 0 {
            let CurveSurface::Meets(points) = answer.value else {
                panic!("one finite wall endpoint expected");
            };
            assert_eq!(points.len(), 1);
            assert_eq!(points[0].contact, Contact::Tangent);
        } else {
            assert_eq!(answer.value, CurveSurface::Disjoint);
        }
    }
}
