//! Finite spline section geometry and pair orientation.

use catallax_intersect::{Component, Contact, SurfacePair, intersect_surfaces};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{BezierExtrusion, Curve, GeomStore, Plane, RPoint3, RVector3, Surface};

fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

#[test]
fn normal_sections_keep_exact_controls_and_reverse_with_pair_order() {
    for height in [-2, 2] {
        for sign in [-1, 1] {
            for swap in [false, true] {
                for z in [-3, -2, -1, 0, 1, 2, 3] {
                    let mut dec = Decider::default();
                    let support = BezierExtrusion::new(
                        vec![p(0, 0, 0), p(0, 4, 0), p(4, 4, 0), p(4, 0, 0)],
                        v(0, 0, height),
                        &mut dec,
                    )
                    .unwrap();
                    let mut geom = GeomStore::new();
                    let wall = geom.add_surface(Surface::BSpline(support.into()));
                    let plane = geom.add_surface(Surface::Plane(
                        Plane::new(
                            p(0, 0, z).to_real(),
                            v(0, 0, sign).to_real(),
                            v(1, 0, 0).to_real(),
                            &mut dec,
                        )
                        .unwrap(),
                    ));
                    let (a, b) = if swap { (wall, plane) } else { (plane, wall) };
                    let answer = intersect_surfaces(&mut geom, a, b, &mut dec).unwrap();
                    assert!(answer.marks.is_empty());
                    if z < height.min(0) || z > height.max(0) {
                        assert_eq!(answer.value, SurfacePair::Disjoint);
                        continue;
                    }
                    let SurfacePair::Meets(parts) = answer.value else {
                        panic!("section expected")
                    };
                    let [
                        Component::Curve {
                            curve,
                            contact: Contact::Crossing,
                            singular,
                        },
                    ] = parts.as_slice()
                    else {
                        panic!("one transverse curve expected")
                    };
                    assert!(singular.is_empty());
                    let Curve::BSpline(curve) = geom.curve(*curve) else {
                        panic!("polynomial section expected")
                    };
                    let reverse = (sign * height < 0) != swap;
                    let expected = if reverse { p(4, 0, z) } else { p(0, 0, z) };
                    let first = curve.point_at(&Rational::zero().into(), &mut dec).unwrap();
                    assert!(
                        dec.sign(&first.minus(&expected.to_real()).norm_squared())
                            .is_coincident()
                    );
                    for point in curve.controls() {
                        assert!(dec.sign(&(&point.z - &z.into())).is_coincident());
                    }
                }
            }
        }
    }
}

#[test]
fn oblique_sections_keep_a_named_refusal() {
    let mut dec = Decider::default();
    let support = BezierExtrusion::new(
        vec![p(0, 0, 0), p(2, 4, 0), p(4, 0, 0)],
        v(0, 0, 2),
        &mut dec,
    )
    .unwrap();
    let mut geom = GeomStore::new();
    let wall = geom.add_surface(Surface::BSpline(support.into()));
    let plane = geom.add_surface(Surface::Plane(
        Plane::new(
            p(0, 0, 0).to_real(),
            v(1, 0, 1).to_real(),
            v(0, 1, 0).to_real(),
            &mut dec,
        )
        .unwrap(),
    ));
    assert_eq!(
        intersect_surfaces(&mut geom, plane, wall, &mut dec)
            .unwrap_err()
            .code,
        "intersect.bezier-plane-oblique"
    );
}

#[test]
fn strict_corner_coefficients_exclude_a_plane_touching_only_the_control_hull() {
    for y in [-4, 4] {
        let mut dec = Decider::default();
        let support = BezierExtrusion::new(
            vec![p(0, 0, 0), p(0, 4, 0), p(4, 4, 0), p(4, 0, 0)],
            v(0, 0, 2),
            &mut dec,
        )
        .unwrap();
        let mut geom = GeomStore::new();
        let wall = geom.add_surface(Surface::BSpline(support.into()));
        let plane = geom.add_surface(Surface::Plane(
            Plane::new(
                p(0, y, 0).to_real(),
                v(0, 1, 0).to_real(),
                v(1, 0, 0).to_real(),
                &mut dec,
            )
            .unwrap(),
        ));
        assert_eq!(
            intersect_surfaces(&mut geom, plane, wall, &mut dec)
                .unwrap()
                .value,
            SurfacePair::Disjoint
        );
    }
}
