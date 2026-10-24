#![allow(missing_docs)]

use catallax_numbers::{Decider, Real};
use catallax_shapes::{Cone, Cylinder, Plane, Point3, Sphere, Surface, Torus, Vector3};

fn p(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(x.into(), y.into(), z.into())
}
fn v(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(x.into(), y.into(), z.into())
}
fn zero(value: Real, dec: &mut Decider) {
    assert!(dec.sign(&value).is_coincident());
}

#[test]
fn nonunit_plane_normals_move_by_the_metric_distance() {
    let mut dec = Decider::default();
    let source = Surface::Plane(Plane::new(p(2, 0, 0), v(3, 4, 0), v(0, 0, 2), &mut dec).unwrap());
    let Some(Surface::Plane(parallel)) = source
        .parallel_support(&Real::from_i64(5), &mut dec)
        .unwrap()
    else {
        panic!("plane vanished");
    };
    zero(parallel.origin.x - Real::from_i64(5), &mut dec);
    zero(parallel.origin.y - Real::from_i64(4), &mut dec);
    zero(parallel.origin.z, &mut dec);
}

#[test]
fn cylinder_and_sphere_collapses_do_not_reappear_with_a_squared_negative_radius() {
    let mut dec = Decider::default();
    let surfaces = [
        Surface::Cylinder(
            Cylinder::new(
                p(0, 0, 0),
                v(0, 0, 2),
                Real::from_i64(2),
                v(1, 0, 0),
                &mut dec,
            )
            .unwrap(),
        ),
        Surface::Sphere(
            Sphere::new(
                p(0, 0, 0),
                Real::from_i64(2),
                v(0, 0, 2),
                v(1, 0, 0),
                &mut dec,
            )
            .unwrap(),
        ),
    ];
    for surface in surfaces {
        for distance in [-2, -3] {
            assert!(
                surface
                    .parallel_support(&Real::from_i64(distance), &mut dec)
                    .unwrap()
                    .is_none()
            );
        }
    }
}

#[test]
fn cone_parallel_is_an_exact_translated_nappe_with_the_original_slope() {
    let mut dec = Decider::default();
    let source = Surface::Cone(
        Cone::new(
            p(0, 0, 0),
            v(0, 0, 7),
            Real::from_i64(3) / Real::from_i64(4),
            v(1, 0, 0),
            &mut dec,
        )
        .unwrap(),
    );
    let Some(Surface::Cone(parallel)) = source
        .parallel_support(&Real::from_i64(3), &mut dec)
        .unwrap()
    else {
        panic!("cone vanished");
    };
    // sqrt(1 + (3/4)^2) / (3/4) = 5/3, so its apex moves -5.
    zero(parallel.apex.z + Real::from_i64(5), &mut dec);
    zero(
        parallel.slope - Real::from_i64(3) / Real::from_i64(4),
        &mut dec,
    );
}

#[test]
fn ring_torus_parallel_preserves_the_core_and_refuses_a_false_ring_carrier() {
    let mut dec = Decider::default();
    let source = Surface::Torus(
        Torus::new(
            p(1, 2, 3),
            v(0, 0, 2),
            Real::from_i64(5),
            Real::one(),
            v(1, 0, 0),
            &mut dec,
        )
        .unwrap(),
    );
    let Some(Surface::Torus(parallel)) = source
        .parallel_support(&Real::from_i64(2), &mut dec)
        .unwrap()
    else {
        panic!("torus vanished");
    };
    zero(parallel.major - Real::from_i64(5), &mut dec);
    zero(parallel.minor - Real::from_i64(3), &mut dec);
    assert_eq!(
        source
            .parallel_support(&Real::from_i64(4), &mut dec)
            .unwrap_err()
            .code,
        "offset.spindle-carrier"
    );
}
