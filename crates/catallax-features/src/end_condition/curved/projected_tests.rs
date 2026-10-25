use super::*;

#[test]
fn quartic_torus_ray_cannot_be_truncated_to_a_quadratic() {
    let mut dec = Decider::default();
    let torus = Surface::Torus(
        catallax_shapes::Torus::new(
            Point3::origin(),
            axis(2),
            Real::from_i64(3),
            Real::one(),
            axis(0),
            &mut dec,
        )
        .unwrap(),
    );
    let frame = Frame {
        origin: Point3::new(Real::from_i64(3), Real::zero(), Real::from_i64(-5)),
        u: axis(0),
        v: axis(1),
        direction: axis(2),
    };
    let error = coefficients(&torus, &frame).unwrap_err();
    assert_eq!(error.code, "features.unsupported-curved-limit");
}
