//! Certified algebraic dimensions, rigid motions and domain-validation precedence.
use catallax_numbers::{
    Decider, Outcome, Rational, Real, StrictOrder, StrictSign, sin_cos_degrees,
};
use catallax_shapes::{
    Circle, Curve, Cylinder, GeomStore, Point3, Pose, Rotation, Surface, Vector3,
};

fn real(n: i64) -> Real {
    Real::from_i64(n)
}
fn point([x, y, z]: [i64; 3]) -> Point3 {
    Point3::new(real(x), real(y), real(z))
}
fn vector([x, y, z]: [i64; 3]) -> Vector3 {
    Vector3::new(real(x), real(y), real(z))
}

#[test]
fn algebraic_dimensions_preserve_incidence_gradient_and_implicit_scale() {
    let mut dec = Decider::default();
    let root = real(2).sqrt();
    let cylinder = Surface::Cylinder(
        Cylinder::new(
            Point3::origin(),
            Vector3::new(real(0), real(0), root.clone()),
            root.clone(),
            vector([1, 0, 0]),
            &mut dec,
        )
        .unwrap(),
    );
    let implicit = cylinder.implicit();
    // Axis length squared is 2, so the unscaled carrier is 2(x²+y²)-4.
    for (p, expected) in [([1, 1, 5], 0), ([1, 0, 5], -2), ([2, 0, 5], 4)] {
        let p = p.map(real);
        assert!(
            dec.sign(&(&implicit.eval_real(&p) - &real(expected)))
                .is_coincident()
        );
    }
    let gradient = cylinder.gradient(&point([1, 1, 5]));
    for (got, expected) in [&gradient.x, &gradient.y, &gradient.z]
        .into_iter()
        .zip([4, 4, 0])
    {
        assert!(dec.compare(got, &real(expected)).is_coincident());
    }
    // The squared radius itself may now be irrational, not just its square root.
    let circle = Curve::Circle(
        Circle::new(
            Point3::origin(),
            vector([0, 0, 1]),
            root.clone(),
            vector([1, 0, 0]),
            &mut dec,
        )
        .unwrap(),
    );
    let p = [root.sqrt(), real(0), real(0)];
    for equation in circle.equations(&GeomStore::new()) {
        assert!(dec.sign(&equation.eval_real(&p)).is_coincident());
    }
}

#[test]
fn nonquadratic_rational_degree_rotation_has_the_selected_exact_axis() {
    let mut dec = Decider::default();
    let [sine, cosine] = sin_cos_degrees(&Rational::from_i64(20), &mut dec).unwrap();
    let rotation =
        Rotation::from_real_quaternion(cosine, real(0), sine, real(0), &mut dec).unwrap();
    let pose = Pose::new(rotation, Vector3::zero(), &mut dec).unwrap();
    let axis = pose.apply_vector(&vector([0, 0, 1]));
    // z = cos(40°) is the unique root of 8z³-6z+1 in (3/4,4/5).
    let cubic = &(&(&axis.z.pow_u32(3) * &real(8)) - &(&axis.z * &real(6))) + &real(1);
    assert!(dec.sign(&cubic).is_coincident());
    assert!(matches!(
        dec.compare(&axis.z, &Real::from(Rational::from_ratio(3, 4).unwrap())),
        Outcome::Proven(StrictOrder::Greater)
    ));
    assert!(matches!(
        dec.compare(&axis.z, &Real::from(Rational::from_ratio(4, 5).unwrap())),
        Outcome::Proven(StrictOrder::Less)
    ));
    assert!(matches!(
        dec.sign(&axis.x),
        Outcome::Proven(StrictSign::Positive)
    ));
    assert!(dec.sign(&axis.y).is_coincident());
}

#[test]
fn defining_expression_domains_precede_geometric_invariant_failures() {
    let mut dec = Decider::default();
    let bad_origin = Point3::new(&real(1) / &real(0), real(0), real(0));
    let error = Cylinder::new(
        bad_origin,
        Vector3::zero(),
        real(-1),
        Vector3::zero(),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "numbers.division-by-zero");
    let error = Circle::new(
        Point3::origin(),
        Vector3::zero(),
        real(-1).sqrt(),
        Vector3::zero(),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "numbers.sqrt-of-negative");
}
