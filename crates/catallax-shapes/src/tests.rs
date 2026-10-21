use super::*;
use catallax_numbers::{
    Ball, BoxRoots, BudgetPolicy, Coincidence, Decider, Fact, FactKind, FactOrigin, MPoly, Outcome,
    PolySystem, Rational, Real, RootBox, Sign, StrictSign, ZeroProof,
};
use std::{
    hash::{DefaultHasher, Hash, Hasher},
    sync::Arc,
    time::{Duration, Instant},
};

fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}
fn b(x: i64, precision: u32) -> Ball {
    Ball::from_rational(&x.into(), precision)
}
fn quaternion(q: [i64; 4]) -> Rotation {
    let [w, x, y, z] = q.map(Rational::from_i64);
    Rotation::from_quaternion(w, x, y, z).unwrap()
}
fn shapes() -> [Surface; 5] {
    [
        Surface::Plane(Plane::new(p(0, 0, 1), v(0, 0, 2), v(1, 0, 0)).unwrap()),
        Surface::Cylinder(Cylinder::new(p(0, 0, 0), v(0, 0, 2), 5.into(), v(3, 0, 0)).unwrap()),
        Surface::Cone(Cone::new(p(0, 0, 0), v(0, 0, 1), r("1/2"), v(1, 0, 0)).unwrap()),
        Surface::Sphere(Sphere::new(p(0, 0, 0), 3.into(), v(0, 0, 1), v(1, 0, 0)).unwrap()),
        Surface::Torus(Torus::new(p(0, 0, 0), v(0, 0, 1), 3.into(), 1.into(), v(1, 0, 0)).unwrap()),
    ]
}
fn explicit(store: &mut GeomStore, point: RPoint3) -> PointId {
    store.add_point(PointDef::Explicit(point))
}
fn contains(actual: &[Ball; 3], expected: &RPoint3) {
    for (a, e) in actual.iter().zip([&expected.x, &expected.y, &expected.z]) {
        assert!(a.contains(e), "{a:?} excludes {e}");
    }
}
fn coincident<T: std::fmt::Debug>(outcome: Outcome<T>) {
    assert!(outcome.is_coincident(), "{outcome:?}");
}

#[test]
fn quaternion_matrix_composition_and_pose_order() {
    for n in [1, 2, -3] {
        assert_eq!(quaternion([n, 0, 0, 0]), Rotation::identity());
    }
    let turn = quaternion([1, 0, 0, 1]);
    assert_eq!(turn.apply(&v(1, 0, 0)), v(0, 1, 0));
    assert_eq!(quaternion([1, 1, 1, 1]).apply(&v(0, 1, 0)), v(0, 0, 1));
    let q = quaternion([1, 2, 3, 4]);
    assert_eq!(
        q.matrix(),
        &[
            [r("-2/3"), r("2/15"), r("11/15")],
            [r("2/3"), r("-1/3"), r("2/3")],
            [r("1/3"), r("14/15"), r("2/15")]
        ]
    );
    assert_eq!(
        q.then(&quaternion([2, -1, 0, 5])),
        quaternion([-16, -12, 20, 10])
    );
    let hash = |a: &Rotation| {
        let mut h = DefaultHasher::new();
        a.hash(&mut h);
        h.finish()
    };
    for scale in [r("-1"), r("5"), r("1/7")] {
        let [w, x, y, z] = [1, 2, 3, 4].map(|x| Rational::from_i64(x) * scale.clone());
        let scaled = Rotation::from_quaternion(w, x, y, z).unwrap();
        assert_eq!(scaled, q);
        assert_eq!(hash(&scaled), hash(&q));
    }
    let a = Pose::new(turn, v(1, 2, 3));
    let next = Pose::new(quaternion([2, 1, 0, 0]), v(-7, 2, 4));
    assert_eq!(a.apply_point(&p(1, 0, 0)), p(1, 3, 3));
    assert_eq!(a.inverse().apply_point(&p(1, 3, 3)), p(1, 0, 0));
    assert_eq!(
        a.then(&next).apply_point(&p(2, 5, 7)),
        next.apply_point(&a.apply_point(&p(2, 5, 7)))
    );
    assert_eq!(a.then(&a.inverse()), Pose::identity());
    assert_eq!(a.apply_vector(&v(1, 0, 0)), v(0, 1, 0));
    assert_eq!(
        Rotation::from_quaternion(0.into(), 0.into(), 0.into(), 0.into())
            .unwrap_err()
            .code,
        "shapes.zero-quaternion"
    );
    for i in 1..=100 {
        let rotation = quaternion([i, i - 3, 2 * i, 7 - i]);
        assert_eq!(rotation.then(&rotation.inverse()), Rotation::identity());
        let x = rotation.apply(&v(1, 0, 0));
        let y = rotation.apply(&v(0, 1, 0));
        let z = rotation.apply(&v(0, 0, 1));
        assert_eq!(x.dot(&y.cross(&z)), Rational::one());
    }
}

#[test]
fn constructor_refusal_precedence() {
    assert_eq!(
        Plane::new(p(0, 0, 0), v(0, 0, 0), v(1, 0, 1))
            .unwrap_err()
            .code,
        "shapes.zero-vector"
    );
    assert_eq!(
        Plane::new(p(0, 0, 0), v(0, 0, 1), v(1, 0, 1))
            .unwrap_err()
            .code,
        "shapes.not-perpendicular"
    );
    assert_eq!(
        Cylinder::new(p(0, 0, 0), v(0, 0, 0), (-1).into(), v(1, 0, 1))
            .unwrap_err()
            .code,
        "shapes.zero-vector"
    );
    for radius in [0, -1] {
        assert_eq!(
            Cylinder::new(p(0, 0, 0), v(0, 0, 1), radius.into(), v(1, 0, 1))
                .unwrap_err()
                .code,
            "shapes.non-positive-radius"
        );
    }
    assert_eq!(
        Cone::new(p(0, 0, 0), v(0, 0, 1), 0.into(), v(1, 0, 0))
            .unwrap_err()
            .code,
        "shapes.non-positive-slope"
    );
    assert_eq!(
        Sphere::new(p(0, 0, 0), 0.into(), v(0, 0, 1), v(1, 0, 0))
            .unwrap_err()
            .code,
        "shapes.non-positive-radius"
    );
    assert_eq!(
        Torus::new(p(0, 0, 0), v(0, 0, 1), 0.into(), 1.into(), v(1, 0, 0))
            .unwrap_err()
            .code,
        "shapes.non-positive-radius"
    );
    for minor in [3, 4] {
        assert_eq!(
            Torus::new(p(0, 0, 0), v(0, 0, 1), 3.into(), minor.into(), v(1, 0, 1))
                .unwrap_err()
                .code,
            "shapes.torus-self-intersecting"
        );
    }
    assert_eq!(
        Line::new(p(0, 0, 0), v(0, 0, 0)).unwrap_err().code,
        "shapes.zero-vector"
    );
    assert_eq!(
        Circle::new(p(0, 0, 0), v(0, 0, 1), 0.into(), v(1, 0, 0))
            .unwrap_err()
            .code,
        "shapes.non-positive-radius"
    );
    for minor in [0, 25, 26] {
        assert_eq!(
            Ellipse::new(p(0, 0, 0), v(0, 0, 1), v(2, 0, 0), 25.into(), minor.into())
                .unwrap_err()
                .code,
            "shapes.ellipse-axes"
        );
    }
}

#[test]
fn exact_polynomials_and_gradients_under_poses() {
    let surfaces = shapes();
    let points = [p(5, 7, 1), p(3, 4, 7), p(3, 4, -10), p(1, 2, 2), p(4, 0, 0)];
    for (surface, point) in surfaces.iter().zip(points.iter()) {
        assert_eq!(
            surface
                .implicit()
                .eval(&[point.x.clone(), point.y.clone(), point.z.clone()]),
            Rational::zero()
        );
    }
    assert_eq!(
        surfaces[1].implicit().eval(&[1.into(), 1.into(), 1.into()]),
        Rational::from_i64(-92)
    );
    assert_eq!(
        surfaces[4].implicit().eval(&[0.into(), 0.into(), 0.into()]),
        Rational::from_i64(64)
    );
    assert_eq!(
        surfaces[4].implicit().eval(&[3.into(), 0.into(), 0.into()]),
        Rational::from_i64(-35)
    );
    for i in 1..=20 {
        let pose = Pose::new(quaternion([i, 2, 3, 4]), v(1_000_000, -7, i));
        let point = p(i - 5, 2 * i - 7, 3 - i);
        for surface in &surfaces {
            let moved = surface.transformed(&pose);
            let image = pose.apply_point(&point);
            assert_eq!(
                moved.implicit().eval(&[image.x, image.y, image.z]),
                surface
                    .implicit()
                    .eval(&[point.x.clone(), point.y.clone(), point.z.clone()])
            );
            let gradient = surface.gradient(&point.to_real());
            for (axis, coordinate) in [gradient.x, gradient.y, gradient.z].into_iter().enumerate() {
                assert_eq!(
                    coordinate.as_rational().unwrap(),
                    surface.implicit().partial(axis).eval(&[
                        point.x.clone(),
                        point.y.clone(),
                        point.z.clone()
                    ])
                );
            }
        }
    }
    let g = surfaces[3].gradient(&p(1, 2, 2).to_real());
    assert_eq!(
        [g.x.as_rational(), g.y.as_rational(), g.z.as_rational()],
        [Some(2.into()), Some(4.into()), Some(4.into())]
    );
}

#[test]
fn sides_near_carriers_and_scaled_placed_families() {
    let surfaces = shapes();
    let on = [
        p(5, 7, 1),
        p(3, 4, 7),
        p(3, 4, 10),
        p(1, 2, 2),
        RPoint3::new(r("54/25"), r("72/25"), r("4/5")),
    ];
    let mut store = GeomStore::new();
    let mut dec = Decider::default();
    for (surface, point) in surfaces.iter().zip(on.iter()) {
        let s = store.add_surface(surface.clone());
        coincident(store.side_of(s, &point.to_real(), &mut dec));
    }
    // A fixed rational family covers 30,000 exact signs, including scale and placement.
    for i in 0..6000_i64 {
        let scale = r(&format!("1e{}", i % 61 - 30));
        let e = r(if i % 2 == 0 { "1e-30" } else { "-1e-30" });
        let pose = Pose::new(
            quaternion([1 + i % 11, 2, 3, 4]),
            v(i % 17, -i % 19, i % 23),
        );
        for (kind, surface) in surfaces.iter().enumerate() {
            let mut surface = surface.clone();
            match &mut surface {
                Surface::Plane(s) => {
                    s.origin =
                        RPoint3::origin().plus(&s.origin.minus(&RPoint3::origin()).scaled(&scale))
                }
                Surface::Cylinder(s) => s.radius = &s.radius * &scale,
                Surface::Cone(_) => {}
                Surface::Sphere(s) => s.radius = &s.radius * &scale,
                Surface::Torus(s) => {
                    s.major = &s.major * &scale;
                    s.minor = &s.minor * &scale;
                }
            }
            let local = match kind {
                0 => RPoint3::new(5.into(), 7.into(), Rational::one() + e.clone()),
                1 => RPoint3::new(
                    Rational::from_i64(3) * (Rational::one() + e.clone()),
                    Rational::from_i64(4) * (Rational::one() + e.clone()),
                    7.into(),
                ),
                2 => RPoint3::new(
                    3.into(),
                    4.into(),
                    Rational::from_i64(10) * (Rational::one() - e.clone()),
                ),
                3 => RPoint3::new(
                    Rational::one() + e.clone(),
                    Rational::from_i64(2) * (Rational::one() + e.clone()),
                    Rational::from_i64(2) * (Rational::one() + e.clone()),
                ),
                _ => RPoint3::new(
                    r("54/25"),
                    r("72/25"),
                    r("4/5") * (Rational::one() + e.clone()),
                ),
            };
            let point = RPoint3::origin().plus(&local.minus(&RPoint3::origin()).scaled(&scale));
            let mut temporary = GeomStore::new();
            let s = temporary.add_surface(surface.transformed(&pose));
            let expected = if e.sign() == Sign::Positive {
                Side::Front
            } else {
                Side::Back
            };
            assert_eq!(
                temporary.side_of(s, &pose.apply_point(&point).to_real(), &mut dec),
                Outcome::Proven(expected),
                "family {i} kind {kind}"
            );
        }
    }
}

#[test]
fn curve_equations_origins_and_membership() {
    let mut store = GeomStore::new();
    let line = Curve::Line(Line::new(p(1, 2, 3), v(2, -1, 0)).unwrap());
    let x = MPoly::var(3, 0);
    let y = MPoly::var(3, 1);
    let z = MPoly::var(3, 2);
    assert_eq!(
        line.equations(&store),
        vec![
            &z - &MPoly::constant(3, 3.into()),
            z.scale(&2.into()) - MPoly::constant(3, 6.into()),
            -x - y.scale(&2.into()) + MPoly::constant(3, 5.into())
        ]
    );
    let curves = [
        line,
        Curve::Circle(Circle::new(p(0, 0, 0), v(0, 0, 1), 25.into(), v(1, 0, 0)).unwrap()),
        Curve::Ellipse(
            Ellipse::new(p(0, 0, 0), v(0, 0, 1), v(2, 0, 0), 25.into(), 9.into()).unwrap(),
        ),
    ];
    let on = [
        p(5, 0, 3),
        p(3, 4, 0),
        RPoint3::new(3.into(), r("12/5"), 0.into()),
    ];
    let mut dec = Decider::default();
    for (curve, point) in curves.into_iter().zip(on) {
        let origin = curve.parameter_origin(&store);
        let c = store.add_curve(curve);
        let plane = store.add_surface(Surface::Plane(
            Plane::new(
                RPoint3::new(0.into(), 0.into(), point.z.clone()),
                v(0, 0, 1),
                v(1, 0, 0),
            )
            .unwrap(),
        ));
        let o = store.add_point(PointDef::Meet(Meet {
            of: MeetOf::CurveSurface(c, plane),
            position: Position::Closed(origin),
        }));
        // Use a distinct explicit point to exercise numeric membership, not its construction fact.
        let q = explicit(&mut store, point.clone());
        coincident(store.point_on_curve(c, q, &mut dec));
        let off = explicit(
            &mut store,
            point.plus(&RVector3::new(r("1e-30"), r("1e-30"), r("1e-30"))),
        );
        assert_eq!(
            store.point_on_curve(c, off, &mut dec),
            Outcome::Proven(Apart)
        );
        coincident(store.point_on_curve(c, o, &mut dec));
    }
}

#[test]
fn closed_and_isolated_positions_and_refusal_propagation() {
    let mut store = GeomStore::new();
    let s = store.add_surface(Surface::Sphere(
        Sphere::new(p(0, 0, 0), 5.into(), v(0, 0, 1), v(1, 0, 0)).unwrap(),
    ));
    let s1 = store.add_surface(Surface::Sphere(
        Sphere::new(p(6, 0, 0), 5.into(), v(0, 0, 1), v(1, 0, 0)).unwrap(),
    ));
    let s2 = store.add_surface(Surface::Sphere(
        Sphere::new(p(0, 6, 0), 5.into(), v(0, 0, 1), v(1, 0, 0)).unwrap(),
    ));
    let system = Arc::new(
        PolySystem::new(vec![
            store.implicit(s),
            store.implicit(s1),
            store.implicit(s2),
        ])
        .unwrap(),
    );
    let mut dec = Decider::default();
    let root = match dec.certify_root(
        &system,
        &RootBox::new(
            vec![2.into(), 2.into(), 2.into()],
            vec![4.into(), 4.into(), 3.into()],
        )
        .unwrap(),
    ) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        other => panic!("{other:?}"),
    };
    let q = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![s, s1, s2]),
        position: Position::Isolated(root),
    }));
    for surface in [s, s1, s2] {
        coincident(store.side_of(surface, &store.position(q), &mut dec));
    }
    let balls = store
        .enclose_point(q, &Rational::one().mul_pow2(-100), &mut dec)
        .unwrap();
    assert!(balls[0].contains(&3.into()) && balls[1].contains(&3.into()));
    assert!(balls[2].square().contains(&7.into()));
    for b in balls {
        assert!(b.rad() * Rational::from_i64(2) <= Rational::one().mul_pow2(-100));
    }
    let closed = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![s, s1]),
        position: Position::Closed(Point3::new(3.into(), 3.into(), Real::from_i64(7).sqrt())),
    }));
    coincident(store.same_point(q, closed, &mut dec));
    let cylinder = store.add_surface(Surface::Cylinder(
        Cylinder::new(p(0, 0, 0), v(0, 0, 1), 1.into(), v(1, 0, 0)).unwrap(),
    ));
    let root = Real::from(r("1/2")).sqrt();
    coincident(store.side_of(
        cylinder,
        &Point3::new(root.clone(), root.clone(), Real::zero()),
        &mut dec,
    ));
    let delta = Real::from(r("1e-30"));
    assert_eq!(
        store.side_of(
            cylinder,
            &Point3::new(&root + &delta, root.clone(), Real::zero()),
            &mut dec
        ),
        Outcome::Proven(Side::Front)
    );
    let undefined =
        Real::one() / (Real::from_i64(2).sqrt() * Real::from_i64(2).sqrt() - Real::from_i64(2));
    assert!(
        matches!(store.side_of(s,&Point3::new(undefined,Real::zero(),Real::zero()),&mut dec),Outcome::Refused(e) if e.code=="numbers.division-by-zero")
    );
}

#[test]
fn fact_closure_deadlines_and_numeric_short_circuit() {
    let mut store = GeomStore::new();
    let s = store.add_surface(shapes()[3].clone());
    let plane = store.add_surface(Surface::Plane(
        Plane::new(p(0, 0, 2), v(0, 0, 1), v(1, 0, 0)).unwrap(),
    ));
    let witness = explicit(&mut store, p(1, 2, 2));
    let c = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [s, plane],
        witness,
    }));
    let meet = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::CurveSurface(c, plane),
        position: Position::Closed(p(1, 2, 2).to_real()),
    }));
    let facts: Vec<_> = store.facts().iter().collect();
    assert_eq!(facts.len(), 4);
    assert_eq!(
        (facts[0].a, facts[0].b),
        (GeomKey::Curve(c), GeomKey::Surface(s))
    );
    assert_eq!(
        (facts[1].a, facts[1].b),
        (GeomKey::Curve(c), GeomKey::Surface(plane))
    );
    assert_eq!(
        (facts[2].a, facts[2].b),
        (GeomKey::Point(meet), GeomKey::Curve(c))
    );
    for fact in facts {
        assert_eq!(fact.origin, FactOrigin::Construction("shapes.definition"));
        assert_eq!(fact.at, None);
    }
    let mut dec = Decider::default();
    dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    assert_eq!(
        store.side(s, meet, &mut dec),
        Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ConstructionFact
        })
    );
    coincident(store.point_on_curve(c, meet, &mut dec));
    coincident(store.same_point(meet, meet, &mut dec));
    coincident(store.same_surface(s, s, &mut dec));
    assert_eq!(dec.stats().decisions, 4);
    assert_eq!(dec.stats().coincident_fact, 4);
    assert_eq!(dec.stats().work, 0);
    assert!(
        matches!(store.side(s,witness,&mut dec),Outcome::Refused(e) if e.code=="numbers.deadline")
    );
    let a = explicit(&mut store, p(0, 0, 0));
    let z = explicit(&mut store, p(1, 2, 3));
    let mut dec = Decider::default();
    assert_eq!(store.same_point(a, z, &mut dec), Outcome::Proven(Apart));
    assert_eq!(dec.stats().decisions, 1);
    let copy = store.add_surface(store.surface(s).clone());
    store.record_fact(Fact {
        kind: FactKind::Coincident,
        a: GeomKey::Surface(s),
        b: GeomKey::Surface(copy),
        at: None,
        origin: FactOrigin::Proof(ZeroProof::ExactArithmetic),
    });
    coincident(store.side(copy, meet, &mut dec));
}

#[test]
fn equality_is_geometric_but_cones_keep_the_nappe() {
    let mut store = GeomStore::new();
    let mut dec = Decider::default();
    let a = store.add_surface(Surface::Plane(
        Plane::new(p(0, 0, 1), v(0, 0, 1), v(1, 0, 0)).unwrap(),
    ));
    let b = store.add_surface(Surface::Plane(
        Plane::new(p(5, 7, 1), v(0, 0, -4), v(0, 1, 0)).unwrap(),
    ));
    coincident(store.same_surface(a, b, &mut dec));
    assert_eq!(dec.stats().decisions, 0);
    let at = explicit(&mut store, p(0, 0, 1));
    assert_eq!(
        store.normals_agree(a, b, at, &mut dec),
        Outcome::Proven(StrictSign::Negative)
    );
    let cone = store.add_surface(shapes()[2].clone());
    let opposite = store.add_surface(Surface::Cone(
        Cone::new(p(0, 0, 0), v(0, 0, -1), r("1/2"), v(1, 0, 0)).unwrap(),
    ));
    assert_eq!(
        store.same_surface(cone, opposite, &mut dec),
        Outcome::Proven(Apart)
    );
    let apex = explicit(&mut store, p(0, 0, 0));
    coincident(store.normals_agree(cone, a, apex, &mut dec));
}

#[test]
fn enclosures_derivatives_and_inverse_angles() {
    let surfaces = shapes();
    for precision in [1, 64, 1000] {
        let zero = b(0, precision);
        let half_pi = Ball::pi(precision) * Ball::from_rational(&r("1/2"), precision);
        contains(&surfaces[1].evaluate(&zero, &b(3, precision)), &p(5, 0, 3));
        contains(
            &surfaces[1].evaluate(&half_pi, &b(3, precision)),
            &p(0, 5, 3),
        );
        contains(&surfaces[3].evaluate(&zero, &half_pi), &p(0, 0, 3));
        contains(&surfaces[4].evaluate(&zero, &zero), &p(4, 0, 0));
        contains(
            &surfaces[2].evaluate(&zero, &b(-10, precision)),
            &p(-5, 0, -10),
        );
        let expected = [
            [p(1, 0, 0), p(0, 1, 0)],
            [p(0, 5, 0), p(0, 0, 1)],
            [p(0, 0, 0), RPoint3::new(r("1/2"), 0.into(), 1.into())],
            [p(0, 3, 0), p(0, 0, 3)],
            [p(0, 4, 0), p(0, 0, 1)],
        ];
        for (s, derivatives) in surfaces.iter().zip(expected) {
            let result = s.derivatives(&zero, &zero);
            contains(&result[0], &derivatives[0]);
            contains(&result[1], &derivatives[1]);
        }
    }
    let precision = 128;
    for (x, y) in [(3, 4), (-5, 0), (-3, -4), (3, -4)] {
        let parameters = surfaces[1]
            .parameters_of(&[b(x, precision), b(y, precision), b(7, precision)])
            .unwrap();
        assert!(parameters[1].contains(&7.into()));
        contains(
            &surfaces[1].evaluate(&parameters[0], &parameters[1]),
            &p(x, y, 7),
        );
    }
    let uv = surfaces[4]
        .parameters_of(&[b(0, precision), b(2, precision), b(0, precision)])
        .unwrap();
    contains(&surfaces[4].evaluate(&uv[0], &uv[1]), &p(0, 2, 0));
    let plane = Surface::Plane(Plane::new(p(1, 1, 1), v(0, 0, 1), v(0, 2, 0)).unwrap());
    contains(
        &plane.evaluate(&b(1, precision), &b(1, precision)),
        &p(0, 2, 1),
    );
    assert!(
        plane
            .parameters_of(&[b(1, precision), b(1, precision), b(1, precision)])
            .is_some()
    );
    let circle = Curve::Circle(Circle::new(p(0, 0, 0), v(0, 0, 1), 4.into(), v(0, 2, 0)).unwrap());
    contains(&circle.evaluate(&b(0, precision)).unwrap(), &p(0, 2, 0));
    contains(
        &circle
            .evaluate(&(Ball::pi(precision) * Ball::from_rational(&r("1/2"), precision)))
            .unwrap(),
        &p(-2, 0, 0),
    );
    let ellipse = Curve::Ellipse(
        Ellipse::new(p(0, 0, 0), v(0, 0, 1), v(2, 0, 0), 25.into(), 9.into()).unwrap(),
    );
    contains(&ellipse.evaluate(&b(0, precision)).unwrap(), &p(5, 0, 0));
}

#[test]
fn bounded_uncertainty_is_in_model_units() {
    let mut store = GeomStore::new();
    let cylinder = store.add_surface(Surface::Cylinder(
        Cylinder::new(p(0, 0, 0), v(0, 0, 2), 5.into(), v(1, 0, 0)).unwrap(),
    ));
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    });
    let eps = r("1e-30");
    let point = RPoint3::new(
        Rational::from_i64(3) * (Rational::one() + eps.clone()),
        Rational::from_i64(4) * (Rational::one() + eps.clone()),
        7.into(),
    );
    match store.side_of(cylinder, &point.to_real(), &mut dec) {
        Outcome::NotExact(n) => assert!(n.bound >= Rational::from_i64(5) * eps),
        Outcome::Proven(Side::Front) => {}
        other => panic!("{other:?}"),
    }
    // Exercise uncertainty with an algebraic zero hidden from the rational path.
    let half = Real::from(r("1/2")).sqrt();
    let unit = store.add_surface(Surface::Cylinder(
        Cylinder::new(p(0, 0, 0), v(0, 0, 1), 1.into(), v(1, 0, 0)).unwrap(),
    ));
    let point = Point3::new(&half + &Real::from(r("1e-30")), half, Real::zero());
    match store.side_of(unit, &point, &mut dec) {
        Outcome::NotExact(n) => assert!(n.bound >= r("7e-31")),
        Outcome::Proven(Side::Front) => {}
        other => panic!("{other:?}"),
    }
}

#[test]
fn zero_work_distance_bounds_cover_each_object() {
    let mut store = GeomStore::new();
    let surfaces = shapes().map(|s| store.add_surface(s));
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 0,
        max_work: 0,
    });
    // These points have rational, independently known distances to the carriers.
    let points = [p(0, 0, 3), p(7, 0, 2), p(0, 0, 0), p(5, 0, 0), p(6, 0, 0)];
    let distances = [2, 2, 0, 2, 2];
    for ((s, point), distance) in surfaces.into_iter().zip(points).zip(distances) {
        match store.side_of(s, &point.to_real(), &mut dec) {
            Outcome::NotExact(n) => assert!(n.bound >= Rational::from_i64(distance)),
            other => panic!("zero-work decision: {other:?}"),
        }
    }
    let curves = [
        Curve::Line(Line::new(p(0, 0, 0), v(2, 0, 0)).unwrap()),
        Curve::Circle(Circle::new(p(0, 0, 0), v(0, 0, 2), 25.into(), v(1, 0, 0)).unwrap()),
        Curve::Ellipse(
            Ellipse::new(p(0, 0, 0), v(0, 0, 2), v(2, 0, 0), 25.into(), 9.into()).unwrap(),
        ),
    ]
    .map(|c| store.add_curve(c));
    let point = explicit(&mut store, p(0, 0, 2));
    // Distances are 2, sqrt(29), sqrt(13), respectively.
    for (c, squared_distance) in curves.into_iter().zip([4, 29, 13]) {
        match store.point_on_curve(c, point, &mut dec) {
            Outcome::NotExact(n) => {
                assert!(&n.bound * &n.bound >= Rational::from_i64(squared_distance))
            }
            other => panic!("zero-work curve decision: {other:?}"),
        }
    }
    let other = explicit(&mut store, p(3, 4, 2));
    match store.same_point(point, other, &mut dec) {
        Outcome::NotExact(n) => assert!(n.bound >= 5.into()),
        other => panic!("{other:?}"),
    }
}

#[test]
fn stable_arenas_enclosure_short_circuit_and_determinism() {
    let mut arena = Arena::new();
    let a = arena.push("first");
    let z = arena.push("second");
    assert_eq!(arena.handles(), vec![a, z]);
    assert_eq!((a.index(), z.index()), (0, 1));
    assert_eq!(arena.get(a), &"first");
    let mut store = GeomStore::new();
    let s = store.add_surface(shapes()[3].clone());
    let copy = store.add_surface(shapes()[3].clone());
    assert_ne!(s, copy);
    let rational = explicit(&mut store, RPoint3::new(r("1/3"), 0.into(), 0.into()));
    assert_eq!(
        store.approx_point(rational)[0].to_bits(),
        r("1/3").to_f64_nearest().to_bits()
    );
    let undefined = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![s, copy]),
        position: Position::Closed(Point3::new(
            Real::one() / Real::zero(),
            Real::zero(),
            Real::zero(),
        )),
    }));
    let mut dec = Decider::default();
    assert!(
        matches!(store.enclose_point(undefined,&r("1e-20"),&mut dec),Err(e) if e.code=="numbers.division-by-zero")
    );
    assert_eq!(dec.stats().decisions, 1);
    let run = move |store: GeomStore| {
        let mut dec = Decider::default();
        let outcome = store.side(s, rational, &mut dec);
        (outcome, dec.take_stats())
    };
    let expected = run(store.clone());
    let workers: Vec<_> = (0..8)
        .map(|_| {
            let store = store.clone();
            std::thread::spawn(move || run(store))
        })
        .collect();
    for worker in workers {
        assert_eq!(worker.join().unwrap(), expected);
    }
}

#[test]
fn inverse_parameters_enclose_boxes_across_both_angle_cuts() {
    let cylinder = &shapes()[1];
    let precision = 128;
    let between = |a: i64, z: i64| b(a, precision).hull(&b(z, precision));
    for x in [-5, 5] {
        let input = [b(x, precision), between(-1, 2), b(7, precision)];
        let uv = cylinder.parameters_of(&input).unwrap();
        for y in [-1, 0, 1, 2] {
            let expected = Ball::atan2(&b(y, precision), &b(x, precision)).unwrap();
            let expected = if y < 0 {
                expected + Ball::pi(precision) * b(2, precision)
            } else {
                expected
            };
            assert!(uv[0].lo() <= expected.hi() && uv[0].hi() >= expected.lo());
        }
    }
    // Exact local subtraction must survive cancellation of an enormous translation.
    let origin = RPoint3::new(r("1e1000"), r("-1e1000"), r("1e1000"));
    let translated =
        Surface::Cylinder(Cylinder::new(origin.clone(), v(0, 0, 1), 5.into(), v(1, 0, 0)).unwrap());
    let point = origin.plus(&v(3, 4, 7));
    let coordinates = [point.x, point.y, point.z].map(|x| Ball::from_rational(&x, 4096));
    let uv = translated.parameters_of(&coordinates).unwrap();
    assert!(uv[1].contains(&7.into()));
}
