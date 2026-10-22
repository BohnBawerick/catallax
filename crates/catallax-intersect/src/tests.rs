use super::*;
use catallax_numbers::{BudgetPolicy, Real, StrictSign};
use catallax_shapes::{
    Circle, Curve, Cylinder, GeomKey, Line, Plane, PointDef, RPoint3, RVector3, Sphere, Surface,
};
fn r(x: i64) -> Rational {
    Rational::from_i64(x)
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
fn dec() -> Decider {
    Decider::new(BudgetPolicy::default())
}
fn plane(g: &mut GeomStore, origin: RPoint3, normal: RVector3, reference: RVector3) -> SurfaceId {
    g.add_surface(Surface::Plane(
        Plane::new(
            (origin).into(),
            (normal).into(),
            (reference).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ))
}
fn cyl(
    g: &mut GeomStore,
    origin: RPoint3,
    axis: RVector3,
    radius: i64,
    reference: RVector3,
) -> SurfaceId {
    g.add_surface(Surface::Cylinder(
        Cylinder::new(
            (origin).into(),
            (axis).into(),
            (r(radius)).into(),
            (reference).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ))
}
fn sphere(g: &mut GeomStore, center: RPoint3, radius: i64) -> SurfaceId {
    g.add_surface(Surface::Sphere(
        Sphere::new(
            (center).into(),
            (r(radius)).into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ))
}
fn components(g: &mut GeomStore, a: SurfaceId, b: SurfaceId, d: &mut Decider) -> Vec<Component> {
    match intersect_surfaces(g, a, b, d).unwrap().value {
        SurfacePair::Meets(c) => c,
        other => panic!("expected components: {other:?}"),
    }
}
fn curve(c: &Component) -> CurveId {
    match c {
        Component::Curve { curve, .. } => *curve,
        _ => panic!("expected curve"),
    }
}
fn point_equal(a: &Point3, b: &Point3, d: &mut Decider) {
    assert!(d.sign(&a.minus(b).norm_squared()).is_coincident());
}
fn hit_positions(
    g: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    d: &mut Decider,
) -> Vec<(Point3, Contact)> {
    let answer = intersect_curve_surface(g, c, s, d).unwrap();
    assert!(answer.is_exact());
    match answer.value {
        CurveSurface::Meets(h) => h
            .into_iter()
            .map(|h| (g.position(h.point), h.contact))
            .collect(),
        other => panic!("expected hits: {other:?}"),
    }
}

#[test]
fn constructed_tangent_locus_skips_gradient_zero_proofs() {
    let mut g = GeomStore::new();
    let a = plane(&mut g, p(0, 5, 0), v(0, 1, 0), v(1, 0, 0));
    let b = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    let parts = components(&mut g, a, b, &mut dec());
    let locus = curve(&parts[0]);
    assert!(matches!(
        parts[0],
        Component::Curve {
            contact: Contact::Tangent,
            ..
        }
    ));
    let witness = g.add_point(PointDef::Meet(catallax_shapes::Meet {
        of: catallax_shapes::MeetOf::CurveSurface(locus, a),
        position: catallax_shapes::Position::Closed(p(0, 5, 1).to_real()),
    }));
    let definition = g.add_curve(Curve::Intersection(catallax_shapes::IntersectionCurve {
        surfaces: [a, b],
        witness,
    }));
    let mut d = dec();
    let vectors = tangent_at(&g, definition, witness, &mut d).unwrap();
    assert_eq!(d.stats().coincident_fact, 1);
    assert_eq!(d.stats().coincident_exact, 0);
    for tangent in vectors {
        point_equal(
            &Point3::origin().plus(&tangent),
            &p(0, 0, 1).to_real(),
            &mut dec(),
        );
    }
}
#[test]
fn rational_configuration_fields_and_no_decisions() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let a = plane(&mut g, p(0, 0, 1), v(0, 0, 2), v(1, 0, 0));
    let b = plane(&mut g, p(3, 0, 0), v(1, 0, 0), v(0, 1, 0));
    let c = curve(&components(&mut g, a, b, &mut d)[0]);
    assert_eq!(d.stats().decisions, 0);
    assert!(
        matches!(g.curve(c),Curve::Line(l) if l.origin.as_rational().unwrap()==p(3,0,1)&&l.dir.as_rational().unwrap()==v(0,2,0))
    );
    assert!(g.facts().holds(
        catallax_numbers::FactKind::Incident,
        GeomKey::Curve(c),
        GeomKey::Surface(a)
    ));
    let z = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    let tilted = plane(&mut g, p(0, 0, 0), v(0, 3, 4), v(1, 0, 0));
    let c = curve(&components(&mut g, tilted, z, &mut d)[0]);
    assert!(
        matches!(g.curve(c),Curve::Ellipse(e) if e.major_sq.rational_ref()==Some(&Rational::from_ratio(625,16).unwrap())&&e.minor_sq.rational_ref()==Some(&r(25)))
    );
    let vertex = Point3::new(
        Real::zero(),
        Real::from_i64(-5),
        Real::from_rational(Rational::from_ratio(15, 4).unwrap()),
    );
    assert!(g.side_of(z, &vertex, &mut d).is_coincident());
    assert!(g.side_of(tilted, &vertex, &mut d).is_coincident());
}
#[test]
fn surface_definitions_and_argument_order_survive_repeated_queries() {
    // Every fresh store reuses the same handles. Different radii, centres and
    // argument orders must still give their own exact oriented circle.
    for offset in (0..40).chain([0, 17, 39]) {
        let mut g = GeomStore::new();
        let mut d = dec();
        let radius = 5 + offset % 3;
        let a = sphere(&mut g, p(offset, 0, 0), radius);
        let b = sphere(&mut g, p(offset + 6, 0, 0), radius);
        for (left, right, direction) in [(a, b, 6), (b, a, -6), (a, b, 6)] {
            let c = curve(&components(&mut g, left, right, &mut d)[0]);
            let Curve::Circle(circle) = g.curve(c) else {
                panic!("sphere section must be a rational circle");
            };
            assert_eq!(circle.center.as_rational().unwrap(), p(offset + 3, 0, 0));
            assert_eq!(
                circle.radius_sq.rational_ref(),
                Some(&r(radius * radius - 9))
            );
            assert_eq!(circle.normal.as_rational().unwrap(), v(direction, 0, 0));
            assert_eq!(circle.ref_dir.as_rational().unwrap(), v(0, 0, direction));
        }
    }
}

#[test]
fn irrational_lines_and_circles_keep_exact_definitions() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let z = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    let x = plane(&mut g, p(1, 0, 0), v(1, 0, 0), v(0, 1, 0));
    let cs = components(&mut g, x, z, &mut d);
    assert_eq!(g.point_count(), 2);
    for (c, expected_dir) in cs.iter().zip([v(0, 0, -1), v(0, 0, 1)]) {
        let c = curve(c);
        let Curve::Intersection(i) = g.curve(c) else {
            panic!()
        };
        let point = g.position(i.witness);
        assert!(g.side_of(x, &point, &mut d).is_coincident());
        assert!(g.side_of(z, &point, &mut d).is_coincident());
        assert!(
            matches!(curve_form(&g,c,&mut d).unwrap(),CurveForm::Line{dir,..} if dir.as_rational().unwrap()==expected_dir)
        );
    }
    assert!(same_curve(&g, curve(&cs[0]), curve(&cs[1]), &mut d).is_proven());
    let s = sphere(&mut g, p(0, 0, 1), 6);
    let cs = components(&mut g, z, s, &mut d);
    for c in &cs {
        let c = curve(c);
        let Curve::Intersection(i) = g.curve(c) else {
            panic!()
        };
        let point = g.position(i.witness);
        assert!(g.side_of(z, &point, &mut d).is_coincident());
        assert!(g.side_of(s, &point, &mut d).is_coincident());
        assert!(
            matches!(curve_form(&g,c,&mut d).unwrap(),CurveForm::Circle{radius_sq,..}if radius_sq.as_rational()==Some(r(25)))
        );
    }
}
#[test]
fn cylinder_nodes_and_steinmetz_refusal_are_transactional() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let a = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 2, v(1, 0, 0));
    let b = cyl(&mut g, p(0, 1, 0), v(1, 0, 0), 1, v(0, 0, 1));
    let cs = components(&mut g, a, b, &mut d);
    let singular: Vec<_> = cs
        .iter()
        .map(|c| match c {
            Component::Curve { singular, .. } => singular[0],
            _ => panic!(),
        })
        .collect();
    assert_eq!(singular[0], singular[1]);
    point_equal(&g.position(singular[0]), &p(0, 2, 0).to_real(), &mut d);
    let a = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 1, v(1, 0, 0));
    let b = cyl(&mut g, p(0, 0, 0), v(1, 1, 0), 1, v(0, 0, 1));
    let count = (g.curve_count(), g.point_count());
    assert_eq!(
        intersect_surfaces(&mut g, a, b, &mut d).unwrap_err().code,
        "intersect.unsupported-pair"
    );
    assert_eq!(count, (g.curve_count(), g.point_count()));
}
#[test]
fn inside_tangency_near_misses_follow_exact_distance() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let c = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 3, v(1, 0, 0));
    let eps = Rational::parse("1e-30").unwrap();
    let s = sphere(&mut g, RPoint3::new(r(1) + eps.clone(), r(0), r(0)), 2);
    assert_eq!(components(&mut g, c, s, &mut d).len(), 1);
    let s = sphere(&mut g, RPoint3::new(r(1) - eps, r(0), r(0)), 2);
    assert_eq!(
        intersect_surfaces(&mut g, c, s, &mut d).unwrap().value,
        SurfacePair::Disjoint
    );
}
#[test]
fn line_meetings_sort_and_classify_multiplicity() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let s = sphere(&mut g, p(0, 0, 0), 5);
    let c = g.add_curve(Curve::Line(
        Line::new(
            (p(-10, 0, 0)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let hits = hit_positions(&mut g, c, s, &mut d);
    assert_eq!(hits.len(), 2);
    for ((point, contact), x) in hits.iter().zip([-5, 5]) {
        point_equal(point, &p(x, 0, 0).to_real(), &mut d);
        assert_eq!(*contact, Contact::Crossing);
    }
    let c = g.add_curve(Curve::Line(
        Line::new(
            (p(-10, 5, 0)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let hits = hit_positions(&mut g, c, s, &mut d);
    assert_eq!(hits.len(), 1);
    point_equal(&hits[0].0, &p(0, 5, 0).to_real(), &mut d);
    assert_eq!(hits[0].1, Contact::Tangent);
}
#[test]
fn circle_quadratic_meetings_and_containment() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let c = g.add_curve(Curve::Circle(
        Circle::new(
            (p(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            (r(25)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let s = cyl(&mut g, p(0, 0, 0), v(1, 0, 0), 3, v(0, 0, 1));
    let hits = hit_positions(&mut g, c, s, &mut d);
    assert_eq!(hits.len(), 4);
    for ((point, contact), (x, y)) in hits.iter().zip([(4, 3), (-4, 3), (-4, -3), (4, -3)]) {
        point_equal(point, &p(x, y, 0).to_real(), &mut d);
        assert_eq!(*contact, Contact::Crossing);
    }
    let s = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    assert!(matches!(
        intersect_curve_surface(&mut g, c, s, &mut d).unwrap().value,
        CurveSurface::Contained
    ));
    let s = cyl(&mut g, p(3, 0, 0), v(0, 1, 1), 2, v(1, 0, 0));
    let hits = hit_positions(&mut g, c, s, &mut d);
    assert_eq!(hits.len(), 1);
    point_equal(&hits[0].0, &p(5, 0, 0).to_real(), &mut d);
    assert_eq!(hits[0].1, Contact::Tangent);
}
fn viviani(g: &mut GeomStore, d: &mut Decider) -> (SurfaceId, SurfaceId, Vec<Component>) {
    let c = cyl(g, p(0, 0, 0), v(0, 0, 1), 1, v(1, 0, 0));
    let s = sphere(g, p(1, 0, 0), 2);
    let cs = components(g, c, s, d);
    (c, s, cs)
}
#[test]
fn viviani_branch_hits_order_node_contact_and_tangents() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let (_, _, cs) = viviani(&mut g, &mut d);
    let upper = curve(&cs[0]);
    let z = plane(&mut g, p(0, 0, 1), v(0, 0, 1), v(1, 0, 0));
    let hits = hit_positions(&mut g, upper, z, &mut d);
    assert_eq!(hits.len(), 2);
    assert!(matches!(
        d.sign(&hits[0].0.y),
        Outcome::Proven(StrictSign::Negative)
    ));
    assert!(matches!(
        d.sign(&hits[1].0.y),
        Outcome::Proven(StrictSign::Positive)
    ));
    let z = plane(&mut g, p(0, 0, 0), v(0, 0, 1), v(1, 0, 0));
    let hits = hit_positions(&mut g, upper, z, &mut d);
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].1, Contact::Tangent);
    point_equal(&hits[0].0, &p(-1, 0, 0).to_real(), &mut d);
    let node = match &cs[0] {
        Component::Curve { singular, .. } => singular[0],
        _ => panic!(),
    };
    let ts = tangent_at(&g, upper, node, &mut d).unwrap();
    for (t, expected) in ts.iter().zip([v(0, 1, -1), v(0, 1, 1)]) {
        assert!(
            d.sign(&t.cross(&expected.to_real()).norm_squared())
                .is_coincident()
        );
        assert!(matches!(
            d.sign(&t.dot(&expected.to_real())),
            Outcome::Proven(StrictSign::Positive)
        ));
    }
}
#[test]
fn circle_trace_taylor_coefficients_and_arc_point() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let c = g.add_curve(Curve::Circle(
        Circle::new(
            (p(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            (r(25)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let a = g.add_point(PointDef::Explicit((p(5, 0, 0)).into()));
    let b = g.add_point(PointDef::Explicit((p(0, 5, 0)).into()));
    let trace = trace(&g, c, Some(a), Some(b), &mut d).unwrap().value;
    assert_eq!(trace.pieces().len(), 1);
    let t = Ball::from_rational(&r(0), 128);
    let j = trace.pieces()[0]
        .taylor(&t, 4, &Decider::default())
        .unwrap();
    assert!(j[0][0].contains(&r(5)));
    assert!(j[1][1].contains(&r(10)));
    assert!(j[2][0].contains(&r(-10)));
    assert!(j[3][1].contains(&r(-10)));
    assert!(j[4][0].contains(&r(10)));
    let mut expired = Decider::default();
    expired.set_deadline(Some(std::time::Instant::now()));
    assert_eq!(
        trace.pieces()[0].taylor(&t, 4, &expired).unwrap_err().code,
        "numbers.deadline",
    );
    let mid = arc_point(&g, c, Some(a), Some(b), &mut d).unwrap().value;
    assert!(
        d.sign(&(mid.x.square() + mid.y.square() - Real::from_i64(25)))
            .is_coincident()
    );
    assert!(matches!(
        d.sign(&mid.x),
        Outcome::Proven(StrictSign::Positive)
    ));
    assert!(matches!(
        d.sign(&mid.y),
        Outcome::Proven(StrictSign::Positive)
    ));
}
#[test]
fn nodal_trace_is_regular_at_both_ends() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let (c, s, cs) = viviani(&mut g, &mut d);
    let curve = curve(&cs[0]);
    let node = match &cs[0] {
        Component::Curve { singular, .. } => singular[0],
        _ => panic!(),
    };
    let trace = trace(&g, curve, Some(node), Some(node), &mut d)
        .unwrap()
        .value;
    for piece in trace.pieces() {
        for t in [r(0), Rational::from_ratio(1, 2).unwrap(), r(1)] {
            let t = Ball::from_rational(&t, 128);
            let p = piece.evaluate(&t, &Decider::default()).unwrap();
            let j = piece.taylor(&t, 4, &Decider::default()).unwrap();
            for i in 0..3 {
                assert!(j[0][i].lo() <= p[i].hi() && p[i].lo() <= j[0][i].hi());
            }
            assert!(
                g.implicit(c)
                    .eval_ball(&p, &catallax_numbers::Decider::default())
                    .unwrap()
                    .contains_zero()
            );
            assert!(
                g.implicit(s)
                    .eval_ball(&p, &catallax_numbers::Decider::default())
                    .unwrap()
                    .contains_zero()
            );
            let derivative = piece.derivative(&t, &Decider::default()).unwrap();
            assert!(derivative.iter().any(|v| v.sign().is_some()));
        }
    }
    let middle = arc_point(&g, curve, None, None, &mut d).unwrap().value;
    assert!(g.side_of(c, &middle, &mut d).is_coincident());
    assert!(g.side_of(s, &middle, &mut d).is_coincident());
    assert!(!matches!(
        d.sign(&middle.z),
        Outcome::Proven(StrictSign::Negative)
    ));
}
#[test]
fn fact_answers_ignore_expired_deadline() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let (c, _, cs) = viviani(&mut g, &mut d);
    let curve = curve(&cs[0]);
    let node = match &cs[0] {
        Component::Curve { singular, .. } => singular[0],
        _ => panic!(),
    };
    d.set_deadline(Some(std::time::Instant::now()));
    assert!(matches!(
        intersect_curve_surface(&mut g, curve, c, &mut d)
            .unwrap()
            .value,
        CurveSurface::Contained
    ));
    assert!(on_curve(&g, curve, node, &mut d).is_coincident());
    assert!(same_curve(&g, curve, curve, &mut d).is_coincident());
}

#[test]
fn quartic_charts_cover_every_regular_branch_family() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let mut pairs = vec![];
    for offset in [0, 1, 2] {
        let a = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 2, v(1, 0, 0));
        let b = cyl(&mut g, p(0, offset, 0), v(1, 0, 0), 1, v(0, 0, 1));
        pairs.push((a, b));
    }
    for (offset, radius) in [(2, 2), (1, 2), (1, 3)] {
        let a = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 1, v(1, 0, 0));
        let b = sphere(&mut g, p(offset, 0, 0), radius);
        pairs.push((a, b));
    }
    for (a, b) in pairs {
        let cs = components(&mut g, a, b, &mut d);
        for component in cs {
            let c = curve(&component);
            let charts = crate::query::charts(&g, c, &mut d).unwrap();
            for ch in &charts {
                let t = Real::from_rational(Rational::from_ratio(1, 3).unwrap());
                let point = ch.position(&t);
                assert!(g.side_of(a, &point, &mut d).is_coincident());
                assert!(g.side_of(b, &point, &mut d).is_coincident());
                let j = ch.exact_jets(&t, 1);
                let derivative = Vector3::new(j[1][0].clone(), j[1][1].clone(), j[1][2].clone());
                let tangent = g
                    .surface(a)
                    .gradient(&point)
                    .cross(&g.surface(b).gradient(&point));
                assert!(matches!(
                    d.sign(&derivative.dot(&tangent)),
                    Outcome::Proven(StrictSign::Positive)
                ));
            }
            for i in 0..charts.len() {
                point_equal(
                    &charts[i].position(&Real::one()),
                    &charts[(i + 1) % charts.len()].position(&Real::zero()),
                    &mut d,
                );
            }
        }
    }
}

#[test]
fn intersection_line_order_uses_one_decision() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let z = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    let x = plane(&mut g, p(1, 0, 0), v(1, 0, 0), v(0, 1, 0));
    let cs = components(&mut g, x, z, &mut d);
    let c = curve(&cs[0]);
    let y = -Real::from_i64(24).sqrt();
    let a = g.add_point(PointDef::Meet(catallax_shapes::Meet {
        of: catallax_shapes::MeetOf::Surfaces(vec![x, z]),
        position: catallax_shapes::Position::Closed(Point3::new(
            Real::one(),
            y.clone(),
            Real::zero(),
        )),
    }));
    let b = g.add_point(PointDef::Meet(catallax_shapes::Meet {
        of: catallax_shapes::MeetOf::Surfaces(vec![x, z]),
        position: catallax_shapes::Position::Closed(Point3::new(
            Real::one(),
            y,
            Real::from_i64(-1),
        )),
    }));
    let before = d.stats().decisions;
    assert!(matches!(
        order_on_curve(&g, c, a, b, &mut d),
        Outcome::Proven(StrictOrder::Less)
    ));
    assert_eq!(d.stats().decisions - before, 1);
}

#[test]
fn irrational_conic_dimensions_preserve_exact_points_and_taylor_bounds() {
    let radius = Real::from_i64(3) + Real::from_i64(3).sqrt() / Real::from_i64(2);
    let major = Real::from_i64(2) + Real::from_i64(2).sqrt();
    let forms = [
        (
            CurveForm::Circle {
                center: p(0, 0, 0).to_real(),
                normal: v(0, 0, 2).into(),
                radius_sq: radius.square(),
                ref_dir: v(3, 0, 0).into(),
            },
            radius.clone(),
            radius,
        ),
        (
            CurveForm::Ellipse {
                center: p(0, 0, 0).to_real(),
                normal: v(0, 0, 2).into(),
                major_dir: v(3, 0, 0).into(),
                major_sq: major.square(),
                minor_sq: Real::from_i64(2),
            },
            major,
            Real::from_i64(2).sqrt(),
        ),
    ];
    let mut dec = dec();
    for (form, a, b) in forms {
        let charts = crate::chart::form_charts(&form);
        let parameter = Real::from(Rational::from_ratio(1, 3).unwrap());
        let point = charts[0].position(&parameter);
        assert!(
            dec.sign(
                &(point.x.square() / a.square() + point.y.square() / b.square() - Real::one())
            )
            .is_coincident()
        );
        let coords = charts[0]
            .restricted(Real::zero(), Real::one(), &mut dec)
            .unwrap();
        let piece = TracePiece::new(coords);
        let coefficients = piece
            .taylor(&Ball::from_rational(&r(0), 128), 4, &dec)
            .unwrap();
        for (order, axis, expected) in [
            (0, 0, a.clone()),
            (1, 1, &Real::from_i64(2) * &b),
            (2, 0, &Real::from_i64(-2) * &a),
            (3, 1, &Real::from_i64(-2) * &b),
            (4, 0, &Real::from_i64(2) * &a),
        ] {
            let bound = &coefficients[order][axis];
            let expected = dec.enclose_at(&expected, 256).unwrap();
            assert!(bound.lo() <= expected.lo() && bound.hi() >= expected.hi());
            assert!(bound.rad() < Rational::one().mul_pow2(-110));
        }
    }
}

#[test]
fn cyclic_system_keeps_both_roots_under_every_generator_order() {
    use crate::solve::{Solutions, intersections};
    use catallax_numbers::MPoly;

    // xy = yz = zx = 1 forces x = y = z = +/-1. This cyclic ideal
    // needs critical pairs whose necessity depends on previously resolved pairs.
    let x = MPoly::var(3, 0);
    let y = MPoly::var(3, 1);
    let z = MPoly::var(3, 2);
    let one = MPoly::constant(3, r(1));
    let equations = [&(&x * &y) - &one, &(&y * &z) - &one, &(&z * &x) - &one]
        .map(catallax_numbers::RealPoly::from);
    for order in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let mut dec = dec();
        let Solutions::Points(points) = intersections(
            &[equations[order[0]].clone(), equations[order[1]].clone()],
            &equations[order[2]],
            &mut dec,
        )
        .unwrap() else {
            panic!("the cyclic system has exactly two isolated real solutions");
        };
        let mut counts = [0, 0];
        for point in points {
            let (index, expected) = match dec.sign(&point.x) {
                Outcome::Proven(StrictSign::Negative) => (0, Real::from_i64(-1)),
                Outcome::Proven(StrictSign::Positive) => (1, Real::one()),
                other => panic!("nonzero cyclic coordinate: {other:?}"),
            };
            counts[index] += 1;
            for coordinate in [point.x, point.y, point.z] {
                assert!(matches!(
                    dec.sign(&(&coordinate - &expected)),
                    Outcome::ProvenCoincident(_)
                ));
            }
        }
        assert_eq!(counts, [1, 1]);
    }
}
