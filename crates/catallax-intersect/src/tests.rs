use super::*;
use catallax_numbers::{BudgetPolicy, Real, StrictSign};
use catallax_shapes::{
    Circle, Curve, Cylinder, GeomKey, Line, Plane, PointDef, RPoint3, Sphere, Surface,
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
        Plane::new(origin, normal, reference).unwrap(),
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
        Cylinder::new(origin, axis, r(radius), reference).unwrap(),
    ))
}
fn sphere(g: &mut GeomStore, center: RPoint3, radius: i64) -> SurfaceId {
    g.add_surface(Surface::Sphere(
        Sphere::new(center, r(radius), v(0, 0, 1), v(1, 0, 0)).unwrap(),
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
fn rational_configuration_fields_and_no_decisions() {
    let mut g = GeomStore::new();
    let mut d = dec();
    let a = plane(&mut g, p(0, 0, 1), v(0, 0, 2), v(1, 0, 0));
    let b = plane(&mut g, p(3, 0, 0), v(1, 0, 0), v(0, 1, 0));
    let c = curve(&components(&mut g, a, b, &mut d)[0]);
    assert_eq!(d.stats().decisions, 0);
    assert!(matches!(g.curve(c),Curve::Line(l) if l.origin==p(3,0,1)&&l.dir==v(0,2,0)));
    assert!(g.facts().holds(
        catallax_numbers::FactKind::Incident,
        GeomKey::Curve(c),
        GeomKey::Surface(a)
    ));
    let z = cyl(&mut g, p(0, 0, 0), v(0, 0, 1), 5, v(1, 0, 0));
    let tilted = plane(&mut g, p(0, 0, 0), v(0, 3, 4), v(1, 0, 0));
    let c = curve(&components(&mut g, tilted, z, &mut d)[0]);
    assert!(
        matches!(g.curve(c),Curve::Ellipse(e) if e.major_sq==Rational::from_ratio(625,16).unwrap()&&e.minor_sq==r(25))
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
            matches!(curve_form(&g,c,&mut d).unwrap(),CurveForm::Line{dir,..} if dir==expected_dir)
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
            matches!(curve_form(&g,c,&mut d).unwrap(),CurveForm::Circle{radius_sq,..}if radius_sq==r(25))
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
    let c = g.add_curve(Curve::Line(Line::new(p(-10, 0, 0), v(1, 0, 0)).unwrap()));
    let hits = hit_positions(&mut g, c, s, &mut d);
    assert_eq!(hits.len(), 2);
    for ((point, contact), x) in hits.iter().zip([-5, 5]) {
        point_equal(point, &p(x, 0, 0).to_real(), &mut d);
        assert_eq!(*contact, Contact::Crossing);
    }
    let c = g.add_curve(Curve::Line(Line::new(p(-10, 5, 0), v(1, 0, 0)).unwrap()));
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
        Circle::new(p(0, 0, 0), v(0, 0, 1), r(25), v(1, 0, 0)).unwrap(),
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
        Circle::new(p(0, 0, 0), v(0, 0, 1), r(25), v(1, 0, 0)).unwrap(),
    ));
    let a = g.add_point(PointDef::Explicit(p(5, 0, 0)));
    let b = g.add_point(PointDef::Explicit(p(0, 5, 0)));
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
            assert!(g.implicit(c).eval_ball(&p).contains_zero());
            assert!(g.implicit(s).eval_ball(&p).contains_zero());
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
