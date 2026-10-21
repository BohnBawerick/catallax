//! Geometry definitions, fact closure and carried marks at the topology boundary.
mod common;
use common::*;

blind_test!(
    intersection_edges_and_meeting_vertices_need_no_coordinates,
    {
        for with_vertex in [false, true] {
            let mut m = Model::new();
            let sphere = m.geom_mut().add_surface(sphere_surface(r(2)));
            let equator = m
                .geom_mut()
                .add_surface(plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)));
            let meridian = m
                .geom_mut()
                .add_surface(plane(p(0, 0, 0), v(0, 1, 0), v(1, 0, 0)));
            // x²+y²+z²=4, z=0, y=0; select the positive-x solution.
            let point = m.geom_mut().add_point(PointDef::Meet(Meet {
                of: MeetOf::Surfaces(vec![sphere, equator, meridian]),
                position: Position::Closed(p(2, 0, 0).to_real()),
            }));
            let curve = m
                .geom_mut()
                .add_curve(Curve::Intersection(IntersectionCurve {
                    surfaces: [sphere, equator],
                    witness: point,
                }));
            assert_eq!(m.geom().curve(curve).is_closed(), None);
            let end = if with_vertex {
                Some(m.add_vertex(Vertex { point }))
            } else {
                None
            };
            record(
                &mut m,
                FactKind::Incident,
                GeomKey::Point(point),
                GeomKey::Curve(curve),
            );
            // grad(sphere) cross grad(z) is clockwise at positive x, so the edge
            // reverses it to use the same equator orientation as numbered case 18.
            let edge = m.add_edge(Edge {
                curve,
                sense: Sense::Reversed,
                start: end,
                end,
            });
            let north = cycle(&mut m, &[(edge, Sense::Same)]);
            let south = cycle(&mut m, &[(edge, Sense::Reversed)]);
            let north = face_on(&mut m, sphere, Sense::Same, vec![north], false);
            let south = face_on(&mut m, sphere, Sense::Same, vec![south], false);
            let b = body_from_faces(&mut m, vec![north, south]);
            sound(&m, b);
            assert_eq!(m.vertices_of_body(b), end.into_iter().collect::<Vec<_>>());
            let c = m.counts(b);
            assert_eq!(c.edges_by_kind, [(CurveKind::Intersection, 1)].into());
            assert_eq!(c.ring_edges, usize::from(!with_vertex));
            assert_eq!(euler(&m, shell_of(&m, b)), 2);
        }
    }
);

blind_test!(check_uses_factbook_closure_not_only_direct_records, {
    let mut m = Model::new();
    let s0 = m.geom_mut().add_surface(sphere_surface(r(2)));
    let s1 = m.geom_mut().add_surface(sphere_surface(r(2)));
    let point0 = m.geom_mut().add_point(PointDef::Explicit(p(2, 0, 0)));
    let point1 = m.geom_mut().add_point(PointDef::Explicit(p(2, 0, 0)));
    let curve = m.geom_mut().add_curve(circle(p(0, 0, 0), r(2)));
    record(
        &mut m,
        FactKind::Coincident,
        GeomKey::Surface(s0),
        GeomKey::Surface(s1),
    );
    record(
        &mut m,
        FactKind::Coincident,
        GeomKey::Point(point0),
        GeomKey::Point(point1),
    );
    record(
        &mut m,
        FactKind::Incident,
        GeomKey::Point(point1),
        GeomKey::Curve(curve),
    );
    record(
        &mut m,
        FactKind::Incident,
        GeomKey::Curve(curve),
        GeomKey::Surface(s0),
    );
    let vertex = m.add_vertex(Vertex { point: point0 });
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(vertex),
        end: Some(vertex),
    });
    let a = cycle(&mut m, &[(edge, Sense::Same)]);
    let b = cycle(&mut m, &[(edge, Sense::Reversed)]);
    let a = face_on(&mut m, s1, Sense::Same, vec![a], false);
    let b = face_on(&mut m, s1, Sense::Same, vec![b], false);
    let body = body_from_faces(&mut m, vec![a, b]);
    sound(&m, body);
    assert!(m.geom().facts().holds(
        FactKind::Incident,
        GeomKey::Point(point0),
        GeomKey::Curve(curve)
    ));
    assert!(m.geom().facts().holds(
        FactKind::Incident,
        GeomKey::Curve(curve),
        GeomKey::Surface(s1)
    ));
    assert_eq!(m.geom().facts().len(), 4);
});

blind_test!(body_marks_survive_sharing_cloning_and_queries, {
    let mut m = Model::new();
    let b = assemble(
        &mut m,
        &expected(&Primitive::Sphere(r(2)), &Pose::identity()),
    );
    let mut marks = MarkLog::new();
    marks.push(catallax_numbers::Mark::new(
        "fixture.bound",
        catallax_numbers::BoundKind::Distance,
        q(1, 7),
        "carried operand bound",
    ));
    let body = Body {
        kind: BodyKind::Solid,
        lumps: m.body(b).lumps.clone(),
        marks: marks.clone(),
    };
    let marked = m.add_body(body.clone());
    sound(&m, marked);
    assert_eq!(m.body(marked), &body);
    assert_eq!(m.counts(marked), m.counts(b));
    assert_eq!(m.body(marked).marks, marks);
    assert!(m.body(b).marks.is_empty());
    let mut clone = m.clone();
    empty(&mut clone);
    assert_eq!(clone.body(marked), &body);
    assert_eq!(m.body(marked), &body);
});

blind_test!(third_distinct_coedge_use_is_an_edge_use_fault, {
    let mut m = Model::new();
    let edge = edge(&mut m, circle(p(0, 0, 0), r(2)), None, None);
    let surface = m.geom_mut().add_surface(sphere_surface(r(2)));
    let mut faces = vec![];
    for sense in [Sense::Same, Sense::Reversed, Sense::Same] {
        let l = cycle(&mut m, &[(edge, sense)]);
        faces.push(face_on(&mut m, surface, Sense::Same, vec![l], true));
    }
    let b = body_from_faces(&mut m, faces);
    faults(&m, b, &["topology.edge-use"]);
});
