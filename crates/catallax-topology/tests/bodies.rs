//! Hand-built sound bodies and ordered topology queries.
mod common;
use common::*;

blind_test!(case_13_backlinks_before_ownership_and_first_face_wins, {
    let mut m = Model::new();
    let e = edge(&mut m, circle(p(0, 0, 0), r(2)), None, None);
    let c0 = m.add_coedge(Coedge {
        edge: e,
        sense: Sense::Same,
    });
    let c1 = m.add_coedge(Coedge {
        edge: e,
        sense: Sense::Reversed,
    });
    assert_eq!(m.coedges_of_edge(e), vec![c0, c1]);
    assert_eq!(m.owner_of_coedge(c0), None);
    let l = m.add_loop(Loop::Coedges(vec![c0]));
    assert_eq!(m.owner_of_coedge(c0), None);
    let s = m.geom_mut().add_surface(sphere_surface(r(2)));
    let f0 = face_on(&mut m, s, Sense::Same, vec![l], true);
    let f1 = face_on(&mut m, s, Sense::Same, vec![l], true);
    assert_eq!(m.owner_of_coedge(c0), Some((l, f0)));
    assert_eq!(m.owner_of_coedge(c1), None);
    let orphan = m.add_coedge(Coedge {
        edge: e,
        sense: Sense::Same,
    });
    assert_eq!(m.coedges_of_edge(e), vec![c0, c1, orphan]);
    assert_eq!(m.owner_of_coedge(orphan), None);
    let b = body_from_faces(&mut m, vec![f1, f0]);
    faults(&m, b, &["topology.shared-entity"]);
});

blind_test!(case_14_walk_order_deduplicates_at_first_sight, {
    let mut m = Model::new();
    let a = assemble(&mut m, &expected(&primitives()[0], &Pose::identity()));
    let pose = Pose::new(
        Rotation::identity(),
        v(10, 20, 30).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let b = assemble(&mut m, &expected(&primitives()[0], &pose));
    let mut fa = m.shell(shell_of(&m, a)).faces.clone();
    let fb = m.shell(shell_of(&m, b)).faces.clone();
    fa.reverse();
    fa.push(fa[0]);
    let s0 = m.add_shell(Shell { faces: fa.clone() });
    let s1 = m.add_shell(Shell { faces: fb.clone() });
    let l0 = m.add_lump(Lump {
        shells: vec![s1, s0],
    });
    let l1 = m.add_lump(Lump { shells: vec![s0] });
    let body = body_from_lumps(&mut m, vec![l0, l1, l0]);
    let (fs, es, vs) = walk(&m, body);
    assert_eq!(fs, [fb, fa[..6].to_vec()].concat());
    assert_eq!(m.faces_of_body(body), fs);
    assert_eq!(m.edges_of_body(body), es);
    assert_eq!(m.vertices_of_body(body), vs);
    let counts = m.counts(body);
    assert_eq!(
        (
            counts.lumps,
            counts.shells,
            counts.faces,
            counts.loops,
            counts.edges,
            counts.vertices
        ),
        (3, 5, 12, 12, 24, 16)
    );
    let c = cone(&mut m, r(3), r(0), r(4), &Pose::identity()).unwrap();
    let apex = m.vertices_of_body(c);
    assert_eq!(apex.len(), 1);
    assert_eq!(explicit(&m, apex[0]), p(0, 0, 4));
});

blind_test!(case_14_vertex_order_uses_coedge_sense_not_edge_sense, {
    let mut m = Model::new();
    let a = vertex(&mut m, p(0, 0, 0));
    let b = vertex(&mut m, p(1, 0, 0));
    let curve = m.geom_mut().add_curve(Curve::Line(Line {
        origin: p(0, 0, 0).into(),
        dir: v(1, 0, 0).into(),
    }));
    let e = m.add_edge(Edge {
        curve,
        sense: Sense::Reversed,
        start: Some(a),
        end: Some(b),
    });
    let l = cycle(&mut m, &[(e, Sense::Reversed), (e, Sense::Same)]);
    let f = face(
        &mut m,
        plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)),
        Sense::Same,
        vec![l],
    );
    let body = body_from_faces(&mut m, vec![f]);
    assert_eq!(m.vertices_of_body(body), vec![b, a]);
    assert_eq!(m.edges_of_body(body), vec![e]);
});

blind_test!(cases_15_16_two_lumps_touching_with_distinct_topology, {
    let mut m = Model::new();
    let shape = Primitive::Block([r(1), r(1), r(1)]);
    let a = expected(&shape, &Pose::identity());
    let b = expected(
        &shape,
        &Pose::new(
            Rotation::identity(),
            v(1, 1, 0).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let ba = assemble(&mut m, &a);
    let bb = assemble(&mut m, &b);
    let la = m.body(ba).lumps[0];
    let lb = m.body(bb).lumps[0];
    let body = body_from_lumps(&mut m, vec![la, lb]);
    sound(&m, body);
    let c = m.counts(body);
    assert_eq!(
        (c.lumps, c.shells, c.faces, c.edges, c.vertices),
        (2, 2, 12, 24, 16)
    );
    let seam = |m: &Model, b| {
        walk(m, b)
            .1
            .into_iter()
            .find(|&e| {
                let e = m.edge(e);
                explicit(m, e.start.unwrap()) == p(1, 1, 0)
                    && explicit(m, e.end.unwrap()) == p(1, 1, 1)
            })
            .unwrap()
    };
    let ea = seam(&m, ba);
    let eb = seam(&m, bb);
    assert_ne!(ea, eb);
    assert_ne!(m.edge(ea).start, m.edge(eb).start);
    assert!(
        m.geom()
            .curve(m.edge(ea).curve)
            .same_definition(m.geom().curve(m.edge(eb).curve))
    );
    // Also share the actual line handle, not just equal line definitions.
    let shared = m.edge(ea).curve;
    let old = m.edge(eb).clone();
    let replacement = m.add_edge(Edge {
        curve: shared,
        ..old
    });
    for vertex in [m.edge(replacement).start, m.edge(replacement).end]
        .into_iter()
        .flatten()
    {
        let point = m.vertex(vertex).point;
        record(
            &mut m,
            FactKind::Incident,
            GeomKey::Point(point),
            GeomKey::Curve(shared),
        );
    }
    let mut faces = vec![];
    for f in walk(&m, bb).0 {
        let f = m.face(f).clone();
        let loops = f
            .loops
            .iter()
            .map(|&l| {
                let uses = cs(&m, l)
                    .iter()
                    .map(|&c| {
                        let c = m.coedge(c);
                        (if c.edge == eb { replacement } else { c.edge }, c.sense)
                    })
                    .collect::<Vec<_>>();
                cycle(&mut m, &uses)
            })
            .collect();
        faces.push(face_on(&mut m, f.surface, f.sense, loops, true));
    }
    let rebuilt = body_from_faces(&mut m, faces);
    let lb = m.body(rebuilt).lumps[0];
    let shared_body = body_from_lumps(&mut m, vec![la, lb]);
    sound(&m, shared_body);
    assert_eq!(m.counts(shared_body), c);
    assert_eq!(euler(&m, shell_of(&m, ba)), a.chi);
    assert_eq!(euler(&m, shell_of(&m, bb)), b.chi);
});

blind_test!(case_17_void_shell_reverses_faces_and_boundary_travel, {
    let mut m = Model::new();
    let outer = assemble(
        &mut m,
        &expected(&Primitive::Block([r(10), r(10), r(10)]), &Pose::identity()),
    );
    let mut inner = expected(
        &Primitive::Block([r(2), r(2), r(2)]),
        &Pose::new(
            Rotation::identity(),
            v(4, 4, 4).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    for (_, sense, loops) in &mut inner.faces {
        *sense = Sense::Reversed;
        for l in loops {
            if let Boundary::Cycle(cs) = l {
                cs.reverse();
                for (_, sense) in cs {
                    *sense = flip(*sense);
                }
            }
        }
    }
    let hole = assemble(&mut m, &inner);
    let shells = vec![shell_of(&m, outer), shell_of(&m, hole)];
    let body = body_from_shells(&mut m, shells.clone());
    sound(&m, body);
    let c = m.counts(body);
    assert_eq!(
        (c.lumps, c.shells, c.faces, c.loops, c.edges, c.vertices),
        (1, 2, 12, 12, 24, 16)
    );
    for s in shells {
        assert_eq!(euler(&m, s), 2);
    }
});

blind_test!(case_18_sphere_split_at_equator, {
    let mut m = Model::new();
    let s = m.geom_mut().add_surface(sphere_surface(r(2)));
    let e = edge(&mut m, circle(p(0, 0, 0), r(2)), None, None);
    let north = cycle(&mut m, &[(e, Sense::Same)]);
    let south = cycle(&mut m, &[(e, Sense::Reversed)]);
    let n = face_on(&mut m, s, Sense::Same, vec![north], true);
    let s = face_on(&mut m, s, Sense::Same, vec![south], true);
    let body = body_from_faces(&mut m, vec![n, s]);
    sound(&m, body);
    assert_eq!(euler(&m, shell_of(&m, body)), 2);
    let c = m.counts(body);
    assert_eq!(
        (c.faces, c.loops, c.edges, c.ring_edges, c.vertices),
        (2, 2, 1, 1, 0)
    );
    assert_eq!(c.faces_by_kind, [(SurfaceKind::Sphere, 2)].into());
});

blind_test!(case_19_torus_cut_along_meridian, {
    let mut m = Model::new();
    let e = edge(
        &mut m,
        Curve::Circle(Circle {
            center: p(3, 0, 0).into(),
            normal: v(0, 1, 0).into(),
            radius_sq: r(1).into(),
            ref_dir: v(1, 0, 0).into(),
        }),
        None,
        None,
    );
    let a = cycle(&mut m, &[(e, Sense::Same)]);
    let b = cycle(&mut m, &[(e, Sense::Reversed)]);
    let f = face(&mut m, torus_surface(), Sense::Same, vec![a, b]);
    let body = body_from_faces(&mut m, vec![f]);
    sound(&m, body);
    assert_eq!(euler(&m, shell_of(&m, body)), 0);
    let c = m.counts(body);
    assert_eq!(
        (c.faces, c.loops, c.edges, c.ring_edges, c.vertices),
        (1, 2, 1, 1, 0)
    );
});

blind_test!(case_20_block_with_cylindrical_hole, {
    let mut d = expected(&Primitive::Block([r(4), r(4), r(4)]), &Pose::identity());
    let bottom = d.edges.len();
    d.edges.push((circle(p(2, 2, 0), r(1)), None, None));
    let top = d.edges.len();
    d.edges.push((circle(p(2, 2, 4), r(1)), None, None));
    d.faces[4]
        .2
        .push(Boundary::Cycle(vec![(bottom, Sense::Same)]));
    d.faces[5]
        .2
        .push(Boundary::Cycle(vec![(top, Sense::Reversed)]));
    d.faces.push((
        Surface::Cylinder(Cylinder {
            origin: p(2, 2, 0).into(),
            axis: v(0, 0, 1).into(),
            radius_sq: r(1).into(),
            ref_dir: v(1, 0, 0).into(),
        }),
        Sense::Reversed,
        vec![
            Boundary::Cycle(vec![(bottom, Sense::Reversed)]),
            Boundary::Cycle(vec![(top, Sense::Same)]),
        ],
    ));
    let mut m = Model::new();
    let body = assemble(&mut m, &d);
    sound(&m, body);
    assert_eq!(
        m.counts(body),
        EntityCounts {
            lumps: 1,
            shells: 1,
            faces: 7,
            loops: 10,
            edges: 14,
            ring_edges: 2,
            vertices: 8,
            faces_by_kind: [(SurfaceKind::Plane, 6), (SurfaceKind::Cylinder, 1)].into(),
            edges_by_kind: [(CurveKind::Line, 12), (CurveKind::Circle, 2)].into()
        }
    );
    // Four side disks contribute 4, the two punctured caps and cylinder each 0.
    assert_eq!(euler(&m, shell_of(&m, body)), 8 + 2 - 14 + 4);
});

blind_test!(closed_curve_with_one_vertex_is_sound, {
    let mut m = Model::new();
    let v = vertex(&mut m, p(2, 0, 0));
    let e = edge(&mut m, circle(p(0, 0, 0), r(2)), Some(v), Some(v));
    let s = m.geom_mut().add_surface(sphere_surface(r(2)));
    let a = cycle(&mut m, &[(e, Sense::Same)]);
    let b = cycle(&mut m, &[(e, Sense::Reversed)]);
    let a = face_on(&mut m, s, Sense::Same, vec![a], true);
    let b = face_on(&mut m, s, Sense::Same, vec![b], true);
    let body = body_from_faces(&mut m, vec![a, b]);
    sound(&m, body);
    assert_eq!(m.counts(body).ring_edges, 0);
    assert_eq!(m.vertices_of_body(body), vec![v]);
    assert_eq!(euler(&m, shell_of(&m, body)), 2);
});
