//! Planted faults checked against the architect's exact fault codes.
mod common;
use common::*;

fn box_model() -> (Model, BodyId) {
    let mut m = Model::new();
    let b = assemble(&mut m, &expected(&primitives()[0], &Pose::identity()));
    sound(&m, b);
    (m, b)
}
fn replace_face(m: &mut Model, b: BodyId, index: usize, loops: Vec<LoopId>) -> BodyId {
    let mut fs = m.shell(shell_of(m, b)).faces.clone();
    let old = m.face(fs[index]).clone();
    fs[index] = m.add_face(Face { loops, ..old });
    body_from_faces(m, fs)
}

blind_test!(case_21_missing_block_face, {
    let (mut m, b) = box_model();
    let mut fs = m.shell(shell_of(&m, b)).faces.clone();
    fs.pop();
    let bad = body_from_faces(&mut m, fs);
    assert_eq!(euler(&m, shell_of(&m, bad)), 1);
    faults(&m, bad, &["topology.edge-use", "topology.euler"]);
    sound(&m, b);
});
blind_test!(case_22_flipped_coedge, {
    let (mut m, b) = box_model();
    let f = m.shell(shell_of(&m, b)).faces[0];
    let mut co = cs(&m, m.face(f).loops[0]);
    let old = m.coedge(co[0]).clone();
    co[0] = m.add_coedge(Coedge {
        sense: flip(old.sense),
        ..old
    });
    let l = m.add_loop(Loop::Coedges(co));
    let bad = replace_face(&mut m, b, 0, vec![l]);
    faults(&m, bad, &["topology.loop-not-closed", "topology.edge-use"]);
    sound(&m, b);
});
blind_test!(case_23_missing_coedge, {
    let (mut m, b) = box_model();
    let f = m.shell(shell_of(&m, b)).faces[0];
    let mut co = cs(&m, m.face(f).loops[0]);
    co.pop();
    let l = m.add_loop(Loop::Coedges(co));
    let bad = replace_face(&mut m, b, 0, vec![l]);
    faults(&m, bad, &["topology.loop-not-closed"]);
});
blind_test!(case_24_duplicate_face_and_shared_loop, {
    for new_face in [false, true] {
        let (mut m, b) = box_model();
        let mut fs = m.shell(shell_of(&m, b)).faces.clone();
        let f = if new_face {
            m.add_face(m.face(fs[0]).clone())
        } else {
            fs[0]
        };
        fs.push(f);
        let bad = body_from_faces(&mut m, fs);
        faults(&m, bad, &["topology.shared-entity"]);
    }
});
blind_test!(case_25_invalid_edge_ends, {
    // Both one-ended variants, a ring on a line, and a one-vertex line edge.
    for ends in 0..4 {
        let mut m = Model::new();
        let vertex = vertex(&mut m, p(0, 0, 0));
        let (a, b) = match ends {
            0 => (None, None),
            1 => (Some(vertex), None),
            2 => (None, Some(vertex)),
            _ => (Some(vertex), Some(vertex)),
        };
        let e = edge(
            &mut m,
            Curve::Line(Line {
                origin: p(0, 0, 0).into(),
                dir: v(1, 0, 0).into(),
            }),
            a,
            b,
        );
        let l = cycle(&mut m, &[(e, Sense::Same)]);
        let f = face(
            &mut m,
            plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)),
            Sense::Same,
            vec![l],
        );
        let bad = body_from_faces(&mut m, vec![f]);
        faults(&m, bad, &["topology.edge-ends"]);
    }
});
blind_test!(case_26_unbounded_plane_cylinder_and_cone, {
    let surfaces = [
        plane(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0)),
        Surface::Cylinder(Cylinder {
            origin: p(0, 0, 0).into(),
            axis: v(0, 0, 1).into(),
            radius_sq: r(4).into(),
            ref_dir: v(1, 0, 0).into(),
        }),
        Surface::Cone(Cone {
            apex: p(0, 0, 4).into(),
            axis: v(0, 0, -1).into(),
            slope: q(3, 4).into(),
            ref_dir: v(1, 0, 0).into(),
        }),
    ];
    for surface in surfaces {
        let mut m = Model::new();
        let f = face(&mut m, surface, Sense::Same, vec![]);
        let bad = body_from_faces(&mut m, vec![f]);
        faults(&m, bad, &["topology.unbounded-face"]);
    }
});
blind_test!(case_27_vertex_only_boundary_and_repeated_vertex, {
    let mut m = Model::new();
    let apex = vertex(&mut m, p(0, 0, 4));
    let l = m.add_loop(Loop::Vertex(apex));
    let f = face(
        &mut m,
        Surface::Cone(Cone {
            apex: p(0, 0, 4).into(),
            axis: v(0, 0, -1).into(),
            slope: q(3, 4).into(),
            ref_dir: v(1, 0, 0).into(),
        }),
        Sense::Same,
        vec![l],
    );
    let bad = body_from_faces(&mut m, vec![f]);
    faults(&m, bad, &["topology.vertex-loop"]);
    let (mut m, b) = box_model();
    let f = m.shell(shell_of(&m, b)).faces[0];
    let loops = m.face(f).loops.clone();
    let co = cs(&m, loops[0]);
    let vertex = m.edge(m.coedge(co[0]).edge).start.unwrap();
    let point_loop = m.add_loop(Loop::Vertex(vertex));
    let bad = replace_face(&mut m, b, 0, vec![loops[0], point_loop]);
    faults(&m, bad, &["topology.vertex-loop"]);
    let second = m.add_loop(Loop::Vertex(vertex));
    let bad = replace_face(&mut m, b, 0, vec![loops[0], point_loop, second]);
    faults(&m, bad, &["topology.vertex-loop"]);
});
blind_test!(case_28_unrecorded_replacement_circle, {
    let mut m = Model::new();
    let b = assemble(&mut m, &expected(&primitives()[1], &Pose::identity()));
    let fs = m.shell(shell_of(&m, b)).faces.clone();
    let top_loop = m.face(fs[1]).loops[0];
    let old = m.coedge(cs(&m, top_loop)[0]).edge;
    let definition = m.geom().curve(m.edge(old).curve).clone();
    let curve = m.geom_mut().add_curve(definition);
    let replacement = m.add_edge(Edge {
        curve,
        ..m.edge(old).clone()
    });
    let mut rebuilt = vec![];
    for f in fs {
        let f = m.face(f).clone();
        let mut loops = vec![];
        for l in f.loops {
            let uses = cs(&m, l)
                .into_iter()
                .map(|c| {
                    let c = m.coedge(c);
                    (if c.edge == old { replacement } else { c.edge }, c.sense)
                })
                .collect::<Vec<_>>();
            loops.push(cycle(&mut m, &uses));
        }
        rebuilt.push(face_on(&mut m, f.surface, f.sense, loops, false));
    }
    let bad = body_from_faces(&mut m, rebuilt);
    faults(&m, bad, &["topology.missing-incidence"]);
});
blind_test!(case_29_disconnected_sphere_faces_have_chi_four, {
    let mut m = Model::new();
    let s = m.geom_mut().add_surface(sphere_surface(r(2)));
    let a = face_on(&mut m, s, Sense::Same, vec![], true);
    let b = face_on(&mut m, s, Sense::Same, vec![], true);
    let bad = body_from_faces(&mut m, vec![a, b]);
    assert_eq!(euler(&m, shell_of(&m, bad)), 4);
    faults(&m, bad, &["topology.euler", "topology.disconnected-shell"]);
});
blind_test!(case_30_disconnected_blocks_in_one_shell, {
    let (mut m, a) = box_model();
    let b = assemble(
        &mut m,
        &expected(
            &primitives()[0],
            &Pose::new(
                Rotation::identity(),
                v(20, 0, 0).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
    );
    let mut fs = walk(&m, a).0;
    fs.extend(walk(&m, b).0);
    let bad = body_from_faces(&mut m, fs);
    faults(&m, bad, &["topology.disconnected-shell"]);
});
blind_test!(case_31_repeated_shell_and_lump, {
    let (mut m, b) = box_model();
    let lump = m.body(b).lumps[0];
    let shell = shell_of(&m, b);
    let other = m.add_lump(Lump {
        shells: vec![shell],
    });
    let a = body_from_lumps(&mut m, vec![lump, other]);
    faults(&m, a, &["topology.shared-entity"]);
    let a = body_from_lumps(&mut m, vec![lump, lump]);
    faults(&m, a, &["topology.shared-entity"]);
});
blind_test!(case_32_body_kind_and_empty_parts, {
    let mut m = Model::new();
    let b = body_from_lumps(&mut m, vec![]);
    faults(&m, b, &["topology.body-kind"]);
    let sphere = assemble(
        &mut m,
        &expected(&Primitive::Sphere(r(2)), &Pose::identity()),
    );
    let lumps = m.body(sphere).lumps.clone();
    let b = m.add_body(Body {
        kind: BodyKind::Empty,
        lumps,
        marks: MarkLog::new(),
    });
    faults(&m, b, &["topology.body-kind"]);
    let lump = m.add_lump(Lump { shells: vec![] });
    let b = body_from_lumps(&mut m, vec![lump]);
    faults(&m, b, &["topology.empty-part"]);
    let shell = m.add_shell(Shell { faces: vec![] });
    let b = body_from_shells(&mut m, vec![shell]);
    faults(&m, b, &["topology.empty-part"]);
    let l = m.add_loop(Loop::Coedges(vec![]));
    let f = face(&mut m, sphere_surface(r(2)), Sense::Same, vec![l]);
    let b = body_from_faces(&mut m, vec![f]);
    faults(&m, b, &["topology.empty-part"]);
});

blind_test!(ring_coedge_must_be_alone_in_its_loop, {
    let mut m = Model::new();
    let e = edge(&mut m, circle(p(0, 0, 0), r(2)), None, None);
    let l = cycle(&mut m, &[(e, Sense::Same), (e, Sense::Reversed)]);
    let f = face(&mut m, sphere_surface(r(2)), Sense::Same, vec![l]);
    let b = body_from_faces(&mut m, vec![f]);
    faults(&m, b, &["topology.loop-not-closed"]);
});
blind_test!(coedge_repeated_across_distinct_loops, {
    let mut m = Model::new();
    let e = edge(&mut m, circle(p(0, 0, 0), r(2)), None, None);
    let c = m.add_coedge(Coedge {
        edge: e,
        sense: Sense::Same,
    });
    let a = m.add_loop(Loop::Coedges(vec![c]));
    let b = m.add_loop(Loop::Coedges(vec![c]));
    let f = face(&mut m, sphere_surface(r(2)), Sense::Same, vec![a, b]);
    let b = body_from_faces(&mut m, vec![f]);
    faults(&m, b, &["topology.shared-entity"]);
});
blind_test!(distinct_shells_cannot_share_edges_or_vertices, {
    let (mut m, b) = box_model();
    let original = walk(&m, b).0;
    let mut faces = vec![];
    for f in original {
        let f = m.face(f).clone();
        let mut loops = vec![];
        for l in f.loops {
            let uses = cs(&m, l)
                .into_iter()
                .map(|c| {
                    let c = m.coedge(c);
                    (c.edge, c.sense)
                })
                .collect::<Vec<_>>();
            loops.push(cycle(&mut m, &uses));
        }
        faces.push(face_on(&mut m, f.surface, f.sense, loops, true));
    }
    let new = m.add_shell(Shell { faces });
    let old = shell_of(&m, b);
    let bad = body_from_shells(&mut m, vec![old, new]);
    faults(&m, bad, &["topology.shared-entity"]);
    // Isolate shared vertices: clone every edge and every higher entity.
    let mut faces = vec![];
    let mut edges = std::collections::BTreeMap::new();
    for e in walk(&m, b).1 {
        edges.insert(e, m.add_edge(m.edge(e).clone()));
    }
    for f in walk(&m, b).0 {
        let f = m.face(f).clone();
        let mut loops = vec![];
        for l in f.loops {
            let uses = cs(&m, l)
                .into_iter()
                .map(|c| {
                    let c = m.coedge(c);
                    (edges[&c.edge], c.sense)
                })
                .collect::<Vec<_>>();
            loops.push(cycle(&mut m, &uses));
        }
        faces.push(face_on(&mut m, f.surface, f.sense, loops, true));
    }
    let new = m.add_shell(Shell { faces });
    let bad = body_from_shells(&mut m, vec![old, new]);
    faults(&m, bad, &["topology.shared-entity"]);
});
blind_test!(incidence_is_required_for_endpoint_and_vertex_loop, {
    for vertex_loop in [false, true] {
        let mut m = Model::new();
        let b = assemble(&mut m, &expected(&primitives()[5], &Pose::identity()));
        let old = m.vertices_of_body(b)[0];
        let point = m
            .geom_mut()
            .add_point(PointDef::Explicit(p(0, 0, 4).into()));
        let replacement = m.add_vertex(Vertex { point });
        if vertex_loop {
            let fs = m.shell(shell_of(&m, b)).faces.clone();
            let side = m.face(fs[1]).clone();
            let l = m.add_loop(Loop::Vertex(replacement));
            let f = face_on(
                &mut m,
                side.surface,
                side.sense,
                vec![side.loops[0], l],
                false,
            );
            let bad = body_from_faces(&mut m, vec![fs[0], f]);
            faults(&m, bad, &["topology.missing-incidence"]);
        } else {
            let curve = m.geom_mut().add_curve(circle(p(0, 0, 0), r(4)));
            let e = m.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(old),
                end: Some(replacement),
            });
            let l = cycle(&mut m, &[(e, Sense::Same), (e, Sense::Reversed)]);
            let f = face(&mut m, sphere_surface(r(4)), Sense::Same, vec![l]);
            let bad = body_from_faces(&mut m, vec![f]);
            faults(&m, bad, &["topology.missing-incidence"]);
        }
    }
});
blind_test!(equal_positions_do_not_close_distinct_vertex_handles, {
    let (mut m, b) = box_model();
    let f = m.shell(shell_of(&m, b)).faces[0];
    let co = cs(&m, m.face(f).loops[0]);
    let original = m.coedge(co[0]).clone();
    let e = m.edge(original.edge).clone();
    let a = e.start.unwrap();
    let point = explicit(&m, a).clone();
    let clone = vertex(&mut m, point);
    let e = m.add_edge(Edge {
        start: Some(clone),
        ..e
    });
    let mut co = co;
    co[0] = m.add_coedge(Coedge {
        edge: e,
        ..original
    });
    let l = m.add_loop(Loop::Coedges(co));
    let bad = replace_face(&mut m, b, 0, vec![l]);
    faults(&m, bad, &["topology.loop-not-closed"]);
});
blind_test!(check_trusts_facts_and_does_not_validate_coordinates, {
    let mut d = expected(&primitives()[0], &Pose::identity());
    d.points[0] = p(99, 98, 97);
    // Deliberately false fixture facts isolate the documented trust boundary:
    // check is topological, not a second geometric decision engine.
    let mut m = Model::new();
    let b = assemble(&mut m, &d);
    sound(&m, b);
});
