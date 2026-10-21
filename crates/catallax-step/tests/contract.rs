//! Blind STEP text, topology, determinism and refusal contracts.
mod common;
use common::*;
use std::io::{self, Write};
use std::time::{Duration, Instant};

blind_test!(case01_block_fixed_entities_and_exact_counts, {
    let mut m = Model::new();
    let b = block(&mut m, [n(2), n(3), n(5)], &Pose::identity()).unwrap();
    let o = StepOptions::default();
    assert_eq!(o, options(17));
    let (_, a, s) = export(&m, &[b], &o);
    assert_eq!(s.entities.len(), 170);
    s.fixed("catallax", "1.E-7");
    counts(
        &s,
        &[
            ("MANIFOLD_SOLID_BREP", 1),
            ("CLOSED_SHELL", 1),
            ("ADVANCED_FACE", 6),
            ("FACE_BOUND", 6),
            ("EDGE_LOOP", 6),
            ("ORIENTED_EDGE", 24),
            ("EDGE_CURVE", 12),
            ("VERTEX_POINT", 8),
            ("CARTESIAN_POINT", 27),
            ("LINE", 12),
            ("VECTOR", 12),
            ("DIRECTION", 26),
            ("AXIS2_PLACEMENT_3D", 7),
            ("PLANE", 6),
        ],
    );
    report(&a, [1, 6, 12, 8, 0, 0, 0], true);
    let corners = (0..8)
        .map(|i| {
            [
                n(if i & 1 == 0 { 0 } else { 2 }),
                n(if i & 2 == 0 { 0 } else { 3 }),
                n(if i & 4 == 0 { 0 } else { 5 }),
            ]
        })
        .collect();
    assert_eq!(vertex_set(&s), corners);
});
blind_test!(
    case02_depth_first_numbering_in_all_shapes_and_unreachable_store_entries,
    {
        // export() always checks the complete reference graph and numbering.
        let mut m = Model::new();
        let _unselected = block(&mut m, [n(11), n(13), n(17)], &Pose::identity()).unwrap();
        let b = sphere(&mut m, n(5), &Pose::identity()).unwrap();
        let (bytes, a, _) = export(&m, &[b], &options(17));
        report(&a, [1, 1, 1, 2, 2, 1, 0], true);
        let mut clean = Model::new();
        let b = sphere(&mut clean, n(5), &Pose::identity()).unwrap();
        let (other, _, _) = export(&clean, &[b], &options(17));
        assert_eq!(bytes, other);
    }
);
blind_test!(case03_cylinder_ring_vertices_and_seams, {
    let mut m = Model::new();
    let b = cylinder(&mut m, n(2), n(5), &Pose::identity()).unwrap();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [1, 3, 2, 2, 2, 0, 0], true);
    counts(
        &s,
        &[
            ("ADVANCED_FACE", 3),
            ("FACE_BOUND", 4),
            ("EDGE_LOOP", 4),
            ("ORIENTED_EDGE", 4),
            ("CIRCLE", 2),
            ("PLANE", 2),
            ("CYLINDRICAL_SURFACE", 1),
        ],
    );
    assert_eq!(
        vertex_set(&s),
        BTreeSet::from([[n(2), n(0), n(0)], [n(2), n(0), n(5)]])
    );
    for id in s.ids("EDGE_CURVE") {
        let e = &s.entities[&id];
        assert_eq!(e.args[1], e.args[2]);
        assert_eq!(s.e(reference(&e.args[3]), "CIRCLE").args[2], "2.");
    }
});
fn loop_edges(s: &Step) -> Vec<(usize, String)> {
    let id = s.ids("EDGE_LOOP");
    assert_eq!(id.len(), 1);
    list(&s.entities[&id[0]].args[1])
        .iter()
        .map(|v| {
            let e = s.e(reference(v), "ORIENTED_EDGE");
            (reference(&e.args[3]), e.args[4].clone())
        })
        .collect()
}
fn circle_frame(
    s: &Step,
    edge: usize,
    center: [Rational; 3],
    axis: [Rational; 3],
    reference_dir: [Rational; 3],
    radius: &str,
) {
    let c = s.e(reference(&s.e(edge, "EDGE_CURVE").args[3]), "CIRCLE");
    assert_eq!(c.args[2], radius);
    assert_eq!(
        s.placement(reference(&c.args[1])),
        (center, axis, reference_dir)
    );
}
blind_test!(case04_sphere_poles_and_half_meridian, {
    let mut m = Model::new();
    let b = sphere(&mut m, n(5), &Pose::identity()).unwrap();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [1, 1, 1, 2, 2, 1, 0], true);
    counts(
        &s,
        &[
            ("SPHERICAL_SURFACE", 1),
            ("CIRCLE", 1),
            ("ORIENTED_EDGE", 2),
        ],
    );
    let es = loop_edges(&s);
    assert_eq!(es, vec![(es[0].0, ".T.".into()), (es[0].0, ".F.".into())]);
    let e = s.e(es[0].0, "EDGE_CURVE");
    assert_eq!(s.xyz(s.vertex(reference(&e.args[1]))), [n(0), n(0), n(-5)]);
    assert_eq!(s.xyz(s.vertex(reference(&e.args[2]))), [n(0), n(0), n(5)]);
    circle_frame(
        &s,
        es[0].0,
        [n(0), n(0), n(0)],
        [n(0), n(-1), n(0)],
        [n(1), n(0), n(0)],
        "5.",
    );
});
blind_test!(case05_torus_two_closed_seams, {
    let mut m = Model::new();
    let b = torus(&mut m, n(3), n(1), &Pose::identity()).unwrap();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [1, 1, 2, 1, 1, 2, 0], true);
    assert_eq!(vertex_set(&s), BTreeSet::from([[n(4), n(0), n(0)]]));
    let es = loop_edges(&s);
    assert_eq!(es.len(), 4);
    assert_ne!(es[0].0, es[1].0);
    assert_eq!(
        es,
        vec![
            (es[0].0, ".T.".into()),
            (es[1].0, ".T.".into()),
            (es[0].0, ".F.".into()),
            (es[1].0, ".F.".into())
        ]
    );
    circle_frame(
        &s,
        es[0].0,
        [n(3), n(0), n(0)],
        [n(0), n(-1), n(0)],
        [n(1), n(0), n(0)],
        "1.",
    );
    circle_frame(
        &s,
        es[1].0,
        [n(0), n(0), n(0)],
        [n(0), n(0), n(1)],
        [n(1), n(0), n(0)],
        "4.",
    );
    for (id, _) in es {
        let e = s.e(id, "EDGE_CURVE");
        assert_eq!(e.args[1], e.args[2]);
    }
});
blind_test!(case06_cone_apex_and_transcendental_angle, {
    let mut m = Model::new();
    let b = cone(&mut m, n(3), n(0), n(4), &Pose::identity()).unwrap();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [1, 2, 1, 2, 1, 0, 0], false);
    assert!(a.value.max_deviation > n(0));
    assert!(a.value.max_deviation <= r("1e-7"));
    assert!(a.marks.is_empty());
    let id = s.ids("VERTEX_LOOP");
    assert_eq!(id.len(), 1);
    assert_eq!(
        s.xyz(s.vertex(reference(&s.entities[&id[0]].args[1]))),
        [n(0), n(0), n(4)]
    );
    let id = s.ids("CONICAL_SURFACE")[0];
    let e = &s.entities[&id];
    assert_eq!(e.args[2], "0.");
    assert_eq!(
        s.placement(reference(&e.args[1])),
        ([n(0), n(0), n(4)], [n(0), n(0), n(-1)], [n(1), n(0), n(0)])
    );
    // Exact alternating arctan series, independently checked by verify_exact.py.
    assert!(["0.64350110879328438", "0.64350110879328439"].contains(&e.args[3].as_str()));
});
blind_test!(case07_step_strings_and_non_bmp_utf16, {
    for name in [
        "a'b\\c",
        "é",
        "😀",
        "a\n\t\u{7f}z",
        "#21,('x')",
        "",
        "\\X2\\00E9\\X0\\",
    ] {
        let mut o = options(17);
        o.product_name = name.into();
        let (_, _, s) = export(&Model::new(), &[], &o);
        assert_eq!(s.entities.len(), 20);
    }
});
blind_test!(case08_empty_body_lists_and_mixed_empty_bodies, {
    let mut m = Model::new();
    let a = empty(&mut m);
    let b = empty(&mut m);
    let mut previous = None;
    for bs in [vec![], vec![a], vec![a, b]] {
        let (bytes, a, s) = export(&m, &bs, &options(17));
        report(&a, [0, 0, 0, 0, 0, 0, 0], true);
        assert_eq!(s.entities.len(), 20);
        assert_eq!(s.items(), vec![17]);
        if let Some(p) = previous {
            assert_eq!(bytes, p);
        }
        previous = Some(bytes);
    }
    let solid = cylinder(&mut m, n(2), n(5), &Pose::identity()).unwrap();
    assert_eq!(
        export(&m, &[a, solid, b], &options(17)).0,
        export(&m, &[solid], &options(17)).0
    );
});
blind_test!(case09_two_lumps_and_touching_shared_point_definitions, {
    let (m, b) = touching_blocks();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [2, 12, 24, 16, 0, 0, 0], true);
    assert_eq!(s.count("MANIFOLD_SOLID_BREP"), 2);
    assert_eq!(vertex_set(&s).len(), 14);
    // Two common corner PointIds still need four VERTEX_POINT/CARTESIAN_POINT pairs.
    for p in [[n(1), n(1), n(0)], [n(1), n(1), n(1)]] {
        assert_eq!(
            s.ids("VERTEX_POINT")
                .iter()
                .filter(|&&id| s.xyz(s.vertex(id)) == p)
                .count(),
            2
        );
    }
});
blind_test!(case10_void_shell_reverses_model_face_senses, {
    let (m, b) = void_block();
    let (_, a, s) = export(&m, &[b], &options(17));
    report(&a, [1, 12, 24, 16, 0, 0, 0], true);
    counts(
        &s,
        &[
            ("BREP_WITH_VOIDS", 1),
            ("MANIFOLD_SOLID_BREP", 0),
            ("ORIENTED_CLOSED_SHELL", 1),
            ("CLOSED_SHELL", 2),
        ],
    );
    let v = &s.entities[&s.ids("ORIENTED_CLOSED_SHELL")[0]];
    assert_eq!(v.args[3], ".F.");
    let shell = s.e(reference(&v.args[2]), "CLOSED_SHELL");
    for face in list(&shell.args[1]) {
        let f = s.e(reference(face), "ADVANCED_FACE");
        // Architect decision step-void-face-sense: general rule stands; the
        // original case 10's .F. was a typo. Reversed model -> negated -> .T.
        assert_eq!(f.args[3], ".T.");
        for bound in list(&f.args[1]) {
            assert_eq!(s.e(reference(bound), "FACE_BOUND").args[2], ".F.");
        }
    }
});
blind_test!(case11_bodies_keep_lump_order_and_geometry_sharing, {
    let mut m = Model::new();
    let a = block(&mut m, [n(2), n(3), n(5)], &Pose::identity()).unwrap();
    let b = sphere(&mut m, n(5), &translated("20", "0", "0")).unwrap();
    let (_, report, s) = export(&m, &[b, a], &options(17));
    assert_eq!(report.value.faces, 7);
    let items = s.items();
    let first = s.e(items[1], "MANIFOLD_SOLID_BREP");
    let shell = s.e(reference(&first.args[1]), "CLOSED_SHELL");
    assert_eq!(list(&shell.args[1]).len(), 1);
    let (m, b) = split_sphere();
    let (_, a, s) = export(&m, &[b], &options(17));
    report_shared_sphere(&a, &s);
});
fn report_shared_sphere(a: &Answer<StepReport>, s: &Step) {
    report(a, [1, 2, 1, 1, 1, 0, 0], true);
    counts(
        s,
        &[("SPHERICAL_SURFACE", 1), ("CIRCLE", 1), ("EDGE_CURVE", 1)],
    );
}
blind_test!(case19_repeat_and_eight_cloned_models, {
    for primitive in primitives() {
        let mut m = Model::new();
        let b = primitive.build(&mut m, &translated("1000000", "-7", "1/3"));
        let o = options(17);
        let (bytes, a, _) = export(&m, &[b], &o);
        let (again, a2, _) = export(&m, &[b], &o);
        assert_eq!(bytes, again);
        assert_eq!(a.value, a2.value);
        assert_eq!(a.marks, a2.marks);
        let threads = (0..8)
            .map(|_| {
                let m = m.clone();
                let o = o.clone();
                std::thread::spawn(move || {
                    let (b, a, _) = export(&m, &[b], &o);
                    (b, a)
                })
            })
            .collect::<Vec<_>>();
        for thread in threads {
            let (b, other) = thread.join().unwrap();
            assert_eq!(b, bytes);
            assert_eq!(other.value, a.value);
            assert_eq!(other.marks, a.marks);
        }
    }
});
blind_test!(case20_expired_deadline_is_atomic_or_identical, {
    for primitive in primitives() {
        let mut m = Model::new();
        let b = primitive.build(&mut m, &translated("1/3", "-7", "2/3"));
        let o = options(17);
        let (expected, a, _) = export(&m, &[b], &o);
        let mut dec = Decider::default();
        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
        let sentinel = b"already-written".to_vec();
        let mut bytes = sentinel.clone();
        match write_step(&m, &[b], &o, &mut dec, &mut bytes) {
            Ok(actual) => {
                assert_eq!(&bytes[sentinel.len()..], expected);
                assert_eq!(actual.value, a.value);
                assert_eq!(actual.marks, a.marks);
            }
            Err(e) => {
                assert_eq!(e.code, "numbers.deadline");
                assert_eq!(bytes, sentinel);
            }
        }
    }
});
struct BrokenWriter {
    limit: usize,
    written: usize,
}
impl Write for BrokenWriter {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self.written == self.limit {
            return Err(io::Error::other("step-test injected write failure"));
        }
        let count = b.len().min(self.limit - self.written);
        self.written += count;
        Ok(count)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
blind_test!(io_errors_keep_exact_code_and_detail_after_short_writes, {
    for limit in [0, 1, 31, 200] {
        let mut writer = BrokenWriter { limit, written: 0 };
        let e = write_step(
            &Model::new(),
            &[],
            &options(17),
            &mut Decider::default(),
            &mut writer,
        )
        .unwrap_err();
        assert_eq!(e.code, "step.io-error");
        assert_eq!(e.detail, "step-test injected write failure");
        assert_eq!(writer.written, limit);
    }
});
blind_test!(
    nonpositive_digit_count_clamps_to_one_and_uncertainty_rounds_up,
    {
        let mut m = Model::new();
        let b = block(&mut m, [r("1/3"), r("2/3"), r("5/4")], &Pose::identity()).unwrap();
        let mut o = options(0);
        o.curve_tolerance = r("1/3");
        let first = export(&m, &[b], &o).0;
        o.significant_digits = 1;
        let second = export(&m, &[b], &o).0;
        assert_eq!(first, second);
        for digits in [1, 5, 17, 40] {
            let mut o = options(digits);
            o.curve_tolerance = r("1/3");
            let (_, _, s) = export(&Model::new(), &[], &o);
            s.fixed("catallax", &decimal(&r("1/3"), digits, true));
        }
    }
);
blind_test!(body_marks_are_not_export_marks, {
    let mut m = Model::new();
    let b = block(&mut m, [n(2), n(3), n(5)], &Pose::identity()).unwrap();
    let mut copy = m.body(b).clone();
    copy.marks.push(Mark::new(
        "boolean.test-input",
        BoundKind::Distance,
        n(1),
        "input provenance",
    ));
    let b = m.add_body(copy);
    let (_, a, _) = export(&m, &[b], &options(17));
    assert!(a.marks.is_empty());
    assert_eq!(m.body(b).marks.len(), 1);
});
