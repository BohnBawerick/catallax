use crate::{
    BoolOp, Containment, FaceClass, Stage, boolean, classify, classify_point, face_point, imprint,
    intersect,
};
use catallax_numbers::{BudgetPolicy, Decider, DecisionStats, Rational, Real};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use catallax_topology::{BodyId, BodyKind, Model, block, cylinder, empty, sphere};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn pose(x: i64, y: i64, z: i64) -> Pose {
    Pose::new(
        Rotation::identity(),
        (RVector3::new(r(x), r(y), r(z))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
fn counts(model: &Model, body: BodyId) -> (usize, usize, usize, usize, usize, usize) {
    let c = model.counts(body);
    (
        c.lumps,
        c.shells,
        c.faces,
        c.edges,
        c.ring_edges,
        c.vertices,
    )
}
fn decide(model: &mut Model, body: BodyId, p: &RPoint3, inside: bool, dec: &mut Decider) {
    assert_eq!(
        classify_point(model, body, p, dec).proven(),
        Some(&if inside {
            Containment::Inside
        } else {
            Containment::Outside
        }),
        "query {p:?}"
    );
}
fn operation(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    expected: (usize, usize, usize, usize, usize, usize),
    dec: &mut Decider,
) -> BodyId {
    let answer = boolean(model, a, b, op, dec)
        .unwrap_or_else(|e| panic!("{op:?} refused at {:?}: {}", e.stage, e.refusal));
    assert_eq!(counts(model, answer.body), expected, "{op:?}");
    assert!(model.check(answer.body).is_empty());
    assert!(model.check(a).is_empty());
    assert!(model.check(b).is_empty());
    answer.body
}

#[test]
fn prepared_face_queries_keep_opposite_spherical_patches_and_closed_rims_distinct() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let stock = sphere(&mut model, r(2), &Pose::identity()).unwrap();
    let clip = block(&mut model, [r(6), r(6), r(3)], &pose(-3, -3, 0)).unwrap();
    let top = boolean(&mut model, stock, clip, BoolOp::Intersect, &mut dec)
        .unwrap()
        .body;
    let bottom = boolean(&mut model, stock, clip, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let faces = [top, bottom].map(|body| {
        model
            .faces_of_body(body)
            .into_iter()
            .find(|face| {
                matches!(
                    model.geom().surface(model.face(*face).surface),
                    catallax_shapes::Surface::Sphere(_)
                )
            })
            .unwrap()
    });
    let cap = model
        .faces_of_body(top)
        .into_iter()
        .find(|face| {
            matches!(
                model.geom().surface(model.face(*face).surface),
                catallax_shapes::Surface::Plane(_)
            )
        })
        .unwrap();
    let mut queries = crate::PreparedFaceQueries::new(&mut model);
    for (face, north) in faces.into_iter().zip([true, false]) {
        assert_eq!(
            queries.contains(face, &point(0, 0, 2), &mut dec).unwrap(),
            north
        );
        assert_eq!(
            queries.contains(face, &point(0, 0, -2), &mut dec).unwrap(),
            !north
        );
        assert!(queries.contains(face, &point(2, 0, 0), &mut dec).unwrap());
        assert!(!queries.contains(face, &point(0, 0, 0), &mut dec).unwrap());
    }
    assert!(queries.contains(cap, &point(0, 0, 0), &mut dec).unwrap());
    assert!(!queries.contains(cap, &point(3, 0, 0), &mut dec).unwrap());
    dec.set_deadline(Some(std::time::Instant::now()));
    assert_eq!(
        queries
            .contains(faces[0], &point(0, 0, 2), &mut dec)
            .unwrap_err()
            .code,
        "numbers.deadline",
    );
}

#[test]
fn primitive_containment_separates_exact_near_boundary_points() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(2), r(3), r(5)], &Pose::identity()).unwrap();
    for (p, want) in [
        (point(1, 1, 1), Some(Containment::Inside)),
        (point(2, 1, 1), None),
        (point(3, 1, 1), Some(Containment::Outside)),
        (point(2, 3, 5), None),
    ] {
        let answer = classify_point(&mut model, a, &p, &mut dec);
        if let Some(want) = want {
            assert_eq!(answer.proven(), Some(&want));
        } else {
            assert!(answer.is_coincident());
        }
    }
    let c = cylinder(&mut model, r(2), r(5), &Pose::identity()).unwrap();
    decide(&mut model, c, &point(0, 0, 1), true, &mut dec);
    decide(&mut model, c, &point(2, 2, 1), false, &mut dec);
    let rim = RPoint3::new(
        Rational::parse("6/5").unwrap(),
        Rational::parse("8/5").unwrap(),
        r(5),
    );
    assert!(classify_point(&mut model, c, &rim, &mut dec).is_coincident());
    let b = sphere(&mut model, r(5), &Pose::identity()).unwrap();
    let epsilon = Rational::parse("1e-30").unwrap();
    assert!(classify_point(&mut model, b, &point(3, 4, 0), &mut dec).is_coincident());
    decide(
        &mut model,
        b,
        &RPoint3::new(r(3), r(4), epsilon.clone()),
        false,
        &mut dec,
    );
    decide(
        &mut model,
        b,
        &RPoint3::new(r(3), r(4) - epsilon, r(0)),
        true,
        &mut dec,
    );
}

#[test]
fn crossing_blocks_have_exact_material_and_maximal_counts() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
    let b = block(&mut model, [r(2), r(2), r(2)], &pose(1, 1, 1)).unwrap();
    for (op, expected) in [
        (BoolOp::Unite, (1, 1, 12, 30, 0, 20)),
        (BoolOp::Intersect, (1, 1, 6, 12, 0, 8)),
        (BoolOp::Subtract, (1, 1, 9, 21, 0, 14)),
    ] {
        let body = operation(&mut model, a, b, op, expected, &mut dec);
        let mut seed = 0x6915_7392u64;
        for _ in 0..128 {
            let xyz: [Rational; 3] = std::array::from_fn(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                Rational::from_ratio((seed % 10) as i64 * 2 - 3, 4).unwrap()
            });
            let in_a = xyz.iter().all(|x| *x > r(0) && *x < r(2));
            let in_b = xyz.iter().all(|x| *x > r(1) && *x < r(3));
            let want = match op {
                BoolOp::Unite => in_a || in_b,
                BoolOp::Intersect => in_a && in_b,
                BoolOp::Subtract => in_a && !in_b,
            };
            let [x, y, z] = xyz;
            decide(&mut model, body, &RPoint3::new(x, y, z), want, &mut dec);
        }
    }
}

#[test]
fn face_contacts_regularise_and_coplanar_faces_merge() {
    for (offset, union) in [
        ((1, 0, 0), (1, 1, 6, 12, 0, 8)),
        ((1, 1, 0), (2, 2, 12, 24, 0, 16)),
        ((1, 1, 1), (2, 2, 12, 24, 0, 16)),
    ] {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let a = block(&mut model, [r(1), r(1), r(1)], &Pose::identity()).unwrap();
        let b = block(
            &mut model,
            [r(1), r(1), r(1)],
            &pose(offset.0, offset.1, offset.2),
        )
        .unwrap();
        operation(&mut model, a, b, BoolOp::Unite, union, &mut dec);
        operation(
            &mut model,
            a,
            b,
            BoolOp::Intersect,
            (0, 0, 0, 0, 0, 0),
            &mut dec,
        );
        operation(
            &mut model,
            a,
            b,
            BoolOp::Subtract,
            (1, 1, 6, 12, 0, 8),
            &mut dec,
        );
    }
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
    let b = block(&mut model, [r(2), r(2), r(2)], &pose(2, 1, 0)).unwrap();
    operation(
        &mut model,
        a,
        b,
        BoolOp::Unite,
        (1, 1, 10, 24, 0, 16),
        &mut dec,
    );
    operation(
        &mut model,
        a,
        b,
        BoolOp::Subtract,
        (1, 1, 6, 12, 0, 8),
        &mut dec,
    );
}

#[test]
fn through_cylinder_has_hole_rings_and_no_seam_edges() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
    let b = cylinder(&mut model, r(3), r(12), &pose(5, 5, -1)).unwrap();
    let graph = intersect(&mut model, a, b, &mut dec).unwrap().value;
    assert_eq!(graph.edges().len(), 2);
    assert!(graph.points().is_empty());
    for edge in graph.edges() {
        assert!(edge.start.is_none() && edge.end.is_none());
        assert_eq!(edge.faces_of_first.len(), 1);
        assert_eq!(edge.faces_of_second.len(), 1);
    }
    let cut = operation(
        &mut model,
        a,
        b,
        BoolOp::Subtract,
        (1, 1, 7, 14, 2, 8),
        &mut dec,
    );
    decide(&mut model, cut, &point(5, 5, 5), false, &mut dec);
    decide(&mut model, cut, &point(1, 1, 1), true, &mut dec);
    assert!(classify_point(&mut model, cut, &point(5, 8, 5), &mut dec).is_coincident());
    operation(
        &mut model,
        a,
        b,
        BoolOp::Unite,
        (1, 1, 10, 16, 4, 8),
        &mut dec,
    );
    operation(
        &mut model,
        a,
        b,
        BoolOp::Intersect,
        (1, 1, 3, 2, 2, 0),
        &mut dec,
    );
}

#[test]
fn isolated_tangent_face_probe_does_not_decide_a_piece_class() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = sphere(&mut model, r(5), &Pose::identity()).unwrap();
    let face = model.faces_of_body(a)[0];
    let probe = face_point(&mut model, face, &mut dec).unwrap().value;
    let displacement = probe
        .minus(&catallax_shapes::Point3::origin())
        .scaled(&catallax_numbers::Real::from_i64(2));
    let b = sphere(
        &mut model,
        r(5),
        &Pose::new(
            Rotation::identity(),
            displacement,
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let graph = intersect(&mut model, a, b, &mut dec).unwrap().value;
    assert!(graph.edges().is_empty());
    assert_eq!(graph.points().len(), 1);
    let im = imprint(&mut model, &graph, &mut dec).unwrap().value;
    let classes = classify(&mut model, &im, &mut dec).unwrap().value;
    for piece in im
        .pieces_of_first()
        .into_iter()
        .chain(im.pieces_of_second())
    {
        assert_eq!(classes.class_of(piece), Some(FaceClass::Outside));
    }
    operation(
        &mut model,
        a,
        b,
        BoolOp::Unite,
        (2, 2, 2, 0, 0, 0),
        &mut dec,
    );
    operation(
        &mut model,
        a,
        b,
        BoolOp::Intersect,
        (0, 0, 0, 0, 0, 0),
        &mut dec,
    );
}

#[test]
fn spherical_lens_and_void_keep_outward_sense() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = sphere(&mut model, r(5), &Pose::identity()).unwrap();
    let b = sphere(&mut model, r(5), &pose(6, 0, 0)).unwrap();
    for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
        let result = operation(&mut model, a, b, op, (1, 1, 2, 1, 1, 0), &mut dec);
        let want = op != BoolOp::Subtract;
        decide(&mut model, result, &point(3, 0, 0), want, &mut dec);
    }
    let inside = sphere(&mut model, r(3), &Pose::identity()).unwrap();
    let cut = operation(
        &mut model,
        a,
        inside,
        BoolOp::Subtract,
        (1, 2, 2, 0, 0, 0),
        &mut dec,
    );
    decide(&mut model, cut, &point(0, 0, 0), false, &mut dec);
    decide(&mut model, cut, &point(4, 0, 0), true, &mut dec);
    let tangent = sphere(&mut model, r(2), &pose(3, 0, 0)).unwrap();
    operation(
        &mut model,
        a,
        tangent,
        BoolOp::Subtract,
        (1, 2, 2, 0, 0, 0),
        &mut dec,
    );
}

fn sum_stats(records: &[crate::StageRecord]) -> DecisionStats {
    let mut sum = DecisionStats::default();
    for record in records {
        let s = &record.stats;
        sum.decisions += s.decisions;
        sum.filter += s.filter;
        sum.refined += s.refined;
        sum.exact += s.exact;
        sum.coincident_exact += s.coincident_exact;
        sum.coincident_separation += s.coincident_separation;
        sum.coincident_fact += s.coincident_fact;
        sum.not_exact += s.not_exact;
        sum.refused += s.refused;
        sum.work += s.work;
        sum.max_bits = sum.max_bits.max(s.max_bits);
    }
    sum
}
#[test]
fn stage_reports_preserve_a_reused_deciders_counters() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
    let b = block(&mut model, [r(2), r(2), r(2)], &pose(1, 1, 1)).unwrap();
    dec.enclose(
        &Real::from_i64(2).sqrt(),
        &Rational::parse("1e-100").unwrap(),
    )
    .unwrap();
    let before = dec.stats().clone();
    let answer = boolean(&mut model, a, b, BoolOp::Unite, &mut dec).unwrap();
    let after = dec.stats().clone();
    assert_eq!(
        answer
            .report
            .stages
            .iter()
            .map(|r| r.stage)
            .collect::<Vec<_>>(),
        [
            Stage::Intersect,
            Stage::Imprint,
            Stage::Classify,
            Stage::Keep,
            Stage::Sew
        ]
    );
    let sum = sum_stats(&answer.report.stages);
    macro_rules! preserved {
        ($($field:ident),+ $(,)?) => {$(assert_eq!(after.$field, before.$field + sum.$field, stringify!($field));)+};
    }
    preserved!(
        decisions,
        filter,
        refined,
        exact,
        coincident_exact,
        coincident_separation,
        coincident_fact,
        not_exact,
        refused,
        work
    );
    assert!(
        answer
            .report
            .stages
            .windows(2)
            .all(|w| w[0].stats.max_bits <= w[1].stats.max_bits)
    );
    assert!(answer.report.stages[0].stats.max_bits >= before.max_bits);
    assert_eq!(answer.report.stages[4].stats.max_bits, after.max_bits);
}

#[test]
fn same_body_booleans_obey_the_set_rules_without_changing_the_operand() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(2), r(3), r(5)], &Pose::identity()).unwrap();
    for op in [BoolOp::Unite, BoolOp::Subtract, BoolOp::Intersect] {
        let expected = if op == BoolOp::Subtract {
            (0, 0, 0, 0, 0, 0)
        } else {
            (1, 1, 6, 12, 0, 8)
        };
        let result = operation(&mut model, a, a, op, expected, &mut dec);
        decide(
            &mut model,
            result,
            &point(1, 1, 1),
            op != BoolOp::Subtract,
            &mut dec,
        );
        decide(&mut model, result, &point(3, 1, 1), false, &mut dec);
        decide(&mut model, a, &point(1, 1, 1), true, &mut dec);
        assert_eq!(counts(&model, a), (1, 1, 6, 12, 0, 8));
    }
}

#[test]
fn empty_operands_and_expired_deadlines_do_not_corrupt_bodies() {
    let mut model = Model::new();
    let mut dec = Decider::new(BudgetPolicy::default());
    let a = block(&mut model, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
    let e = empty(&mut model);
    for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
        let expected = if op == BoolOp::Intersect {
            (0, 0, 0, 0, 0, 0)
        } else {
            (1, 1, 6, 12, 0, 8)
        };
        operation(&mut model, a, e, op, expected, &mut dec);
    }
    let answer = boolean(&mut model, e, e, BoolOp::Unite, &mut dec).unwrap();
    assert_eq!(model.body(answer.body).kind, BodyKind::Empty);
    assert!(answer.report.stages.is_empty());
    let b = block(&mut model, [r(2), r(2), r(2)], &pose(1, 1, 1)).unwrap();
    dec.set_deadline(Some(std::time::Instant::now()));
    let outcome = boolean(&mut model, a, b, BoolOp::Unite, &mut dec);
    match outcome {
        Err(e) => assert_eq!(e.refusal.code, "numbers.deadline"),
        Ok(answer) => assert_eq!(counts(&model, answer.body), (1, 1, 12, 30, 0, 20)),
    };
    assert!(model.check(a).is_empty());
    assert!(model.check(b).is_empty());
}

#[test]
fn tangent_side_contacts_split_the_section_circles() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
    let b = cylinder(&mut model, r(5), r(12), &pose(5, 5, -1)).unwrap();
    let graph = intersect(&mut model, a, b, &mut dec).unwrap().value;
    assert_eq!(graph.edges().len(), 12);
    assert!(graph.points().is_empty());
    assert!(
        graph
            .edges()
            .iter()
            .all(|e| e.start.is_some() && e.end.is_some())
    );
    let body = operation(
        &mut model,
        a,
        b,
        BoolOp::Intersect,
        (1, 1, 3, 2, 2, 0),
        &mut dec,
    );
    decide(&mut model, body, &point(5, 5, 5), true, &mut dec);
    decide(&mut model, body, &point(1, 1, 5), false, &mut dec);
}

#[test]
fn connected_tangent_corners_refuse_as_a_self_touching_shell() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let a = block(&mut model, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
    let b = cylinder(&mut model, r(5), r(6), &pose(5, 5, 2)).unwrap();
    let error = boolean(&mut model, a, b, BoolOp::Subtract, &mut dec).unwrap_err();
    assert_eq!(error.stage, crate::Stage::Sew);
    assert_eq!(error.refusal.code, "boolean.self-touching-shell");
    assert!(model.check(a).is_empty());
    assert!(model.check(b).is_empty());
}

#[test]
fn skipped_carriers_keep_the_complete_section_order() {
    // A public edge-hunt block/cylinder intersection whose skipped carriers
    // cross the closed section curves. Its STEP bytes depend on which cut each
    // closed curve starts at.
    let q = |n: i64, d: i64| Rational::from_ratio(n, d).unwrap();
    let crossing = |model: &mut Model, dec: &mut Decider| {
        let a = block(model, [r(2), r(2), r(2)], &pose(-1, -1, -1)).unwrap();
        let rotation = Rotation::from_quaternion(q(1, 2), q(1, 3), q(1, 5), q(1, 7)).unwrap();
        let translation = RVector3::new(q(-33491, 37178), q(3161, 1923), q(-34021, 92945));
        let placement = Pose::new(rotation, translation.into(), dec).unwrap();
        (a, cylinder(model, r(2), r(4), &placement).unwrap())
    };
    // Coplanar top and bottom faces overlap, so a carrier could keep a piece
    // inside the common region without a meeting generating pair.
    let overlap = |model: &mut Model, _: &mut Decider| {
        let a = block(model, [r(6), r(4), r(2)], &Pose::identity()).unwrap();
        (
            a,
            block(model, [r(4), r(6), r(2)], &pose(3, -1, 0)).unwrap(),
        )
    };
    type Scene<'a> = &'a dyn Fn(&mut Model, &mut Decider) -> (BodyId, BodyId);
    let scenes: [Scene; 2] = [&crossing, &overlap];
    for scene in scenes {
        let section = |prune: bool, carrier_exclusions: bool| {
            let (answer, work) = catallax_numbers::with_operation_budget(100_000_000, || {
                let mut model = Model::new();
                let mut dec = Decider::default();
                let (a, b) = scene(&mut model, &mut dec);
                let graph =
                    crate::graph::section(&mut model, a, b, prune, carrier_exclusions, &mut dec)
                        .unwrap()
                        .value;
                let at = |point: Option<catallax_shapes::PointId>| {
                    point.map(|point| {
                        let p = model.geom().position(point);
                        [&p.x, &p.y, &p.z].map(|c| format!("{:.6}", c.to_f64()))
                    })
                };
                let edges: Vec<_> = graph
                    .edges()
                    .iter()
                    .map(|edge| {
                        (
                            model.geom().curve(edge.curve).clone(),
                            at(edge.start),
                            at(edge.end),
                            edge.faces_of_first.clone(),
                            edge.faces_of_second.clone(),
                        )
                    })
                    .collect();
                let points: Vec<_> = graph.points().iter().map(|&p| at(Some(p))).collect();
                (edges, points, graph.marks.len())
            });
            assert!(!work.exhausted, "{work:?}");
            (answer, work.used)
        };
        let complete = section(false, false);
        let baseline = section(true, false);
        let bounded = section(true, true);
        let same_graph = |a: &(Vec<_>, Vec<_>, usize), b: &(Vec<_>, Vec<_>, usize)| {
            assert_eq!(a.0.len(), b.0.len());
            assert_eq!(a.1, b.1);
            assert_eq!(a.2, b.2);
            let mut geom = catallax_shapes::GeomStore::new();
            let mut dec = Decider::default();
            for ((ca, a0, a1, af, ag), (cb, b0, b1, bf, bg)) in a.0.iter().zip(&b.0) {
                assert_eq!((a0, a1, af, ag), (b0, b1, bf, bg));
                let ca = geom.add_curve(catallax_shapes::Curve::clone(ca));
                let cb = geom.add_curve(catallax_shapes::Curve::clone(cb));
                assert!(catallax_intersect::same_curve(&geom, ca, cb, &mut dec).is_coincident());
            }
        };
        same_graph(&baseline.0, &complete.0);
        same_graph(&bounded.0, &baseline.0);
        assert!(
            bounded.1 <= baseline.1,
            "bounded {}, baseline {}",
            bounded.1,
            baseline.1
        );
    }
}

/// A copy of `source` whose first edge, and the two vertices it ends at,
/// store a tolerance. Everything else of the copy is exact.
fn with_a_tolerant_edge(model: &mut Model, source: BodyId) -> BodyId {
    use catallax_topology::{Body, Coedge, Edge, Face, Loop, Lump, Shell, Tolerance, Vertex};
    use std::collections::BTreeMap;
    let stored = Tolerance::new(Rational::parse("1/1000000").unwrap()).unwrap();
    let edge = model.edges_of_body(source)[0];
    let ends = [model.edge(edge).start, model.edge(edge).end];
    let mut vertices = BTreeMap::new();
    for v in model.vertices_of_body(source) {
        let copy = Vertex {
            tolerance: ends.contains(&Some(v)).then(|| stored.clone()),
            ..model.vertex(v).clone()
        };
        vertices.insert(v, model.add_vertex(copy));
    }
    let mut edges = BTreeMap::new();
    for e in model.edges_of_body(source) {
        let old = model.edge(e).clone();
        let copy = Edge {
            start: old.start.map(|v| vertices[&v]),
            end: old.end.map(|v| vertices[&v]),
            tolerance: (e == edge).then(|| stored.clone()),
            ..old
        };
        edges.insert(e, model.add_edge(copy));
    }
    let mut faces = Vec::new();
    for f in model.faces_of_body(source) {
        let old = model.face(f).clone();
        let mut loops = Vec::new();
        for l in &old.loops {
            let Loop::Coedges(uses) = model.loop_(*l).clone() else {
                unreachable!("a block face has loops of coedges")
            };
            let uses = uses
                .into_iter()
                .map(|c| {
                    let c = model.coedge(c).clone();
                    model.add_coedge(Coedge {
                        edge: edges[&c.edge],
                        sense: c.sense,
                    })
                })
                .collect();
            loops.push(model.add_loop(Loop::Coedges(uses)));
        }
        faces.push(model.add_face(Face { loops, ..old }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: Default::default(),
    })
}

#[test]
fn a_tolerant_entity_never_reaches_the_exact_intersector() {
    // Checkpoint 2 decisions, section 4.5. One tolerant edge in an otherwise
    // exact block is enough: the boolean refuses by name before it makes a
    // decision, whichever operand holds it and whatever the operation.
    let mut model = Model::new();
    let exact = block(&mut model, [r(4), r(4), r(4)], &pose(0, 0, 0)).unwrap();
    let tool = cylinder(&mut model, r(1), r(10), &pose(2, 2, -3)).unwrap();
    let mixed = with_a_tolerant_edge(&mut model, exact);
    assert_eq!(model.check(mixed), Vec::new());
    assert_eq!(model.tolerant_entity(exact), None);
    assert!(
        model
            .tolerant_entity(mixed)
            .is_some_and(|entity| entity.starts_with("edge "))
    );
    let mut dec = Decider::default();
    for op in [BoolOp::Subtract, BoolOp::Unite, BoolOp::Intersect] {
        for (a, b) in [(mixed, tool), (tool, mixed), (mixed, mixed)] {
            let (refused, work) = catallax_numbers::with_operation_budget(1_000_000, || {
                boolean(&mut model, a, b, op, &mut dec).unwrap_err()
            });
            assert_eq!(refused.refusal.code, crate::TOLERANT_PARTICIPANT);
            assert_eq!(refused.stage, crate::Stage::Intersect);
            // Nothing was decided and nothing was charged.
            assert_eq!((work.used, work.refunded), (0, 0), "{op:?}");
            assert_eq!(dec.stats().decisions, 0);
        }
    }
    // The same cut of the exact block is answered as before.
    assert!(boolean(&mut model, exact, tool, BoolOp::Subtract, &mut dec).is_ok());
}

mod tolerance_publication {
    use super::*;

    #[test]
    fn an_empty_kind_does_not_bypass_publication_topology_checks() {
        let mut model = Model::new();
        let solid = block(&mut model, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
        let mut malformed = model.body(solid).clone();
        malformed.kind = BodyKind::Empty;
        let malformed = model.add_body(malformed);
        let error =
            crate::make_maximal(&mut model, malformed, &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, catallax_numbers::GEOMETRIC_GAP_UNPROVEN);
    }

    #[test]
    fn exact_point_merging_does_not_chain_through_subnanometre_gaps() {
        use catallax_shapes::{Point3, PointDef};
        let mut model = Model::new();
        let points = ["0", "0.0000006", "0.0000012"].map(|x| {
            model.geom_mut().add_point(PointDef::Explicit(Point3::new(
                Real::from_rational(Rational::parse(x).unwrap()),
                Real::zero(),
                Real::zero(),
            )))
        });
        let mut pool = Vec::with_capacity(points.len());
        let mut dec = Decider::default();
        for point in points {
            let represented = crate::geometry::canonical_point(
                &mut model,
                &mut pool,
                point,
                &mut dec,
                "boolean.imprint-ambiguous",
            )
            .unwrap();
            assert_eq!(represented, point);
        }
        assert_eq!(pool, points);
    }
}
