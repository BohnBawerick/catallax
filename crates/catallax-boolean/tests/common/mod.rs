#![allow(dead_code, unused_imports, unused_macros)]
pub use catallax_boolean::*;
pub use catallax_intersect::{CurveForm, curve_form};
pub use catallax_numbers::{
    Answer, Answered, BoundKind, BudgetPolicy, Decider, DecisionStats, Fact, FactKind, FactOrigin,
    Mark, MarkLog, Outcome, Rational, Real, Refusal, Sign,
};
pub use catallax_shapes::*;
pub use catallax_topology::*;
pub use std::collections::{BTreeMap, BTreeSet};
mod cases;
mod oracle;
pub use cases::*;
pub use oracle::*;
mod witness;
pub use witness::*;
mod chains;
pub use chains::*;

// Passing cases run by default; quarantined cases retain their observed reason.
macro_rules! blind_test {
    ($(#[$attr:meta])* $name:ident, $body:block)=>{
        #[test]
        $(#[$attr])*
        fn $name() $body
    };
}
pub(crate) use blind_test;

pub fn exact<T: std::fmt::Debug>(a: Answered<T>) -> T {
    let a = a.unwrap();
    assert!(a.marks.is_empty(), "unexpected marks: {:?}", a.marks);
    a.value
}
pub fn query(m: &mut Model, b: BodyId, p: &RPoint3, want: Location, dec: &mut Decider) {
    let got = classify_point(m, b, p, dec);
    match (want, &got) {
        (Location::Inside, Outcome::Proven(Containment::Inside))
        | (Location::Outside, Outcome::Proven(Containment::Outside))
        | (Location::On, Outcome::ProvenCoincident(_)) => {}
        _ => panic!("point {p:?}: expected {want:?}, got {got:?}"),
    }
}
pub fn counts(m: &Model, b: BodyId) -> Counts {
    let c = m.counts(b);
    let mut curves = [0; 4];
    let mut dec = Decider::default();
    for e in m.edges_of_body(b) {
        let kind = match curve_form(m.geom(), m.edge(e).curve, &mut dec).unwrap() {
            CurveForm::Line { .. } => 0,
            CurveForm::Circle { .. } => 1,
            CurveForm::Ellipse { .. } => 2,
            CurveForm::Branch => 3,
            _ => panic!("new curve form"),
        };
        curves[kind] += 1;
    }
    let surfaces = [
        SurfaceKind::Plane,
        SurfaceKind::Cylinder,
        SurfaceKind::Sphere,
    ]
    .map(|k| c.faces_by_kind.get(&k).copied().unwrap_or(0));
    assert_eq!(
        surfaces.iter().sum::<usize>(),
        c.faces,
        "unexpected surface kind"
    );
    Counts {
        lumps: c.lumps,
        shells: c.shells,
        faces: c.faces,
        edges: c.edges,
        rings: c.ring_edges,
        vertices: c.vertices,
        surfaces,
        curves,
    }
}
pub fn check_counts(m: &Model, b: BodyId, want: &Counts) {
    assert!(m.check(b).is_empty(), "topology faults: {:?}", m.check(b));
    assert_eq!(&counts(m, b), want);
    assert_eq!(
        m.body(b).kind,
        if want.lumps == 0 {
            BodyKind::Empty
        } else {
            BodyKind::Solid
        }
    );
}

pub fn topology_snapshot(m: &Model, bs: &[BodyId]) -> Vec<String> {
    // Debug here serialises immutable public entity values, not an expected
    // textual format. It compares before/after, or identical fresh executions.
    let mut out = Vec::new();
    for &b in bs {
        out.push(format!("{b:?} {:?}", m.body(b)));
        for &l in &m.body(b).lumps {
            out.push(format!("{l:?} {:?}", m.lump(l)));
            for &s in &m.lump(l).shells {
                out.push(format!("{s:?} {:?}", m.shell(s)));
            }
        }
        for f in m.faces_of_body(b) {
            out.push(format!(
                "{f:?} {:?} {:?}",
                m.face(f),
                m.geom().surface(m.face(f).surface)
            ));
            for &l in &m.face(f).loops {
                out.push(format!("{l:?} {:?}", m.loop_(l)));
                if let Loop::Coedges(cs) = m.loop_(l) {
                    for &c in cs {
                        out.push(format!("{c:?} {:?}", m.coedge(c)));
                    }
                }
            }
        }
        for e in m.edges_of_body(b) {
            out.push(format!(
                "{e:?} {:?} {:?}",
                m.edge(e),
                m.geom().curve(m.edge(e).curve)
            ));
        }
        for v in m.vertices_of_body(b) {
            out.push(format!(
                "{v:?} {:?} {:?}",
                m.vertex(v),
                m.geom().point(m.vertex(v).point)
            ));
        }
    }
    out
}

pub fn topology_rules(m: &Model, b: BodyId) {
    assert!(m.check(b).is_empty(), "{:?}", m.check(b));
    let mut owned_faces = BTreeSet::new();
    let mut owned_edges = BTreeSet::new();
    let mut owned_vertices = BTreeSet::new();
    for &l in &m.body(b).lumps {
        let mut lump_faces = BTreeSet::new();
        let mut lump_edges = BTreeSet::new();
        let mut lump_vertices = BTreeSet::new();
        for &s in &m.lump(l).shells {
            let faces = &m.shell(s).faces;
            let mut uses: BTreeMap<EdgeId, Vec<(FaceId, Sense)>> = BTreeMap::new();
            for &f in faces {
                lump_faces.insert(f);
                let mut loop_vertices = BTreeSet::new();
                for &lp in &m.face(f).loops {
                    let mut in_loop = BTreeSet::new();
                    if let Loop::Coedges(cs) = m.loop_(lp) {
                        for &c in cs {
                            let c = m.coedge(c);
                            let e = m.edge(c.edge);
                            uses.entry(c.edge).or_default().push((f, c.sense));
                            lump_edges.insert(c.edge);
                            for v in [e.start, e.end].into_iter().flatten() {
                                in_loop.insert(v);
                                lump_vertices.insert(v);
                            }
                        }
                    }
                    // A node may occur twice in ONE loop, never in two loops.
                    assert!(
                        loop_vertices.is_disjoint(&in_loop),
                        "face has loops sharing a vertex"
                    );
                    loop_vertices.extend(in_loop);
                }
            }
            for (&e, u) in &uses {
                assert_eq!(u.len(), 2, "edge {e:?} must have two uses in this shell");
                assert_ne!(u[0].0, u[1].0, "slit was not removed");
                assert_ne!(u[0].1, u[1].1);
                let mut d = Decider::default();
                assert!(
                    matches!(
                        m.geom().same_surface(
                            m.face(u[0].0).surface,
                            m.face(u[1].0).surface,
                            &mut d
                        ),
                        Outcome::Proven(_)
                    ),
                    "non-maximal coplanar/co-cylindrical face split"
                );
            }
            let mut incident: BTreeMap<VertexId, Vec<EdgeId>> = BTreeMap::new();
            for &e in uses.keys() {
                for v in [m.edge(e).start, m.edge(e).end]
                    .into_iter()
                    .flatten()
                    .collect::<BTreeSet<_>>()
                {
                    incident.entry(v).or_default().push(e);
                }
            }
            for edges in incident.values().filter(|es| es.len() == 2) {
                let face_pair = |e: &EdgeId| uses[e].iter().map(|x| x.0).collect::<BTreeSet<_>>();
                if face_pair(&edges[0]) == face_pair(&edges[1]) {
                    let mut d = Decider::default();
                    assert!(
                        matches!(
                            catallax_intersect::same_curve(
                                m.geom(),
                                m.edge(edges[0]).curve,
                                m.edge(edges[1]).curve,
                                &mut d
                            ),
                            Outcome::Proven(_)
                        ),
                        "redundant degree-two vertex"
                    );
                }
            }
        }
        assert!(owned_faces.is_disjoint(&lump_faces));
        assert!(owned_edges.is_disjoint(&lump_edges));
        assert!(owned_vertices.is_disjoint(&lump_vertices));
        owned_faces.extend(lump_faces);
        owned_edges.extend(lump_edges);
        owned_vertices.extend(lump_vertices);
    }
}

// Counter deltas preserve caller totals. Precision is the running end maximum,
// not a stage-local peak, as specified by boolean-stage-statistics.
pub fn stats_delta(before: &DecisionStats, after: &DecisionStats) -> DecisionStats {
    let delta = |start: u64, end: u64| end.checked_sub(start).expect("caller counter decreased");
    assert!(
        after.max_bits >= before.max_bits,
        "caller precision maximum decreased"
    );
    DecisionStats {
        decisions: delta(before.decisions, after.decisions),
        filter: delta(before.filter, after.filter),
        refined: delta(before.refined, after.refined),
        exact: delta(before.exact, after.exact),
        coincident_exact: delta(before.coincident_exact, after.coincident_exact),
        coincident_separation: delta(before.coincident_separation, after.coincident_separation),
        coincident_fact: delta(before.coincident_fact, after.coincident_fact),
        not_exact: delta(before.not_exact, after.not_exact),
        refused: delta(before.refused, after.refused),
        work: delta(before.work, after.work),
        max_bits: after.max_bits,
    }
}

pub fn sum_stats(records: &[StageRecord]) -> DecisionStats {
    let mut s = DecisionStats::default();
    for x in records {
        let x = &x.stats;
        s.decisions += x.decisions;
        s.filter += x.filter;
        s.refined += x.refined;
        s.exact += x.exact;
        s.coincident_exact += x.coincident_exact;
        s.coincident_separation += x.coincident_separation;
        s.coincident_fact += x.coincident_fact;
        s.not_exact += x.not_exact;
        s.refused += x.refused;
        s.work += x.work;
        s.max_bits = s.max_bits.max(x.max_bits);
    }
    s
}
pub fn check_report(answer: &BooleanAnswer, before: &DecisionStats, after: &DecisionStats) {
    assert_eq!(
        answer
            .report
            .stages
            .iter()
            .map(|s| s.stage)
            .collect::<Vec<_>>(),
        [
            Stage::Intersect,
            Stage::Imprint,
            Stage::Classify,
            Stage::Keep,
            Stage::Sew
        ]
    );
    assert_eq!(sum_stats(&answer.report.stages), stats_delta(before, after));
    let mut previous_max = before.max_bits;
    for record in &answer.report.stages {
        assert!(
            record.stats.max_bits >= previous_max,
            "stage precision maximum decreased"
        );
        previous_max = record.stats.max_bits;
    }
    assert_eq!(previous_max, after.max_bits);
    assert_eq!(
        answer.report.stages.iter().map(|s| s.marks).sum::<usize>(),
        answer.marks.len()
    );
    assert_eq!(
        answer.report.stages[3].stats,
        DecisionStats {
            max_bits: answer.report.stages[2].stats.max_bits,
            ..DecisionStats::default()
        }
    );
    assert_eq!(answer.report.stages[3].marks, 0);
}
pub fn same_real_point(got: &Point3, want: &RPoint3) {
    let mut d = Decider::default();
    for (x, y) in [(&got.x, &want.x), (&got.y, &want.y), (&got.z, &want.z)] {
        assert!(
            matches!(
                d.compare(x, &Real::from(y.clone())),
                Outcome::ProvenCoincident(_)
            ),
            "coordinate mismatch"
        );
    }
}

pub fn check_faces(
    m: &mut Model,
    b: BodyId,
    a: BodyId,
    other: BodyId,
    pa: &Primitive,
    pb: &Primitive,
    op: BoolOp,
) {
    let source: Vec<_> = m
        .faces_of_body(a)
        .into_iter()
        .chain(m.faces_of_body(other))
        .map(|f| m.face(f).surface)
        .collect();
    let mut dec = Decider::default();
    for f in m.faces_of_body(b) {
        let surface = m.face(f).surface;
        assert!(
            source.iter().any(|s| matches!(
                m.geom().same_surface(surface, *s, &mut dec),
                Outcome::ProvenCoincident(_)
            )),
            "result introduced a new supporting surface"
        );
        let point = exact(face_point(m, f, &mut dec))
            .as_rational()
            .expect("rational Boolean fixture");
        let (implicit, _) = implicit_normal(m.geom().surface(surface), &point);
        assert_eq!(implicit, r(0), "face_point is off its supporting surface");
        let point = orientation_witness(m, b, f, [(a, pa), (other, pb)], point);
        let (_, normal) = implicit_normal(m.geom().surface(surface), &point);
        let normal = if m.face(f).sense == Sense::Same {
            normal
        } else {
            normal.scaled(&r(-1))
        };
        let step = safe_step(pa, pb, &point, &normal);
        let outside = point.plus(&normal.scaled(&step));
        let inside = point.plus(&normal.scaled(&(-step)));
        assert_eq!(
            regularized_step(op, pa, pb, &outside),
            Some(Location::Outside),
            "outward normal at {point:?}"
        );
        assert_eq!(
            regularized_step(op, pa, pb, &inside),
            Some(Location::Inside),
            "inward normal at {point:?}"
        );
    }
}

pub fn run_case(
    c: &Case,
    op: BoolOp,
    reverse: bool,
    expected: Option<&Expected>,
    points: usize,
) -> (Counts, MarkLog, DecisionStats, Vec<String>) {
    let (pa, pb) = if reverse { (&c.b, &c.a) } else { (&c.a, &c.b) };
    let mut m = Model::new();
    let a = pa.build(&mut m);
    let b = pb.build(&mut m);
    let before = topology_snapshot(&m, &[a, b]);
    let mut dec = Decider::default();
    let result = boolean(&mut m, a, b, op, &mut dec);
    assert_eq!(before, topology_snapshot(&m, &[a, b]), "operands changed");
    if let Some(Expected::Refused(code, stage)) = expected {
        let refusal = result.unwrap_err();
        assert_eq!(refusal.refusal.code, *code);
        assert_eq!(refusal.stage, *stage);
        return (
            Counts::empty(),
            MarkLog::new(),
            dec.stats().clone(),
            Vec::new(),
        );
    }
    let answer =
        result.unwrap_or_else(|e| panic!("case {} {op:?}, reverse={reverse}: {e:?}", c.id));
    assert_ne!(answer.body, a, "boolean must add a new body");
    assert_ne!(answer.body, b, "boolean must add a new body");
    check_report(&answer, &DecisionStats::default(), dec.stats());
    assert!(answer.marks.is_empty(), "case {} unexpectedly marked", c.id);
    assert_eq!(m.body(answer.body).marks, answer.marks);
    if let Some(Expected::Counts(want)) = expected {
        check_counts(&m, answer.body, want);
    }
    topology_rules(&m, answer.body);
    let stats = dec.stats().clone();
    check_faces(&mut m, answer.body, a, b, pa, pb, op);
    if !reverse {
        for (p, want) in witnesses(c.id, op) {
            query(&mut m, answer.body, &p, want, &mut dec);
        }
    }
    let mut rng = Rng(0x5b001ea5);
    let mut checked = 0;
    while checked < points {
        let point = rng.point();
        if let Some(want) = set_rule(op, pa.locate(&point), pb.locate(&point)) {
            query(&mut m, answer.body, &point, want, &mut dec);
            checked += 1;
        }
    }
    (
        counts(&m, answer.body),
        answer.marks,
        stats,
        topology_snapshot(&m, &[answer.body]),
    )
}
