//! Research only (R10, docs/R10-BODY-CHECKER-AUDIT.md). Read-only: builds each
//! recipe as `catallax-verify` does, then measures, in deterministic work units,
//! what a body checker that never reads the kernel's recorded facts would have
//! to recompute for the finished body. Kernel behaviour is unchanged.
//!
//! The "blind" copy re-adds the same geometry definitions, in the same order,
//! and the result body's topology to a fresh model. Its fact book then holds
//! only the incidences the store records from a definition itself (an
//! intersection curve on its two surfaces, a meeting point on what it meets),
//! none that a builder, the Boolean or a proof recorded.
//!
//! Cases run side by side on `--workers` threads; each case has its own
//! allowance and its own models, so the report is the same at any count.

use catallax_numbers::{
    Decider, FactKind, FactOrigin, Outcome, Rational, Real, Refusal, with_operation_budget,
};
use catallax_shapes::{
    Arena, Curve, CurveId, GeomKey, GeomStore, Line, Point3, PointDef, Surface, SurfaceId, Vector3,
};
use catallax_topology::{
    Body, BodyId, Coedge, Edge, Face, FaceId, Loop, Lump, Model, Sense, Shell,
};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

const USAGE: &str = "usage: r10_audit [--workers N] RECIPE...";
/// The replay allowance, for the build and for each property's pass.
const ALLOWANCE: u64 = 100_000_000;
/// Each surface pair of property 2 gets the application allowance.
const PAIR_ALLOWANCE: u64 = 20_000_000;

fn main() {
    let mut args = std::env::args().skip(1);
    let mut workers = 1;
    let mut paths = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--workers" {
            workers = args.next().and_then(|n| n.parse().ok()).expect(USAGE);
        } else {
            paths.push(arg);
        }
    }
    assert!(!paths.is_empty() && workers > 0, "{USAGE}");
    let next = AtomicUsize::new(0);
    let reports = Mutex::new(vec![Value::Null; paths.len()]);
    std::thread::scope(|scope| {
        for _ in 0..workers.min(paths.len()) {
            scope.spawn(|| {
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    let Some(path) = paths.get(index) else { break };
                    let report = case(path);
                    reports.lock().unwrap()[index] = report;
                }
            });
        }
    });
    let reports = reports.into_inner().unwrap();
    let summary = summary(&reports);
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "format": "catallax-r10-audit/1",
            "allowance": ALLOWANCE,
            "pair_allowance": PAIR_ALLOWANCE,
            "summary": summary,
            "planted_faults": planted(),
            "cases": reports,
        }))
        .unwrap()
    );
}

fn case(path: &str) -> Value {
    let name = std::path::Path::new(path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .into_owned();
    let text = std::fs::read_to_string(path).expect("readable recipe");
    let recipe = catallax_recipe::Recipe::parse(&text).expect("valid recipe");
    let mut dec = Decider::default();
    let (built, work) =
        with_operation_budget(ALLOWANCE, || catallax_features::evaluate(&recipe, &mut dec));
    let built = match built {
        Ok(built) if !work.exhausted => built,
        Ok(_) => {
            return json!({"case": name, "build": {"work": work.used, "refused": "numbers.operation-budget"}});
        }
        Err(error) => {
            let code = match &error {
                catallax_features::BuildError::Geometry { refusal, .. } => refusal.code.to_string(),
                other => format!("{other:?}"),
            };
            return json!({"case": name, "build": {"work": work.used, "refused": code}});
        }
    };
    let stats = dec.stats();
    let model = &built.model;
    let body = built.result;
    let Some(&first) = model.faces_of_body(body).first() else {
        return json!({"case": name, "build": {"work": work.used, "empty_body": true}});
    };
    let surface = model.geom().surface(model.face(first).surface).clone();
    let mut blind = blind_copy(model, body, surface.clone());
    let checks = check_questions(model, body, &blind.model);
    json!({
        "case": name,
        "build": {
            "work": work.used,
            "decisions": stats.decisions,
            "zero_decisions": stats.coincident_exact + stats.coincident_separation + stats.coincident_fact,
            "zero_decisions_settled_by_a_recorded_fact": stats.coincident_fact,
        },
        "body": census(model, body),
        "fact_book": fact_book(model),
        "topology_check": checks,
        "blind_topology_check": blind_check(&blind.model, blind.body),
        "edges_on_faces": edges_on_faces(&mut blind),
        "faces_meet_only_at_edges": {
            "carriers": faces_meet(&mut blind),
            "arrangement": arrangement(&mut blind_copy(model, body, surface.clone())),
            "across_shells": across_shells(&mut blind_copy(model, body, surface.clone())),
        },
        "inside_and_outside": inside_outside(&mut blind),
    })
}

/// The same geometry and the result body's topology, without recorded facts.
struct Blind {
    model: Model,
    body: BodyId,
}

// Handles are positions in an arena. A scratch arena of the same length issues
// the handle of every stored object, so definitions are re-added in order and
// every handle inside a definition keeps its meaning.
fn handles<T: Clone>(sample: T, count: usize) -> Vec<catallax_shapes::Handle<T>> {
    let mut arena = Arena::new();
    for _ in 0..count {
        arena.push(sample.clone());
    }
    arena.handles()
}

/// A fault planted in a copy, for the fault-injection rows.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Plant {
    Nothing,
    /// The first edge's curve becomes a line off every face.
    EdgeOffItsFaces,
    /// The first face's sense is reversed, its loops left as they were.
    FaceTurnedInsideOut,
}

fn blind_copy(model: &Model, body: BodyId, surface: Surface) -> Blind {
    planted_copy(model, body, surface, Plant::Nothing)
}

fn planted_copy(model: &Model, body: BodyId, surface: Surface, plant: Plant) -> Blind {
    let geom = model.geom();
    let line = Line::new(
        Point3::origin(),
        Vector3::new(Real::from_i64(1), Real::zero(), Real::zero()),
        &mut Decider::default(),
    )
    .expect("a unit direction");
    let point = PointDef::Explicit(Point3::origin());
    let mut copy = Model::new();
    let store = copy.geom_mut();
    for id in handles(surface, geom.surface_count()) {
        store.add_surface(geom.surface(id).clone());
    }
    for id in handles(Curve::Line(line), geom.curve_count()) {
        store.add_curve(geom.curve(id).clone());
    }
    for id in handles(point, geom.point_count()) {
        store.add_point(geom.point(id).clone());
    }
    let off = (plant == Plant::EdgeOffItsFaces).then(|| {
        let far = Point3::new(
            Real::from_i64(100),
            Real::from_i64(100),
            Real::from_i64(100),
        );
        let direction = Vector3::new(Real::from_i64(1), Real::zero(), Real::zero());
        let line = Line::new(far, direction, &mut Decider::default()).expect("a unit direction");
        store.add_curve(Curve::Line(line))
    });
    let mut vertex_map = BTreeMap::new();
    for v in model.vertices_of_body(body) {
        vertex_map.insert(v, copy.add_vertex(model.vertex(v).clone()));
    }
    let mut edge_map = BTreeMap::new();
    for (index, e) in model.edges_of_body(body).into_iter().enumerate() {
        let edge = model.edge(e);
        let copied = Edge {
            curve: off.filter(|_| index == 0).unwrap_or(edge.curve),
            sense: edge.sense,
            start: edge.start.map(|v| vertex_map[&v]),
            end: edge.end.map(|v| vertex_map[&v]),
        };
        edge_map.insert(e, copy.add_edge(copied));
    }
    let mut face_map = BTreeMap::new();
    for (index, f) in model.faces_of_body(body).into_iter().enumerate() {
        let face = model.face(f);
        let flip = plant == Plant::FaceTurnedInsideOut && index == 0;
        let loops = face
            .loops
            .iter()
            .map(|l| {
                let copied = match model.loop_(*l) {
                    Loop::Vertex(v) => Loop::Vertex(vertex_map[v]),
                    Loop::Coedges(cs) => Loop::Coedges(
                        cs.iter()
                            .map(|c| {
                                let coedge = model.coedge(*c);
                                copy.add_coedge(Coedge {
                                    edge: edge_map[&coedge.edge],
                                    sense: coedge.sense,
                                })
                            })
                            .collect(),
                    ),
                };
                copy.add_loop(copied)
            })
            .collect();
        let copied = Face {
            surface: face.surface,
            sense: match (flip, face.sense) {
                (false, sense) => sense,
                (true, Sense::Same) => Sense::Reversed,
                (true, Sense::Reversed) => Sense::Same,
            },
            loops,
        };
        face_map.insert(f, copy.add_face(copied));
    }
    let original = model.body(body);
    let lumps = original
        .lumps
        .iter()
        .map(|l| {
            let shells = model
                .lump(*l)
                .shells
                .iter()
                .map(|s| {
                    let faces = model.shell(*s).faces.iter().map(|f| face_map[f]).collect();
                    copy.add_shell(Shell { faces })
                })
                .collect();
            copy.add_lump(Lump { shells })
        })
        .collect();
    let body = copy.add_body(Body {
        kind: original.kind,
        lumps,
        marks: original.marks.clone(),
    });
    Blind { model: copy, body }
}

fn census(model: &Model, body: BodyId) -> Value {
    let counts = model.counts(body);
    let coedges: usize = model
        .faces_of_body(body)
        .iter()
        .flat_map(|f| model.face(*f).loops.iter())
        .map(|l| match model.loop_(*l) {
            Loop::Coedges(cs) => cs.len(),
            Loop::Vertex(_) => 0,
        })
        .sum();
    json!({
        "lumps": counts.lumps,
        "shells": counts.shells,
        "faces": counts.faces,
        "loops": counts.loops,
        "coedges": coedges,
        "edges": counts.edges,
        "ring_edges": counts.ring_edges,
        "vertices": counts.vertices,
        "faces_by_kind": counts.faces_by_kind.iter().map(|(k, n)| (format!("{k:?}"), *n)).collect::<BTreeMap<_, _>>(),
        "edges_by_kind": counts.edges_by_kind.iter().map(|(k, n)| (format!("{k:?}"), *n)).collect::<BTreeMap<_, _>>(),
    })
}

fn origin_name(origin: FactOrigin) -> String {
    match origin {
        FactOrigin::Construction(name) => name.to_string(),
        FactOrigin::Proof(proof) => format!("proof: {proof:?}"),
    }
}

fn fact_book(model: &Model) -> Value {
    let mut by_kind = BTreeMap::<String, BTreeMap<String, usize>>::new();
    for fact in model.geom().facts().iter() {
        *by_kind
            .entry(format!("{:?}", fact.kind))
            .or_default()
            .entry(origin_name(fact.origin))
            .or_default() += 1;
    }
    json!({"facts": model.geom().facts().len(), "by_kind_and_origin": by_kind})
}

/// Every `Incident` question the topology check asks of the body, in its order.
fn questions(model: &Model, body: BodyId) -> Vec<(GeomKey, GeomKey)> {
    let mut asked = Vec::new();
    let mut seen = BTreeSet::new();
    for f in model.faces_of_body(body) {
        let face = model.face(f);
        for l in &face.loops {
            match model.loop_(*l) {
                Loop::Vertex(v) => asked.push((
                    GeomKey::Point(model.vertex(*v).point),
                    GeomKey::Surface(face.surface),
                )),
                Loop::Coedges(cs) => {
                    for c in cs {
                        let e = model.coedge(*c).edge;
                        let edge = model.edge(e);
                        asked.push((GeomKey::Curve(edge.curve), GeomKey::Surface(face.surface)));
                        if seen.insert(e) {
                            for v in [edge.start, edge.end].into_iter().flatten() {
                                asked.push((
                                    GeomKey::Point(model.vertex(v).point),
                                    GeomKey::Curve(edge.curve),
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    asked
}

fn check_questions(model: &Model, body: BodyId, blind: &Model) -> Value {
    let mut direct = BTreeMap::new();
    for fact in model.geom().facts().iter() {
        if fact.kind == FactKind::Incident && fact.at.is_none() {
            direct.entry((fact.a, fact.b)).or_insert(fact.origin);
        }
    }
    let asked = questions(model, body);
    let mut answered_by = BTreeMap::<String, usize>::new();
    let mut by_definition = 0;
    for (a, b) in &asked {
        if blind.geom().facts().holds(FactKind::Incident, *a, *b) {
            by_definition += 1;
        }
        let how = match direct.get(&(*a, *b)) {
            Some(origin) => origin_name(*origin),
            None if model.geom().facts().holds(FactKind::Incident, *a, *b) => {
                "closure of other facts".to_string()
            }
            None => "not recorded".to_string(),
        };
        *answered_by.entry(how).or_default() += 1;
    }
    let faults: Vec<_> = model.check(body).iter().map(|f| f.code).collect();
    json!({
        "faults": faults,
        "incidence_questions": asked.len(),
        "answered_by": answered_by,
        "answered_by_definitions_alone": by_definition,
    })
}

fn blind_check(model: &Model, body: BodyId) -> Value {
    let mut faults = BTreeMap::<&str, usize>::new();
    for fault in model.check(body) {
        *faults.entry(fault.code).or_default() += 1;
    }
    json!({"faults": faults})
}

#[derive(Default)]
struct Tally {
    by_definition: usize,
    proved: usize,
    contradicted: usize,
    not_exact: usize,
    refused: BTreeMap<&'static str, usize>,
}

impl Tally {
    fn outcome<T>(&mut self, outcome: Outcome<T>, holds: impl FnOnce(&T) -> bool) {
        match outcome {
            Outcome::ProvenCoincident(_) => self.proved += 1,
            Outcome::Proven(value) if holds(&value) => self.proved += 1,
            Outcome::Proven(_) => self.contradicted += 1,
            Outcome::NotExact(_) => self.not_exact += 1,
            Outcome::Refused(refusal) => self.refuse(&refusal),
        }
    }
    fn refuse(&mut self, refusal: &Refusal) {
        *self.refused.entry(refusal.code).or_default() += 1;
    }
    fn report(&self) -> Value {
        json!({
            "by_definition": self.by_definition,
            "proved": self.proved,
            "contradicted": self.contradicted,
            "not_exact": self.not_exact,
            "refused": self.refused,
        })
    }
}

/// Property 1: every edge's curve lies on the surface of each face that uses
/// it, and its end points lie on its curve; vertex loops lie on their surface.
fn edges_on_faces(blind: &mut Blind) -> Value {
    let body = blind.body;
    let asked = questions(&blind.model, body);
    let mut curves = Tally::default();
    let mut points = Tally::default();
    let (_, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        for (a, b) in asked {
            let geom = blind.model.geom_mut();
            if geom.facts().holds(FactKind::Incident, a, b) {
                match a {
                    GeomKey::Curve(_) => curves.by_definition += 1,
                    _ => points.by_definition += 1,
                }
                continue;
            }
            match (a, b) {
                (GeomKey::Curve(c), GeomKey::Surface(s)) => {
                    curve_on_surface(geom, c, s, &mut curves, &mut dec)
                }
                (GeomKey::Point(p), GeomKey::Curve(c)) => {
                    let outcome = geom.point_on_curve(c, p, &mut dec);
                    points.outcome(outcome, |_| false);
                }
                (GeomKey::Point(p), GeomKey::Surface(s)) => {
                    let outcome = geom.side_of(s, &geom.position(p), &mut dec);
                    points.outcome(outcome, |_| false);
                }
                _ => unreachable!("the check asks only these three kinds"),
            }
        }
    });
    json!({
        "work": work.used,
        "allowance_exhausted": work.exhausted,
        "curves_on_surfaces": curves.report(),
        "points_on_curves_and_surfaces": points.report(),
    })
}

fn curve_on_surface(
    geom: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    tally: &mut Tally,
    dec: &mut Decider,
) {
    match catallax_intersect::intersect_curve_surface(geom, c, s, dec) {
        Ok(answer) => match answer.value {
            catallax_intersect::CurveSurface::Contained if answer.marks.is_empty() => {
                tally.proved += 1
            }
            catallax_intersect::CurveSurface::Contained => tally.not_exact += 1,
            _ => tally.contradicted += 1,
        },
        Err(refusal) => tally.refuse(&refusal),
    }
}

/// Property 2, its surface stage: which face pairs could meet away from their
/// shared edges, and what proving where their carriers meet costs. Trimming
/// each meeting to both faces' domains is not measured here.
fn faces_meet(blind: &mut Blind) -> Value {
    let body = blind.body;
    let faces = blind.model.faces_of_body(body);
    let mut edges_of = BTreeMap::<FaceId, BTreeSet<_>>::new();
    let mut vertices_of = BTreeMap::<FaceId, BTreeSet<_>>::new();
    for &f in &faces {
        for l in &blind.model.face(f).loops {
            match blind.model.loop_(*l) {
                Loop::Vertex(v) => {
                    vertices_of.entry(f).or_default().insert(*v);
                }
                Loop::Coedges(cs) => {
                    for c in cs {
                        let e = blind.model.coedge(*c).edge;
                        edges_of.entry(f).or_default().insert(e);
                        let edge = blind.model.edge(e);
                        for v in [edge.start, edge.end].into_iter().flatten() {
                            vertices_of.entry(f).or_default().insert(v);
                        }
                    }
                }
            }
        }
    }
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::zero(),
        absolute: Rational::parse("1/1000").unwrap(),
    };
    let mut boxes = BTreeMap::new();
    let mut box_refusals = BTreeMap::<&str, usize>::new();
    let (_, box_work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        for &f in &faces {
            match catallax_measure::face_bounding_box(&blind.model, f, &accuracy, &mut dec) {
                Ok(answer) => {
                    boxes.insert(f, answer.value);
                }
                Err(refusal) => *box_refusals.entry(refusal.code).or_default() += 1,
            }
        }
    });
    let overlap = |a: FaceId, b: FaceId| match (boxes.get(&a), boxes.get(&b)) {
        (Some(x), Some(y)) => (0..3).all(|i| x.min[i] <= y.max[i] && y.min[i] <= x.max[i]),
        _ => true,
    };
    let (mut pairs, mut adjacent, mut touching, mut candidates) = (0, 0, 0, 0);
    let mut same_carrier = 0;
    let mut surface_pairs = BTreeSet::new();
    for (i, &a) in faces.iter().enumerate() {
        for &b in &faces[i + 1..] {
            pairs += 1;
            let shares_edge = shares(&edges_of, a, b);
            let shares_vertex = shares(&vertices_of, a, b);
            adjacent += usize::from(shares_edge);
            touching += usize::from(!shares_edge && shares_vertex);
            if !overlap(a, b) {
                continue;
            }
            candidates += 1;
            let (sa, sb) = (blind.model.face(a).surface, blind.model.face(b).surface);
            if sa == sb {
                same_carrier += 1;
            } else {
                surface_pairs.insert((sa.min(sb), sa.max(sb)));
            }
        }
    }
    let mut outcomes = BTreeMap::<String, usize>::new();
    let mut pair_work = 0;
    let mut largest = 0;
    for (a, b) in &surface_pairs {
        let (outcome, work) = with_operation_budget(PAIR_ALLOWANCE, || {
            let mut dec = Decider::default();
            catallax_intersect::intersect_surfaces(blind.model.geom_mut(), *a, *b, &mut dec)
        });
        pair_work += work.used;
        largest = largest.max(work.used);
        let kind = match outcome {
            _ if work.exhausted => "allowance exhausted".to_string(),
            Ok(answer) => match answer.value {
                catallax_intersect::SurfacePair::Disjoint => "disjoint".into(),
                catallax_intersect::SurfacePair::Coincident { .. } => "coincident".into(),
                catallax_intersect::SurfacePair::Meets(components) => {
                    if answer.marks.is_empty() {
                        format!("meets ({} components)", components.len().min(4))
                    } else {
                        "meets, not exact".into()
                    }
                }
            },
            Err(refusal) => format!("refused: {}", refusal.code),
        };
        *outcomes.entry(kind).or_default() += 1;
    }
    json!({
        "face_pairs": pairs,
        "sharing_an_edge": adjacent,
        "sharing_only_a_vertex": touching,
        "box_overlapping_pairs": candidates,
        "box_overlapping_pairs_on_one_carrier": same_carrier,
        "distinct_carrier_pairs_to_intersect": surface_pairs.len(),
        "face_box_work": box_work.used,
        "face_box_refusals": box_refusals,
        "carrier_pair_work": pair_work,
        "largest_carrier_pair_work": largest,
        "carrier_pair_outcomes": outcomes,
    })
}

/// Property 2, whole: the Boolean's arrangement of the body's own faces
/// against one another, every cell kept. A valid body's faces meet only at
/// their shared boundaries, so no face splits: as many cells as faces.
fn arrangement(blind: &mut Blind) -> Value {
    let faces = blind.model.faces_of_body(blind.body).len();
    let mut cells = 0;
    let (answer, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        catallax_boolean::arrange(
            &mut blind.model,
            blind.body,
            &[],
            |_, _, _, _| {
                cells += 1;
                Ok(Some(catallax_boolean::PatchSelection::Keep(Sense::Same)))
            },
            &mut dec,
        )
    });
    let mut report = json!({
        "faces": faces,
        "cells": cells,
        "work": work.used,
        "allowance_exhausted": work.exhausted,
    });
    match answer {
        Ok(answer) => {
            report["resewn_faces"] = json!(blind.model.faces_of_body(answer.body).len());
            report["work_by_stage"] = answer
                .report
                .stages
                .iter()
                .map(|record| (format!("{:?}", record.stage), json!(record.stats.work)))
                .collect::<serde_json::Map<_, _>>()
                .into();
        }
        Err(refused) => {
            report["refused"] =
                json!({"stage": format!("{:?}", refused.stage), "code": refused.refusal.code});
        }
    }
    report
}

/// Property 2 across shells: faces of different shells share no edge, so
/// they must not meet at all. Each shell becomes a body of its own and the
/// Boolean's intersect stage sections every pair; any section edge or point
/// is a contact. A contact within one shell that splits no face is not seen
/// by this or by the arrangement.
fn across_shells(blind: &mut Blind) -> Value {
    let shells: Vec<_> = blind
        .model
        .body(blind.body)
        .lumps
        .iter()
        .flat_map(|l| blind.model.lump(*l).shells.clone())
        .collect();
    let bodies: Vec<_> = shells
        .iter()
        .map(|&shell| {
            let lump = blind.model.add_lump(Lump {
                shells: vec![shell],
            });
            blind.model.add_body(Body {
                kind: catallax_topology::BodyKind::Solid,
                lumps: vec![lump],
                marks: catallax_numbers::MarkLog::new(),
            })
        })
        .collect();
    let (mut pairs, mut touching) = (0, 0);
    let mut refused = BTreeMap::<&str, usize>::new();
    let (_, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        for (i, &a) in bodies.iter().enumerate() {
            for &b in &bodies[i + 1..] {
                pairs += 1;
                match catallax_boolean::intersect(&mut blind.model, a, b, &mut dec) {
                    Ok(graph) => {
                        let graph = graph.value;
                        touching +=
                            usize::from(!graph.edges().is_empty() || !graph.points().is_empty());
                    }
                    Err(refusal) => *refused.entry(refusal.code).or_default() += 1,
                }
            }
        }
    });
    json!({
        "shells": shells.len(),
        "shell_pairs": pairs,
        "pairs_that_meet": touching,
        "refused": refused,
        "work": work.used,
        "allowance_exhausted": work.exhausted,
    })
}

fn shares<T: Ord>(sets: &BTreeMap<FaceId, BTreeSet<T>>, a: FaceId, b: FaceId) -> bool {
    match (sets.get(&a), sets.get(&b)) {
        (Some(x), Some(y)) => !x.is_disjoint(y),
        _ => false,
    }
}

/// Property 3: for every face, an exact point strictly inside it, moved a
/// little along the face's outward normal, is outside the body, and moved the
/// other way is inside. The step is `STEP` times the carrier's gradient, so
/// it is positive wherever the gradient is; no other face lies that close to
/// these points in the measured parts, which a checker would have to prove.
const STEP: &str = "1/1000000000";

fn inside_outside(blind: &mut Blind) -> Value {
    let body = blind.body;
    let faces = blind.model.faces_of_body(body);
    let mut tally = Tally::default();
    let mut points_refused = BTreeMap::<&str, usize>::new();
    let step = Real::from_rational(Rational::parse(STEP).unwrap());
    let (_, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        for f in faces {
            let point = match catallax_boolean::face_point(&mut blind.model, f, &mut dec) {
                Ok(answer) => answer.value,
                Err(refusal) => {
                    *points_refused.entry(refusal.code).or_default() += 1;
                    continue;
                }
            };
            let face = blind.model.face(f).clone();
            let normal = match blind
                .model
                .geom()
                .surface(face.surface)
                .normal_at(&point, &mut dec)
            {
                Ok(normal) => normal,
                Err(refusal) => {
                    tally.refuse(&refusal);
                    continue;
                }
            };
            let outward = match face.sense {
                Sense::Same => normal.scaled(&step),
                Sense::Reversed => normal.scaled(&Real::from_i64(-1)).scaled(&step),
            };
            let back = outward.scaled(&Real::from_i64(-1));
            let (mut agree, mut undecided) = (true, false);
            for (offset, expected) in [
                (outward, catallax_boolean::Containment::Outside),
                (back, catallax_boolean::Containment::Inside),
            ] {
                let probe = point.plus(&offset);
                match catallax_boolean::classify_point(&mut blind.model, body, &probe, &mut dec) {
                    Outcome::Proven(found) => agree &= found == expected,
                    Outcome::ProvenCoincident(_) => agree = false,
                    other => {
                        tally.outcome(other, |_| false);
                        undecided = true;
                        break;
                    }
                }
            }
            match (undecided, agree) {
                (true, _) => {}
                (false, true) => tally.proved += 1,
                (false, false) => tally.contradicted += 1,
            }
        }
    });
    let mut report = tally.report();
    report["face_points_refused"] = json!(points_refused);
    report["work"] = json!(work.used);
    report["allowance_exhausted"] = json!(work.exhausted);
    report
}

/// Records an `Incident` fact for every question the topology check asks, as
/// the builders' face helpers do for the faces they make.
fn forge(model: &mut Model, body: BodyId) {
    for (a, b) in questions(model, body) {
        model.geom_mut().record_fact(catallax_numbers::Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("r10.planted"),
        });
    }
}

fn rational(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn unit_block(model: &mut Model, at: &str) -> BodyId {
    let pose = catallax_shapes::Pose::from_rational(
        catallax_shapes::Rotation::identity(),
        catallax_shapes::RVector3::new(rational(at), rational(at), rational(at)),
    );
    let size = [rational("2"), rational("2"), rational("2")];
    catallax_topology::block(model, size, &pose).expect("a block")
}

/// Three faults the topology check cannot see, each planted in a 2 mm block,
/// against the facts-blind recomputation of the property it breaks.
fn planted() -> Value {
    let mut model = Model::new();
    let block = unit_block(&mut model, "0");
    let surface = model
        .geom()
        .surface(model.face(model.faces_of_body(block)[0]).surface)
        .clone();
    let mut rows = serde_json::Map::new();
    for (name, plant) in [
        ("edge off its faces", Plant::EdgeOffItsFaces),
        ("face turned inside out", Plant::FaceTurnedInsideOut),
    ] {
        let mut copy = planted_copy(&model, block, surface.clone(), plant);
        forge(&mut copy.model, copy.body);
        let faults: Vec<_> = copy.model.check(copy.body).iter().map(|f| f.code).collect();
        let mut blind = blind_copy(&copy.model, copy.body, surface.clone());
        let found = match plant {
            Plant::EdgeOffItsFaces => edges_on_faces(&mut blind),
            _ => inside_outside(&mut blind),
        };
        rows.insert(
            name.into(),
            json!({"topology_check_faults": faults, "recomputed": found}),
        );
    }
    // Two overlapping blocks as the two lumps of one body: each shell is
    // sound, so the check passes, but their faces cross.
    let other = unit_block(&mut model, "1");
    let lumps = [block, other].map(|b| model.body(b).lumps[0]).to_vec();
    let overlapping = model.add_body(Body {
        kind: catallax_topology::BodyKind::Solid,
        lumps,
        marks: catallax_numbers::MarkLog::new(),
    });
    let faults: Vec<_> = model.check(overlapping).iter().map(|f| f.code).collect();
    let found = arrangement(&mut blind_copy(&model, overlapping, surface));
    rows.insert(
        "faces crossing away from edges".into(),
        json!({"topology_check_faults": faults, "recomputed": found}),
    );
    Value::Object(rows)
}

/// Denominators for the report: recipes attempted, bodies built, and for each
/// property how many bodies completed it and how many it established.
fn summary(cases: &[Value]) -> Value {
    let built: Vec<_> = cases.iter().filter(|c| c.get("body").is_some()).collect();
    let empty = |v: &Value| v.as_object().is_none_or(|m| m.is_empty());
    let names = |keep: &dyn Fn(&Value) -> bool| -> Vec<Value> {
        built
            .iter()
            .filter(|c| keep(c))
            .map(|c| c["case"].clone())
            .collect()
    };
    let edges = |c: &Value| {
        let e = &c["edges_on_faces"];
        ["curves_on_surfaces", "points_on_curves_and_surfaces"]
            .iter()
            .all(|k| e[k]["contradicted"] == 0 && e[k]["not_exact"] == 0 && empty(&e[k]["refused"]))
            && e["allowance_exhausted"] == false
    };
    let arranged = |c: &Value| {
        c["faces_meet_only_at_edges"]["arrangement"]
            .get("refused")
            .is_none()
    };
    let no_split = |c: &Value| {
        let a = &c["faces_meet_only_at_edges"]["arrangement"];
        arranged(c) && a["cells"] == a["faces"]
    };
    let shells_apart = |c: &Value| {
        let s = &c["faces_meet_only_at_edges"]["across_shells"];
        s["pairs_that_meet"] == 0 && empty(&s["refused"]) && s["allowance_exhausted"] == false
    };
    let oriented = |c: &Value| {
        let o = &c["inside_and_outside"];
        o["proved"] == c["body"]["faces"] && o["allowance_exhausted"] == false
    };
    json!({
        "recipes_attempted": cases.len(),
        "bodies_built": built.len(),
        "edges_on_faces_established": names(&edges),
        "arrangement_completed": names(&arranged),
        "arrangement_completed_without_a_split_face": names(&no_split),
        "no_contact_across_shells": names(&shells_apart),
        "inside_and_outside_agree_on_every_face": names(&oriented),
        "all_measured_checks_pass": names(&|c| edges(c) && no_split(c) && shells_apart(c) && oriented(c)),
    })
}
