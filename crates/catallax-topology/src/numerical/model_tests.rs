//! A numerical body's model is a body of the one model (section 4.5 of the
//! Checkpoint 2 decisions): its vertices and edges store the tolerances
//! admission measured, an edge between two exact surfaces is re-tightened to
//! exact, and `Model::check` is its checker.
use super::test_fixture::{q, quarter};
use super::*;

fn admitted() -> NumericalBody {
    build(
        quarter(),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

/// The tolerance on each vertex and each edge of the body's model, in
/// definition order, and how many of its faces store one.
fn stored(body: &NumericalBody) -> (Vec<Option<Rational>>, Vec<Option<Rational>>, usize) {
    let mm = |tolerance: &Option<Tolerance>| tolerance.as_ref().map(|t| t.mm().clone());
    (
        body.model
            .vertex_handles()
            .into_iter()
            .map(|v| mm(&body.model.vertex(v).tolerance))
            .collect(),
        body.model
            .edge_handles()
            .into_iter()
            .map(|e| mm(&body.model.edge(e).tolerance))
            .collect(),
        body.model
            .faces_of_body(body.body)
            .into_iter()
            .filter(|f| body.model.face(*f).tolerance.is_some())
            .count(),
    )
}

/// What the report states for each vertex and each edge, in definition order.
fn reported(body: &NumericalBody) -> (Vec<Option<Rational>>, Vec<Option<Rational>>) {
    let of = |kind| {
        body.tolerances
            .entities
            .iter()
            .filter(|entity| entity.kind == kind)
            .map(|entity| Some(entity.tolerance_mm.clone()))
            .collect::<Vec<_>>()
    };
    (
        of(tolerance::EntityKind::Vertex),
        of(tolerance::EntityKind::Edge),
    )
}

/// What the model is to store: the report's value, or nothing for an edge
/// that is re-tightened.
fn expected(body: &NumericalBody) -> (Vec<Option<Rational>>, Vec<Option<Rational>>) {
    let (vertices, edges) = reported(body);
    let exact = retighten::exact_edges(&body.definition);
    (
        vertices,
        edges
            .into_iter()
            .zip(exact)
            .map(|(measured, exact)| measured.filter(|_| !exact))
            .collect(),
    )
}

#[test]
fn the_model_stores_what_admission_measured_and_the_one_checker_passes_it() {
    // The wall is 4 thousandths off its edges, so the edges on it keep what
    // admission measured, and the edges between two planes are exact.
    let body = build(
        gapped("1.004"),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let (vertices, edges, faces) = stored(&body);
    assert_eq!((vertices.clone(), edges.clone()), expected(&body));
    assert_eq!(vertices.len(), body.definition.vertices.len());
    assert_eq!(edges.len(), body.definition.edges.len());
    // A face stores none: its surface is the authoritative copy.
    assert_eq!(faces, 0);
    // A stored value is never under rule T1's floor, also where the measured
    // distance is zero. A vertex is a rounded point and always stores one.
    let floor = tolerance::floor_minimum_mm();
    assert!(
        vertices
            .iter()
            .all(|value| value.as_ref().is_some_and(|value| *value >= floor))
    );
    assert!(edges.iter().flatten().all(|value| *value >= floor));
    let kept = edges.iter().flatten().count();
    assert!(0 < kept && kept < edges.len(), "{edges:?}");
    assert_eq!(body.model.check(body.body), Vec::new());
    assert_eq!(body.topology_faults(), Vec::new());
    // By key, as a caller reads it.
    let exact = retighten::exact_edges(&body.definition);
    for entity in &body.tolerances.entities {
        let retightened = entity.kind == tolerance::EntityKind::Edge
            && body
                .definition
                .edges
                .iter()
                .zip(&exact)
                .any(|(edge, exact)| edge.key == entity.key && *exact);
        assert_eq!(
            body.stored_tolerance(entity.kind, &entity.key)
                .map(Tolerance::mm),
            (!retightened).then_some(&entity.tolerance_mm)
        );
    }
    assert_eq!(
        body.stored_tolerance(tolerance::EntityKind::Face, &body.definition.faces[0].key),
        None
    );
    assert_eq!(
        body.stored_tolerance(tolerance::EntityKind::Edge, "no such edge"),
        None
    );
}

#[test]
fn an_edge_whose_curve_lies_on_both_its_surfaces_is_exact_with_its_facts() {
    // Every edge of the quarter cylinder is a line or an arc of its two
    // surfaces, in rationals: all of them are re-tightened.
    let body = admitted();
    let (vertices, edges, _) = stored(&body);
    assert!(edges.iter().all(Option::is_none), "{edges:?}");
    assert!(vertices.iter().all(Option::is_some));
    // The one checker asks an exact edge on an exact face for its fact, and
    // every use has one. The same edges with no fact are not sound, which is
    // the test below.
    assert_eq!(body.model.check(body.body), Vec::new());
    let model = &body.model;
    let mut uses = 0;
    for face in model.faces_of_body(body.body) {
        for boundary in &model.face(face).loops {
            let Loop::Coedges(coedges) = model.loop_(*boundary) else {
                continue;
            };
            for coedge in coedges {
                let curve = model.edge(model.coedge(*coedge).edge).curve;
                assert!(model.geom().facts().holds(
                    catallax_numbers::FactKind::Incident,
                    catallax_shapes::GeomKey::Curve(curve),
                    catallax_shapes::GeomKey::Surface(model.face(face).surface),
                ));
                uses += 1;
            }
        }
    }
    assert_eq!(uses, 2 * body.definition.edges.len());
}

#[test]
fn a_body_whose_tolerances_are_missing_is_not_sound() {
    // The same topology with no tolerance stored claims to be exact, and no
    // fact was ever recorded for it: the checker says so for every vertex on
    // an edge and every edge on a face.
    let body = admitted();
    let (unmeasured, id) =
        ModelGeometry::of(&body.definition).topology(&body.definition, body.kind(), None);
    let faults = unmeasured.check(id);
    let uses: usize = body
        .definition
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .map(Vec::len)
        .sum();
    assert_eq!(faults.len(), 2 * body.definition.edges.len() + uses);
    assert!(
        faults
            .iter()
            .all(|fault| fault.code == "topology.missing-incidence")
    );
}

#[test]
fn a_record_without_tolerances_gets_them_back_from_its_report() {
    // A handover written before format version 3 stored the model with no
    // tolerance on any entity.
    let source = admitted();
    let mut read = source.clone();
    let (model, id) =
        ModelGeometry::of(&read.definition).topology(&read.definition, read.kind(), None);
    (read.model, read.body) = (model, id);
    assert!(!read.topology_faults().is_empty());
    read.restore_model_tolerances().unwrap();
    assert_eq!(stored(&read), stored(&source));
    assert_eq!(read.topology_faults(), Vec::new());
    assert_eq!(read.counts(), source.counts());
    // Entity for entity, it is the model admission built.
    for vertex in source.model.vertex_handles() {
        assert_eq!(read.model.vertex(vertex), source.model.vertex(vertex));
    }
    for edge in source.model.edge_handles() {
        assert_eq!(read.model.edge(edge), source.model.edge(edge));
    }
    let faces = source.model.faces_of_body(source.body);
    assert_eq!(read.model.faces_of_body(read.body), faces);
    for face in faces {
        assert_eq!(read.model.face(face), source.model.face(face));
    }
}

#[test]
fn a_record_older_than_the_report_takes_its_largest_recorded_tolerance() {
    let source = build(
        gapped("1.004"),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let mut read = source.clone();
    read.tolerances = ToleranceReport::default();
    read.restore_model_tolerances().unwrap();
    let largest = tolerance::largest_recorded(&source.evidence);
    assert!(largest > q("0.004"));
    assert_eq!(largest, source.tolerances.largest_mm());
    let (vertices, edges, faces) = stored(&read);
    assert!(
        vertices
            .iter()
            .all(|value| value.as_ref() == Some(&largest))
    );
    // An edge between two planes is exact again. The others take the value.
    assert!(edges.iter().flatten().all(|value| *value == largest));
    let kept = edges.iter().flatten().count();
    assert!(0 < kept && kept < edges.len(), "{edges:?}");
    assert_eq!(faces, 0);
    assert_eq!(read.topology_faults(), Vec::new());

    // Where every recorded distance is zero, rule T1's fixed 1e-7 mm.
    let mut clean = admitted();
    assert!(tolerance::largest_recorded(&clean.evidence) < tolerance::floor_minimum_mm());
    clean.tolerances = ToleranceReport::default();
    clean.restore_model_tolerances().unwrap();
    let (vertices, edges, _) = stored(&clean);
    assert!(
        vertices
            .iter()
            .all(|value| value.as_ref() == Some(&tolerance::floor_minimum_mm()))
    );
    assert!(edges.iter().all(Option::is_none));
}

#[test]
fn a_definition_with_a_missing_edge_is_refused_not_indexed() {
    let mut read = admitted();
    read.definition.faces[0].loops[0][0].edge = read.definition.edges.len();
    assert_eq!(
        read.restore_model_tolerances().unwrap_err().code,
        "handover.invalid-data"
    );
}

/// The quarter cylinder with its wall this far off its edges.
fn gapped(radius: &str) -> NumericalDefinition {
    let mut definition = quarter();
    let Surface::Cylinder(wall) = &mut definition.faces[4].surface else {
        unreachable!()
    };
    wall.radius_sq = (&q(radius) * &q(radius)).into();
    definition
}
