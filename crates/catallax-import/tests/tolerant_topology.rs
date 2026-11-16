//! Tolerance storage does not grant solid validity or exact certification.
mod common;

use catallax_import::{
    ImportRefusal, import_step_bodies, read_step, sha256_hex,
    tolerant::{self, Error, ToleranceEntity, ToleranceOptions, TolerantFile},
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use catallax_topology::Model;
use common::{BoxFixture, box_step, boxes_step, r};

fn read(source: &str) -> TolerantFile {
    tolerant::read_step(
        source.as_bytes(),
        &sha256_hex(source.as_bytes()),
        ToleranceOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn shared_entities_keep_all_source_geometry_and_oriented_uses() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let original = read_step(source.as_bytes()).unwrap();
    let graph = read(&source);
    assert_eq!(graph.definition(&mut Decider::default()).unwrap(), original);
    let body = &graph.bodies()[0];
    assert_eq!(body.vertices().len(), 8);
    assert_eq!(body.edges().len(), 12);
    assert_eq!(body.shells()[0].faces.len(), 6);
    assert!(body.declared_solid());
    assert!(!body.has_tolerance_repairs());
    assert_eq!(graph.source_sha256(), sha256_hex(source.as_bytes()));
    let mut uses = 0;
    for (face, old) in body.shells()[0].faces.iter().zip(&original.solids[0].faces) {
        assert_eq!(
            (&face.name, &face.surface, face.same_sense),
            (&old.name, &old.surface, old.same_sense)
        );
        for (boundary, old) in face.loops.iter().zip(&old.loops) {
            assert_eq!(
                (boundary.outer, boundary.orientation),
                (old.outer, old.orientation)
            );
            for (coedge, old) in boundary.edges.iter().zip(&old.edges) {
                uses += 1;
                let edge = &body.edges()[coedge.edge];
                assert_eq!(
                    (&edge.name, &edge.curve, edge.curve_same_sense),
                    (&old.edge, &old.curve, old.curve_same_sense)
                );
                assert_eq!(coedge.forward, old.forward);
                for (index, old) in edge.vertices.into_iter().zip([&old.start, &old.end]) {
                    assert_eq!(&body.vertices()[index].source, old);
                    assert_eq!(body.vertices()[index].tolerance_mm(), &r("1e-7"));
                }
                assert_eq!(edge.tolerance_mm(), &r("1e-7"));
            }
        }
    }
    assert_eq!(uses, 24);
}

#[test]
fn an_open_inconsistent_source_can_be_stored_without_solid_admission() {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    fixture.omit = Some(5);
    let source = box_step(fixture);
    assert!(
        import_step_bodies(
            &mut Model::new(),
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("source").unwrap()
        )
        .is_err()
    );
    let graph = read(&source);
    assert_eq!(graph.bodies()[0].shells()[0].faces.len(), 5);
    assert!(graph.bodies()[0].declared_solid());
    assert!(graph.bodies()[0].tolerance_changes().is_empty());
}

#[test]
fn edge_growth_propagates_to_both_endpoints_and_keeps_other_entities_local() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let mut graph = read(&source);
    let body = graph.bodies()[0].name().to_owned();
    let edge = graph.bodies()[0].edges()[0].clone();
    assert!(
        graph
            .grow_edge_tolerance(
                &body,
                &edge.name,
                r("0.002"),
                "sewing.edge-gap",
                &mut Decider::default()
            )
            .unwrap()
    );
    let result = &graph.bodies()[0];
    assert_eq!(result.edges()[0].tolerance_mm(), &r("0.002"));
    for (index, vertex) in result.vertices().iter().enumerate() {
        assert_eq!(
            vertex.tolerance_mm(),
            &r(if edge.vertices.contains(&index) {
                "0.002"
            } else {
                "1e-7"
            })
        );
    }
    assert!(
        result.edges()[1..]
            .iter()
            .all(|e| e.tolerance_mm() == &r("1e-7"))
    );
    assert_eq!(result.tolerance_changes().len(), 3);
    assert_eq!(
        result.tolerance_changes()[0].entity,
        ToleranceEntity::Edge(edge.name.clone())
    );
    assert!(
        result.tolerance_changes()[1..]
            .iter()
            .all(|c| c.propagated_from_edge.as_deref() == Some(&edge.name))
    );
    assert!(result.has_tolerance_repairs());
    assert_eq!(result.maximum_tolerance_mm(), &r("0.002"));
    let unchanged = graph.clone();
    assert!(
        !graph
            .grow_edge_tolerance(
                &body,
                &edge.name,
                r("0.001"),
                "sewing.edge-gap",
                &mut Decider::default()
            )
            .unwrap()
    );
    assert_eq!(graph, unchanged);
    let vertex = graph.bodies()[0].vertices()[edge.vertices[0]]
        .source
        .name
        .clone();
    graph
        .grow_vertex_tolerance(
            &body,
            &vertex,
            r("0.003"),
            "wire.vertex-gap",
            &mut Decider::default(),
        )
        .unwrap();
    assert_eq!(graph.bodies()[0].maximum_tolerance_mm(), &r("0.003"));
    assert_eq!(graph.bodies()[0].edges()[0].tolerance_mm(), &r("0.002"));
    assert_eq!(
        graph.definition(&mut Decider::default()).unwrap(),
        read_step(source.as_bytes()).unwrap()
    );
}

#[test]
fn cap_and_invalid_requests_leave_the_entire_graph_unchanged() {
    let source = box_step(BoxFixture::new(["1", "1", "1"]));
    let mut graph = read(&source);
    let body = graph.bodies()[0].name().to_owned();
    let edge = graph.bodies()[0].edges()[0].name.clone();
    let vertex = graph.bodies()[0].vertices()[0].source.name.clone();
    for (entity, value, cause) in [
        (edge.as_str(), "0.01000000001", "sewing"),
        (edge.as_str(), "-0.001", "sewing"),
        (edge.as_str(), "0.001", " "),
        ("#missing", "0.001", "sewing"),
    ] {
        let before = graph.clone();
        assert!(
            graph
                .grow_edge_tolerance(&body, entity, r(value), cause, &mut Decider::default())
                .is_err()
        );
        assert_eq!(graph, before);
    }
    let before = graph.clone();
    assert!(
        graph
            .grow_vertex_tolerance(
                &body,
                &vertex,
                r("0.01000000001"),
                "gap",
                &mut Decider::default()
            )
            .is_err()
    );
    assert_eq!(graph, before);
    assert!(
        graph
            .grow_edge_tolerance(&body, &edge, r("0.01"), "at-cap", &mut Decider::default())
            .unwrap()
    );
    assert_eq!(graph.bodies()[0].maximum_tolerance_mm(), &r("0.01"));
}

#[test]
fn exhaustion_at_every_edge_update_boundary_is_atomic() {
    let source = box_step(BoxFixture::new(["1", "1", "1"]));
    let original = read(&source);
    let body = original.bodies()[0].name().to_owned();
    let edge = original.bodies()[0].edges()[0].name.clone();
    let mut successes = 0;
    let mut failures = 0;
    for allowance in 0..100 {
        let mut graph = original.clone();
        let (result, work) = with_operation_budget(allowance, || {
            graph.grow_edge_tolerance(&body, &edge, r("0.002"), "sewing", &mut Decider::default())
        });
        if result.is_err() {
            assert!(work.exhausted);
            assert_eq!(graph, original);
            failures += 1;
        } else {
            assert!(!work.exhausted);
            assert_eq!(graph.bodies()[0].tolerance_changes().len(), 3);
            successes += 1;
        }
    }
    assert!(successes > 0 && failures > 0);
}

#[test]
fn reopened_history_reproduces_values_and_rejects_tampering_and_changed_sources() {
    let source = box_step(BoxFixture::new(["1", "2", "3"]));
    let mut graph = read(&source);
    let body = graph.bodies()[0].name().to_owned();
    let edge = graph.bodies()[0].edges()[0].name.clone();
    graph
        .grow_edge_tolerance(&body, &edge, r("0.002"), "sewing", &mut Decider::default())
        .unwrap();
    let encoded = serde_json::to_vec(&graph.record()).unwrap();
    let record = serde_json::from_slice(&encoded).unwrap();
    assert_eq!(
        tolerant::reopen(source.as_bytes(), &record, &mut Decider::default()).unwrap(),
        graph
    );
    for tamper in 0..5 {
        let mut bad = graph.record();
        match tamper {
            0 => bad.bodies[0].changes[0].before_mm = "1/2".into(),
            1 => bad.bodies[0].changes[0].after_mm = "1/2".into(),
            2 => {
                bad.bodies[0].changes.pop();
            }
            3 => bad.bodies[0].vertices[0].1 = "1/2".into(),
            _ => bad.bodies[0].changes[0].entity = ToleranceEntity::Edge("#renumbered".into()),
        }
        assert!(tolerant::reopen(source.as_bytes(), &bad, &mut Decider::default()).is_err());
    }
    // Entity renumbering can preserve geometry but must not reuse source bindings.
    let changed = source
        .replace("#1=", "#99999=")
        .replace("#1,", "#99999,")
        .replace("#1)", "#99999)");
    assert_ne!(changed, source);
    assert!(matches!(
        tolerant::reopen(changed.as_bytes(), &record, &mut Decider::default()),
        Err(Error::Source(ImportRefusal::HashMismatch { .. }))
    ));
}

#[test]
fn independent_solids_and_surface_declarations_remain_distinct() {
    let source = boxes_step(&[
        BoxFixture::new(["1", "2", "3"]),
        BoxFixture::new(["4", "5", "6"]),
    ]);
    let mut graph = read(&source);
    let second = graph.bodies()[1].clone();
    let first = graph.bodies()[0].name().to_owned();
    let edge = graph.bodies()[0].edges()[0].name.clone();
    graph
        .grow_edge_tolerance(&first, &edge, r("0.002"), "gap", &mut Decider::default())
        .unwrap();
    assert_eq!(graph.bodies()[1], second);
    let source = box_step(BoxFixture::new(["1", "2", "3"]));
    let surface = source
        .lines()
        .map(|line| {
            if line.contains("MANIFOLD_SOLID_BREP") {
                let (id, _) = line.split_once('=').unwrap();
                let shell = line.split(",#").nth(1).unwrap().trim_end_matches(");");
                format!("{id}=SHELL_BASED_SURFACE_MODEL('',(#{shell}));")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let graph = read(&surface);
    assert!(!graph.bodies()[0].declared_solid());
    assert!(graph.bodies()[0].shells()[0].declared_closed);
    assert_eq!(graph.bodies()[0].edges().len(), 12);
}

#[test]
fn options_hash_and_read_work_are_checked_before_a_graph_is_returned() {
    let source = box_step(BoxFixture::new(["1", "1", "1"]));
    for (initial, maximum) in [("-1", "0.01"), ("0.02", "0.01"), ("0", "0")] {
        assert!(
            tolerant::read_step(
                source.as_bytes(),
                &sha256_hex(source.as_bytes()),
                ToleranceOptions {
                    initial_mm: r(initial),
                    maximum_mm: r(maximum)
                },
                &mut Decider::default()
            )
            .is_err()
        );
    }
    assert!(
        tolerant::read_step(
            source.as_bytes(),
            &"0".repeat(64),
            ToleranceOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
    let (result, work) = with_operation_budget(0, || {
        tolerant::read_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            ToleranceOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(result.is_err() && work.exhausted);
}

#[test]
fn vertex_loops_and_closed_edge_endpoints_keep_their_source_identity() {
    let mut soup = common::Soup::new();
    soup.vertex_loop_sphere(["1", "2", "3"], "2");
    let source = soup.finish();
    let graph = read(&source);
    assert_eq!(
        graph.definition(&mut Decider::default()).unwrap(),
        read_step(source.as_bytes()).unwrap()
    );
    assert_eq!(graph.bodies()[0].vertices().len(), 1);
    assert!(graph.bodies()[0].edges().is_empty());
    assert_eq!(
        graph.bodies()[0].shells()[0].faces[0].loops[0].vertex,
        Some(0)
    );
    let mut soup = common::Soup::new();
    soup.cylinder(["0", "0", "0"], 2, "2", "1", common::Caps::Both);
    let source = soup.finish();
    let mut graph = read(&source);
    let body = graph.bodies()[0].name().to_owned();
    let edge = graph.bodies()[0]
        .edges()
        .iter()
        .find(|e| e.vertices[0] == e.vertices[1])
        .unwrap()
        .name
        .clone();
    graph
        .grow_edge_tolerance(
            &body,
            &edge,
            r("0.001"),
            "closed-edge",
            &mut Decider::default(),
        )
        .unwrap();
    // A ring's single endpoint is propagated once, even though used twice.
    assert_eq!(graph.bodies()[0].tolerance_changes().len(), 2);
    assert_eq!(
        graph.definition(&mut Decider::default()).unwrap(),
        read_step(source.as_bytes()).unwrap()
    );
    assert_eq!(
        tolerant::reopen(source.as_bytes(), &graph.record(), &mut Decider::default()).unwrap(),
        graph
    );
}

#[test]
fn normalized_source_units_and_empty_entity_maxima_are_retained() {
    let mut fixture = BoxFixture::new(["1", "2", "3"]);
    fixture.unit = common::Unit::Metre;
    let source = box_step(fixture);
    let graph = read(&source);
    assert_eq!(
        graph.definition(&mut Decider::default()).unwrap(),
        read_step(source.as_bytes()).unwrap()
    );
    assert_eq!(graph.bodies()[0].maximum_tolerance_mm(), &r("1e-7"));
    let mut soup = common::Soup::new();
    soup.unbounded_plane(["0", "0", "0"], 2);
    let graph = read(&soup.finish());
    assert!(graph.bodies()[0].vertices().is_empty() && graph.bodies()[0].edges().is_empty());
    assert_eq!(graph.bodies()[0].maximum_tolerance_mm(), &r("0"));
}
