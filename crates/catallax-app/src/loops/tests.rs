use super::*;
use crate::curves::{CurveSpan, CurveTrim};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../tests/fixtures/trimmed-curves.json");
const POSE: &str = r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#;

fn parse(steps: &str, result: &str) -> Recipe {
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"loops","steps":[{steps}],"result":"{result}"}}"#
    ))
    .unwrap()
}
fn build(recipe: &Recipe) -> ApplicationBody {
    crate::build(recipe, &Default::default(), &mut Decider::default()).unwrap()
}
fn extrusion(scale: i64) -> Recipe {
    let mut value: Value = serde_json::from_str(FIXTURE).unwrap();
    value["steps"][0]["loops"][0][0]["to"] =
        json!([(3 * scale).to_string(), (4 * scale).to_string()]);
    value["steps"][0]["loops"][0][1]["to"] =
        json!([(3 * scale).to_string(), (-4 * scale).to_string()]);
    Recipe::parse(&value.to_string()).unwrap()
}
fn plate_with_hole() -> Recipe {
    parse(
        &format!(
            r#"{{"id":"base","op":"block","size":["30","8","4"],"pose":{POSE}}},
            {{"id":"drill","op":"hole","target":"base","at":[{{"id":"c","point":["3","4","4"],"direction":["0","0","-1"]}}],"diameter":"2","depth":"through","tip":"flat"}}"#
        ),
        "drill",
    )
}

/// Every coedge cycle closes under the reported travel senses, and the
/// per-edge data agrees with the stored topology and the edge query.
fn check(
    model: &Model,
    body: BodyId,
    queried: &FaceLoops,
    edge_query: &dyn Fn(EdgeId) -> EdgeCurve,
) {
    let stored = model.face(queried.face);
    assert_eq!(queried.surface, stored.surface);
    assert_eq!(queried.sense, stored.sense);
    assert_eq!(queried.origins, model.face_origins(queried.face));
    assert_eq!(queried.loops.len(), stored.loops.len());
    for (loop_, &id) in queried.loops.iter().zip(&stored.loops) {
        match (loop_, model.loop_(id)) {
            (
                FaceLoop::Vertex {
                    id: got,
                    endpoint,
                    name,
                },
                Loop::Vertex(vertex),
            ) => {
                assert_eq!(*got, id);
                assert_eq!(endpoint.vertex, *vertex);
                assert_eq!(*name, model.vertex_name(body, *vertex));
            }
            (FaceLoop::Coedges { id: got, coedges }, Loop::Coedges(cycle)) => {
                assert_eq!(*got, id);
                assert_eq!(
                    coedges.iter().map(|c| c.coedge).collect::<Vec<_>>(),
                    *cycle,
                    "stored order is kept"
                );
                for (index, c) in coedges.iter().enumerate() {
                    let next = &coedges[(index + 1) % coedges.len()];
                    assert_eq!(c.end, next.start, "cycle closes under travel senses");
                    let edge = model.edge(c.edge.edge);
                    assert_eq!(c.sense, model.coedge(c.coedge).sense);
                    match (&c.edge.trim, c.sense) {
                        (CurveTrim::Closed, _) => assert_eq!((c.start, c.end), (None, None)),
                        (CurveTrim::Finite { start, end }, Sense::Same) => {
                            assert_eq!((c.start, c.end), (Some(start.vertex), Some(end.vertex)));
                        }
                        (CurveTrim::Finite { start, end }, Sense::Reversed) => {
                            assert_eq!((c.start, c.end), (Some(end.vertex), Some(start.vertex)));
                        }
                    }
                    assert_eq!((edge.start, edge.end) == (None, None), c.start.is_none());
                    let direct = edge_query(c.edge.edge);
                    assert_eq!(direct.sense, c.edge.sense);
                    assert_eq!(direct.name, c.edge.name);
                    assert_eq!(direct.origins, c.edge.origins);
                    assert_eq!(direct.constructions, c.edge.constructions);
                    assert!(direct.definition.same_definition(&c.edge.definition));
                }
            }
            _ => panic!("loop kind differs from topology"),
        }
    }
}
fn check_all(body: &ApplicationBody) -> Vec<FaceLoops> {
    body.names
        .faces
        .keys()
        .map(|&face| {
            let queried = body.face_loops(face, &mut Decider::default()).unwrap();
            check(&body.built.model, body.built.result, &queried, &|edge| {
                body.edge_curve(edge, &mut Decider::default()).unwrap()
            });
            queried
        })
        .collect()
}

#[test]
fn every_extruded_face_returns_closed_stored_cycles_without_display() {
    let mut body = build(&extrusion(1));
    body.display.value.positions.clear();
    body.display.value.triangles.clear();
    body.display.value.edges.clear();
    let faces = check_all(&body);
    assert_eq!(faces.len(), 4);
    assert!(
        faces
            .iter()
            .all(|f| f.loops.len() == 1 && !f.origins.is_empty())
    );
    let curved = faces
        .iter()
        .flat_map(|f| &f.loops)
        .filter_map(|l| match l {
            FaceLoop::Coedges { coedges, .. } => Some(coedges),
            FaceLoop::Vertex { .. } => None,
        })
        .flatten()
        .filter(|c| matches!(c.edge.span, CurveSpan::CircularArc { .. }))
        .count();
    assert_eq!(
        curved, 4,
        "each of the two caps and the curved wall holds one arc use per end"
    );
}

#[test]
fn a_through_hole_keeps_both_loops_in_stored_order() {
    let body = build(&plate_with_hole());
    let faces = check_all(&body);
    assert!(
        faces.iter().any(|f| f.loops.len() == 2),
        "plate faces with a hole"
    );
    for face in faces.iter().filter(|f| f.loops.len() == 2) {
        let lengths: Vec<_> = face
            .loops
            .iter()
            .map(|l| match l {
                FaceLoop::Coedges { coedges, .. } => coedges.len(),
                FaceLoop::Vertex { .. } => 0,
            })
            .collect();
        assert!(lengths.iter().all(|&n| n >= 1));
    }
    assert!(
        faces.iter().flat_map(|f| &f.loops).any(|l| matches!(
            l,
            FaceLoop::Coedges { coedges, .. } if coedges.iter().any(|c| c.sense == Sense::Reversed)
        )),
        "a reversed coedge is exercised"
    );
}

#[test]
fn ring_edges_whole_faces_and_vertex_loops_are_preserved() {
    let cylinder = build(&parse(
        &format!(r#"{{"id":"s","op":"cylinder","radius":"5","height":"2","pose":{POSE}}}"#),
        "s",
    ));
    let faces = check_all(&cylinder);
    assert!(faces.iter().flat_map(|f| &f.loops).any(|l| matches!(
        l,
        FaceLoop::Coedges { coedges, .. } if coedges.iter().any(|c| c.start.is_none())
    )));
    let sphere = build(&parse(
        &format!(r#"{{"id":"s","op":"sphere","radius":"3","pose":{POSE}}}"#),
        "s",
    ));
    let faces = check_all(&sphere);
    assert_eq!(faces.len(), 1);
    assert!(
        faces[0].loops.is_empty(),
        "a whole closed face has zero loops"
    );
    let cone = build(&parse(
        &format!(
            r#"{{"id":"c","op":"cone","base_radius":"3","top_radius":"0","height":"2","pose":{POSE}}}"#
        ),
        "c",
    ));
    let faces = check_all(&cone);
    let apex = faces
        .iter()
        .flat_map(|f| &f.loops)
        .filter(|l| matches!(l, FaceLoop::Vertex { .. }))
        .count();
    assert_eq!(apex, 1, "the apex is a vertex-only boundary");
}

#[test]
fn a_face_outside_the_result_is_a_stable_typed_refusal() {
    let body = build(&plate_with_hole());
    let base = body.built.bodies[&catallax_recipe::StepId::new("base").unwrap()];
    let outside = body
        .built
        .model
        .faces_of_body(base)
        .into_iter()
        .find(|f| !body.names.faces.contains_key(f))
        .expect("the target block keeps faces the result does not");
    for _ in 0..2 {
        let error = body
            .face_loops(outside, &mut Decider::default())
            .unwrap_err();
        assert_eq!(error.code, "app.loop-face");
    }
    // The valid body is still queryable afterwards.
    check_all(&body);
}

#[test]
fn an_exhausted_allowance_refuses_without_discarding_the_body() {
    let body = build(&extrusion(1));
    let face = *body.names.faces.keys().next().unwrap();
    let (result, work) =
        with_operation_budget(1, || body.face_loops(face, &mut Decider::default()));
    assert!(result.is_err() || !work.exhausted);
    check_all(&body);
}

#[test]
fn edits_keep_cycles_closed_and_scale_the_exact_boundary() {
    let first = crate::build_incremental(
        &extrusion(1),
        None,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let edited = crate::build_incremental(
        &extrusion(2),
        Some(first.checkpoint.clone()),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    for body in [&first, &edited] {
        let built = body.checkpoint.built();
        for &face in body.names.faces.keys() {
            let queried = body.face_loops(face, &mut Decider::default()).unwrap();
            check(&built.model, built.result, &queried, &|edge| {
                body.edge_curve(edge, &mut Decider::default()).unwrap()
            });
        }
    }
    let direct = build(&extrusion(2));
    let (a, b) = (check_all(&direct), edited.names.faces.len());
    assert_eq!(a.len(), b);
    let origins = |faces: Vec<Vec<Origin>>| {
        let mut faces = faces;
        faces.sort();
        faces
    };
    let from_edit = origins(
        edited
            .names
            .faces
            .keys()
            .map(|&f| {
                edited
                    .face_loops(f, &mut Decider::default())
                    .unwrap()
                    .origins
            })
            .collect(),
    );
    assert_eq!(
        from_edit,
        origins(a.into_iter().map(|f| f.origins).collect())
    );
}
