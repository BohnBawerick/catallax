//! Exact draft geometry, material retention and consumed-name regressions.
use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_field::{Location, Plan, Value, Workspace};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value as Json, json};

fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(|n| Rational::parse(n).unwrap());
    RPoint3::new(x, y, z)
}
fn pose() -> Json {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn block() -> Vec<Json> {
    vec![json!({"id":"base","op":"block","size":["8","8","4"],"pose":pose()})]
}
fn cylinder(radius: &str) -> Vec<Json> {
    vec![json!({"id":"base","op":"cylinder","radius":radius,"height":"4","pose":pose()})]
}
fn operation(segment: &str) -> Json {
    json!({"id":"drafted","op":"draft","body":"base","faces":[{"face":{"step":"base","role":{"side":{"segment":segment}}}}],"pull":["0","0","1"],"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"45"})
}
fn recipe(mut steps: Vec<Json>, draft: Json) -> Recipe {
    steps.push(draft);
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"draft regression","steps":steps,"result":"drafted"}).to_string()).unwrap()
}

#[test]
fn disconnected_coplanar_faces_do_not_share_only_their_role_names() {
    let mut left = block()[0].clone();
    left["id"] = json!("left");
    left["size"] = json!(["2", "2", "8"]);
    let mut right = left.clone();
    right["id"] = json!("right");
    right["pose"]["translation"] = json!(["6", "0", "0"]);
    let mut bridge = block()[0].clone();
    bridge["id"] = json!("bridge");
    bridge["size"] = json!(["8", "2", "1"]);
    let steps = vec![
        left,
        right,
        bridge,
        json!({"id":"arms","op":"boolean","kind":"union","target":"left","tool":"right"}),
        json!({"id":"source","op":"boolean","kind":"union","target":"arms","tool":"bridge"}),
    ];
    let mut draft = operation("floor");
    draft["body"] = json!("source");
    draft["faces"] = json!([{"face":{"step":"bridge","role":"end"}}]);
    draft["pull"] = json!(["1", "0", "0"]);
    draft["neutral"]["normal"] = json!(["1", "0", "0"]);
    draft["angle"] = json!("-45");
    let recipe = recipe(steps, draft);
    let field = Plan::compile(&recipe).unwrap();
    for (creator, own, other) in [
        ("left", ["1", "1", "8"], ["7", "1", "8"]),
        ("right", ["7", "1", "8"], ["1", "1", "8"]),
    ] {
        let origin = serde_json::from_value(json!({"step":creator,"role":"end"})).unwrap();
        let id = field.select(&origin).unwrap();
        assert_eq!(
            field
                .face_query(id, &point(own), &mut Workspace::default())
                .unwrap()
                .state,
            catallax_field::FaceState::OnFace
        );
        assert_eq!(
            field
                .face_query(id, &point(other), &mut Workspace::default())
                .unwrap()
                .state,
            catallax_field::FaceState::NotOnFace
        );
    }
}

fn both(recipe: &Recipe, cases: &[([&str; 3], Location)]) {
    let mut brep = evaluate(recipe, &mut Decider::default()).unwrap();
    assert!(brep.model.check(brep.result).is_empty());
    let field = Plan::compile(recipe).unwrap();
    for (xyz, expected) in cases {
        let p = point(*xyz);
        let query = field.proven(&p, &mut Workspace::default());
        assert_eq!(query.location, *expected, "field at {xyz:?}: {query:?}");
        let answer = classify_point(&mut brep.model, brep.result, &p, &mut Decider::default());
        let agrees = match expected {
            Location::Inside => matches!(answer, Outcome::Proven(Containment::Inside)),
            Location::Outside => matches!(answer, Outcome::Proven(Containment::Outside)),
            Location::Boundary => matches!(answer, Outcome::ProvenCoincident(_)),
            _ => panic!("test needs a decided location"),
        };
        assert!(agrees, "B-rep at {xyz:?}: {answer:?}");
    }
}

fn split_u(keep: Option<Json>) -> Recipe {
    let vertices = [
        (0, 0),
        (8, 0),
        (8, 8),
        (6, 8),
        (6, 1),
        (2, 1),
        (2, 8),
        (0, 8),
    ];
    let edges:Vec<_>=(0..vertices.len()).map(|i| {
        let to=vertices[(i+1)%vertices.len()];
        json!({"kind":"line","id":if i==4 {"floor".to_owned()}else{format!("e{i}")},"to":[to.0.to_string(),to.1.to_string()]})
    }).collect();
    let steps = vec![
        json!({"id":"profile","op":"sketch","plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},"loops":[edges]}),
        json!({"id":"base","op":"extrude","profile":"profile","distance":"2"}),
    ];
    let mut draft = operation("floor");
    draft["pull"] = json!(["1", "0", "0"]);
    draft["neutral"] = json!({"point":["0","0","0"],"normal":["1","0","0"]});
    if let Some(keep) = keep {
        draft["keep_regions"] = keep;
    }
    recipe(steps, draft)
}
#[test]
fn redrafting_a_mirrored_nonconvex_definition_keeps_its_copy_binding() {
    let mut document = serde_json::to_value(split_u(None)).unwrap();
    document["steps"][2]["angle"] = json!("-45");
    let source: catallax_recipe::Origin =
        serde_json::from_value(document["steps"][2]["faces"][0]["face"].clone()).unwrap();
    let copied = catallax_features::copy_origin(
        &catallax_recipe::StepId::new("reflected").unwrap(),
        1,
        &source,
    );
    document["steps"].as_array_mut().unwrap().push(json!({
        "id":"reflected","op":"mirror","body":"drafted","keep_original":false,
        "plane":{"point":["0","0","0"],"normal":["1","0","0"]}
    }));
    document["steps"].as_array_mut().unwrap().push(json!({
        "id":"redrafted","op":"draft","body":"reflected","faces":[{"face":copied}],
        "pull":["-1","0","0"],"neutral":{"point":["0","0","0"],"normal":["1","0","0"]},
        "angle":"0"
    }));
    document["result"] = json!("redrafted");
    let recipe = Recipe::parse(&document.to_string()).unwrap();
    both(
        &recipe,
        &[
            (["-4", "-1", "1/2"], Location::Inside),
            (["-4", "-1", "4"], Location::Outside),
            (["-1", "-1", "4"], Location::Inside),
            (["-7", "-1", "4"], Location::Inside),
        ],
    );
    let field = Plan::compile(&recipe).unwrap();
    let query = field
        .face_query(
            field.select(&copied).unwrap(),
            &point(["-4", "-1", "1"]),
            &mut Workspace::default(),
        )
        .unwrap();
    assert_eq!(query.state, catallax_field::FaceState::OnFace);
}

#[test]
fn coincident_cylindrical_sources_draft_once_and_keep_both_origin_bindings() {
    let mut steps = cylinder("8");
    let mut alias = steps[0].clone();
    alias["id"] = json!("alias");
    steps.push(alias);
    steps.push(json!({"id":"source","op":"boolean","kind":"union","target":"base","tool":"alias"}));
    let mut draft = operation("wall");
    draft["body"] = json!("source");
    let recipe = recipe(steps, draft);
    both(
        &recipe,
        &[
            (["4", "0", "3"], Location::Inside),
            (["5", "0", "3"], Location::Boundary),
            (["6", "0", "3"], Location::Outside),
        ],
    );
    let field = Plan::compile(&recipe).unwrap();
    for creator in ["base", "alias"] {
        let origin = serde_json::from_value(json!({
            "step":creator,"role":{"side":{"segment":"wall"}}
        }))
        .unwrap();
        let query = field
            .face_query(
                field.select(&origin).unwrap(),
                &point(["5", "0", "3"]),
                &mut Workspace::default(),
            )
            .unwrap();
        assert_eq!(query.state, catallax_field::FaceState::OnFace);
    }
}

fn boolean_cut_u(through: bool) -> Recipe {
    let mut base = block()[0].clone();
    base["size"] = json!(["8", "2", "8"]);
    base["pose"]["translation"] = json!(["0", "-2", "0"]);
    let mut steps = vec![base];
    let mut draft = operation("floor");
    if through {
        steps.push(json!({
            "id":"notch-profile","op":"sketch",
            "plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
            "loops":[[
                {"kind":"line","id":"floor","to":["6","1"]},
                {"kind":"line","id":"a","to":["6","9"]},
                {"kind":"line","id":"b","to":["2","9"]},
                {"kind":"line","id":"c","to":["2","1"]}
            ]]
        }));
        steps.push(
            json!({"id":"source","op":"extrude_through","target":"base","profile":"notch-profile"}),
        );
        draft["faces"][0]["face"]["step"] = json!("source");
    } else {
        let mut notch = block()[0].clone();
        notch["id"] = json!("notch");
        notch["size"] = json!(["4", "4", "8"]);
        notch["pose"]["translation"] = json!(["2", "-3", "1"]);
        steps.push(notch);
        steps.push(
            json!({"id":"source","op":"boolean","kind":"subtract","target":"base","tool":"notch"}),
        );
        draft["faces"][0]["face"] = json!({"step":"notch","role":"start"});
    }
    draft["body"] = json!("source");
    draft["pull"] = json!(["1", "0", "0"]);
    draft["neutral"]["normal"] = json!(["1", "0", "0"]);
    recipe(steps, draft)
}

#[test]
fn coactive_origins_do_not_move_a_disconnected_piece_with_only_the_alias() {
    let mut seed = block()[0].clone();
    seed["id"] = json!("seed");
    seed["size"] = json!(["14", "8", "4"]);
    let mut left = block()[0].clone();
    left["id"] = json!("left");
    left["size"] = json!(["6", "8", "4"]);
    let mut tool = block()[0].clone();
    tool["id"] = json!("tool");
    tool["size"] = json!(["4", "10", "6"]);
    tool["pose"]["translation"] = json!(["6", "-1", "-1"]);
    let steps = vec![
        seed,
        left,
        tool,
        json!({"id":"split","op":"boolean","kind":"subtract","target":"seed","tool":"tool"}),
        json!({"id":"source","op":"boolean","kind":"union","target":"left","tool":"split"}),
    ];
    let mut draft = operation("x_max");
    draft["body"] = json!("source");
    draft["faces"] = json!([{"face":{"step":"left","role":"end"}}]);
    draft["pull"] = json!(["1", "0", "0"]);
    draft["neutral"]["normal"] = json!(["1", "0", "0"]);
    draft["keep_regions"] = json!([["1", "4", "1"], ["11", "4", "1"]]);
    let recipe = recipe(steps, draft);
    both(
        &recipe,
        &[
            (["1", "4", "2"], Location::Inside),
            (["1", "4", "7/2"], Location::Outside),
            (["11", "4", "7/2"], Location::Inside),
            (["11", "4", "4"], Location::Boundary),
        ],
    );
    let field = Plan::compile(&recipe).unwrap();
    let query = field.proven(&point(["11", "4", "5"]), &mut Workspace::default());
    assert!(matches!(query.value, Value::TrueDistance(interval) if interval.contains(1.0)));
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let expected: std::collections::BTreeSet<_> = built
        .consumed
        .iter()
        .map(|event| (&event.source, &event.by))
        .collect();
    let actual: std::collections::BTreeSet<_> = field
        .consumed_faces()
        .iter()
        .map(|event| (&event.source, &event.by))
        .collect();
    assert_eq!(
        actual, expected,
        "consumption retains the complete physical source origin set"
    );
}

#[test]
fn boolean_boundary_draft_preserves_named_creators_despite_shared_geometry() {
    for reflected in [false, true] {
        let mut masked = block()[0].clone();
        masked["id"] = json!("masked");
        let mut gate = block()[0].clone();
        gate["id"] = json!("gate");
        gate["size"] = json!(["1", "1", "1"]);
        gate["pose"]["translation"] = json!(["40", "40", "40"]);
        let mut steps = block();
        steps.extend([
            masked,
            gate,
            json!({"id":"masked-result","op":"boolean","kind":"intersect","target":"masked","tool":"gate"}),
            json!({"id":"source","op":"boolean","kind":"union","target":"base","tool":"masked-result"}),
        ]);
        let mut draft = operation("x_max");
        draft["body"] = json!("source");
        let base: catallax_recipe::Origin =
            serde_json::from_value(draft["faces"][0]["face"].clone()).unwrap();
        let mut masked = base.clone();
        masked.step = catallax_recipe::StepId::new("masked").unwrap();
        let (base, masked) = if reflected {
            steps.push(json!({
                "id":"reflected","op":"mirror","body":"source","keep_original":false,
                "plane":{"point":["0","0","0"],"normal":["1","0","0"]}
            }));
            draft["body"] = json!("reflected");
            let id = catallax_recipe::StepId::new("reflected").unwrap();
            (
                catallax_features::copy_origin(&id, 1, &base),
                catallax_features::copy_origin(&id, 1, &masked),
            )
        } else {
            (base, masked)
        };
        draft["faces"][0]["face"] = serde_json::to_value(&base).unwrap();
        let recipe = recipe(steps, draft);
        let (inside, boundary, outside) = if reflected {
            (["-4", "4", "3"], ["-5", "4", "3"], ["-6", "4", "3"])
        } else {
            (["4", "4", "3"], ["5", "4", "3"], ["6", "4", "3"])
        };
        both(
            &recipe,
            &[
                (inside, Location::Inside),
                (boundary, Location::Boundary),
                (outside, Location::Outside),
            ],
        );
        let field = Plan::compile(&recipe).unwrap();
        for (origin, expected) in [
            (&base, catallax_field::FaceState::OnFace),
            (&masked, catallax_field::FaceState::NotOnFace),
        ] {
            let query = field
                .face_query(
                    field.select(origin).unwrap(),
                    &point(boundary),
                    &mut Workspace::default(),
                )
                .unwrap();
            assert_eq!(query.state, expected, "{origin:?}");
        }
        assert!(
            field.consumed_faces().is_empty(),
            "already-absent masked operands must not become draft consumption"
        );
    }
}

#[test]
fn a_boolean_cut_input_rebuilds_and_drafts_its_actual_trimmed_boundary() {
    let recipe = boolean_cut_u(false);
    both(
        &recipe,
        &[
            (["1", "-1", "4"], Location::Inside),
            (["7", "-1", "4"], Location::Outside),
            (["4", "-1", "1/2"], Location::Outside),
        ],
    );
    let field = Plan::compile(&recipe).unwrap();
    let query = field.proven(&point(["3", "-1", "-1"]), &mut Workspace::default());
    assert!(
        matches!(query.value, Value::TrueDistance(interval) if interval.contains(2_f64.sqrt()))
    );
    let notch = catallax_recipe::StepId::new("notch").unwrap();
    assert!(field.consumed_faces().iter().any(|event| {
        event
            .source
            .iter()
            .any(|origin| origin.step == notch && origin.role == Role::Start)
    }));
}

#[test]
fn through_cut_source_clipping_never_publishes_construction_caps() {
    let recipe = boolean_cut_u(true);
    both(
        &recipe,
        &[
            (["1", "-1", "4"], Location::Inside),
            (["7", "-1", "4"], Location::Outside),
            (["1", "-3", "4"], Location::Outside),
        ],
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 6);
    let field = Plan::compile(&recipe).unwrap();
    let query = field.proven(&point(["1", "-3", "4"]), &mut Workspace::default());
    assert!(matches!(query.value, Value::TrueDistance(interval) if interval.contains(1.0)));
    assert!(field.consumed_faces().iter().any(|event| {
        event.source.iter().any(
            |origin| matches!(&origin.role, Role::Side{segment,..} if segment.to_string()=="floor"),
        )
    }));
}

#[test]
fn neutral_retention_keeps_remote_material_connected_through_the_drafted_bridge() {
    let mut document = serde_json::to_value(split_u(None)).unwrap();
    document["steps"][2]["angle"] = json!("-45");
    let recipe = Recipe::parse(&document.to_string()).unwrap();
    both(
        &recipe,
        &[
            (["7", "-1", "4"], Location::Inside),
            (["4", "-1", "4"], Location::Inside),
            (["4", "-1", "6"], Location::Outside),
            (["4", "-1", "5"], Location::Boundary),
        ],
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(built.model.body(built.result).lumps.len(), 1);
    let field = Plan::compile(&recipe).unwrap();
    let q = field.proven(&point(["4", "-1", "6"]), &mut Workspace::default());
    assert!(
        matches!(q.value,Value::TrueDistance(interval) if interval.contains(0.5_f64.sqrt())),
        "{q:?}"
    );
}

#[test]
fn swallowed_floor_separates_components_without_merging_coplanar_names() {
    for (keep, left) in [(None, true), (Some(json!([["7", "-1", "4"]])), false)] {
        let recipe = split_u(keep);
        let (inside, outside, corner) = if left {
            (["1", "-1", "4"], ["7", "-1", "4"], ["3", "-1", "-1"])
        } else {
            (["7", "-1", "4"], ["1", "-1", "4"], ["5", "-1", "-1"])
        };
        both(
            &recipe,
            &[
                (inside, Location::Inside),
                (outside, Location::Outside),
                (["4", "-1", "1/2"], Location::Outside),
            ],
        );
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert_eq!(built.model.body(built.result).lumps.len(), 1);
        assert_eq!(built.model.counts(built.result).faces, 6);
        let names: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .collect();
        let has = |id: &str| {
            names
                .iter()
                .any(|o| matches!(&o.role,Role::Side{segment,..} if segment.to_string()==id))
        };
        assert_eq!(has("e6"), left);
        assert_eq!(has("e2"), !left);
        assert!(!has("floor"));
        let field = Plan::compile(&recipe).unwrap();
        for history in [&built.consumed[..], field.consumed_faces()] {
            assert!(!history.iter().any(|r| {
                r.source
                    .iter()
                    .any(|o| matches!(o.role, Role::Start | Role::End))
            }));
            assert!(history.iter().any(|r| r.source.iter().any(
                |o| matches!(&o.role,Role::Side{segment,..} if segment.to_string()=="floor")
            )));
        }
        let q = field.proven(&point(corner), &mut Workspace::default());
        assert!(
            matches!(q.value,Value::TrueDistance(interval) if interval.contains(2_f64.sqrt())),
            "{q:?}"
        );
    }
}

#[test]
fn retained_picks_accept_internal_subdivision_seams_and_deduplicate_components() {
    let recipe = split_u(Some(json!([
        ["1/2", "-1", "1/2"],
        ["1", "-1", "4"],
        ["7", "-1", "4"]
    ])));
    both(
        &recipe,
        &[
            (["1/2", "-1", "1/2"], Location::Inside),
            (["7", "-1", "4"], Location::Inside),
            (["4", "-1", "1/2"], Location::Outside),
            (["2", "-1", "1"], Location::Boundary),
        ],
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(built.model.body(built.result).lumps.len(), 2);
    assert_eq!(built.model.counts(built.result).faces, 12);
    let field = Plan::compile(&recipe).unwrap();
    let q = field.proven(&point(["1/2", "-1", "1/2"]), &mut Workspace::default());
    assert!(
        matches!(q.value,Value::TrueDistance(interval) if interval.contains(-0.5)),
        "{q:?}"
    );
}

#[test]
fn reflected_retained_components_preserve_trim_normals_and_distance_proofs() {
    let mut document =
        serde_json::to_value(split_u(Some(json!([["1", "-1", "4"], ["7", "-1", "4"]])))).unwrap();
    document["steps"].as_array_mut().unwrap().push(json!({
        "id":"reflected","op":"mirror","body":"drafted","keep_original":false,
        "plane":{"point":["0","0","0"],"normal":["1","0","0"]}
    }));
    document["result"] = json!("reflected");
    let recipe = Recipe::parse(&document.to_string()).unwrap();
    both(
        &recipe,
        &[
            (["-1/2", "-1", "1/2"], Location::Inside),
            (["-7", "-1", "4"], Location::Inside),
            (["-4", "-1", "1/2"], Location::Outside),
            (["-2", "-1", "4"], Location::Boundary),
        ],
    );
    let field = Plan::compile(&recipe).unwrap();
    for (xyz, expected) in [(["-1/2", "-1", "1/2"], -0.5), (["-4", "-1", "1/2"], 2.0)] {
        let q = field.proven(&point(xyz), &mut Workspace::default());
        assert!(
            matches!(q.value,Value::TrueDistance(interval) if interval.contains(expected)),
            "{q:?}"
        );
        assert!(q.distance.millimetres.hi() - q.distance.millimetres.lo() <= 1e-9);
    }
}

#[test]
fn an_oblique_neutral_does_not_reverse_the_material_sense() {
    let mut steps = block();
    steps[0]["pose"]["translation"] = json!(["0", "0", "-4"]);
    let mut draft = operation("x_max");
    draft["neutral"] = json!({"point":["8","0","-2"],"normal":["2","0","1"]});
    let recipe = recipe(steps, draft);
    both(
        &recipe,
        &[
            (["9", "4", "-7/2"], Location::Inside),
            (["8", "4", "-1"], Location::Outside),
            (["7", "4", "-1"], Location::Boundary),
        ],
    );
}
#[test]
fn planar_draft_changes_material_and_measures_true_trimmed_distance() {
    let recipe = recipe(block(), operation("x_max"));
    both(
        &recipe,
        &[
            (["4", "4", "2"], Location::Inside),
            (["7", "4", "3"], Location::Outside),
            (["5", "4", "3"], Location::Boundary),
        ],
    );
    let plan = Plan::compile(&recipe).unwrap();
    let right = plan
        .faces()
        .iter()
        .find(|f| matches!(&f.origin.role,Role::Side{segment,..} if segment.to_string()=="x_max"))
        .unwrap();
    let selected = plan.select(&right.origin).unwrap();
    let face = plan
        .face_query(selected, &point(["5", "4", "3"]), &mut Workspace::default())
        .unwrap();
    assert_eq!(face.state, catallax_field::FaceState::OnFace);
    for (xyz, distance) in [
        (["4", "4", "2"], -2_f64.sqrt()),
        (["7", "4", "3"], 2_f64.sqrt()),
    ] {
        let q = plan.proven(&point(xyz), &mut Workspace::default());
        let Value::TrueDistance(interval) = q.value else {
            panic!("draft cell must not claim a residual is distance")
        };
        assert!(interval.contains(distance), "{q:?} excludes {distance}");
        assert!(q.distance.millimetres.contains(distance.abs()));
        assert!(q.distance.millimetres.hi() - q.distance.millimetres.lo() <= 1e-9);
    }
}
#[test]
fn axial_cylinder_becomes_the_neutral_circle_cone() {
    let recipe = recipe(cylinder("8"), operation("wall"));
    both(
        &recipe,
        &[
            (["1", "0", "3"], Location::Inside),
            (["6", "0", "3"], Location::Outside),
            (["5", "0", "3"], Location::Boundary),
        ],
    );
    let plan = Plan::compile(&recipe).unwrap();
    let q = plan.proven(&point(["6", "0", "3"]), &mut Workspace::default());
    assert!(
        matches!(q.value,Value::TrueDistance(interval) if interval.contains(0.5_f64.sqrt())),
        "{q:?}"
    );
}
#[test]
fn apex_swallows_the_end_without_an_invented_cap() {
    let recipe = recipe(cylinder("2"), operation("wall"));
    both(
        &recipe,
        &[
            (["0", "0", "1"], Location::Inside),
            (["0", "0", "3"], Location::Outside),
            (["0", "0", "2"], Location::Boundary),
        ],
    );
    let brep = evaluate(&recipe, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    assert_eq!(brep.model.counts(brep.result).faces, 2);
    for history in [&brep.consumed[..], plan.consumed_faces()] {
        assert!(
            history
                .iter()
                .any(|record| record.source.iter().any(|o| o.role == Role::End))
        );
        assert!(history.iter().all(|record| {
            record
                .by
                .iter()
                .all(|o| matches!(o.role, Role::Side { .. }))
        }));
    }
    let end = plan
        .faces()
        .iter()
        .find(|f| f.origin.role == Role::End)
        .unwrap();
    let query = plan
        .face_query(
            plan.select(&end.origin).unwrap(),
            &point(["0", "0", "4"]),
            &mut Workspace::default(),
        )
        .unwrap();
    assert_eq!(query.state, catallax_field::FaceState::NotOnFace);
}
#[test]
fn a_fully_swallowed_planar_neighbour_loses_its_selectable_name() {
    let loops = json!([[{"kind":"line","id":"bottom","to":["8","0"]},{"kind":"line","id":"right","to":["8","7"]},{"kind":"line","id":"bevel","to":["7","8"]},{"kind":"line","id":"top","to":["0","8"]},{"kind":"line","id":"left","to":["0","0"]}]]);
    let steps = vec![
        json!({"id":"profile","op":"sketch","plane":pose(),"loops":loops}),
        json!({"id":"base","op":"extrude","profile":"profile","distance":"4"}),
    ];
    let mut draft = operation("right");
    draft["neutral"]["point"] = json!(["0", "0", "-2"]);
    draft["keep_regions"] = json!([["1", "1", "1"]]);
    let recipe = recipe(steps, draft);
    both(
        &recipe,
        &[
            (["1", "1", "3"], Location::Inside),
            (["7", "1", "3"], Location::Outside),
        ],
    );
    let brep = evaluate(&recipe, &mut Decider::default()).unwrap();
    let field = Plan::compile(&recipe).unwrap();
    let bevel = |role: &Role| matches!(role,Role::Side{segment,..} if segment.to_string()=="bevel");
    for history in [&brep.consumed[..], field.consumed_faces()] {
        assert!(
            history
                .iter()
                .any(|record| record.source.iter().any(|o| bevel(&o.role)))
        );
    }
    assert!(
        brep.model.faces_of_body(brep.result).iter().all(|f| brep
            .model
            .face_origins(*f)
            .iter()
            .all(|o| !bevel(&o.role)))
    );
    let bevel = field
        .faces()
        .iter()
        .find(|f| bevel(&f.origin.role))
        .unwrap();
    let query = field
        .face_query(
            field.select(&bevel.origin).unwrap(),
            &point(["15/2", "15/2", "2"]),
            &mut Workspace::default(),
        )
        .unwrap();
    assert_eq!(query.state, catallax_field::FaceState::NotOnFace);
}
#[test]
fn retention_distinguishes_facets_from_edge_and_vertex_contacts() {
    for (normal, neutral, keep) in [
        (["1", "0", "0"], ["8", "0", "0"], true),
        (["1", "1", "0"], ["8", "8", "0"], false),
        (["1", "1", "1"], ["8", "8", "4"], false),
        (["1", "0", "0"], ["-1", "0", "0"], false),
    ] {
        let mut draft = operation("x_max");
        draft["faces"] = json!([]);
        draft["neutral"] = json!({"point":neutral,"normal":normal});
        let recipe = recipe(block(), draft);
        let brep = evaluate(&recipe, &mut Decider::default()).unwrap();
        let field = Plan::compile(&recipe).unwrap();
        assert_eq!(
            brep.model.body(brep.result).kind == catallax_topology::BodyKind::Empty,
            !keep
        );
        assert_eq!(field.is_proven_empty(), !keep);
    }
}
#[test]
fn explicit_region_picks_override_default_and_reject_boundary_or_exterior() {
    for (pick, valid, empty) in [
        (json!([["1", "1", "1"]]), true, false),
        (json!([]), true, true),
        (json!([["0", "1", "1"]]), false, false),
        (json!([["20", "1", "1"]]), false, false),
    ] {
        let mut draft = operation("x_max");
        draft["faces"] = json!([]);
        draft["neutral"]["point"] = json!(["0", "0", "-2"]);
        draft["keep_regions"] = pick;
        let recipe = recipe(block(), draft);
        if valid {
            let brep = evaluate(&recipe, &mut Decider::default()).unwrap();
            assert_eq!(
                brep.model.body(brep.result).kind == catallax_topology::BodyKind::Empty,
                empty
            );
            assert_eq!(Plan::compile(&recipe).unwrap().is_proven_empty(), empty);
        } else {
            assert!(
                matches!(evaluate(&recipe,&mut Decider::default()),Err(BuildError::Geometry{refusal,..}) if refusal.code=="draft.invalid-keep-region")
            );
            assert_eq!(
                Plan::compile(&recipe).unwrap_err().code,
                "draft.invalid-keep-region"
            );
        }
    }
}
#[test]
fn an_oblique_cylindrical_neutral_is_unfinished_not_an_impossible_input() {
    let mut draft = operation("wall");
    draft["neutral"]["normal"] = json!(["-1", "0", "10"]);
    let recipe = recipe(cylinder("8"), draft);
    assert!(
        matches!(evaluate(&recipe,&mut Decider::default()),Err(BuildError::Geometry{refusal,..}) if refusal.code=="draft.unsupported-oblique-neutral")
    );
    assert_eq!(
        Plan::compile(&recipe).unwrap_err().code,
        "draft.unsupported-oblique-neutral"
    );
}

#[test]
fn drafted_flange_unites_with_plate_inside_scoreboard_work_allowance() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P1-drafted-bracket.json")).unwrap();
    let base = recipe
        .steps
        .iter()
        .position(|step| step.id.to_string() == "base")
        .unwrap();
    recipe.steps.truncate(base + 1);
    recipe.result = catallax_recipe::StepId::new("base").unwrap();
    let (_, work) = catallax_numbers::with_operation_budget(1_000_000, || {
        both(
            &recipe,
            &[
                (["50", "2", "10"], Location::Inside),
                (["58", "2", "10"], Location::Outside),
                (["59", "20", "2"], Location::Inside),
            ],
        );
    });
    assert!(!work.exhausted);
}
