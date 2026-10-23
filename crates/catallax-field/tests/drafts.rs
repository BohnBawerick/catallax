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
