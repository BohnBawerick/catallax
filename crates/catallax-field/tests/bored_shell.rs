//! Independent material, finite-face distance and source-name proofs for bored shells.

use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role};
use catallax_shapes::{RPoint3, Surface, distance::Options};
use serde_json::{Value, json};

const WORK: u64 = 20_000_000;

fn document() -> Value {
    serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/bored-shell.json"
    ))
    .unwrap()
}

fn recipe(source: &Value) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}

fn compile(source: &Value) -> Plan {
    let (plan, work) = with_operation_budget(WORK, || Plan::compile(&recipe(source)).unwrap());
    assert!(!work.exhausted, "{work:?}");
    plan
}

fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(|x| Rational::parse(x).unwrap());
    RPoint3::new(x, y, z)
}

fn material(plan: &Plan, cases: &[([&str; 3], Location)]) {
    for (xyz, want) in cases {
        let query = plan.proven(&point(*xyz), &mut Workspace::default());
        assert_eq!(query.location, *want, "{xyz:?}: {query:?}");
    }
}

fn origin(value: Value) -> Origin {
    serde_json::from_value(value).unwrap()
}

fn inner_bore(plan: &Plan) -> &catallax_field::Face {
    plan.faces()
        .iter()
        .find(|face| {
            matches!(&face.origin.role,
        Role::Inner { of } if of.step.as_str() == "bored")
        })
        .unwrap()
}

#[test]
fn shell_compiles_independently_and_keeps_inner_bore_names_after_thickness_edits() {
    let mut previous = None;
    for thickness in ["2", "3/2"] {
        let mut source = document();
        source["steps"][2]["thickness"] = json!(thickness);
        let plan = compile(&source);
        material(
            &plan,
            &[
                (["20", "15", "1"], Location::Outside),
                (["24", "15", "6"], Location::Inside),
                (["26", "15", "6"], Location::Outside),
                (["10", "10", "1"], Location::Inside),
                (["10", "10", "6"], Location::Outside),
                (["1", "10", "6"], Location::Inside),
                (["10", "10", "13"], Location::Outside),
            ],
        );
        let mut names: Vec<_> = plan
            .faces()
            .iter()
            .filter_map(|face| {
                matches!(face.origin.role, Role::Inner { .. }).then_some(face.origin.clone())
            })
            .collect();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), 6);
        if let Some(previous) = previous {
            assert_eq!(names, previous);
        }
        previous = Some(names);
        let inner = inner_bore(&plan);
        let Surface::Cylinder(cylinder) = &inner.surface else {
            panic!("the bore's inner wall must stay cylindrical");
        };
        let radius = Rational::from_i64(3) + Rational::parse(thickness).unwrap();
        assert_eq!(
            Decider::default()
                .rational_value(&cylinder.radius_sq)
                .unwrap(),
            Some(&radius * &radius)
        );
        let xyz = if thickness == "2" {
            ["25", "15", "6"]
        } else {
            ["49/2", "15", "6"]
        };
        let query = plan
            .face_query(
                plan.select(&inner.origin).unwrap(),
                &point(xyz),
                &mut Workspace::default(),
            )
            .unwrap();
        assert_eq!(query.state, FaceState::OnFace);
        let Direction::Unique(normal) = query.direction else {
            panic!("{query:?}")
        };
        assert!(normal[0].contains(1.0));
        assert!(normal[0].lo() > 0.99);
        assert!(normal[1].contains(0.0) && normal[2].contains(0.0));
        let removed = origin(json!({"step":"stock","role":"end"}));
        assert_eq!(
            plan.face_query(
                plan.select(&removed).unwrap(),
                &point(["10", "10", "12"]),
                &mut Workspace::default()
            )
            .unwrap()
            .state,
            FaceState::NotOnFace
        );
    }
}

#[test]
fn distances_use_finite_bore_walls_and_caps_with_holes_not_supports() {
    let plan = compile(&document());
    let outer = origin(json!({"step":"bored","role":{"hole":{"at":"bore","part":"wall"}}}));
    let floor =
        origin(json!({"step":"shell","role":{"inner":{"of":{"step":"stock","role":"start"}}}}));
    let bottom = origin(json!({"step":"stock","role":"start"}));
    let controls = Options {
        width: Rational::parse("1/100").unwrap(),
        max_cells: 16384,
    };
    let mut dec = Decider::default();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    for (name, xyz, expected, max_cells) in [
        (
            inner_bore(&plan).origin.clone(),
            ["26", "15", "6"],
            Real::one(),
            16384,
        ),
        (
            inner_bore(&plan).origin.clone(),
            ["25", "15", "13"],
            Real::one(),
            16384,
        ),
        (outer, ["24", "15", "6"], Real::one(), 16384),
        (floor, ["20", "15", "3"], Real::from_i64(26).sqrt(), 1),
        (bottom, ["20", "15", "-1"], Real::from_i64(10).sqrt(), 1),
    ] {
        let controls = Options {
            max_cells,
            width: controls.width.clone(),
        };
        let (answer, work) = with_operation_budget(WORK, || {
            prepared
                .point_to_face(
                    plan.select(&name).unwrap(),
                    &point(xyz).to_real(),
                    &controls,
                    &mut dec,
                )
                .unwrap()
        });
        assert!(!work.exhausted, "{name:?}: {work:?}");
        let expected = dec.enclose_at(&expected, 128).unwrap();
        let upper = answer.upper.as_ref().expect("finite physical witness");
        assert!(
            answer.lower <= expected.lo() && upper >= &expected.hi(),
            "{name:?}: {answer:?}"
        );
        assert!(
            answer.converged && upper - &answer.lower <= controls.width,
            "{name:?}: {answer:?}"
        );
    }
}

#[test]
fn closed_bottom_open_both_open_and_consumed_cavities_follow_the_definition() {
    let mut source = document();
    source["steps"][2]["open"] = json!([]);
    material(
        &compile(&source),
        &[
            (["10", "10", "11"], Location::Inside),
            (["10", "10", "6"], Location::Outside),
            (["24", "15", "6"], Location::Inside),
        ],
    );
    source["steps"][2]["open"] = json!([{"face":{"step":"stock","role":"start"}}]);
    material(
        &compile(&source),
        &[
            (["10", "10", "1"], Location::Outside),
            (["10", "10", "11"], Location::Inside),
        ],
    );
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}}, {"face":{"step":"stock","role":"end"}}
    ]);
    material(
        &compile(&source),
        &[
            (["10", "10", "1"], Location::Outside),
            (["24", "15", "1"], Location::Inside),
            (["10", "10", "11"], Location::Outside),
        ],
    );
    source["steps"][2]["open"] = json!([]);
    source["steps"][2]["thickness"] = json!("6");
    let plan = compile(&source);
    material(
        &plan,
        &[
            (["10", "10", "6"], Location::Inside),
            (["20", "15", "6"], Location::Outside),
        ],
    );
    assert!(
        !plan
            .faces()
            .iter()
            .any(|face| matches!(face.origin.role, Role::Inner { .. }))
    );
}

#[test]
fn multiple_bores_keep_distinct_inner_names_and_material_columns() {
    let mut source = document();
    source["steps"][1]["at"] = json!([
        {"id":"left","point":["12","15","12"],"direction":["0","0","-1"]},
        {"id":"right","point":["28","15","12"],"direction":["0","0","-1"]}
    ]);
    let plan = compile(&source);
    material(
        &plan,
        &[
            (["12", "15", "6"], Location::Outside),
            (["28", "15", "6"], Location::Outside),
            (["16", "15", "6"], Location::Inside),
            (["24", "15", "6"], Location::Inside),
            (["20", "15", "6"], Location::Outside),
        ],
    );
    for (at, own, other) in [("left", "17", "33"), ("right", "33", "17")] {
        let name = origin(json!({"step":"shell","role":{"inner":{"of":{
            "step":"bored","role":{"hole":{"at":at,"part":"wall"}}
        }}}}));
        let id = plan.select(&name).unwrap();
        for (x, want) in [(own, FaceState::OnFace), (other, FaceState::NotOnFace)] {
            assert_eq!(
                plan.face_query(id, &point([x, "15", "6"]), &mut Workspace::default())
                    .unwrap()
                    .state,
                want
            );
        }
    }
}

#[test]
fn rational_rigid_placement_keeps_the_same_material_and_finite_bore() {
    let mut source = document();
    source["steps"][0]["pose"] =
        json!({"rotation":["3/5","4/5","0","0"],"translation":["7","-5","11"]});
    source["steps"][1]["at"][0]["point"] = json!(["27", "-518/25", "551/25"]);
    source["steps"][1]["at"][0]["direction"] = json!(["0", "24/25", "7/25"]);
    let plan = compile(&source);
    material(
        &plan,
        &[
            (["27", "-374/25", "593/25"], Location::Outside),
            (["31", "-374/25", "593/25"], Location::Inside),
            (["33", "-374/25", "593/25"], Location::Outside),
        ],
    );
    assert_eq!(
        plan.face_query(
            plan.select(&inner_bore(&plan).origin).unwrap(),
            &point(["32", "-374/25", "593/25"]),
            &mut Workspace::default()
        )
        .unwrap()
        .state,
        FaceState::OnFace
    );
}

#[test]
fn unsupported_openings_and_touching_or_crossing_envelopes_still_refuse() {
    let unsupported = |source: Value| {
        assert_eq!(
            Plan::compile(&recipe(&source)).unwrap_err().code,
            "distance.unsupported-envelope"
        );
    };
    for open in [
        json!([{"face":{"step":"bored"}}]),
        json!([{"face":{"step":"stock"}}]),
        json!([{"face":{"step":"stock"},"near":["0","0","12"]}]),
        json!([{"face":{"step":"stock","role":"end"},"touching":[{"step":"bored"}]}]),
    ] {
        let mut source = document();
        source["steps"][2]["open"] = open;
        unsupported(source);
    }
    for x in ["7", "6"] {
        let mut source = document();
        source["steps"][1]["at"][0]["point"] = json!([x, "15", "12"]);
        unsupported(source);
    }
    for right in ["25", "24"] {
        let mut source = document();
        source["steps"][1]["at"] = json!([
            {"id":"left","point":["15","15","12"],"direction":["0","0","-1"]},
            {"id":"right","point":[right,"15","12"],"direction":["0","0","-1"]}
        ]);
        unsupported(source);
    }
    let mut blind = document();
    blind["steps"][1]["depth"] = json!({"blind":"6"});
    unsupported(blind);
}

#[test]
fn boolean_cylinder_tools_and_circular_through_cuts_keep_their_own_inner_names() {
    let base = document()["steps"][0].clone();
    let shell = document()["steps"][2].clone();
    let tool = json!({"id":"tool","op":"cylinder","radius":"3","height":"14",
        "pose":{"rotation":["1","0","0","0"],"translation":["20","15","-1"]}});
    let profile = json!({"id":"profile","op":"sketch",
        "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
        "loops":[[{"kind":"circle","id":"circle","centre":["20","15"],"radius":"3"}]]});
    for (steps, wall) in [
        (
            vec![
                base.clone(),
                tool,
                json!({"id":"bored","op":"boolean","kind":"subtract","target":"stock","tool":"tool"}),
                shell.clone(),
            ],
            json!({"step":"tool","role":{"side":{"segment":"wall"}}}),
        ),
        (
            vec![
                base,
                profile,
                json!({"id":"bored","op":"extrude_through","target":"stock","profile":"profile"}),
                shell,
            ],
            json!({"step":"bored","role":{"side":{"segment":"circle"}}}),
        ),
    ] {
        let plan = compile(
            &json!({"format":"catallax-recipe/4","name":"tool shell","steps":steps,"result":"shell"}),
        );
        material(
            &plan,
            &[
                (["20", "15", "6"], Location::Outside),
                (["24", "15", "6"], Location::Inside),
                (["26", "15", "6"], Location::Outside),
            ],
        );
        let inner = origin(json!({"step":"shell","role":{"inner":{"of":wall}}}));
        assert_eq!(
            plan.face_query(
                plan.select(&inner).unwrap(),
                &point(["25", "15", "6"]),
                &mut Workspace::default()
            )
            .unwrap()
            .state,
            FaceState::OnFace
        );
        cap_distance(&plan, ["20", "15", "3"], Real::from_i64(26).sqrt(), 1);
    }
}

fn cap_distance(plan: &Plan, xyz: [&str; 3], expected: Real, max_cells: usize) {
    let floor =
        origin(json!({"step":"shell","role":{"inner":{"of":{"step":"stock","role":"start"}}}}));
    let controls = Options {
        width: Rational::parse("1/100").unwrap(),
        max_cells,
    };
    let mut dec = Decider::default();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    let (answer, work) = with_operation_budget(WORK, || {
        prepared
            .point_to_face(
                plan.select(&floor).unwrap(),
                &point(xyz).to_real(),
                &controls,
                &mut dec,
            )
            .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let expected = dec.enclose_at(&expected, 128).unwrap();
    let upper = answer.upper.as_ref().expect("admitted finite-cap witness");
    assert!(
        answer.lower <= expected.lo() && upper >= &expected.hi(),
        "{answer:?}"
    );
    assert!(
        answer.converged && upper - &answer.lower <= controls.width,
        "{answer:?}"
    );
}

#[test]
fn circular_cap_bounds_cover_edits_multiple_holes_rotations_and_small_clearances() {
    let mut source = document();
    source["steps"][2]["thickness"] = json!("3/2");
    cap_distance(
        &compile(&source),
        ["21", "15", "5/2"],
        Real::from(Rational::parse("53/4").unwrap()).sqrt(),
        1,
    );

    source = document();
    source["steps"][1]["at"] = json!([
        {"id":"left","point":["12","15","12"],"direction":["0","0","-1"]},
        {"id":"right","point":["28","15","12"],"direction":["0","0","-1"]}
    ]);
    let plan = compile(&source);
    cap_distance(&plan, ["12", "15", "3"], Real::from_i64(26).sqrt(), 1);
    cap_distance(&plan, ["29", "15", "3"], Real::from_i64(17).sqrt(), 1);

    source = document();
    source["steps"][0]["pose"] =
        json!({"rotation":["3/5","4/5","0","0"],"translation":["7","-5","11"]});
    source["steps"][1]["at"][0]["point"] = json!(["27", "-518/25", "551/25"]);
    source["steps"][1]["at"][0]["direction"] = json!(["0", "24/25", "7/25"]);
    cap_distance(
        &compile(&source),
        ["27", "-302/25", "614/25"],
        Real::from_i64(26).sqrt(),
        1,
    );

    source = document();
    source["steps"][1]["at"][0]["point"] = json!(["700001/100000", "15", "12"]);
    cap_distance(
        &compile(&source),
        ["600001/100000", "15", "3"],
        Real::from_i64(17).sqrt(),
        32,
    );
}
