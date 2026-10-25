//! Single-plane thickening of nonconvex faces: exact prisms, rim origins and both kernels.

use catallax_features::{Built, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use catallax_topology::BodyKind;
use serde_json::{Value, json};

fn pose() -> Value {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn polygon(points: &[[i64; 2]], prefix: &str) -> Value {
    Value::Array(
        points
            .iter()
            .enumerate()
            .map(|(i, p)| {
                json!({"kind":"line","id":format!("{prefix}{i}"),"to":[p[0].to_string(), p[1].to_string()]})
            })
            .collect(),
    )
}
/// An L with arms of width 2 and length 4, 3 tall.
fn l_shape() -> Vec<Value> {
    vec![
        json!({"id":"sketch","op":"sketch","plane":pose(),"loops":[polygon(&[[4,0],[4,2],[2,2],[2,4],[0,4],[0,0]], "l")]}),
        json!({"id":"stock","op":"extrude","profile":"sketch","distance":"3"}),
    ]
}
/// A 6 by 6 plate with a 2 by 2 hole, 3 tall.
fn holed() -> Vec<Value> {
    vec![
        json!({"id":"sketch","op":"sketch","plane":pose(),"loops":[
            polygon(&[[6,0],[6,6],[0,6],[0,0]], "o"),
            polygon(&[[2,2],[4,2],[4,4],[2,4]], "h")
        ]}),
        json!({"id":"stock","op":"extrude","profile":"sketch","distance":"3"}),
    ]
}
fn recipe(mut steps: Vec<Value>, role: Value, thickness: &str, side: &str) -> Recipe {
    steps.push(json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":role}}],"thickness":thickness,"side":side}));
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"nonconvex thicken","steps":steps,"result":"result"})
            .to_string(),
    )
    .unwrap()
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut Decider::default()).unwrap();
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    built
}
fn volume(built: &Built, expected: &str) {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let expected = Rational::parse(expected).unwrap();
    assert_eq!(mass.volume.lo, expected);
    assert_eq!(mass.volume.hi, expected);
}
fn roles(built: &Built) -> Vec<Role> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).to_vec())
        .map(|o| o.role)
        .collect()
}
/// Checks the same points in the B-rep and the independently compiled field.
fn agree(recipe: &Recipe, points: &[([&str; 3], Location)]) {
    let mut built = build(recipe);
    let plan = Plan::compile(recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, want) in points {
        let [x, y, z] = point.map(|s| Rational::parse(s).unwrap());
        let p = RPoint3::new(x, y, z);
        let field = plan
            .proven_checked(&p, &mut workspace, &mut Decider::default())
            .unwrap()
            .location;
        let brep = match catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &p,
            &mut Decider::default(),
        ) {
            Outcome::Proven(catallax_boolean::Containment::Inside) => Location::Inside,
            Outcome::Proven(catallax_boolean::Containment::Outside) => Location::Outside,
            Outcome::ProvenCoincident(_) => Location::Boundary,
            other => panic!("unresolved {point:?}: {other:?}"),
        };
        assert_eq!(field, *want, "field at {point:?}");
        assert_eq!(brep, *want, "brep at {point:?}");
    }
}

#[test]
fn l_shaped_end_face_thickens_front_and_back_with_exact_prism_sides() {
    for (side, z_in, z_out) in [("front", "7/2", "9/2"), ("back", "5/2", "3/2")] {
        let recipe = recipe(l_shape(), json!("end"), "1", side);
        let built = build(&recipe);
        volume(&built, "12");
        let roles = roles(&built);
        let rims = roles
            .iter()
            .filter(|r| matches!(r, Role::Rim { .. }))
            .count();
        assert_eq!(rims, 6, "one prism side per L edge and no partition faces");
        assert_eq!(built.model.faces_of_body(built.result).len(), 8);
        assert!(roles.contains(&Role::End));
        assert!(roles.iter().any(|r| matches!(r, Role::Offset { .. })));
        let (low, high) = if side == "front" {
            ("3", "4")
        } else {
            ("2", "3")
        };
        agree(
            &recipe,
            &[
                (["1", "1", z_in], Location::Inside),
                (["3", "1", z_in], Location::Inside),
                (["1", "3", z_in], Location::Inside),
                (["3", "3", z_in], Location::Outside),
                (["1", "1", z_out], Location::Outside),
                (["3", "1", low], Location::Boundary),
                (["3", "1", high], Location::Boundary),
                (["2", "3", z_in], Location::Boundary),
                (["5", "1", z_in], Location::Outside),
            ],
        );
    }
}

#[test]
fn start_face_and_a_face_with_a_hole_thicken_with_straight_hole_walls() {
    let recipe = recipe(holed(), json!("start"), "1", "front");
    let built = build(&recipe);
    volume(&built, "32");
    let rims = roles(&built)
        .iter()
        .filter(|r| matches!(r, Role::Rim { .. }))
        .count();
    assert_eq!(rims, 8, "four outer and four hole walls");
    agree(
        &recipe,
        &[
            (["1", "1", "-1/2"], Location::Inside),
            (["3", "3", "-1/2"], Location::Outside),
            (["3", "1", "-1/2"], Location::Inside),
            (["2", "3", "-1/2"], Location::Boundary),
            (["3", "3", "1"], Location::Outside),
            (["1", "1", "1/2"], Location::Outside),
        ],
    );
}

#[test]
fn two_sided_union_of_a_nonconvex_face_removes_the_shared_base() {
    let mut steps = l_shape();
    for side in ["front", "back"] {
        steps.push(json!({"id":side,"op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":"end"}}],"thickness":"1","side":side}));
    }
    steps.push(json!({"id":"result","op":"boolean","kind":"union","target":"front","tool":"back"}));
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"symmetric","steps":steps,"result":"result"})
            .to_string(),
    )
    .unwrap();
    let built = build(&recipe);
    volume(&built, "24");
    agree(
        &recipe,
        &[
            (["1", "1", "3"], Location::Inside),
            (["1", "1", "7/2"], Location::Inside),
            (["1", "1", "5/2"], Location::Inside),
            (["3", "3", "3"], Location::Outside),
            (["1", "1", "9/2"], Location::Outside),
        ],
    );
}

#[test]
fn reversed_face_sense_sweeps_toward_the_same_outward_normal() {
    // A downward extrusion reverses the sense of both caps.
    let steps = vec![
        l_shape().remove(0),
        json!({"id":"stock","op":"extrude","profile":"sketch","distance":"-3"}),
    ];
    for (side, z_in, z_out) in [("front", "-7/2", "-9/2"), ("back", "-5/2", "-3/2")] {
        let recipe = recipe(steps.clone(), json!("end"), "1", side);
        volume(&build(&recipe), "12");
        agree(
            &recipe,
            &[
                (["1", "1", z_in], Location::Inside),
                (["1", "3", z_in], Location::Inside),
                (["3", "3", z_in], Location::Outside),
                (["1", "1", z_out], Location::Outside),
            ],
        );
    }
}

#[test]
fn narrow_valid_regions_keep_their_exact_prisms() {
    let narrow = vec![
        json!({"id":"sketch","op":"sketch","plane":pose(),"loops":[polygon(&[[100,0],[100,1],[1,1],[1,100],[0,100],[0,0]], "n")]}),
        json!({"id":"stock","op":"extrude","profile":"sketch","distance":"3"}),
    ];
    let recipe = recipe(narrow, json!("end"), "1/1000", "front");
    let built = build(&recipe);
    volume(&built, "199/1000");
    assert_eq!(built.model.body(built.result).kind, BodyKind::Solid);
    agree(
        &recipe,
        &[
            (["50", "1/2", "6001/2000"], Location::Inside),
            (["1/2", "50", "6001/2000"], Location::Inside),
            (["50", "50", "6001/2000"], Location::Outside),
            (["50", "1/2", "3001/1000"], Location::Boundary),
            (["50", "1/2", "3002/1000"], Location::Outside),
        ],
    );
}

fn result_origins(plan: &Plan) -> std::collections::BTreeSet<catallax_recipe::Origin> {
    (0..plan.faces().len())
        .filter(|i| plan.result_face(*i).is_some())
        .map(|i| plan.faces()[i].origin.clone())
        .collect()
}
fn brep_origins(built: &Built) -> std::collections::BTreeSet<catallax_recipe::Origin> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).to_vec())
        .filter(|o| o.step.as_str() == "result" || o.step.as_str() == "stock")
        .collect()
}

#[test]
fn both_kernels_name_the_same_rims_for_every_selection() {
    for (steps, role) in [(l_shape(), "end"), (l_shape(), "start"), (holed(), "end")] {
        for side in ["front", "back"] {
            let recipe = recipe(steps.clone(), json!(role), "1/2", side);
            let built = build(&recipe);
            let plan = Plan::compile(&recipe).unwrap();
            let field: std::collections::BTreeSet<_> = result_origins(&plan)
                .into_iter()
                .filter(|o| o.step.as_str() == "result")
                .collect();
            let brep: std::collections::BTreeSet<_> = brep_origins(&built)
                .into_iter()
                .filter(|o| o.step.as_str() == "result")
                .collect();
            assert!(!brep.is_empty());
            assert_eq!(field, brep, "{role} {side}");
        }
    }
}

#[test]
fn display_follows_the_notch_and_step_export_succeeds() {
    let recipe = recipe(l_shape(), json!("end"), "1", "front");
    let plan = Plan::compile(&recipe).unwrap();
    let mesh = catallax_field::display::mesh_plan(
        &plan,
        &catallax_field::display::Options {
            tolerance: Rational::from_ratio(1, 10).unwrap(),
            ..Default::default()
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(!mesh.triangles.is_empty());
    for p in &mesh.positions {
        assert!((-0.1..=4.1).contains(&p[0]) && (-0.1..=4.1).contains(&p[1]));
        assert!(
            (2.9..=4.1).contains(&p[2]),
            "only the new prism is displayed"
        );
    }
    for t in &mesh.triangles {
        let c: Vec<f64> = (0..3)
            .map(|k| {
                t.vertices
                    .iter()
                    .map(|v| mesh.positions[*v][k])
                    .sum::<f64>()
                    / 3.0
            })
            .collect();
        assert!(
            !(c[0] > 2.1 && c[1] > 2.1),
            "no hull is invented over the notch: {c:?}"
        );
    }
    let built = build(&recipe);
    let mut out = Vec::new();
    let report = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut out,
    )
    .unwrap();
    assert!(report.marks.is_empty());
    assert!(!out.is_empty());
}

#[test]
fn editing_the_thickness_replays_to_the_fresh_build() {
    let recipe = recipe(l_shape(), json!("end"), "1", "front");
    let checkpoint =
        catallax_features::evaluate_incremental(&recipe, None, &mut Decider::default()).unwrap();
    let mut edited = recipe.clone();
    if let catallax_recipe::Operation::Thicken {
        thickness, side, ..
    } = &mut edited.steps.last_mut().unwrap().op
    {
        *thickness = catallax_recipe::Number::new(Rational::parse("3/2").unwrap());
        *side = catallax_recipe::ThickenSide::Back;
    } else {
        panic!("thicken step");
    }
    let resumed =
        catallax_features::evaluate_incremental(&edited, Some(checkpoint), &mut Decider::default())
            .unwrap();
    let resumed = resumed.built();
    let fresh = build(&edited);
    volume(&fresh, "18");
    assert_eq!(
        resumed.model.counts(resumed.result),
        fresh.model.counts(fresh.result)
    );
    for (a, b) in resumed
        .model
        .faces_of_body(resumed.result)
        .into_iter()
        .zip(fresh.model.faces_of_body(fresh.result))
    {
        assert_eq!(resumed.model.face_origins(a), fresh.model.face_origins(b));
    }
}

#[test]
fn multi_face_and_side_selections_of_nonconvex_bodies_still_refuse() {
    for faces in [
        json!([{"face":{"step":"stock","role":"end"}},{"face":{"step":"stock","role":"start"}}]),
        json!([{"face":{"step":"stock","role":{"side":{"segment":"l0"}}}}]),
        json!([{"face":{"step":"stock"}}]),
    ] {
        let mut steps = l_shape();
        steps.push(json!({"id":"result","op":"thicken","body":"stock","faces":faces,"thickness":"1","side":"front"}));
        let recipe = Recipe::parse(
            &json!({"format":"catallax-recipe/4","name":"refusal","steps":steps,"result":"result"})
                .to_string(),
        )
        .unwrap();
        let brep = evaluate(&recipe, &mut Decider::default()).err().unwrap();
        let field = Plan::compile(&recipe).err().unwrap();
        assert!(brep.to_string().contains("distance."), "{brep}");
        assert!(field.code.starts_with("distance."), "{field:?}");
    }
}
