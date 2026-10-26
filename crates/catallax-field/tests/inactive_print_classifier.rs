//! Retained source definitions must not change live print or mass consumers.
use catallax_field::{FaceState, Location, Plan, Workspace, mass, printing};
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn q(value: &str) -> Rational {
    value.parse().unwrap()
}
fn pose(x: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,"0","0"]})
}
fn feature(id: &str) -> Vec<Value> {
    let profile = format!("{id}_profile");
    let x = match id {
        "unused_between" => "2",
        "unused_after" => "3",
        _ => "1",
    };
    vec![
        json!({"id":profile,"op":"sketch","plane":pose(x),"loops":[[
            {"kind":"line","id":"a","to":["1","0"]},
            {"kind":"line","id":"b","to":["0","1"]},
            {"kind":"line","id":"c","to":["0","0"]}
        ]]}),
        json!({"id":id,"op":"extrude","profile":profile,"distance":"1"}),
    ]
}
fn recipe(kind: &str, retained: bool) -> Recipe {
    let mut steps = Vec::new();
    if retained {
        steps.extend(feature("unused_before"));
    }
    steps.push(json!({"id":"left","op":"block","size":["1/64","1/64","1/64"],"pose":pose("0")}));
    if retained {
        steps.extend(feature("unused_between"));
    }
    steps.push(
        json!({"id":"right","op":"block","size":["1/64","1/64","1/64"],"pose":pose("1/128")}),
    );
    steps.push(json!({"id":"result","op":"boolean","kind":kind,"target":"left","tool":"right"}));
    if retained {
        steps.extend(feature("unused_after"));
    }
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"retained print sources","steps":steps,"result":"result"}).to_string()).unwrap()
}
fn plan(kind: &str, retained: bool) -> Plan {
    Plan::compile(&recipe(kind, retained)).unwrap()
}
fn print(plan: &Plan) -> (Result<printing::Mesh, Refusal>, u64) {
    let (result, work) = with_operation_budget(100_000_000, || {
        printing::mesh_plan(plan, &printing::Options::default(), &mut Decider::default())
    });
    assert!(!work.exhausted);
    (result, work.used)
}
fn bytes(mesh: &printing::Mesh) -> (Vec<u8>, Vec<u8>) {
    let (mut stl, mut obj) = (Vec::new(), Vec::new());
    mesh.write_stl(&mut stl).unwrap();
    mesh.write_obj(&mut obj).unwrap();
    (stl, obj)
}

#[test]
fn inactive_features_preserve_print_bytes_work_and_boolean_material() {
    for kind in ["union", "intersect", "subtract"] {
        let compact = plan(kind, false);
        let retained = plan(kind, true);
        let (expected, expected_work) = print(&compact);
        let (actual, actual_work) = print(&retained);
        match (actual, expected) {
            (Ok(actual), Ok(expected)) => {
                assert!(!actual.triangles.is_empty());
                assert_eq!(bytes(&actual), bytes(&expected), "{kind}");
                assert_eq!(actual.cells, expected.cells);
                assert!(actual.max_deviation <= printing::Options::default().tolerance);
            }
            (Err(actual), Err(expected)) => {
                assert_eq!(actual.code, expected.code, "{kind}");
                assert_eq!(actual.detail, expected.detail, "{kind}");
            }
            (actual, expected) => panic!("{kind}: retained {actual:?}; compact {expected:?}"),
        }
        assert_eq!(actual_work, expected_work);
        // Exact points cover the overlap, both ends, ties and exterior.
        for (x, union, intersect, subtract) in [
            (
                "1/256",
                Location::Inside,
                Location::Outside,
                Location::Inside,
            ),
            (
                "1/128",
                Location::Inside,
                Location::Boundary,
                Location::Boundary,
            ),
            (
                "3/256",
                Location::Inside,
                Location::Inside,
                Location::Outside,
            ),
            (
                "5/256",
                Location::Inside,
                Location::Outside,
                Location::Outside,
            ),
            (
                "1/32",
                Location::Outside,
                Location::Outside,
                Location::Outside,
            ),
        ] {
            let point = RPoint3::new(q(x), q("1/128"), q("1/128"));
            let expected = match kind {
                "union" => union,
                "intersect" => intersect,
                _ => subtract,
            };
            assert_eq!(
                retained.proven(&point, &mut Workspace::default()).location,
                expected
            );
            assert_eq!(
                compact.proven(&point, &mut Workspace::default()).location,
                expected
            );
        }
        let origin = Origin {
            step: StepId::new("left").unwrap(),
            role: Role::Start,
        };
        for (x, state) in [("3/256", FaceState::OnFace), ("1/32", FaceState::NotOnFace)] {
            // The subtract result has no material at the overlap point.
            let x = if kind == "subtract" && state == FaceState::OnFace {
                "1/256"
            } else {
                x
            };
            let point = RPoint3::new(q(x), q("1/128"), q("0"));
            for plan in [&compact, &retained] {
                assert_eq!(
                    plan.face_query(
                        plan.select(&origin).unwrap(),
                        &point,
                        &mut Workspace::default()
                    )
                    .unwrap()
                    .state,
                    state
                );
            }
        }
    }
}

#[test]
fn inactive_features_preserve_an_actual_printable_output() {
    let compile = |retained| {
        let mut recipe = recipe("union", retained);
        recipe.result = StepId::new("left").unwrap();
        Plan::compile(&recipe).unwrap()
    };
    let (expected, expected_work) = print(&compile(false));
    let (actual, actual_work) = print(&compile(true));
    let (expected, actual) = (expected.unwrap(), actual.unwrap());
    assert!(!actual.triangles.is_empty());
    assert_eq!(bytes(&actual), bytes(&expected));
    assert_eq!(actual.cells, expected.cells);
    assert_eq!(actual_work, expected_work);
    assert!(actual.max_deviation <= printing::Options::default().tolerance);
}

#[test]
fn ordinary_mass_ignores_inactive_features_at_default_accuracy_and_caps() {
    for (kind, cells) in [("union", 12), ("intersect", 4), ("subtract", 4)] {
        let measure = |plan: &Plan| {
            let (result, work) = with_operation_budget(100_000_000, || {
                mass::properties(
                    plan,
                    &mass::Accuracy::default(),
                    &mass::Options::default(),
                    &mut Decider::default(),
                )
            });
            assert!(!work.exhausted);
            (result.unwrap(), work.used)
        };
        let (expected, expected_work) = measure(&plan(kind, false));
        let (actual, actual_work) = measure(&plan(kind, true));
        // Independent inclusion/exclusion of two overlapping 2s by 2s by 2s boxes, s=1/128.
        let volume = Rational::from_i64(cells)
            .checked_div(&q("2097152"))
            .unwrap();
        assert!(actual.volume.is_exact());
        assert_eq!(actual.volume.lo, volume);
        assert_eq!(actual.volume, expected.volume);
        assert_eq!(actual.centroid, expected.centroid);
        assert_eq!(actual.inertia, expected.inertia);
        assert_eq!(actual.cells, expected.cells);
        assert_eq!(actual_work, expected_work);
    }
}

#[test]
fn an_active_unsupported_feature_still_refuses_print_and_mass() {
    let recipe = Recipe::parse(&json!({"format":"catallax-recipe/3","name":"active triangle","steps":feature("triangle"),"result":"triangle"}).to_string()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    assert_eq!(
        printing::mesh_plan(
            &plan,
            &printing::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.print-coverage"
    );
    assert_eq!(
        mass::properties(
            &plan,
            &mass::Accuracy::default(),
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.print-coverage"
    );
}

fn handed_off_plan() -> Plan {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/boolean-face-thicken/collinear-neighbours.json"
    ))
    .unwrap();
    Plan::compile(&recipe).unwrap()
}

#[test]
fn handed_off_triangle_history_has_exact_box_mass_at_original_defaults() {
    // Preserve the original ordinary library call: no operation-budget wrapper.
    let result = mass::properties(
        &handed_off_plan(),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(result.volume.is_exact());
    assert_eq!(result.volume.lo, q("16"));
    // Independent moments of [0,4] x [0,4] x [3,4], not a mesh or B-rep query.
    for (value, expected) in result
        .centroid
        .as_ref()
        .unwrap()
        .iter()
        .zip(["2", "2", "7/2"])
    {
        assert!(value.is_exact());
        assert_eq!(value.lo, q(expected));
    }
    for (i, row) in result.inertia.as_ref().unwrap().iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert!(value.is_exact());
            assert_eq!(
                value.lo,
                if i == j {
                    q(["68/3", "68/3", "128/3"][i])
                } else {
                    q("0")
                }
            );
        }
    }
}
