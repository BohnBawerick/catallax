//! Public synthetic planar Boolean faces, checked against independent exact oracles.

use catallax_features::{Built, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role};
use catallax_shapes::{Curve, RPoint3, Surface};
use serde_json::{Value, json};
use std::collections::BTreeSet;

const L_FRONT: &str = include_str!("fixtures/boolean-face-thicken/l-front.json");
const HOLE_FRONT: &str = include_str!("fixtures/boolean-face-thicken/hole-front.json");
const SPLIT: &str = include_str!("fixtures/boolean-face-thicken/split-near.json");

fn parse(text: &str) -> Recipe {
    Recipe::parse(text).unwrap()
}
fn edited(text: &str, edit: impl FnOnce(&mut Value)) -> Recipe {
    let mut value: Value = serde_json::from_str(text).unwrap();
    edit(&mut value);
    parse(&value.to_string())
}
fn build(recipe: &Recipe) -> Built {
    let (built, work) =
        with_operation_budget(100_000_000, || evaluate(recipe, &mut Decider::default()));
    let built = built.unwrap();
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    built
}
fn names(built: &Built) -> BTreeSet<Origin> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).to_vec())
        .collect()
}
fn exact_volume(built: &Built, want: &str) {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(mass.volume.lo, Rational::parse(want).unwrap());
    assert_eq!(mass.volume.lo, mass.volume.hi);
}
fn agree(recipe: &Recipe, points: &[([&str; 3], Location)]) {
    let mut built = build(recipe);
    let (plan, work) = with_operation_budget(100_000_000, || Plan::compile(recipe));
    let plan = plan.unwrap();
    assert!(!work.exhausted);
    let field: BTreeSet<_> = (0..plan.faces().len())
        .filter(|i| plan.result_face(*i).is_some())
        .map(|i| plan.faces()[i].origin.clone())
        .collect();
    let expected = names(&built);
    if recipe.name == "Two sided Boolean L" {
        // The registry contains candidates. Query physical base survival.
        let extra: BTreeSet<_> = field.difference(&expected).cloned().collect();
        assert_eq!(extra.len(), 1);
        assert!(
            extra
                .iter()
                .all(|o| o.step.as_str() == "stock" && o.role == Role::End)
        );
        for i in 0..plan.faces().len() {
            if !extra.contains(&plan.faces()[i].origin) {
                continue;
            }
            if let Some(id) = plan.result_face(i) {
                for p in [
                    RPoint3::new(1.into(), 1.into(), 3.into()),
                    RPoint3::new(3.into(), 1.into(), 3.into()),
                ] {
                    assert_eq!(
                        plan.face_query_checked(
                            id,
                            &p,
                            &mut Workspace::default(),
                            &mut Decider::default()
                        )
                        .unwrap()
                        .state,
                        catallax_field::FaceState::NotOnFace
                    );
                }
            }
        }
        assert!(expected.is_subset(&field));
    } else {
        assert_eq!(field, expected, "complete source, offset and rim origins");
    }
    for (point, expected) in points {
        let [x, y, z] = point.map(|s| Rational::parse(s).unwrap());
        let p = RPoint3::new(x, y, z);
        let actual = match catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &p,
            &mut Decider::default(),
        ) {
            Outcome::Proven(catallax_boolean::Containment::Inside) => Location::Inside,
            Outcome::Proven(catallax_boolean::Containment::Outside) => Location::Outside,
            Outcome::ProvenCoincident(_) => Location::Boundary,
            other => panic!("{point:?}: {other:?}"),
        };
        assert_eq!(actual, *expected, "B-rep at {point:?}");
        assert_eq!(
            plan.proven_checked(&p, &mut Workspace::default(), &mut Decider::default())
                .unwrap()
                .location,
            *expected,
            "field at {point:?}"
        );
    }
}

#[test]
fn retained_l_and_through_hole_have_exact_straight_prisms_front_and_back() {
    for (text, volume, faces, hole) in [
        (L_FRONT, "12", 8, false),
        (
            include_str!("fixtures/boolean-face-thicken/l-back.json"),
            "12",
            8,
            false,
        ),
        (HOLE_FRONT, "32", 10, true),
        (
            include_str!("fixtures/boolean-face-thicken/hole-back.json"),
            "32",
            10,
            true,
        ),
    ] {
        let recipe = parse(text);
        let built = build(&recipe);
        exact_volume(&built, volume);
        assert_eq!(
            built.model.faces_of_body(built.result).len(),
            faces,
            "no partition faces"
        );
        assert!(
            built
                .model
                .edges_of_body(built.result)
                .iter()
                .all(|e| matches!(
                    built.model.geom().curve(built.model.edge(*e).curve),
                    Curve::Line(_)
                ))
        );
        assert!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .all(|f| matches!(
                    built.model.geom().surface(built.model.face(*f).surface),
                    Surface::Plane(_)
                ))
        );
        assert_eq!(
            names(&built)
                .iter()
                .filter(|o| matches!(o.role, Role::Rim { .. }))
                .count(),
            faces - 2
        );
        let inside = if recipe.name.ends_with("front") {
            "7/2"
        } else {
            "5/2"
        };
        let outside = if recipe.name.ends_with("front") {
            "5/2"
        } else {
            "7/2"
        };
        let mut probes = vec![
            (["1", "1", inside], Location::Inside),
            (["2", "1", inside], Location::Inside),
            (["3", "3", inside], Location::Outside),
            (["1", "1", outside], Location::Outside),
            (["4", "4", inside], Location::Outside),
        ];
        if hole {
            probes.pop();
            probes.push((["5", "5", inside], Location::Inside));
            probes.push((["2", "3", inside], Location::Boundary));
        }
        agree(&recipe, &probes);
    }
}

#[test]
fn explicit_two_sided_union_removes_its_internal_base() {
    let recipe = parse(include_str!("fixtures/boolean-face-thicken/symmetric.json"));
    let built = build(&recipe);
    exact_volume(&built, "24");
    assert_eq!(built.model.faces_of_body(built.result).len(), 8);
    assert!(names(&built).iter().all(|o| o.step.as_str() != "stock"));
    agree(
        &recipe,
        &[
            (["1", "1", "3"], Location::Inside),
            (["1", "1", "7/2"], Location::Inside),
            (["1", "1", "5/2"], Location::Inside),
            (["3", "3", "3"], Location::Outside),
        ],
    );
}

#[test]
fn subtracted_tool_face_uses_its_reversed_material_sense() {
    let recipe = parse(include_str!(
        "fixtures/boolean-face-thicken/reversed-wall.json"
    ));
    exact_volume(&build(&recipe), "6");
    agree(
        &recipe,
        &[
            (["5/2", "3", "1"], Location::Inside),
            (["3/2", "3", "1"], Location::Outside),
            (["7/2", "3", "1"], Location::Outside),
            (["5/2", "1", "1"], Location::Outside),
        ],
    );
}

#[test]
fn closed_brep_display_preserves_the_notch_and_hole() {
    for (text, hole) in [(L_FRONT, false), (HOLE_FRONT, true)] {
        let built = build(&parse(text));
        let mesh = catallax_mesh::mesh_body(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        let mut uses = std::collections::BTreeMap::<_, (usize, i32)>::new();
        for t in &mesh.triangles {
            for i in 0..3 {
                let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
                let entry = uses.entry((a.min(b), a.max(b))).or_default();
                entry.0 += 1;
                entry.1 += if a < b { 1 } else { -1 };
            }
            let centre: [f64; 3] = std::array::from_fn(|axis| {
                t.vertices
                    .iter()
                    .map(|v| mesh.positions[*v][axis])
                    .sum::<f64>()
                    / 3.0
            });
            if hole {
                assert!(
                    !(centre[0] > 2.0 && centre[0] < 4.0 && centre[1] > 2.0 && centre[1] < 4.0)
                );
            } else {
                assert!(!(centre[0] > 2.0 && centre[1] > 2.0));
            }
        }
        assert!(!uses.is_empty());
        assert!(
            uses.values().all(|v| *v == (2, 0)),
            "closed, consistently oriented display"
        );
        let mut output = Vec::new();
        let report = catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut output,
        )
        .unwrap();
        assert!(report.marks.is_empty());
        assert!(!output.is_empty());
    }
}

#[test]
fn independently_integrated_field_volumes_contain_the_exact_oracles() {
    for (text, want) in [
        (L_FRONT, "12"),
        (HOLE_FRONT, "32"),
        (
            include_str!("fixtures/boolean-face-thicken/symmetric.json"),
            "24",
        ),
    ] {
        let plan = Plan::compile(&parse(text)).unwrap();
        let (mass, work) = with_operation_budget(100_000_000, || {
            catallax_field::mass::properties(
                &plan,
                &Default::default(),
                &Default::default(),
                &mut Decider::default(),
            )
        });
        let mass = mass.unwrap();
        assert!(!work.exhausted);
        assert_eq!(mass.volume.lo, Rational::parse(want).unwrap(), "{mass:?}");
        assert_eq!(
            mass.volume.lo, mass.volume.hi,
            "exact independent field volume"
        );
    }
}

#[test]
fn split_origin_near_pick_tracks_the_closed_trim_after_an_edit() {
    let original = parse(SPLIT);
    let checkpoint =
        catallax_features::evaluate_incremental(&original, None, &mut Decider::default()).unwrap();
    exact_volume(&build(&original), "8");
    agree(
        &original,
        &[
            (["1", "1", "7/2"], Location::Inside),
            (["5", "1", "7/2"], Location::Outside),
        ],
    );
    let changed = edited(SPLIT, |v| {
        v["steps"][1]["pose"]["translation"][0] = json!("1")
    });
    let resumed = catallax_features::evaluate_incremental(
        &changed,
        Some(checkpoint),
        &mut Decider::default(),
    )
    .unwrap();
    let fresh = build(&changed);
    exact_volume(&fresh, "4");
    assert_eq!(names(resumed.built()), names(&fresh));
    assert_eq!(names(&build(&original)), names(&fresh));
    agree(
        &changed,
        &[
            (["1/2", "1", "7/2"], Location::Inside),
            (["3/2", "1", "7/2"], Location::Outside),
            (["5", "1", "7/2"], Location::Outside),
        ],
    );
}

#[test]
fn lost_ambiguous_multiface_and_curved_requests_remain_typed_refusals() {
    for (recipe, expected) in [
        (
            edited(SPLIT, |v| {
                v["steps"][3]["faces"][0]["near"] = json!(["3", "1", "3"])
            }),
            "reference.ambiguous-face",
        ),
        (
            edited(SPLIT, |v| {
                v["steps"][3]["faces"][0]
                    .as_object_mut()
                    .unwrap()
                    .remove("near");
            }),
            "distance.",
        ),
        (
            edited(
                L_FRONT,
                |v| v["steps"][3]["faces"] = json!([{"face":{"step":"stock","role":"end"}},{"face":{"step":"stock","role":"start"}}]),
            ),
            "distance.",
        ),
    ] {
        let err = evaluate(&recipe, &mut Decider::default()).unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } = err else {
            panic!("geometry refusal")
        };
        let field = Plan::compile(&recipe).unwrap_err();
        assert!(refusal.code.starts_with(expected), "{refusal:?}");
        assert!(field.code.starts_with(expected), "{field:?}");
    }
    let lost = edited(L_FRONT, |v| {
        v["steps"][3]["faces"][0]["face"] = json!({"step":"tool","role":"end"})
    });
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&lost, &mut Decider::default()).unwrap_err()
    else {
        panic!("geometry refusal")
    };
    assert_eq!(refusal.code, "reference.lost-face");
    assert_eq!(
        Plan::compile(&lost).unwrap_err().code,
        "reference.lost-face"
    );
    let curved = edited(
        HOLE_FRONT,
        |v| v["steps"][1] = json!({"id":"tool","op":"cylinder","radius":"1","height":"5","pose":{"rotation":["1","0","0","0"],"translation":["3","3","-1"]}}),
    );
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&curved, &mut Decider::default()).unwrap_err()
    else {
        panic!("geometry refusal")
    };
    assert!(refusal.code.starts_with("distance."));
    assert!(
        Plan::compile(&curved)
            .unwrap_err()
            .code
            .starts_with("distance.")
    );
}

#[test]
fn merged_source_aliases_have_complete_offset_and_rim_products() {
    let recipe = parse(include_str!(
        "fixtures/boolean-face-thicken/l-merged-origins.json"
    ));
    let built = build(&recipe);
    exact_volume(&built, "12");
    assert_eq!(built.model.faces_of_body(built.result).len(), 8);
    for face in built.model.faces_of_body(built.result) {
        let origins = built.model.face_origins(face);
        if origins.iter().any(|o| matches!(o.role, Role::Rim { .. })) {
            let caps: BTreeSet<_> = origins
                .iter()
                .filter_map(|o| {
                    if let Role::Rim { edge } = &o.role {
                        edge.between().iter().find(|o| o.role == Role::End).cloned()
                    } else {
                        None
                    }
                })
                .collect();
            assert_eq!(
                caps.len(),
                2,
                "every rim retains both selected-face source aliases"
            );
        }
    }
    agree(
        &recipe,
        &[
            (["1", "1", "7/2"], Location::Inside),
            (["3", "3", "7/2"], Location::Outside),
        ],
    );
}

#[test]
fn non_axis_aligned_retained_face_uses_its_exact_analytic_terms() {
    let recipe = parse(include_str!("fixtures/boolean-face-thicken/l-tilted.json"));
    exact_volume(&build(&recipe), "12");
    agree(
        &recipe,
        &[
            (["6", "-21/5", "39/10"], Location::Inside),
            (["8", "-3", "11/2"], Location::Outside),
            (["6", "-17/5", "33/10"], Location::Outside),
        ],
    );
}

#[test]
fn physical_adjacency_disambiguates_a_split_origin_without_a_near_point() {
    let recipe = edited(SPLIT, |v| {
        let pick = v["steps"][3]["faces"][0].as_object_mut().unwrap();
        pick.remove("near");
        pick.insert(
            "touching".to_owned(),
            json!([{"step":"stock","role":{"side":{"segment":"x_min"}}}]),
        );
    });
    exact_volume(&build(&recipe), "8");
    agree(
        &recipe,
        &[
            (["1", "1", "7/2"], Location::Inside),
            (["5", "1", "7/2"], Location::Outside),
        ],
    );
}

#[test]
fn authored_plain_cap_keeps_its_existing_collinear_rim_alias_choice() {
    let points = [[2, 0], [4, 0], [4, 2], [2, 2], [2, 4], [0, 4], [0, 0]];
    let loop_:Vec<_>=points.iter().enumerate().map(|(i,p)|json!({"kind":"line","id":format!("edge{i}"),"to":[p[0].to_string(),p[1].to_string()]})).collect();
    let recipe=parse(&json!({"format":"catallax-recipe/4","name":"retained plain cap aliases","steps":[
        {"id":"sketch","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[loop_]},
        {"id":"stock","op":"extrude","profile":"sketch","distance":"3"},
        {"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":"end"}}],"thickness":"1","side":"front"}
    ],"result":"result"}).to_string());
    exact_volume(&build(&recipe), "12");
    agree(
        &recipe,
        &[
            (["1", "1", "7/2"], Location::Inside),
            (["3", "3", "7/2"], Location::Outside),
        ],
    );
}

#[test]
fn changing_neighbours_on_one_straight_rim_do_not_partition_the_sweep() {
    let recipe = parse(include_str!(
        "fixtures/boolean-face-thicken/collinear-neighbours.json"
    ));
    let built = build(&recipe);
    exact_volume(&built, "16");
    assert_eq!(built.model.faces_of_body(built.result).len(), 6);
    assert_eq!(built.model.edges_of_body(built.result).len(), 12);
    assert_eq!(built.model.vertices_of_body(built.result).len(), 8);
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .any(|f| {
                let origins = built.model.face_origins(f);
                origins.len() > 1 && origins.iter().all(|o| matches!(o.role, Role::Rim { .. }))
            }),
        "one merged wall retains every adjacent source origin"
    );
    agree(
        &recipe,
        &[
            (["1", "1", "7/2"], Location::Inside),
            (["7/2", "3", "7/2"], Location::Inside),
            (["1", "1", "5/2"], Location::Outside),
            (["4", "1", "7/2"], Location::Boundary),
            (["4", "3", "7/2"], Location::Boundary),
            (["4", "2", "7/2"], Location::Boundary),
        ],
    );
    // The authored triangular source is inactive in the final box. Its retained
    // history must not prevent ordinary field mass at the original defaults.
    let mass = catallax_field::mass::properties(
        &Plan::compile(&recipe).unwrap(),
        &Default::default(),
        &Default::default(),
        &mut Decider::default(),
    );
    let mass = mass.unwrap();
    assert!(mass.volume.is_exact());
    assert_eq!(mass.volume.lo, Rational::from_i64(16));
}
