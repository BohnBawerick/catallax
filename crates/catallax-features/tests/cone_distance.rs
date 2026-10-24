//! Finite conical domains, normal-image rims and independent material signs.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{Built, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, Surface};
use catallax_topology::BodyKind;
use serde_json::{Value, json};

fn pose() -> Value {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn cone(top: &str) -> Value {
    json!({"id":"stock","op":"cone","base_radius":"4","top_radius":top,"height":"4","pose":pose()})
}
fn recipe(stock: Value, operation: Value) -> Recipe {
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"conical distance","steps":[stock,operation],"result":"result"}).to_string()).unwrap()
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
fn locations(recipe: &Recipe, points: &[([&str; 3], Location)]) {
    let mut built = build(recipe);
    let field = Plan::compile(recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in points {
        let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
        let point = RPoint3::new(x, y, z);
        let got = field
            .proven_checked(&point, &mut workspace, &mut Decider::default())
            .unwrap();
        assert_eq!(got.location, *want, "field at {p:?}");
        if *want != Location::Boundary {
            let got = classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert_eq!(
                got.proven(),
                Some(&if *want == Location::Inside {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "brep at {p:?}: {got:?}"
            );
        }
    }
}
#[test]
fn inward_frustum_meets_both_shifted_caps_at_true_distance() {
    let recipe = recipe(
        cone("1"),
        json!({"id":"result","op":"offset","body":"stock","distance":"-1/5"}),
    );
    locations(
        &recipe,
        &[
            (["2", "0", "2"], Location::Inside),
            (["9/4", "0", "2"], Location::Boundary),
            (["23/10", "0", "2"], Location::Outside),
            (["0", "0", "1/5"], Location::Boundary),
            (["0", "0", "1/10"], Location::Outside),
            (["0", "0", "39/10"], Location::Outside),
            (["0", "0", "-10"], Location::Outside),
        ],
    );
}
#[test]
fn cone_tip_moves_without_inventing_a_planar_cap_and_collapse_is_empty() {
    let pointed = recipe(
        cone("0"),
        json!({"id":"result","op":"offset","body":"stock","distance":"-1/5"}),
    );
    let built = build(&pointed);
    assert_eq!(built.model.vertices_of_body(built.result).len(), 1);
    assert_eq!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .filter(|f| matches!(
                built.model.geom().surface(built.model.face(**f).surface),
                Surface::Plane(_)
            ))
            .count(),
        1
    );
    locations(
        &pointed,
        &[
            (["0", "0", "3"], Location::Inside),
            (["0", "0", "4"], Location::Outside),
        ],
    );
    let consumed = recipe(
        cone("1"),
        json!({"id":"result","op":"offset","body":"stock","distance":"-2"}),
    );
    let built = build(&consumed);
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
    locations(&consumed, &[(["0", "0", "2"], Location::Outside)]);
}
#[test]
fn conical_wall_sweeps_keep_sloping_normal_rims_on_both_sides() {
    for (side, r, z, out) in [
        ("front", "129/50", "103/50", "14/5"),
        ("back", "121/50", "97/50", "11/5"),
    ] {
        let recipe = recipe(
            cone("1"),
            json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":{"side":{"segment":"wall"}}}}],"thickness":"1/5","side":side}),
        );
        locations(
            &recipe,
            &[
                ([r, "0", z], Location::Inside),
                ([out, "0", "2"], Location::Outside),
                (["0", "0", "2"], Location::Outside),
                (["4", "0", "-1/10"], Location::Outside),
            ],
        );
        let built = build(&recipe);
        assert_eq!(built.model.faces_of_body(built.result).len(), 4);
        assert!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .all(|f| matches!(
                    built.model.geom().surface(built.model.face(*f).surface),
                    Surface::Cone(_)
                ))
        );
    }
}
#[test]
fn open_wide_cap_shell_has_a_real_rim_and_no_opening_lip() {
    let recipe = recipe(
        cone("1"),
        json!({"id":"result","op":"shell","body":"stock","open":[{"face":{"step":"stock","role":"start"}}],"thickness":"1/5"}),
    );
    locations(
        &recipe,
        &[
            (["0", "0", "1"], Location::Outside),
            (["12/5", "0", "2"], Location::Inside),
            (["0", "0", "39/10"], Location::Inside),
            (["0", "0", "0"], Location::Outside),
        ],
    );
}
#[test]
fn unsupported_tip_sweeps_keep_typed_refusals() {
    for operation in [
        json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":{"side":{"segment":"wall"}}}}],"thickness":"1/5","side":"front"}),
        json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":{"side":{"segment":"wall"}}}}],"thickness":"1/5","side":"back"}),
    ] {
        let recipe = recipe(cone("0"), operation);
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected a geometric refusal")
        };
        assert_eq!(refusal.code, "distance.cone-domain");
        assert_eq!(
            Plan::compile(&recipe).unwrap_err().code,
            "distance.cone-domain"
        );
    }
}

#[test]
fn coaxial_cone_cylinder_offsets_share_the_analytic_transition_circle() {
    let recipe = Recipe::parse(include_str!("fixtures/cone-distance/coaxial-inward.json")).unwrap();
    locations(
        &recipe,
        &[
            (["1", "0", "1"], Location::Inside),
            (["3/2", "0", "1"], Location::Boundary),
            (["23/10", "0", "3"], Location::Boundary),
            (["12/5", "0", "3"], Location::Outside),
            (["0", "0", "1/10"], Location::Outside),
        ],
    );
    let built = build(&recipe);
    assert_eq!(built.model.faces_of_body(built.result).len(), 4);
    assert_eq!(built.model.edges_of_body(built.result).len(), 3);
}

#[test]
fn oblique_cone_plane_offsets_keep_a_native_ellipse_and_finite_nappe() {
    let recipe = Recipe::parse(include_str!("fixtures/cone-distance/oblique-inward.json")).unwrap();
    locations(
        &recipe,
        &[
            (["0", "0", "5/2"], Location::Inside),
            (["1", "0", "3"], Location::Inside),
            (["2", "0", "3"], Location::Outside),
            (["0", "0", "2"], Location::Outside),
            (["0", "0", "-3"], Location::Outside),
        ],
    );
    let built = build(&recipe);
    assert!(
        built
            .model
            .edges_of_body(built.result)
            .iter()
            .any(|e| matches!(
                built.model.geom().curve(built.model.edge(*e).curve),
                catallax_shapes::Curve::Ellipse(_)
            ))
    );
}

#[test]
fn outward_frustum_rounds_both_circular_edges_at_true_distance() {
    let recipe =
        Recipe::parse(include_str!("fixtures/cone-distance/frustum-outward.json")).unwrap();
    locations(
        &recipe,
        &[
            (["27/10", "0", "2"], Location::Inside),
            (["11/4", "0", "2"], Location::Boundary),
            (["14/5", "0", "2"], Location::Outside),
            (["103/25", "0", "-4/25"], Location::Boundary),
            (["411/100", "0", "-4/25"], Location::Inside),
            (["413/100", "0", "-4/25"], Location::Outside),
            (["0", "0", "-1/5"], Location::Boundary),
            (["0", "0", "21/5"], Location::Boundary),
        ],
    );
    let built = build(&recipe);
    assert_eq!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .filter(|f| matches!(
                built.model.geom().surface(built.model.face(**f).surface),
                Surface::Torus(_)
            ))
            .count(),
        2
    );
}

#[test]
fn outward_cone_tips_keep_a_typed_coverage_refusal() {
    let recipe = Recipe::parse(include_str!("fixtures/cone-distance/tip-outward.json")).unwrap();
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected a geometric coverage refusal");
    };
    assert_eq!(refusal.code, "distance.cone-domain");
    assert_eq!(
        Plan::compile(&recipe).unwrap_err().code,
        "distance.cone-domain"
    );
}

#[test]
fn conical_cap_thickening_sweeps_only_the_selected_disk() {
    for (side, z) in [("front", "41/10"), ("back", "39/10")] {
        let recipe = recipe(
            cone("1"),
            json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock","role":"end"}}],"thickness":"1/5","side":side}),
        );
        locations(
            &recipe,
            &[
                (["9/10", "0", z], Location::Inside),
                (["101/100", "0", z], Location::Outside),
                (["0", "0", "7/2"], Location::Outside),
            ],
        );
    }
}

#[test]
fn inward_offsets_remove_swallowed_small_caps_without_an_invented_face() {
    for distance in ["-1/2", "-3/5"] {
        let stock = json!({"id":"stock","op":"cone","base_radius":"13/4","top_radius":"1/4","height":"4","pose":pose()});
        let recipe = recipe(
            stock,
            json!({"id":"result","op":"offset","body":"stock","distance":distance}),
        );
        locations(
            &recipe,
            &[
                (["1/10", "0", "3"], Location::Inside),
                (["0", "0", "18/5"], Location::Outside),
                (["0", "0", "9/2"], Location::Outside),
            ],
        );
        let built = build(&recipe);
        assert_eq!(built.model.faces_of_body(built.result).len(), 2);
        assert_eq!(built.model.vertices_of_body(built.result).len(), 1);
        assert!(built.model.faces_of_body(built.result).iter().all(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .all(|o| o.role != catallax_recipe::Role::End)
        }));
    }
}

#[test]
fn a_cap_on_the_coaxial_transition_drops_the_zero_height_cylinder_band() {
    let mut value: Value =
        serde_json::from_str(include_str!("fixtures/cone-distance/coaxial-inward.json")).unwrap();
    value["steps"].as_array_mut().unwrap().last_mut().unwrap()["distance"] = json!("-3/2");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    locations(
        &recipe,
        &[
            (["0", "0", "2"], Location::Inside),
            (["1", "0", "5/2"], Location::Boundary),
            (["0", "0", "13/5"], Location::Outside),
        ],
    );
    let built = build(&recipe);
    assert_eq!(built.model.faces_of_body(built.result).len(), 3);
    assert_eq!(built.model.edges_of_body(built.result).len(), 2);
}
