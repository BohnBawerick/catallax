//! Consumer-visible distance semantics, regularisation, names and independent field signs.

use catallax_features::{Built, evaluate};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use catallax_topology::BodyKind;
use serde_json::{Value, json};

fn pose() -> Value {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn recipe(stock: Value, operation: Value) -> Recipe {
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"distance feature","steps":[stock,operation],"result":"result"}).to_string()).unwrap()
}
fn block() -> Value {
    json!({"id":"stock","op":"block","size":["10","8","6"],"pose":pose()})
}
fn sphere() -> Value {
    json!({"id":"stock","op":"sphere","radius":"2","pose":pose()})
}
fn face(role: Value) -> Value {
    json!({"face":{"step":"stock","role":role}})
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
fn volume(built: &Built, expected: i64) {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(mass.volume.lo, Rational::from_i64(expected));
    assert_eq!(mass.volume.hi, Rational::from_i64(expected));
}
fn locations(recipe: &Recipe, points: &[([&str; 3], Location)]) {
    let plan = Plan::compile(recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, want) in points {
        let [x, y, z] = point.map(|p| Rational::parse(p).unwrap());
        let actual = plan
            .proven_checked(
                &RPoint3::new(x, y, z),
                &mut workspace,
                &mut Decider::default(),
            )
            .unwrap();
        assert_eq!(actual.location, *want, "at {point:?}: {actual:?}");
    }
}

#[test]
fn open_tray_keeps_a_rim_and_true_inward_walls() {
    let recipe = recipe(
        block(),
        json!({"id":"result","op":"shell","body":"stock","open":[face(json!("end"))],"thickness":"1"}),
    );
    let built = build(&recipe);
    volume(&built, 240);
    let origins: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f))
        .collect();
    assert!(
        origins
            .iter()
            .any(|o| o.step.as_str() == "stock" && o.role == Role::End),
        "the open face's rim retains its origin"
    );
    assert!(
        origins
            .iter()
            .any(|o| matches!(&o.role, Role::Inner {of} if of.role == Role::Start)),
        "the cavity floor names its kept generator"
    );
    locations(
        &recipe,
        &[
            (["5", "4", "1/2"], Location::Inside),
            (["5", "4", "1"], Location::Boundary),
            (["5", "4", "5"], Location::Outside),
            (["1/2", "4", "6"], Location::Boundary),
            (["5", "4", "6"], Location::Outside),
        ],
    );
}

#[test]
fn consuming_a_cavity_preserves_the_solid_and_opening_every_face_removes_it() {
    let consumed = recipe(
        block(),
        json!({"id":"result","op":"shell","body":"stock","open":[],"thickness":"3"}),
    );
    volume(&build(&consumed), 480);
    locations(&consumed, &[(["5", "4", "3"], Location::Inside)]);
    let empty = recipe(
        block(),
        json!({"id":"result","op":"shell","body":"stock","open":[{"face":{"step":"stock"}}],"thickness":"1"}),
    );
    let built = build(&empty);
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
    locations(&empty, &[(["1", "1", "1"], Location::Outside)]);
}

#[test]
fn inward_plane_offsets_regularise_collapsed_cells() {
    let shrunk = recipe(
        block(),
        json!({"id":"result","op":"offset","body":"stock","distance":"-1"}),
    );
    volume(&build(&shrunk), 192);
    locations(
        &shrunk,
        &[
            (["5", "4", "3"], Location::Inside),
            (["1", "4", "3"], Location::Boundary),
            (["0", "4", "3"], Location::Outside),
        ],
    );
    let vanished = recipe(
        block(),
        json!({"id":"result","op":"offset","body":"stock","distance":"-3"}),
    );
    let built = build(&vanished);
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
    locations(&vanished, &[(["5", "4", "3"], Location::Outside)]);
}

#[test]
fn planar_thicken_has_normal_sweep_rims_without_a_distance_lip() {
    for (side, z, outside) in [("front", "13/2", "11/2"), ("back", "11/2", "13/2")] {
        let recipe = recipe(
            block(),
            json!({"id":"result","op":"thicken","body":"stock","faces":[face(json!("end"))],"thickness":"1","side":side}),
        );
        let built = build(&recipe);
        volume(&built, 80);
        locations(
            &recipe,
            &[
                (["5", "4", z], Location::Inside),
                (["5", "4", outside], Location::Outside),
                (["-1/10", "4", z], Location::Outside),
            ],
        );
        assert!(built.model.faces_of_body(built.result).iter().any(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| matches!(o.role, Role::Rim { .. }))
        }));
    }
}

#[test]
fn radial_back_sweep_regularises_folds_instead_of_rejecting_or_clipping_them() {
    for (t, radius, centre) in [
        ("1", "3/2", Location::Outside),
        ("3", "1", Location::Inside),
        ("5", "5/2", Location::Inside),
    ] {
        let recipe = recipe(
            sphere(),
            json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock"}}],"thickness":t,"side":"back"}),
        );
        let built = build(&recipe);
        locations(
            &recipe,
            &[
                ([radius, "0", "0"], Location::Inside),
                (["0", "0", "0"], centre),
                (["4", "0", "0"], Location::Outside),
            ],
        );
        let far = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .any(|o| matches!(o.role, Role::Offset { .. }));
        assert_eq!(far, t != "3", "buried far sheets have no boundary origin");
    }
}

#[test]
fn target_tray_boolean_preserves_the_open_cavity() {
    let recipe = Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap();
    let built = build(&recipe);
    assert!(built.model.faces_of_body(built.result).iter().any(|f| {
        built
            .model
            .face_origins(*f)
            .iter()
            .any(|o| o.step.as_str() == "rounding")
    }));
    locations(
        &recipe,
        &[
            (["0", "0", "1"], Location::Inside),
            (["0", "0", "10"], Location::Outside),
            (["19", "0", "10"], Location::Inside),
            (["19", "14", "1"], Location::Outside),
            (["0", "0", "12"], Location::Outside),
        ],
    );
}

#[test]
fn outward_polyhedral_offset_rounds_edges_and_vertices_at_true_distance() {
    let recipe = recipe(
        block(),
        json!({"id":"result","op":"offset","body":"stock","distance":"1"}),
    );
    let built = build(&recipe);
    let mut planes = 0;
    let mut cylinders = 0;
    let mut spheres = 0;
    for f in built.model.faces_of_body(built.result) {
        match built.model.geom().surface(built.model.face(f).surface) {
            catallax_shapes::Surface::Plane(_) => planes += 1,
            catallax_shapes::Surface::Cylinder(_) => cylinders += 1,
            catallax_shapes::Surface::Sphere(_) => spheres += 1,
            other => panic!("unexpected distance support {other:?}"),
        }
    }
    assert_eq!((planes, cylinders, spheres), (6, 12, 8));
    locations(
        &recipe,
        &[
            (["-1/2", "-1/2", "-1/2"], Location::Inside),
            (["-3/4", "-3/4", "-3/4"], Location::Outside),
            (["-3/5", "-4/5", "3"], Location::Boundary),
            (["-1", "4", "3"], Location::Boundary),
        ],
    );
    let measured = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(measured.marks.is_empty(), "{:?}", measured.marks);
    let expected = catallax_numbers::Ball::from_rational(&Rational::from_i64(856), 512)
        + catallax_numbers::Ball::pi(512)
            * catallax_numbers::Ball::from_rational(&Rational::from_ratio(76, 3).unwrap(), 512);
    assert!(measured.value.volume.lo <= expected.lo());
    assert!(measured.value.volume.hi >= expected.hi());
}
