//! Independent field trim and inner-wall material controls.
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn source() -> Value {
    serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/draft-annular/both.json"
    ))
    .unwrap()
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(
        Rational::from_i64(x),
        Rational::from_i64(y),
        Rational::from_i64(z),
    )
}
fn origin(role: Value) -> Origin {
    serde_json::from_value(json!({"step":"base","role":role})).unwrap()
}

#[test]
fn inner_material_normal_and_finite_annular_cap_queries_are_exact() {
    let plan = Plan::compile(&Recipe::parse(&source().to_string()).unwrap()).unwrap();
    for (segment, point, radial_sign) in [("inner", p(3, 0, 1), -1.0), ("outer", p(7, 0, 1), 1.0)] {
        let name = origin(json!({"side":{"segment":segment}}));
        let id = plan.select(&name).unwrap();
        let query = plan
            .face_query(id, &point, &mut Workspace::default())
            .unwrap();
        assert_eq!(query.state, FaceState::OnFace);
        let Direction::Unique(normal) = query.direction else {
            panic!("{query:?}")
        };
        assert!(normal[0].lo() * radial_sign > 0.0 && normal[0].hi() * radial_sign > 0.0);
        assert!(normal[2].lo() > 0.0);
    }
    for (role, z, inner) in [(json!("start"), 0, 2), (json!("end"), 2, 4)] {
        let id = plan.select(&origin(role)).unwrap();
        for (radius, expected) in [
            (0, FaceState::NotOnFace),
            (inner + 1, FaceState::OnFace),
            (9, FaceState::NotOnFace),
        ] {
            assert_eq!(
                plan.face_query(id, &p(radius, 0, z), &mut Workspace::default())
                    .unwrap()
                    .state,
                expected
            );
        }
    }
    let id = plan
        .select(&origin(json!({"side":{"segment":"inner"}})))
        .unwrap();
    assert_eq!(
        plan.face_query(id, &p(1, 0, -1), &mut Workspace::default())
            .unwrap()
            .state,
        FaceState::NotOnFace
    );
}

#[test]
fn independent_carriers_retain_hole_after_reflection_and_scalar_rebuild() {
    for reflect in [false, true] {
        let mut source = source();
        source["steps"][0]["loops"][1][0]["radius"] = json!("3");
        if reflect {
            source["steps"].as_array_mut().unwrap().push(json!({"id":"copy","op":"mirror","body":"drafted","keep_original":false,"plane":{"point":["0","0","0"],"normal":["1","0","0"]}}));
            source["result"] = json!("copy");
        }
        let plan = Plan::compile(&Recipe::parse(&source.to_string()).unwrap()).unwrap();
        for (radius, expected) in [
            (0, Location::Outside),
            (3, Location::Outside),
            (4, Location::Boundary),
            (5, Location::Inside),
            (7, Location::Boundary),
            (8, Location::Outside),
        ] {
            let point = p(if reflect { -radius } else { radius }, 0, 1);
            assert_eq!(
                plan.proven(&point, &mut Workspace::default()).location,
                expected
            );
        }
    }
}

#[test]
fn nearest_selection_distinguishes_annular_walls_and_ties_stay_ambiguous() {
    let mut source = source();
    // The hole axis side of the inner wall: the inner wall alone is drafted.
    source["steps"][2]["faces"] = json!([{"face":{"step":"base"},"near":["1","0","1"]}]);
    let plan = Plan::compile(&Recipe::parse(&source.to_string()).unwrap()).unwrap();
    let mut storage = Workspace::default();
    // Inner radius 2+z at height 1; the outer wall stays at radius 8.
    assert_eq!(
        plan.proven(&p(3, 0, 1), &mut storage).location,
        Location::Boundary
    );
    assert_eq!(
        plan.proven(&p(4, 0, 1), &mut storage).location,
        Location::Inside
    );
    assert_eq!(
        plan.proven(&p(7, 0, 1), &mut storage).location,
        Location::Inside
    );
    assert_eq!(
        plan.proven(&p(2, 0, 1), &mut storage).location,
        Location::Outside
    );
    source["steps"][2]["faces"] = json!([{"face":{"step":"base"},"near":["5","0","1"]}]);
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    assert_eq!(
        Plan::compile(&recipe).unwrap_err().code,
        "recipe.ambiguous-reference"
    );
}

#[test]
fn retained_nearest_selection_accepts_aliases_of_one_physical_cylinder() {
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    let document = json!({"format":"catallax-recipe/3","name":"retained cylinder aliases","steps":[
        {"id":"left","op":"cylinder","radius":"8","height":"2","pose":pose},
        {"id":"right","op":"cylinder","radius":"8","height":"2","pose":pose},
        {"id":"source","op":"boolean","kind":"union","target":"left","tool":"right"},
        {"id":"drafted","op":"draft","body":"source","faces":[{"face":{"step":"left"},"near":["8","0","1"]}],"pull":["0","0","1"],"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"45"}
    ],"result":"drafted"});
    let recipe = Recipe::parse(&document.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let plan = Plan::compile(&recipe).unwrap();
    assert_eq!(
        plan.proven(&p(0, 0, 1), &mut Workspace::default()).location,
        Location::Inside
    );
    for creator in ["left", "right"] {
        let name: Origin =
            serde_json::from_value(json!({"step":creator,"role":{"side":{"segment":"wall"}}}))
                .unwrap();
        let id = plan.select(&name).unwrap();
        assert_eq!(
            plan.face_query(id, &p(7, 0, 1), &mut Workspace::default())
                .unwrap()
                .state,
            FaceState::OnFace
        );
    }
}

#[test]
fn true_distance_reaches_finite_annuli_and_both_wall_segments() {
    let plan = Plan::compile(&Recipe::parse(&source().to_string()).unwrap()).unwrap();
    for (point, oracle) in [
        (p(0, 0, 0), 2.0),
        (p(0, 0, -1), 5.0_f64.sqrt()),
        (p(0, 0, 3), 12.5_f64.sqrt()),
        (p(5, 0, 1), -1.0),
    ] {
        let query = plan
            .proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
            .unwrap();
        let catallax_field::Value::TrueDistance(distance) = query.value else {
            panic!("{query:?}")
        };
        assert!(
            distance.contains(oracle),
            "{point:?}: {query:?} vs {oracle}"
        );
        assert!(distance.hi() - distance.lo() < 1e-9);
        let real = point.to_real();
        let fast = plan
            .fast(
                [real.x.to_f64(), real.y.to_f64(), real.z.to_f64()],
                &mut Workspace::default(),
            )
            .unwrap();
        assert!(fast.true_distance_formula);
        assert!((fast.value - oracle).abs() < 1e-12);
    }
}
