//! Complete target copies preserve names, exact output and work across threads.
use catallax_app::pattern::{FeatureInstance, expand_feature};
use catallax_features::evaluate;
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{Isometry, RPoint3};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::num::NonZeroUsize;

fn id(value: &str) -> StepId {
    StepId::new(value).unwrap()
}

// These fixtures have only rational dimensions; rotations and angles stay unchanged.
fn scale(value: &mut Value, dimension: bool) {
    match value {
        Value::String(text) if dimension => {
            *text = (Rational::parse(text).unwrap() * Rational::from_i64(2)).to_string();
        }
        Value::Array(values) => {
            for value in values {
                scale(value, dimension);
            }
        }
        Value::Object(values) => {
            for (key, value) in values {
                scale(
                    value,
                    matches!(
                        key.as_str(),
                        "translation" | "to" | "radius" | "height" | "distance" | "size"
                    ),
                );
            }
        }
        _ => (),
    }
}

fn placement(kind: &str) -> (Value, Isometry) {
    let operation = match kind {
        "linear" => json!({"id":"copies","op":"pattern","body":"unused",
            "layout":{"linear":{"vector":["500","0","0"],"count":2}}}),
        "circular" => json!({"id":"copies","op":"pattern","body":"unused",
            "layout":{"circular":{"axis":{"point":["-500","0","0"],
                "direction":["0","0","7"]},"angle":"90","count":2}}}),
        "mirror" => json!({"id":"copies","op":"mirror","body":"unused",
            "plane":{"point":["-500","0","0"],"normal":["2","3","0"]},
            "keep_original":false}),
        _ => unreachable!(),
    };
    let step: catallax_recipe::Step = serde_json::from_value(operation.clone()).unwrap();
    let map = match step.op {
        catallax_recipe::Operation::Pattern { layout, .. } => {
            catallax_features::pattern_placement(&layout, 1, &mut Decider::default()).unwrap()
        }
        catallax_recipe::Operation::Mirror { plane, .. } => {
            catallax_features::mirror_placement(&plane).unwrap()
        }
        _ => unreachable!(),
    };
    (operation, map)
}

fn copy(source: &Recipe, kind: &str, feature: bool) -> (Recipe, Isometry) {
    let (mut operation, map) = placement(kind);
    let recipe = if feature {
        let graph = source
            .steps
            .iter()
            .map(|s| s.id.clone())
            .collect::<Vec<_>>();
        let instance = FeatureInstance {
            placement: map.clone(),
            steps: graph
                .iter()
                .map(|s| (s.clone(), id(&format!("copy-{s}"))))
                .collect(),
            inputs: BTreeMap::new(),
            origins: BTreeMap::new(),
        };
        expand_feature(
            source,
            &graph,
            &source.result,
            &[instance],
            &mut Decider::default(),
        )
        .unwrap()
    } else {
        operation["body"] = json!(source.result);
        let mut recipe = source.clone();
        recipe
            .steps
            .push(serde_json::from_value(operation).unwrap());
        recipe.result = id("copies");
        recipe
    };
    let saved = Recipe::parse(&recipe.to_text()).unwrap();
    assert_eq!(saved, recipe);
    (saved, map)
}

struct Snapshot {
    names: Vec<String>,
    step: Vec<u8>,
    field: String,
    work: u64,
}

fn snapshot(recipe: &Recipe, map: &Isometry, inside: &RPoint3, on_face: &RPoint3) -> Snapshot {
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate(recipe, &mut dec).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let (seed, copies) = match &recipe.steps.last().unwrap().op {
            catallax_recipe::Operation::Pattern { body, layout } => (body.clone(), layout.count()),
            catallax_recipe::Operation::Mirror { body, .. } => (body.clone(), 1),
            _ => (id(recipe.result.as_str().strip_prefix("copy-").unwrap()), 1),
        };
        let source_faces = built.model.faces_of_body(built.bodies[&seed]).len();
        assert_eq!(
            built.model.faces_of_body(built.result).len(),
            source_faces * copies as usize
        );
        let mut names = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .map(|face| {
                let mut origins = built.model.face_origins(face).to_vec();
                origins.sort();
                serde_json::to_string(&origins).unwrap()
            })
            .collect::<Vec<_>>();
        names.sort();
        assert!(!names.is_empty());
        assert!(
            names
                .iter()
                .any(|name| name.contains("copy-") || name.contains("\"copy\""))
        );

        let plan = Plan::compile_checked(recipe, &mut dec).unwrap();
        let inside = map.apply_rational_point(inside).unwrap();
        let on_face = map.apply_rational_point(on_face).unwrap();
        let mut workspace = Workspace::default();
        let query = plan
            .proven_checked(&inside, &mut workspace, &mut dec)
            .unwrap();
        assert_eq!(query.location, Location::Inside);
        // A registry hit alone cannot establish that a copied face survived trimming.
        let mut boundary = vec![];
        for index in 0..plan.faces().len() {
            if let Some(face) = plan.result_face(index)
                && plan
                    .face_query_checked(face, &on_face, &mut workspace, &mut dec)
                    .unwrap()
                    .state
                    == FaceState::OnFace
            {
                boundary.push(serde_json::to_string(&plan.faces()[index].origin).unwrap());
            }
        }
        boundary.sort();
        assert!(!boundary.is_empty());
        let mut step = vec![];
        let export = catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut dec,
            &mut step,
        )
        .unwrap();
        assert!(export.marks.is_empty());
        assert!(export.value.max_deviation <= Rational::parse("1e-7").unwrap());
        (names, step, format!("{query:?}:{boundary:?}"))
    });
    assert!(!work.exhausted);
    Snapshot {
        names: result.0,
        step: result.1,
        field: result.2,
        work: work.used,
    }
}

enum Target {
    Bracket,
    Collar,
}

fn check_target(target: Target, kind: &str, feature: bool, threads: usize) {
    // Cargo can run cases in one process; the kernel thread setting is process-wide.
    static THREADS: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _guard = THREADS
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let (fixture, edited_step, parameter, edit, inside, on_face) = match target {
        Target::Bracket => (
            include_str!("../../../parts/P1-bracket.json"),
            "boss",
            "radius",
            "9",
            ["36", "25", "10"],
            ["60", "20", "2"],
        ),
        Target::Collar => (
            include_str!("../../../parts/P2-collar.json"),
            "collar",
            "distance",
            "0.75",
            ["9", "6", "0"],
            ["12", "6", "0"],
        ),
    };
    let mut source: Value = serde_json::from_str(fixture).unwrap();
    // Exercise doubled dimensions as well as moved and turned placements.
    if kind == "circular" {
        scale(&mut source, false);
    }
    let factor = if kind == "circular" { 2 } else { 1 };
    let point = |coordinates: [&str; 3]| {
        let [x, y, z] =
            coordinates.map(|s| Rational::parse(s).unwrap() * Rational::from_i64(factor));
        RPoint3::new(x, y, z)
    };
    let inside = point(inside);
    let before = Recipe::parse(&source.to_string()).unwrap();
    let changed = source["steps"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|s| s["id"] == edited_step)
        .unwrap();
    changed[parameter] =
        json!((Rational::parse(edit).unwrap() * Rational::from_i64(factor)).to_string());
    let after = Recipe::parse(&source.to_string()).unwrap();
    let (before, map) = copy(&before, kind, feature);
    let (after, _) = copy(&after, kind, feature);
    // Use a face unchanged by either edit, at a point away from its edges.
    let on_face = point(on_face);
    set_threads(NonZeroUsize::new(1));
    let baseline = [
        snapshot(&before, &map, &inside, &on_face),
        snapshot(&after, &map, &inside, &on_face),
    ];
    assert_eq!(
        baseline[0].names, baseline[1].names,
        "{edited_step}/{kind}/feature={feature}"
    );
    assert!(
        baseline[0].step != baseline[1].step,
        "the edit must change geometry"
    );
    set_threads(NonZeroUsize::new(threads));
    for (recipe, expected) in [(&before, &baseline[0]), (&after, &baseline[1])] {
        let actual = snapshot(recipe, &map, &inside, &on_face);
        let context = format!("{edited_step}/{kind}/feature={feature}/threads={threads}");
        assert_eq!(actual.names, expected.names, "names: {context}");
        assert!(actual.step == expected.step, "STEP bytes: {context}");
        assert_eq!(actual.field, expected.field, "field: {context}");
        assert_eq!(actual.work, expected.work, "work: {context}");
    }
    println!(
        "{edited_step}/{kind}/feature={feature}: {} faces; work {}/{}; STEP bytes {}/{}; threads 1,{threads}",
        baseline[0].names.len(),
        baseline[0].work,
        baseline[1].work,
        baseline[0].step.len(),
        baseline[1].step.len()
    );
    set_threads(None);
}

macro_rules! matrix_case {
    ($name:ident, $target:ident, $kind:literal, $feature:literal) => {
        mod $name {
            use super::*;

            #[test]
            fn two_threads() {
                check_target(Target::$target, $kind, $feature, 2);
            }

            #[test]
            fn four_threads() {
                check_target(Target::$target, $kind, $feature, 4);
            }

            #[test]
            fn all_threads() {
                check_target(
                    Target::$target,
                    $kind,
                    $feature,
                    std::thread::available_parallelism().unwrap().get(),
                );
            }
        }
    };
}

matrix_case!(p1_linear_body, Bracket, "linear", false);
matrix_case!(p1_linear_feature, Bracket, "linear", true);
matrix_case!(p1_circular_body, Bracket, "circular", false);
matrix_case!(p1_circular_feature, Bracket, "circular", true);
matrix_case!(p1_mirror_body, Bracket, "mirror", false);
matrix_case!(p1_mirror_feature, Bracket, "mirror", true);
matrix_case!(p2_linear_body, Collar, "linear", false);
matrix_case!(p2_linear_feature, Collar, "linear", true);
matrix_case!(p2_circular_body, Collar, "circular", false);
matrix_case!(p2_circular_feature, Collar, "circular", true);
matrix_case!(p2_mirror_body, Collar, "mirror", false);
matrix_case!(p2_mirror_feature, Collar, "mirror", true);
