//! Bend and twist (roadmap 3.4, RECIPE.md section 4.20) through the build the
//! application calls, `evaluate_incremental_outcome`: build, save, reopen,
//! edit and STEP.
//!
//! Set `CATALLAX_DEFORM_PROOF` to a directory to keep the STEP files for the
//! raw OpenCASCADE check in `tools/step-regressions/test_bend_twist.py`.
use catallax_features::{BuildError, EvaluationOutcome, evaluate_incremental_outcome};
use catallax_numbers::{Decider, OperationWork, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_topology::numerical::NumericalBody;
use serde_json::{Value, json};

/// A 10 by 10 bar, 100 long on the z axis, and a twist about that axis.
fn twisted_bar(angle: &str) -> Recipe {
    recipe(
        json!({"id": "bar", "op": "block", "size": ["10", "10", "100"],
            "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["-5", "-5", "0"]}}),
        json!({"id": "twisted", "op": "twist", "body": "bar",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "0", "1"]},
            "angle": angle, "length": "100"}),
    )
}

/// A plate 30 long and 2 thick whose mid-plane lies 20 from the z axis,
/// bent round that axis with the mid-plane as the neutral plane.
fn bent_plate(radius: &str, toward: [&str; 3]) -> Recipe {
    recipe(
        json!({"id": "plate", "op": "block", "size": ["30", "2", "20"],
            "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["-15", "19", "0"]}}),
        json!({"id": "bent", "op": "bend", "body": "plate",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "0", "1"]},
            "toward": toward, "radius": radius}),
    )
}

fn recipe(stock: Value, warp: Value) -> Recipe {
    let result = warp["id"].clone();
    Recipe::parse(
        &json!({"format": "catallax-recipe/5", "name": "bend and twist",
            "result": result, "steps": [stock, warp]})
        .to_string(),
    )
    .unwrap()
}

/// The application's 100M allowance, as the ordinary no-display build gives.
fn build(recipe: &Recipe) -> (Result<EvaluationOutcome, BuildError>, OperationWork) {
    with_operation_budget(100_000_000, || {
        evaluate_incremental_outcome(recipe, None, &Default::default(), &mut Decider::default())
    })
}

/// The admitted body, its STEP bytes, and the work of the build.
fn built(recipe: &Recipe) -> (NumericalBody, Vec<u8>, OperationWork) {
    let (outcome, work) = build(recipe);
    assert!(!work.exhausted, "{work:?}");
    let EvaluationOutcome::Numerical(result) = outcome.unwrap() else {
        panic!("bend and twist build only as numerical geometry");
    };
    // The faces keep the stock's names: a warp makes no faces of its own.
    for origins in result.face_origins().values() {
        assert!(
            origins
                .iter()
                .all(|o| o.step.as_str() != recipe.result.as_str()),
            "{origins:?}"
        );
        assert!(!origins.is_empty());
    }
    let body = result.body().clone();
    let (step, export) = with_operation_budget(100_000_000, || {
        catallax_step::numerical::build(&body, &Default::default(), &mut Decider::default())
    });
    assert!(!export.exhausted, "{export:?}");
    let step = step.unwrap();
    step.check_written(step.bytes(), &mut Decider::default())
        .unwrap();
    println!(
        "{}: build {work:?}, export {export:?}, largest tolerance {} mm, cap {} mm",
        recipe.result,
        body.tolerances().largest_mm().to_f64_nearest(),
        body.tolerances().cap_mm.to_f64_nearest(),
    );
    (body, step.bytes().to_vec(), work)
}

fn keep(name: &str, bytes: &[u8]) {
    if let Some(directory) = std::env::var_os("CATALLAX_DEFORM_PROOF") {
        std::fs::write(std::path::Path::new(&directory).join(name), bytes).unwrap();
    }
}

fn positions(body: &NumericalBody) -> Vec<[f64; 3]> {
    body.definition()
        .vertices
        .iter()
        .map(|v| v.position.each_ref().map(Rational::to_f64_nearest))
        .collect()
}

/// Every expected corner is a vertex of the body, within `within` mm.
fn corners_match(body: &NumericalBody, expected: &[[f64; 3]], within: f64) {
    let actual = positions(body);
    assert_eq!(actual.len(), expected.len());
    for e in expected {
        let nearest = actual
            .iter()
            .map(|a| (0..3).map(|i| (a[i] - e[i]).powi(2)).sum::<f64>().sqrt())
            .fold(f64::INFINITY, f64::min);
        assert!(nearest <= within, "corner {e:?} is {nearest:e} mm away");
    }
}

fn block_corners(low: [f64; 3], size: [f64; 3]) -> Vec<[f64; 3]> {
    (0..8)
        .map(|i| std::array::from_fn(|axis| low[axis] + size[axis] * ((i >> axis) & 1) as f64))
        .collect()
}

fn twist(p: [f64; 3], degrees: f64) -> [f64; 3] {
    let (sin, cos) = (degrees.to_radians() * p[2] / 100.0).sin_cos();
    [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1], p[2]]
}

/// The bend of `bent_plate` with `toward` +y: binormal -x, angle -x / 20.
fn bend(p: [f64; 3]) -> [f64; 3] {
    let (u, v) = (p[1], -p[0]);
    let (sin, cos) = (v / 20.0).sin_cos();
    [-u * sin, u * cos, p[2]]
}

#[test]
fn a_twisted_bar_builds_saves_reopens_edits_and_exports() {
    let recipe = twisted_bar("90");
    let (body, step, work) = built(&recipe);
    assert_eq!(body.definition().faces.len(), 6);
    assert_eq!(body.definition().edges.len(), 12);
    assert!(work.used <= 20_000_000, "{work:?}");
    let cap = body.tolerances().cap_mm.to_f64_nearest();
    let bar = block_corners([-5.0, -5.0, 0.0], [10.0, 10.0, 100.0]);
    corners_match(
        &body,
        &bar.iter().map(|p| twist(*p, 90.0)).collect::<Vec<_>>(),
        cap,
    );
    keep("twisted.step", &step);

    // Saving writes the recipe; reopening rebuilds it to the same file.
    let reopened = Recipe::parse(&recipe.to_text()).unwrap();
    assert_eq!(reopened, recipe);
    let (_, again, rework) = built(&reopened);
    assert_eq!(
        again, step,
        "a reopened twist must write the same STEP bytes"
    );
    assert_eq!(rework.used, work.used);

    // A parameter edit rebuilds the whole twist.
    let (edited, edited_step, _) = built(&twisted_bar("120"));
    corners_match(
        &edited,
        &bar.iter().map(|p| twist(*p, 120.0)).collect::<Vec<_>>(),
        cap,
    );
    assert_ne!(edited_step, step);
    keep("twisted-edited.step", &edited_step);
}

#[test]
fn a_bent_plate_builds_saves_reopens_and_exports() {
    let recipe = bent_plate("20", ["0", "1", "0"]);
    let (body, step, work) = built(&recipe);
    assert_eq!(body.definition().faces.len(), 6);
    assert!(work.used <= 20_000_000, "{work:?}");
    let cap = body.tolerances().cap_mm.to_f64_nearest();
    let plate = block_corners([-15.0, 19.0, 0.0], [30.0, 2.0, 20.0]);
    corners_match(
        &body,
        &plate.iter().map(|p| bend(*p)).collect::<Vec<_>>(),
        cap,
    );
    keep("bent.step", &step);
    let reopened = Recipe::parse(&recipe.to_text()).unwrap();
    let (_, again, _) = built(&reopened);
    assert_eq!(
        again, step,
        "a reopened bend must write the same STEP bytes"
    );
}

#[test]
fn bends_outside_smlibs_domain_and_curved_input_refuse_by_name() {
    let refusal = |recipe: &Recipe| {
        let (outcome, work) = build(recipe);
        assert!(!work.exhausted, "{work:?}");
        let BuildError::Geometry { refusal, .. } = outcome.unwrap_err() else {
            panic!("expected a geometric refusal");
        };
        refusal.code
    };
    // The plate lies on the far side of the axis from `toward`.
    assert_eq!(
        refusal(&bent_plate("20", ["0", "-1", "0"])),
        "deform.bend-axis"
    );
    // `toward` parallel to an oblique axis has no part across it. In binary64
    // the cross-axis part of [3, 5, 7] against itself is a rounding residue,
    // so the decision is exact.
    let mut oblique: Value = serde_json::from_str(&bent_plate("20", ["3", "5", "7"]).to_text())
        .expect("canonical recipe text");
    oblique["steps"][1]["axis"]["direction"] = json!(["3", "5", "7"]);
    let oblique = Recipe::parse(&oblique.to_string()).expect("a valid recipe");
    assert_eq!(refusal(&oblique), "deform.definition");
    // Thirty millimetres round a 4 mm radius is more than one turn.
    assert_eq!(
        refusal(&bent_plate("4", ["0", "1", "0"])),
        "deform.bend-overlap"
    );
    // Curved faces are coverage follow-on item 3, not this first build.
    let cylinder = recipe(
        json!({"id": "rod", "op": "cylinder", "radius": "5", "height": "100",
            "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}}),
        json!({"id": "twisted", "op": "twist", "body": "rod",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "0", "1"]},
            "angle": "90", "length": "100"}),
    );
    assert!(refusal(&cylinder).starts_with("deform.unsupported-"));
    // The exact-only build names the numerical path instead of guessing.
    let (exact, _) = with_operation_budget(100_000_000, || {
        catallax_features::evaluate_incremental(&twisted_bar("90"), None, &mut Decider::default())
    });
    let BuildError::Geometry { refusal, .. } = exact.unwrap_err() else {
        panic!("expected a geometric refusal");
    };
    assert_eq!(refusal.code, "deform.requires-numerical-evaluation");
}

#[test]
fn a_toward_nearly_along_the_axis_keeps_its_exact_part_across() {
    // The part of `toward` across the z axis is 1e-400 along +y. It underflows
    // in binary64, so only the exact inputs keep it, and the bend is the one
    // `toward` +y gives.
    let (_, near, _) = built(&bent_plate("20", ["0", "1e-400", "1"]));
    let (_, plain, _) = built(&bent_plate("20", ["0", "1", "0"]));
    assert_eq!(near, plain, "the same frame must write the same STEP bytes");
}

#[test]
fn bend_and_twist_need_recipe_format_five() {
    let mut text: Value = serde_json::from_str(&twisted_bar("90").to_text()).unwrap();
    text["format"] = json!("catallax-recipe/4");
    let error = Recipe::parse(&text.to_string()).unwrap_err();
    assert!(error.to_string().contains("catallax-recipe/5"), "{error}");
}
