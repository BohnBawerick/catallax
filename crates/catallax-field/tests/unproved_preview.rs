//! The unproved processor preview: its fence and its sign agreement with
//! the deterministic evaluator.

use catallax_field::{Plan, PreviewHardware, UnprovedBatch, Workspace};
use catallax_recipe::Recipe;
use std::path::{Path, PathBuf};

const NAMES: [&str; 4] = [
    "unproved_preview",
    "UnprovedValue",
    "PreviewHardware",
    "UnprovedBatch",
];

// Library source that may name the hardware preview. Everything else in the
// workspace's library and binary sources may not, so no export, proven
// measurement or proof route can reach it. A display caller added later is
// listed here in its own reviewed change; `FORBIDDEN` paths never are.
const ALLOWED: [&str; 5] = [
    "crates/catallax-field/src/lib.rs",
    "crates/catallax-field/src/unproved_batch.rs",
    "crates/catallax-field/src/unproved_wgsl.rs",
    // GPU batch evaluator with the processor preview as fallback. It returns
    // unproved `Samples` to the preview crate's first paint and checking
    // binary only; no other crate depends on catallax-preview, and it has no
    // writer, proof or proven-measurement output.
    "crates/catallax-preview/src/batch.rs",
    "crates/catallax-preview/src/main.rs",
];
const FORBIDDEN: [&str; 9] = [
    "crates/catallax-field/src/printing",
    "crates/catallax-field/src/measurement",
    "crates/catallax-field/src/mass",
    "crates/catallax-field/src/distance",
    "crates/catallax-measure/",
    "crates/catallax-step/",
    "crates/catallax-iges/",
    "crates/catallax-compare/",
    "crates/catallax-cli/src/printing",
];

fn sources(dir: &Path, found: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            sources(&path, found);
        } else if path.extension().is_some_and(|e| e == "rs") {
            found.push(path);
        }
    }
}

#[test]
fn hardware_preview_is_fenced_from_proof_and_export_routes() {
    assert!(
        ALLOWED
            .iter()
            .all(|a| !FORBIDDEN.iter().any(|f| a.starts_with(f)))
    );
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut files = Vec::new();
    for crate_dir in std::fs::read_dir(root.join("crates")).unwrap() {
        let src = crate_dir.unwrap().path().join("src");
        if src.is_dir() {
            sources(&src, &mut files);
        }
    }
    assert!(
        files.len() > 100,
        "workspace sources found: {}",
        files.len()
    );
    let mut definition = false;
    let mut offenders = Vec::new();
    for file in files {
        let relative = file
            .strip_prefix(&root)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        let text = std::fs::read_to_string(&file).unwrap();
        if !NAMES.iter().any(|n| text.contains(n)) {
            continue;
        }
        if ALLOWED.contains(&relative.as_str()) {
            definition = true;
        } else {
            offenders.push(relative);
        }
    }
    assert!(definition, "the allowed definition file names the preview");
    assert!(
        offenders.is_empty(),
        "unreviewed callers of the hardware preview: {offenders:?}"
    );
}

fn sample(n: usize, size: [f64; 3], mut state: u64) -> Vec<[f64; 3]> {
    (0..n)
        .map(|_| {
            std::array::from_fn(|i| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
                size[i] * (1.2 * unit - 0.1)
            })
        })
        .collect()
}

#[test]
fn processor_preview_keeps_the_deterministic_sign_on_p1_and_p3() {
    for (recipe, size) in [
        (
            include_str!("../../../parts/P1-bracket.json"),
            [60.0, 40.0, 30.0],
        ),
        (
            include_str!("../../../parts/P3-filleted-boss.json"),
            [40.0, 30.0, 22.0],
        ),
    ] {
        let plan = Plan::compile(&Recipe::parse(recipe).unwrap()).unwrap();
        let mut storage = Workspace::default();
        let (mut inside, mut outside) = (0, 0);
        for point in sample(20_000, size, 7) {
            let fast = plan.fast(point, &mut storage).unwrap();
            let preview = plan.unproved_preview(point, &mut storage).unwrap();
            assert_eq!(preview.hardware, PreviewHardware::Processor);
            assert_eq!(preview.true_distance_formula, fast.true_distance_formula);
            assert_eq!(
                (preview.value < 0.0, preview.value > 0.0),
                (fast.value < 0.0, fast.value > 0.0),
                "sign differs at {point:?}: {} against {}",
                preview.value,
                fast.value
            );
            assert!((preview.value - fast.value).abs() <= 1e-9 * (1.0 + fast.value.abs()));
            if fast.value < 0.0 {
                inside += 1;
            } else {
                outside += 1;
            }
        }
        assert!(inside > 1_000 && outside > 1_000, "{inside} {outside}");
    }
}

#[test]
fn processor_preview_refuses_what_fast_refuses() {
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    let mut storage = Workspace::default();
    for point in [[f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]] {
        let fast = plan.fast(point, &mut storage).unwrap_err();
        let preview = plan.unproved_preview(point, &mut storage).unwrap_err();
        assert_eq!(preview.code, fast.code);
    }
}

#[test]
fn threaded_preview_has_scalar_bits_at_every_worker_count() {
    for (recipe, size) in [
        (
            include_str!("../../../parts/P1-bracket.json"),
            [60.0, 40.0, 30.0],
        ),
        (
            include_str!("../../../parts/P3-filleted-boss.json"),
            [40.0, 30.0, 22.0],
        ),
    ] {
        let plan = Plan::compile(&Recipe::parse(recipe).unwrap()).unwrap();
        let points = sample(20_003, size, 7);
        let mut w = Workspace::default();
        let scalar: Vec<_> = points
            .iter()
            .map(|p| plan.unproved_preview(*p, &mut w).unwrap())
            .collect();
        let mut batch = UnprovedBatch::default();
        // Reuse storage while changing worker count and input length.
        for workers in [1, 2, 4, 8, 4, 1] {
            for len in [0, 1, 3, 17, points.len()] {
                let values = plan
                    .unproved_preview_batch(&points[..len], workers, &mut batch)
                    .unwrap();
                assert_eq!(values.len(), len);
                for (i, (a, b)) in values.iter().zip(&scalar).enumerate() {
                    assert_eq!(
                        a.value.to_bits(),
                        b.value.to_bits(),
                        "worker count {workers}, point {i}"
                    );
                    assert_eq!(a.hardware, b.hardware);
                    assert_eq!(a.true_distance_formula, b.true_distance_formula);
                }
            }
        }
        for point in [
            [f64::NAN, 0.0, 0.0],
            [0.0, f64::NEG_INFINITY, 0.0],
            [f64::MAX; 3],
        ] {
            let scalar = plan.unproved_preview(point, &mut w).unwrap_err();
            for workers in [1, 2, 4, 8] {
                let error = plan
                    .unproved_preview_batch(&[point; 17], workers, &mut batch)
                    .unwrap_err();
                assert_eq!(error.code, scalar.code);
                assert_eq!(error.detail, scalar.detail);
            }
        }
    }
}

#[test]
fn threaded_preview_keeps_boundary_zero_bits() {
    let plan =
        Plan::compile(&Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap())
            .unwrap();
    let mut batch = UnprovedBatch::default();
    let mut w = Workspace::default();
    let points = [[0.0; 3], [-0.0; 3], [60.0, 0.0, 0.0], [1e-300; 3]];
    for workers in [1, 2, 4, 8] {
        for (point, result) in points.iter().zip(
            plan.unproved_preview_batch(&points, workers, &mut batch)
                .unwrap(),
        ) {
            let scalar = plan.unproved_preview(*point, &mut w).unwrap();
            assert_eq!(result.value.to_bits(), scalar.value.to_bits());
            if result.value == 0.0 {
                assert_eq!(result.value.to_bits(), 0.0f64.to_bits());
            }
        }
    }
}

#[test]
fn threaded_preview_leaves_the_callers_proof_allowance_and_answer_unchanged() {
    use catallax_numbers::{Decider, Rational, current_operation_work, with_operation_budget};
    use catallax_shapes::RPoint3;

    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    let points = sample(2_003, [40.0, 30.0, 22.0], 7);
    let point = RPoint3::new(
        Rational::from_i64(10),
        Rational::from_i64(10),
        Rational::from_i64(5),
    );
    let query = || {
        plan.proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
            .unwrap()
    };
    let (reference, baseline) = with_operation_budget(100_000, query);
    let mut batch = UnprovedBatch::default();
    for workers in [1, 2, 4, 8] {
        let (answer, work) = with_operation_budget(100_000, || {
            plan.unproved_preview_batch(&points, workers, &mut batch)
                .unwrap();
            assert_eq!(current_operation_work().used, 0);
            assert!(!current_operation_work().exhausted);
            query()
        });
        assert_eq!(answer, reference);
        assert_eq!(work.used, baseline.used);
        assert_eq!(work.exhausted, baseline.exhausted);
        let (_, zero) = with_operation_budget(0, || {
            plan.unproved_preview_batch(&points, workers, &mut batch)
                .unwrap();
        });
        assert_eq!(zero.used, 0);
        assert!(!zero.exhausted);
    }
}

#[test]
fn wgsl_generation_does_no_proof_work_and_keeps_the_processor_answer() {
    use catallax_numbers::with_operation_budget;
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    let point = [20.0, 15.0, 14.0];
    let before = plan.fast(point, &mut Workspace::default()).unwrap();
    let (source, work) = with_operation_budget(0, || plan.unproved_preview_wgsl().unwrap());
    assert!(!work.exhausted);
    assert_eq!(work.used, 0);
    assert_eq!(source, plan.unproved_preview_wgsl().unwrap());
    assert_eq!(
        before.value.to_bits(),
        plan.fast(point, &mut Workspace::default())
            .unwrap()
            .value
            .to_bits()
    );
}

#[test]
fn wgsl_refuses_constants_outside_binary32() {
    let recipe = serde_json::json!({"format":"catallax-recipe/3","name":"large block","steps":[{"id":"base","op":"block","size":["100000000000000000000000000000000000000000000000000","8","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"base"});
    let plan = Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap();
    assert_eq!(
        plan.unproved_preview_wgsl().unwrap_err().code,
        "field.preview-wgsl-range"
    );
}

#[test]
fn wgsl_refuses_an_unsupported_leaf_without_a_partial_program() {
    let recipe = serde_json::json!({"format":"catallax-recipe/3","name":"draft preview refusal","steps":[{"id":"base","op":"block","size":["8","8","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"drafted","op":"draft","body":"base","faces":[{"face":{"step":"base","role":{"side":{"segment":"x_max"}}}}],"pull":["0","0","1"],"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"45"}],"result":"drafted"});
    let plan = Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap();
    assert_eq!(
        plan.unproved_preview_wgsl().unwrap_err().code,
        "field.preview-wgsl-unsupported"
    );
}
