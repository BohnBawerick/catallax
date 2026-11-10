//! The unproved processor preview: its fence and its sign agreement with
//! the deterministic evaluator.

use catallax_field::{Plan, PreviewHardware, Workspace};
use catallax_recipe::Recipe;
use std::path::{Path, PathBuf};

const NAMES: [&str; 3] = ["unproved_preview", "UnprovedValue", "PreviewHardware"];

// Library source that may name the hardware preview. Everything else in the
// workspace's library and binary sources may not, so no export, proven
// measurement or proof route can reach it. A display caller added later is
// listed here in its own reviewed change; `FORBIDDEN` paths never are.
const ALLOWED: [&str; 1] = ["crates/catallax-field/src/lib.rs"];
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
