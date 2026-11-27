use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use catallax_features::{EvaluationOutcome, ImportFiles};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use serde_json::Value;

const EXACT: &str = "
    P0-crossing-pipes P1-bracket P1-drafted-bracket loft-ruled-fillet-r1_4
    fillet-partial-oblique fillet-sphere-cylinder loft-ruled-fillet-r1_10
    fillet-l-mixed-all fillet-vertex-4-all fillet-cross-equal P0-block-sphere
    fillet-cross-unequal fillet-tangent-pocket fillet-box-all P3-filleted-boss
    m8-k12-far-disjoint P5-flanged-housing fillet-tangent-boss
    fillet-cylinder-collapse fillet-cylinder-past-collapse P0-seam
    fillet-roll-over-mitres fillet-sphere-dome fillet-narrow-mitres
    fillet-vertex-4-all-unequal-slopes m10-top-cut fillet-oblique-cone
    fillet-mixed-ordered fillet-box-mitre fillet-unequal-mitres
    fillet-oblique-cylinder fillet-vertex-4-two fillet-unequal-ordered
    fillet-acute-halfdeg fillet-acute-5deg fillet-l-concave fillet-vertex-4-one
    P4-shelled-tray fillet-r12-cube-past-contact fillet-r12-cube-at-contact
    fillet-past-face-limit fillet-at-face-limit P6-lattice-bracket
    fillet-swallowed-face fillet-explicit-setbacks fillet-near-face-limit
    fillet-box-single holes-P1-bolt-holes holes-P5-countersinks
    fillet-torus-intersection fillet-obtuse-174deg fillet-obtuse-179deg
    fillet-torus-slice P0-line-tangent fillet-boss-seam-rotated fillet-hole-mouth
    fillet-boss-rim P0-near-miss fillet-cone-base fillet-cone-top
    fillet-narrow-cone P0-equal-cylinder-sphere P2-collar fillet-variable-linear
    fillet-variable-zero-end fillet-cylinder-rim
";
const FITTED: &str = "
    rim-five-r7.99 rim-five-r7.9 rim-five-r5 rim-six-r3.99 rim-six-r2
    rim-six-r1_2 rim-six-r1 loft-smooth-fillet-r1_10 rim-five-r1
";
const HELD_OUT: &str = "
    abc-00000858 abc-00001883 abc-00001130 abc-00001664 abc-00001285
    abc-00001344 abc-00001237 abc-00001453 abc-00001534 abc-00001427
    abc-00001260 abc-00001533 abc-00001006 abc-00000938 abc-00001339
";

pub fn bodies() -> impl Iterator<Item = (&'static str, EvaluationOutcome)> {
    let root = std::env::var_os("CARGO_MANIFEST_DIR")
        .map_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")), PathBuf::from)
        .join("../..");
    let manifest: Value =
        serde_json::from_str(include_str!("../../../../tools/predictor/corpus.json")).unwrap();
    let sources = match (
        std::env::var_os("P5_SOURCES"),
        std::env::var_os("P5_SELECTION"),
    ) {
        (Some(sources), Some(selection)) => {
            let members: Vec<(String, PathBuf)> =
                serde_json::from_slice(&fs::read(selection).unwrap()).unwrap();
            Some((
                PathBuf::from(sources),
                members.into_iter().collect::<BTreeMap<_, _>>(),
            ))
        }
        (None, None) => None,
        _ => panic!("set both P5_SOURCES and P5_SELECTION"),
    };
    let public = EXACT.split_whitespace().count() + FITTED.split_whitespace().count();
    assert_eq!(public, 75);
    assert_eq!(HELD_OUT.split_whitespace().count(), 15);
    let enabled = sources.is_some();
    let held_out = HELD_OUT.split_whitespace().filter(move |_| enabled);
    EXACT
        .split_whitespace()
        .chain(held_out)
        .map(|id| (id, false))
        .chain(FITTED.split_whitespace().map(|id| (id, true)))
        .map(move |(id, fitted)| {
            let row = manifest["rows"]
                .as_array()
                .unwrap()
                .iter()
                .find(|row| row["id"] == id)
                .unwrap();
            let source = row["import"]["model"].as_str().map(|model| {
                let (root, files) = sources.as_ref().unwrap();
                let path = files
                    .get(model)
                    .unwrap_or_else(|| panic!("{id}: no selected STEP source"));
                let bytes = fs::read(root.join(path))
                    .unwrap_or_else(|error| panic!("{id}: cannot read STEP source: {error}"));
                (row["import"]["sha256"].as_str().unwrap(), bytes)
            });
            let recipe = recipe(&root, row, source.as_ref().map(|(hash, _)| *hash));
            let files = source
                .iter()
                .map(|(hash, bytes)| ((*hash).to_owned(), bytes.as_slice()))
                .collect::<ImportFiles<'_>>();
            let built = std::thread::scope(|scope| {
                std::thread::Builder::new()
                    .stack_size(8 * 1024 * 1024)
                    .spawn_scoped(scope, || {
                        catallax_features::evaluate_incremental_outcome(
                            &recipe,
                            None,
                            &files,
                            &mut Decider::default(),
                        )
                    })
                    .unwrap()
                    .join()
                    .unwrap()
            })
            .unwrap_or_else(|error| panic!("{id}: {error:?}"));
            assert_eq!(
                matches!(&built, EvaluationOutcome::Numerical(_)),
                fitted,
                "{id}: corpus build grade changed"
            );
            (id, built)
        })
}

fn recipe(root: &Path, row: &Value, hash: Option<&str>) -> Recipe {
    let mut value = if let Some(hash) = hash {
        serde_json::json!({
            "format": "catallax-recipe/3",
            "name": row["id"],
            "steps": [{"id": "imported", "op": "import_step",
                "file": {"name": "part.step", "sha256": hash},
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}}],
            "result": "imported"
        })
    } else if let Some(path) = row["recipe"].as_str() {
        serde_json::from_slice(&fs::read(root.join(path)).unwrap()).unwrap()
    } else if let Some(path) = row["cases"].as_str() {
        fs::read_to_string(root.join(path))
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| serde_json::from_str::<Value>(line).unwrap())
            .find(|case| case["id"] == row["case"])
            .unwrap()["recipe"]
            .clone()
    } else {
        assert!(row["inline"].is_object());
        row["inline"].clone()
    };
    if let Some(steps) = row["append"].as_array() {
        value["steps"]
            .as_array_mut()
            .unwrap()
            .extend(steps.iter().cloned());
        value["result"] = row["result"].clone();
    }
    if let Some(format) = row["format"].as_str() {
        value["format"] = Value::from(format);
    }
    Recipe::parse(&value.to_string()).unwrap()
}
