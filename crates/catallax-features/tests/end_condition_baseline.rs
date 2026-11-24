//! Differential oracle for the first-hit recipes shipped before face starts.
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use catallax_features::evaluate;
use catallax_import::sha256_hex;
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;

const ALLOWANCE: u64 = 100_000_000;

fn cases() -> BTreeMap<String, PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut cases = BTreeMap::new();
    for group in ["curved-limits", "polygon-limits", "circular-sphere-limits"] {
        let directory = root.join("tests/fixtures").join(group);
        for entry in directory.read_dir().unwrap() {
            let path = entry.unwrap().path();
            if path.extension().is_some_and(|ext| ext == "json") {
                let name = format!("{group}/{}", path.file_name().unwrap().to_str().unwrap());
                assert!(cases.insert(name, path).is_none());
            }
        }
    }
    cases.insert(
        "stepped-limit".into(),
        root.join("tests/fixtures/stepped-limit.json"),
    );
    cases.insert(
        "application-next".into(),
        root.join("../catallax-app/tests/fixtures/face-extrusion-next.json"),
    );
    cases
}

fn signature(recipe: &Recipe) -> String {
    let (built, work) =
        with_operation_budget(ALLOWANCE, || evaluate(recipe, &mut Decider::default()));
    let used = work.used;
    let built = match built {
        Ok(built) => built,
        Err(catallax_features::BuildError::Geometry { refusal, .. }) => {
            return format!("{used} refusal {}", refusal.code);
        }
        Err(error) => panic!("{}: {error}", recipe.name),
    };
    let model = &built.model;
    let body = built.result;
    assert!(model.check(body).is_empty(), "{}", recipe.name);
    let mut faces: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .map(|id| {
            let face = model.face(id);
            format!(
                "{:?} {:?}",
                model.geom().surface(face.surface).kind(),
                model.face_origins(id)
            )
        })
        .collect();
    faces.sort_unstable();
    let topology = sha256_hex(format!("{:?} {faces:?}", model.body(body).kind).as_bytes());
    let mut step = Vec::new();
    let (export, work) = with_operation_budget(ALLOWANCE, || {
        catallax_step::write_step(
            model,
            &[body],
            &Default::default(),
            &mut Decider::default(),
            &mut step,
        )
    });
    let exchange = match export {
        Ok(_) => sha256_hex(&step),
        Err(error) if !work.exhausted => error.code.to_owned(),
        Err(_) => "numbers.operation-budget".to_owned(),
    };
    format!("{used} {topology} {exchange}")
}

fn axial_wall_refusal() -> Recipe {
    let mut value: serde_json::Value = serde_json::from_str(include_str!(
        "fixtures/circular-sphere-limits/sphere-entry-to.json"
    ))
    .unwrap();
    value["steps"][1]["op"] = serde_json::json!("cylinder");
    value["steps"][1]["height"] = serde_json::json!("3");
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn first_hit_matches_main() {
    let cases = cases();
    let expected = include_str!("fixtures/end-condition-baseline.txt");
    for line in expected.lines().filter(|line| !line.starts_with('#')) {
        let (name, expected) = line.split_once(' ').unwrap();
        let recipe = if name == "axial-cylinder-wall" {
            axial_wall_refusal()
        } else {
            Recipe::parse(&std::fs::read_to_string(&cases[name]).unwrap()).unwrap()
        };
        assert_eq!(signature(&recipe), expected, "{name}");
    }
}
