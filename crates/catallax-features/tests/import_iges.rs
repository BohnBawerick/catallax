//! `import_iges` recipe steps: IGES solids enter a recipe as STEP files do,
//! with the same face names, selections, field compilation and replacement.

use std::collections::BTreeMap;

use catallax_features::{
    ImportFiles, evaluate_with_files,
    import_edit::{StepReplacement, replace_step},
};
use catallax_iges::{IgesOptions, write_iges};
use catallax_import::{Surface, read_iges, sha256_hex};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{FileRef, Recipe, Role, StepId};
use catallax_shapes::Pose;
use catallax_topology::{Model, block};
use serde_json::{Value, json};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}

/// A 2 by 3 by 5 block written by Catallax's IGES writer.
fn catallax_box() -> Vec<u8> {
    let mut model = Model::new();
    let body = block(&mut model, [2, 3, 5].map(r), &Pose::identity()).unwrap();
    let mut bytes = Vec::new();
    write_iges(
        &model,
        &[body],
        &IgesOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
}

/// The same block written by OpenCASCADE 8.0 in B-rep mode.
fn occt_box() -> Vec<u8> {
    include_bytes!("../../catallax-import/tests/fixtures/iges/occt-box.igs").to_vec()
}

fn recipe(bytes: &[u8], extra: Vec<Value>, result: &str) -> Recipe {
    let mut steps = vec![json!({"id": "source", "op": "import_iges",
        "file": {"name": "part.igs", "sha256": sha256_hex(bytes)},
        "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["7", "-2", "3"]}})];
    steps.extend(extra);
    Recipe::parse(
        &json!({"format": "catallax-recipe/5", "name": "IGES import", "steps": steps,
            "result": result})
        .to_string(),
    )
    .unwrap()
}

/// The face key of the face whose outward normal is `+z`.
fn top_key(bytes: &[u8]) -> String {
    read_iges(bytes).unwrap().solids[0]
        .faces
        .iter()
        .find(|f| match &f.surface {
            Surface::Plane(p) => {
                let up = p.axis[2].sign() == catallax_numbers::Sign::Positive;
                p.axis[0].is_zero() && p.axis[1].is_zero() && (up == f.same_sense)
            }
            _ => false,
        })
        .unwrap()
        .name
        .clone()
}

fn shelled(bytes: &[u8]) -> Recipe {
    let key = top_key(bytes);
    recipe(
        bytes,
        vec![
            json!({"id": "hollow", "op": "shell", "body": "source", "thickness": "1/2",
            "open": [{"face": {"step": "source", "role": {"imported": {"key": key}}}}]}),
        ],
        "hollow",
    )
}

fn measured(built: &catallax_features::Built) -> (Rational, Rational) {
    let v = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    (v.lo, v.hi)
}

#[test]
fn an_iges_block_builds_placed_with_its_directory_names() {
    for bytes in [catallax_box(), occt_box()] {
        let files = ImportFiles::from([(sha256_hex(&bytes), &bytes[..])]);
        let built = evaluate_with_files(
            &recipe(&bytes, vec![], "source"),
            &files,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(measured(&built), (r(30), r(30)));
        let source = read_iges(&bytes).unwrap();
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 6);
        for face in faces {
            let Role::Imported { key } = &built.model.face_origins(face)[0].role else {
                panic!("an imported face")
            };
            assert!(key.starts_with('D'));
            assert!(source.solids[0].faces.iter().any(|f| &f.name == key));
        }
        // A downstream selection names an imported face by its key.
        let hollow =
            evaluate_with_files(&shelled(&bytes), &files, &mut Decider::default()).unwrap();
        assert!(measured(&hollow).1 < r(30));
    }
}

#[test]
fn the_field_compiles_an_iges_import() {
    let bytes = catallax_box();
    let files = ImportFiles::from([(sha256_hex(&bytes), &bytes[..])]);
    catallax_field::Plan::compile_with_files_checked(
        &recipe(&bytes, vec![], "source"),
        &files,
        &mut Decider::default(),
    )
    .unwrap();
}

#[test]
fn replacing_an_iges_file_rewrites_every_reference() {
    let (old, new) = (catallax_box(), occt_box());
    let recipe = shelled(&old);
    // Pair faces by outward normal: the editorial relation, stated explicitly.
    let normals = |bytes: &[u8]| -> BTreeMap<String, String> {
        read_iges(bytes).unwrap().solids[0]
            .faces
            .iter()
            .map(|f| {
                let Surface::Plane(p) = &f.surface else {
                    unreachable!()
                };
                let sign = if f.same_sense { 1 } else { -1 };
                let n = p.axis.clone().map(|c| {
                    let unit = match c.sign() {
                        catallax_numbers::Sign::Positive => 1,
                        catallax_numbers::Sign::Negative => -1,
                        _ => 0,
                    };
                    (unit * sign).to_string()
                });
                (n.join(","), f.name.clone())
            })
            .collect()
    };
    let (before, after) = (normals(&old), normals(&new));
    let faces = before
        .iter()
        .map(|(normal, key)| (key.clone(), vec![after[normal].clone()]))
        .collect();
    let edit = StepReplacement {
        step: StepId::new("source").unwrap(),
        old_sha256: sha256_hex(&old),
        file: FileRef {
            name: "replacement.igs".into(),
            sha256: sha256_hex(&new),
        },
        decisions: BTreeMap::new(),
        faces,
    };
    let files = ImportFiles::from([(sha256_hex(&old), &old[..]), (sha256_hex(&new), &new[..])]);
    let previous = recipe.to_text();
    let built = replace_step(&recipe, &edit, &files, &mut Decider::default()).unwrap();
    assert_eq!(recipe.to_text(), previous);
    let text = built.recipe().to_text();
    assert!(text.contains(&sha256_hex(&new)));
    assert!(text.contains(&format!("\"{}\"", top_key(&new))));
    assert!(text.contains("import_iges"));
    let original = evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    assert_eq!(measured(&built), measured(&original));
    // An unmapped key fails loudly.
    let mut partial = edit.clone();
    partial.faces.remove(&top_key(&old));
    let error = replace_step(&recipe, &partial, &files, &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "import.replace-unmapped");
}
