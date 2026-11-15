//! Qualify a complete replacement on actual immutable STEP sources.
use std::collections::{BTreeMap, BTreeSet};

use catallax_features::{
    Built, ImportFiles, evaluate_with_files,
    import_replace::{FileReplacementOutcome, StepFileReplacement, replace_step_file},
};
use catallax_import::{correspondence::EntityKind, read_step, sha256_hex};
use catallax_numbers::{Decider, Refusal, current_operation_work, with_operation_budget};
use catallax_recipe::{EdgeName, FileRef, Origin, Recipe, Role, VertexName};
use serde_json::{Value, json};

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("import.replacement-control", detail)
}
fn origins(built: &Built) -> BTreeMap<String, Vec<Origin>> {
    built
        .bodies
        .iter()
        .map(|(step, body)| {
            let mut origins = built
                .model
                .faces_of_body(*body)
                .into_iter()
                .flat_map(|face| built.model.face_origins(face).iter().cloned())
                .collect::<Vec<_>>();
            origins.sort();
            (step.as_str().to_owned(), origins)
        })
        .collect()
}
fn run(old: &[u8], new: &[u8], dec: &mut Decider) -> Result<Value, Refusal> {
    let source = read_step(old).map_err(|e| bad(format!("{e:?}")))?;
    let edit = StepFileReplacement {
        old_sha256: sha256_hex(old),
        file: FileRef {
            name: "replacement.step".into(),
            sha256: sha256_hex(new),
        },
        decisions: BTreeMap::new(),
    };
    let steps = source
        .solids
        .iter()
        .enumerate()
        .map(|(i, solid)| {
            json!({"id":format!("source-{i}"),"op":"import_step",
            "file":{"name":"source.step","sha256":edit.old_sha256},"source_solid":solid.name,
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}})
        })
        .collect::<Vec<_>>();
    let recipe = Recipe::parse(&json!({"format":"catallax-recipe/5","name":"file replacement control",
        "result":format!("source-{}",steps.len().checked_sub(1).ok_or_else(||bad("source has no solid"))?),"steps":steps}).to_string())
        .map_err(|e|bad(e.to_string()))?;
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old),
        (edit.file.sha256.clone(), new),
    ]);
    let old_text = recipe.to_text();
    let before = evaluate_with_files(&recipe, &files, dec).map_err(|e| bad(e.to_string()))?;
    let old_origins = origins(&before);
    let before_work = current_operation_work().used;
    let outcome = replace_step_file(&recipe, &edit, &files, dec)?;
    let replacement_work = current_operation_work().used - before_work;
    let FileReplacementOutcome::Ready(result) = outcome else {
        return Err(bad("fixed control requires unresolved selections"));
    };
    if recipe.to_text() != old_text || origins(&before) != old_origins {
        return Err(bad("replacement changed the old revision"));
    }
    let mut checked_faces = 0;
    let mut checked_edges = 0;
    let mut checked_vertices = 0;
    let mut unnamed_edges = 0;
    let mut unnamed_vertices = 0;
    for import in &result.analysis.imports {
        let a = before.bodies[&import.step];
        let b = result.built.bodies[&import.step];
        let faces = import
            .correspondence
            .verified
            .iter()
            .filter(|m| m.old.kind == EntityKind::Face)
            .map(|m| (m.old.key.as_str(), m.new.key.as_str()))
            .collect::<BTreeMap<_, _>>();
        let remap = |origin: &Origin| -> Result<Origin, Refusal> {
            let Role::Imported { key } = &origin.role else {
                return Err(bad("expected a direct imported origin"));
            };
            Ok(Origin {
                step: origin.step.clone(),
                role: Role::Imported {
                    key: faces
                        .get(key.as_str())
                        .ok_or_else(|| bad("unverified actual face origin"))?
                        .to_string(),
                },
            })
        };
        for face in before.model.faces_of_body(a) {
            for origin in before.model.face_origins(face) {
                dec.charge_work(1)?;
                if result
                    .built
                    .model
                    .faces_with_origin(b, &remap(origin)?)
                    .len()
                    != 1
                {
                    return Err(bad(
                        "mapped actual face reference does not resolve uniquely",
                    ));
                }
                checked_faces += 1;
            }
        }
        for edge in before.model.edges_of_body(a) {
            dec.charge_work(1)?;
            let Some(name) = before.model.edge_name(a, edge) else {
                unnamed_edges += 1;
                continue;
            };
            let [x, y] = name.between();
            let mapped = EdgeName::new(remap(x)?, remap(y)?).map_err(|e| bad(e.to_string()))?;
            let expected = before.model.edges_with_name(a, &name).len();
            if result.built.model.edges_with_name(b, &mapped).len() != expected {
                return Err(bad("mapped actual edge references differ"));
            }
            checked_edges += 1;
        }
        let new_names = result
            .built
            .model
            .vertices_of_body(b)
            .into_iter()
            .filter_map(|v| result.built.model.vertex_name(b, v))
            .map(|n| serde_json::to_string(&n).unwrap())
            .collect::<BTreeSet<_>>();
        for vertex in before.model.vertices_of_body(a) {
            dec.charge_work(1)?;
            let Some(name) = before.model.vertex_name(a, vertex) else {
                unnamed_vertices += 1;
                continue;
            };
            let mapped = VertexName::new(name.faces().iter().map(remap).collect::<Result<_, _>>()?)
                .map_err(|e| bad(e.to_string()))?;
            if !new_names.contains(&serde_json::to_string(&mapped).unwrap()) {
                return Err(bad("mapped actual vertex reference is absent"));
            }
            checked_vertices += 1;
        }
    }
    let saved = result.built.recipe().to_text();
    let reopened = Recipe::parse(&saved).map_err(|e| bad(e.to_string()))?;
    let again = evaluate_with_files(&reopened, &files, dec).map_err(|e| bad(e.to_string()))?;
    if origins(&again) != origins(&result.built) {
        return Err(bad("replacement origins differ after reopening"));
    }
    let restored = evaluate_with_files(&recipe, &files, dec).map_err(|e| bad(e.to_string()))?;
    if origins(&restored) != old_origins {
        return Err(bad("retained old revision no longer replays"));
    }
    Ok(
        json!({"status":"replaced","source_sha256":edit.old_sha256,"replacement_sha256":edit.file.sha256,
        "import_steps":result.analysis.imports.len(),"unresolved_references":result.analysis.unresolved_references.len(),
        "checked_model_face_references":checked_faces,"checked_model_edge_names":checked_edges,"checked_model_vertex_names":checked_vertices,
        "unnamed_model_edges":unnamed_edges,"unnamed_model_vertices":unnamed_vertices,
        "replacement_work":replacement_work,"old_revision_preserved":true,"save_reopen_checked":true,"old_revision_replay_checked":true}),
    )
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<_> = std::env::args_os().skip(1).collect();
    if paths.len() != 2 {
        return Err("expected original STEP and changed STEP paths".into());
    }
    let old = std::fs::read(&paths[0])?;
    let new = std::fs::read(&paths[1])?;
    let (result, work) =
        with_operation_budget(20_000_000, || run(&old, &new, &mut Decider::default()));
    let mut value = match result {
        Ok(report) => report,
        Err(e) => json!({"status":"refused","code":e.code,"detail":e.detail}),
    };
    if work.exhausted {
        value = json!({"status":"refused","code":"numbers.operation-budget"});
    }
    value["work_used"] = json!(work.used);
    value["work_limit"] = json!(20_000_000);
    println!("{value}");
    Ok(())
}
