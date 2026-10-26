//! Development-only phase measurements. Never linked into an application entry.
use catallax_features::Built;
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Decider, with_operation_budget, with_operation_log};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{env, fs, io, path::Path, time::Instant};

fn enclosure(value: &Enclosure) -> Value {
    json!([value.lo.to_string(), value.hi.to_string()])
}

fn observe<T, E: std::fmt::Display>(
    limit: u64,
    run: impl FnOnce() -> Result<T, E>,
) -> (Option<T>, Value) {
    let start = Instant::now();
    let ((result, work), operations) = with_operation_log(|| with_operation_budget(limit, run));
    let wall_ms = start.elapsed().as_secs_f64() * 1000.0;
    let error = result.as_ref().err().map(ToString::to_string);
    let operations: Vec<_> = operations
        .iter()
        .map(|op| {
            json!({
                "kernel":op.kernel,"id":op.id,"operation":op.operation,
                "work":op.work_after-op.work_before,"wall_ms":op.wall_ms,
                "refusal":op.refusal.as_ref().map(|r| json!({"code":r.code,"detail":r.detail}))
            })
        })
        .collect();
    let status = if work.exhausted {
        "work-exhausted"
    } else if error.is_some() {
        "refused"
    } else {
        "answered"
    };
    let report = json!({"status":status,"wall_ms":wall_ms,"work":work.used,
        "work_limit":limit,"error":error,"operations":operations,
        "refusal_ms":if error.is_some() || work.exhausted {Some(wall_ms)} else {None}});
    (if work.exhausted { None } else { result.ok() }, report)
}

fn names(built: &Built) -> Value {
    let faces: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|face| json!({"origins":built.model.face_origins(face)}))
        .collect();
    let edges: Vec<_> = built
        .model
        .edges_of_body(built.result)
        .into_iter()
        .map(|edge| built.model.edge_name(built.result, edge))
        .collect();
    json!({"faces":faces,"edges":edges})
}

fn survival(before: Option<&Built>, after: Option<&Built>) -> Value {
    let (Some(before), Some(after)) = (before, after) else {
        return json!({"status":"not-run","reason":"both complete builds required"});
    };
    let mut missing_faces = Vec::new();
    let mut missing_edges = Vec::new();
    let mut unnamed = 0;
    for built in [before, after] {
        unnamed += built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .filter(|&f| built.model.face_origins(f).is_empty())
            .count();
        unnamed += built
            .model
            .edges_of_body(built.result)
            .into_iter()
            .filter(|&e| built.model.edge_name(built.result, e).is_none())
            .count();
    }
    for face in before.model.faces_of_body(before.result) {
        for origin in before.model.face_origins(face) {
            if after
                .model
                .faces_with_origin(after.result, origin)
                .is_empty()
            {
                missing_faces.push(origin);
            }
        }
    }
    for edge in before.model.edges_of_body(before.result) {
        if let Some(name) = before.model.edge_name(before.result, edge)
            && after.model.edges_with_name(after.result, &name).is_empty()
        {
            missing_edges.push(name);
        }
    }
    json!({"status":if unnamed == 0 && missing_faces.is_empty() && missing_edges.is_empty() {"passed"} else {"failed"},
        "coverage":"every original result face origin and physical-edge name resolves after edit; not a semantic correspondence proof",
        "missing_faces":missing_faces,"missing_edges":missing_edges,"unnamed":unnamed})
}

fn consumers(built: Option<&Built>, out: &Path, limit: u64) -> io::Result<Value> {
    let Some(built) = built else {
        return Ok(
            json!({"preview":{"status":"not-run","reason":"build required"},
            "display":{"status":"not-run","reason":"build required"},
            "measurement":{"status":"not-run","reason":"build required"},"step":{"status":"not-run","reason":"build required"}}),
        );
    };
    // The unproved preview carries the selection names of this revision. It is
    // timed on its own allowance and never stands in for the certified display.
    let (preview, mut preview_report) = observe(limit, || {
        catallax_app::preview_body(built, &Default::default(), &mut Decider::default())
    });
    preview_report["proof_state"] = json!("unproved");
    if let Some(preview) = &preview {
        preview_report["vertices"] = json!(preview.mesh.positions.len());
        preview_report["triangles"] = json!(preview.mesh.triangles.len());
        preview_report["named_faces"] = json!(preview.names.faces.len());
    }
    let (mesh, mut display) = observe(limit, || {
        catallax_mesh::mesh_body(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
    });
    if let Some(mesh) = mesh {
        display["vertices"] = json!(mesh.value.positions.len());
        display["triangles"] = json!(mesh.value.triangles.len());
        mesh.value
            .write_obj(&mut fs::File::create(out.join("display.obj"))?)?;
    }
    let (properties, mut measurement) = observe(limit, || {
        catallax_measure::volume_area_centroid(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
    });
    if let Some(properties) = properties {
        let accuracy_reached = properties
            .marks
            .iter()
            .all(|mark| mark.code != "measure.accuracy-not-reached");
        measurement["accuracy_reached"] = json!(accuracy_reached);
        measurement["values"] = json!({"volume":enclosure(&properties.value.volume),
            "area":enclosure(&properties.value.area),
            "centroid":properties.value.centroid.as_ref().map(|c| c.iter().map(enclosure).collect::<Vec<_>>())});
        measurement["marks"] = json!(
            properties
                .marks
                .iter()
                .map(|m| json!({"code":m.code,"detail":m.detail}))
                .collect::<Vec<_>>()
        );
    }
    let mut step_bytes = Vec::new();
    let (step, mut export) = observe(limit, || {
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut step_bytes,
        )
    });
    export["expected_faces"] = json!(step.as_ref().map(|answer| answer.value.faces));
    if step.is_some() {
        fs::write(out.join("answer.step"), step_bytes)?;
    }
    Ok(json!({"preview":preview_report,"display":display,"measurement":measurement,"step":export}))
}

const PREVIEW_TARGET_MS: f64 = 1000.0;
const PROOF_TARGET_MS: f64 = 2000.0;

fn answered_ms(phase: &Value) -> Option<f64> {
    (phase["status"] == "answered")
        .then(|| phase["wall_ms"].as_f64())
        .flatten()
}

/// Separates the wait for the unproved preview from background proof completion.
/// A phase that did not answer yields no wait, so a fast refusal is never a pass.
/// The certified display, requested measurement and STEP export are proof phases.
fn acceptance(variant: &Value) -> Value {
    let preview_wait = answered_ms(&variant["build"])
        .zip(answered_ms(&variant["preview"]))
        .map(|(build, preview)| build + preview);
    let proof_parts: Vec<_> = ["display", "measurement", "step"]
        .iter()
        .map(|phase| answered_ms(&variant[*phase]))
        .collect();
    let proof_wait = answered_ms(&variant["build"]).and_then(|build| {
        proof_parts
            .iter()
            .try_fold(build, |sum, part| part.map(|ms| sum + ms))
    });
    json!({
        "preview":{"state":"unproved","wait_ms":preview_wait,"target_ms":PREVIEW_TARGET_MS,
            "wait":"build plus preview display and selection names",
            "met":preview_wait.map(|ms| ms <= PREVIEW_TARGET_MS)},
        "proof":{"wait_ms":proof_wait,"target_ms":PROOF_TARGET_MS,
            "wait":"build plus certified display, requested measurement and STEP export, run in sequence",
            "met":proof_wait.map(|ms| ms <= PROOF_TARGET_MS)},
        "built":proof_wait.is_some_and(|ms| ms <= PROOF_TARGET_MS),
        "built_rule":"a part counts as built only after its proof phases all answer within the proof target; a preview alone never does",
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err(io::Error::other("usage: head_to_head REQUEST.json OUT_DIR").into());
    }
    let input: Value = serde_json::from_slice(&fs::read(&args[0])?)?;
    let limit = input["max_work"].as_u64().ok_or("max_work required")?;
    let original = Recipe::parse(&serde_json::to_string(&input["recipe"])?)?;
    let edited = Recipe::parse(&serde_json::to_string(&input["edited_recipe"])?)?;
    let root = Path::new(&args[1]);
    fs::create_dir_all(root)?;
    let mut report = json!({"format":"catallax-head-to-head-phases/1","variants":{},
        "edit_mode":"complete recipe rebuild, no incremental reuse","measurement_accuracy":{"relative":"1e-12","absolute":"1e-30"}});
    let mut builds = Vec::new();
    for (label, recipe) in [("construction", &original), ("edit", &edited)] {
        let out = root.join(label);
        fs::create_dir_all(&out)?;
        let (built, mut build) = observe(limit, || {
            catallax_features::evaluate(recipe, &mut Decider::default())
        });
        build["operations_rebuilt"] = json!(
            build["operations"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|op| {
                    op["refusal"].is_null()
                        && recipe.steps.iter().any(|step| {
                            step.id.as_str() == op["id"].as_str().unwrap_or("")
                                && step.op.name() == op["operation"].as_str().unwrap_or("")
                        })
                })
                .count()
        );
        let build_names = built.as_ref().map(names);
        report["variants"][label] = json!({"build":build,"names":build_names});
        fs::write(
            root.join("phases.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        // Each consumer gets its own fixed work allowance, not construction's remainder.
        let values = consumers(built.as_ref(), &out, limit)?;
        for (key, value) in values.as_object().unwrap() {
            report["variants"][label][key] = value.clone();
        }
        report["variants"][label]["acceptance"] = acceptance(&report["variants"][label]);
        fs::write(
            root.join("phases.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        builds.push(built);
    }
    report["edit_survival"] = survival(builds[0].as_ref(), builds[1].as_ref());
    fs::write(
        root.join("phases.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!("{}", serde_json::to_string(&report)?);
    Ok(())
}
