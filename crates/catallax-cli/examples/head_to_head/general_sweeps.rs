//! Ordinary helical recipes use numerical application results, never exact checkpoints.
use super::*;
use catallax_app::numerical_sweeps::SweepBody;
use catallax_numbers::Rational;
use catallax_topology::numerical::NumericalOptions;
use std::collections::BTreeSet;

fn consumers(body: Option<&SweepBody>, out: &Path, limit: u64) -> io::Result<Value> {
    let Some(body) = body else {
        return Ok(
            json!({"preview":{"status":"not-run","reason":"build required"},
            "display":{"status":"not-run","reason":"build required"},
            "measurement":{"status":"not-run","reason":"build required"},
            "step":{"status":"not-run","reason":"build required"}}),
        );
    };
    let (preview, mut preview_report) = observe(limit, || {
        body.preview(&Default::default(), &mut Decider::default())
    });
    preview_report["proof_state"] = json!("unproved");
    if let Some(preview) = preview {
        preview_report["vertices"] = json!(preview.positions.len());
        preview_report["triangles"] = json!(preview.triangles.len());
    }
    let (mesh, mut display) = observe(limit, || {
        body.display(&Default::default(), &mut Decider::default())
    });
    display["proof"] = json!("numerical");
    display["scope"] =
        json!("display bound concerns the fitted representation, not the authored sweep");
    if let Some(mesh) = mesh {
        display["vertices"] = json!(mesh.fitted_mesh.positions.len());
        display["triangles"] = json!(mesh.fitted_mesh.triangles.len());
        mesh.fitted_mesh
            .write_obj(&mut fs::File::create(out.join("display.obj"))?)?;
    }
    let (properties, mut measurement) = observe(limit, || {
        body.properties(&Accuracy::default(), &mut Decider::default())
    });
    measurement["proof"] = json!("numerical");
    measurement["scope"] = json!("integration enclosures concern the fitted representation only");
    measurement["original_error_bound"] = json!("unestablished");
    if let Some(properties) = properties {
        measurement["accuracy_reached"] = json!(
            properties
                .marks
                .iter()
                .all(|m| m.code != "measure.accuracy-not-reached")
        );
        measurement["values"] = json!({"volume":enclosure(&properties.fitted.volume),
            "area":enclosure(&properties.fitted.area),
            "centroid":properties.fitted.centroid.as_ref().map(|c|c.iter().map(enclosure).collect::<Vec<_>>())});
        measurement["marks"] = json!(
            properties
                .marks
                .iter()
                .map(|m| json!({"code":m.code,"detail":m.detail}))
                .collect::<Vec<_>>()
        );
    }
    let (artifact, mut export) = observe(limit, || {
        catallax_step::numerical::build(body.body(), &Default::default(), &mut Decider::default())
    });
    export["proof"] = json!("numerical");
    export["exact"] = json!(false);
    export["expected_faces"] = json!(body.body().counts().faces);
    export["tolerance_mm"] = json!(body.body().options().tolerance_mm.to_string());
    if let Some(artifact) = artifact {
        let path = out.join("answer.step");
        fs::write(&path, artifact.bytes())?;
        let bytes = fs::read(&path)?;
        if let Err(error) = artifact.check_written(&bytes, &mut Decider::default()) {
            export["status"] = json!("refused");
            export["error"] = json!(error.to_string());
            fs::remove_file(path)?;
        }
    }
    Ok(json!({"preview":preview_report,"display":display,"measurement":measurement,"step":export}))
}

pub(super) fn run(
    original: &Recipe,
    edited: &Recipe,
    input: &Value,
    root: &Path,
    limit: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let tolerance = input["numerical_tolerance_mm"]
        .as_str()
        .map(Rational::parse)
        .transpose()?
        .unwrap_or_else(|| NumericalOptions::default().tolerance_mm);
    let mut report = json!({"format":"catallax-head-to-head-phases/1","variants":{},
        "proof":"numerical","exact":false,"tolerance_mm":tolerance.to_string(),
        "edit_mode":"complete numerical recipe rebuild, no incremental reuse",
        "measurement_accuracy":{"relative":"1e-12","absolute":"1e-30","scope":"fitted representation only"}});
    let mut origins = Vec::new();
    for (label, recipe) in [("construction", original), ("edit", edited)] {
        let out = root.join(label);
        fs::create_dir_all(&out)?;
        fs::write(out.join("recipe.json"), serde_json::to_vec_pretty(recipe)?)?;
        let (body, mut build) = observe(limit, || {
            catallax_app::numerical_sweeps::build(
                recipe,
                NumericalOptions {
                    tolerance_mm: tolerance.clone(),
                    vertex_tolerance_mm: None,
                },
                &mut Decider::default(),
            )
        });
        build["proof"] = json!("numerical");
        build["exact"] = json!(false);
        let names = body.as_ref().map(|b| {
            b.face_origins()
                .values()
                .map(|origins| json!({"origins":origins}))
                .collect::<Vec<_>>()
        });
        let checks = body
            .as_ref()
            .map(|b| serde_json::to_value(b.checks()).expect("finite sweep checks"));
        report["variants"][label] =
            json!({"build":build,"names":{"faces":names},"sweep_fitting":checks});
        fs::write(
            root.join("phases.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        let values = consumers(body.as_ref(), &out, limit)?;
        for (key, value) in values.as_object().expect("consumer object") {
            report["variants"][label][key] = value.clone();
        }
        report["variants"][label]["acceptance"] = acceptance(&report["variants"][label]);
        origins.push(body.as_ref().map(|b| {
            b.face_origins()
                .values()
                .flatten()
                .map(|o| o.key())
                .collect::<BTreeSet<_>>()
        }));
        fs::write(
            root.join("phases.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
    }
    report["edit_survival"] = if let (Some(before), Some(after)) = (&origins[0], &origins[1]) {
        let missing: Vec<_> = before.difference(after).collect();
        json!({"status":if missing.is_empty(){"passed"}else{"failed"},
            "coverage":"every original face origin resolves after edit; internal chart edges are not persistent edge names",
            "missing_faces":missing})
    } else {
        json!({"status":"not-run","reason":"both complete numerical builds required"})
    };
    fs::write(
        root.join("phases.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!("{}", serde_json::to_string(&report)?);
    Ok(())
}
