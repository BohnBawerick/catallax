//! Fixed material-edit qualification; never extends the numerical edit family.
use catallax_features::{
    evaluate_numerical,
    numerical::{NumericalEdit, NumericalOperation, NumericalRecipe},
};
use catallax_import::{
    numerical::{self, StepImport},
    sha256_hex,
};
use catallax_numbers::{
    Decider, Rational, Real, Refusal, current_operation_work, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Point3, Pose, RVector3, Rotation, Surface};
use catallax_step::numerical::{NumericalStepOptions, build as export};
use catallax_topology::numerical::{NumericalBody, NumericalOptions};
use serde_json::{Value, json};
use std::{io::Write, path::Path};

fn problem(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

fn state(body: &NumericalBody, dec: &mut Decider) -> Result<String, Refusal> {
    let value = format!(
        "{:?}\n{:?}\n{:?}\n{:?}",
        body.definition(),
        body.source_definition(),
        body.history(),
        body.evidence()
    );
    dec.charge_work(value.len() as u64)?;
    Ok(sha256_hex(value.as_bytes()))
}

fn posed(body: &NumericalBody, dec: &mut Decider) -> Result<NumericalBody, Refusal> {
    let pose = Pose::from_rational(
        Rotation::identity(),
        RVector3::new(3.into(), 5.into(), 7.into()),
    );
    let mut definition = body.definition().clone();
    for vertex in &mut definition.vertices {
        dec.charge_work(1)?;
        let [x, y, z] = vertex.position.clone().map(Real::from);
        let p = pose.apply_point(&Point3::new(x, y, z));
        vertex.position = [&p.x, &p.y, &p.z].map(|q| q.as_rational().expect("rational pose"));
    }
    for edge in &mut definition.edges {
        dec.charge_work(1)?;
        let Curve::Line(line) = &mut edge.curve else {
            return Err(problem(
                "import.material-control",
                "straight boundaries required",
            ));
        };
        line.origin = pose.apply_point(&line.origin);
    }
    for face in &mut definition.faces {
        dec.charge_work(1)?;
        let Surface::Plane(plane) = &mut face.surface else {
            return Err(problem("import.material-control", "planar source required"));
        };
        plane.origin = pose.apply_point(&plane.origin);
    }
    body.rebuild(definition, body.options().clone(), dec)
}

fn resize(
    body: &NumericalBody,
    axis: usize,
    length: Rational,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    let source = NumericalRecipe {
        source: body.definition().clone(),
        sew: Vec::new(),
        options: body.options().clone(),
        edits: Vec::new(),
    };
    let old = evaluate_numerical(&source, dec)?;
    let next = old.edit(
        NumericalEdit {
            id: "material-dimension".into(),
            options: body.options().clone(),
            operation: NumericalOperation::Resize {
                axis,
                length_mm: length,
            },
        },
        dec,
    )?;
    // Preserve the actual protected input's prior checked states rather than
    // treating a freshly evaluated edit recipe as its complete provenance.
    body.rebuild(
        next.body().definition().clone(),
        next.body().options().clone(),
        dec,
    )
}

fn save(body: &NumericalBody, path: &Path, dec: &mut Decider) -> Result<Value, Refusal> {
    let file = export(body, &NumericalStepOptions::default(), dec)?;
    let mut output = std::fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(path)
        .map_err(|e| problem("import.material-io", e.to_string()))?;
    output
        .write_all(file.bytes())
        .map_err(|e| problem("import.material-io", e.to_string()))?;
    drop(output);
    let bytes = std::fs::read(path).map_err(|e| problem("import.material-io", e.to_string()))?;
    file.check_written(&bytes, dec)?;
    Ok(
        json!({"file":path.file_name().unwrap().to_string_lossy(),"sha256":sha256_hex(&bytes),"category":"numerical approximate",
        "history_stages":body.history().len(),"repairs":body.evidence().repairs().len(),"checks":body.evidence().checks().len(),
        "source_history_stages":file.report().source_history.len(),"complete_body_bound":"unestablished"}),
    )
}

fn positive(out: &Path, dec: &mut Decider) -> Result<Value, Refusal> {
    let mut recipe = NumericalRecipe::parse(
        include_str!("../../../evidence/phase1/N04/recipe.json"),
        dec,
    )?;
    recipe.edits.clear();
    let before = evaluate_numerical(&recipe, dec)?;
    let before_file = save(before.body(), &out.join("control-before.step"), dec)?;
    let after = before.edit(
        NumericalEdit {
            id: "double-height".into(),
            options: before.body().options().clone(),
            operation: NumericalOperation::Resize {
                axis: 2,
                length_mm: 2.into(),
            },
        },
        dec,
    )?;
    if after.body().history().len() != before.body().history().len() + 1 {
        return Err(problem(
            "import.material-history",
            "positive edit lost a prior stage",
        ));
    }
    let after_file = save(after.body(), &out.join("control-after.step"), dec)?;
    let saved = after.recipe().to_json()?;
    let replayed = catallax_features::replay_numerical(&saved, dec)?;
    if state(replayed.body(), dec)? != state(after.body(), dec)? {
        return Err(problem(
            "import.material-replay",
            "saved material edit differs on replay",
        ));
    }
    Ok(
        json!({"status":"numerical","exact":false,"expected_before_volume_mm3":"12","expected_after_volume_mm3":"24","expected_delta_mm3":"12",
        "history_retained":true,"replay_checked":true,"before":before_file,"after":after_file}),
    )
}

fn run(bytes: &[u8], out: &Path, dec: &mut Decider) -> Result<Value, Refusal> {
    let start = current_operation_work().used;
    let controls: Value = serde_json::from_str(include_str!(
        "../../../docs/phase1/import-oblique-step/real-controls.json"
    ))
    .expect("committed controls");
    let expected = controls["rows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["model_id"] == "00000022")
        .unwrap();
    let sha = sha256_hex(bytes);
    if sha != expected["source_sha256"].as_str().unwrap() {
        return Err(problem(
            "import.material-source",
            "expected fixed census source 00000022",
        ));
    }
    let imported = numerical::import_step(
        bytes,
        &sha,
        &StepId::new("source").unwrap(),
        &NumericalOptions::default(),
        dec,
    )
    .map_err(|e| problem("import.material-admission", format!("{e:?}")))?;
    let import_work = current_operation_work().used - start;
    let StepImport::Numerical(imported) = imported else {
        return Err(problem(
            "import.material-category",
            "fixed source must remain numerical",
        ));
    };
    let derivative = imported.record().derivative.as_ref().ok_or_else(|| {
        problem(
            "import.material-source",
            "accepted repair lacks a derivative",
        )
    })?;
    if derivative.sha256 != expected["repaired_definition_sha256"].as_str().unwrap() {
        return Err(problem(
            "import.material-source",
            "accepted source definition differs from writer qualification",
        ));
    }
    let before_work = current_operation_work().used;
    let protected = posed(&imported.bodies()[0].1, dec)?;
    let pose_work = current_operation_work().used - before_work;
    let fingerprint = state(&protected, dec)?;
    let before_work = current_operation_work().used;
    let before_file = save(&protected, &out.join("real-before.step"), dec)?;
    let before_step_work = current_operation_work().used - before_work;
    if before_file["sha256"] != expected["edited_step"]["sha256"] {
        return Err(problem(
            "import.material-source",
            "posed input differs from writer qualification",
        ));
    }
    let lo = protected
        .definition()
        .vertices
        .iter()
        .map(|v| v.position[0].clone())
        .min()
        .unwrap();
    let hi = protected
        .definition()
        .vertices
        .iter()
        .map(|v| v.position[0].clone())
        .max()
        .unwrap();
    dec.reserve_exact_arithmetic(&hi, &lo, false)?;
    let width = hi - lo;
    dec.reserve_exact_arithmetic(&width, &Rational::from_i64(2), true)?;
    let requested = &width * &Rational::from_i64(2);
    let before_work = current_operation_work().used;
    let attempt = resize(&protected, 0, requested.clone(), dec);
    let edit_work = current_operation_work().used - before_work;
    let unchanged = state(&protected, dec)? == fingerprint;
    if !unchanged {
        return Err(problem(
            "import.material-history",
            "edit attempt changed the protected input",
        ));
    }
    let before_work = current_operation_work().used;
    let real = match attempt {
        Ok(after) => {
            if after.history().len() != protected.history().len() + 1 {
                return Err(problem(
                    "import.material-history",
                    "material edit lost input history",
                ));
            }
            let after_file = save(&after, &out.join("real-after.step"), dec)?;
            json!({"status":"numerical","exact":false,"before":before_file,"after":after_file,"history_retained":true})
        }
        Err(refusal) => {
            if refusal.code != "features.numerical-unsupported" {
                return Err(refusal);
            }
            let preserved = save(&protected, &out.join("real-after-refusal.step"), dec)?;
            if preserved["sha256"] != before_file["sha256"] {
                return Err(problem(
                    "import.material-history",
                    "refusal changed the prior STEP artifact",
                ));
            }
            json!({"status":"refused","exact":false,"refusal":{"code":refusal.code,"detail":refusal.detail},"before":before_file,"preserved_after_refusal":preserved,
                "prior_geometry_and_history_unchanged":true,"material_delta_applied":false})
        }
    };
    let after_step_work = current_operation_work().used - before_work;
    let before_work = current_operation_work().used;
    let positive = positive(out, dec)?;
    let positive_work = current_operation_work().used - before_work;
    Ok(
        json!({"status":"observed","proof":"numerical","exact":false,"model_id":"00000022","source_sha256":sha,"repaired_definition_sha256":derivative.sha256,
        "requested_edit":{"kind":"resize","axis":0,"length_mm":requested.to_string(),"expected_volume_ratio":"2","expected_delta_volume_ratio":"1"},
        "protected_history_stages":protected.history().len(),"real":real,"positive_control":positive,
        "stages":{"import":import_work,"pose":pose_work,"before_step":before_step_work,"material_edit":edit_work,"after_step":after_step_work,"positive_control":positive_work}}),
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: import_material_edit SOURCE.step NEW_OUTPUT_DIRECTORY".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let out = Path::new(&args[1]);
    // A new directory and create_new files preserve existing sources/artifacts.
    std::fs::create_dir(out)?;
    let (result, work) =
        with_operation_budget(20_000_000, || run(&bytes, out, &mut Decider::default()));
    let mut report = match result {
        Ok(v) => v,
        Err(e) => json!({"status":"failed","code":e.code,"detail":e.detail}),
    };
    if work.exhausted {
        report = json!({"status":"failed","code":"numbers.operation-budget"});
    }
    report["work_used"] = json!(work.used);
    report["work_limit"] = json!(20_000_000);
    std::fs::write(out.join("report.json"), format!("{report}\n"))?;
    println!("{report}");
    if report["status"] == "failed" {
        return Err("material edit qualification failed".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn supported_material_control_preserves_history_and_replays() {
        let out = std::env::temp_dir().join(format!("cx-material-positive-{}", std::process::id()));
        std::fs::create_dir(&out).unwrap();
        let (result, work) =
            with_operation_budget(20_000_000, || positive(&out, &mut Decider::default()));
        std::fs::remove_dir_all(&out).unwrap();
        let result = result.unwrap();
        assert!(!work.exhausted);
        assert_eq!(result["expected_delta_mm3"], "12");
        assert_eq!(
            result["after"]["history_stages"].as_u64().unwrap(),
            result["before"]["history_stages"].as_u64().unwrap() + 1
        );
    }

    #[test]
    fn unsupported_oblique_dimension_keeps_prior_geometry_and_history() {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut recipe = NumericalRecipe::parse(
                include_str!("../../../evidence/phase1/N04/recipe.json"),
                &mut dec,
            )?;
            recipe.edits.clear();
            let built = evaluate_numerical(&recipe, &mut dec)?;
            let pose = Pose::from_rational(
                Rotation::from_quaternion(4.into(), 0.into(), 3.into(), 0.into())?,
                RVector3::zero(),
            );
            let mut definition = built.body().definition().clone();
            for vertex in &mut definition.vertices {
                let [x, y, z] = vertex.position.clone().map(Real::from);
                let p = pose.apply_point(&Point3::new(x, y, z));
                vertex.position = [&p.x, &p.y, &p.z].map(|v| v.as_rational().unwrap());
            }
            for edge in &mut definition.edges {
                let Curve::Line(line) = &mut edge.curve else {
                    unreachable!()
                };
                line.origin = pose.apply_point(&line.origin);
                line.dir = pose.apply_vector(&line.dir);
            }
            for face in &mut definition.faces {
                let Surface::Plane(plane) = &mut face.surface else {
                    unreachable!()
                };
                plane.origin = pose.apply_point(&plane.origin);
                plane.normal = pose.apply_vector(&plane.normal);
                plane.ref_dir = pose.apply_vector(&plane.ref_dir);
            }
            let protected =
                built
                    .body()
                    .rebuild(definition, built.body().options().clone(), &mut dec)?;
            assert_eq!(protected.history().len(), built.body().history().len() + 1);
            let before = state(&protected, &mut dec)?;
            let refusal = resize(&protected, 0, 8.into(), &mut dec).unwrap_err();
            assert_eq!(refusal.code, "features.numerical-unsupported");
            assert_eq!(before, state(&protected, &mut dec)?);
            Ok::<_, Refusal>(())
        });
        result.unwrap();
        assert!(!work.exhausted);
    }
}
