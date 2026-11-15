//! Actual source -> numerical B-rep -> pose rebuild -> approximate disk STEP.

use std::path::Path;

use catallax_import::{
    numerical::{self, ImportError, StepImport},
    sha256_hex,
};
use catallax_numbers::{
    Decider, NumericalCheck, Rational, Real, Refusal, current_operation_work, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Point3, Surface};
use catallax_step::numerical::{NumericalStepOptions, build as export};
use catallax_topology::numerical as n;
use serde_json::{Value, json};

fn error(stage: &str, error: Refusal) -> Value {
    json!({"stage":stage,"code":error.code,"detail":error.detail})
}
fn evidence(body: &n::NumericalBody) -> Value {
    let checks = body.evidence().checks().iter().map(|check| match check {
        NumericalCheck::Passed { name, scope } => json!({"kind":"passed","name":name,"scope":scope}),
        NumericalCheck::Samples { name, scope, samples } => json!({"kind":"samples","name":name,"scope":scope,
            "samples":samples.iter().map(|sample|json!({"parameters":sample.parameters.iter().map(ToString::to_string).collect::<Vec<_>>(),"distance_upper_mm":sample.distance_upper_mm.to_string()})).collect::<Vec<_>>() }),
    }).collect::<Vec<_>>();
    json!({"proof":"numerical","exact":false,"tolerance_mm":body.options().tolerance_mm.to_string(),"checks":checks,
        "counts":{"faces":body.counts().faces,"edges":body.counts().edges,"vertices":body.counts().vertices}})
}
fn shifted(p: &Point3, shift: &[Real; 3]) -> Point3 {
    Point3::new(&p.x + &shift[0], &p.y + &shift[1], &p.z + &shift[2])
}
fn pose_edit(body: &n::NumericalBody, dec: &mut Decider) -> Result<n::NumericalBody, Refusal> {
    let mut definition = body.definition().clone();
    let shift = [3, 5, 7].map(Rational::from_i64);
    let delta = shift.clone().map(Real::from);
    for vertex in &mut definition.vertices {
        for i in 0..3 {
            dec.charge_work(1)?;
            vertex.position[i] = &vertex.position[i] + &shift[i];
        }
    }
    for edge in &mut definition.edges {
        dec.charge_work(1)?;
        match &mut edge.curve {
            Curve::Line(line) => line.origin = shifted(&line.origin, &delta),
            Curve::Circle(circle) => circle.center = shifted(&circle.center, &delta),
            _ => {
                return Err(Refusal::new(
                    "import.example-pose",
                    "pose control requires analytic lines or circles",
                ));
            }
        }
    }
    for face in &mut definition.faces {
        dec.charge_work(1)?;
        match &mut face.surface {
            Surface::Plane(plane) => plane.origin = shifted(&plane.origin, &delta),
            Surface::Cylinder(cylinder) => cylinder.origin = shifted(&cylinder.origin, &delta),
            _ => {
                return Err(Refusal::new(
                    "import.example-pose",
                    "pose control requires planes or cylinders",
                ));
            }
        }
    }
    body.rebuild(definition, body.options().clone(), dec)
}
fn file(body: &n::NumericalBody, path: &Path, dec: &mut Decider) -> Result<Value, Value> {
    let file = export(body, &NumericalStepOptions::default(), dec).map_err(|e| error("step", e))?;
    std::fs::write(path, file.bytes())
        .map_err(|e| json!({"stage":"write","detail":e.to_string()}))?;
    let bytes =
        std::fs::read(path).map_err(|e| json!({"stage":"reread","detail":e.to_string()}))?;
    file.check_written(&bytes, dec)
        .map_err(|e| error("reread-check", e))?;
    let phase_bounds = file
        .report()
        .emitted_evidence
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, samples, .. }
                if name == "file.closed-circle-zero-phase-tangent-bound" =>
            {
                Some(samples)
            }
            _ => None,
        })
        .flat_map(|samples| samples.iter().map(|sample| &sample.distance_upper_mm))
        .collect::<Vec<_>>();
    Ok(
        json!({"sha256":sha256_hex(&bytes),"bytes":bytes.len(),"category":"numerical approximate",
        "checks":file.report().checks,"scalar_rounding_limit_mm":file.report().scalar_rounding_limit_mm.to_string(),
        "max_scalar_rounding_error_mm":file.report().max_scalar_rounding_error_mm.to_string(),
        "max_direction_ratio_rounding_error":file.report().max_direction_ratio_rounding_error.to_string(),
        "closed_circle_phase_checks":{"count":phase_bounds.len(),"max_tangent_bound_mm":phase_bounds.iter().max().map(ToString::to_string),"angular_limit_rad":"1/1000000000000","distance_limit_mm":"1/1000000000"},
        "source_history_stages":file.report().source_history.len(),"complete_body_bound":"unestablished"}),
    )
}
fn run(bytes: &[u8], expected_sha256: &str, out: &Path) -> Result<Value, Value> {
    let mut dec = Decider::default();
    let imported = numerical::import_step(
        bytes,
        expected_sha256,
        &StepId::new("source").unwrap(),
        &n::NumericalOptions::default(),
        &mut dec,
    )
    .map_err(|e| match e {
        ImportError::Source(e) => json!({"stage":"import","source_refusal":e}),
        ImportError::Numerical(e) => error("import", e),
    })?;
    let import_work = current_operation_work().used;
    let StepImport::Numerical(imported) = imported else {
        return Ok(json!({"status":"certified","work":{"import":import_work}}));
    };
    let admission = imported.source_admission().map(|record| {
        json!({
        "method":record.method,"tolerance_mm":record.tolerance_mm.to_string(),
        "strict_refusal":record.strict_refusal,"intersection_checks":record.intersection_checks,
        "topology_repairs":record.topology_repairs,
        "optional_exact_work_reserved":record.optional_exact_work_reserved,
        "optional_exact_work_used":record.optional_exact_work_used,
        "cross_solid_pairs":"not_compared_by_import_policy","shared_3d_geometry":"retained",
        "parameter_trims":"rebuilt","complete_body_bound":"unestablished"})
    });
    let mut solids = Vec::new();
    for (index, (key, body)) in imported.bodies().iter().enumerate() {
        let before = current_operation_work().used;
        let original = file(body, &out.join(format!("solid-{index}.step")), &mut dec)
            .map_err(|refusal| json!({"stage":"source-step","cause":refusal,
                "source_body":evidence(body),"source_record":imported.record(),"source_admission":admission,"import_work":import_work}))?;
        let export_work = current_operation_work().used - before;
        let before = current_operation_work().used;
        let edited = pose_edit(body, &mut dec).map_err(|e| error("pose-edit", e))?;
        let edit_work = current_operation_work().used - before;
        let before = current_operation_work().used;
        let edited_file = file(
            &edited,
            &out.join(format!("solid-{index}-edited.step")),
            &mut dec,
        )?;
        let edited_export_work = current_operation_work().used - before;
        solids.push(json!({"solid":key,"source_body":evidence(body),"source_step":original,
            "pose_translation_mm":["3","5","7"],"edited_body":evidence(&edited),"edited_step":edited_file,
            "work":{"source_export":export_work,"pose_rebuild":edit_work,"edited_export":edited_export_work}}));
    }
    Ok(
        json!({"status":"numerical","proof":"numerical","exact":false,"source_record":imported.record(),"source_admission":admission,
        "tolerance_mm":"0.01","solids":solids,"work":{"import":import_work}}),
    )
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err("expected source STEP path and output directory".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let out = Path::new(&args[1]);
    std::fs::create_dir_all(out)?;
    let source_sha256 = sha256_hex(&bytes);
    let (result, work) =
        with_operation_budget(n::DEFAULT_WORK, || run(&bytes, &source_sha256, out));
    let mut response = match result {
        Ok(v) => v,
        Err(e) => json!({"status":"refused","refusal":e}),
    };
    if work.exhausted {
        response = json!({"status":"refused","refusal":{"code":"resource_limit"}});
    }
    response["source_sha256"] = json!(source_sha256);
    response["work_used"] = json!(work.used);
    response["work_limit"] = json!(n::DEFAULT_WORK);
    std::fs::write(out.join("response.json"), format!("{response}\n"))?;
    println!("{response}");
    Ok(())
}
