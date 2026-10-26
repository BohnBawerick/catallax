//! Explicit IGES output for the deterministic verification interface.
use catallax_iges::{IgesError, IgesOptions, IgesReport, Units, write_iges};
use catallax_numbers::{Decider, Refusal};
use catallax_recipe::{Number, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{fs, io, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    format: String,
    recipe: Value,
    #[serde(default)]
    tolerance: Option<Number>,
    #[serde(default)]
    scale_digits: Option<u32>,
    #[serde(default)]
    units: Option<String>,
}

fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/3","status":"refused","exact":false,"result":null,
        "step":null,"marks":[],"refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}
fn invalid(code: &str, detail: impl ToString) -> Value {
    let mut value = refused(code, detail);
    value["status"] = json!("invalid-question");
    value
}

/// The stable refusal for a typed writer failure. An I/O failure is not here: it is
/// returned to the caller as an infrastructure error.
fn writer_refusal(error: &IgesError) -> Refusal {
    match error {
        IgesError::InvalidOptions => Refusal::new(
            "iges.invalid-options",
            "tolerance must not be negative and scale_digits must be 1 to 100",
        ),
        IgesError::InvalidTopology(detail) => Refusal::new("iges.invalid-topology", detail.clone()),
        IgesError::UncertifiedModel => Refusal::new(
            "iges.uncertified-model",
            "the model lacks the exact-definition correspondence the writer needs",
        ),
        IgesError::UnsupportedSurface(kind) => Refusal::new(
            "iges.unsupported-surface",
            format!("surface kind {kind:?} is not written"),
        ),
        IgesError::UnsupportedCurve(kind) => Refusal::new(
            "iges.unsupported-curve",
            format!("curve kind {kind:?} is not written"),
        ),
        IgesError::UnrepresentableDefinition(detail) => {
            Refusal::new("iges.unrepresentable-definition", *detail)
        }
        IgesError::ToleranceExceeded { proven, requested } => Refusal::new(
            "iges.tolerance-exceeded",
            format!("proven bound {proven} mm exceeds the requested {requested} mm"),
        ),
        IgesError::Decision(refusal) => refusal.clone(),
        IgesError::Io(error) => Refusal::new("iges.io", error.to_string()),
    }
}

pub(crate) fn evaluate(problem: &Value, out: &Path) -> io::Result<Value> {
    evaluate_with(problem, out, &mut |path, bytes| fs::write(path, bytes))
}

/// `store` is the artifact sink, so a failing one can be exercised without a
/// platform-specific full device.
fn evaluate_with(
    problem: &Value,
    out: &Path,
    store: &mut dyn FnMut(&Path, &[u8]) -> io::Result<()>,
) -> io::Result<Value> {
    let request: Request = match serde_json::from_value(problem.clone()) {
        Ok(request) => request,
        Err(error) => return Ok(invalid("cli.invalid-iges", error)),
    };
    if request.format != "catallax-iges/1" {
        return Ok(invalid("cli.invalid-iges", "expected catallax-iges/1"));
    }
    if request.units.as_deref().is_some_and(|units| units != "mm") {
        return Ok(invalid("cli.invalid-iges", "units must be \"mm\""));
    }
    let defaults = IgesOptions::default();
    let options = IgesOptions {
        tolerance: request
            .tolerance
            .map_or(defaults.tolerance, |n| n.value().clone()),
        scale_digits: request.scale_digits.unwrap_or(defaults.scale_digits),
        units: Units::Millimetres,
    };
    let recipe = match Recipe::parse(&request.recipe.to_string()) {
        Ok(recipe) => recipe,
        Err(error) => return Ok(invalid("cli.invalid-recipe", error)),
    };
    let mut dec = Decider::default();
    let built = match catallax_features::evaluate(&recipe, &mut dec) {
        Ok(built) => built,
        Err(catallax_features::BuildError::Geometry { step, refusal }) => {
            let mut value = refused(refusal.code, refusal.detail);
            value["refusal"]["step"] = json!(step.to_string());
            return Ok(value);
        }
        Err(error) => return Ok(invalid("cli.invalid-recipe", error)),
    };
    if !built.stand_ins.is_empty() {
        return Ok(refused(
            "field.stand-in-export",
            "the lattice envelope alone is not a mass-equivalent stand-in",
        ));
    }
    let marks: Vec<_> = built
        .model
        .body(built.result)
        .marks
        .iter()
        .map(|mark| crate::mark_entry(mark, None))
        .collect();
    let span = catallax_numbers::OperationSpan::begin("brep", &recipe.result, "iges", &dec);
    let mut bytes = Vec::new();
    let written = write_iges(
        &built.model,
        &[built.result],
        &options,
        &mut dec,
        &mut bytes,
    );
    let refusal = written.as_ref().err().map(writer_refusal);
    span.finish(&dec, refusal.as_ref(), None);
    let report = match written {
        Ok(report) => report,
        Err(IgesError::Io(error)) => return Err(error),
        Err(_) => {
            let refusal = refusal.expect("a typed failure has a refusal");
            return Ok(refused(refusal.code, refusal.detail));
        }
    };
    let path = out.join("answer.igs");
    if let Err(error) = store(&path, &bytes) {
        // A prefix must not survive as an artifact of an unanswered request.
        let _ = fs::remove_file(&path);
        return Err(error);
    }
    Ok(json!({"format":"catallax-manifest/3","status":"answered",
        "exact":marks.is_empty() && report.max_deviation == catallax_numbers::Rational::zero(),
        "result":null,"properties":null,"step":null,"marks":marks,"refusal":null,
        "iges":published(&options, &report, &bytes)}))
}

fn published(options: &IgesOptions, report: &IgesReport, bytes: &[u8]) -> Value {
    json!({"file":"answer.igs","bytes":bytes.len(),
        "sha256":Sha256::digest(bytes).iter().map(|b| format!("{b:02x}")).collect::<String>(),
        "units":"mm","tolerance_mm":options.tolerance.to_string(),"scale_digits":options.scale_digits,
        "max_deviation_mm":report.max_deviation.to_string(),
        "trim_deviation_mm":report.trim_deviation.to_string(),
        "scale":report.scale.to_string(),
        "coordinate_bound_mm":report.coordinate_bound.to_string(),
        "solids":report.solids,"faces":report.faces,"edges":report.edges,
        "added_edges":report.added_edges,"vertices":report.vertices,
        "added_vertices":report.added_vertices,
        "coverage":"whole-domain file-to-definition boundary displacement, all vertices, curves, faces and trims"})
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_failing_artifact_sink_is_an_io_error_and_leaves_no_file() {
        let recipe: Value = serde_json::from_str(include_str!(
            "../../../.claude/skills/verify-catallax/references/examples/recipe.json"
        ))
        .unwrap();
        let dir = std::env::temp_dir().join(format!("catallax-iges-sink-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let problem = json!({"format":"catallax-iges/1","recipe":recipe});
        let error = evaluate_with(&problem, &dir, &mut |path, bytes| {
            fs::write(path, &bytes[..bytes.len() / 2])?;
            Err(io::Error::other("device full"))
        })
        .unwrap_err();
        assert_eq!(error.to_string(), "device full");
        assert!(!dir.join("answer.igs").exists());
    }
}
