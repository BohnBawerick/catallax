//! Explicit massless-envelope STEP outputs backed by field moment enclosures.
use catallax_field::{Plan, mass};
use catallax_measure::{Accuracy as QuantityAccuracy, Enclosure};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Number, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{fs, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct QuantityRequest {
    relative: Number,
    absolute: Number,
}
impl QuantityRequest {
    fn accuracy(&self) -> QuantityAccuracy {
        QuantityAccuracy {
            relative: self.relative.value().clone(),
            absolute: self.absolute.value().clone(),
        }
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AccuracyRequest {
    volume: QuantityRequest,
    centroid: QuantityRequest,
    inertia: QuantityRequest,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    format: String,
    recipe: Value,
    accuracy: AccuracyRequest,
    #[serde(default)]
    max_cells: Option<usize>,
}
fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/3","status":"refused","exact":false,"result":null,"properties":null,"step":null,"marks":[],"refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}
fn bound(value: &Enclosure) -> Value {
    json!({"lo":value.lo.to_string(),"hi":value.hi.to_string()})
}
fn scalars(values: &[Rational; 3]) -> [String; 3] {
    values.each_ref().map(ToString::to_string)
}
fn accuracy(value: &QuantityAccuracy) -> Value {
    json!({"relative":value.relative.to_string(),"absolute":value.absolute.to_string()})
}

pub(crate) fn evaluate(problem: &Value, out: &Path) -> std::io::Result<Value> {
    let request: Request = match serde_json::from_value(problem.clone()) {
        Ok(request) => request,
        Err(error) => return Ok(refused("cli.invalid-stand-in", error)),
    };
    if request.format != "catallax-stand-in/1" {
        return Ok(refused(
            "cli.invalid-stand-in",
            "expected catallax-stand-in/1",
        ));
    }
    let recipe = match Recipe::parse(&request.recipe.to_string()) {
        Ok(recipe) => recipe,
        Err(error) => return Ok(refused("cli.invalid-recipe", error)),
    };
    let requested = mass::Accuracy {
        volume: request.accuracy.volume.accuracy(),
        centroid: request.accuracy.centroid.accuracy(),
        inertia: request.accuracy.inertia.accuracy(),
    };
    let mut dec = Decider::default();
    let mut built = match catallax_features::evaluate(&recipe, &mut dec) {
        Ok(built) => built,
        Err(catallax_features::BuildError::Geometry { step, refusal }) => {
            let mut response = refused(refusal.code, refusal.detail);
            response["refusal"]["step"] = json!(step.to_string());
            return Ok(response);
        }
        Err(error) => return Ok(refused("cli.invalid-recipe", error)),
    };
    let plan = match Plan::compile_checked(&recipe, &mut dec) {
        Ok(plan) => plan,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    if let Err(error) = mass::check_region_envelope(&plan) {
        return Ok(refused(error.code, error.detail));
    }
    let span =
        catallax_numbers::OperationSpan::begin("field", &recipe.result, "mass-properties", &dec);
    let measured = mass::properties(
        &plan,
        &requested,
        &mass::Options {
            max_cells: request
                .max_cells
                .unwrap_or(mass::Options::default().max_cells),
        },
        &mut dec,
    );
    span.finish(&dec, measured.as_ref().err(), None);
    let measured = match measured {
        Ok(properties) => properties,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let field_properties = json!({"format":"catallax-field-mass/1","scope":"original-field-definition","source":"field-interval-integration",
        "units":{"volume":"mm3","centroid":"mm","inertia":"mm5"},"volume":bound(&measured.volume),
        "centroid":measured.centroid.as_ref().map(|c|c.each_ref().map(bound)),
        "inertia":measured.inertia.as_ref().map(|m|m.each_ref().map(|row|row.each_ref().map(bound))),"cells":measured.cells});
    let Some(target) = measured.estimate.as_ref() else {
        return Ok(
            json!({"format":"catallax-manifest/3","status":"answered","exact":true,"result":null,"properties":null,"step":null,"marks":[],"refusal":null,
            "stand_in":{"representation":"empty","field_properties":field_properties,"mass_bodies":0}}),
        );
    };
    let span =
        catallax_numbers::OperationSpan::begin("brep", &recipe.result, "mass-surrogate", &dec);
    let surrogate = catallax_features::stand_in::equivalent_blocks(&mut built.model, target);
    span.finish(&dec, surrogate.as_ref().err(), None);
    let surrogate = match surrogate {
        Ok(surrogate) => surrogate,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let inertia_limit = Rational::from_ratio(1, 1_000_000_000).expect("nonzero denominator");
    if surrogate.inertia_residual >= inertia_limit {
        return Ok(refused(
            "stand-in.inertia-not-matched",
            format!(
                "closed-form surrogate inertia residual {} exceeds the strict 1e-9 mm5 comparison",
                surrogate.inertia_residual
            ),
        ));
    }
    if !measured.contains(&surrogate.achieved) {
        return Ok(refused(
            "stand-in.properties-not-matched",
            format!(
                "surrogate properties exceed the field enclosures; exact inertia residual {}",
                surrogate.inertia_residual
            ),
        ));
    }
    let options = catallax_step::StepOptions {
        product_name: format!("{} field stand-in", recipe.name),
        ..Default::default()
    };
    let span =
        catallax_numbers::OperationSpan::begin("brep", &recipe.result, "stand-in-step", &dec);
    let exported = catallax_step::write_stand_in_step(
        &built.model,
        built.result,
        &surrogate.bodies,
        &options,
        &mut dec,
        &mut fs::File::create(out.join("answer.step"))?,
    );
    span.finish(
        &dec,
        exported.as_ref().err(),
        exported.as_ref().ok().map(|a| &a.marks),
    );
    let exported = match exported {
        Ok(exported) => exported,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let exchanged = match mass::surrogate_properties(&surrogate, &exported.value.max_deviation) {
        Ok(properties) => properties,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    if !measured.encloses(&exchanged) {
        return Ok(refused(
            "stand-in.exchange-not-enclosed",
            "STEP-rounded surrogate property bounds exceed the original field enclosures",
        ));
    }
    let marks: Vec<_> = exported
        .marks
        .iter()
        .map(|mark| crate::mark_entry(mark, None))
        .collect();
    Ok(
        json!({"format":"catallax-manifest/3","status":"answered","exact":false,"result":null,"properties":null,"marks":marks,"refusal":null,
        "step":{"file":"answer.step","max_deviation":exported.value.max_deviation.to_string(),"solids":exported.value.solids,
            "faces":exported.value.faces,"surface_envelopes":1},
        "stand_in":{"representation":"massless-surface-envelope-plus-unit-density-solids","scope":"field-mass-equivalence-not-boundary-conversion",
            "density_metadata":false,"material_metadata":false,"within_field_enclosures":true,"mass_bodies":surrogate.bodies.len(),
            "field_accuracy":{"volume":accuracy(&requested.volume),"centroid":accuracy(&requested.centroid),"inertia":accuracy(&requested.inertia)},
            "field_properties":field_properties,
            "surrogate_properties":{"volume":surrogate.achieved.volume.to_string(),"centroid":scalars(&surrogate.achieved.centroid),
                "inertia":surrogate.achieved.inertia.each_ref().map(scalars)},
            "exported_property_bounds":{"volume":bound(&exchanged.volume),
                "centroid":exchanged.centroid.each_ref().map(bound),
                "inertia":exchanged.inertia.each_ref().map(|row|row.each_ref().map(bound))},
            "inertia_residual":surrogate.inertia_residual.to_string()}}),
    )
}
