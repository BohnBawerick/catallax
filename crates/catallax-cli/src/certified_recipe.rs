//! Opt-in certified geometry. Existing request adapters retain their contracts.

use std::{fs, path::Path};

use catallax_app::CertifiedApplicationBody;
use catallax_field::{CertifiedPrismField, Direction, Location};
use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign};
use catallax_recipe::{Number, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};

use crate::recipe_edit::ScalarEdit;
mod certificate;

pub(crate) const FORMAT: &str = "catallax-certified-recipe/1";
pub(crate) fn is_problem(problem: &Value) -> bool {
    problem["format"] == FORMAT
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    format: String,
    recipe: Value,
    outputs: Outputs,
    properties: Option<PropertyOptions>,
    #[serde(default)]
    field_points: Vec<[Number; 3]>,
    edit: Option<ScalarEdit>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Outputs {
    mesh: Option<MeshOptions>,
    step: Option<StepOptions>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MeshOptions {
    chord_tolerance: Number,
    max_vertices: usize,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StepOptions {
    file_tolerance: Number,
    decimal_places: u32,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PropertyOptions {
    relative: Number,
    absolute: Number,
    #[serde(default)]
    inertia: bool,
}

pub(crate) fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/4","status":"refused","exact":false,
        "result":null,"certificate":null,"display":null,"properties":null,"field":null,"step":null,"edit":null,"marks":[],
        "refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}
fn consumer_refusal(error: &Refusal) -> Value {
    json!({"status":"refused","refusal":{"code":error.code,"detail":error.detail}})
}
fn invalid(detail: impl ToString) -> Value {
    let mut result = refused("cli.invalid-certified-recipe", detail);
    result["status"] = json!("invalid-question");
    result
}
fn app_refusal(error: &catallax_app::Error) -> Value {
    let mut result = refused(error.code(), error);
    if let catallax_app::Error::Build(catallax_features::BuildError::Geometry { step, refusal }) =
        error
    {
        result["refusal"]["detail"] = json!(refusal.detail);
        result["refusal"]["step"] = json!(step.to_string());
    }
    result
}
fn composed(app: &CertifiedApplicationBody, output_error: &Rational) -> Option<String> {
    app.built()
        .result_record()
        .uniform_bound
        .as_ref()
        .map(|model| (&model.upper + output_error).to_string())
}
fn names(app: &CertifiedApplicationBody) -> Value {
    json!({"faces":app.names().faces.iter().map(|(face, origins)| json!({"face":face.index(),"origins":origins})).collect::<Vec<_>>(),
        "edges":app.names().edges.iter().map(|(edge, name)| json!({"edge":edge.index(),"name":name})).collect::<Vec<_>>(),
        "vertices":app.names().vertices.iter().map(|(vertex, name)| json!({"vertex":vertex.index(),"name":name})).collect::<Vec<_>>()})
}
fn display(
    app: &CertifiedApplicationBody,
    options: Option<&MeshOptions>,
    out: &Path,
) -> std::io::Result<Value> {
    let Some(display) = app.display() else {
        return Ok(Value::Null);
    };
    let display = match display {
        Ok(display) => display,
        Err(error) => return Ok(consumer_refusal(error)),
    };
    let mesh = &display.value;
    mesh.write_obj(&mut fs::File::create(out.join("display.obj"))?)?;
    Ok(
        json!({"status":"answered","file":"display.obj","scope":"represented-model","vertices":mesh.positions.len(),"triangles":mesh.triangles.len(),
        "requested_chord_tolerance_mm":options.unwrap().chord_tolerance.value().to_string(),
        "faces":mesh.faces.iter().map(|face| json!({"face":face.face.index(),"group":format!("face_{}",face.face.index()),
            "triangle_range":[face.triangles.start,face.triangles.end],"model_chord_bound_mm":face.max_chord_error.to_string(),
            "definition_bound_mm":composed(app,&face.max_chord_error)})).collect::<Vec<_>>(),
        "edges":mesh.edges.iter().map(|edge| json!({"edge":edge.edge.index(),"vertices":edge.vertices,
            "model_chord_bound_mm":edge.max_chord_error.to_string(),"definition_bound_mm":composed(app,&edge.max_chord_error)})).collect::<Vec<_>>() }),
    )
}
fn properties(
    app: &CertifiedApplicationBody,
    options: Option<&PropertyOptions>,
    dec: &mut Decider,
) -> Value {
    let Some(options) = options else {
        return Value::Null;
    };
    let accuracy = catallax_measure::Accuracy {
        relative: options.relative.value().clone(),
        absolute: options.absolute.value().clone(),
    };
    let (volume, area, center, inertia, measured_marks) = if options.inertia {
        match app.mass_properties(&accuracy, dec) {
            Ok(answer) => (
                answer.value.volume,
                answer.value.area,
                answer.value.centroid,
                Some(answer.value.inertia),
                answer.marks,
            ),
            Err(error) => return consumer_refusal(&error),
        }
    } else {
        match app.volume_area_centroid(&accuracy, dec) {
            Ok(answer) => (
                answer.value.volume,
                answer.value.area,
                answer.value.centroid,
                None,
                answer.marks,
            ),
            Err(error) => return consumer_refusal(&error),
        }
    };
    let values = crate::Properties {
        volume: Some(crate::bound(&volume, &accuracy)),
        area: Some(crate::bound(&area, &accuracy)),
        centroid: Some(center.as_ref().map(|center| {
            center
                .each_ref()
                .map(|value| crate::bound(value, &accuracy))
        })),
        inertia: inertia.as_ref().map(|tensor| {
            tensor.as_ref().map(|tensor| {
                tensor
                    .each_ref()
                    .map(|row| row.each_ref().map(|value| crate::bound(value, &accuracy)))
            })
        }),
    };
    let mut marks: Vec<_> = measured_marks
        .iter()
        .map(|mark| crate::mark_entry(mark, None))
        .collect();
    crate::mark_decimal_accuracy(&values, &accuracy, &mut marks);
    if !marks.is_empty() {
        return json!({"status":"refused","scope":"original-definition","marks":marks,
            "refusal":{"code":"measure.accuracy-not-reached","detail":"the original-target enclosure did not reach the requested width"}});
    }
    json!({"status":"answered","scope":"original-definition","relative":options.relative,"absolute":options.absolute,
        "units":{"volume":"mm3","area":"mm2","centroid":"mm","inertia":"mm5"},"values":values})
}
fn scalar(value: f64) -> String {
    if value == f64::INFINITY {
        "infinity".into()
    } else if value == f64::NEG_INFINITY {
        "-infinity".into()
    } else {
        Rational::from_f64_exact(value)
            .expect("finite proven interval endpoint")
            .to_string()
    }
}
fn interval(value: Interval) -> Value {
    json!({"lo":scalar(value.lo()),"hi":scalar(value.hi())})
}
fn field(
    app: &CertifiedApplicationBody,
    points: &[[Number; 3]],
    dec: &mut Decider,
) -> Result<Value, Refusal> {
    if points.is_empty() {
        return Ok(Value::Null);
    }
    let plan = CertifiedPrismField::new(app.built(), dec)?;
    let mut results = Vec::with_capacity(points.len());
    for point in points {
        let exact = catallax_shapes::RPoint3::new(
            point[0].value().clone(),
            point[1].value().clone(),
            point[2].value().clone(),
        );
        let query = plan.proven_checked(&exact, dec)?;
        let location = match query.location {
            Location::Inside => "inside",
            Location::Outside => "outside",
            Location::Boundary => "boundary",
            Location::Indeterminate => "indeterminate",
            _ => {
                return Err(Refusal::new(
                    catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                    "unknown field classification",
                ));
            }
        };
        let direction = match query.direction {
            Direction::Unique(normal) => json!({"kind":"unique","components":normal.map(interval)}),
            Direction::Ambiguous => json!({"kind":"ambiguous"}),
            Direction::Undefined => json!({"kind":"undefined"}),
            _ => {
                return Err(Refusal::new(
                    catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                    "unknown field direction",
                ));
            }
        };
        let catallax_field::Value::TrueDistance(value) = query.value else {
            return Err(Refusal::new(
                catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                "nominal field did not provide true distance",
            ));
        };
        results.push(json!({"point":point,"location":location,"value":{"kind":"true-distance","mm":interval(value)},
            "unsigned_distance_mm":interval(query.distance.millimetres),"direction":direction}));
    }
    Ok(json!({"status":"answered","scope":"original-definition","points":results}))
}
fn step(
    app: &CertifiedApplicationBody,
    options: Option<&StepOptions>,
    out: &Path,
    dec: &mut Decider,
) -> std::io::Result<Value> {
    let Some(options) = options else {
        return Ok(Value::Null);
    };
    let certified = app.built();
    let built = certified.built();
    if let Some(junction) = certified.junction() {
        return junction_step(app, junction, options, out, dec);
    }
    let certificate = certified.result_certificate();
    let mut bytes = Vec::new();
    let report = if certificate.target().is_empty() {
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &catallax_step::StepOptions {
                product_name: built.recipe().name.clone(),
                curve_tolerance: options.file_tolerance.value().clone(),
                significant_digits: 17,
            },
            dec,
            &mut bytes,
        )
    } else if certificate.spline().is_some() {
        catallax_step::write_spline_prism_step(
            &built.model,
            built.result,
            &catallax_step::CertifiedSplineStepOptions {
                product_name: built.recipe().name.clone(),
                file_tolerance: options.file_tolerance.value().clone(),
                decimal_places: options.decimal_places,
            },
            dec,
            &mut bytes,
        )
    } else {
        return Ok(consumer_refusal(&Refusal::new(
            catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
            "certified spline export requires a genuine spline body; analytic export retains its existing adapter",
        )));
    };
    let report = match report {
        Ok(answer) => answer.value,
        Err(error) => return Ok(consumer_refusal(&error)),
    };
    fs::write(out.join("answer.step"), bytes)?;
    let regions: Vec<_> = certificate.faces().into_iter().flatten().map(|face| {
        let region = certificate.region_for_face(face).expect("complete geometry-bound region");
        json!({"model_region":certificate::key(&region.identity.representation),"face":face.index(),"file_region":format!("face_{}",face.index()),
            "file_to_model_bound_mm":report.max_deviation.to_string(),"file_to_definition_bound_mm":composed(app,&report.max_deviation)})
    }).collect();
    Ok(
        json!({"status":"answered","file":"answer.step","scope":"represented-model","requested_file_to_model_mm":options.file_tolerance,
        "decimal_places":options.decimal_places,"file_to_model_bound_mm":report.max_deviation.to_string(),
        "file_to_definition_bound_mm":composed(app,&report.max_deviation),"regions":regions,
        "added_vertices":report.added_vertices,"added_edges":report.added_edges}),
    )
}
// Every region composes its own bound with the file's. Nothing reaches the
// answer directory until the whole file is proved.
fn junction_step(
    app: &CertifiedApplicationBody,
    junction: &catallax_features::JunctionCertificate,
    options: &StepOptions,
    out: &Path,
    dec: &mut Decider,
) -> std::io::Result<Value> {
    let built = app.built().built();
    let mut bytes = Vec::new();
    let report = match catallax_step::write_rolling_junction_step(
        &built.model,
        built.result,
        &catallax_step::CertifiedSplineStepOptions {
            product_name: built.recipe().name.clone(),
            file_tolerance: options.file_tolerance.value().clone(),
            decimal_places: options.decimal_places,
        },
        dec,
        &mut bytes,
    ) {
        Ok(answer) => answer.value,
        Err(error) => return Ok(consumer_refusal(&error)),
    };
    fs::write(out.join("answer.step"), bytes)?;
    let regions: Vec<_> = junction.faces().iter().map(|&face| {
        let region = junction.region_for_face(face).expect("complete geometry-bound region");
        json!({"model_region":certificate::key(&region.identity.representation),"face":face.index(),"file_region":format!("face_{}",face.index()),
            "file_to_model_bound_mm":report.max_deviation.to_string(),
            "file_to_definition_bound_mm":(junction.face_bound(face) + report.max_deviation.clone()).to_string()})
    }).collect();
    Ok(
        json!({"status":"answered","file":"answer.step","scope":"represented-model","requested_file_to_model_mm":options.file_tolerance,
        "decimal_places":options.decimal_places,"file_to_model_bound_mm":report.max_deviation.to_string(),
        "file_to_definition_bound_mm":composed(app,&report.max_deviation),"regions":regions,
        "added_vertices":report.added_vertices,"added_edges":report.added_edges}),
    )
}
fn edit(
    app: &CertifiedApplicationBody,
    recipe: &Recipe,
    edit: Option<&ScalarEdit>,
    dec: &mut Decider,
) -> Result<Value, Refusal> {
    let Some(edit) = edit else {
        return Ok(Value::Null);
    };
    let recipe = edit
        .apply(recipe)
        .map_err(|error| Refusal::new("cli.certified-edit", error))?;
    let edited = catallax_app::build_certified(&recipe, None, dec)
        .map_err(|error| Refusal::new(error.code(), error.to_string()))?;
    let built = edited.built().built();
    let summary = crate::result_summary(
        &built.model,
        built.result,
        dec,
        crate::SummaryExchange::Certified,
    )?;
    let mut resolved = Vec::new();
    for (face, origins) in &app.names().faces {
        let matches: Vec<_> = edited
            .names()
            .faces
            .iter()
            .filter(|(_, candidate)| origins.iter().any(|origin| candidate.contains(origin)))
            .map(|(face, origins)| json!({"face":face.index(),"origins":origins}))
            .collect();
        resolved.push(json!({"kind":"face","old_handle":face.index(),"reference":origins,"status":if matches.is_empty(){"missing"}else{"resolved"},"matches":matches}));
    }
    for (edge, name) in &app.names().edges {
        if let Some(name) = name {
            let matches: Vec<_> = edited
                .names()
                .edges
                .iter()
                .filter(|(_, candidate)| candidate.as_ref() == Some(name))
                .map(|(edge, _)| edge.index())
                .collect();
            resolved.push(json!({"kind":"edge","old_handle":edge.index(),"reference":name,"status":if matches.is_empty(){"missing"}else{"resolved"},"matches":matches}));
        }
    }
    Ok(
        json!({"status":"answered","step":edit.step,"parameter":edit.parameter,"value":edit.value,
        "result":summary,"certificate":certificate::record(edited.built()),"names":names(&edited),"resolved":resolved}),
    )
}

pub(crate) fn evaluate(problem: &Value, out: &Path) -> std::io::Result<Value> {
    let request: Request = match Request::deserialize(problem) {
        Ok(request) => request,
        Err(error) => return Ok(invalid(error)),
    };
    if request.format != FORMAT
        || !matches!(
            request.recipe["format"].as_str(),
            Some("catallax-recipe/3" | "catallax-recipe/4")
        )
    {
        return Ok(invalid(
            "expected certified-recipe/1 containing unchanged recipe/3 or recipe/4",
        ));
    }
    if request.properties.as_ref().is_some_and(|p| {
        p.relative.value().sign() == Sign::Negative || p.absolute.value().sign() == Sign::Negative
    }) {
        return Ok(invalid("property accuracy must be nonnegative"));
    }
    if request
        .outputs
        .step
        .as_ref()
        .is_some_and(|options| options.file_tolerance.value().sign() == Sign::Negative)
    {
        return Ok(invalid("file tolerance must be nonnegative"));
    }
    let recipe = match Recipe::parse(&request.recipe.to_string()) {
        Ok(recipe) => recipe,
        Err(error) => return Ok(invalid(error)),
    };
    let mesh_options = request
        .outputs
        .mesh
        .as_ref()
        .map(|options| catallax_mesh::MeshOptions {
            chord_tolerance: options.chord_tolerance.value().clone(),
            max_vertices: options.max_vertices,
            ..Default::default()
        });
    let mut dec = Decider::default();
    let app = match catallax_app::build_certified(&recipe, mesh_options.as_ref(), &mut dec) {
        Ok(app) => app,
        Err(error) => return Ok(app_refusal(&error)),
    };
    respond(&app, &recipe, &request, out, &mut dec)
}

/// The ordinary driver's answer for a body with fitted faces: the same manifest/4
/// record as the opt-in wrapper, because manifest/3 promises exact geometry.
pub(crate) fn fitted_recipe(
    recipe: &Recipe,
    built: catallax_features::Built,
    out: &Path,
    dec: &mut Decider,
) -> std::io::Result<Value> {
    let app = match catallax_app::certify_built(built, None, dec) {
        Ok(app) => app,
        Err(error) => return Ok(app_refusal(&error)),
    };
    let request: Request = serde_json::from_value(json!({
        "format":FORMAT,"recipe":Value::Null,
        "outputs":{"step":{"file_tolerance":"1/10000000","decimal_places":15}},
        "properties":{"relative":"1/1000000000","absolute":"1/1000000000"}}))
    .expect("fixed ordinary outputs");
    respond(&app, recipe, &request, out, dec)
}

fn respond(
    app: &CertifiedApplicationBody,
    recipe: &Recipe,
    request: &Request,
    out: &Path,
    dec: &mut Decider,
) -> std::io::Result<Value> {
    let built = app.built().built();
    let summary = match crate::result_summary(
        &built.model,
        built.result,
        dec,
        crate::SummaryExchange::Certified,
    ) {
        Ok(summary) => summary,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let record = app.built().result_record();
    let exact = record
        .regions
        .values()
        .all(|region| region.exception.is_none())
        && record.decisions.is_empty();
    let display = display(app, request.outputs.mesh.as_ref(), out)?;
    let properties = properties(app, request.properties.as_ref(), dec);
    let field =
        field(app, &request.field_points, dec).unwrap_or_else(|error| consumer_refusal(&error));
    let step = step(app, request.outputs.step.as_ref(), out, dec)?;
    let edit = edit(app, recipe, request.edit.as_ref(), dec)
        .unwrap_or_else(|error| consumer_refusal(&error));
    Ok(
        json!({"format":"catallax-manifest/4","status":"answered","exact":exact,"result":summary,
        "certificate":certificate::record(app.built()),"names":names(app),"display":display,"properties":properties,"field":field,"step":step,"edit":edit,
        "marks":built.model.body(built.result).marks.iter().map(|mark|crate::mark_entry(mark,None)).collect::<Vec<_>>(),
        "refusal":null,"replay":"question.json"}),
    )
}
