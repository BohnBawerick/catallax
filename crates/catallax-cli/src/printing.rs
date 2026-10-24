//! Explicit field print outputs for the deterministic verification interface.
use catallax_field::{
    Plan,
    printing::{Mesh, Options},
};
use catallax_numbers::Decider;
use catallax_recipe::{Number, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{fs, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    format: String,
    recipe: Value,
    tolerance: Number,
    #[serde(default)]
    max_cells: Option<usize>,
    #[serde(default)]
    max_vertices: Option<usize>,
}
fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/3","status":"refused","exact":false,"result":null,"step":null,"marks":[],"refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}
pub(crate) fn evaluate(problem: &Value, out: &Path) -> std::io::Result<Value> {
    let (recipe, options) = if problem["format"] == "catallax-print/1" {
        let request: Request = match serde_json::from_value(problem.clone()) {
            Ok(request) => request,
            Err(error) => return Ok(refused("cli.invalid-print", error)),
        };
        if request.format != "catallax-print/1" {
            return Ok(refused("cli.invalid-print", "expected catallax-print/1"));
        }
        let defaults = Options::default();
        (
            request.recipe,
            Options {
                tolerance: request.tolerance.value().clone(),
                max_cells: request.max_cells.unwrap_or(defaults.max_cells),
                max_vertices: request.max_vertices.unwrap_or(defaults.max_vertices),
            },
        )
    } else {
        (problem.clone(), Options::default())
    };
    let recipe = match Recipe::parse(&recipe.to_string()) {
        Ok(recipe) => recipe,
        Err(error) => return Ok(refused("cli.invalid-recipe", error)),
    };
    let mut dec = Decider::default();
    let plan = match Plan::compile_checked(&recipe, &mut dec) {
        Ok(plan) => plan,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let span = catallax_numbers::OperationSpan::begin("field", &recipe.result, "print", &dec);
    let mesh = catallax_field::printing::mesh_plan(&plan, &options, &mut dec);
    span.finish(&dec, mesh.as_ref().err(), None);
    let mesh = match mesh {
        Ok(mesh) => mesh,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    publish(&plan, &options, out, &mesh)
}
fn publish(plan: &Plan, options: &Options, out: &Path, mesh: &Mesh) -> std::io::Result<Value> {
    write(mesh, out)?;
    Ok(
        json!({"format":"catallax-manifest/3","status":"answered","exact":false,
        "field":{"faces":plan.faces().len(),"proven_empty":plan.is_proven_empty()},
        "printing":{"file":"print.stl","format":"binary-stl","units":"mm","watertight":true,
            "vertices":mesh.positions.len(),"triangles":mesh.triangles.len(),"cells":mesh.cells,
            "tolerance_mm":options.tolerance.to_string(),"deviation_bound_mm":mesh.max_deviation.to_string(),
            "coverage":"two-sided result-boundary distance; interval-excluded cells and witnessed neighbouring boundary cells"},
        "face_measurements":{"format":"catallax-face-measurements/1","faces":[],"edit":null,
            "display":{"file":"display.obj","vertices":mesh.positions.len(),"triangles":mesh.triangles.len(),
                "chord_tolerance_mm":options.tolerance.to_string(),"deviation_bound_mm":mesh.max_deviation.to_string(),
                "coverage":"same watertight field mesh as print output","cells":mesh.cells}},
        "marks":[],"step":null,"refusal":null}),
    )
}
fn write(mesh: &Mesh, out: &Path) -> std::io::Result<()> {
    mesh.write_stl(&mut fs::File::create(out.join("print.stl"))?)?;
    mesh.write_obj(&mut fs::File::create(out.join("display.obj"))?)
}
