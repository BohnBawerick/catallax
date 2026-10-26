//! Qualified print outputs. Route selection never retries another representation.
use catallax_field::{
    Plan,
    printing::{Mesh, Options},
};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Number, Operation, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufWriter, Write},
    path::Path,
};
mod brep;
mod nesting;
mod qualification;

#[derive(Clone, Copy, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
enum Route {
    Brep,
    Field,
}
impl Route {
    fn name(self) -> &'static str {
        match self {
            Self::Brep => "brep",
            Self::Field => "field",
        }
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    format: String,
    recipe: Value,
    tolerance: Number,
    max_cells: Option<usize>,
    max_vertices: Option<usize>,
    route: Option<Route>,
    #[serde(default)]
    certified: bool,
}
fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/3","status":"refused","exact":false,"result":null,"step":null,"marks":[],"refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}
pub(crate) fn request_refused(
    problem: &Value,
    force_field: bool,
    code: &str,
    detail: impl ToString,
) -> Value {
    let mut options = Options::default();
    let mut explicit = None;
    let recipe = if problem["format"] == "catallax-print/1" {
        serde_json::from_value::<Request>(problem.clone())
            .ok()
            .map(|request| {
                options.tolerance = request.tolerance.value().clone();
                options.max_cells = request.max_cells.unwrap_or(options.max_cells);
                options.max_vertices = request.max_vertices.unwrap_or(options.max_vertices);
                explicit = request.route;
                request.recipe
            })
    } else {
        Some(problem.clone())
    };
    let route = if force_field {
        Route::Field
    } else {
        explicit.unwrap_or_else(|| {
            if recipe
                .and_then(|value| Recipe::parse(&value.to_string()).ok())
                .as_ref()
                .is_some_and(native_field)
            {
                Route::Field
            } else {
                Route::Brep
            }
        })
    };
    let mut result = refused(code, detail);
    result["printing"] = settings(route, &options);
    result
}
fn native_field(recipe: &Recipe) -> bool {
    let steps: BTreeMap<_, _> = recipe.steps.iter().map(|s| (&s.id, &s.op)).collect();
    let mut pending = vec![&recipe.result];
    let mut seen = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !seen.insert(id) {
            continue;
        }
        let op = steps[id];
        if matches!(
            op,
            Operation::Lattice { .. }
                | Operation::Gyroid { .. }
                | Operation::SmoothBlend { .. }
                | Operation::ImportMesh { .. }
        ) {
            return true;
        }
        pending.extend(op.inputs().into_iter().map(|(id, _)| id));
    }
    false
}
fn settings(route: Route, options: &Options) -> Value {
    json!({"route":route.name(),"file":null,"format":"binary-stl","metadata":"catallax-print-stl/1","units":"mm",
        "coordinate_encoding":"ieee754-binary32-le","tolerance_mm":options.tolerance.to_string(),
        "max_cells":options.max_cells,"cell_cap_applicable":route == Route::Field,"max_vertices":options.max_vertices,
        "cells":if route == Route::Brep {json!(0)} else {Value::Null},"vertices":null,"triangles":null,
        "encoding_rounding_bound_mm":null,"deviation_bound_mm":null,
        "component_preservation":"all components and cavities; no filtering or simplification",
        "completeness":"unproven","topology":"unproven","orientation":"unproven"})
}
pub(crate) fn evaluate(problem: &Value, out: &Path, force_field: bool) -> std::io::Result<Value> {
    let mut options = Options::default();
    let mut explicit = None;
    let mut certified = false;
    let recipe_value = if problem["format"] == "catallax-print/1" {
        let request: Request = match serde_json::from_value(problem.clone()) {
            Ok(request) => request,
            Err(error) => {
                let mut result = refused("cli.invalid-print", error);
                result["printing"] = settings(
                    if force_field {
                        Route::Field
                    } else {
                        Route::Brep
                    },
                    &options,
                );
                return Ok(result);
            }
        };
        debug_assert_eq!(request.format, "catallax-print/1");
        options.tolerance = request.tolerance.value().clone();
        options.max_cells = request.max_cells.unwrap_or(options.max_cells);
        options.max_vertices = request.max_vertices.unwrap_or(options.max_vertices);
        explicit = request.route;
        certified = request.certified;
        request.recipe
    } else {
        problem.clone()
    };
    let parsed = Recipe::parse(&recipe_value.to_string());
    let route = if force_field {
        Route::Field
    } else {
        explicit.unwrap_or_else(|| {
            if parsed.as_ref().is_ok_and(native_field) {
                Route::Field
            } else {
                Route::Brep
            }
        })
    };
    let mut record = settings(route, &options);
    let mut result = (|| -> std::io::Result<Value> {
        if force_field && explicit == Some(Route::Brep) {
            return Ok(refused(
                "cli.print-route",
                "source field cannot request the B-rep route",
            ));
        }
        let recipe = match parsed {
            Ok(recipe) => recipe,
            Err(error) => return Ok(refused("cli.invalid-recipe", error)),
        };
        if options.tolerance <= Rational::zero()
            || options.max_vertices < 4
            || (route == Route::Field && options.max_cells == 0)
        {
            return Ok(refused(
                "print.invalid-options",
                "positive accuracy, at least four vertices, and a positive field cell cap required",
            ));
        }
        let mut dec = Decider::default();
        let span =
            catallax_numbers::OperationSpan::begin(route.name(), &recipe.result, "print", &dec);
        let mut field = Value::Null;
        let output = match route {
            Route::Brep => brep::mesh(&recipe, &options, certified, &mut dec),
            Route::Field => {
                if certified {
                    Err(Refusal::new(
                        "cli.print-route",
                        "certified model bounds require the B-rep route",
                    ))
                } else {
                    Plan::compile_checked(&recipe,&mut dec).and_then(|plan| {
                    field = json!({"faces":plan.faces().len(),"proven_empty":plan.is_proven_empty()});
                    catallax_field::printing::mesh_plan(&plan,&options,&mut dec)
                }).map(|mesh| {
                    let tessellation_bound = mesh.max_deviation.clone();
                    brep::Output {mesh, model_bound:Rational::zero(),tessellation_bound,rounding_bound:Rational::zero()}
                })
                }
            }
        };
        span.finish(&dec, output.as_ref().err(), None);
        let mut output = match output {
            Ok(output) => output,
            Err(error) => return Ok(refused(error.code, error.detail)),
        };
        let header = match header(route, &options, &mut output) {
            Ok(header) => header,
            Err(error) => return Ok(refused(error.code, error.detail)),
        };
        record["file"] = json!("print.stl");
        record["watertight"] = json!(true);
        record["vertices"] = json!(output.mesh.positions.len());
        record["triangles"] = json!(output.mesh.triangles.len());
        record["cells"] = json!(output.mesh.cells);
        record["model_bound_mm"] = json!(output.model_bound.to_string());
        record["tessellation_bound_mm"] = json!(output.tessellation_bound.to_string());
        record["encoding_rounding_bound_mm"] = json!(output.rounding_bound.to_string());
        record["deviation_bound_mm"] = json!(output.mesh.max_deviation.to_string());
        record["completeness"] = json!("proven-two-sided");
        record["topology"] = json!("closed-opposite-edges-and-cyclic-vertex-links");
        record["orientation"] = json!("material-outward-including-cavities");
        record["coverage"] = json!(if route == Route::Brep {
            "complete certified finite face/trim correspondence; affine encoding isotopy"
        } else {
            "interval-excluded cells and witnessed neighbouring boundary cells"
        });
        let files = match publish(out, &output.mesh, &header, &record) {
            Ok(files) => files,
            Err(error) => {
                crate::prepare_output(out)?;
                return Err(error);
            }
        };
        record["sha256"] = files["print.stl"]["sha256"].clone();
        Ok(
            json!({"format":"catallax-manifest/3","status":"answered","exact":false,"marks":[],"step":null,"refusal":null,
            "field":field,"files":files,"face_measurements":{"format":"catallax-face-measurements/1","faces":[],"edit":null,
            "display":{"file":"display.obj","vertices":output.mesh.positions.len(),"triangles":output.mesh.triangles.len(),
                "chord_tolerance_mm":options.tolerance.to_string(),"deviation_bound_mm":output.mesh.max_deviation.to_string(),
                "coverage":"same encoded qualified print mesh","cells":output.mesh.cells}}}),
        )
    })()?;
    result["printing"] = record;
    Ok(result)
}
fn upward(bound: &Rational) -> Result<(f64, Rational), Refusal> {
    let mut value = bound.to_f64_nearest();
    if !value.is_finite() {
        return Err(Refusal::new(
            "print.metadata-range",
            "certificate bound exceeds binary64 metadata range",
        ));
    }
    if Rational::from_f64_exact(value).unwrap() < *bound {
        value = value.next_up();
    }
    Rational::from_f64_exact(value)
        .map(|q| (value, q))
        .ok_or_else(|| {
            Refusal::new(
                "print.metadata-range",
                "certificate bound exceeds binary64 metadata range",
            )
        })
}
fn header(route: Route, options: &Options, output: &mut brep::Output) -> Result<[u8; 80], Refusal> {
    let accuracy = options.tolerance.to_string();
    let (n, d) = accuracy.split_once('/').unwrap_or((&accuracy, "1"));
    let n: i64 = n.parse().map_err(|_| Refusal::new("print.metadata-range","accuracy numerator does not fit STL metadata; use an exactly representable requested accuracy"))?;
    let d: u64 = d.parse().map_err(|_| {
        Refusal::new(
            "print.metadata-range",
            "accuracy denominator does not fit STL metadata; requested accuracy was not truncated",
        )
    })?;
    let (u, model) = upward(&output.model_bound)?;
    let (v, tessellation) = upward(&output.tessellation_bound)?;
    let (r, rounding) = upward(&output.rounding_bound)?;
    let (w, total) = upward(&(&(&model + &tessellation) + &rounding))?;
    if total > options.tolerance {
        return Err(Refusal::new(
            "print.metadata-accuracy",
            "lossless conservative header bounds exceed requested accuracy",
        ));
    }
    output.model_bound = model;
    output.tessellation_bound = tessellation;
    output.rounding_bound = rounding;
    output.mesh.max_deviation = total;
    let mut bytes = [0u8; 80];
    bytes[..8].copy_from_slice(b"CXPRINT1");
    bytes[8] = if route == Route::Brep { 1 } else { 2 };
    bytes[9] = 1; // millimetres
    bytes[10] = 1; // binary32 little endian
    bytes[11] = if route == Route::Field { 1 } else { 0 }; // cell cap applicability
    bytes[12..20].copy_from_slice(&n.to_le_bytes());
    bytes[20..28].copy_from_slice(&d.to_le_bytes());
    bytes[28..36].copy_from_slice(&(options.max_cells as u64).to_le_bytes());
    bytes[36..44].copy_from_slice(&(options.max_vertices as u64).to_le_bytes());
    for (offset, value) in [(44, u), (52, v), (60, r), (68, w)] {
        bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    }
    Ok(bytes)
}
struct Hashed<W> {
    out: W,
    hash: Sha256,
    bytes: u64,
}
impl<W: Write> Write for Hashed<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let n = self.out.write(bytes)?;
        self.hash.update(&bytes[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.out.flush()
    }
}
struct SkipHeader<'a, W> {
    out: &'a mut W,
    remaining: usize,
}
impl<W: Write> Write for SkipHeader<'_, W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let skipped = bytes.len().min(self.remaining);
        self.remaining -= skipped;
        if skipped == bytes.len() {
            Ok(skipped)
        } else {
            self.out.write(&bytes[skipped..]).map(|n| n + skipped)
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.out.flush()
    }
}
fn publish(out: &Path, mesh: &Mesh, header: &[u8; 80], record: &Value) -> std::io::Result<Value> {
    let mut files = json!({});
    for name in ["print.stl", "display.obj"] {
        let mut sink = Hashed {
            out: BufWriter::new(fs::File::create(out.join(name))?),
            hash: Sha256::new(),
            bytes: 0,
        };
        if name == "print.stl" {
            sink.write_all(header)?;
            mesh.write_stl(&mut SkipHeader {
                out: &mut sink,
                remaining: 80,
            })?;
        } else {
            writeln!(
                sink,
                "# catallax-print/1 {}",
                serde_json::to_string(record)?
            )?;
            for point in &mesh.positions {
                writeln!(
                    sink,
                    "v {} {} {}",
                    f64::from(point[0]),
                    f64::from(point[1]),
                    f64::from(point[2])
                )?;
            }
            for &[a, b, c] in &mesh.triangles {
                writeln!(sink, "f {} {} {}", a + 1, b + 1, c + 1)?;
            }
        }
        sink.flush()?;
        files[name] = json!({"sha256":format!("{:x}",sink.hash.finalize()),"bytes":sink.bytes});
    }
    Ok(files)
}
