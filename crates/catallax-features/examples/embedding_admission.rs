//! Report exact construction, native embedding and numerical admission per recipe.
//! Usage: cargo run -p catallax-features --example embedding_admission -- [--inspect] RECIPE...

use catallax_numbers::{
    Ball, BudgetPolicy, Decider, Rational, current_operation_work, with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, Surface};
use catallax_topology::numerical::{NumericalOptions, build_with_native_trims};
use catallax_topology::{BodyId, Loop, Model};
use serde_json::json;
use std::{error::Error, fs, process::ExitCode};

fn source_topology(path: &str, model: &Model, body: BodyId) {
    for id in model.faces_of_body(body) {
        let face = model.face(id);
        let boundaries: Vec<_> = face.loops.iter().map(|lp| {
            let Loop::Coedges(uses) = model.loop_(*lp) else { return json!({"point_loop": true}) };
            json!(uses.iter().map(|coedge| {
                let coedge = model.coedge(*coedge);
                let edge = model.edge(coedge.edge);
                let kind = match model.geom().curve(edge.curve) {
                    Curve::Circle(_) => "circle", Curve::Line(_) => "line",
                    Curve::BSpline(_) => "spline", _ => "traced",
                };
                let geometry = if let Curve::Circle(circle) = model.geom().curve(edge.curve) {
                    let mut dec = Decider::default();
                    Some(json!({
                        "center": ([&circle.center.x, &circle.center.y, &circle.center.z]
                            .map(|value| dec.enclose_at(value, 128).map(|v| v.mid().to_f64_nearest()).ok())),
                        "radius_sq": dec.enclose_at(&circle.radius_sq, 128)
                            .map(|v| v.mid().to_f64_nearest()).ok(),
                    }))
                } else {
                    None
                };
                json!({"edge": format!("{:?}", coedge.edge), "curve": kind,
                    "sense": format!("{:?}", coedge.sense), "edge_sense": format!("{:?}", edge.sense),
                    "start": format!("{:?}", edge.start), "end": format!("{:?}", edge.end),
                    "geometry": geometry})
            }).collect::<Vec<_>>())
        }).collect();
        let surface = model.geom().surface(face.surface);
        let spline = match surface {
            Surface::BSpline(source) => {
                let domain = source.carrier().domain();
                let mut corners = Vec::with_capacity(4);
                for u in domain[0] {
                    for v in domain[1] {
                        let derivatives = source.carrier().derivatives(
                            &Ball::from_rational(u, 128), &Ball::from_rational(v, 128), &Decider::default(),
                        );
                        let derivatives = match derivatives {
                            Ok(jet) => json!(jet.map(|axis| axis.each_ref().map(|value| value.mid().to_f64_nearest()))),
                            Err(error) => json!({"error": error.to_string()}),
                        };
                        corners.push(json!({"uv": [u.to_string(), v.to_string()], "derivatives": derivatives}));
                    }
                }
                Some(json!({
                    "family_tag": format!("{:?}", std::mem::discriminant(source)),
                    "degrees": source.carrier().degrees(), "control_count": source.carrier().control_count(),
                    "knots": source.carrier().knots().map(|row| row.iter().map(ToString::to_string).collect::<Vec<_>>()),
                    "corners": corners,
                }))
            }
            _ => None,
        };
        println!("{}", json!({"recipe": path, "stage": "source-topology", "face": format!("{id:?}"),
            "surface": format!("{:?}", surface.kind()), "spline": spline, "boundaries": boundaries}));
    }
}

fn run(path: &str, inspect: bool) -> Result<bool, Box<dyn Error>> {
    let recipe = Recipe::parse(&fs::read_to_string(path)?)?;
    let (passed, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 100_000_000,
            ..Default::default()
        });
        let built = match catallax_features::evaluate(&recipe, &mut dec) {
            Ok(built) => built,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "construction", "error": error.to_string()}));
                return false;
            }
        };
        let construction_work = current_operation_work().used;
        if inspect {
            source_topology(path, &built.model, built.result);
            println!("{}", json!({"recipe": path, "stage": "inspected", "faces": built.model.faces_of_body(built.result).len(),
                "edges": built.model.edges_of_body(built.result).len(), "construction_work": construction_work}));
            return true;
        }
        let source_edges = built.model.edges_of_body(built.result);
        let tolerance = Rational::parse("0.0000001").unwrap();
        let (definition, native) = match catallax_features::numerical_embedding::definition_from_exact_model(
            &built.model, built.result, recipe.name.clone(), &tolerance, &mut dec,
        ) {
            Ok(definition) => definition,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "embedding", "error": error.to_string(), "construction_work": construction_work, "embedding_work": current_operation_work().used - construction_work}));
                source_topology(path, &built.model, built.result);
                return false;
            }
        };
        let embedding_end = current_operation_work().used;
        let edge_rows: Vec<_> = definition.edges.iter().zip(&source_edges).map(|(edge, source)| {
            let family = match built.model.geom().curve(built.model.edge(*source).curve) {
                Curve::Line(_) => "line",
                Curve::Circle(_) => "circle",
                Curve::BSpline(_) => "spline",
                _ => "traced",
            };
            json!({"edge": edge.key, "source_curve": family, "native_trims": native.iter().filter(|trim| trim.edge == edge.key).count()})
        }).collect();
        println!("{}", json!({
            "recipe": path, "stage": "embedded",
            "faces": definition.faces.len(), "edges": definition.edges.len(),
            "source_edges": source_edges.len(), "native_trims": native.len(),
            "construction_work": construction_work,
            "embedding_work": embedding_end - construction_work,
        }));
        let body = match build_with_native_trims(
            definition, native, NumericalOptions { tolerance_mm: tolerance }, &mut dec,
        ) {
            Ok(body) => body,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "admission", "error": error.to_string(), "construction_work": construction_work, "embedding_work": embedding_end - construction_work, "admission_work": current_operation_work().used - embedding_end}));
                return false;
            }
        };
        println!("{}", json!({
            "recipe": path, "stage": "admitted", "tolerance_mm": "0.0000001",
            "faces": body.definition().faces.len(), "edges": body.definition().edges.len(),
            "construction_work": construction_work, "embedding_work": embedding_end - construction_work,
            "admission_work": current_operation_work().used - embedding_end,
        }));
        for mut row in edge_rows {
            row["recipe"] = json!(path);
            row["stage"] = json!("admitted");
            println!("{row}");
        }
        true
    });
    println!("{}", json!({"recipe": path, "work": work.used, "exhausted": work.exhausted}));
    Ok(passed && !work.exhausted)
}

fn main() -> ExitCode {
    let mut paths: Vec<_> = std::env::args().skip(1).collect();
    let inspect = paths.first().is_some_and(|arg| arg == "--inspect");
    if inspect { paths.remove(0); }
    if paths.is_empty() {
        eprintln!("usage: embedding_admission [--inspect] RECIPE...");
        return ExitCode::FAILURE;
    }
    let mut passed = true;
    for path in paths {
        match run(&path, inspect) {
            Ok(result) => passed &= result,
            Err(error) => {
                eprintln!("{path}: {error}");
                passed = false;
            }
        }
    }
    if passed { ExitCode::SUCCESS } else { ExitCode::FAILURE }
}
