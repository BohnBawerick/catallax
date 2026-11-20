//! Reproduce two retained loft cuts and their numerical STEP evidence.

use catallax_numbers::{
    Decider, NumericalCheck, NumericalEvidence, Rational, current_operation_work,
    with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineCurve, Curve, Point3, SplineSupport, Surface};
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::{
    Sense,
    numerical::{self, NumericalOptions},
};
use serde_json::{Value, json};
use std::{path::Path, process::ExitCode};

const WORK: u64 = 20_000_000;

fn evidence(source: &NumericalEvidence) -> Value {
    let checks: Vec<_> = source.checks().iter().map(|check| match check {
        NumericalCheck::Passed {name,scope} => json!({"kind":"passed","name":name,"scope":scope}),
        NumericalCheck::Samples {name,scope,samples} => json!({"kind":"samples","name":name,"scope":scope,"samples":samples.iter().map(|sample| json!({"parameters":sample.parameters.iter().map(ToString::to_string).collect::<Vec<_>>(),"distance_upper_mm":sample.distance_upper_mm.to_string()})).collect::<Vec<_>>()}),
    }).collect();
    json!({"tolerance_mm":source.tolerance_mm().to_string(),"checks":checks,"repair_count":source.repairs().len()})
}

fn nonuniform_edge(
    definition: &mut numerical::NumericalDefinition,
    dec: &mut Decider,
) -> Result<(Value, Vec<numerical::NumericalNativeTrim>), String> {
    let (face, trim) = definition
        .faces
        .iter()
        .filter(|f| matches!(f.surface, Surface::BSpline(SplineSupport::RationalPatch(_))))
        .find_map(|f| {
            f.loops
                .iter()
                .flatten()
                .find(|t| {
                    t.uv[0][0] == t.uv[1][0]
                        && matches!(definition.edges[t.edge].curve, Curve::Line(_))
                })
                .map(|t| (f.key.clone(), t.clone()))
        })
        .ok_or("part has no linear generatrix on a curved face")?;
    let edge = &mut definition.edges[trim.edge];
    let endpoints = edge
        .vertices
        .map(|v| definition.vertices[v].position.clone());
    let coefficients = ["2", "-7/2", "9/2", "-1"].map(|q| Rational::parse(q).unwrap());
    let controls = coefficients
        .iter()
        .map(|c| {
            let p: [Rational; 3] = std::array::from_fn(|a| {
                &endpoints[0][a] + &(c * &(&endpoints[1][a] - &endpoints[0][a]))
            });
            let [x, y, z] = p;
            Point3::new(x.into(), y.into(), z.into())
        })
        .collect();
    edge.curve = Curve::BSpline(
        BSplineCurve::new(
            3,
            [-1, -1, -1, -1, 2, 2, 2, 2]
                .map(Rational::from_i64)
                .to_vec(),
            controls,
            vec![1.into(); 4],
            dec,
        )
        .map_err(|e| format!("{e:?}"))?,
    );
    edge.parameters = [0.into(), 1.into()];
    let key = edge.key.clone();
    let mut bindings = Vec::new();
    for face in &definition.faces {
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, used) in row.iter().enumerate() {
                if used.edge != trim.edge {
                    continue;
                }
                let mut uv = used.uv.clone();
                if used.sense == Sense::Reversed {
                    uv.reverse();
                }
                let controls = coefficients
                    .iter()
                    .map(|c| {
                        let p: [Rational; 2] =
                            std::array::from_fn(|a| &uv[0][a] + &(c * &(&uv[1][a] - &uv[0][a])));
                        Point3::new(p[0].clone().into(), p[1].clone().into(), 0.into())
                    })
                    .collect();
                let curve = BSplineCurve::new(
                    3,
                    [-1, -1, -1, -1, 2, 2, 2, 2]
                        .map(Rational::from_i64)
                        .to_vec(),
                    controls,
                    vec![1.into(); 4],
                    dec,
                )
                .map_err(|e| format!("{e:?}"))?;
                bindings.push(numerical::NumericalNativeTrim {
                    tag: catallax_topology::numerical::CopyTag::Independent,
                    face: face.key.clone(),
                    boundary,
                    use_index,
                    edge: key.clone(),
                    curve,
                });
            }
        }
    }
    Ok((
        json!({"face":face,"edge":key,"native_carrier_domain":["-1","2"],"finite_interval":["0","1"],"map":"t/2+3*t^2/2-t^3"}),
        bindings,
    ))
}

fn native_control(
    bindings: &mut Vec<numerical::NumericalNativeTrim>,
    control: &str,
    dec: &mut Decider,
) -> Result<(), String> {
    match control {
        "native-uv" => {}
        "native-uv-missing" => bindings.clear(),
        "native-uv-stale" => bindings[0].edge = "stale-native-edge".into(),
        "native-uv-duplicate" => bindings.push(bindings[0].clone()),
        "native-uv-reversed" | "native-uv-partial" | "native-uv-wrong" => {
            let source = &bindings[0].curve;
            let mut points = source.controls().cloned().collect::<Vec<_>>();
            let mut knots = source.knots().to_vec();
            if control == "native-uv-reversed" {
                points.reverse();
            }
            if control == "native-uv-wrong" {
                points[1].y =
                    &points[1].y + &catallax_numbers::Real::from(Rational::parse("1/100").unwrap());
            }
            if control == "native-uv-partial" {
                for t in &mut knots[source.degree() + 1..] {
                    *t = Rational::parse("1/2").unwrap();
                }
            }
            bindings[0].curve = BSplineCurve::new(
                source.degree(),
                knots,
                points,
                source.weights().cloned().collect(),
                dec,
            )
            .map_err(|e| format!("{e:?}"))?;
        }
        _ => return Err("unknown native UV control".into()),
    }
    Ok(())
}

fn native_summary(bindings: &[numerical::NumericalNativeTrim]) -> Value {
    json!(bindings.iter().map(|b|json!({"face":b.face,"boundary":b.boundary,"use_index":b.use_index,"edge":b.edge,
        "degree":b.curve.degree(),"knots":b.curve.knots().iter().map(ToString::to_string).collect::<Vec<_>>(),
        "controls":b.curve.controls().map(|p|[&p.x,&p.y].map(|q|q.as_rational().unwrap().to_string())).collect::<Vec<_>>(),
        "weights":b.curve.weights().map(|q|q.as_rational().unwrap().to_string()).collect::<Vec<_>>()
    })).collect::<Vec<_>>())
}

fn run(out: &Path, case: &str, control: &str) -> Result<Value, String> {
    let (name, size, at) = match case {
        "tiny" => (
            "loft crossing",
            ["10", "10", "1"],
            ["-3", "1999999/1000000", "1"],
        ),
        "m10-008" => ("m7-cut", ["2", "4", "1"], ["1", "1", "1"]),
        _ => return Err("unknown part".into()),
    };
    let mut input: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/loft-spline-prism.json"
    ))
    .map_err(|e| e.to_string())?;
    input["name"] = json!(name);
    let steps = input["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"box","op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":at}}));
    steps.push(json!({"id":"cut","op":"boolean","kind":"subtract","target":"loft","tool":"box"}));
    input["result"] = json!("cut");
    std::fs::write(
        out.join("recipe.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    let recipe = Recipe::parse(&input.to_string()).map_err(|e| format!("{e:?}"))?;
    let mut dec = Decider::default();
    let built = catallax_features::evaluate(&recipe, &mut dec).map_err(|e| format!("{e:?}"))?;
    let recipe_work = current_operation_work().used;
    let mut definition =
        numerical::definition_from_model(&built.model, built.result, name.into(), &mut dec)
            .map_err(|e| format!("{e:?}"))?;
    let mut native_trims = Vec::new();
    if control.starts_with("native-uv") {
        let (control_data, bindings) = nonuniform_edge(&mut definition, &mut dec)?;
        native_trims = bindings;
        native_control(&mut native_trims, control, &mut dec)?;
        std::fs::write(
            out.join("native-uv-control.json"),
            serde_json::to_vec_pretty(&control_data).unwrap(),
        )
        .map_err(|e| e.to_string())?;
    } else if control != "source" {
        let face = definition
            .faces
            .iter_mut()
            .find(|face| face.loops.len() == 2)
            .ok_or("part has no physical hole")?;
        match control {
            "reverse-hole" => {
                let inner = &mut face.loops[1];
                inner.reverse();
                for trim in inner {
                    trim.uv.reverse();
                    trim.sense = if trim.sense == Sense::Same {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    };
                }
            }
            "missing-hole" => {
                face.loops.pop();
            }
            _ => return Err("unknown control".into()),
        }
    }
    let body = numerical::build_with_native_trims(
        definition,
        native_trims,
        NumericalOptions {
            tolerance_mm: Rational::parse("1/10000000").unwrap(),
            vertex_tolerance_mm: None,
        },
        &mut dec,
    )
    .map_err(|e| format!("{e:?}"))?;
    let admitted_work = current_operation_work().used;
    let file = step::build(&body, &NumericalStepOptions::default(), &mut dec)
        .map_err(|e| format!("{e:?}"))?;
    std::fs::write(out.join("answer.step"), file.bytes()).map_err(|e| e.to_string())?;
    let reread = std::fs::read(out.join("answer.step")).map_err(|e| e.to_string())?;
    file.check_written(&reread, &mut dec)
        .map_err(|e| format!("{e:?}"))?;
    let report = file.report();
    let rounding = report.spline_rounding.as_ref().unwrap();
    Ok(json!({
        "status":"numerical","category":"approximate","exact":false,"case":case,"control":control,
        "counts":report.counts,"two_loop_faces":report.source_definition.faces.iter().filter(|face| face.loops.len()==2).count(),
        "bytes":file.bytes().len(),"source_evidence":evidence(&report.source_evidence),"emitted_evidence":evidence(&report.emitted_evidence),
        "native_trims":{"original":native_summary(&report.original_native_trims),"source":native_summary(&report.source_native_trims),"emitted":native_summary(&report.emitted_native_trims)},
        "file_checks":report.checks,"max_scalar_rounding_error_mm":report.max_scalar_rounding_error_mm.to_string(),
        "max_direction_ratio_rounding_error":report.max_direction_ratio_rounding_error.to_string(),
        "max_native_parameter_error":rounding.max_parameter_error.to_string(),"max_weight_error":rounding.max_weight_error.to_string(),
        "complete_body_bound":"unestablished",
        "work_stages":{"recipe":recipe_work,"numerical_admission":admitted_work-recipe_work,"step_and_reread":current_operation_work().used-admitted_work}
    }))
}

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(2..=4).contains(&args.len())
        || !matches!(args[1].as_str(), "tiny" | "m10-008")
        || args.get(2).is_some_and(|s| {
            !matches!(
                s.as_str(),
                "source"
                    | "reverse-hole"
                    | "missing-hole"
                    | "native-uv-missing"
                    | "native-uv"
                    | "native-uv-reversed"
                    | "native-uv-partial"
                    | "native-uv-wrong"
                    | "native-uv-stale"
                    | "native-uv-duplicate"
            )
        })
    {
        eprintln!(
            "usage: numerical_loft_cut OUTPUT_DIRECTORY tiny|m10-008 [source|reverse-hole|missing-hole|native-uv[-missing|-reversed|-partial|-wrong|-stale|-duplicate]] [1|2|4]"
        );
        return ExitCode::from(2);
    }
    if let Some(workers) = args.get(3) {
        let Ok(workers) = workers.parse::<usize>() else {
            return ExitCode::from(2);
        };
        if ![1, 2, 4].contains(&workers) {
            return ExitCode::from(2);
        }
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(workers));
    }
    let out = Path::new(&args[0]);
    let control = args.get(2).map_or("source", String::as_str);
    if let Err(error) = std::fs::create_dir_all(out) {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    if let Err(error) = std::fs::remove_file(out.join("answer.step"))
        && error.kind() != std::io::ErrorKind::NotFound
    {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    let (result, work) = with_operation_budget(WORK, || run(out, &args[1], control));
    let mut response = match result {
        Ok(value) if !work.exhausted => value,
        other => {
            let _ = std::fs::remove_file(out.join("answer.step"));
            let reason = match other {
                Err(reason) => reason,
                Ok(_) => "operation work allowance exhausted".into(),
            };
            json!({"status":"refused","exact":false,"case":args[1],"control":control,"reason":reason})
        }
    };
    response["work"] =
        json!({"used":work.used,"limit":WORK,"target":WORK,"exhausted":work.exhausted});
    if let Err(error) = std::fs::write(
        out.join("response.json"),
        serde_json::to_vec_pretty(&response).unwrap(),
    ) {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    println!(
        "{}",
        json!({"status":response["status"],"case":args[1],"control":control,"work":response["work"]})
    );
    if response["status"] == "numerical" {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
