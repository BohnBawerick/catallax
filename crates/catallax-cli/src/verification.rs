//! Versioned deterministic verification exchange. Legacy commands do not use it.
use std::{fs, path::PathBuf, time::Instant};

use catallax_numbers::{Decider, with_operation_budget};
use serde::Deserialize;
use serde_json::{Value, json};

/// Default complete-request allowance, in arithmetic and checkpoint units.
pub const DEFAULT_WORK: u64 = 100_000_000;

/// One synchronous stream request. Workers run independent requests.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    /// Must be `catallax-request/1`.
    pub format: String,
    /// `question`, `recipe`, `field`, `compare`, or `name-edit`.
    pub source: String,
    /// Original question or recipe object.
    pub problem: Value,
    /// Owned directory for this request's manifest and STEP artifacts.
    pub out: PathBuf,
    /// Complete-operation work allowance, independent of machine speed.
    #[serde(default = "default_work")]
    pub max_work: u64,
}
fn default_work() -> u64 {
    DEFAULT_WORK
}

/// Portable application-loop question, without artifact paths.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExportedQuestion {
    /// Must be `catallax-question/1`.
    pub format: String,
    /// Verification adapter.
    pub source: String,
    /// Complete question or recipe.
    pub problem: Value,
    /// Original deterministic allowance.
    pub max_work: u64,
}

impl ExportedQuestion {
    /// Binds the portable question to an owned artifact directory.
    pub fn request(self, out: PathBuf) -> std::io::Result<Request> {
        if self.format != "catallax-question/1" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "expected catallax-question/1",
            ));
        }
        Ok(Request {
            format: "catallax-request/1".into(),
            source: self.source,
            problem: self.problem,
            out,
            max_work: self.max_work,
        })
    }
}

fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":"catallax-manifest/3", "status":"refused", "exact":false,
        "result":null, "properties":null, "step":null, "marks":[],
        "refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}

fn request_refused(request: &Request, code: &str, detail: impl ToString) -> Value {
    if request.format == "catallax-request/1"
        && request.source == "recipe"
        && crate::certified_recipe::is_problem(&request.problem)
    {
        crate::certified_recipe::refused(code, detail)
    } else {
        refused(code, detail)
    }
}

/// Answers one request and publishes a deterministic manifest plus timing.
/// A panic remains an internal error, never a geometric refusal.
pub fn answer(request: &Request) -> std::io::Result<Value> {
    let started = Instant::now();
    // Intersection plans retain shared expression nodes in thread-local caches.
    // A fresh thread prevents previous requests from changing memoized work.
    let mut response = std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("catallax-request".into())
            .stack_size(8 * 1024 * 1024)
            .spawn_scoped(scope, || answer_isolated(request))?
            .join()
            .map_err(|_| std::io::Error::other("request thread panicked"))?
    })?;
    response["wall_ms"] = json!(started.elapsed().as_secs_f64() * 1000.0);
    Ok(response)
}

fn answer_isolated(request: &Request) -> std::io::Result<Value> {
    crate::prepare_output(&request.out)?;
    let ((result, work), entries) = catallax_numbers::with_operation_log(|| {
        with_operation_budget(request.max_work, || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| evaluate(request)))
        })
    });
    let mut manifest = match result {
        Ok(result) => result?,
        Err(_) => {
            let mut value = request_refused(request, "cli.internal-error", "kernel panicked");
            value["status"] = json!("internal-error");
            value
        }
    };
    if work.exhausted && manifest["status"] != "internal-error" {
        let display_diagnostics = manifest.get("display_diagnostics").cloned();
        manifest = request_refused(
            request,
            "numbers.operation-budget",
            "operation work allowance exhausted",
        );
        if let Some(diagnostics) = display_diagnostics {
            manifest["display_diagnostics"] = diagnostics;
        }
    }
    if manifest["status"] != "answered" {
        crate::prepare_output(&request.out)?;
    }
    manifest["kernel"] = serde_json::to_value(crate::kernel_info())?;
    let mut measured_faces = manifest["face_measurements"]["faces"]
        .as_array()
        .into_iter()
        .flatten();
    let operations: Vec<_> = entries.iter().enumerate().map(|(sequence, entry)| {
        let decisions = entry.after.as_ref().map(|after| {
            let before = &entry.before;
            json!({
                "total":after.decisions.saturating_sub(before.decisions),
                "filter":after.filter.saturating_sub(before.filter),
                "refined":after.refined.saturating_sub(before.refined),
                "exact":after.exact.saturating_sub(before.exact),
                "coincident_exact":after.coincident_exact.saturating_sub(before.coincident_exact),
                "coincident_separation":after.coincident_separation.saturating_sub(before.coincident_separation),
                "coincident_fact":after.coincident_fact.saturating_sub(before.coincident_fact),
                "not_exact":after.not_exact.saturating_sub(before.not_exact),
                "refused":after.refused.saturating_sub(before.refused),
                "max_bits_high_water":after.max_bits
            })
        });
        let marks: Vec<_> = entry.marks.iter().map(|mark| crate::mark_entry(mark, Some(&entry.id))).collect();
        json!({"sequence":sequence,"kernel":entry.kernel,"id":entry.id,"operation":entry.operation,
            "status":if entry.after.is_none() {"interrupted"} else if entry.refusal.is_some() {"refused"} else {"answered"},
            "decisions":decisions,"marks":marks,"local_tolerances":if entry.operation == "face_measure" {
                measured_faces.next().map(|face| json!([face]))
            } else { None },
            "refusal":entry.refusal.as_ref().map(|r| json!({"code":r.code,"detail":r.detail})),
            "work":{"before":entry.work_before,"after":entry.work_after,"used":entry.work_after-entry.work_before},
            "wall_ms":entry.wall_ms})
    }).collect();
    let log = json!({"format":"catallax-operation-log/1","status":manifest["status"],
        "refusal":manifest["refusal"],"marks":manifest["marks"],"operations":operations});
    fs::write(
        request.out.join("operations.json"),
        serde_json::to_vec_pretty(&log)?,
    )?;
    fs::write(
        request.out.join("question.json"),
        serde_json::to_vec_pretty(&json!({
            "format":"catallax-question/1","source":request.source,"problem":request.problem,
            "max_work":request.max_work
        }))?,
    )?;
    let response = json!({"format":"catallax-response/1", "manifest":manifest,
        "operation_log":log,
        "work":{"used":work.used,"limit":request.max_work,"exhausted":work.exhausted},
        "wall_ms":0});
    fs::write(
        request.out.join("manifest.json"),
        serde_json::to_vec_pretty(&response["manifest"])?,
    )?;
    Ok(response)
}

fn evaluate(request: &Request) -> std::io::Result<Value> {
    if request.format != "catallax-request/1" {
        let mut value = refused("cli.invalid-request", "expected catallax-request/1");
        value["status"] = json!("invalid-question");
        return Ok(value);
    }
    match request.source.as_str() {
        "name-edit" => {
            let report = crate::name_checks::check(&request.problem);
            Ok(json!({"format":"catallax-manifest/3",
                "status":if report["status"] == "invalid" {"invalid-question"} else {"answered"},
                "exact":false,"marks":[],"refusal":null,"step":null,"name_edit":report}))
        }
        "question" => {
            let text = serde_json::to_string(&request.problem)?;
            let manifest = match crate::parse_question(&text) {
                Ok(question) => crate::answer_deterministic(&question, &request.out)?,
                Err(error) => crate::failure(
                    error.question_id,
                    crate::Status::InvalidQuestion,
                    "cli.invalid-question",
                    error.detail,
                ),
            };
            let mut value: Value = serde_json::from_str(&crate::render_manifest(&manifest))?;
            value["format"] = json!("catallax-manifest/3");
            Ok(value)
        }
        "recipe" if crate::certified_recipe::is_problem(&request.problem) => {
            crate::certified_recipe::evaluate(&request.problem, &request.out)
        }
        "recipe" | "compare" | "field" => {
            let display_started = Instant::now();
            let checks = if request.problem["format"] == "catallax-face-check/1" {
                if request.source != "recipe" {
                    return Ok(refused(
                        "cli.face-definition",
                        "face checks require the recipe adapter",
                    ));
                }
                match crate::face_checks::FaceChecks::parse(&request.problem) {
                    Ok(checks) => Some(checks),
                    Err(error) => return Ok(refused(error.code, error.detail)),
                }
            } else {
                None
            };
            let problem = checks.as_ref().map_or(&request.problem, |c| &c.recipe);
            let recipe = match catallax_recipe::Recipe::parse(&problem.to_string()) {
                Ok(recipe) => recipe,
                Err(error) => {
                    let mut value = refused("cli.invalid-recipe", error);
                    value["status"] = json!("invalid-question");
                    return Ok(value);
                }
            };
            if request.source == "field" {
                let mut dec = Decider::default();
                let plan = match catallax_field::Plan::compile_checked(&recipe, &mut dec) {
                    Ok(plan) => plan,
                    Err(error) => return Ok(refused(error.code, error.detail)),
                };
                let options = catallax_field::display::Options::default();
                let span = catallax_numbers::OperationSpan::begin(
                    "field",
                    &recipe.result,
                    "display",
                    &dec,
                );
                let mut diagnostics = catallax_field::display::Diagnostics::default();
                let mesh = catallax_field::display::mesh_plan_reported(
                    &plan,
                    &options,
                    &mut dec,
                    &mut diagnostics,
                );
                span.finish(&dec, mesh.as_ref().err(), None);
                let first_mesh_ms = display_started.elapsed().as_secs_f64() * 1000.0;
                let diagnostics = display_diagnostics(&plan, &diagnostics);
                let mesh = match mesh {
                    Ok(mesh) => mesh,
                    Err(error) => {
                        let mut manifest = refused(error.code, error.detail);
                        manifest["display_diagnostics"] = diagnostics;
                        return Ok(manifest);
                    }
                };
                let mut file = fs::File::create(request.out.join("display.obj"))?;
                mesh.write_obj(&mut file)?;
                return Ok(json!({
                    "format":"catallax-manifest/3","status":"answered","exact":true,
                    "field":{"faces":plan.faces().len(),"proven_empty":plan.is_proven_empty()},
                    "face_measurements":{
                        "format":"catallax-face-measurements/1","faces":[],"edit":null,
                        "display":{"first_mesh_ms":first_mesh_ms,"file":"display.obj",
                            "vertices":mesh.positions.len(),"triangles":mesh.triangles.len(),
                            "chord_tolerance_mm":options.tolerance.to_string(),
                            "deviation_bound_mm":mesh.max_deviation.to_string(),
                            "coverage":"two-sided distance to the defined result boundary",
                            "cells":mesh.cells}
                    },
                    "marks":[],"refusal":null,"step":null,"display_diagnostics":diagnostics
                }));
            }
            if request.source == "compare" {
                let controller = catallax_compare::Controller::default();
                controller
                    .begin()
                    .check(&recipe, &catallax_compare::Options::default());
                let report = controller
                    .latest()
                    .expect("synchronous comparison publishes");
                if let Some(refusal) = &report.refusal {
                    return Ok(refused(refusal.code, &refusal.detail));
                }
                let verdict = match report.verdict {
                    catallax_compare::Verdict::Agree => "agree",
                    catallax_compare::Verdict::Disagree => "disagree",
                    catallax_compare::Verdict::Unresolved => "unresolved",
                };
                let findings: Vec<_> = report.findings.iter().map(|finding| {
                    let (code, detail) = match &finding.reason {
                        catallax_compare::Reason::Gap => ("compare.gap", String::new()),
                        catallax_compare::Reason::Origins => ("compare.origins", String::new()),
                        catallax_compare::Reason::Direction => ("compare.direction", String::new()),
                        catallax_compare::Reason::Unresolved(detail) => ("compare.unresolved", detail.clone()),
                    };
                    json!({"code":code,"detail":detail,"pass":match finding.pass {
                        catallax_compare::Pass::Forward=>"forward",catallax_compare::Pass::Reverse=>"reverse"},
                        "witness":format!("{:?}",finding.point),
                        "brep_faces":format!("{:?}",finding.brep_faces),
                        "field_faces":format!("{:?}",finding.field_faces),
                        "gap":format!("{:?}",finding.gap)})
                }).collect();
                let unsupported: Vec<_> = report.not_comparable.iter().map(|entry| json!({"step":entry.step.to_string(),"operation":entry.operation,"reason":entry.reason})).collect();
                return Ok(
                    json!({"format":"catallax-manifest/3","status":"answered", "exact":false,
                    "comparison":{"verdict":verdict,"forward_vertices":report.forward_vertices,
                    "reverse_samples":report.reverse_samples,
                    "sampled_distance_upper_mm":report.sampled_distance_upper.as_ref().map(ToString::to_string),
                    "findings":findings,"not_comparable":unsupported},
                    "marks":[],"refusal":null,"step":null}),
                );
            }
            let mut dec = Decider::default();
            let mut built = match catallax_features::evaluate(&recipe, &mut dec) {
                Ok(built) => built,
                Err(catallax_features::BuildError::Geometry { step, refusal }) => {
                    let mut value = refused(refusal.code, refusal.detail);
                    value["refusal"]["step"] = json!(step.to_string());
                    return Ok(value);
                }
                Err(error) => return Ok(refused("cli.invalid-recipe", error)),
            };
            let measurements = if let Some(checks) = checks {
                match checks.report(&recipe, &mut built, display_started, &mut dec) {
                    Ok(report) => Some(report),
                    Err(error) => return Ok(refused(error.code, error.detail)),
                }
            } else {
                None
            };
            let span =
                catallax_numbers::OperationSpan::begin("brep", &recipe.result, "summarize", &dec);
            let summary = crate::result_summary(
                &built.model,
                built.result,
                &mut dec,
                crate::SummaryExchange::Analytic,
            );
            span.finish(&dec, summary.as_ref().err(), None);
            let summary = match summary {
                Ok(summary) => summary,
                Err(error) => return Ok(refused(error.code, error.detail)),
            };
            let mut step = Value::Null;
            let mut marks: Vec<_> = built
                .model
                .body(built.result)
                .marks
                .iter()
                .map(|mark| crate::mark_entry(mark, None))
                .collect();
            if summary.kind != "empty" {
                let mut file = fs::File::create(request.out.join("answer.step"))?;
                let span =
                    catallax_numbers::OperationSpan::begin("brep", &recipe.result, "step", &dec);
                let exported = catallax_step::write_step(
                    &built.model,
                    &[built.result],
                    &Default::default(),
                    &mut dec,
                    &mut file,
                );
                span.finish(
                    &dec,
                    exported.as_ref().err(),
                    exported.as_ref().ok().map(|answer| &answer.marks),
                );
                match exported {
                    Ok(answer) => {
                        marks.extend(
                            answer
                                .marks
                                .iter()
                                .map(|mark| crate::mark_entry(mark, None)),
                        );
                        let report = answer.value;
                        step = json!({"file":"answer.step","max_deviation":report.max_deviation.to_decimal_bounds(40).1,
                            "approximated_curves":report.approximated_curves,"added_vertices":report.added_vertices,"added_edges":report.added_edges});
                    }
                    Err(error) => return Ok(refused(error.code, error.detail)),
                }
            }
            let mut manifest = json!({"format":"catallax-manifest/3","status":"answered","exact":marks.is_empty(),
                "result":summary,"marks":marks,"step":step,"refusal":null});
            if let Some(measurements) = measurements {
                manifest["face_measurements"] = measurements;
            }
            Ok(manifest)
        }
        _ => {
            let mut value = refused(
                "cli.invalid-source",
                "use question, recipe, field, compare, or name-edit",
            );
            value["status"] = json!("invalid-question");
            Ok(value)
        }
    }
}

fn display_diagnostics(
    plan: &catallax_field::Plan,
    report: &catallax_field::display::Diagnostics,
) -> Value {
    let cost = |c: &catallax_field::display::Cost| {
        json!({
            "cells":c.cells,"max_depth":c.max_depth,"decisions":c.decisions,"work":c.work
        })
    };
    let faces:Vec<_> = report.faces.iter().map(|f|json!({
        "face":f.face.index(),"origin":plan.faces()[f.face.index()].origin,
        "candidate_cells":f.candidate_cells,"max_depth":f.max_depth,"queries":cost(&f.queries)
    })).collect();
    let supports:Vec<_> = report.supports.iter().enumerate().map(|(i,s)|json!({
        "support":i,"faces":s.faces.iter().map(|f|f.index()).collect::<Vec<_>>(),
        "origins":s.faces.iter().map(|f|&plan.faces()[f.index()].origin).collect::<Vec<_>>(),
        "total":cost(&s.total),"classification":cost(&s.classification),
        "queries":cost(&s.queries),"emission":cost(&s.emission)
    })).collect();
    json!({"format":"catallax-display-diagnostics/1","total":cost(&report.total),"common":cost(&report.common),
        "faces":faces,"supports":supports})
}
