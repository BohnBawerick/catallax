//! Development-only sampled observer shared by all three B-rep producers.
use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, distance::Options};
use serde_json::{Value, json};
use std::{env, fs};

fn observe(input: &Value) -> Result<Value, Box<dyn std::error::Error>> {
    let recipe = Recipe::parse(&serde_json::to_string(&input["recipe"])?)?;
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec)?;
    if input["action"] == "samples" {
        let mut points = Vec::new();
        let mut unresolved = 0usize;
        let mut workspace = Workspace::default();
        let mut failure = None;
        let completed = plan.visit_candidates(4, |point, face| {
            match plan.face_query_checked(face, &point, &mut workspace, &mut dec) {
                Ok(query) => match query.state {
                    FaceState::OnFace => points.push(json!([
                        point.x.to_string(),
                        point.y.to_string(),
                        point.z.to_string()
                    ])),
                    FaceState::NotOnFace => (),
                    _ => unresolved += 1,
                },
                Err(error) => {
                    failure = Some(error.to_string());
                    return false;
                }
            }
            true
        })?;
        return Ok(
            json!({"status":"answered","points":points,"completed":completed,
            "unresolved_candidates":unresolved,"failure":failure,
            "scope":"trim-checked definition-native rational candidates, four divisions; finite sampling, not complete boundary coverage"}),
        );
    }
    let options = Options {
        width: Rational::parse(
            input["distance_width_mm"]
                .as_str()
                .ok_or("width required")?,
        )?,
        max_cells: 256,
    };
    let prepared = plan.prepare_boundary_distances(&mut dec)?;
    let mut distances = Vec::new();
    for coordinates in input["points"].as_array().ok_or("points required")? {
        let values = coordinates.as_array().ok_or("point array required")?;
        if values.len() != 3 {
            return Err("three coordinates required".into());
        }
        let coordinate = |index: usize| -> Result<Rational, Box<dyn std::error::Error>> {
            Ok(Rational::parse(
                values[index].as_str().ok_or("exact coordinate required")?,
            )?)
        };
        let point = RPoint3::new(coordinate(0)?, coordinate(1)?, coordinate(2)?);
        let distance = prepared.point_to_boundary(&point.to_real(), &options, &mut dec)?;
        distances.push(json!({"lower_mm":distance.lower.to_string(),
            "upper_mm":distance.upper.map(|q|q.to_string()),"converged":distance.converged}));
    }
    Ok(json!({"status":"answered","distances":distances,
        "scope":"certified point-to-definition distance intervals at supplied samples; not a global distance certificate"}))
}

fn refusal_report(error: &(dyn std::error::Error + 'static)) -> Value {
    let mut report = json!({"status":"unresolved","detail":error.to_string()});
    if let Some(refusal) = error.downcast_ref::<catallax_numbers::Refusal>() {
        report["code"] = json!(refusal.code);
        if let Some(undecided) = refusal.undecided() {
            report["undecided"] = json!({"code":undecided.code,
                "detail":undecided.detail,"bound":undecided.bound.to_string()});
        }
    }
    report
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: kernel_compare INPUT.json OUTPUT.json".into());
    }
    let input: Value = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ((result, work), elapsed) = {
        let start = std::time::Instant::now();
        let result = with_operation_budget(100_000_000, || observe(&input));
        (result, start.elapsed())
    };
    let mut report = match result {
        Ok(value) if !work.exhausted => value,
        Ok(_) => json!({"status":"unresolved","detail":"observer work allowance exhausted"}),
        Err(error) => refusal_report(error.as_ref()),
    };
    report["work"] = json!(work.used);
    report["wall_ms"] = json!(elapsed.as_secs_f64() * 1000.0);
    fs::write(&args[1], serde_json::to_vec_pretty(&report)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request(result: &str, extra: Vec<Value>) -> Value {
        let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
        let mut steps = vec![json!({"id":"ball","op":"sphere","radius":"2","pose":pose})];
        steps.extend(extra);
        json!({"action":"distances","distance_width_mm":"1/1000","points":[["1","0","0"]],
            "recipe":{"format":"catallax-recipe/3","name":"observer coverage control","steps":steps,"result":result}})
    }

    fn sheet() -> Value {
        json!({"id":"sheet","op":"gyroid","region":"ball","period":"4","level":{"constant":"1"},"kind":"sheet",
            "frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}})
    }

    #[test]
    fn sheet_boundary_is_not_a_certified_named_sphere_distance() {
        let input = request("sheet", vec![sheet()]);
        // g(1,0,0)=1 and dg/dy=pi/2: this is a genuine sheet boundary,
        // whereas the surviving named sphere witness is one millimetre away.
        let recipe = Recipe::parse(&input["recipe"].to_string()).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let point = RPoint3::new(Rational::one(), Rational::zero(), Rational::zero());
        assert_eq!(
            plan.proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
                .unwrap()
                .location,
            catallax_field::Location::Boundary
        );
        let (answer, work) = with_operation_budget(100_000_000, || observe(&input));
        assert!(!work.exhausted);
        assert!(
            answer.is_err(),
            "observer falsely certified a whole-boundary distance: {answer:?}"
        );
        let error = answer.unwrap_err();
        let refusal = error
            .downcast_ref::<catallax_numbers::Refusal>()
            .expect("typed refusal");
        assert_eq!(refusal.code, "solver.distance-coverage");
        let report = refusal_report(error.as_ref());
        assert_eq!(report["status"], "unresolved");
        assert_eq!(report["code"], "solver.distance-coverage");
        assert!(report.get("distances").is_none());
    }

    #[test]
    fn covered_controls_and_compacted_implicit_leaves_remain_supported() {
        let sphere = observe(&request("ball", vec![])).unwrap();
        assert_eq!(sphere["status"], "answered");
        assert_eq!(sphere["distances"][0]["lower_mm"], "1");
        assert_eq!(sphere["distances"][0]["upper_mm"], "1");
        assert_eq!(sphere["distances"][0]["converged"], true);
        assert_eq!(observe(&request("ball", vec![sheet()])).unwrap(), sphere);
        let alias = json!({"id":"alias","op":"sphere","radius":"2",
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}});
        assert_eq!(
            observe(&request("alias", vec![sheet(), alias])).unwrap(),
            sphere
        );
        let folded = vec![
            sheet(),
            json!({"id":"empty","op":"boolean","kind":"subtract","target":"sheet","tool":"sheet"}),
            json!({"id":"result","op":"boolean","kind":"union","target":"ball","tool":"empty"}),
        ];
        assert_eq!(observe(&request("result", folded)).unwrap(), sphere);
        let block = json!({"id":"block","op":"block","size":["4","4","4"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-2","-2","-2"]}});
        assert_eq!(observe(&request("block", vec![block])).unwrap(), sphere);
    }

    #[test]
    fn other_live_unnamed_boundaries_refuse() {
        let mut network = sheet();
        network["kind"] = json!("network");
        let lattice = json!({"id":"sheet","op":"lattice","region":"ball","cell":"cubic","size":"2",
            "radius":{"constant":"0.25"},"frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}});
        for implicit in [network, lattice] {
            let error = observe(&request("sheet", vec![implicit])).unwrap_err();
            assert_eq!(
                refusal_report(error.as_ref())["code"],
                "solver.distance-coverage"
            );
        }
    }

    #[test]
    fn coverage_check_preserves_named_distance_work_and_refusal_evidence() {
        let (exhausted, work) = with_operation_budget(0, || observe(&request("ball", vec![])));
        assert!(work.exhausted);
        let error = exhausted.expect_err("exhaustion must not expose distances");
        assert_eq!(refusal_report(error.as_ref())["status"], "unresolved");
        let recipe = Recipe::parse(&request("ball", vec![])["recipe"].to_string()).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let query = |covered: bool| {
            with_operation_budget(100_000_000, || {
                let mut dec = Decider::default();
                let prepared = if covered {
                    plan.prepare_boundary_distances(&mut dec)
                } else {
                    plan.prepare_distances(&mut dec)
                }
                .unwrap();
                prepared
                    .point_to_boundary(
                        &RPoint3::new(Rational::one(), Rational::zero(), Rational::zero())
                            .to_real(),
                        &Options {
                            width: Rational::parse("1/1000").unwrap(),
                            max_cells: 256,
                        },
                        &mut dec,
                    )
                    .unwrap()
            })
        };
        let (original, old_work) = query(false);
        let (guarded, new_work) = query(true);
        assert_eq!(original.lower, guarded.lower);
        assert_eq!(original.upper, guarded.upper);
        assert_eq!(original.converged, guarded.converged);
        assert_eq!(old_work.used, new_work.used);
        assert_eq!(old_work.exhausted, new_work.exhausted);
        let error = catallax_numbers::Refusal::new("observer.test", "test refusal").keeping(
            catallax_numbers::NotExact {
                bound: Rational::one(),
                code: "numbers.budget-exhausted",
                detail: "test evidence".into(),
            },
        );
        assert_eq!(
            refusal_report(&error)["undecided"],
            json!({"code":"numbers.budget-exhausted","detail":"test evidence","bound":"1"})
        );
    }
}
