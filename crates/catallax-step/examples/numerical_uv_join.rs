//! Actual native pcurve joins and body handover controls on a 3 mm cube.
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Point3, Pose, RationalPatch, SplineSupport, Surface,
};
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::{Model, Sense, numerical as n};
use serde_json::json;
use std::{path::Path, process::ExitCode};
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
/// Builds a native join control and writes its numerical artifacts.
pub fn run(out: &Path, control: &str) -> Result<serde_json::Value, String> {
    let tolerance = q("1/100");
    let mut dec = Decider::default();
    let mut model = Model::new();
    let body = catallax_topology::block(&mut model, [q("3"), q("3"), q("3")], &Pose::identity())
        .map_err(|e| format!("{e:?}"))?;
    let mut definition =
        n::definition_from_model(&model, body, "native join cube".into(), &mut dec)
            .map_err(|e| format!("{e:?}"))?;
    let face = &mut definition.faces[0];
    let rounded = control == "rounded-gap";
    if rounded {
        let Surface::Plane(plane) = face.surface.clone() else {
            panic!("plane")
        };
        let v = plane.normal.cross(&plane.ref_dir);
        let mut points = Vec::new();
        for u in [-2, 2] {
            for t in [-2, 2] {
                points.push(
                    plane
                        .origin
                        .plus(&plane.ref_dir.scaled(&Real::from(3 * u)))
                        .plus(&v.scaled(&Real::from(3 * t))),
                );
            }
        }
        let source = BSplineSurface::new(
            [1, 1],
            [
                [-2, -2, 2, 2].map(Rational::from_i64).to_vec(),
                [-2, -2, 2, 2].map(Rational::from_i64).to_vec(),
            ],
            [2, 2],
            points,
            vec![Real::one(); 4],
            &mut dec,
        )
        .map_err(|e| format!("{e:?}"))?;
        face.surface = Surface::BSpline(SplineSupport::RationalPatch(
            RationalPatch::new(source, &mut dec).map_err(|e| format!("{e:?}"))?,
        ));
        for trim in face.loops.iter_mut().flatten() {
            for uv in &mut trim.uv {
                for x in uv {
                    *x = &*x / &q("3");
                }
            }
        }
    }
    let row = &face.loops[0];
    let (first, axis) = row
        .iter()
        .enumerate()
        .find_map(|(i, t)| (0..2).find(|&a| t.uv[1][a].is_zero()).map(|a| (i, a)))
        .ok_or("need zero-valued join coordinate")?;
    let uses = [first, (first + 1) % row.len()];
    let shift = if rounded {
        q("1/625")
    } else if control == "source-gap" {
        q("9/1000")
    } else {
        q("0")
    };
    let mut bindings = Vec::new();
    for (k, use_index) in uses.into_iter().enumerate() {
        let trim = &row[use_index];
        let edge = &definition.edges[trim.edge];
        let mut ends = trim.uv.clone();
        if trim.sense == Sense::Reversed {
            ends.reverse();
        }
        let mut domain = edge.parameters.clone();
        if domain[0] > domain[1] {
            domain.reverse();
            ends.reverse();
        }
        let delta = if k == 0 { shift.clone() } else { -&shift };
        for at in &mut ends {
            at[axis] = &at[axis] + &delta;
        }
        let curve = BSplineCurve::new(
            1,
            vec![
                domain[0].clone(),
                domain[0].clone(),
                domain[1].clone(),
                domain[1].clone(),
            ],
            ends.into_iter()
                .map(|[u, v]| Point3::new(u.into(), v.into(), 0.into()))
                .collect(),
            vec![Real::one(); 2],
            &mut dec,
        )
        .map_err(|e| format!("{e:?}"))?;
        bindings.push(n::NumericalNativeTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            face: face.key.clone(),
            boundary: 0,
            use_index,
            edge: edge.key.clone(),
            curve,
        });
    }
    std::fs::write(out.join("control.json"),serde_json::to_vec_pretty(&json!({"face":face.key,"uses":uses,"axis":axis,"shift_parameter":shift.to_string(),"tolerance_mm":tolerance.to_string(),"source_join_gap_mm":if rounded{"0.0096"}else if control=="source-gap"{"0.018"}else{"0"},"emitted_join_gap_mm":if rounded{"0.012"}else if control=="source-gap"{"0.018"}else{"0"}})).unwrap()).map_err(|e|e.to_string())?;
    let body = n::build_with_native_trims(
        definition,
        bindings,
        n::NumericalOptions {
            tolerance_mm: tolerance,
            vertex_tolerance_mm: None,
        },
        &mut dec,
    )
    .map_err(|e| format!("source admission: {e:?}"))?;
    if control == "valid" {
        let previous = body
            .rebuild(body.definition().clone(), body.options().clone(), &mut dec)
            .map_err(|e| format!("{e:?}"))?;
        std::fs::write(
            out.join("body.handover"),
            catallax_numbers::transfer::encode(&previous, &mut dec)
                .map_err(|e| format!("{e:?}"))?,
        )
        .map_err(|e| e.to_string())?;
        std::fs::write(
            out.join("history.handover"),
            catallax_numbers::transfer::encode(&previous.history()[0], &mut dec)
                .map_err(|e| format!("{e:?}"))?,
        )
        .map_err(|e| e.to_string())?;
    }
    let options = NumericalStepOptions {
        significant_digits: if rounded { 1 } else { 17 },
        ..Default::default()
    };
    let file =
        step::build(&body, &options, &mut dec).map_err(|e| format!("decimal admission: {e:?}"))?;
    std::fs::write(out.join("answer.step"), file.bytes()).map_err(|e| e.to_string())?;
    Ok(
        json!({"status":"numerical","exact":false,"control":control,"tolerance_mm":body.options().tolerance_mm.to_string(),"source_checks":format!("{:?}",body.evidence()),"emitted_checks":format!("{:?}",file.report().emitted_evidence)}),
    )
}
fn main() -> ExitCode {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 || !matches!(args[1].as_str(), "valid" | "source-gap" | "rounded-gap") {
        eprintln!("usage: numerical_uv_join OUT valid|source-gap|rounded-gap");
        return ExitCode::from(2);
    }
    let out = Path::new(&args[0]);
    std::fs::create_dir_all(out).unwrap();
    match std::fs::remove_file(out.join("answer.step")) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => panic!("{e}"),
    }
    let (result, work) = with_operation_budget(20_000_000, || run(out, &args[1]));
    let ok = result.is_ok() && !work.exhausted;
    let response = match result {
        Ok(v) => v,
        Err(e) => json!({"status":"refused","reason":e}),
    };
    let response = json!({"result":response,"work":{"used":work.used,"exhausted":work.exhausted,"limit":20_000_000}});
    std::fs::write(
        out.join("response.json"),
        serde_json::to_vec_pretty(&response).unwrap(),
    )
    .unwrap();
    println!("{}", response["result"]["status"]);
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
