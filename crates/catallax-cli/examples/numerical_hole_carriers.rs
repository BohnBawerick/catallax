//! Retain numerical cone/hole-carrier fits from actual recipe prefixes.
//! Usage: numerical_hole_carriers OUT WORKERS RECIPE.json...

use catallax_numbers::{
    Decider, Rational, Real, Sign, current_operation_work, with_operation_budget,
};
use catallax_shapes::{BSplineCurve, BSplineSurface, Cylinder, Point3, Surface, Vector3};
use serde_json::{Value, json};
use std::{path::Path, process::ExitCode};

fn rational(value: &Value) -> Result<Rational, String> {
    Rational::parse(value.as_str().ok_or("expected a rational string")?)
        .map_err(|_| "invalid rational string".into())
}

fn triple(value: &Value) -> Result<[Rational; 3], String> {
    let values = value
        .as_array()
        .ok_or("expected three rational coordinates")?;
    if values.len() != 3 {
        return Err("expected three rational coordinates".into());
    }
    Ok([
        rational(&values[0])?,
        rational(&values[1])?,
        rational(&values[2])?,
    ])
}

fn position(point: &Point3) -> [String; 3] {
    [&point.x, &point.y, &point.z].map(|q| {
        q.as_rational()
            .expect("retained rational control")
            .to_string()
    })
}

fn source(surface: &BSplineSurface) -> Value {
    json!({"degrees":surface.degrees(),"control_count":surface.control_count(),
        "knots":surface.knots().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()),
        "controls":surface.controls().iter().map(position).collect::<Vec<_>>(),
        "weights":surface.weights().iter().map(|q|q.as_rational().expect("retained rational weight").to_string()).collect::<Vec<_>>()})
}

fn curve(curve: &BSplineCurve) -> Value {
    json!({"degree":curve.degree(),"knots":curve.knots().iter().map(ToString::to_string).collect::<Vec<_>>(),
        "controls":curve.controls().map(position).collect::<Vec<_>>(),
        "weights":curve.weights().map(|q|q.as_rational().expect("retained rational weight").to_string()).collect::<Vec<_>>()})
}

fn query(recipe: &Value, face_index: usize, reverse: bool) -> Result<Value, String> {
    let mut prefix = recipe.clone();
    let hole = prefix["steps"]
        .as_array_mut()
        .ok_or("recipe needs steps")?
        .pop()
        .ok_or("empty recipe")?;
    if hole["op"] != "hole" {
        return Err("the last recipe step must be a hole".into());
    }
    prefix["result"] = hole["target"].clone();
    let parsed = catallax_recipe::Recipe::parse(&prefix.to_string()).map_err(|e| e.to_string())?;
    let mut dec = Decider::default();
    let built = catallax_features::evaluate(&parsed, &mut dec).map_err(|e| e.to_string())?;
    let prefix_work = current_operation_work().used;
    let cones: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter_map(
            |id| match built.model.geom().surface(built.model.face(id).surface) {
                Surface::Cone(cone) => Some((id, cone.clone())),
                _ => None,
            },
        )
        .collect();
    let (face, cone) = cones.get(face_index).ok_or("prefix needs two cone faces")?;
    let at = hole["at"]
        .as_array()
        .filter(|at| at.len() == 1)
        .ok_or("example expects one hole placement")?;
    let point = triple(&at[0]["point"])?;
    let direction = triple(&at[0]["direction"])?;
    let diameter = rational(&hole["diameter"])?;
    if diameter.sign() != Sign::Positive {
        return Err("hole diameter must be positive".into());
    }
    let reference = (0..3).min_by_key(|&i| direction[i].abs()).unwrap();
    let reference: [Real; 3] = std::array::from_fn(|i| Real::from(i64::from(i == reference)));
    let [x, y, z] = point.map(Real::from);
    let [dx, dy, dz] = direction.map(Real::from);
    let [rx, ry, rz] = reference;
    let radius = diameter / Rational::from_i64(2);
    let cylinder = Surface::Cylinder(Cylinder {
        origin: Point3::new(x, y, z),
        axis: Vector3::new(dx, dy, dz),
        radius_sq: Real::from(&radius * &radius),
        ref_dir: Vector3::new(rx, ry, rz),
    });
    let cone = Surface::Cone(cone.clone());
    let ordered = if reverse {
        [&cylinder, &cone]
    } else {
        [&cone, &cylinder]
    };
    let result = catallax_intersect::numerical_native_intersection(ordered, &mut dec)
        .map_err(|e| format!("{}: {}", e.code, e.detail))?
        .ok_or("unsupported carrier pair")?;
    let windows = result.windows();
    let stages = result.work_stages();
    Ok(
        json!({"format":"catallax-numerical-native-pair/1","proof":"numerical","exact":false,
        "category":"approximate","cone_face":format!("{face:?}"),"reverse":reverse,
        "sources":windows.sources().iter().map(source).collect::<Vec<_>>(),
        "cone_heights":windows.cone_heights().each_ref().map(ToString::to_string),
        "cylinder_parameters":windows.cylinder_parameters().each_ref().map(ToString::to_string),
        "pieces":result.pieces().iter().map(|p|json!({"kind":format!("{:?}",p.kind()),
            "curves":std::iter::once(p.curve()).chain(p.uv().iter()).map(curve).collect::<Vec<_>>()})).collect::<Vec<_>>(),
        "evidence":catallax_cli::numerical_evidence_summary(result.evidence()),
        "complete_body_bound":"unestablished","trace_checkpoints":format!("{:?}",result.trace_work()),
        "work_stages":{"recipe":prefix_work,"windows":stages.windows,"setup":stages.setup,
            "tracing":stages.tracing,"fit_checks":stages.fit_checks}}),
    )
}

fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 3 {
        return Err("usage: numerical_hole_carriers OUT WORKERS RECIPE.json...".into());
    }
    let out = Path::new(&args[0]);
    let workers: usize = args[1].parse().map_err(|_| "workers must be 1, 2 or 4")?;
    if ![1, 2, 4].contains(&workers) {
        return Err("workers must be 1, 2 or 4".into());
    }
    if out.exists() {
        return Err("choose a fresh output directory".into());
    }
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let mut jobs = Vec::new();
    for path in &args[2..] {
        let recipe: Value =
            serde_json::from_slice(&std::fs::read(path).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        let name = Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or("recipe needs a filename")?;
        for face in 0..2 {
            for reverse in [false, true] {
                jobs.push((
                    format!("{name}-{face}-{reverse}"),
                    recipe.clone(),
                    face,
                    reverse,
                ));
            }
        }
    }
    let mut rows = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers).map(|rank| {
            let jobs = &jobs;
            scope.spawn(move || jobs.iter().enumerate().filter(|(i, _)| i % workers == rank).map(|(index, (name, recipe, face, reverse))| {
                let (result, work) = with_operation_budget(20_000_000, || query(recipe, *face, *reverse));
                let mut response = match result {
                    Ok(response) if !work.exhausted => response,
                    Ok(_) => json!({"proof":"refused","refusal":"operation allowance exhausted"}),
                    Err(error) => json!({"proof":"refused","refusal":error}),
                };
                response["part"] = json!(name);
                response["work"] = json!({"used":work.used,"limit":20_000_000,"exhausted":work.exhausted});
                (index, name.clone(), response)
            }).collect::<Vec<_>>())
        }).collect();
        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .map_err(|_| "query worker panicked".to_owned())
            })
            .collect::<Result<Vec<_>, _>>()
            .map(|rows| rows.into_iter().flatten().collect::<Vec<_>>())
    })?;
    rows.sort_by_key(|row| row.0);
    let mut failed = false;
    for (_, name, response) in &rows {
        failed |= response["proof"] != "numerical";
        std::fs::write(
            out.join(format!("{name}.json")),
            serde_json::to_vec_pretty(response).unwrap(),
        )
        .map_err(|e| e.to_string())?;
    }
    if failed {
        return Err("a carrier query refused; inspect its retained report".into());
    }
    println!(
        "{} numerical carrier queries completed at {workers} workers",
        rows.len()
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}
