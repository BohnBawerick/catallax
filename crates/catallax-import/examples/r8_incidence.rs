//! Private source audit for saved R8 incidence refusals.

use catallax_import::{Curve, ImportRefusal, Placement, Surface, Vec3, read_step};
use serde_json::{Value, json};

fn point(p: &Vec3) -> Value {
    json!(p.each_ref().map(ToString::to_string))
}

fn placement(p: &Placement) -> Value {
    json!({"origin":point(&p.origin),"axis":point(&p.axis),"reference":point(&p.ref_dir)})
}

fn curve(c: &Curve) -> Value {
    match c {
        Curve::Line(o, d) => json!({"kind":"line","origin":point(o),"direction":point(d)}),
        Curve::Circle(p, r) => {
            json!({"kind":"circle","placement":placement(p),"radius":r.to_string()})
        }
        Curve::Ellipse(p, a, b) => {
            json!({"kind":"ellipse","placement":placement(p),"radii":[a.to_string(),b.to_string()]})
        }
        Curve::BSpline(s) => json!({"kind":"spline","degree":s.geometry.degree(),
            "knots":s.geometry.knots().iter().map(ToString::to_string).collect::<Vec<_>>(),
            "controls":s.geometry.controls().iter().map(point).collect::<Vec<_>>(),
            "weights":s.geometry.weights().iter().map(ToString::to_string).collect::<Vec<_>>()}),
    }
}

fn surface(s: &Surface) -> Value {
    match s {
        Surface::Plane(p) => json!({"kind":"plane","placement":placement(p)}),
        Surface::Cylinder(p, r) => {
            json!({"kind":"cylinder","placement":placement(p),"radius":r.to_string()})
        }
        Surface::Sphere(p, r) => {
            json!({"kind":"sphere","placement":placement(p),"radius":r.to_string()})
        }
        Surface::Cone45(p, r) => {
            json!({"kind":"cone45","placement":placement(p),"radius":r.to_string()})
        }
        Surface::Cone(p, r, a) => {
            json!({"kind":"cone","placement":placement(p),"radius":r.to_string(),"angle":a.value.to_string(),"degrees":a.degrees})
        }
        Surface::Torus(p, a, b) => {
            json!({"kind":"torus","placement":placement(p),"radii":[a.to_string(),b.to_string()]})
        }
        Surface::BSpline(s) => {
            json!({"kind":"spline","degrees":s.geometry.degrees(),"counts":s.geometry.control_count()})
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let bytes = std::fs::read(args.next().ok_or("expected STEP path")?)?;
    let saved: ImportRefusal = serde_json::from_slice(&std::fs::read(
        args.next().ok_or("expected saved refusal")?,
    )?)?;
    let mode = args.next();
    if mode.as_deref() == Some("--measure") {
        let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || {
            catallax_import::measure_incidence(&bytes, &saved)
        });
        let result = if work.exhausted {
            Err(ImportRefusal::ResourceLimit)
        } else {
            result
        };
        println!(
            "{}",
            json!({"measurement":result,"work_used":work.used,"work_limit":20_000_000})
        );
        return Ok(());
    }
    let file = read_step(&bytes).map_err(|e| format!("{e:?}"))?;
    if mode.as_deref() == Some("--identities") {
        let mut pairs = std::collections::BTreeMap::new();
        let mut seen = std::collections::BTreeSet::new();
        let mut duplicates = Vec::new();
        for edge in file
            .solids
            .iter()
            .flat_map(|s| &s.faces)
            .flat_map(|f| &f.loops)
            .flat_map(|b| &b.edges)
        {
            let pair = if edge.start.point <= edge.end.point {
                (edge.start.point.clone(), edge.end.point.clone())
            } else {
                (edge.end.point.clone(), edge.start.point.clone())
            };
            if let Some(old) = pairs.get(&pair) {
                let old: &catallax_import::OrientedEdge = old;
                if old.edge != edge.edge && seen.insert((old.edge.clone(), edge.edge.clone())) {
                    duplicates.push(json!({"first":old.edge,"second":edge.edge,"curve_equal":old.curve==edge.curve,"first_kind":curve(&old.curve)["kind"],"second_kind":curve(&edge.curve)["kind"]}));
                }
            } else {
                pairs.insert(pair, edge.clone());
            }
        }
        println!("{}", json!({"duplicate_endpoint_pairs":duplicates}));
        return Ok(());
    }
    let ImportRefusal::NotSolid { detail } = saved else {
        return Err("expected incidence refusal".into());
    };
    let names = detail
        .split_whitespace()
        .filter(|s| s.starts_with('#'))
        .collect::<Vec<_>>();
    let edge_name = names.first().ok_or("missing edge name")?;
    let mut uses = Vec::new();
    for solid in &file.solids {
        for face in &solid.faces {
            for edge in face
                .loops
                .iter()
                .flat_map(|b| &b.edges)
                .filter(|e| e.edge == *edge_name)
            {
                uses.push(
                    json!({"solid":solid.name,"face":face.name,"surface":surface(&face.surface),
                    "edge":edge.edge,"start":point(&edge.start.point),"end":point(&edge.end.point),
                    "same_sense":edge.curve_same_sense,"curve":curve(&edge.curve)}),
                );
            }
        }
    }
    if uses.is_empty() {
        return Err("saved edge absent from source".into());
    }
    println!("{}", json!({"detail":detail,"uses":uses}));
    Ok(())
}
