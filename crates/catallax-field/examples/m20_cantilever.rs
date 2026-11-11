//! UNPROVED M20 aligned-grid experiment. See docs/M20-CANTILEVER.md.
#[path = "m20/elasticity.rs"]
mod elasticity;
use catallax_field::{
    Location, Workspace,
    printing::{BoxState, BoxStorage},
    solver::{Condition, Definition, MaterialAt, Materials},
};
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};
const RECIPE: &str = include_str!("m20/recipe.json");
const ALLOWANCE: u64 = 100_000_000;
fn q(v: i64) -> Rational {
    Rational::from(v)
}
fn fail(s: &str) -> Refusal {
    Refusal::new("m20.coverage", s)
}
fn refusal_json(r: Refusal) -> Value {
    json!({"code":r.code,"detail":r.detail,"undecided":r.undecided().map(|u|json!({"code":u.code,"detail":u.detail,"bound":u.bound.to_string()}))})
}
fn geometry(n: [usize; 3], dec: &mut Decider) -> Result<Value, Refusal> {
    let recipe = Recipe::parse(RECIPE).map_err(|_| fail("recipe parse"))?;
    let def = Definition::compile(
        &recipe,
        Materials::Homogeneous {
            region: "beam".into(),
            material: "E200000-nu0.3".into(),
        },
        dec,
    )?;
    let queries = def.prepare(dec)?;
    let half = Rational::from_ratio(5, 2).unwrap();
    let expected = [
        [q(0), q(100)],
        [-half.clone(), half.clone()],
        [-half.clone(), half.clone()],
    ];
    if queries.bounds() != &expected {
        return Err(fail("unexpected geometry bounds"));
    }
    let centre = RPoint3::new(q(50), q(0), q(0));
    if queries
        .membership(&centre, &mut Workspace::default(), dec)?
        .location
        != Location::Inside
    {
        return Err(fail("membership"));
    }
    if !matches!(queries.material(&centre,&mut Workspace::default(),dec)?, MaterialAt::Region{region,material} if region=="beam" && material=="E200000-nu0.3")
    {
        return Err(fail("material"));
    }
    let mut volume = q(0);
    let mut unknown = 0;
    let mut inside = 0;
    let mut storage = BoxStorage::default();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let indices = [k, i, j];
                let counts = [n[2], n[0], n[1]];
                let b = std::array::from_fn(|a| {
                    let h = Rational::from_ratio([100, 5, 5][a], counts[a] as i64).unwrap();
                    std::array::from_fn(|side| {
                        &expected[a][0] + &(&h * &q((indices[a] + side) as i64))
                    })
                });
                let cell = queries.integration(&b, &mut storage, dec)?;
                match cell.state {
                    BoxState::Inside => inside += 1,
                    BoxState::Unknown => unknown += 1,
                    BoxState::Outside => return Err(fail("outside cell")),
                }
                if !cell.connected_box
                    || !cell.fraction.is_exact()
                    || cell.fraction.lo != q(1)
                    || cell.pieces.len() != 1
                    || !cell.pieces[0].certain
                    || cell.pieces[0].bounds != b
                    || cell.pieces[0].region != "beam"
                    || cell.pieces[0].material != "E200000-nu0.3"
                {
                    return Err(fail("uncertified full homogeneous cell"));
                }
                // Moment 000 is physical volume; all moments are exact for this certificate.
                let moments = cell.moments("beam", dec)?;
                if moments.iter().any(|m| !m.is_exact()) {
                    return Err(fail("inexact moment"));
                }
                volume = volume + moments[0].lo.clone();
            }
        }
    }
    if volume != q(2500) {
        return Err(fail("volume"));
    }
    let names = [Role::Start, Role::End].map(|role| Origin {
        step: StepId::new("beam").unwrap(),
        role,
    });
    let fixed = queries.box_patch(
        &names[0],
        &RPoint3::new(q(0), q(0), q(0)),
        Condition::Fixed,
        dec,
    )?;
    let loaded = queries.box_patch(
        &names[1],
        &RPoint3::new(q(100), q(0), q(0)),
        Condition::Traction([q(0), q(0), Rational::from_ratio(-1, 25).unwrap()]),
        dec,
    )?;
    if fixed.area != q(25)
        || loaded.area != q(25)
        || fixed.centroid != [q(0), q(0), q(0)]
        || loaded.centroid != [q(100), q(0), q(0)]
    {
        return Err(fail("named patches"));
    }
    Ok(
        json!({"revision":def.revision(),"closed_inside":inside,"closed_unknown_exact_box_certificate":unknown,"volume_mm3":2500,"support":"beam/Start","load":"beam/End","material":"E200000-nu0.3","patch_area_mm2":25,"force_N":[0,0,-1],"moment_about_origin_N_mm":[0,100,0]}),
    )
}
fn beam(x: f64) -> f64 {
    x * x * (300. - x) / (6. * 200000. * (625. / 12.))
}
fn traces(n: [usize; 3], u: &[f64], x: f64, z: f64) -> Value {
    let internal = [2.5, z + 2.5, x];
    let size = [5., 5., 100.];
    let mut samples = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let c = [i, j, k];
                let lo: [f64; 3] = std::array::from_fn(|a| c[a] as f64 * size[a] / n[a] as f64);
                let hi: [f64; 3] =
                    std::array::from_fn(|a| (c[a] + 1) as f64 * size[a] / n[a] as f64);
                if (0..3).any(|a| internal[a] < lo[a] || internal[a] > hi[a]) {
                    continue;
                }
                let local =
                    std::array::from_fn(|a| 2. * (internal[a] - lo[a]) / (hi[a] - lo[a]) - 1.);
                let tensor = elasticity::tensor(n, u, c, local);
                let vm = elasticity::vm(tensor);
                let down = -elasticity::uz(n, u, c, local);
                let bending = (100. - x) * z.abs() / (625. / 12.);
                samples.push(json!({"cell_yzx":c,"tensor_xx_yy_zz_xy_yz_xz":tensor,"von_mises":vm,"down_mm":down,"deflection_signed_error_mm":down-beam(x),"vm_signed_error":vm-bending,"vm_absolute_error":(vm-bending).abs()}));
            }
        }
    }
    let vals: Vec<f64> = samples
        .iter()
        .map(|s| s["von_mises"].as_f64().unwrap())
        .collect();
    json!({"x_mm":x,"z_mm":z,"beam_down_mm":beam(x),"bending_reference":(100.-x)*z.abs()/(625./12.),"vm_range":[vals.iter().copied().fold(f64::INFINITY,f64::min),vals.iter().copied().fold(0.,f64::max)],"one_sided_traces":samples})
}
fn run(m: usize, workers: usize, allowance: u64, iterations: usize) -> Value {
    let n = [m, m, 10 * m];
    let (answer, work) = with_operation_budget(allowance, || -> Result<Value, Refusal> {
        if ![2, 4, 8].contains(&m) || workers == 0 {
            return Err(fail("unsupported run policy"));
        }
        let mut dec = Decider::default();
        let certificate = geometry(n, &mut dec)?;
        let solution = elasticity::solve(
            n,
            [5., 5., 100.],
            200000.,
            0.3,
            -1.,
            1e-9,
            iterations,
            workers,
            &mut dec,
        )?;
        // Fixed postprocessing reservation: one cell-centre evaluation and up to
        // eight incident-cell traces at each of three paths per axial station.
        dec.charge_work((n.iter().product::<usize>() + (n[2] + 1) * 3 * 8) as u64)?;
        let mut raw = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let c = [i, j, k];
                    let s = elasticity::tensor(n, &solution.displacement, c, [0.; 3]);
                    if s.iter().any(|v| !v.is_finite()) || !elasticity::vm(s).is_finite() {
                        return Err(fail("nonfinite stress"));
                    }
                    raw.push(json!({"cell_yzx":c,"centre_xyz":[(k as f64+0.5)*100./n[2] as f64,(i as f64+0.5)*5./n[0] as f64-2.5,(j as f64+0.5)*5./n[1] as f64-2.5],"tensor_xx_yy_zz_xy_yz_xz":s,"von_mises":elasticity::vm(s)}));
                }
            }
        }
        let profiles: Vec<_> = (0..=n[2])
            .flat_map(|k| {
                [-2.5, 0., 2.5]
                    .map(|z| traces(n, &solution.displacement, k as f64 * 100. / n[2] as f64, z))
            })
            .collect();
        let delta = -solution.tip;
        let r = solution.reaction;
        let moment = solution.reaction_moment;
        Ok(
            json!({"status":"UNPROVED","grid_yzx":n,"certificate":certificate,"tip_down_mm":delta,"tip_reference_mm":0.032,"tip_signed_error_mm":delta-0.032,"tip_absolute_error_mm":(delta-0.032).abs(),"tip_relative_error":delta/0.032-1.,"tip_within_2_percent":(delta/0.032-1.).abs()<=0.02,"reaction_xyz_N":[r[2],r[0],r[1]],"reaction_moment_xyz_N_mm":[moment[2],moment[0],moment[1]],"force_balance_xyz_N":[r[2],r[0],r[1]-1.],"moment_balance_xyz_N_mm":[moment[2],moment[0]+100.,moment[1]],"iterations":solution.iterations,"actual_relative_residual":solution.residual,"matrix_applications":solution.applications,"matrix_entry_products":solution.products,"displacement_internal_yzx":solution.displacement,"profiles":profiles,"raw_cell_centres":raw}),
        )
    });
    let answer = match answer {
        Ok(a) if !work.exhausted => a,
        Ok(_) => {
            json!({"refusal":refusal_json(Refusal::new("m20.exhausted","request allowance exhausted"))})
        }
        Err(r) => json!({"refusal":refusal_json(r)}),
    };
    json!({"answer":answer,"work_used":work.used,"work_exhausted":work.exhausted,"allowance":allowance})
}
fn picture(report: &Value) -> String {
    let a = &report["answer"];
    let mut svg = String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"650\" viewBox=\"0 0 1100 650\"><rect width=\"1100\" height=\"650\" fill=\"white\"/><g font-family=\"sans-serif\" fill=\"black\"><text x=\"35\" y=\"28\" font-size=\"20\">UNPROVED — M20 pipeline demonstration</text><text x=\"35\" y=\"52\">Computed cell-centre von Mises (N/mm²), one cell slice beside y=0; unsmoothed</text>",
    );
    if let Some(raw) = a["raw_cell_centres"].as_array() {
        let n: Vec<usize> = a["grid_yzx"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as usize)
            .collect();
        for cell in raw {
            let c = cell["cell_yzx"].as_array().unwrap();
            if c[0].as_u64().unwrap() as usize != n[0] / 2 {
                continue;
            }
            let k = c[2].as_u64().unwrap() as f64;
            let j = c[1].as_u64().unwrap() as f64;
            let vm = cell["von_mises"].as_f64().unwrap();
            let f = (vm / 6.).clamp(0., 1.);
            let red = (255. * f) as u8;
            let blue = (255. * (1. - f)) as u8;
            svg.push_str(&format!("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"rgb({red},40,{blue})\"><title>cell {k},{j}: {vm:.9} N/mm2</title></rect>",50.+k*1000./n[2] as f64,90.+(n[1] as f64-1.-j)*140./n[1] as f64,1000./n[2] as f64,140./n[1] as f64));
        }
        svg.push_str(&format!("<text x=\"50\" y=\"258\">x: 0 → 100 mm; z: -2.5 → 2.5 mm (vertical scale exaggerated). Blue=0, red≥6 N/mm².</text><text x=\"50\" y=\"285\">Grid {:?}; tip {:.9} mm; EB 0.032 mm; signed error {:.4}%</text>",n,a["tip_down_mm"].as_f64().unwrap(),100.*a["tip_relative_error"].as_f64().unwrap()));
    }
    if let Some(profiles) = a["profiles"].as_array() {
        let mut computed = String::new();
        let mut reference = String::new();
        for p in profiles.iter().filter(|p| p["z_mm"].as_f64() == Some(0.)) {
            let x = 50. + 10. * p["x_mm"].as_f64().unwrap();
            let down = p["one_sided_traces"][0]["down_mm"].as_f64().unwrap();
            computed.push_str(&format!("{x},{} ", 600. - 3000. * down));
            reference.push_str(&format!(
                "{x},{} ",
                600. - 3000. * p["beam_down_mm"].as_f64().unwrap()
            ));
        }
        svg.push_str(&format!("<text x=\"50\" y=\"435\">Centreline downward deflection: blue=computed, red=Euler–Bernoulli comparison only.</text><text x=\"50\" y=\"458\">x=0 to 100 mm; vertical displacement scale 3000 px/mm; zero at bottom.</text><polyline points=\"{reference}\" fill=\"none\" stroke=\"red\"/><polyline points=\"{computed}\" fill=\"none\" stroke=\"blue\"/><text x=\"50\" y=\"630\">Converged discrete solve; 2% tip mark: {}. Geometry revision: {}</text>",if a["tip_within_2_percent"].as_bool()==Some(true){"MET"}else{"NOT MET"},a["certificate"]["revision"].as_str().unwrap()));
    }
    svg.push_str("<text x=\"50\" y=\"315\">E=200000 N/mm², nu=0.3; beam/Start fixed; beam/End traction -0.04 N/mm² in z.</text><text x=\"50\" y=\"342\">Clamp/load end regions and corners are not beam-theory stress acceptance locations.</text><text x=\"50\" y=\"369\">Raw tensors, both outer-fibre one-sided trace ranges and centreline: companion JSON.</text><text x=\"50\" y=\"396\">Not a certified stress bound, yield assessment or safety sign-off.</text></g></svg>");
    svg
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let m = args.get(1).map_or(Ok(2), |v| v.parse::<usize>())?;
    let workers = args
        .get(2)
        .filter(|v| v.as_str() != "default")
        .map_or_else(
            || Ok(std::thread::available_parallelism().map_or(1, usize::from)),
            |v| v.parse::<usize>(),
        )?;
    let output = args
        .get(3)
        .ok_or("usage: m20_cantilever 2|4|8 workers|default output-stem")?;
    let report = run(m, workers, ALLOWANCE, 10000);
    std::fs::write(format!("{output}.json"), serde_json::to_vec(&report)?)?;
    std::fs::write(format!("{output}.svg"), picture(&report))?;
    eprintln!(
        "workers={workers} grid_m={m} work={} tip={} error={} iterations={} residual={} refusal={}",
        report["work_used"],
        report["answer"]["tip_down_mm"],
        report["answer"]["tip_relative_error"],
        report["answer"]["iterations"],
        report["answer"]["actual_relative_residual"],
        report["answer"]["refusal"]
    );
    if report["answer"]["refusal"].is_object() {
        return Err("typed refusal retained in output".into());
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn worker_answers_and_refusals() {
        for (cap, iterations) in [
            (100_000_000, 10000),
            (1, 10000),
            (100_000, 10000),
            (100_000_000, 1),
        ] {
            let baseline = run(2, 1, cap, iterations);
            for workers in [2, 3, 4, 17] {
                assert_eq!(run(2, workers, cap, iterations), baseline);
            }
            if cap == 100_000_000 && iterations == 10000 {
                assert!(baseline["answer"]["tip_down_mm"].is_number());
            } else {
                assert!(baseline["answer"]["refusal"].is_object());
            }
        }
    }
    #[test]
    fn recipe_also_builds_brep() {
        let r = Recipe::parse(RECIPE).unwrap();
        let built = catallax_features::evaluate(&r, &mut Decider::default())
            .expect("same native recipe builds B-rep");
        let measured = catallax_measure::volume_area_centroid(
            &built.model,
            built.result,
            &catallax_measure::Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        assert!(measured.volume.is_exact());
        assert_eq!(measured.volume.lo, q(2500));
        for (v, expect) in measured.centroid.unwrap().iter().zip([50, 0, 0]) {
            assert!(v.is_exact());
            assert_eq!(v.lo, q(expect));
        }
        for role in [Role::Start, Role::End] {
            let name = Origin {
                step: StepId::new("beam").unwrap(),
                role,
            };
            assert_eq!(built.model.faces_with_origin(built.result, &name).len(), 1);
        }
    }
}
