//! Agreed M20 cantilever. All numerical elasticity output is UNPROVED.
#[path = "beam/solver.rs"]
mod elasticity;

use catallax_field::{
    printing::BoxStorage,
    solver::{Condition, Definition, Materials},
};
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use std::{fmt::Write, path::Path};

const RECIPE: &str = include_str!("beam/agreed.json");
const THEORY: f64 = -0.032;
const GRIDS: [[usize; 3]; 3] = [[2, 2, 20], [4, 4, 40], [8, 8, 80]];
fn q(n: i64) -> Rational {
    Rational::from(n)
}
fn half(n: i64) -> Rational {
    Rational::from_ratio(n, 2).unwrap()
}
fn fail(code: &'static str, detail: &str) -> Refusal {
    Refusal::new(code, detail)
}

#[derive(Debug, PartialEq)]
struct Run {
    estimate: elasticity::Estimate,
    geometry_work: u64,
    revision: String,
}

fn run(n: [usize; 3]) -> Result<Run, Refusal> {
    let (geometry, work) = with_operation_budget(100_000_000, || {
        let recipe = Recipe::parse(RECIPE).unwrap();
        let mut dec = Decider::default();
        let definition = Definition::compile(
            &recipe,
            Materials::Homogeneous {
                region: "beam".into(),
                material: "isotropic-E200000-nu0.3".into(),
            },
            &mut dec,
        )?;
        let queries = definition.prepare(&mut dec)?;
        let mut volume = q(0);
        // The solver's local (x,y,z) maps to world (y,z,x), a proper rotation.
        // Every assembled full cell must first have exact field integration.
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let index = [k, i, j];
                    let counts = [n[2], n[0], n[1]];
                    let bounds = std::array::from_fn(|axis| {
                        let h =
                            Rational::from_ratio([100, 5, 5][axis], counts[axis] as i64).unwrap();
                        let origin = if axis == 0 { q(0) } else { half(-5) };
                        [
                            &origin + &(&h * &q(index[axis] as i64)),
                            &origin + &(&h * &q(index[axis] as i64 + 1)),
                        ]
                    });
                    let cell =
                        queries.integration(&bounds, &mut BoxStorage::default(), &mut dec)?;
                    if !cell.connected_box || !cell.fraction.is_exact() || cell.fraction.lo != q(1)
                    {
                        return Err(fail(
                            "fea.cell",
                            &format!(
                                "cell={index:?} fraction=[{},{}] connected={}",
                                cell.fraction.lo, cell.fraction.hi, cell.connected_box
                            ),
                        ));
                    }
                    let moments = cell.moments("beam", &mut dec)?;
                    if !moments.iter().all(|m| m.is_exact()) {
                        return Err(fail(
                            "fea.moments",
                            &format!(
                                "cell={index:?} exact_moments={}/27",
                                moments.iter().filter(|m| m.is_exact()).count()
                            ),
                        ));
                    }
                    volume = volume + moments[0].lo.clone();
                }
            }
        }
        if volume != q(2500) {
            return Err(fail(
                "fea.volume",
                &format!("volume_mm3={volume} expected=2500"),
            ));
        }
        for (segment, x, condition) in [
            ("x_min", 0, Condition::Fixed),
            (
                "x_max",
                100,
                Condition::Traction([q(0), q(0), Rational::from_ratio(-1, 25).unwrap()]),
            ),
        ] {
            let origin = Origin {
                step: StepId::new("beam").unwrap(),
                role: Role::Side {
                    segment: ElementId::new(segment).unwrap(),
                    path: None,
                },
            };
            let patch = queries.box_patch(
                &origin,
                &RPoint3::new(q(x), q(0), q(0)),
                condition,
                &mut dec,
            )?;
            if patch.area != q(25) {
                return Err(fail(
                    "fea.face-area",
                    &format!("face={segment} area_mm2={} expected=25", patch.area),
                ));
            }
            if let Condition::Traction(t) = patch.condition
                && &t[2] * &patch.area != q(-1)
            {
                return Err(fail("fea.load", "expected resultant_z_N=-1"));
            }
        }
        Ok(definition.revision().to_owned())
    });
    let revision = geometry.map_err(|mut refusal: Refusal| {
        write!(
            refusal.detail,
            " geometry_work={} geometry_limit=100000000",
            work.used
        )
        .unwrap();
        refusal
    })?;
    if work.exhausted {
        return Err(fail(
            "fea.geometry-budget",
            &format!("used={} limit=100000000", work.used),
        ));
    }
    let estimate = elasticity::solve(n, [5., 5., 100.], 200_000., 0.3, -1., 1e-9, 10_000)?;
    Ok(Run {
        estimate,
        geometry_work: work.used,
        revision,
    })
}

// Tensor order in world axes: xx, yy, zz, xy, yz, xz.
fn world_tensor(s: [f64; 6]) -> [f64; 6] {
    [s[2], s[0], s[1], s[5], s[3], s[4]]
}
fn world_vector(v: [f64; 3]) -> [f64; 3] {
    [v[2], v[0], v[1]]
}
fn bending(x: f64) -> f64 {
    (100. - x) * 2.5 / (625. / 12.)
}

fn artifacts(dir: &Path, n: [usize; 3], result: &elasticity::Estimate) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let mut profile = String::from(
        "evidence,x_mm,y_mm,z_mm,cell_local_i,cell_local_j,cell_local_k,ux_mm,uy_mm,uz_mm,sxx,syy,szz,txy,tyz,txz,vm_N_per_mm2,EB_uz_mm,EB_bending_abs,uz_signed_error_mm,uz_abs_error_mm,vm_signed_error,vm_abs_error\n",
    );
    for station in 0..=80 {
        let x = station as f64 * 1.25;
        for z in [-2.5, 0., 2.5] {
            let traces = result.field.traces([2.5, z + 2.5, x]);
            let min = traces
                .iter()
                .map(|s| s.von_mises)
                .fold(f64::INFINITY, f64::min);
            let max = traces.iter().map(|s| s.von_mises).fold(0., f64::max);
            if [20, 40, 60].contains(&station) {
                println!(
                    "UNPROVED station_x_mm={x} z_mm={z} traces={} vm_range=[{min:.12e},{max:.12e}] beam_bending_abs={:.12e}",
                    traces.len(),
                    if z == 0. { 0. } else { bending(x) }
                );
            }
            let eb = -x * x * (300. - x) / (6. * 200_000. * (625. / 12.));
            let reference = if z == 0. { 0. } else { bending(x) };
            for trace in traces {
                let u = world_vector(trace.displacement);
                let s = world_tensor(trace.tensor);
                writeln!(profile,"UNPROVED,{x},0,{z},{},{},{},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{eb:.16e},{reference:.16e},{:.16e},{:.16e},{:.16e},{:.16e}",trace.cell[0],trace.cell[1],trace.cell[2],u[0],u[1],u[2],s[0],s[1],s[2],s[3],s[4],s[5],trace.von_mises,u[2]-eb,(u[2]-eb).abs(),trace.von_mises-reference,(trace.von_mises-reference).abs()).unwrap();
            }
        }
    }
    std::fs::write(dir.join("profiles.csv"), profile)?;
    let mut raw = String::from("evidence,x_mm,y_mm,z_mm,sxx,syy,szz,txy,tyz,txz,vm_N_per_mm2\n");
    let mut values = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let p = [
                    (i as f64 + 0.5) * 5. / n[0] as f64,
                    (j as f64 + 0.5) * 5. / n[1] as f64,
                    (k as f64 + 0.5) * 100. / n[2] as f64,
                ];
                let traces = result.field.traces(p);
                assert_eq!(traces.len(), 1);
                let trace = &traces[0];
                let s = world_tensor(trace.tensor);
                writeln!(
                    raw,
                    "UNPROVED,{},{},{},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e}",
                    p[2],
                    p[0] - 2.5,
                    p[1] - 2.5,
                    s[0],
                    s[1],
                    s[2],
                    s[3],
                    s[4],
                    s[5],
                    trace.von_mises
                )
                .unwrap();
                if i == n[0] / 2 {
                    values.push((k, j, trace.von_mises));
                }
            }
        }
    }
    std::fs::write(dir.join("cell-stresses.csv"), raw)?;
    let peak = values.iter().map(|v| v.2).fold(0., f64::max);
    let mut svg = String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"650\" viewBox=\"0 0 1000 650\"><rect width=\"1000\" height=\"650\" fill=\"white\"/><g font-family=\"sans-serif\" fill=\"#15212c\"><text x=\"50\" y=\"35\" font-size=\"22\">UNPROVED numerical FEA: agreed cantilever</text><text x=\"50\" y=\"65\">100 x 5 x 5 mm; E=200000 N/mm²; nu=0.3; fixed beam/side.x_min</text><text x=\"50\" y=\"90\">beam/side.x_max: traction (0,0,-0.04) N/mm²; force (0,0,-1) N</text>",
    );
    writeln!(svg,"<text x=\"50\" y=\"120\">Cell-centre von Mises at y={} mm; z/x scale 4/3; colour scale 0 to {peak:.5} N/mm²</text>",2.5/n[0] as f64).unwrap();
    for (k, j, vm) in values {
        let t = (vm / peak).clamp(0., 1.);
        let red = (255. * t) as u8;
        let blue = (255. * (1. - t)) as u8;
        writeln!(svg,"<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"rgb({red},50,{blue})\"><title>x-cell {k}, z-cell {j}: {vm:.12e} N/mm²</title></rect>",50.+900.*k as f64/n[2] as f64,145.+60.*(n[1]-1-j) as f64/n[1] as f64,900./n[2] as f64,60./n[1] as f64).unwrap();
    }
    svg.push_str("<text x=\"50\" y=\"235\">Centreline downward displacement, mm: field blue; Euler-Bernoulli grey</text><path d=\"M50 280V510H950\" fill=\"none\" stroke=\"black\"/>");
    for (theory, color) in [(true, "#777"), (false, "#175baf")] {
        svg.push_str(&format!(
            "<polyline fill=\"none\" stroke=\"{color}\" points=\""
        ));
        for station in 0..=80 {
            let x = station as f64 * 1.25;
            let w = if theory {
                x * x * (300. - x) / (6. * 200_000. * (625. / 12.))
            } else {
                -result.field.traces([2.5, 2.5, x])[0].displacement[1]
            };
            write!(svg, "{},{} ", 50. + 9. * x, 510. - w * 6500.).unwrap();
        }
        svg.push_str("\"/>");
    }
    writeln!(svg,"<text x=\"50\" y=\"550\">x: 0 to 100 mm; downward displacement: 0 to 0.035 mm</text><text x=\"50\" y=\"580\">Area-mean tip = {:.12e} mm; EB = -0.032 mm; error = {:.6}%</text><text x=\"50\" y=\"610\">Raw tensors and separate incident-cell traces: cell-stresses.csv and profiles.csv</text></g></svg>",result.tip,(result.tip/THEORY-1.).abs()*100.).unwrap();
    std::fs::write(dir.join("beam.svg"), svg)
}

fn accuracy(tip: f64) -> Result<(), Refusal> {
    let error = (tip / THEORY - 1.).abs();
    if !error.is_finite() || error > 0.02 {
        Err(fail(
            "fea.deflection-accuracy",
            &format!(
                "tip_mm={tip:.16e} reference_mm={THEORY} relative_error={error:.12e} allowed=0.02"
            ),
        ))
    } else {
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 2 {
        return Err("usage: m20_beam <output-directory> [cross-section-cells: 2|4|8]".into());
    }
    let selected: Option<usize> = args.get(1).map(|value| value.parse()).transpose()?;
    if selected.is_some_and(|cells| ![2, 4, 8].contains(&cells)) {
        return Err(fail("fea.grid", "supported cross-section grids: 2, 4, 8").into());
    }
    let recipe = Recipe::parse(RECIPE)?;
    let (reference, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = catallax_features::evaluate(&recipe, &mut dec)
            .map_err(|error| fail("fea.brep-reference", &error.to_string()))?;
        catallax_measure::volume(&built.model, built.result, &Default::default(), &mut dec)
    });
    let reference = reference?;
    assert!(!work.exhausted);
    assert_eq!(reference.value.lo, q(2500));
    assert_eq!(reference.value.hi, q(2500));
    println!(
        "UNPROVED elasticity; certified geometry: same recipe as independent B-rep volume=[2500,2500] mm3 reference_work={}",
        work.used
    );
    println!(
        "policy: Q1 full 2x2x2 Gauss; serial fixed-order Jacobi CG; residual<=1e-9; iterations<=10000; cells<=5120; axis_cells<=80; geometry_limit=100000000 per solve"
    );
    for n in GRIDS {
        if selected.is_some_and(|cells| cells != n[0]) {
            continue;
        }
        let mut expected = None;
        for workers in [1, 2, 4] {
            let answers = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers).map(|_| scope.spawn(move || run(n))).collect();
                handles
                    .into_iter()
                    .map(|h| h.join().unwrap())
                    .collect::<Vec<_>>()
            });
            for answer in answers {
                let answer = answer?;
                if let Some(prior) = &expected {
                    assert_eq!(&answer, prior);
                } else {
                    expected = Some(answer);
                }
            }
            println!(
                "grid_local={n:?} independent_workers={workers} complete_answer_identical=true"
            );
        }
        let answer = expected.unwrap();
        let result = &answer.estimate;
        println!("geometry_material_revision={}", answer.revision);
        let error = (result.tip / THEORY - 1.).abs();
        println!(
            "UNPROVED grid_world={:?} cells={} unknowns={} tip_mm={:.16e} signed_error_mm={:.16e} absolute_error_mm={:.16e} relative_error={error:.12e} iterations={} residual={:.12e} geometry_work={}",
            [n[2], n[0], n[1]],
            n.iter().product::<usize>(),
            3 * (n[0] + 1) * (n[1] + 1) * n[2],
            result.tip,
            result.tip - THEORY,
            (result.tip - THEORY).abs(),
            result.iterations,
            result.residual,
            answer.geometry_work
        );
        println!(
            "UNPROVED reaction_N={:?} reaction_moment_root_N_mm={:?}",
            world_vector(result.reaction),
            world_vector(result.reaction_moment)
        );
        let r = world_vector(result.reaction);
        let m = world_vector(result.reaction_moment);
        println!(
            "UNPROVED force_balance_N={:?} moment_balance_root_N_mm={:?}",
            [r[0], r[1], r[2] - 1.],
            [m[0], m[1] + 100., m[2]]
        );
        match accuracy(result.tip) {
            Ok(()) => println!("UNPROVED verdict=beam_deflection_mark_met"),
            Err(refusal) => println!(
                "UNPROVED refusal={} detail={} cells={} axis_limit=80",
                refusal.code,
                refusal.detail,
                n.iter().product::<usize>()
            ),
        }
        println!(
            "UNPROVED stiffness_entries={} matrix_vector_products={} scalar_matrix_products={}",
            result.matrix_entries,
            result.matvecs,
            result.matrix_entries * result.matvecs
        );
        artifacts(
            &Path::new(&args[0]).join(format!("{}-{}-{}", n[2], n[0], n[1])),
            n,
            result,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn accuracy_refusal_keeps_the_attained_numbers() {
        assert!(accuracy(THEORY).is_ok());
        for tip in [0., -0.03, f64::NAN] {
            let refusal = accuracy(tip).unwrap_err();
            assert_eq!(refusal.code, "fea.deflection-accuracy");
            assert!(refusal.detail.contains("relative_error="));
            assert!(refusal.detail.contains("allowed=0.02"));
        }
    }
    #[test]
    fn agreed_recipe_drives_geometry_load_and_solve() {
        let answer = run([2, 2, 20]).unwrap();
        assert!(answer.estimate.tip < 0.);
        let reaction = world_vector(answer.estimate.reaction);
        let moment = world_vector(answer.estimate.reaction_moment);
        assert!((reaction[2] - 1.).abs() < 1e-6);
        assert!((moment[1] + 100.).abs() < 1e-5);
        assert_eq!(
            world_tensor([1., 2., 3., 4., 5., 6.]),
            [3., 1., 2., 6., 4., 5.]
        );
    }
}
