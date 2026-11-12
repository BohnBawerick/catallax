//! Bounded, unproved 3D Q1 elasticity experiment; see docs/FEA-FIRST-CUT.md.
#[path = "beam/solver.rs"]
mod solver;

use catallax_field::{Batch, FaceState, Location, Plan};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;

const RECIPE: &str = include_str!("beam/recipe.json");
const SIZE: [f64; 3] = [10., 10., 100.];
const E: f64 = 200_000.;
const NU: f64 = 0.3;
const FORCE: f64 = -10.;

fn refuse(code: &'static str, reason: &str) -> Refusal {
    Refusal::new(code, reason)
}
fn point(p: [f64; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|x| Rational::from_f64_exact(x).expect("finite grid coordinate"));
    RPoint3::new(x, y, z)
}

// This is deliberately a whitelist, not a claim that point samples certify an
// arbitrary shape. Recipe identity establishes the aligned box; field queries
// supply every occupied-cell decision and every boundary patch selection.
fn geometry(recipe: &Recipe, n: [usize; 3]) -> Result<(), Refusal> {
    if recipe != &Recipe::parse(RECIPE).expect("public fixture") {
        return Err(refuse(
            "fea.unsupported-geometry",
            "only the public aligned beam recipe is supported",
        ));
    }
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(recipe, &mut dec)?;
    let mut probes = Vec::new();
    let faces = [Role::Start, Role::End].map(|role| {
        plan.select(&Origin {
            step: StepId::new("beam").unwrap(),
            role,
        })
        .ok_or_else(|| refuse("fea.face", "required named end face is absent"))
    });
    let [fixed, loaded] = faces;
    let (fixed, loaded) = (fixed?, loaded?);
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let ijk = [i, j, k];
                let p = point(std::array::from_fn(|a| {
                    (ijk[a] as f64 + 0.5) * SIZE[a] / n[a] as f64
                }));
                probes.push((p, Location::Inside));
            }
        }
    }
    // Explicit void halo probes; no boundary point is silently treated as inside.
    for axis in 0..3 {
        for side in [-1., 1.] {
            let mut p = SIZE.map(|v| v / 2.);
            p[axis] = if side < 0. {
                -SIZE[axis] / n[axis] as f64 / 2.
            } else {
                SIZE[axis] + SIZE[axis] / n[axis] as f64 / 2.
            };
            probes.push((point(p), Location::Outside));
        }
    }
    let points = probes
        .iter()
        .map(|(point, _)| point.clone())
        .collect::<Vec<_>>();
    for (query, (_, expected)) in plan
        .proven_batch_checked(&points, &mut Batch::default())?
        .iter()
        .zip(probes)
    {
        if query.location != expected {
            return Err(refuse(
                "fea.cell",
                "beam cell or void halo classification disagrees with the supported domain",
            ));
        }
    }
    let mut face_points = Vec::new();
    for j in 0..n[1] {
        for i in 0..n[0] {
            for (id, z) in [(fixed, 0.), (loaded, SIZE[2])] {
                let p = point([
                    (i as f64 + 0.5) * SIZE[0] / n[0] as f64,
                    (j as f64 + 0.5) * SIZE[1] / n[1] as f64,
                    z,
                ]);
                face_points.push((id, p));
            }
        }
    }
    for query in plan.face_batch_checked(&face_points)? {
        if query.state != FaceState::OnFace {
            return Err(refuse(
                "fea.face",
                "end patch is not on the requested surviving named face",
            ));
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err("usage: beam_first_cut [cross-section-cells: 2|4|8] [recipe.json]".into());
    }
    let m: usize = args.first().map(|s| s.parse()).transpose()?.unwrap_or(4);
    if ![2, 4, 8].contains(&m) {
        return Err(refuse("fea.grid", "supported cross-section grids: 2, 4, 8").into());
    }
    let text = args
        .get(1)
        .map(std::fs::read_to_string)
        .transpose()?
        .unwrap_or_else(|| RECIPE.into());
    let recipe = Recipe::parse(&text)?;
    let n = [m, m, 10 * m];
    geometry(&recipe, n)?;
    let result = solver::solve(n, SIZE, E, NU, FORCE, 1e-9, 10_000)?;
    let inertia = SIZE[0] * SIZE[1].powi(3) / 12.;
    let theory = FORCE * SIZE[2].powi(3) / (3. * E * inertia);
    println!("UNPROVED ONLY: mm, N, N/mm^2; clamp beam/start; uniform -y traction on beam/end");
    println!(
        "UNPROVED grid={n:?} iterations={} relative_residual={:.12e}",
        result.iterations, result.residual
    );
    println!(
        "UNPROVED tip_area_mean_y={:.12e} Euler_Bernoulli={theory:.12e} relative_difference={:.8}",
        result.tip,
        (result.tip / theory - 1.).abs()
    );
    println!(
        "UNPROVED reaction_N={:?} reaction_moment_about_root_centre_N_mm={:?}",
        result.reaction, result.reaction_moment
    );
    for (z, y, vm) in result.stress {
        let bending = (FORCE * (SIZE[2] - z) * (y - SIZE[1] / 2.) / inertia).abs();
        println!("UNPROVED z={z:.8} y={y:.8} von_Mises={vm:.12e} beam_bending_abs={bending:.12e}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn geometry_and_refusals() {
        let recipe = Recipe::parse(RECIPE).unwrap();
        geometry(&recipe, [2, 2, 20]).unwrap();
        let other = Recipe::parse(&RECIPE.replace("100", "101")).unwrap();
        assert!(geometry(&other, [2, 2, 20]).is_err());
    }

    #[test]
    fn field_grid_queries_keep_their_work_at_every_worker_count() {
        use catallax_field::workers;
        use catallax_numbers::with_operation_budget;
        let run = || {
            let (result, work) = with_operation_budget(100_000, || {
                geometry(&Recipe::parse(RECIPE).unwrap(), [2, 2, 20])
            });
            (
                result.map_err(|r| (r.code, r.detail)),
                work.used,
                work.exhausted,
            )
        };
        let expected = workers::with_workers(1, run).unwrap();
        assert!(expected.0.is_ok());
        for count in [2, 4, 8, workers::worker_count()] {
            assert_eq!(workers::with_workers(count, run).unwrap(), expected);
        }
    }
}
