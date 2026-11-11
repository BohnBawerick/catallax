//! F7 geometry/region/load adapter feeding the existing UNPROVED aligned Q1 solve.
#[path = "beam/solver.rs"]
mod elasticity;
use catallax_field::{
    printing::BoxStorage,
    solver::{Condition, Definition, Materials},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
fn q(n: i64) -> Rational {
    Rational::from(n)
}
fn main() {
    let (result, work) = with_operation_budget(100_000_000, || {
        let recipe = Recipe::parse(include_str!("beam/recipe.json")).unwrap();
        let mut dec = Decider::default();
        let definition = Definition::compile(
            &recipe,
            Materials::Homogeneous {
                region: "beam".into(),
                material: "isotropic-E200000-nu0.3".into(),
            },
            &mut dec,
        )
        .unwrap();
        let queries = definition.prepare(&mut dec).unwrap();
        assert!(
            matches!(queries.material(&RPoint3::new(q(5), q(5), q(50)), &mut catallax_field::Workspace::default(), &mut dec).unwrap(), catallax_field::solver::MaterialAt::Region {region, material} if region=="beam" && material=="isotropic-E200000-nu0.3")
        );
        let n = [4usize, 4, 40];
        let mut volume = q(0);
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let index = [i, j, k];
                    let b = std::array::from_fn(|a| {
                        let h = Rational::from_ratio([10, 10, 100][a], n[a] as i64).unwrap();
                        [&h * &q(index[a] as i64), &h * &q(index[a] as i64 + 1)]
                    });
                    let cell = queries
                        .integration(&b, &mut BoxStorage::default(), &mut dec)
                        .unwrap();
                    assert!(
                        cell.connected_box && cell.fraction.is_exact() && cell.fraction.lo == q(1)
                    );
                    let moments = cell.moments("beam", &mut dec).unwrap();
                    assert!(moments.iter().all(|m| m.is_exact()));
                    volume = volume + moments[0].lo.clone();
                }
            }
        }
        assert_eq!(volume, q(10000));
        let names = [Role::Start, Role::End].map(|role| Origin {
            step: StepId::new("beam").unwrap(),
            role,
        });
        let fixed = queries
            .box_patch(
                &names[0],
                &RPoint3::new(q(5), q(5), q(0)),
                Condition::Fixed,
                &mut dec,
            )
            .unwrap();
        let loaded = queries
            .box_patch(
                &names[1],
                &RPoint3::new(q(5), q(5), q(100)),
                Condition::Traction([q(0), Rational::from_ratio(-1, 10).unwrap(), q(0)]),
                &mut dec,
            )
            .unwrap();
        assert_eq!(fixed.area, q(100));
        let Condition::Traction(t) = loaded.condition else {
            unreachable!()
        };
        let force = &loaded.area * &t[1];
        assert_eq!(force, q(-10));
        // Only this aligned, full-cell fixture is sent to the unchanged example
        // solver. No cut-cell alpha-scaled stiffness or material mixing occurs.
        elasticity::solve(
            n,
            [10., 10., 100.],
            200_000.,
            0.3,
            force.to_f64_nearest(),
            1e-9,
            10_000,
        )
        .unwrap()
    });
    assert!(!work.exhausted);
    let theory = -0.02;
    println!(
        "UNPROVED F7 beam: tip_mm={:.12e} Euler_Bernoulli_mm={theory:.12e} relative_difference={:.8} iterations={} residual={:.12e} geometry_work={}",
        result.tip,
        (result.tip / theory - 1.).abs(),
        result.iterations,
        result.residual,
        work.used
    );
    println!(
        "UNPROVED reactions_N={:?} reaction_moment_N_mm={:?}",
        result.reaction, result.reaction_moment
    );
    for (z, y, vm) in result.stress {
        let bending = (10. * (100. - z) * (y - 5.) / (10000. / 12.)).abs();
        println!("UNPROVED z={z:.8} y={y:.8} von_Mises={vm:.12e} beam_bending_abs={bending:.12e}");
    }
}
