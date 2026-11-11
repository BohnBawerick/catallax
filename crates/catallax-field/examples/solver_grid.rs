//! Public P1/P4 64-cubed enclosure and complete semantic worker comparison.
use catallax_field::solver::{Definition, GridOptions, Materials};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
fn main() {
    for name in ["P1-bracket", "P4-shelled-tray"] {
        let text = std::fs::read_to_string(format!("parts/{name}.json")).unwrap();
        let recipe = Recipe::parse(&text).unwrap();
        let (definition, preparation) = with_operation_budget(100_000_000, || {
            Definition::compile(
                &recipe,
                Materials::Homogeneous {
                    region: "body".into(),
                    material: "explicit-control-material".into(),
                },
                &mut Decider::default(),
            )
        });
        assert!(!preparation.exhausted);
        let definition = definition.unwrap();
        let (queries, prepared) =
            with_operation_budget(100_000_000, || definition.prepare(&mut Decider::default()));
        assert!(!prepared.exhausted);
        let queries = queries.unwrap();
        let (reference, brep_work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
            catallax_measure::volume(&built.model, built.result, &Default::default(), &mut dec)
                .unwrap()
                .value
        });
        assert!(!brep_work.exhausted);
        println!(
            "{name} compile_work={} prepare_work={} independent_brep_work={} brep_mm3=[{},{}]",
            preparation.used, prepared.used, brep_work.used, reference.lo, reference.hi
        );
        let (distance_probe, probe_work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let d = queries.distances(&mut dec)?;
            let b = queries.bounds();
            let p = catallax_shapes::RPoint3::new(
                &b[0][1] + &Rational::from(1),
                (&b[1][0] + &b[1][1])
                    .checked_div(&Rational::from(2))
                    .unwrap(),
                (&b[2][0] + &b[2][1])
                    .checked_div(&Rational::from(2))
                    .unwrap(),
            )
            .to_real();
            let distance = d.point_to_boundary(&p, &Default::default(), &mut dec)?;
            let owner = d.nearest_face(&p, &Default::default(), &mut dec)?;
            Ok::<_, catallax_numbers::Refusal>((distance, owner))
        });
        if probe_work.exhausted {
            let error = distance_probe.err().unwrap_or_else(|| {
                catallax_numbers::Refusal::new(
                    "solver.probe-budget",
                    "probe exhausted its operation; any partial answer discarded",
                )
            });
            println!(
                "{name} trimmed_distance_probe_refusal={} detail={} work={} exhausted=true",
                error.code, error.detail, probe_work.used
            );
        } else {
            match distance_probe {
                Ok((d, owner)) => {
                    println!(
                        "{name} trimmed_distance_mm=[{}, {:?}] converged={} probe_work={} exhausted=false",
                        d.lower, d.upper, d.converged, probe_work.used
                    );
                    match owner {
                        catallax_field::distance::Ownership::Proved(o) => println!(
                            "{name} owner={:?} owner_distance_mm=[{},{}] witness_direction={:?} margin_lower_mm={:?}",
                            o.origin, o.lower, o.upper, o.direction, o.margin
                        ),
                        catallax_field::distance::Ownership::Refused(a) => println!(
                            "{name} ownership_refusal={} contenders={}",
                            a.refusal().code,
                            a.contenders.len()
                        ),
                    }
                }
                Err(e) => println!(
                    "{name} trimmed_distance_probe_refusal={} detail={} work={} exhausted=false",
                    e.code, e.detail, probe_work.used
                ),
            }
        }
        let mut prior = None;
        for workers in [Some(1), Some(2), Some(4), None] {
            let now = std::time::Instant::now();
            let (run, work) = with_operation_budget(100_000_000, || {
                queries.grid(
                    &GridOptions {
                        workers,
                        ..Default::default()
                    },
                    &mut Decider::default(),
                )
            });
            println!(
                "{name} requested_workers={workers:?} spawned={} actual_work={} charged={} exhausted={} elapsed_s={:.3}",
                run.spawned_workers,
                run.actual_work,
                work.used,
                work.exhausted,
                now.elapsed().as_secs_f64()
            );
            assert!(!work.exhausted);
            let grid = run.answer.unwrap();
            assert!(grid.volume.lo <= reference.lo && reference.hi <= grid.volume.hi);
            println!(
                "{name} volume_mm3=[{},{}] width_mm3={} width_fraction={:.8} cells={}",
                grid.volume.lo,
                grid.volume.hi,
                grid.volume.width(),
                grid.volume.width().to_f64_nearest() / reference.lo.to_f64_nearest(),
                grid.cells.len()
            );
            if let Some(prior) = &prior {
                assert_eq!(&grid, prior);
            } else {
                let (checked, moment_work) = with_operation_budget(20_000_000, || {
                    let mut checked = 0;
                    let mut dec = Decider::default();
                    for index in (0..grid.cells.len()).step_by(4093) {
                        let b = grid.cell_bounds(index).unwrap();
                        let cell = queries
                            .integration(
                                &b,
                                &mut catallax_field::printing::BoxStorage::default(),
                                &mut dec,
                            )
                            .unwrap();
                        assert_eq!(cell.state, grid.cells[index].state);
                        assert_eq!(cell.fraction, grid.cells[index].fraction);
                        if cell.pieces.is_empty() {
                            continue;
                        }
                        let moments = cell.moments("body", &mut dec).unwrap();
                        let v = b.iter().fold(Rational::from(1), |v, [a, b]| v * (b - a));
                        assert_eq!(moments[0].lo, &cell.fraction.lo * &v);
                        assert_eq!(moments[0].hi, &cell.fraction.hi * &v);
                        assert!(moments.iter().all(|m| m.lo <= m.hi));
                        checked += 1;
                    }
                    checked
                });
                assert!(!moment_work.exhausted);
                println!(
                    "{name} attributed_moment_cells={checked} moment_work={}",
                    moment_work.used
                );
                prior = Some(grid);
            }
        }
    }
}
