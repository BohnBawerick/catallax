//! Ordered field batches, mesh artifacts and deterministic request exhaustion.

use catallax_field::{Batch, Direction, Location, Plan, Query, Value, display, printing, workers};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;

const PARTS: [&str; 4] = [
    include_str!("../../../parts/P1-bracket.json"),
    include_str!("../../../parts/P2-collar.json"),
    include_str!("../../../parts/P3-filleted-boss.json"),
    include_str!("../../../parts/P4-shelled-tray.json"),
];

fn bits(q: Query) -> Vec<u64> {
    let loc = match q.location {
        Location::Inside => 0,
        Location::Outside => 1,
        Location::Boundary => 2,
        Location::Indeterminate => 3,
        _ => panic!("new classification needs a fingerprint"),
    };
    let value = q.value.enclosure();
    let mut bits = vec![
        loc,
        u64::from(matches!(q.value, Value::TrueDistance(_))),
        value.lo().to_bits(),
        value.hi().to_bits(),
        q.distance.millimetres.lo().to_bits(),
        q.distance.millimetres.hi().to_bits(),
    ];
    match q.direction {
        Direction::Undefined => bits.push(0),
        Direction::Ambiguous => bits.push(1),
        Direction::Unique(n) => {
            bits.push(2);
            for i in n {
                bits.extend([i.lo().to_bits(), i.hi().to_bits()]);
            }
        }
        _ => panic!("new direction needs a fingerprint"),
    }
    bits
}

fn queries(
    recipe: &'static str,
    count: usize,
    limit: u64,
) -> (Result<Vec<Vec<u64>>, String>, u64, bool) {
    std::thread::spawn(move || {
        workers::with_workers(count, || {
            let plan = Plan::compile(&Recipe::parse(recipe).unwrap()).unwrap();
            let points = (0..257)
                .map(|i| {
                    RPoint3::new(
                        Rational::from_i64(i % 17 - 8),
                        Rational::from_i64(i % 11 - 5),
                        Rational::from_i64(i % 7 - 3),
                    )
                })
                .collect::<Vec<_>>();
            let (result, work) = with_operation_budget(limit, || {
                plan.proven_batch_checked(&points, &mut Batch::default())
                    .map(|values| values.iter().copied().map(bits).collect())
                    .map_err(|error| format!("{}: {}", error.code, error.detail))
            });
            (result, work.used, work.exhausted)
        })
        .unwrap()
    })
    .join()
    .unwrap()
}

#[test]
fn proven_points_have_identical_bits_work_and_exhausted_prefixes() {
    for recipe in PARTS {
        let full = queries(recipe, 1, 100_000_000);
        assert!(full.0.is_ok());
        assert!(full.1 >= 257);
        for limit in [0, 1, full.1 / 2, full.1.saturating_sub(1), full.1] {
            let expected = queries(recipe, 1, limit);
            if limit == 0 {
                assert!(
                    expected
                        .0
                        .as_ref()
                        .unwrap_err()
                        .contains("numbers.operation-budget")
                );
                assert!(expected.2);
            }
            for count in [2, 4, 8, workers::worker_count()] {
                assert_eq!(
                    queries(recipe, count, limit),
                    expected,
                    "workers={count}, limit={limit}"
                );
            }
        }
    }
}

fn meshes(
    recipe: &'static str,
    count: usize,
    work_limit: u64,
    cell_limit: usize,
) -> Vec<(Result<Vec<u8>, String>, u64, bool)> {
    std::thread::spawn(move || {
        workers::with_workers(count, || {
            let mut rows = Vec::new();
            for print in [false, true] {
                let (result, work) = with_operation_budget(work_limit, || {
                    let mut dec = Decider::default();
                    let plan = Plan::compile_checked(&Recipe::parse(recipe).unwrap(), &mut dec)?;
                    let mut bytes = Vec::new();
                    if print {
                        let mesh = printing::mesh_plan(
                            &plan,
                            &printing::Options {
                                max_cells: cell_limit,
                                ..Default::default()
                            },
                            &mut dec,
                        )?;
                        mesh.write_stl(&mut bytes).unwrap();
                        mesh.write_obj(&mut bytes).unwrap();
                        bytes.extend_from_slice(&(mesh.cells as u64).to_le_bytes());
                    } else {
                        let mesh = display::mesh_plan(
                            &plan,
                            &display::Options {
                                max_cells: cell_limit,
                                ..Default::default()
                            },
                            &mut dec,
                        )?;
                        mesh.write_obj(&mut bytes).unwrap();
                        for point in mesh.positions {
                            for coordinate in point {
                                bytes.extend_from_slice(&coordinate.to_bits().to_le_bytes());
                            }
                        }
                        bytes.extend_from_slice(&(mesh.cells as u64).to_le_bytes());
                    }
                    Ok::<_, catallax_numbers::Refusal>(bytes)
                });
                rows.push((
                    result.map_err(|error| format!("{}: {}", error.code, error.detail)),
                    work.used,
                    work.exhausted,
                ));
            }
            rows
        })
        .unwrap()
    })
    .join()
    .unwrap()
}

#[test]
fn display_and_print_refuse_at_the_same_work_and_cell_limits() {
    for recipe in PARTS {
        for (work, cells) in [(0, 2_000_000), (100, 2_000_000), (100_000_000, 8)] {
            let expected = meshes(recipe, 1, work, cells);
            assert!(expected.iter().all(|row| row.0.is_err()));
            for count in [2, 4, 8, workers::worker_count()] {
                assert_eq!(meshes(recipe, count, work, cells), expected);
            }
        }
    }
}

#[test]
#[ignore = "one measured R4 round at the default display and print accuracy"]
fn target_display_and_print_bytes_and_work_match_at_every_count() {
    for (index, recipe) in PARTS.into_iter().enumerate() {
        let start = std::time::Instant::now();
        let expected = meshes(recipe, 1, 1_000_000, 2_000_000);
        eprintln!(
            "P{} workers=1 wall_ms={} rows={:?}",
            index + 1,
            start.elapsed().as_millis(),
            expected
                .iter()
                .map(|r| (r.0.as_ref().map(Vec::len), r.1, r.2))
                .collect::<Vec<_>>()
        );
        for count in [2, 4, 8, workers::worker_count()] {
            let start = std::time::Instant::now();
            let actual = meshes(recipe, count, 1_000_000, 2_000_000);
            eprintln!(
                "P{} workers={} wall_ms={} rows={:?}",
                index + 1,
                count,
                start.elapsed().as_millis(),
                actual
                    .iter()
                    .map(|r| (r.0.as_ref().map(Vec::len), r.1, r.2))
                    .collect::<Vec<_>>()
            );
            assert_eq!(actual, expected);
        }
    }
}
