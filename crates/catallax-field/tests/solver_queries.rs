//! F7 exact integration, named boundary, material and parallel allowance controls.
use catallax_field::{Workspace, printing::BoxStorage, solver::*};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn material() -> Materials {
    Materials::Homogeneous {
        region: "beam".into(),
        material: "steel".into(),
    }
}
fn recipe() -> Recipe {
    Recipe::parse(include_str!("../examples/beam/recipe.json")).unwrap()
}
fn p(x: &str, y: &str, z: &str) -> RPoint3 {
    RPoint3::new(q(x), q(y), q(z))
}
#[test]
fn box_moments_slivers_and_material_interfaces() {
    let mut dec = Decider::default();
    let materials = Materials::Slabs {
        axis: 2,
        regions: vec![
            Region {
                id: "lower".into(),
                material: "A".into(),
                lo: q("0"),
                hi: q("50"),
            },
            Region {
                id: "upper".into(),
                material: "B".into(),
                lo: q("50"),
                hi: q("100"),
            },
        ],
    };
    let def = Definition::compile(&recipe(), materials, &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    let cell = queries
        .integration(
            &[[q("0"), q("10")], [q("0"), q("10")], [q("0"), q("100")]],
            &mut BoxStorage::default(),
            &mut dec,
        )
        .unwrap();
    assert_eq!(cell.fraction.lo, q("1"));
    assert_eq!(cell.fraction.hi, q("1"));
    let a = cell.moments("lower", &mut dec).unwrap();
    let b = cell.moments("upper", &mut dec).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let m = &a[9 * i + 3 * j + k];
                let n = &b[9 * i + 3 * j + k];
                assert!(m.is_exact() && n.is_exact());
                assert_eq!(
                    &m.lo + &n.lo,
                    q("10000")
                        .checked_div(&Rational::from(((i + 1) * (j + 1) * (k + 1)) as i64))
                        .unwrap()
                );
            }
        }
    }
    assert!(
        matches!(queries.material(&p("5","5","50"),&mut Workspace::default(),&mut dec).unwrap(),MaterialAt::Interface(v) if v.len()==2)
    );
    for denominator in [2, 16, 256, 4096] {
        let eps = Rational::from_ratio(1, denominator).unwrap();
        for shift in [q("0"), q("1/7")] {
            let b = [
                [&eps - &q("1"), eps.clone()],
                [q("1") + shift.clone(), q("2") + shift],
                [q("10"), q("11")],
            ];
            let c = queries
                .integration(&b, &mut BoxStorage::default(), &mut dec)
                .unwrap();
            assert_eq!(c.fraction.lo, eps);
            assert!(c.fraction.is_exact() && c.connected_box);
            let m = c.moments("lower", &mut dec).unwrap();
            assert!(m.iter().all(|m| m.is_exact()));
            assert_eq!(
                m[9].lo,
                &eps * &(&q("1") - &eps.checked_div(&q("2")).unwrap())
            );
        }
    }
}
#[test]
fn material_gap_refuses() {
    let m = Materials::Slabs {
        axis: 0,
        regions: vec![Region {
            id: "bad".into(),
            material: "A".into(),
            lo: q("1"),
            hi: q("10"),
        }],
    };
    assert!(Definition::compile(&recipe(), m, &mut Decider::default()).is_err());
}
#[test]
fn distance_ties_zero_and_named_loads() {
    let mut dec = Decider::default();
    let def = Definition::compile(&recipe(), material(), &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    let d = queries
        .box_distance(&p("13", "14", "50"), &mut dec)
        .unwrap();
    assert_eq!(d.squared, q("25"));
    assert!(d.distance(&q("0"), &mut dec).is_err());
    assert_eq!(d.displacements, vec![[q("-3"), q("-4"), q("0")]]);
    assert!(
        d.distance(&q("1/1000000"), &mut dec)
            .unwrap()
            .contains(&q("5"))
    );
    assert_eq!(
        queries
            .box_distance(&p("5", "5", "50"), &mut dec)
            .unwrap()
            .displacements
            .len(),
        4
    );
    assert_eq!(
        queries
            .box_distance(&p("5", "5", "0"), &mut dec)
            .unwrap()
            .squared,
        q("0")
    );
    for (role, z, condition) in [
        (Role::Start, "0", Condition::Fixed),
        (
            Role::End,
            "100",
            Condition::Traction([q("0"), q("-1/10"), q("0")]),
        ),
    ] {
        let name = Origin {
            step: StepId::new("beam").unwrap(),
            role,
        };
        let patch = queries
            .box_patch(&name, &p("5", "5", z), condition.clone(), &mut dec)
            .unwrap();
        assert_eq!(patch.area, q("100"));
        assert_eq!(patch.condition, condition);
        assert!(
            queries
                .box_patch(&name, &p("0", "0", z), Condition::Fixed, &mut dec)
                .is_err()
        );
    }
    let distances = queries.distances(&mut dec).unwrap();
    assert!(matches!(
        distances
            .nearest_face(&p("5", "5", "-1").to_real(), &Default::default(), &mut dec)
            .unwrap(),
        catallax_field::distance::Ownership::Proved(_)
    ));
}
#[test]
fn parallel_answers_and_both_exhaustion_paths_are_equal() {
    let mut dec = Decider::default();
    let def = Definition::compile(&recipe(), material(), &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    assert!(
        queries
            .grid(&GridOptions::default(), &mut dec)
            .answer
            .is_err()
    );
    for limit in [1_000_000, 20_000_000] {
        let mut previous = None;
        for workers in [Some(1), Some(2), Some(3), Some(4), None] {
            let (run, work) = with_operation_budget(limit, || {
                queries.grid(
                    &GridOptions {
                        workers,
                        ..Default::default()
                    },
                    &mut dec,
                )
            });
            assert!(run.answer.is_err() && work.exhausted);
            assert_eq!(run.spawned_workers, 0);
            let value = (run.answer, run.actual_work, work.used, work.exhausted);
            if let Some(old) = &previous {
                assert_eq!(old, &value);
            } else {
                previous = Some(value);
            }
        }
    }
    for (limit, work_per_cell) in [(100_000_000, 256), (1_000_000, 256), (100_000_000, 0)] {
        let mut prior = None;
        for workers in [Some(1), Some(2), Some(3), Some(4), None] {
            let (run, work) = with_operation_budget(limit, || {
                queries.grid(
                    &GridOptions {
                        divisions: 16,
                        workers,
                        work_per_cell,
                    },
                    &mut dec,
                )
            });
            if run.answer.is_ok() {
                assert_eq!(
                    run.spawned_workers,
                    workers
                        .unwrap_or_else(|| std::thread::available_parallelism().unwrap().get())
                        .min(32)
                );
            }
            let semantic = (run.answer, run.actual_work, work.used, work.exhausted);
            if let Some(p) = &prior {
                assert_eq!(&semantic, p);
            } else {
                prior = Some(semantic);
            }
        }
    }
    // A curved field incurs at least nine charged classifications in an unknown
    // cell, so one unit per cell deterministically exhausts a canonical job.
    let r=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere","steps":[{"id":"s","op":"sphere","radius":"1","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"s"}"#).unwrap();
    let def = Definition::compile(&r, material(), &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    let mut prior = None;
    for workers in [Some(1), Some(2), Some(3), Some(4), None] {
        let (run, work) = with_operation_budget(20_000_000, || {
            queries.grid(
                &GridOptions {
                    divisions: 8,
                    workers,
                    work_per_cell: 1,
                },
                &mut dec,
            )
        });
        assert!(run.answer.is_err());
        let v = (run.answer, run.actual_work, work.used, work.exhausted);
        if let Some(p) = &prior {
            assert_eq!(p, &v);
        } else {
            prior = Some(v);
        }
    }
}

#[test]
fn unknown_moments_cover_sphere_and_unsupported_fast_distance_refuses() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere","steps":[{"id":"s","op":"sphere","radius":"1","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"s"}"#).unwrap();
    let mut dec = Decider::default();
    let definition = Definition::compile(&recipe, material(), &mut dec).unwrap();
    let queries = definition.prepare(&mut dec).unwrap();
    assert!(queries.box_distance(&p("0", "0", "0"), &mut dec).is_err());
    let cell = queries
        .integration(queries.bounds(), &mut BoxStorage::default(), &mut dec)
        .unwrap();
    assert!(!cell.connected_box);
    let moments = cell.moments("beam", &mut dec).unwrap();
    // Unit ball volume is 4pi/3, and local xi_x=(x+1)/2.
    // Integral xi_x = V/2, integral xi_x² = (V+4pi/15)/4 = 2pi/5.
    let pi_lo = q("3.1415926535897932384626433832795028");
    let pi_hi = q("3.1415926535897932384626433832795029");
    for (index, factor) in [(0, q("4/3")), (9, q("2/3")), (18, q("2/5"))] {
        assert!(moments[index].lo <= &pi_lo * &factor && moments[index].hi >= &pi_hi * &factor);
    }
    let changed = Definition::compile(
        &recipe,
        Materials::Homogeneous {
            region: "beam".into(),
            material: "different".into(),
        },
        &mut dec,
    )
    .unwrap();
    assert_ne!(definition.revision(), changed.revision());
}

#[test]
fn material_partition_work_cannot_escape_reserved_jobs() {
    let mut dec = Decider::default();
    let materials = Materials::Slabs {
        axis: 2,
        regions: (0..300)
            .map(|i| Region {
                id: format!("layer-{i}"),
                material: "A".into(),
                lo: Rational::from_ratio(i, 3).unwrap(),
                hi: Rational::from_ratio(i + 1, 3).unwrap(),
            })
            .collect(),
    };
    let definition = Definition::compile(&recipe(), materials, &mut dec).unwrap();
    let queries = definition.prepare(&mut dec).unwrap();
    let mut previous = None;
    for workers in [Some(1), Some(2), Some(3), Some(4), None] {
        let (run, work) = with_operation_budget(100_000_000, || {
            queries.grid(
                &GridOptions {
                    divisions: 8,
                    workers,
                    ..Default::default()
                },
                &mut dec,
            )
        });
        assert!(run.answer.is_err());
        assert!(run.actual_work > 0 && !work.exhausted); // reserved child jobs exhausted, caller did not
        let value = (run.answer, run.actual_work, work.used, work.exhausted);
        if let Some(old) = &previous {
            assert_eq!(old, &value);
        } else {
            previous = Some(value);
        }
    }
}
