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
        for workers in [Some(1), Some(2), Some(3), Some(4), Some(8), None] {
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
        for workers in [Some(1), Some(2), Some(3), Some(4), Some(8), None] {
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
                        .unwrap_or_else(catallax_field::workers::worker_count)
                        .min(catallax_field::workers::worker_count())
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
    for workers in [Some(1), Some(2), Some(3), Some(4), Some(8), None] {
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
    for workers in [Some(1), Some(2), Some(3), Some(4), Some(8), None] {
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

/// P1 witness: a sheet gyroid band inside a native sphere has a regular level-set
/// boundary that no named face registers. At (1,0,0), g = 1 and dg/dy = pi/2, so
/// the point lies on the sheet boundary strictly inside the radius-2 sphere. The
/// general whole-boundary distance must refuse rather than report the named
/// sphere-face distance (about 1 mm) as the distance to the whole boundary.
#[test]
fn sheet_gyroid_boundary_refuses_named_face_only_distance() {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere sheet","steps":[{"id":"ball","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"sheet","op":"gyroid","region":"ball","period":"4","level":{"constant":"1"},"kind":"sheet","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"sheet"}"#).unwrap();
    let mut dec = Decider::default();
    let definition = Definition::compile(&recipe, material(), &mut dec).unwrap();
    let queries = definition.prepare(&mut dec).unwrap();
    let witness = p("1", "0", "0");
    let at = queries
        .membership(&witness, &mut Workspace::default(), &mut dec)
        .unwrap();
    assert_eq!(at.location, catallax_field::Location::Boundary);
    let Err(refusal) = queries.distances(&mut dec) else {
        panic!("whole-boundary distance must refuse an unnamed sheet boundary");
    };
    assert_eq!(refusal.code, "solver.distance-coverage");
}

#[test]
fn p1_and_p4_distance_probes_and_general_boundary_distance_converge() {
    for (name, text, expect_owner) in [
        (
            "P1-bracket",
            include_str!("../../../parts/P1-bracket.json"),
            false,
        ),
        (
            "P4-shelled-tray",
            include_str!("../../../parts/P4-shelled-tray.json"),
            true,
        ),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        let mut dec = Decider::default();
        let def = Definition::compile(&recipe, material(), &mut dec).unwrap();
        let queries = def.prepare(&mut dec).unwrap();
        let mut prior_work = None;
        for workers in [1usize, 2, 4] {
            let (probe, work) = catallax_field::workers::with_workers(workers, || {
                with_operation_budget(100_000_000, || {
                    let mut dec = Decider::default();
                    let d = queries.distances(&mut dec)?;
                    let b = queries.bounds();
                    let probe_pt = RPoint3::new(
                        &b[0][1] + &Rational::from(1),
                        (&b[1][0] + &b[1][1])
                            .checked_div(&Rational::from(2))
                            .unwrap(),
                        (&b[2][0] + &b[2][1])
                            .checked_div(&Rational::from(2))
                            .unwrap(),
                    );
                    let real_pt = probe_pt.clone().to_real();
                    let distance = d.point_to_boundary(&real_pt, &Default::default(), &mut dec)?;
                    let owner = d.nearest_face(&real_pt, &Default::default(), &mut dec)?;
                    let q4 = queries.boundary_distance_with(
                        &d,
                        &probe_pt,
                        &Default::default(),
                        &mut Workspace::default(),
                        &mut dec,
                    )?;
                    Ok::<_, catallax_numbers::Refusal>((distance, owner, q4))
                })
            })
            .unwrap();
            let (distance, _owner, q4) = probe.unwrap();
            assert!(!work.exhausted, "{name} exhausted 100M");
            assert!(distance.converged, "{name} did not converge: {distance:?}");
            assert!(q4.converged, "{name} Q4 did not converge: {q4:?}");
            assert_eq!(q4.location, Some(catallax_field::Location::Outside));
            if expect_owner {
                assert!(matches!(
                    &q4.direction,
                    BoundaryDirection::UniqueNearest {
                        owner,
                        margin: Some(m),
                        ..
                    } if owner.step.as_str() == "blank" && *m > q("0")
                ));
            } else {
                assert!(matches!(
                    &q4.direction,
                    BoundaryDirection::WitnessOnly {
                        contenders,
                        ..
                    } if contenders.len() == 2
                ));
            }
            if let Some(prev) = prior_work {
                assert_eq!(prev, work.used, "{name} work differed at {workers} workers");
            } else {
                prior_work = Some(work.used);
            }
        }
    }
    // Also check boundary_distance at a smooth boundary point and at a box edge.
    let mut dec = Decider::default();
    let def = Definition::compile(&recipe(), material(), &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    let on_face = queries
        .boundary_distance(
            &p("5", "5", "0"),
            &Default::default(),
            &mut Workspace::default(),
            &mut dec,
        )
        .unwrap();
    assert_eq!(on_face.lower, q("0"));
    assert_eq!(on_face.upper, Some(q("0")));
    assert!(matches!(
        on_face.direction,
        BoundaryDirection::BoundaryNormal(catallax_field::Direction::Unique(_))
    ));
    let on_edge = queries
        .boundary_distance(
            &p("0", "0", "0"),
            &Default::default(),
            &mut Workspace::default(),
            &mut dec,
        )
        .unwrap();
    assert_eq!(on_edge.lower, q("0"));
    assert_eq!(
        on_edge.direction,
        BoundaryDirection::BoundaryNormal(catallax_field::Direction::Ambiguous)
    );
}

#[test]
fn named_face_loads_and_supports_integrate_planar_holed_slotted_and_curved_faces() {
    let mut dec = Decider::default();
    let width = q("1/1000000");
    let pi_lo = q("3.1415926535897932384626433832795028");
    let pi_hi = q("3.1415926535897932384626433832795029");

    // 1. P1 bracket: holed plate base, slotted flange wall, cylindrical boss/bore/slot-end, and buried flange/start.
    let p1_recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let p1_def = Definition::compile(&p1_recipe, material(), &mut dec).unwrap();
    let p1_q = p1_def.prepare(&mut dec).unwrap();
    let p1_dist = p1_q.distances(&mut dec).unwrap();

    // plate/start: 60x40 minus r=4 bore and four r=2 bolt holes => 2400 - 32*pi.
    let plate_start = Origin {
        step: StepId::new("plate").unwrap(),
        role: Role::Start,
    };
    let patch_plate = p1_q
        .face_patch_with(
            &p1_dist,
            &plate_start,
            Condition::Pressure(q("2")),
            &width,
            &mut dec,
        )
        .unwrap();
    let expected_plate_lo = &q("2400") - &(&q("32") * &pi_hi);
    let expected_plate_hi = &q("2400") - &(&q("32") * &pi_lo);
    assert!(
        patch_plate.area.lo <= expected_plate_hi && patch_plate.area.hi >= expected_plate_lo,
        "plate/start area: {:?}",
        patch_plate.area
    );
    // Outward normal on plate/start is (0, 0, -1), so inward pressure +2 pushes in (0, 0, +2*area).
    assert!(
        patch_plate.resultant_force[2].lo <= &expected_plate_hi * &q("2")
            && patch_plate.resultant_force[2].hi >= &expected_plate_lo * &q("2")
    );

    // flange/side.top (y = 5): 60x25 = 1500 minus stadium slot (14x6 + pi*3^2 = 84 + 9*pi) => 1416 - 9*pi.
    let flange_front = Origin {
        step: StepId::new("flange").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("top").unwrap(),
            path: None,
        },
    };
    let patch_flange = p1_q
        .face_patch_with(
            &p1_dist,
            &flange_front,
            Condition::TotalForce([q("0"), q("-100"), q("0")]),
            &width,
            &mut dec,
        )
        .unwrap();
    let expected_flange_lo = &q("1416") - &(&q("9") * &pi_hi);
    let expected_flange_hi = &q("1416") - &(&q("9") * &pi_lo);
    assert!(
        patch_flange.area.lo <= expected_flange_hi && patch_flange.area.hi >= expected_flange_lo,
        "flange/top area: {:?}",
        patch_flange.area
    );
    assert!(patch_flange.resultant_force[1].contains(&q("-100")));

    // boss/side.wall: cylinder r=8, z=5..15 (h=10) => area = 160*pi, centroid = (30, 25, 10).
    let boss_wall = Origin {
        step: StepId::new("boss").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("wall").unwrap(),
            path: None,
        },
    };
    let patch_boss = p1_q
        .face_patch_with(
            &p1_dist,
            &boss_wall,
            Condition::Traction([q("0"), q("0"), q("-1")]),
            &width,
            &mut dec,
        )
        .unwrap();
    assert_eq!(patch_boss.pieces, 1);
    assert!(
        patch_boss.area.lo <= &q("160") * &pi_hi && patch_boss.area.hi >= &q("160") * &pi_lo,
        "boss/wall area: {:?}",
        patch_boss.area
    );
    assert!(patch_boss.centroid[0].contains(&q("30")));
    assert!(patch_boss.centroid[1].contains(&q("25")));
    assert!(patch_boss.centroid[2].contains(&q("10")));
    // Full cylinder normal integral is zero.
    for k in 0..3 {
        assert!(patch_boss.normal_integral[k].contains(&q("0")));
    }

    // slot-tool-0/side.wall: half-cylinder r=3, y=0..5 (h=5), dtheta=pi => area = 15*pi, normal integral along +x is 2*r*h = 30.
    let slot_left = Origin {
        step: StepId::new("slot-tool-0").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("wall").unwrap(),
            path: None,
        },
    };
    let patch_slot = p1_q
        .face_patch_with(
            &p1_dist,
            &slot_left,
            Condition::Pressure(q("1")),
            &width,
            &mut dec,
        )
        .unwrap();
    assert_eq!(patch_slot.pieces, 1);
    assert!(
        patch_slot.area.lo <= &q("15") * &pi_hi && patch_slot.area.hi >= &q("15") * &pi_lo,
        "slot-tool-0/wall area: {:?}",
        patch_slot.area
    );
    // Outward normal of the subtracted left semicircle points into the slot (+x), so inward pressure pushes in -x with force -30 N.
    assert!(patch_slot.normal_integral[0].contains(&q("30")));
    assert!(patch_slot.resultant_force[0].contains(&q("-30")));

    // Buried face flange/start refuses with solver.patch-membership.
    let flange_start = Origin {
        step: StepId::new("flange").unwrap(),
        role: Role::Start,
    };
    let buried_err = p1_q
        .face_patch_with(&p1_dist, &flange_start, Condition::Fixed, &width, &mut dec)
        .unwrap_err();
    assert_eq!(buried_err.code, "solver.patch-membership");

    // 2. P4 shelled tray: rounding/side.wall has 4 disconnected curved corner patches.
    let p4_recipe = Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap();
    let p4_def = Definition::compile(&p4_recipe, material(), &mut dec).unwrap();
    let p4_q = p4_def.prepare(&mut dec).unwrap();
    let p4_dist = p4_q.distances(&mut dec).unwrap();
    let rounding_wall = Origin {
        step: StepId::new("rounding").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("wall").unwrap(),
            path: None,
        },
    };
    let patch_rounding = p4_q
        .face_patch_with(&p4_dist, &rounding_wall, Condition::Fixed, &width, &mut dec)
        .unwrap();
    assert_eq!(patch_rounding.pieces, 4);
    assert!(patch_rounding.centroid[0].contains(&q("0")));
    assert!(patch_rounding.centroid[1].contains(&q("0")));
    assert!(patch_rounding.centroid[2].contains(&q("6")));

    // blank/side.x_max on P4: trimmed to y in [-sqrt(129), sqrt(129)], z in [0, 12], area = 24*sqrt(129).
    let blank_xmax = Origin {
        step: StepId::new("blank").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("x_max").unwrap(),
            path: None,
        },
    };
    let patch_xmax = p4_q
        .face_patch_with(&p4_dist, &blank_xmax, Condition::Fixed, &width, &mut dec)
        .unwrap();
    let expected_xmax = dec
        .enclose(
            &(catallax_numbers::Real::from_i64(24) * catallax_numbers::Real::from_i64(129).sqrt()),
            &q("1/1000000000"),
        )
        .unwrap();
    assert!(
        patch_xmax.area.lo <= expected_xmax.hi() && patch_xmax.area.hi >= expected_xmax.lo(),
        "blank/x_max area: {:?}",
        patch_xmax.area
    );
}

#[test]
fn recipe_material_regions_support_non_slab_and_curved_partitions() {
    let mut dec = Decider::default();
    // 1. Curved cylindrical core ("copper") inside a 10x10x10 block with outer shell ("steel").
    let curved_recipe = Recipe::parse(
        r#"{
  "format": "catallax-recipe/4",
  "name": "cored block",
  "steps": [
    {
      "id": "whole",
      "op": "block",
      "size": ["10", "10", "10"],
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["-5", "-5", "0"] }
    },
    {
      "id": "core",
      "op": "cylinder",
      "radius": "2",
      "height": "10",
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }
    },
    {
      "id": "shell",
      "op": "boolean",
      "kind": "subtract",
      "target": "whole",
      "tool": "core"
    },
    {
      "id": "part",
      "op": "boolean",
      "kind": "union",
      "target": "shell",
      "tool": "core"
    }
  ],
  "result": "part",
  "materials": [
    { "id": "core", "region": "core", "material": "copper" },
    { "id": "shell", "region": "shell", "material": "steel" }
  ]
}"#,
    )
    .unwrap();
    let def = Definition::from_recipe(&curved_recipe, &mut dec).unwrap();
    let queries = def.prepare(&mut dec).unwrap();
    let mut ws = Workspace::default();
    assert_eq!(
        queries
            .material(&p("0", "0", "5"), &mut ws, &mut dec)
            .unwrap(),
        MaterialAt::Region {
            region: "core".into(),
            material: "copper".into(),
        }
    );
    assert_eq!(
        queries
            .material(&p("4", "0", "5"), &mut ws, &mut dec)
            .unwrap(),
        MaterialAt::Region {
            region: "shell".into(),
            material: "steel".into(),
        }
    );
    assert_eq!(
        queries
            .material(&p("2", "0", "5"), &mut ws, &mut dec)
            .unwrap(),
        MaterialAt::Interface(vec![
            ("core".into(), "copper".into()),
            ("shell".into(), "steel".into()),
        ])
    );
    assert_eq!(
        queries
            .material(&p("5", "0", "5"), &mut ws, &mut dec)
            .unwrap(),
        MaterialAt::Boundary
    );
    assert_eq!(
        queries
            .material(&p("6", "0", "5"), &mut ws, &mut dec)
            .unwrap(),
        MaterialAt::Void
    );

    // Verify cell integration on a cell straddling the curved interface keeps whole-cell fraction == [1, 1]
    // while enclosing both regions' moments.
    let interface_cell = queries
        .integration(
            &[[q("-3"), q("3")], [q("-3"), q("3")], [q("0"), q("10")]],
            &mut BoxStorage::default(),
            &mut dec,
        )
        .unwrap();
    assert_eq!(interface_cell.fraction.lo, q("1"));
    assert_eq!(interface_cell.fraction.hi, q("1"));
    let core_moments = interface_cell.moments("core", &mut dec).unwrap();
    let shell_moments = interface_cell.moments("shell", &mut dec).unwrap();
    let pi_lo = q("3.1415926535897932384626433832795028");
    let pi_hi = q("3.1415926535897932384626433832795029");
    // True core volume is pi * 2^2 * 10 = 40 * pi; shell volume in [-3,3]^2x[0,10] is 360 - 40*pi.
    assert!(
        core_moments[0].lo <= &q("40") * &pi_lo && core_moments[0].hi >= &q("40") * &pi_hi,
        "core M000: {:?}",
        core_moments[0]
    );
    assert!(
        shell_moments[0].lo <= &q("360") - &(&q("40") * &pi_hi)
            && shell_moments[0].hi >= &q("360") - &(&q("40") * &pi_lo),
        "shell M000: {:?}",
        shell_moments[0]
    );

    // Parallel grid determinism at 1, 2, 4 workers.
    let mut prior = None;
    for workers in [Some(1), Some(2), Some(4)] {
        let (run, work) = with_operation_budget(20_000_000, || {
            queries.grid(
                &GridOptions {
                    divisions: 4,
                    workers,
                    work_per_cell: 512,
                },
                &mut dec,
            )
        });
        let grid = run.answer.unwrap();
        assert_eq!(grid.volume.lo, q("1000"));
        assert_eq!(grid.volume.hi, q("1000"));
        let record = (grid, run.actual_work, work.used, work.exhausted);
        if let Some(prev) = &prior {
            assert_eq!(prev, &record);
        } else {
            prior = Some(record);
        }
    }

    // Overlapping recipe material regions refuse at compile.
    let overlap_recipe = Recipe::parse(
        r#"{
  "format": "catallax-recipe/4",
  "name": "overlap",
  "steps": [
    {
      "id": "whole",
      "op": "block",
      "size": ["10", "10", "10"],
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["-5", "-5", "0"] }
    },
    {
      "id": "core",
      "op": "cylinder",
      "radius": "2",
      "height": "10",
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }
    }
  ],
  "result": "whole",
  "materials": [
    { "id": "core", "region": "core", "material": "copper" },
    { "id": "all", "region": "whole", "material": "steel" }
  ]
}"#,
    )
    .unwrap();
    assert!(Definition::from_recipe(&overlap_recipe, &mut dec).is_err());

    // Gapped recipe material regions refuse at compile.
    let gap_recipe = Recipe::parse(
        r#"{
  "format": "catallax-recipe/4",
  "name": "gap",
  "steps": [
    {
      "id": "whole",
      "op": "block",
      "size": ["10", "10", "10"],
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["-5", "-5", "0"] }
    },
    {
      "id": "core",
      "op": "cylinder",
      "radius": "2",
      "height": "10",
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }
    }
  ],
  "result": "whole",
  "materials": [
    { "id": "core", "region": "core", "material": "copper" }
  ]
}"#,
    )
    .unwrap();
    assert!(Definition::from_recipe(&gap_recipe, &mut dec).is_err());
}
