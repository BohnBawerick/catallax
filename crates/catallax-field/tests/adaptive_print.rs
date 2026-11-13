//! Public analytic controls for the opt-in M19 route.
use catallax_field::{
    Location, Plan, Workspace,
    printing::{Options, adaptive},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
fn q(s: &str) -> Rational {
    s.parse().unwrap()
}
fn recipe(extra: &str, result: &str) -> Recipe {
    Recipe::parse(&format!(r#"{{"format":"catallax-recipe/3","name":"M19 public control","steps":[{{"id":"box","op":"block","size":["4","1","1"],"pose":{{"rotation":["1","0","0","0"],"translation":["-0.5","-0.5","-0.5"]}}}}{extra}],"result":"{result}"}}"#)).unwrap()
}
fn gyroid(kind: &str, level: &str) -> Recipe {
    recipe(
        &format!(
            r#",{{"id":"g","op":"gyroid","region":"box","period":"4","level":{level},"kind":"{kind}","frame":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}}}"#
        ),
        "g",
    )
}
fn options() -> Options {
    Options {
        tolerance: q("0.125"),
        max_cells: 100_000,
        max_vertices: 100_000,
    }
}
fn location(plan: &Plan, x: &str, y: &str, z: &str) -> Location {
    plan.proven_checked(
        &RPoint3::new(q(x), q(y), q(z)),
        &mut Workspace::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .location
}
#[test]
fn sheet_is_a_band_not_a_network_and_grades_retain_both_sides() {
    let band = Plan::compile(&gyroid("sheet", r#"{"constant":"1"}"#)).unwrap();
    for x in ["1", "3"] {
        assert_eq!(location(&band, x, "0", "0"), Location::Boundary);
    }
    for level in [
        r#"{"constant":"0.25"}"#,
        r#"{"linear":{"from":["0","0","-1"],"to":["0","0","1"],"start":"0.125","end":"0.25"}}"#,
    ] {
        let sheet = Plan::compile(&gyroid("sheet", level)).unwrap();
        let network = Plan::compile(&gyroid("network", level)).unwrap();
        assert_eq!(location(&sheet, "0", "0", "0"), Location::Inside);
        for x in ["-0.25", "0.25"] {
            assert_eq!(location(&sheet, x, "0", "0"), Location::Outside);
        }
        assert_eq!(location(&network, "-0.25", "0", "0"), Location::Inside);
        assert_eq!(location(&network, "0.25", "0", "0"), Location::Outside);
        assert_eq!(location(&sheet, "0", "0", "1"), Location::Outside);
    }
    for level in [
        r#"{"constant":"0"}"#,
        r#"{"constant":"-0.1"}"#,
        r#"{"linear":{"from":["0","0","0"],"to":["1","0","0"],"start":"0.1","end":"0"}}"#,
    ] {
        assert_eq!(
            Plan::compile(&gyroid("sheet", level)).unwrap_err().code,
            "gyroid.impossible"
        );
    }
}
#[test]
fn adaptive_cells_are_bound_to_immutable_complete_output_and_workers_agree() {
    let plan = Plan::compile(&recipe("", "box")).unwrap();
    let mut reference = None;
    for workers in [1, 2, 4, adaptive::default_workers()] {
        let (result, work) = with_operation_budget(20_000_000, || {
            adaptive::mesh_with_workers(&plan, &options(), workers, &mut Decider::default())
        });
        let output = result.unwrap();
        assert!(!work.exhausted);
        assert!(!output.cells().is_empty());
        let mut end = 0;
        for cell in output.cells() {
            assert_eq!(cell.triangles.start, end);
            assert!(cell.triangles.end > end);
            assert!(cell.max_deviation <= options().tolerance);
            for &t in &output.mesh().triangles[cell.triangles.clone()] {
                for v in t {
                    for a in 0..3 {
                        let p =
                            Rational::from_f64_exact(output.mesh().positions[v][a] as f64).unwrap();
                        assert!(p >= cell.bounds[a][0] && p <= cell.bounds[a][1]);
                    }
                }
            }
            end = cell.triangles.end;
        }
        assert_eq!(end, output.mesh().triangles.len());
        let widths: std::collections::BTreeSet<_> = output
            .cells()
            .iter()
            .flat_map(|cell| cell.bounds.iter().map(|b| &b[1] - &b[0]))
            .collect();
        assert!(
            widths.len() > 1,
            "the control must exercise unequal cell sizes"
        );
        if workers == 1 {
            let mut bytes = Vec::new();
            output.mesh().write_stl(&mut bytes).unwrap();
            let mesh =
                catallax_import::mesh::ClosedMesh::read(&bytes, &q("1"), &mut Decider::default())
                    .unwrap();
            assert_eq!(mesh.shells().len(), 1);
            assert!(mesh.shells()[0].outward);
            assert_eq!(
                mesh.winding(&[q("0"), q("0"), q("0")], &mut Decider::default())
                    .unwrap(),
                Some(1)
            );
            assert_eq!(
                mesh.winding(&[q("1"), q("1"), q("1")], &mut Decider::default())
                    .unwrap(),
                Some(0)
            );
        }
        let signature = (
            output.mesh().positions.clone(),
            output.mesh().triangles.clone(),
            output.mesh().max_deviation.clone(),
            output.mesh().cells,
            output.cells().to_vec(),
            work.used,
        );
        eprintln!(
            "cube workers={workers} cells={} facets={} work={}",
            output.mesh().cells,
            end,
            work.used
        );
        if let Some(old) = &reference {
            assert_eq!(&signature, old);
        } else {
            reference = Some(signature);
        }
    }
    let required = reference.unwrap().5;
    let mut stopped = None;
    for workers in [1, 2, 4, adaptive::default_workers()] {
        let (result, work) = with_operation_budget(required - 1, || {
            adaptive::mesh_with_workers(&plan, &options(), workers, &mut Decider::default())
        });
        let failure = result.unwrap_err();
        assert_eq!(failure.stage, adaptive::Stage::Representation);
        let signature = (*failure, work.used, work.exhausted);
        if let Some(old) = &stopped {
            assert_eq!(&signature, old);
        } else {
            stopped = Some(signature);
        }
    }
}
#[test]
fn refusal_keeps_original_request_evidence_and_separate_choice() {
    let plan = Plan::compile(&recipe("", "box")).unwrap();
    for limit in [1, 100, 1000] {
        let mut previous = None;
        for workers in [1, 2, 4, adaptive::default_workers()] {
            let (result, work) = with_operation_budget(limit, || {
                adaptive::mesh_with_workers(&plan, &options(), workers, &mut Decider::default())
            });
            let failure = result.unwrap_err();
            assert_eq!(failure.requested, options());
            let alternative = failure.coarser.as_ref().unwrap();
            assert!(alternative.label.starts_with("UNVALIDATED"));
            assert_eq!(
                alternative.options.tolerance,
                options().tolerance.mul_pow2(1)
            );
            assert_eq!(alternative.options.max_cells, options().max_cells);
            assert_eq!(failure.request_work, Some(work.used));
            let signature = (*failure, work.used, work.exhausted);
            // Refusal's Eq includes attached undecided evidence, unlike Debug.
            if let Some(old) = &previous {
                assert_eq!(&signature, old);
            } else {
                previous = Some(signature);
            }
        }
    }
    let tiny = Options {
        max_cells: 1,
        ..options()
    };
    let failure = adaptive::mesh(&plan, &tiny, &mut Decider::default()).unwrap_err();
    assert_eq!(failure.stage, adaptive::Stage::Preparation);
    assert_eq!(failure.cause.code, "field.print-budget");
    let vertices = Options {
        max_vertices: 4,
        ..options()
    };
    let failure = adaptive::mesh(&plan, &vertices, &mut Decider::default()).unwrap_err();
    assert_eq!(failure.stage, adaptive::Stage::Candidate);
    assert_eq!(failure.cause.code, "field.print-budget");
    assert_eq!(
        adaptive::mesh_with_workers(&plan, &options(), 0, &mut Decider::default())
            .unwrap_err()
            .cause
            .code,
        "field.print-options"
    );
    for workers in [1, 2, 4, adaptive::default_workers()] {
        let mut dec = Decider::default();
        dec.set_deadline(Some(std::time::Instant::now()));
        let failure =
            adaptive::mesh_with_workers(&plan, &options(), workers, &mut dec).unwrap_err();
        assert_eq!(failure.cause.code, "numbers.deadline");
        assert_eq!(failure.stage, adaptive::Stage::Preparation);
    }
}
#[test]
fn sheet_and_thin_cavity_never_publish_an_uncertified_candidate() {
    let cavity = recipe(
        r#",{"id":"void","op":"sphere","radius":"0.03125","pose":{"rotation":["1","0","0","0"],"translation":["0.13","0.11","0.09"]}},{"id":"part","op":"boolean","kind":"subtract","target":"box","tool":"void"}"#,
        "part",
    );
    for (name, recipe) in [
        ("sheet", gyroid("sheet", r#"{"constant":"0.25"}"#)),
        ("thin-cavity", cavity),
    ] {
        let plan = Plan::compile(&recipe).unwrap();
        let mut reference = None;
        for workers in [1, 2, 4] {
            let (result, work) = with_operation_budget(20_000_000, || {
                adaptive::mesh_with_workers(&plan, &options(), workers, &mut Decider::default())
            });
            let failure = result.unwrap_err();
            assert!(!work.exhausted);
            assert_eq!(failure.stage, adaptive::Stage::Refinement);
            assert_eq!(failure.cause.code, "field.print-unresolved");
            assert_eq!(failure.request_work, Some(work.used));
            assert_eq!(failure.requested, options());
            eprintln!(
                "{name} workers={workers}: stage={:?} refusal={} {} work={}",
                failure.stage, failure.cause.code, failure.cause.detail, work.used
            );
            let signature = (*failure, work.used, work.exhausted);
            if let Some(expected) = &reference {
                assert_eq!(expected, &signature);
            } else {
                reference = Some(signature);
            }
        }
    }
}

#[test]
fn full_p6_adaptive_request_retains_its_recorded_defaults_and_work() {
    let recipe = Recipe::parse(include_str!("../../../parts/P6-lattice-bracket.json")).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let request = Options::default();
    let mut reference = None;
    for workers in [1, 2, 4] {
        let (result, work) = with_operation_budget(100_000_000, || {
            adaptive::mesh_with_workers(&plan, &request, workers, &mut Decider::default())
        });
        let failure = result.unwrap_err();
        eprintln!(
            "P6-default workers={workers}: stage={:?} refusal={} {} work={} exhausted={}",
            failure.stage, failure.cause.code, failure.cause.detail, work.used, work.exhausted
        );
        assert_eq!(failure.requested, request);
        assert_eq!(failure.request_work, Some(work.used));
        let choice = failure.coarser.as_ref().unwrap();
        assert_eq!(choice.options.tolerance, q("0.02"));
        assert_eq!(choice.options.max_cells, request.max_cells);
        assert_eq!(choice.options.max_vertices, request.max_vertices);
        assert!(choice.label.starts_with("UNVALIDATED"));
        assert!(choice.estimate.contains("calibrated time unavailable"));
        let signature = (*failure, work.used, work.exhausted);
        if let Some(expected) = &reference {
            assert_eq!(expected, &signature);
        } else {
            reference = Some(signature);
        }
    }
}

#[test]
fn two_disjoint_blocks_reach_the_material_topology_guard() {
    let plan = Plan::compile(&recipe(
        r#",{"id":"other","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["5","-0.5","-0.5"]}},{"id":"part","op":"boolean","kind":"union","target":"box","tool":"other"}"#,
        "part",
    )).unwrap();
    let mut reference = None;
    for workers in [1, 2, 4] {
        let (result, work) = with_operation_budget(20_000_000, || {
            adaptive::mesh_with_workers(&plan, &options(), workers, &mut Decider::default())
        });
        let failure = result.unwrap_err();
        eprintln!(
            "disjoint-blocks workers={workers}: stage={:?} refusal={} {} work={}",
            failure.stage, failure.cause.code, failure.cause.detail, work.used
        );
        assert_eq!(failure.stage, adaptive::Stage::MaterialTopology);
        assert_eq!(failure.cause.code, "field.print-detail-topology");
        assert!(!work.exhausted);
        let signature = (*failure, work.used, work.exhausted);
        if let Some(expected) = &reference {
            assert_eq!(expected, &signature);
        } else {
            reference = Some(signature);
        }
    }
}

#[test]
fn distance_coverage_alone_can_hide_a_real_cavity() {
    let filled = Plan::compile(&recipe("", "box")).unwrap();
    let hollow = Plan::compile(&recipe(
        r#",{"id":"void","op":"sphere","radius":"0.03125","pose":{"rotation":["1","0","0","0"],"translation":["0.13","0.11","0.44"]}},{"id":"part","op":"boolean","kind":"subtract","target":"box","tool":"void"}"#,
        "part",
    )).unwrap();
    assert_eq!(location(&filled, "0.13", "0.11", "0.44"), Location::Inside);
    assert_eq!(location(&hollow, "0.13", "0.11", "0.44"), Location::Outside);
    for z in ["0.40875", "0.47125"] {
        assert_eq!(location(&hollow, "0.13", "0.11", z), Location::Boundary);
        assert_eq!(location(&filled, "0.13", "0.11", z), Location::Inside);
    }
    // The unchanged outer boundary covers every cavity point within this bound.
    let maximum_gap_to_top_face = &(&q("0.5") - &q("0.44")) + &q("0.03125");
    assert_eq!(maximum_gap_to_top_face, q("0.09125"));
    assert!(maximum_gap_to_top_face < options().tolerance);
    assert!(&q("0.44") + &q("0.03125") < q("0.5"));
}

#[test]
fn public_p6_and_network_caps_are_reported_without_coarsening() {
    let p6 = Recipe::parse(include_str!("../../../parts/P6-lattice-bracket.json")).unwrap();
    // Same public V7c region/period/level as the terminal network control.
    let network=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"V7c public control","steps":[{"id":"region","op":"block","size":["20","20","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-10","-10","-5"]}},{"id":"gyroid","op":"gyroid","region":"region","period":"10","level":{"constant":"0"},"kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"gyroid"}"#).unwrap();
    for (name, recipe) in [("P6", p6), ("V7c", network)] {
        let plan = Plan::compile(&recipe).unwrap();
        let request = Options {
            max_cells: 256,
            ..Options::default()
        };
        let mut reference = None;
        for workers in [1, 2, 4, adaptive::default_workers()] {
            let (result, work) = with_operation_budget(20_000_000, || {
                adaptive::mesh_with_workers(&plan, &request, workers, &mut Decider::default())
            });
            let failure = result.unwrap_err();
            assert_eq!(failure.requested.tolerance, q("0.01"));
            assert_eq!(failure.coarser.as_ref().unwrap().options.max_cells, 256);
            eprintln!(
                "{name} workers={workers} stage={:?} code={} work={}",
                failure.stage, failure.cause.code, work.used
            );
            assert_eq!(failure.request_work, Some(work.used));
            let signature = (*failure, work.used, work.exhausted);
            if let Some(old) = &reference {
                assert_eq!(&signature, old);
            } else {
                reference = Some(signature);
            }
        }
    }
}

#[test]
fn sheet_stage_builds_but_legacy_export_cannot_bypass_topology_refusal() {
    use catallax_field::lattice_stage::{GyroidSpec, Settings, Stage};
    use catallax_recipe::{Grade, GyroidKind, Number, Pose};
    let spec = GyroidSpec {
        period: Number::new(q("4")),
        level: Grade::Constant(Number::new(q("0.25"))),
        kind: GyroidKind::Sheet,
        frame: Pose::identity(),
    };
    let stage = Stage::new(
        &recipe("", "box"),
        spec,
        Settings {
            accuracy: q("0.125"),
            minimum_detail: q("0.25"),
            max_cells: 100_000,
            max_vertices: 100_000,
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        stage.export(&mut Decider::default()).unwrap_err().code,
        "field.print-detail-topology"
    );
    assert_eq!(
        catallax_field::printing::mesh_plan(stage.plan(), &options(), &mut Decider::default())
            .unwrap_err()
            .code,
        "field.print-detail-topology"
    );
}
