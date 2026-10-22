use super::*;
use catallax_features::evaluate;
use catallax_numbers::Decider;
use serde_json::{Value, json};
fn pose(x: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,"0","0"]})
}
fn recipe(steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(&json!({"format":catallax_recipe::FORMAT,"name":"comparison","steps":steps,"result":result}).to_string()).unwrap()
}
fn block(id: &str, size: &str, x: &str) -> Value {
    json!({"id":id,"op":"block","size":[size,size,size],"pose":pose(x)})
}
fn options() -> Options {
    Options {
        divisions: 2,
        ..Default::default()
    }
}
fn check(recipe: &Recipe) -> Report {
    check_with_options(recipe, &options())
}
fn check_with_options(recipe: &Recipe, options: &Options) -> Report {
    let controller = Controller::default();
    assert_eq!(
        controller.begin().check(recipe, options),
        Completion::Published
    );
    controller.latest().unwrap()
}
fn compare(brep: &Recipe, field: &Recipe, options: &Options) -> Report {
    let built = evaluate(brep, &mut Decider::default()).unwrap();
    let controller = Controller::default();
    assert_eq!(
        controller.begin().check_built(&built, field, options),
        Completion::Published
    );
    controller.latest().unwrap()
}
#[test]
fn agreeing_stock_across_poses_and_scales() {
    // Periodic chart copies count toward the mesher's work allowance. Keep the
    // same display accuracy at every scale, with enough work for the largest case.
    let mut controls = options();
    controls.mesh.max_vertices = 4_000_000;
    for size in ["1/10", "1", "10"] {
        for rotation in [["1", "0", "0", "0"], ["1", "1", "2", "3"]] {
            for op in ["block", "cylinder", "sphere", "cone", "torus"] {
                eprintln!("checking {op} size={size} rotation={rotation:?}");
                let step = match op {
                    "block" => json!({"id":"stock","op":op,"size":[size,size,size]}),
                    "cylinder" => json!({"id":"stock","op":op,"radius":size,"height":size}),
                    "sphere" => json!({"id":"stock","op":op,"radius":size}),
                    "cone" => {
                        json!({"id":"stock","op":op,"base_radius":size,"top_radius":"0","height":size})
                    }
                    _ => json!({"id":"stock","op":op,"major":"2","minor":size}),
                };
                // Ring tori require major > minor at every scale.
                let mut step = step;
                if op == "torus" {
                    step["major"] =
                        json!((Rational::parse(size).unwrap() * Rational::from_i64(3)).to_string());
                }
                step["pose"] = json!({"rotation":rotation,"translation":["7/3","-5/7","2"]});
                let report = check_with_options(&recipe(vec![step]), &controls);
                assert_eq!(
                    report.verdict,
                    Verdict::Agree,
                    "{op} size={size} rotation={rotation:?}: {report:?}"
                );
            }
        }
    }
}
#[test]
fn booleans_and_merged_origin_sets() {
    for kind in ["union", "subtract", "intersect"] {
        let mut tool = block("b", "2", "1");
        if kind == "subtract" {
            tool["pose"]["translation"] = json!(["1", "1", "1"]);
        }
        let report = check(&recipe(vec![
            block("a", "2", "0"),
            tool,
            json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"b"}),
        ]));
        assert_eq!(report.verdict, Verdict::Agree, "{kind}: {report:?}");
    }
}
#[test]
fn independent_perturbations_prove_gaps_in_both_directions() {
    let a = recipe(vec![block("a", "2", "0")]);
    let b = recipe(vec![block("a", "2", "1/10")]);
    let report = compare(&a, &b, &options());
    assert_eq!(report.verdict, Verdict::Disagree);
    for pass in [Pass::Forward, Pass::Reverse] {
        assert!(
            report
                .findings
                .iter()
                .any(|f| f.pass == pass && f.reason == Reason::Gap && f.gap.unwrap().lo() > 0.01),
            "{report:?}"
        );
    }
}
#[test]
fn reverse_check_catches_boundary_missing_from_the_brep() {
    let a = recipe(vec![block("a", "2", "0")]);
    let b = recipe(vec![
        block("a", "2", "0"),
        block("b", "1", "5"),
        json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let report = compare(&a, &b, &options());
    assert_eq!(report.verdict, Verdict::Disagree);
    assert!(
        !report
            .findings
            .iter()
            .any(|f| f.pass == Pass::Forward && f.reason == Reason::Gap)
    );
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.pass == Pass::Reverse && f.reason == Reason::Gap)
    );
}
#[test]
fn provenance_changes_are_disagreements_even_at_zero_gap() {
    let report = compare(
        &recipe(vec![block("a", "2", "0")]),
        &recipe(vec![block("renamed", "2", "0")]),
        &options(),
    );
    assert_eq!(report.verdict, Verdict::Disagree);
    assert!(report.findings.iter().any(|f| f.reason == Reason::Origins
        && !f.brep_faces.is_empty()
        && !f.field_faces.is_empty()));
}
#[test]
fn unbounded_boolean_upper_distance_never_agrees() {
    let a = recipe(vec![
        block("a", "2", "0"),
        block("b", "1", "5"),
        json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let b = recipe(vec![
        block("a", "2", "1/10000"),
        block("b", "1", "5"),
        json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let report = compare(&a, &b, &options());
    assert_eq!(report.verdict, Verdict::Unresolved, "{report:?}");
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.pass == Pass::Forward && matches!(f.reason, Reason::Unresolved(_)))
    );
}
#[test]
fn region_filters_both_directions_and_empty_scope_is_unresolved() {
    let a = recipe(vec![
        block("a", "2", "0"),
        block("b", "1", "5"),
        json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let b = recipe(vec![
        block("a", "2", "0"),
        block("b", "1", "6"),
        json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
    ]);
    let region = Region {
        min: std::array::from_fn(|_| (-1).into()),
        max: std::array::from_fn(|_| 3.into()),
    };
    let options = Options {
        region: Some(region),
        ..options()
    };
    assert_eq!(compare(&a, &b, &options).verdict, Verdict::Agree);
    let options = Options {
        region: Some(Region {
            min: std::array::from_fn(|_| 100.into()),
            max: std::array::from_fn(|_| 101.into()),
        }),
        ..options
    };
    assert_eq!(compare(&a, &b, &options).verdict, Verdict::Unresolved);
}

#[test]
fn a_remote_origin_change_on_a_merged_face_is_not_a_local_disagreement() {
    let part = |tool: &str| {
        recipe(vec![
            block("a", "2", "0"),
            block(tool, "2", "2"),
            json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":tool}),
        ])
    };
    let options = Options {
        region: Some(Region {
            min: std::array::from_fn(|_| (-1).into()),
            max: [1.into(), 3.into(), 3.into()],
        }),
        ..options()
    };
    let report = compare(&part("b"), &part("renamed"), &options);
    assert_eq!(report.verdict, Verdict::Unresolved, "{report:?}");
}
#[test]
fn unsupported_operations_are_all_reported() {
    let recipe = recipe(vec![
        block("a", "2", "0"),
        json!({"id":"offset","op":"offset","body":"a","distance":"1"}),
        json!({"id":"mirror","op":"mirror","body":"offset","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}),
    ]);
    let report = check(&recipe);
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(
        report
            .not_comparable
            .iter()
            .map(|o| o.operation.as_str())
            .collect::<Vec<_>>(),
        ["offset", "mirror"]
    );
}
#[test]
fn complete_reports_are_identical_across_one_two_and_four_workers() {
    let parts = [
        recipe(vec![
            block("a", "2", "0"),
            block("b", "2", "1"),
            json!({"id":"result","op":"boolean","kind":"subtract","target":"a","tool":"b"}),
        ]),
        recipe(vec![
            sketch_rectangle("profile", "0", "0", "4", "3"),
            json!({"id":"feature","op":"extrude","profile":"profile","distance":"5","direction":["3","0","4"]}),
        ]),
    ];
    for recipe in parts {
        let expected = check(&recipe);
        for workers in [1, 2, 4] {
            std::thread::scope(|scope| {
                let jobs: Vec<_> = (0..workers)
                    .map(|_| scope.spawn(|| check(&recipe)))
                    .collect();
                for job in jobs {
                    assert_eq!(job.join().unwrap(), expected);
                }
            });
        }
    }
}
#[test]
fn cancelled_and_superseded_runs_cannot_publish_or_clear_newer_results() {
    let recipe = recipe(vec![block("a", "2", "0")]);
    let controller = Controller::default();
    let old = controller.begin();
    let current = controller.begin();
    assert_eq!(old.check(&recipe, &options()), Completion::Cancelled);
    assert!(controller.latest().is_none());
    assert_eq!(current.check(&recipe, &options()), Completion::Published);
    let report = controller.latest().unwrap();
    old.cancel();
    assert_eq!(controller.latest(), Some(report));
    current.cancel();
    assert!(controller.latest().is_none());
}
#[test]
fn cancellation_at_publication_discards_a_complete_computation() {
    let controller = Controller::default();
    let old = controller.begin();
    let mut report = check(&recipe(vec![block("a", "2", "0")]));
    report.revision = old.id();
    controller.begin();
    assert_eq!(old.publish(report), Completion::Cancelled);
    assert!(controller.latest().is_none());
}

#[test]
fn removed_coplanar_boolean_zeros_stay_unresolved() {
    let report = check(&recipe(vec![
        block("a", "2", "0"),
        block("b", "2", "1"),
        json!({"id":"result","op":"boolean","kind":"subtract","target":"a","tool":"b"}),
    ]));
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.pass == Pass::Reverse && matches!(f.reason, Reason::Unresolved(_)))
    );
}

#[test]
fn certified_empty_results_and_empty_nonempty_mismatches() {
    let empty = recipe(vec![
        block("a", "2", "0"),
        block("b", "2", "0"),
        json!({"id":"result","op":"boolean","kind":"subtract","target":"a","tool":"b"}),
    ]);
    assert_eq!(check(&empty).verdict, Verdict::Agree);
    let solid = recipe(vec![block("a", "2", "0")]);
    for (a, b, pass) in [
        (&empty, &solid, Pass::Reverse),
        (&solid, &empty, Pass::Forward),
    ] {
        let report = compare(a, b, &options());
        assert_eq!(report.verdict, Verdict::Disagree);
        assert!(report.findings.iter().any(|f| f.pass == pass
            && f.reason == Reason::Gap
            && f.gap.unwrap().lo() == f64::INFINITY));
    }
}

#[test]
fn repeated_body_operand_is_reported_without_entering_the_kernel_assertion() {
    let report = check(&recipe(vec![
        block("a", "2", "0"),
        json!({"id":"result","op":"boolean","kind":"subtract","target":"a","tool":"a"}),
    ]));
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(report.not_comparable[0].step.as_str(), "result");
    assert_eq!(report.not_comparable[0].operation, "boolean");
}

#[test]
fn cancellation_after_forward_pass_never_exposes_partial_or_stale_results() {
    use std::sync::{Arc, Barrier};
    let recipe = recipe(vec![block("a", "2", "0")]);
    let controller = Controller::default();
    let mut old = controller.begin();
    let barriers = Arc::new((Barrier::new(2), Barrier::new(2)));
    old.checkpoint = Some(barriers.clone());
    std::thread::scope(|scope| {
        let job = scope.spawn(|| old.check(&recipe, &options()));
        barriers.0.wait();
        assert!(controller.latest().is_none());
        let current = controller.begin();
        assert_eq!(current.check(&recipe, &options()), Completion::Published);
        let complete = controller.latest().unwrap();
        barriers.1.wait();
        assert_eq!(job.join().unwrap(), Completion::Cancelled);
        assert_eq!(controller.latest(), Some(complete));
    });
}

#[test]
fn reversed_exact_face_normal_is_a_disagreement() {
    use catallax_topology::{Body, Lump, Sense, Shell};
    let recipe = recipe(vec![
        json!({"id":"ball","op":"sphere","radius":"1","pose":pose("0")}),
    ]);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let old = built.model.faces_of_body(built.result)[0];
    let mut face = built.model.face(old).clone();
    face.sense = Sense::Reversed;
    let origins = built.model.face_origins(old).to_vec();
    let face = built.model.add_named_face(face, origins);
    let shell = built.model.add_shell(Shell { faces: vec![face] });
    let lump = built.model.add_lump(Lump {
        shells: vec![shell],
    });
    built.result = built.model.add_body(Body {
        kind: catallax_topology::BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    let controller = Controller::default();
    assert_eq!(
        controller.begin().check_built(&built, &recipe, &options()),
        Completion::Published
    );
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Disagree);
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.reason == Reason::Direction)
    );
}

#[test]
fn reverse_projection_catches_curved_perturbations_above_display_accuracy() {
    for op in ["sphere", "cylinder", "cone", "torus"] {
        let step = |radius: &str| match op {
            "sphere" => json!({"id":"stock","op":op,"radius":radius,"pose":pose("0")}),
            "cylinder" => {
                json!({"id":"stock","op":op,"radius":radius,"height":"1/5","pose":pose("0")})
            }
            "cone" => {
                json!({"id":"stock","op":op,"base_radius":radius,"top_radius":"0","height":"1/5","pose":pose("0")})
            }
            _ => json!({"id":"stock","op":op,"major":"3/10","minor":radius,"pose":pose("0")}),
        };
        let report = compare(
            &recipe(vec![step("1/10")]),
            &recipe(vec![step("1/5")]),
            &options(),
        );
        assert_eq!(report.verdict, Verdict::Disagree, "{op}: {report:?}");
        // Enlarged cylinder caps still contain the old cap rims. A vertex-only
        // forward pass need not detect the changed cylindrical wall.
        assert!(
            report.findings.iter().any(|f| f.pass == Pass::Reverse
                && f.reason == Reason::Gap
                && f.gap.unwrap().lo() > 0.01),
            "{op}: {report:?}"
        );
    }
}

#[test]
fn extra_brep_origin_cannot_be_mistaken_for_insufficient_sampling() {
    use catallax_topology::{Body, Lump, Shell};
    let recipe = recipe(vec![
        json!({"id":"ball","op":"sphere","radius":"1/10","pose":pose("0")}),
    ]);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let old = built.model.faces_of_body(built.result)[0];
    let mut origins = built.model.face_origins(old).to_vec();
    let mut unknown = origins[0].clone();
    unknown.step = StepId::new("ghost").unwrap();
    origins.push(unknown.clone());
    let face = built
        .model
        .add_named_face(built.model.face(old).clone(), origins);
    let shell = built.model.add_shell(Shell { faces: vec![face] });
    let lump = built.model.add_lump(Lump {
        shells: vec![shell],
    });
    built.result = built.model.add_body(Body {
        kind: catallax_topology::BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    let controller = Controller::default();
    controller.begin().check_built(&built, &recipe, &options());
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Disagree);
    assert!(report.findings.iter().any(|f| f.reason == Reason::Origins
        && f.brep_faces.iter().any(|names| names.contains(&unknown))));
    assert!(!report.findings.iter().any(|f| f.reason == Reason::Gap));
}

#[test]
fn round_boolean_checks_keep_unbounded_vertex_distances_unresolved() {
    for op in ["sphere", "cylinder"] {
        for kind in ["union", "subtract", "intersect"] {
            let solid = |id: &str, x: &str| {
                if op == "sphere" {
                    json!({"id":id,"op":op,"radius":"1/5","pose":pose(x)})
                } else {
                    json!({"id":id,"op":op,"radius":"1/5","height":"2/5","pose":pose(x)})
                }
            };
            let report = check(&recipe(vec![
                solid("a", "0"),
                solid("b", "1/5"),
                json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"b"}),
            ]));
            assert_eq!(
                report.verdict,
                Verdict::Unresolved,
                "{op} {kind}: {report:?}"
            );
            assert!(
                report.findings.iter().any(|f| f.pass == Pass::Forward
                    && matches!(f.reason, Reason::Unresolved(_))
                    && f.brep_faces
                        .iter()
                        .flatten()
                        .any(|o| o.step.as_str() == "a")),
                "{op} {kind}: {report:?}"
            );
        }
    }
}

fn sketch_rectangle(id: &str, x: &str, y: &str, width: &str, height: &str) -> Value {
    let x = Rational::parse(x).unwrap();
    let y = Rational::parse(y).unwrap();
    let right = &x + &Rational::parse(width).unwrap();
    let top = &y + &Rational::parse(height).unwrap();
    json!({"id":id,"op":"sketch","plane":pose("0"),"loops":[[
        {"kind":"line","id":"bottom","to":[right.to_string(),y.to_string()]},
        {"kind":"line","id":"right","to":[right.to_string(),top.to_string()]},
        {"kind":"line","id":"top","to":[x.to_string(),top.to_string()]},
        {"kind":"line","id":"left","to":[x.to_string(),y.to_string()]}
    ]]})
}

#[test]
fn feature_extrusions_compare_signed_and_oblique_sweeps() {
    for distance in ["2", "-2"] {
        let part = recipe(vec![
            sketch_rectangle("profile", "0", "0", "4", "3"),
            json!({"id":"feature","op":"extrude","profile":"profile","distance":distance}),
        ]);
        assert_eq!(check(&part).verdict, Verdict::Agree);
    }
    let part = recipe(vec![
        sketch_rectangle("profile", "0", "0", "4", "3"),
        json!({"id":"feature","op":"extrude","profile":"profile","distance":"5","direction":["3","0","4"]}),
    ]);
    let report = check(&part);
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn feature_through_slot_uses_the_field_target_bound() {
    let part = recipe(vec![
        block("stock", "4", "0"),
        sketch_rectangle("profile", "1", "1", "2", "2"),
        json!({"id":"cut","op":"extrude_through","target":"stock","profile":"profile"}),
    ]);
    let controller = Controller::default();
    let options = Options {
        divisions: 3,
        ..options()
    };
    controller.begin().check(&part, &options);
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn feature_reverse_check_catches_an_extended_extrusion() {
    let part = |height: &str| {
        recipe(vec![
            sketch_rectangle("profile", "0", "0", "4", "3"),
            json!({"id":"feature","op":"extrude","profile":"profile","distance":height}),
        ])
    };
    let report = compare(&part("2"), &part("3"), &options());
    assert_eq!(report.verdict, Verdict::Disagree, "{report:?}");
    assert!(
        report.findings.iter().any(|f| f.pass == Pass::Reverse
            && f.reason == Reason::Gap
            && f.gap.unwrap().lo() > 0.01)
    );
}

#[test]
fn feature_curved_extrusions_and_full_partial_revolutions_keep_unbounded_uppers() {
    let extrusion = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"1/5"}]]}),
        json!({"id":"feature","op":"extrude","profile":"profile","distance":"2/5"}),
    ]);
    let mut parts = vec![extrusion];
    for angle in ["90", "180", "270", "360"] {
        parts.push(recipe(vec![
            sketch_rectangle("profile","1/5","0","1/5","2/5"),
            json!({"id":"feature","op":"revolve","profile":"profile","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":angle}),
        ]));
    }
    for part in parts {
        let report = check(&part);
        assert_eq!(report.verdict, Verdict::Unresolved, "{report:?}");
        assert!(
            report.findings.iter().any(|f| f.pass == Pass::Forward
                && matches!(f.reason, Reason::Unresolved(_))
                && f.brep_faces
                    .iter()
                    .flatten()
                    .any(|o| o.step.as_str() == "feature")),
            "{report:?}"
        );
        assert!(
            !report
                .findings
                .iter()
                .any(|f| matches!(f.reason, Reason::Gap | Reason::Origins | Reason::Direction)),
            "{report:?}"
        );
    }
}

#[test]
fn feature_empty_profile_is_proven_empty_without_samples() {
    let part = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[]}),
        json!({"id":"feature","op":"extrude","profile":"profile","distance":"2"}),
    ]);
    assert_eq!(check(&part).verdict, Verdict::Agree);
}

#[test]
fn feature_planar_agreement_across_dyadic_poses_and_scales() {
    for scale in ["1/2", "2", "8"] {
        let scale = Rational::parse(scale).unwrap();
        for rotation in [["1", "0", "0", "0"], ["1", "0", "1", "0"]] {
            let mut profile = sketch_rectangle(
                "profile",
                "0",
                "0",
                &(&scale * &4.into()).to_string(),
                &(&scale * &3.into()).to_string(),
            );
            profile["plane"] = json!({"rotation":rotation,"translation":["7/2","-3/4","1/2"]});
            let part = recipe(vec![
                profile,
                json!({"id":"feature","op":"extrude","profile":"profile","distance":(&scale*&2.into()).to_string()}),
            ]);
            let report = check(&part);
            assert_eq!(
                report.verdict,
                Verdict::Agree,
                "scale={scale} rotation={rotation:?}: {report:?}"
            );
        }
    }
}

#[test]
fn an_exhausted_mesh_allowance_stays_unresolved_without_partial_queries() {
    let part = recipe(vec![
        json!({"id":"stock","op":"sphere","radius":"1/10","pose":pose("0")}),
    ]);
    let mut controls = options();
    controls.mesh.max_vertices = 4;
    let report = check_with_options(&part, &controls);
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(report.forward_vertices, 0);
    assert_eq!(report.reverse_samples, 0);
    assert!(
        report
            .findings
            .iter()
            .all(|f| matches!(f.reason, Reason::Unresolved(_)))
    );
}
