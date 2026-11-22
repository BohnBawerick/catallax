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

fn m11_loft([a, b]: [usize; 2], planar: bool, heights: [&str; 2]) -> Recipe {
    let polygon = |prefix: &str, points: [[i64; 2]; 4]| {
        points
            .into_iter()
            .enumerate()
            .map(|(i, [x, y])| {
                json!({"kind":"line","id":format!("{prefix}_{i}"),
                "to":[x.to_string(),y.to_string()]})
            })
            .collect::<Vec<_>>()
    };
    let side = |i: usize| json!({"step":"loft","role":{"side":{"segment":format!("a_{i}")}}});
    let section = |i: usize, height: &str, loop_: &[_]| {
        json!({"id":format!("section{i}"),"op":"sketch","loops":[loop_],
          "plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]}})
    };
    let bottom = polygon("a", [[1, -1], [1, 1], [-1, 1], [-1, -1]]);
    let top = if planar {
        [[1, 0], [1, 2], [-2, 0], [0, -1]]
    } else {
        [[1, 0], [0, 2], [-2, 0], [0, -1]]
    };
    let top = polygon("b", top);
    Recipe::parse(
        &json!({
            "format":"catallax-recipe/5","name":"M11 loft comparison","steps":[
                section(0,heights[0],&bottom), section(1,heights[1],&top),
                {"id":"loft","op":"loft","shape":"ruled","sections":[
                    {"profile":"section0","start":["a_0"]},{"profile":"section1","start":["b_0"]}]},
                {"id":"round","op":"fillet","body":"loft","radius":"1/10",
                    "edges":[{"between":[side(a),side(b)]}]}
            ],"result":"round"
        })
        .to_string(),
    )
    .unwrap()
}

#[test]
fn loft_spline_rounds_report_display_coverage_or_a_named_brep_refusal() {
    for ([a, b], planar) in [([0, 1], false), ([1, 2], false), ([0, 1], true)] {
        let recipe = m11_loft([a, b], planar, ["0", "2"]);
        let (report, work) =
            catallax_numbers::with_operation_budget(100_000_000, || check(&recipe));
        println!(
            "loft edges {a}/{b} planar {planar}: {} forward vertices, {} work, upper {:?}, refusal {:?}",
            report.forward_vertices,
            work.used,
            report
                .sampled_distance_upper
                .as_ref()
                .map(Rational::to_f64_nearest),
            report.refusal
        );
        assert!(!work.exhausted, "{report:?}");
        if planar {
            assert_eq!(report.forward_vertices, 0);
            assert_eq!(
                report.refusal.as_ref().map(|e| e.code),
                Some("fillet.unsupported-contact")
            );
            assert!(report.findings.iter().any(|f| {
                matches!(&f.reason, Reason::Unresolved(reason) if reason.starts_with("B-rep build"))
            }));
            continue;
        }
        assert_eq!(report.forward_vertices, if a == 0 { 245 } else { 335 });
        assert!(report.refusal.is_none(), "{report:?}");
        let findings: Vec<_> = report
            .findings
            .iter()
            .filter(|f| f.pass == Pass::Forward)
            .collect();
        assert!(findings.is_empty(), "{findings:?}");
    }
}

#[test]
fn analytic_named_boss_loop_keeps_the_base_refusal() {
    let sides = ["x_min", "y_min", "x_max", "y_max"];
    let side = |segment| json!({"step":"block","role":{"side":{"segment":segment}}});
    let cap = json!({"step":"block","role":"end"});
    let upright_edges = (0..4)
        .map(|i| json!({"between":[side(sides[i]), side(sides[(i + 1) % 4])]}))
        .collect::<Vec<_>>();
    let mut top_edges = sides
        .map(|segment| json!({"between":[cap, side(segment)]}))
        .to_vec();
    top_edges.push(json!({"between":[cap, {"step":"upright"}]}));
    let mut stock = block("block", "20", "0");
    stock["size"][2] = json!("10");
    let boss = Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"fillet over fillet","result":"top","steps":[
            stock,
            {"id":"upright","op":"fillet","body":"block","radius":"1","edges":upright_edges},
            {"id":"top","op":"fillet","body":"upright","radius":"3","edges":top_edges}
        ]})
        .to_string(),
    )
    .unwrap();
    let report = check(&boss);
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(
        report.refusal.as_ref().map(|error| error.code),
        Some("field.rolling-path-reach"),
    );
    assert_eq!(report.sampled_distance_upper, Some(Rational::zero()));
}

#[test]
fn analytic_forward_gap_leaves_a_refused_comparison_unresolved() {
    use catallax_numbers::with_operation_budget;

    let brep = recipe(vec![block("stock", "4", "0")]);
    let field = recipe(vec![block("stock", "4", "5")]);
    let built = evaluate(&brep, &mut Decider::default()).unwrap();
    let compare = || {
        let controller = Controller::default();
        assert_eq!(
            controller.begin().check_built(&built, &field, &options()),
            Completion::Published,
        );
        controller.latest().unwrap()
    };
    let (complete, work) = with_operation_budget(1_000_000, compare);
    assert_eq!(complete.verdict, Verdict::Disagree);
    assert!(complete.refusal.is_none());
    let (partial, spent) = with_operation_budget(work.used - 1, compare);
    assert!(spent.exhausted);
    assert!(
        partial
            .findings
            .iter()
            .any(|finding| { finding.pass == Pass::Forward && finding.reason == Reason::Gap })
    );
    assert_eq!(partial.verdict, Verdict::Unresolved);
    assert_eq!(
        partial.refusal.as_ref().map(|error| error.code),
        Some("numbers.operation-budget"),
    );
    assert_eq!(partial.sampled_distance_upper, None);
}

#[test]
fn analytic_build_and_rounded_copy_refusals_remain_typed() {
    let analytic = recipe(vec![
        block("stock", "4", "0"),
        json!({"id":"round","op":"fillet","body":"stock","radius":"1/10",
            "edges":[{"between":[{"step":"stock","role":{"side":{"segment":"x_min"}}},
                {"step":"stock","role":{"side":{"segment":"x_max"}}}]}]}),
    ]);
    let report = check(&analytic);
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(
        report.refusal.as_ref().map(|error| error.code),
        Some("fillet.lost-edge"),
    );
    assert_eq!(report.sampled_distance_upper, Some(Rational::zero()));
    assert!(report.findings.iter().any(|finding| {
        matches!(&finding.reason, Reason::Unresolved(reason) if reason.starts_with("field build "))
    }));

    let original = m11_loft([0, 1], false, ["0", "2"]);
    for copy in [
        json!({"id":"copy","op":"mirror","body":"round","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}),
        json!({"id":"copy","op":"mirror","body":"round","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":true}),
        json!({"id":"copy","op":"pattern","body":"round","layout":{"linear":{"vector":["4","0","0"],"count":2}}}),
    ] {
        let mut value = serde_json::to_value(&original).unwrap();
        value["steps"].as_array_mut().unwrap().push(copy);
        value["result"] = json!("copy");
        let report = check(&Recipe::parse(&value.to_string()).unwrap());
        assert_eq!(
            report.refusal.as_ref().map(|e| e.code),
            Some("field.spline-copy"),
            "{report:?}"
        );
        assert_eq!(report.sampled_distance_upper, None);
    }
}

#[test]
fn spline_build_refusals_have_no_sampled_distance_bound() {
    use catallax_numbers::with_operation_budget;

    let missing_edge = m11_loft([0, 2], false, ["0", "2"]);
    let valid_edge = m11_loft([0, 1], false, ["0", "2"]);
    let budget = 100_000;
    let (unrounded, _) =
        with_operation_budget(budget, || catallax_field::Plan::compile(&valid_edge));
    assert_eq!(unrounded.unwrap_err().code, "fillet.lost-edge");

    let (prerequisite_refusal, no_work) = with_operation_budget(0, || check(&valid_edge));
    assert!(no_work.exhausted);
    assert_eq!(
        prerequisite_refusal.refusal.as_ref().unwrap().code,
        "numbers.operation-budget"
    );
    assert_eq!(
        prerequisite_refusal.sampled_distance_upper,
        Some(Rational::zero())
    );
    let mut translated = serde_json::to_value(&valid_edge).unwrap();
    translated["steps"][1]["loops"] = translated["steps"][0]["loops"].clone();
    translated["steps"][2]["sections"][1]["start"] = json!(["a_0"]);
    translated["steps"][1]["plane"]["rotation"] = json!(["2", "0", "0", "0"]);
    let translated = Recipe::parse(&translated.to_string()).unwrap();
    let (analytic, _) = with_operation_budget(0, || check(&translated));
    assert_eq!(
        analytic.refusal.as_ref().map(|error| error.code),
        Some("numbers.operation-budget")
    );
    assert_eq!(analytic.sampled_distance_upper, Some(Rational::zero()));

    let mut union = serde_json::to_value(&valid_edge).unwrap();
    let steps = union["steps"].as_array_mut().unwrap();
    steps[3]["body"] = json!("same");
    steps.insert(
        3,
        json!({"id":"same","op":"boolean","kind":"union","target":"loft","tool":"loft"}),
    );
    let union = Recipe::parse(&union.to_string()).unwrap();
    let (union_refusal, work) = with_operation_budget(budget, || check(&union));
    assert!(work.exhausted);

    let (budget_refusal, work) = with_operation_budget(budget, || check(&valid_edge));
    assert!(work.exhausted);
    for (report, code) in [
        (check(&missing_edge), "fillet.lost-edge"),
        (union_refusal, "numbers.operation-budget"),
        (budget_refusal, "numbers.operation-budget"),
    ] {
        assert_eq!(report.refusal.as_ref().map(|error| error.code), Some(code));
        assert_eq!(report.verdict, Verdict::Unresolved);
        assert_eq!(report.forward_vertices, 0);
        assert_eq!(report.sampled_distance_upper, None);
    }
}

#[test]
fn invalid_spline_prerequisites_keep_the_sketch_refusal() {
    let mut value = serde_json::to_value(m11_loft([0, 1], false, ["0", "2"])).unwrap();
    for (section, prefix, scale) in [(0, 'a', 1), (1, 'b', 2)] {
        let edges = [(0, 4), (3, -3), (-4, 1), (4, 1), (-3, -3)]
            .into_iter()
            .enumerate()
            .map(|(i, (x, y))| {
                json!({"kind":"line","id":format!("{prefix}_{i}"),
                    "to":[(scale * x).to_string(),y.to_string()]})
            })
            .collect::<Vec<_>>();
        value["steps"][section]["loops"] = json!([edges]);
    }

    let report = check(&Recipe::parse(&value.to_string()).unwrap());
    assert_eq!(
        report.refusal.as_ref().unwrap().code,
        "sketch.self-intersection"
    );
    assert_eq!(report.sampled_distance_upper, Some(Rational::zero()));
}

#[test]
fn spline_prerequisite_refusal_does_not_expand_shared_unions() {
    let mut value = serde_json::to_value(m11_loft([0, 1], false, ["0", "2"])).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    let mut fillet = steps.pop().unwrap();
    let mut inputs = ["loft".to_owned(), "loft".to_owned()];
    for depth in 0..30 {
        let outputs = [format!("left{depth}"), format!("right{depth}")];
        for id in &outputs {
            steps.push(json!({"id":id,"op":"boolean","kind":"union",
                "target":inputs[0],"tool":inputs[1]}));
        }
        inputs = outputs;
    }
    fillet["body"] = json!(inputs[0]);
    steps.push(fillet);

    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let (report, work) = catallax_numbers::with_operation_budget(0, || check(&recipe));
    assert!(work.exhausted);
    assert_eq!(
        report.refusal.as_ref().unwrap().code,
        "numbers.operation-budget"
    );
    assert_eq!(report.sampled_distance_upper, Some(Rational::zero()));
}

#[test]
fn a_mirrored_loft_can_be_queried_but_not_filleted_with_stale_origins() {
    use catallax_field::{Plan, Workspace};

    let mut value = serde_json::to_value(m11_loft([0, 1], false, ["0", "2"])).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    let mut fillet = steps.pop().unwrap();
    steps.push(json!({"id":"copy","op":"mirror","body":"loft",
        "plane":{"point":["5","0","0"],"normal":["1","0","0"]},"keep_original":false}));
    value["result"] = json!("copy");
    let mirrored = Recipe::parse(&value.to_string()).unwrap();
    let plan = Plan::compile(&mirrored).unwrap();
    let mut workspace = Workspace::default();
    assert!(plan.fast([10.0, 0.0, 1.0], &mut workspace).unwrap().value < 0.0);
    assert!(plan.fast([0.0, 0.0, 1.0], &mut workspace).unwrap().value > 0.0);

    fillet["body"] = json!("copy");
    value["steps"].as_array_mut().unwrap().push(fillet);
    value["result"] = json!("round");
    let report = check(&Recipe::parse(&value.to_string()).unwrap());
    assert_eq!(
        report.refusal.as_ref().map(|e| e.code),
        Some("fillet.unsupported-spline")
    );
    assert_eq!(report.sampled_distance_upper, None);
}

#[test]
fn scoped_spline_display_selects_either_blend_in_a_disjoint_body() {
    use catallax_shapes::{Isometry, RVector3, SplineSupport, Surface};
    let recipe = m11_loft([0, 1], false, ["0", "2"]);
    let ((), work) = catallax_numbers::with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let mut built = evaluate(&recipe, &mut dec).unwrap();
        let map = Isometry::translation(RVector3::new(
            Rational::from_i64(8),
            Rational::zero(),
            Rational::zero(),
        ));
        let copied =
            catallax_topology::transformed_body(&mut built.model, built.result, &map, Clone::clone);
        let mut body = built.model.body(built.result).clone();
        body.lumps
            .extend_from_slice(&built.model.body(copied).lumps);
        let body = built.model.add_body(body);
        let mut extents = Vec::new();
        for face in built.model.faces_of_body(body) {
            if !matches!(
                built.model.geom().surface(built.model.face(face).surface),
                Surface::BSpline(SplineSupport::SplineFillet(_))
            ) {
                continue;
            }
            let display = catallax_mesh::display_spline_fillet_face(
                &built.model,
                body,
                face,
                &Options::default().mesh,
                &mut dec,
            )
            .unwrap();
            assert_eq!(display.measurement.covered_faces, 1);
            assert_eq!(display.geometry.positions.len(), 245);
            let xs: Vec<_> = display.geometry.positions.iter().map(|p| p[0]).collect();
            extents.push((
                xs.iter().copied().fold(f64::INFINITY, f64::min),
                xs.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            ));
        }
        extents.sort_by(|a, b| a.0.total_cmp(&b.0));
        let [original, translated] = extents.as_slice() else {
            panic!("{extents:?}")
        };
        assert!((translated.0 - original.0 - 8.0).abs() < 1e-9);
        assert!((translated.1 - original.1 - 8.0).abs() < 1e-9);
    });
    assert!(!work.exhausted);
}

#[test]
fn loft_vertices_must_fit_the_brep_bound_even_below_display_accuracy() {
    let field = m11_loft([0, 1], false, ["0", "2"]);
    let brep = m11_loft([0, 1], false, ["3/250", "503/250"]);
    let mut options = options();
    options.region = Some(Region {
        min: ["-2", "-2", "503/250"].map(|x| Rational::parse(x).unwrap()),
        max: ["2", "2", "503/250"].map(|x| Rational::parse(x).unwrap()),
    });
    let (report, work) =
        catallax_numbers::with_operation_budget(100_000_000, || compare(&brep, &field, &options));
    assert!(!work.exhausted && report.refusal.is_none(), "{report:?}");
    assert_eq!(report.verdict, Verdict::Disagree, "{report:?}");
    let shifted_cap = Rational::parse("503/250").unwrap().to_f64_nearest();
    assert!(
        report.findings.iter().any(|f| {
            f.pass == Pass::Forward
                && f.reason == Reason::Gap
                && (f.point.z.to_f64_nearest() - shifted_cap).abs() < 1e-9
                && f.gap.unwrap().lo() < options.accuracy.to_f64_nearest()
        }),
        "{report:?}"
    );
    let (partial, spent) =
        catallax_numbers::with_operation_budget(work.used - 1, || compare(&brep, &field, &options));
    assert!(spent.exhausted && partial.refusal.is_some(), "{partial:?}");
    assert_eq!(partial.verdict, Verdict::Disagree, "{partial:?}");
    assert_eq!(partial.forward_vertices + 1, report.forward_vertices);
    assert!(partial.sampled_distance_upper.is_none());
    assert!(partial.findings.iter().any(|f| f.reason == Reason::Gap));
}

#[test]
fn an_empty_spline_region_explains_missing_coverage() {
    let mut options = options();
    options.region = Some(Region {
        min: [10; 3].map(Rational::from_i64),
        max: [11; 3].map(Rational::from_i64),
    });
    let (report, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        check_with_options(&m11_loft([0, 1], false, ["0", "2"]), &options)
    });
    assert!(!work.exhausted && report.refusal.is_none(), "{report:?}");
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert_eq!(report.forward_vertices, 0);
    assert!(report.sampled_distance_upper.is_none());
    assert!(report.findings.iter().any(|f| {
        matches!(&f.reason, Reason::Unresolved(reason) if reason == "no spline display vertices in the requested scope")
    }));
}

fn biquadratic_supports(dec: &mut Decider) -> [catallax_shapes::RationalPatch; 2] {
    use catallax_numbers::Real;
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch};
    let q = Rational::from_i64;
    let quarter = |n| q(n).mul_pow2(-2);
    [false, true].map(|second| {
        let mut controls = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                let (x, y) = if second {
                    (quarter(4 - i * j), quarter(-4 + 4 * i + 2 * j))
                } else {
                    (quarter(-4 + 4 * i + 2 * j - i * j), quarter(-4 + i * j))
                };
                let z = if i == 1 && j == 1 {
                    q(j) + Rational::from_ratio(if second { -1 } else { 1 }, 20).unwrap()
                } else {
                    q(j)
                };
                controls.push(Point3::new(x.into(), y.into(), z.into()));
            }
        }
        let patch = BSplineSurface::new(
            [2, 2],
            std::array::from_fn(|_| [0, 0, 0, 1, 1, 1].map(q).to_vec()),
            [3, 3],
            controls,
            vec![Real::one(); 9],
            dec,
        )
        .unwrap();
        RationalPatch::new(patch, dec).unwrap()
    })
}

fn biquadratic_spine(
    dec: &mut Decider,
) -> Result<catallax_field::spline_spine::Spine, catallax_numbers::Refusal> {
    use catallax_field::spline_spine::{Spine, SpineDefinition};
    Spine::prove(
        &SpineDefinition {
            supports: biquadratic_supports(dec),
            radius: Rational::parse("1/10").unwrap(),
            sides: [catallax_numbers::Sign::Negative; 2],
            axis: 1,
            seed: Rational::parse("1/2").unwrap(),
            tolerance: Rational::parse("1/1000").unwrap(),
        },
        dec,
    )
}

#[test]
fn spline_comparison_field_builds_agree_at_one_two_and_four_workers() {
    use catallax_field::{Plan, spline_spine::sweep::Sweep};
    use catallax_numbers::{set_threads, threads, with_operation_budget};
    use std::num::NonZeroUsize;

    let mut baseline = None;
    for count in [1, 2, 4] {
        set_threads(NonZeroUsize::new(count));
        assert_eq!(threads(), count);
        let mut records = Vec::new();
        for (name, edges, planar) in [
            ("M27", [0, 1], false),
            ("next loft", [1, 2], false),
            ("spline-plane", [0, 1], true),
        ] {
            let recipe = m11_loft(edges, planar, ["0", "2"]);
            let (plan, work) = with_operation_budget(20_000_000, || {
                Plan::compile_for_comparison(
                    &recipe,
                    &catallax_features::ImportFiles::new(),
                    &mut Decider::default(),
                )
            });
            let plan = plan.unwrap();
            assert!(plan.spline_blend_origin().is_some());
            assert!(!work.exhausted, "{name} at {count} workers");
            records.push((name, work.used, format!("{plan:?}")));
        }
        let (bare, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let spine = biquadratic_spine(&mut dec)?;
            let sweep = Sweep::new(&spine, &mut dec)?;
            Ok::<_, catallax_numbers::Refusal>(format!("{sweep:?}"))
        });
        assert!(!work.exhausted, "biquadratic at {count} workers");
        records.push(("biquadratic", work.used, bare.unwrap()));
        set_threads(None);
        for (name, used, _) in &records {
            println!("{name} field build at {count} workers: {used} work");
        }
        if let Some(expected) = &baseline {
            assert!(
                &records == expected,
                "field build or work changed at {count} workers"
            );
        } else {
            baseline = Some(records);
        }
    }
}

#[test]
fn biquadratic_carrier_vertices_fit_the_field_within_one_request() {
    use catallax_field::spline_spine::sweep::Sweep;
    use catallax_numbers::{Refusal, Sign, current_operation_work, with_operation_budget};
    use catallax_shapes::{SplineFillet, SplineFilletDefinition};
    let options = options();
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let radius = Rational::parse("1/10").unwrap();
        let sides = [Sign::Negative; 2];
        let spine = biquadratic_spine(&mut dec)?;
        let sweep = Sweep::new(&spine, &mut dec)?;
        let field_work = current_operation_work().used;
        assert!(field_work <= 20_000_000);

        let definition = SplineFilletDefinition::new(
            biquadratic_supports(&mut dec),
            radius,
            sides,
            1,
            [Rational::parse("1/32").unwrap(), Rational::one()],
        )?;
        let fit = SplineFillet::fit(definition, &mut dec)?;
        let fitted_work = current_operation_work().used;
        let fit_work = fitted_work - field_work;
        assert!(fit_work <= 20_000_000);
        let grid = catallax_mesh::display_spline_fillet_carrier(&fit, &options.mesh, &mut dec)?;
        let display_work = current_operation_work().used - fitted_work;
        assert_eq!(grid.dimensions, [4, 26]);
        assert!(grid.maximum_deviation <= options.mesh.chord_tolerance);
        let error = &fit.bound().upper + &grid.maximum_deviation;
        let mut report = Report {
            import_agreement: None,
            revision: 0,
            region: None,
            verdict: Verdict::Unresolved,
            forward_vertices: 0,
            reverse_samples: 0,
            sampled_distance_upper: Some(Rational::zero()),
            findings: vec![],
            not_comparable: vec![],
            refusal: None,
        };
        for point in grid.positions {
            let [x, y, z] = point.map(|x| Rational::from_f64_exact(x).unwrap());
            let p = RPoint3::new(x, y, z);
            let range = sweep.distance_at_most(&p, &error, &mut dec)?;
            check::spline_vertex(&mut report, &p, None, range, &error, &options.accuracy);
        }
        assert_eq!(report.forward_vertices, 104);
        assert!(report.findings.is_empty(), "{report:?}");
        println!(
            "biquadratic carrier: {} vertices, {:.9e} mm display bound, {:.9e} mm widened upper, {field_work} field work, {fit_work} fit work, {display_work} display work",
            report.forward_vertices,
            grid.maximum_deviation.to_f64_nearest(),
            report.sampled_distance_upper.unwrap().to_f64_nearest(),
        );
        Ok::<_, Refusal>(())
    });
    println!("biquadratic comparison: {} work", work.used);
    assert!(!work.exhausted, "{result:?}");
    result.unwrap();
}

#[test]
fn thin_boolean_faces_receive_trim_checked_reverse_samples() {
    let part = recipe(vec![
        block("stock", "4", "0"),
        block("tool", "4", "3999/1000"),
        json!({"id":"result","op":"boolean","kind":"intersect","target":"stock","tool":"tool"}),
    ]);
    let report = check_with_options(&part, &Options::default());
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
    assert!(
        report
            .sampled_distance_upper
            .as_ref()
            .is_some_and(|upper| { upper <= &Options::default().accuracy }),
        "{report:?}"
    );
}

#[test]
fn empty_comparison_requires_a_whole_field_exclusion_not_missing_samples() {
    let part = |x: &str| {
        recipe(vec![
            block("stock", "8", "0"),
            block("cut", "4", "1"),
            json!({"id":"small","op":"cylinder","radius":"1/4","height":"1",
            "pose":{"rotation":["1","0","0","0"],"translation":[x,"2","2"]}}),
            json!({"id":"cut_stock","op":"boolean","kind":"subtract","target":"stock","tool":"cut"}),
            json!({"id":"result","op":"boolean","kind":"intersect","target":"cut_stock","tool":"small"}),
        ])
    };
    let empty = part("2");
    assert!(
        !catallax_field::Plan::compile(&empty)
            .unwrap()
            .is_proven_empty()
    );
    let report = check_with_options(&empty, &Options::default());
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
    assert_eq!((report.forward_vertices, report.reverse_samples), (0, 0));
    let protruding = catallax_field::Plan::compile(&part("4750001/1000000")).unwrap();
    assert!(!protruding.prove_empty(&mut Decider::default()).unwrap());
}

#[test]
fn drafted_stock_cut_samples_the_current_definition_not_its_old_cutter_grid() {
    let mut base = block("base", "8", "0");
    base["size"] = json!(["8", "2", "8"]);
    base["pose"]["translation"] = json!(["0", "-2", "0"]);
    let mut notch = block("notch", "4", "2");
    notch["size"] = json!(["4", "4", "8"]);
    notch["pose"]["translation"] = json!(["2", "-3", "1"]);
    let part = recipe(vec![
        base,
        notch,
        json!({"id":"source","op":"boolean","kind":"subtract","target":"base","tool":"notch"}),
        json!({
            "id":"drafted","op":"draft","body":"source",
            "faces":[{"face":{"step":"notch","role":"start"}}],"pull":["1","0","0"],
            "neutral":{"point":["0","0","0"],"normal":["1","0","0"]},"angle":"45"
        }),
    ]);
    let report = check_with_options(&part, &Options::default());
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn drafted_split_origins_cover_each_remaining_geometric_support() {
    let mut seed = block("seed", "8", "0");
    seed["size"] = json!(["14", "8", "4"]);
    let mut left = block("left", "8", "0");
    left["size"] = json!(["6", "8", "4"]);
    let mut tool = block("tool", "4", "6");
    tool["size"] = json!(["4", "10", "6"]);
    tool["pose"]["translation"] = json!(["6", "-1", "-1"]);
    let part = recipe(vec![
        seed,
        left,
        tool,
        json!({"id":"split","op":"boolean","kind":"subtract","target":"seed","tool":"tool"}),
        json!({"id":"source","op":"boolean","kind":"union","target":"left","tool":"split"}),
        json!({
            "id":"drafted","op":"draft","body":"source",
            "faces":[{"face":{"step":"left","role":"end"}}],"pull":["1","0","0"],
            "neutral":{"point":["0","0","0"],"normal":["1","0","0"]},"angle":"45",
            "keep_regions":[["1","4","1"],["11","4","1"]]
        }),
    ]);
    let report = check_with_options(&part, &Options::default());
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn rational_bisector_arc_with_irrational_radius_agrees() {
    let part = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[
            {"kind":"arc","id":"round","centre":["1","1/2"],"to":["0","0"],"turn":"ccw"},
            {"kind":"line","id":"chord","to":["2","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"3"}),
    ]);
    let report = check(&part);
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn reflected_irrational_arc_agrees_at_rounded_boundary_points() {
    let part = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[
            {"kind":"arc","id":"round","centre":["1","1/2"],"to":["0","0"],"turn":"ccw"},
            {"kind":"line","id":"chord","to":["2","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"3"}),
        json!({"id":"reflected","op":"mirror","body":"body","plane":{
            "point":["0","0","0"],"normal":["1","2","0"]
        },"keep_original":false}),
    ]);
    let report = check(&part);
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
}

#[test]
fn arc_just_past_a_chart_seam_keeps_its_shared_display_boundary() {
    let n = Rational::parse("10000000000000000").unwrap();
    let n2 = &n * &n;
    let denominator = &n2 + &Rational::one();
    let x = (&Rational::one() - &n2) / denominator.clone();
    let y = (-Rational::from_i64(2) * n) / denominator;
    let part = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[
            {"kind":"arc","id":"round","centre":["0","0"],"to":[x.to_string(),y.to_string()],"turn":"ccw"},
            {"kind":"line","id":"chord","to":["1","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"2"}),
    ]);
    let report = check(&part);
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
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
fn small_boolean_perturbations_agree_only_with_finite_trimmed_distance_bounds() {
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
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
    let upper = report
        .sampled_distance_upper
        .as_ref()
        .expect("finite sampled bound");
    assert!(upper >= &Rational::parse("1/10000").unwrap(), "{report:?}");
    assert!(upper <= &options().accuracy, "{report:?}");
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
fn coplanar_boolean_supports_agree_without_unresolved_observations() {
    for kind in ["union", "subtract", "intersect"] {
        for x in ["1", "2"] {
            let report = check(&recipe(vec![
                block("a", "2", "0"),
                block("b", "2", x),
                json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"b"}),
            ]));
            assert_eq!(report.verdict, Verdict::Agree, "{kind} x={x}: {report:?}");
            assert!(report.findings.is_empty(), "{kind} x={x}: {report:?}");
        }
    }
}

#[test]
fn joined_blocks_remove_the_meeting_plane_without_unresolved_observations() {
    let report = check_with_options(
        &recipe(vec![
            json!({"id":"a","op":"block","size":["20","10","5"],"pose":pose("0")}),
            json!({"id":"b","op":"block","size":["5","10","5"],"pose":pose("20")}),
            json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
        ]),
        &Options::default(),
    );
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
    assert!(report.findings.is_empty(), "{report:?}");
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
fn repeated_body_operands_reach_the_set_rules() {
    for kind in ["union", "subtract", "intersect"] {
        let report = check(&recipe(vec![
            block("a", "2", "0"),
            json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"a"}),
        ]));
        assert_eq!(report.verdict, Verdict::Agree, "{kind}: {report:?}");
        assert!(report.findings.is_empty(), "{kind}: {report:?}");
        assert!(report.not_comparable.is_empty(), "{kind}: {report:?}");
    }
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
        wires: Vec::new(),
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
        wires: Vec::new(),
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
fn round_boolean_self_checks_have_finite_trimmed_distance_bounds() {
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
            assert_eq!(report.verdict, Verdict::Agree, "{op} {kind}: {report:?}");
            assert!(
                report
                    .sampled_distance_upper
                    .as_ref()
                    .is_some_and(|u| u <= &options().accuracy),
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
fn feature_curved_extrusions_and_full_partial_revolutions_agree_with_finite_bounds() {
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
        assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
        assert!(
            report
                .sampled_distance_upper
                .as_ref()
                .is_some_and(|u| u <= &options().accuracy),
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

#[test]
fn copied_and_reflected_features_preserve_comparison_certainty() {
    let feature = vec![
        json!({"id":"profile","op":"sketch","plane":pose("2"),"loops":[[
            {"kind":"line","id":"a","to":["1","0"]},{"kind":"line","id":"b","to":["1","1"]},{"kind":"line","id":"c","to":["0","1"]},{"kind":"line","id":"d","to":["0","0"]}
        ]]}),
        json!({"id":"seed","op":"extrude","profile":"profile","distance":"1"}),
    ];
    let cylinder = vec![
        json!({"id":"profile","op":"sketch","plane":pose("2"),"loops":[[
            {"kind":"circle","id":"wall","centre":["0","0"],"radius":"1/2"}
        ]]}),
        json!({"id":"seed","op":"extrude","profile":"profile","distance":"1"}),
    ];
    let revolution = vec![
        json!({"id":"profile","op":"sketch","plane":pose("2"),"loops":[[
            {"kind":"line","id":"a","to":["2","0"]},{"kind":"line","id":"b","to":["2","1"]},{"kind":"line","id":"c","to":["1","1"]},{"kind":"line","id":"d","to":["1","0"]}
        ]]}),
        json!({"id":"seed","op":"revolve","profile":"profile","axis":{"point":["2","0","0"],"direction":["0","1","0"]},"angle":"90"}),
    ];
    for steps in [vec![block("seed", "1", "2")], feature, cylinder, revolution] {
        for operation in [
            json!({"id":"copies","op":"pattern","body":"seed","layout":{"linear":{"vector":["3","0","0"],"count":3}}}),
            json!({"id":"copies","op":"pattern","body":"seed","layout":{"circular":{"axis":{"point":["0","0","0"],"direction":["0","0","5"]},"angle":"90","count":4}}}),
            json!({"id":"copies","op":"mirror","body":"seed","plane":{"point":["0","0","0"],"normal":["2","3","0"]},"keep_original":false}),
        ] {
            let mut steps = steps.clone();
            steps.push(operation);
            let report = check(&recipe(steps));
            assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
            assert!(
                report
                    .sampled_distance_upper
                    .as_ref()
                    .is_some_and(|upper| upper <= &options().accuracy),
                "{report:?}"
            );
        }
    }
}

#[test]
fn a_numeric_refusal_during_copy_comparison_is_not_a_kernel_panic() {
    let part = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0"),"loops":[[
            {"kind":"line","id":"radius","to":["1","0"]},
            {"kind":"line","id":"wall","to":["1","4"]},
            {"kind":"line","id":"top","to":["0","4"]},
            {"kind":"line","id":"axis","to":["0","0"]}
        ]]}),
        json!({"id":"seed","op":"revolve","profile":"profile","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"360"}),
        json!({"id":"copies","op":"pattern","body":"seed","layout":{"linear":{"vector":["4","0","0"],"count":3}}}),
    ]);
    let (report, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        check_with_options(&part, &Options::default())
    });
    assert!(matches!(
        report.verdict,
        Verdict::Agree | Verdict::Unresolved
    ));
    if work.exhausted {
        assert_eq!(report.verdict, Verdict::Unresolved);
        assert!(
            report
                .findings
                .iter()
                .any(|finding| matches!(&finding.reason,
            Reason::Unresolved(reason) if reason.contains("numbers.operation-budget")))
        );
    }
}

#[test]
fn separated_curved_features_keep_finite_gap_bounds_after_near_samples() {
    let part = |radius: &str, z: &str| {
        recipe(vec![
            json!({"id":"pa","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","1"]},"loops":[[
                {"kind":"circle","id":"wall","centre":["0","0"],"radius":"1/5"}
            ]]}),
            json!({"id":"a","op":"extrude","profile":"pa","distance":"2/5"}),
            json!({"id":"pb","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["1","0",z]},"loops":[[
                {"kind":"circle","id":"wall","centre":["0","0"],"radius":radius}
            ]]}),
            json!({"id":"b","op":"extrude","profile":"pb","distance":"2/5"}),
            json!({"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}),
        ])
    };
    // One component supplies genuine near witnesses. The other is both
    // narrower and higher, so none of its old bottom rims is near the new
    // side or cap supports. A failed seed must retain the gap proof path.
    let report = compare(&part("1/5", "0"), &part("1/10", "1/10"), &options());
    assert_eq!(report.verdict, Verdict::Disagree, "{report:?}");
    assert!(
        report.findings.iter().any(|finding| {
            finding.pass == Pass::Forward
                && finding.reason == Reason::Gap
                && finding
                    .gap
                    .is_some_and(|gap| gap.lo() > 0.01 && gap.hi().is_finite())
        }),
        "{report:?}"
    );
}

#[test]
fn motion_shell_comparison_checks_spline_sweep_vertices_and_reports_reverse_gap() {
    for source in [
        include_str!("../../catallax-features/tests/fixtures/spline-sweep-planar.json"),
        include_str!("../../catallax-features/tests/fixtures/spline-sweep-spatial.json"),
    ] {
        let recipe = Recipe::parse(source).unwrap();
        let options = Options::default();
        let (report, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            check_with_options(&recipe, &options)
        });
        assert!(!work.exhausted, "{work:?}");
        assert!(report.refusal.is_none(), "{:?}", report.refusal);
        assert!(report.forward_vertices > 0);
        assert_eq!(report.reverse_samples, 0);
        assert_eq!(report.verdict, Verdict::Unresolved);
        assert_eq!(report.findings.len(), 1, "{:?}", report.findings.first());
        assert_eq!(report.findings[0].pass, Pass::Reverse);
        assert!(report.sampled_distance_upper.unwrap() <= options.accuracy);
    }
}
