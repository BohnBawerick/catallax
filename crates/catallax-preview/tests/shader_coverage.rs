//! Target-part and sweep coverage through the actual display batch evaluator.
use catallax_field::{Plan, Workspace};
use catallax_numbers::with_operation_budget;
use catallax_preview::{
    batch::{Evaluator, Hardware, MAX_DISPATCH, Options},
    gpu::Gpu,
};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn compile(value: Value) -> Plan {
    Plan::compile(&Recipe::parse(&value.to_string()).unwrap()).unwrap()
}

fn cases() -> Vec<(String, Plan)> {
    let mut cases = Vec::new();
    for (name, text) in [
        ("P1", include_str!("../../../parts/P1-bracket.json")),
        (
            "P1-drafted",
            include_str!("../../../parts/P1-drafted-bracket.json"),
        ),
        ("P2", include_str!("../../../parts/P2-collar.json")),
        ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
        ("P4", include_str!("../../../parts/P4-shelled-tray.json")),
    ] {
        cases.push((
            name.to_string(),
            Plan::compile(&Recipe::parse(text).unwrap()).unwrap(),
        ));
    }
    for (name, text) in [
        (
            "annular-draft",
            include_str!("../../catallax-features/tests/fixtures/draft-annular/both.json"),
        ),
        (
            "concave-draft",
            include_str!("../../catallax-features/tests/fixtures/concave-draft/R-P45.json"),
        ),
        (
            "oblique-draft-fixture",
            include_str!("../../catallax-features/tests/fixtures/draft-oblique.json"),
        ),
    ] {
        cases.push((
            name.into(),
            Plan::compile(&Recipe::parse(text).unwrap()).unwrap(),
        ));
    }
    let mut p5: Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    p5["steps"].as_array_mut().unwrap().last_mut().unwrap()["radius"] = json!("1/2");
    cases.push(("P5-radius-half-edit".into(), compile(p5)));
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    for (name, base, wall, pull) in [
        (
            "planar-draft",
            json!({"id":"base","op":"block","size":["8","8","4"],"pose":pose}),
            "x_max",
            ["0", "0", "1"],
        ),
        (
            "round-draft",
            json!({"id":"base","op":"cylinder","radius":"8","height":"4","pose":pose}),
            "wall",
            ["0", "0", "1"],
        ),
    ] {
        cases.push((name.into(), compile(json!({"format":"catallax-recipe/3","name":name,"steps":[base,
            {"id":"drafted","op":"draft","body":"base","faces":[{"face":{"step":"base","role":{"side":{"segment":wall}}}}],
             "pull":pull,"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"45"}
        ],"result":"drafted"}))));
    }
    for angle in ["90", "120", "180", "240", "270"] {
        let mut ring: Value =
            serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
        ring["steps"].as_array_mut().unwrap().truncate(2);
        ring["steps"][1]["angle"] = json!(angle);
        ring["result"] = json!("ring");
        cases.push((format!("revolve-{angle}"), compile(ring)));
    }
    cases.push(("through-cut".into(), compile(json!({"format":"catallax-recipe/4","name":"through cut","steps":[
        {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","-5"]}},
        {"id":"profile","op":"sketch","plane":pose,"loops":[[{"kind":"circle","id":"circle","centre":["5","5"],"radius":"1"}]]},
        {"id":"cut","op":"extrude_through","target":"stock","profile":"profile"}
    ],"result":"cut"}))));
    cases
}

#[test]
fn authored_p5_keeps_its_compiler_refusal_before_shader_translation() {
    let recipe = Recipe::parse(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let (answer, work) = with_operation_budget(100_000_000, || Plan::compile(&recipe));
    let refusal = answer.unwrap_err();
    assert_eq!(refusal.code, "fillet.unsupported-overflow");
    assert!(!work.exhausted);
    eprintln!(
        "authored P5 radius=2: {}: {}; work={}",
        refusal.code, refusal.detail, work.used
    );
}

fn points() -> Vec<[f64; 3]> {
    let mut state = 0x15_2026_1113u64;
    (0..20_000)
        .map(|_| {
            std::array::from_fn(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                // Binary32 inputs keep GPU and processor coordinates identical.
                f64::from(((state >> 40) as f32 / (1u32 << 24) as f32) * 80.0 - 30.0)
            })
        })
        .collect()
}

#[test]
fn target_and_sweep_shaders_validate_without_a_device_or_proof_work() {
    for (name, plan) in cases() {
        let before = plan
            .fast([2.0, 3.0, 1.0], &mut Workspace::default())
            .unwrap();
        let (source, work) = with_operation_budget(0, || plan.unproved_preview_wgsl().unwrap());
        assert_eq!(work.used, 0, "{name}");
        assert!(!work.exhausted, "{name}");
        assert_eq!(source, plan.unproved_preview_wgsl().unwrap(), "{name}");
        assert_eq!(
            before.value.to_bits(),
            plan.fast([2.0, 3.0, 1.0], &mut Workspace::default())
                .unwrap()
                .value
                .to_bits()
        );
        if let Some(dir) = std::env::var_os("CATALLAX_SHADER_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&dir).join(format!("{name}.wgsl")),
                &source,
            )
            .unwrap();
        }
        let source = format!("{source}\n{}", include_str!("../src/batch.wgsl"));
        let module = naga::front::wgsl::parse_str(&source)
            .unwrap_or_else(|e| panic!("{name}: {}", e.emit_to_string(&source)));
        naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        )
        .validate(&module)
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
    }
}

#[test]
fn target_gpu_values_repeat_and_match_processor_signs_when_available() {
    let gpu = match Gpu::new(false) {
        Ok(gpu) => gpu,
        Err(reason) => {
            assert!(
                std::env::var_os("CATALLAX_REQUIRE_GPU").is_none(),
                "{reason}"
            );
            eprintln!("GPU coverage skipped: {reason}");
            return;
        }
    };
    eprintln!("shader coverage adapter: {}", gpu.name());
    let points = points();
    for (name, plan) in cases() {
        let evaluate = |gpu, workers, dispatch| {
            let evaluator = Evaluator::new(&plan, gpu, Options { workers, dispatch }).unwrap();
            if gpu.is_some() {
                assert!(
                    matches!(evaluator.hardware(), Hardware::GraphicsCard { .. }),
                    "{name}: {:?}",
                    evaluator.hardware()
                );
            }
            evaluator.points(&points).unwrap().values
        };
        let processor = evaluate(None, 1, MAX_DISPATCH);
        let bits = |values: &[f64]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        for workers in [2, 4] {
            assert_eq!(
                bits(&processor),
                bits(&evaluate(None, workers, MAX_DISPATCH)),
                "{name} workers={workers}"
            );
        }
        let values = evaluate(Some(&gpu), 1, MAX_DISPATCH);
        assert_eq!(
            bits(&values),
            bits(&evaluate(Some(&gpu), 1, MAX_DISPATCH)),
            "{name} repeat"
        );
        assert_eq!(
            bits(&values),
            bits(&evaluate(Some(&gpu), 4, 1003)),
            "{name} dispatch"
        );
        let mut near = 0;
        let mut max_difference = 0.0f64;
        let mut farthest = 0;
        for (i, (&value, &reference)) in values.iter().zip(&processor).enumerate() {
            if (value - reference).abs() > max_difference {
                max_difference = (value - reference).abs();
                farthest = i;
            }
            if (value < 0.0, value > 0.0) != (reference < 0.0, reference > 0.0) {
                assert!(
                    reference.abs() <= 1e-4,
                    "{name} point {:?}: GPU={value}, processor={reference}",
                    points[i]
                );
                near += 1;
            }
        }
        eprintln!(
            "{name}: points={} near-boundary sign differences={near} max value difference={max_difference}; repeat, dispatch and 1/2/4 processor bits agree",
            points.len()
        );
        eprintln!(
            "{name} largest difference at {:?}: GPU={}, processor={}",
            points[farthest], values[farthest], processor[farthest]
        );
    }
}

#[test]
fn drafted_p1_keeps_the_nearest_planar_foot_on_the_gpu() {
    let gpu = match Gpu::new(false) {
        Ok(gpu) => gpu,
        Err(reason) => {
            assert!(
                std::env::var_os("CATALLAX_REQUIRE_GPU").is_none(),
                "{reason}"
            );
            eprintln!("GPU planar-foot control skipped: {reason}");
            return;
        }
    };
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P1-drafted-bracket.json")).unwrap(),
    )
    .unwrap();
    let points = [
        [9.944385528564453, 11.361885070800781, 16.847496032714844],
        [32.99707794189453, 11.38400650024414, 22.54297637939453],
    ];
    // Each nearest material point is on the flange's y=5 wall, inside every
    // other trim. Retesting that plane rounded the foot outside it.
    let values = Evaluator::new(&plan, Some(&gpu), Options::default())
        .unwrap()
        .points(&points)
        .unwrap()
        .values;
    for (point, value) in points.into_iter().zip(values) {
        assert!(
            (value - (point[1] - 5.0)).abs() < 1e-5,
            "point {point:?}: GPU={value}"
        );
    }
}
