//! Fixed binary64 inputs and complete output bits for the R7 hosted comparison.

use catallax_field::{Plan, PreviewHardware, UnprovedBatch, Workspace};
use catallax_recipe::{BooleanKind, Operation, Pose, Recipe, Step, StepId};
use std::fmt::Write;

fn points() -> Vec<[f64; 3]> {
    let mut points = vec![[0.0; 3], [-0.0; 3], [40.0, 30.0, 22.0]];
    let mut state = 7u64;
    for _ in 0..256 {
        points.push(std::array::from_fn(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            // Exactly representable dyadic coordinates, with no library maths.
            ((state >> 48) as f64 - 32768.0) / 512.0
        }));
    }
    for x in [
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1e-150,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(0x7ff8_0000_0000_0042),
        f64::from_bits(0xfff0_0000_0000_0001),
    ] {
        for axis in 0..3 {
            let mut point = [1.0; 3];
            point[axis] = x;
            points.push(point);
        }
    }
    points
}

#[test]
fn fixed_processor_bits_and_refusals() {
    let points = points();
    assert_eq!(points.len(), 286);
    let mut transcript = String::from("catallax-r7-bits/1\n");
    for (name, source) in [
        ("P1", include_str!("../../../parts/P1-bracket.json")),
        ("P2", include_str!("../../../parts/P2-collar.json")),
        ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
        ("P4", include_str!("../../../parts/P4-shelled-tray.json")),
    ] {
        let plan = Plan::compile(&Recipe::parse(source).unwrap()).unwrap();
        let mut workspace = Workspace::default();
        let mut batch = UnprovedBatch::default();
        for (i, &point) in points.iter().enumerate() {
            let input = point.map(f64::to_bits);
            let scalar = plan.unproved_preview(point, &mut workspace);
            let fast = plan.fast(point, &mut workspace);
            let record = match &scalar {
                Ok(value) => {
                    assert_eq!(value.hardware, PreviewHardware::Processor);
                    assert!(!value.value.is_nan());
                    format!(
                        "ok {:016x} {}",
                        value.value.to_bits(),
                        value.true_distance_formula
                    )
                }
                Err(error) => format!("refused {} {}", error.code, error.detail),
            };
            let fast_record = match fast {
                Ok(value) => format!(
                    "ok {:016x} {}",
                    value.value.to_bits(),
                    value.true_distance_formula
                ),
                Err(error) => format!("refused {} {}", error.code, error.detail),
            };
            writeln!(
                transcript,
                "{name} {i} {input:016x?} preview {record} fast {fast_record}"
            )
            .unwrap();
            for workers in [1, 2, 4] {
                // Include successes before and after a failing input. Reuse storage.
                let result =
                    plan.unproved_preview_batch(&[[1.0; 3], point, [2.0; 3]], workers, &mut batch);
                match (&scalar, result) {
                    (Ok(a), Ok(values)) => {
                        assert_eq!(values.len(), 3);
                        assert_eq!(values[1].value.to_bits(), a.value.to_bits());
                        assert_eq!(values[1].true_distance_formula, a.true_distance_formula);
                        assert_eq!(values[1].hardware, a.hardware);
                    }
                    (Err(a), Err(b)) => {
                        assert_eq!(a.code, b.code);
                        assert_eq!(a.detail, b.detail);
                    }
                    pair => panic!("{name} point {i}, {workers} workers: {pair:?}"),
                }
            }
        }
    }
    assert_eq!(transcript.lines().count(), 1 + 4 * points.len());
    if let Some(path) = std::env::var_os("CATALLAX_R7_BITS") {
        std::fs::write(path, transcript).unwrap();
    }
}

#[test]
fn nonfinite_policy_has_typed_boundaries_and_an_empty_sentinel() {
    let plan = Plan::compile(
        &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
    )
    .unwrap();
    let mut workspace = Workspace::default();
    for point in points()
        .into_iter()
        .filter(|p| p.iter().any(|x| !x.is_finite()))
    {
        assert_eq!(
            plan.unproved_preview(point, &mut workspace)
                .unwrap_err()
                .code,
            "field.non-finite-point"
        );
    }
    assert_eq!(
        plan.unproved_preview([f64::MAX; 3], &mut workspace)
            .unwrap_err()
            .code,
        "field.fast-range"
    );

    let id = StepId::new("base").unwrap();
    let empty = StepId::new("empty").unwrap();
    let plan = Plan::compile(&Recipe::new(
        "empty",
        vec![
            Step {
                id: id.clone(),
                op: Operation::Block {
                    size: [8.into(), 8.into(), 8.into()],
                    pose: Pose::identity(),
                },
            },
            Step {
                id: empty.clone(),
                op: Operation::Boolean {
                    kind: BooleanKind::Subtract,
                    target: id.clone(),
                    tool: id,
                },
            },
        ],
        empty,
    ))
    .unwrap();
    for point in [[0.0; 3], [f64::MAX; 3]] {
        assert_eq!(
            plan.unproved_preview(point, &mut workspace)
                .unwrap()
                .value
                .to_bits(),
            f64::INFINITY.to_bits()
        );
    }
    assert_eq!(
        plan.unproved_preview([f64::NAN; 3], &mut workspace)
            .unwrap_err()
            .code,
        "field.non-finite-point"
    );
}
