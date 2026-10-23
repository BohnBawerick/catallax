//! Definition-native display coverage, boolean cancellation and resource refusal.
use catallax_field::{
    Plan,
    display::{Options, mesh_plan},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{BooleanKind, Operation, Pose, Recipe, Step, StepId};

fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn block() -> Operation {
    Operation::Block {
        size: [2.into(), 3.into(), 4.into()],
        pose: Pose::identity(),
    }
}
fn plan(steps: Vec<Step>, result: &str) -> Plan {
    Plan::compile(&Recipe::new("display", steps, id(result))).unwrap()
}
fn options() -> Options {
    Options {
        tolerance: Rational::from_ratio(1, 10).unwrap(),
        ..Options::default()
    }
}

#[test]
fn displayed_box_stays_on_its_boundary_and_covers_every_side() {
    let plan = plan(
        vec![Step {
            id: id("box"),
            op: block(),
        }],
        "box",
    );
    let mesh = mesh_plan(&plan, &options(), &mut Decider::default()).unwrap();
    let mut sides = [false; 6];
    for t in &mesh.triangles {
        let p = t.vertices.map(|i| mesh.positions[i]);
        let centre = std::array::from_fn::<_, 3, _>(|i| (p[0][i] + p[1][i] + p[2][i]) / 3.0);
        for (axis, extent) in [2.0, 3.0, 4.0].into_iter().enumerate() {
            assert!(centre[axis] >= -0.1 && centre[axis] <= extent + 0.1);
            if centre[axis].abs() < 0.1 {
                sides[axis * 2] = true;
            }
            if (centre[axis] - extent).abs() < 0.1 {
                sides[axis * 2 + 1] = true;
            }
        }
        assert!(
            [
                centre[0].abs(),
                (centre[0] - 2.0).abs(),
                centre[1].abs(),
                (centre[1] - 3.0).abs(),
                centre[2].abs(),
                (centre[2] - 4.0).abs()
            ]
            .into_iter()
            .any(|d| d <= 0.1)
        );
    }
    assert_eq!(
        sides, [true; 6],
        "a missing box side is not a complete display"
    );
}
#[test]
fn cancelled_solid_has_no_display_and_budget_failure_publishes_no_mesh() {
    let plan = plan(
        vec![
            Step {
                id: id("box"),
                op: block(),
            },
            Step {
                id: id("empty"),
                op: Operation::Boolean {
                    kind: BooleanKind::Subtract,
                    target: id("box"),
                    tool: id("box"),
                },
            },
        ],
        "empty",
    );
    let mesh = mesh_plan(&plan, &options(), &mut Decider::default()).unwrap();
    assert!(mesh.positions.is_empty() && mesh.triangles.is_empty());
    assert_eq!(mesh.max_deviation, Rational::zero());
    let plan = Plan::compile(&Recipe::new(
        "solid",
        vec![Step {
            id: id("box"),
            op: block(),
        }],
        id("box"),
    ))
    .unwrap();
    let limited = Options {
        max_cells: 1,
        ..options()
    };
    assert_eq!(
        mesh_plan(&plan, &limited, &mut Decider::default())
            .unwrap_err()
            .code,
        "field.display-budget"
    );
}
#[test]
fn native_mesh_is_identical_across_independent_workers() {
    let plan = plan(
        vec![Step {
            id: id("box"),
            op: block(),
        }],
        "box",
    );
    let expected = mesh_plan(&plan, &options(), &mut Decider::default()).unwrap();
    for workers in [1, 2, 4] {
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for _ in 0..workers {
                let plan = &plan;
                handles.push(
                    scope.spawn(move || {
                        mesh_plan(plan, &options(), &mut Decider::default()).unwrap()
                    }),
                );
            }
            for handle in handles {
                assert_eq!(handle.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn p1_plate_bottom_aligned_edges_preserve_source_face_area() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    // Isolate the unchanged creating face; full P1's work refusal is measured
    // through the CLI, not pinned as a permanent expected outcome.
    let recipe = Recipe::new(
        "P1 plate source",
        recipe
            .steps
            .into_iter()
            .filter(|s| s.id == id("plate-profile") || s.id == id("plate"))
            .collect(),
        id("plate"),
    );
    let plan = Plan::compile(&recipe).unwrap();
    let mesh = mesh_plan(&plan, &Options::default(), &mut Decider::default()).unwrap();
    let selected = plan
        .select(&catallax_recipe::Origin {
            step: id("plate"),
            role: catallax_recipe::Role::Side {
                segment: catallax_recipe::ElementId::new("bottom").unwrap(),
                path: None,
            },
        })
        .unwrap();
    let mut area = 0.0;
    for triangle in &mesh.triangles {
        if triangle.face != selected {
            continue;
        }
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        for p in [a, b, c] {
            assert!(p[1].abs() <= 0.01);
            assert!(p[0] >= -0.01 && p[0] <= 60.01 && p[2] >= -0.01 && p[2] <= 5.01);
        }
        area += ((b[0] - a[0]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[0] - a[0])).abs() * 0.5;
    }
    assert!(
        (area - 300.0).abs() <= 2.0,
        "lost or double-counted plate-bottom area: {area}"
    );
}

#[test]
fn observations_preserve_the_mesh_at_its_exact_work_allowance() {
    use catallax_field::display::{Diagnostics, mesh_plan_reported};
    use catallax_numbers::with_operation_budget;
    let plan = plan(
        vec![Step {
            id: id("box"),
            op: block(),
        }],
        "box",
    );
    let (baseline, work) = with_operation_budget(100_000_000, || {
        mesh_plan(&plan, &options(), &mut Decider::default())
    });
    let baseline = baseline.unwrap();
    let mut diagnostics = Diagnostics::default();
    let (observed, observed_work) = with_operation_budget(work.used, || {
        mesh_plan_reported(&plan, &options(), &mut Decider::default(), &mut diagnostics)
    });
    assert!(
        !observed_work.exhausted,
        "observing work must not consume geometry's allowance"
    );
    assert_eq!(observed.unwrap(), baseline);
    assert_eq!(observed_work.used, work.used);
    assert_eq!(diagnostics.total.work, work.used);
    assert_eq!(
        diagnostics.common.work
            + diagnostics
                .supports
                .iter()
                .map(|s| s.total.work)
                .sum::<u64>(),
        work.used
    );
    assert_eq!(diagnostics.total.cells, baseline.cells);
}

#[test]
fn translated_cylinder_displays_its_boundary_without_quadric_cancellation() {
    use catallax_numbers::with_operation_budget;
    let mut pose = Pose::identity();
    pose.translation.0 = [1_000_000.into(), 1_000_000.into(), 0.into()];
    let plan = plan(
        vec![Step {
            id: id("cylinder"),
            op: Operation::Cylinder {
                radius: 3.into(),
                height: 7.into(),
                pose,
            },
        }],
        "cylinder",
    );
    let (result, work) = with_operation_budget(100_000_000, || {
        mesh_plan(&plan, &Options::default(), &mut Decider::default())
    });
    let mesh = result.unwrap();
    assert!(!work.exhausted);
    let bound = mesh.max_deviation.to_f64_nearest() + 1e-8;
    for p in &mesh.positions {
        let rho = ((p[0] - 1_000_000.0).powi(2) + (p[1] - 1_000_000.0).powi(2)).sqrt();
        let z = p[2];
        let distance = if rho <= 3.0 && (0.0..=7.0).contains(&z) {
            (3.0 - rho).min(z).min(7.0 - z)
        } else {
            ((rho - 3.0).max(0.0).powi(2) + (-z).max(z - 7.0).max(0.0).powi(2)).sqrt()
        };
        assert!(
            distance <= bound,
            "display point outside certified cylinder bound: {p:?}, {distance}"
        );
    }
    let area = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.vertices.map(|i| mesh.positions[i]);
            let u = std::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
            let v = std::array::from_fn::<_, 3, _>(|i| c[i] - a[i]);
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            cross.into_iter().map(|x| x * x).sum::<f64>().sqrt() * 0.5
        })
        .sum::<f64>();
    assert!(
        (area - 60.0 * std::f64::consts::PI).abs() < 2.0,
        "missing or duplicated cylinder boundary: {area}"
    );
}
