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
fn p1_plate_bottom_keeps_its_area_within_the_native_display_budget() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
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
