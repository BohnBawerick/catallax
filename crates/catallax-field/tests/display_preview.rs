//! Unproved field preview beside the unchanged certified display request.

use catallax_field::{
    Plan,
    display::{Options, PreviewOptions, mesh_plan, preview_plan},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Operation, Pose, Recipe, Step, StepId};

fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn plan(op: Operation) -> Plan {
    Plan::compile(&Recipe::new(
        "preview",
        vec![Step { id: id("part"), op }],
        id("part"),
    ))
    .unwrap()
}
fn block() -> Plan {
    plan(Operation::Block {
        size: [2.into(), 3.into(), 4.into()],
        pose: Pose::identity(),
    })
}
fn cylinder() -> Plan {
    plan(Operation::Cylinder {
        radius: 2.into(),
        height: 5.into(),
        pose: Pose::identity(),
    })
}
fn options(tolerance: Rational) -> Options {
    Options {
        tolerance,
        ..Options::default()
    }
}
fn tenth() -> Rational {
    Rational::from_ratio(1, 10).unwrap()
}

#[test]
fn preview_answers_where_the_certified_request_refuses_on_cells() {
    let plan = block();
    let limited = Options {
        max_cells: 4,
        ..options(tenth())
    };
    let refusal = mesh_plan(&plan, &limited, &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.code, "field.display-budget");
    let preview = preview_plan(&plan, &(&limited).into(), &mut Decider::default()).unwrap();
    assert!(
        preview.unresolved_cells > 0,
        "the unvisited cells are counted"
    );
    assert!(
        preview.triangles.len()
            < mesh_plan(&plan, &options(tenth()), &mut Decider::default())
                .unwrap()
                .triangles
                .len(),
        "a gap stays a gap: nothing fills the unvisited cells"
    );
    // The certified request is unchanged by the preview call.
    assert_eq!(
        mesh_plan(&plan, &limited, &mut Decider::default())
            .unwrap_err()
            .code,
        "field.display-budget"
    );
}

#[test]
fn preview_never_invents_geometry_outside_the_defined_boundary() {
    let plan = block();
    let limited = Options {
        max_cells: 7,
        ..options(tenth())
    };
    let preview = preview_plan(&plan, &(&limited).into(), &mut Decider::default()).unwrap();
    for t in &preview.triangles {
        let p = t.vertices.map(|i| preview.positions[i]);
        let c: Vec<f64> = (0..3)
            .map(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.0)
            .collect();
        let on_boundary = [2.0, 3.0, 4.0]
            .iter()
            .zip(&c)
            .any(|(extent, v)| v.abs() < 1e-9 || (v - extent).abs() < 1e-9);
        assert!(on_boundary, "triangle centre {c:?} is off the boundary");
    }
}

#[test]
fn preview_matches_certified_output_when_nothing_is_unresolved() {
    let plan = block();
    let certified = mesh_plan(&plan, &options(tenth()), &mut Decider::default()).unwrap();
    let preview =
        preview_plan(&plan, &(&options(tenth())).into(), &mut Decider::default()).unwrap();
    assert_eq!(preview.unresolved_cells, 0);
    assert_eq!(preview.positions, certified.positions);
    assert_eq!(preview.triangles, certified.triangles);
    assert_eq!(preview.cells, certified.cells);
}

#[test]
fn a_coarser_preview_target_differs_from_the_certified_mesh() {
    let plan = cylinder();
    let certified = mesh_plan(
        &plan,
        &options(Rational::from_ratio(1, 100).unwrap()),
        &mut Decider::default(),
    )
    .unwrap();
    let preview = preview_plan(
        &plan,
        &PreviewOptions {
            tolerance: Rational::from_ratio(1, 2).unwrap(),
            ..PreviewOptions::default()
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(!preview.triangles.is_empty());
    assert!(preview.triangles.len() < certified.triangles.len());
}

#[test]
fn an_exhausted_allowance_gives_no_partial_preview() {
    let plan = cylinder();
    let (result, work) = with_operation_budget(50, || {
        preview_plan(&plan, &PreviewOptions::default(), &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn the_vertex_limit_still_refuses_the_preview() {
    let preview = preview_plan(
        &block(),
        &PreviewOptions {
            max_vertices: 4,
            ..(&options(tenth())).into()
        },
        &mut Decider::default(),
    );
    assert_eq!(preview.unwrap_err().code, "field.display-budget");
}

/// Diagnostic only: `cargo test -p catallax-field --test display_preview -- --ignored --nocapture`.
#[test]
#[ignore = "diagnostic timing, not a speed qualification"]
fn preview_against_certified_time_on_public_fixtures() {
    let p1 = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let plate = Recipe::new(
        "plate",
        p1.steps
            .into_iter()
            .filter(|s| s.id == id("plate-profile") || s.id == id("plate"))
            .collect(),
        id("plate"),
    );
    for (name, plan) in [
        ("block", block()),
        ("cylinder", cylinder()),
        ("p1-plate", Plan::compile(&plate).unwrap()),
    ] {
        let start = std::time::Instant::now();
        let certified = mesh_plan(&plan, &Options::default(), &mut Decider::default());
        let certified_time = start.elapsed();
        let start = std::time::Instant::now();
        let preview = preview_plan(&plan, &PreviewOptions::default(), &mut Decider::default());
        println!(
            "{name}: certified {certified_time:?} ({} triangles), preview {:?} ({} triangles, {} gaps)",
            certified.map_or(0, |m| m.triangles.len()),
            start.elapsed(),
            preview.as_ref().map_or(0, |m| m.triangles.len()),
            preview.as_ref().map_or(0, |m| m.unresolved_cells),
        );
    }
}
