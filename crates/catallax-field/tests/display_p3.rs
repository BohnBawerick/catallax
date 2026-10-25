//! The unchanged public P3 recipe: full native display, independent domains,
//! outward winding and surviving origin sets at the scoreboard allowance.
use std::{collections::BTreeSet, time::Instant};

use catallax_field::{
    Plan,
    display::{Diagnostics, Options, mesh_plan_reported},
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role};

fn radius(p: [f64; 3]) -> f64 {
    ((p[0] - 20.0).powi(2) + (p[1] - 15.0).powi(2)).sqrt()
}

// Independent analytic controls in the public recipe's model frame. The field
// mesher's whole-cell certificates remain responsible for complete coverage.
fn control(origin: &Origin, p: [f64; 3], tolerance: f64) -> [f64; 3] {
    let [x, y, z] = p;
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-8, "{origin:?}: {p:?}: {a} != {b}");
    assert!(x >= -tolerance && x <= 40.0 + tolerance && y >= -tolerance && y <= 30.0 + tolerance);
    match (origin.step.as_str(), &origin.role) {
        ("block", Role::Start) => {
            near(z, 0.0);
            [0.0, 0.0, -1.0]
        }
        ("block", Role::End) => {
            near(z, 10.0);
            assert!(
                x >= 1.0 - tolerance
                    && x <= 39.0 + tolerance
                    && y >= 1.0 - tolerance
                    && y <= 29.0 + tolerance
            );
            assert!(radius(p) >= 10.0 - tolerance);
            [0.0, 0.0, 1.0]
        }
        ("block", Role::Side { segment, .. }) => {
            assert!(z >= -tolerance && z <= 9.0 + tolerance);
            match segment.as_str() {
                "x_min" => {
                    near(x, 0.0);
                    [-1.0, 0.0, 0.0]
                }
                "x_max" => {
                    near(x, 40.0);
                    [1.0, 0.0, 0.0]
                }
                "y_min" => {
                    near(y, 0.0);
                    [0.0, -1.0, 0.0]
                }
                "y_max" => {
                    near(y, 30.0);
                    [0.0, 1.0, 0.0]
                }
                _ => panic!("unexpected origin {origin:?}"),
            }
        }
        ("boss", Role::End) => {
            near(z, 22.0);
            assert!(radius(p) <= 8.0 + tolerance);
            [0.0, 0.0, 1.0]
        }
        ("boss", Role::Side { .. }) => {
            near(radius(p), 8.0);
            assert!(z >= 12.0 - tolerance && z <= 22.0 + tolerance);
            [x - 20.0, y - 15.0, 0.0]
        }
        ("junction", Role::Blend { .. }) => {
            let rho = radius(p);
            near((rho - 10.0).powi(2) + (z - 12.0).powi(2), 4.0);
            assert!(
                rho >= 8.0 - tolerance
                    && rho <= 10.0 + tolerance
                    && z >= 10.0 - tolerance
                    && z <= 12.0 + tolerance
            );
            [
                (10.0 - rho) * (x - 20.0) / rho,
                (10.0 - rho) * (y - 15.0) / rho,
                12.0 - z,
            ]
        }
        ("rounded", Role::Blend { edge }) => {
            let side = edge
                .between()
                .into_iter()
                .find_map(|o| {
                    if let Role::Side { segment, .. } = &o.role {
                        Some(segment.as_str())
                    } else {
                        None
                    }
                })
                .unwrap();
            let (axis, centre, distance, other, extent) = match side {
                "x_min" => (0, 1.0, x, y, 30.0),
                "x_max" => (0, 39.0, 40.0 - x, y, 30.0),
                "y_min" => (1, 1.0, y, x, 40.0),
                "y_max" => (1, 29.0, 30.0 - y, x, 40.0),
                _ => panic!("unexpected blend {origin:?}"),
            };
            near((p[axis] - centre).powi(2) + (z - 9.0).powi(2), 1.0);
            assert!(
                distance >= -tolerance
                    && distance <= 1.0 + tolerance
                    && z >= 9.0 - tolerance
                    && z <= 10.0 + tolerance
            );
            assert!(other >= distance - tolerance && other <= extent - distance + tolerance);
            let mut normal = [0.0, 0.0, z - 9.0];
            normal[axis] = p[axis] - centre;
            normal
        }
        _ => panic!("buried or unexpected origin {origin:?}"),
    }
}

#[test]
fn complete_p3_native_display_fits_scoreboard_allowance() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let ((plan, answer, diagnostics), work) = with_operation_budget(1_000_000, || {
        let mut dec = Decider::default();
        let started = Instant::now();
        let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
        let build = started.elapsed();
        let started = Instant::now();
        let mut diagnostics = Diagnostics::default();
        let answer = mesh_plan_reported(&plan, &Options::default(), &mut dec, &mut diagnostics);
        eprintln!(
            "P3 field build={build:?}, display={:?}, display work={}",
            started.elapsed(),
            diagnostics.total.work
        );
        (plan, answer, diagnostics)
    });
    let mesh = answer.unwrap();
    assert!(!work.exhausted);
    assert!(mesh.max_deviation <= Options::default().tolerance);
    assert_eq!(diagnostics.total.cells, mesh.cells);
    let (built, brep_work) = with_operation_budget(1_000_000, || {
        catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap()
    });
    assert!(!brep_work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    let expected: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).iter().cloned())
        .collect();
    let actual: BTreeSet<_> = mesh
        .triangles
        .iter()
        .map(|t| plan.faces()[t.face.index()].origin.clone())
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(actual.len(), 13);
    let tolerance = mesh.max_deviation.to_f64_nearest() + 1e-8;
    let mut bottom_area = 0.0;
    let mut top_area = 0.0;
    for triangle in &mesh.triangles {
        let origin = &plan.faces()[triangle.face.index()].origin;
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        for p in [a, b, c] {
            control(origin, p, tolerance);
        }
        let u: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
        let v: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let normals = [a, b, c].map(|p| control(origin, p, tolerance));
        assert!(
            (0..3)
                .map(|i| n[i] * (normals[0][i] + normals[1][i] + normals[2][i]))
                .sum::<f64>()
                > 0.0,
            "inward or degenerate triangle {origin:?}: {a:?} {b:?} {c:?}"
        );
        if origin.step.as_str() == "block" {
            if origin.role == Role::Start {
                bottom_area += n[2].abs() * 0.5;
            }
            if origin.role == Role::End {
                top_area += n[2].abs() * 0.5;
            }
        }
    }
    assert!((bottom_area - 1200.0).abs() < 1e-8);
    assert!(
        (top_area - (38.0 * 28.0 - 100.0 * std::f64::consts::PI)).abs() < 2.0,
        "top area={top_area}"
    );
}
