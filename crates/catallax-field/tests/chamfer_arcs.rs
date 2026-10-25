//! Field-native finite arc material and origin queries, with analytic controls.
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn source(name: &str) -> Value {
    serde_json::from_str(match name {
        "convex" => include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/convex.json"),
        "concave" => {
            include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/concave.json")
        }
        "chord" => include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/chord.json"),
        "seam" => include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/seam.json"),
        "major" => include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/major.json"),
        "major-seam" => {
            include_str!("../../catallax-features/tests/fixtures/chamfer-arcs/major-seam.json")
        }
        _ => unreachable!(),
    })
    .unwrap()
}
fn recipe(name: &str, distance: &str, other: &str, mirror: bool) -> Recipe {
    let mut value = source(name);
    let bevel = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
    bevel["distance"] = json!(distance);
    bevel["other_distance"] = json!(other);
    if mirror {
        value["steps"]
            .as_array_mut()
            .unwrap()
            .push(json!({"id":"reflected","op":"mirror","body":"bevel",
            "plane":{"point":["0","0","0"],"normal":["3","0","0"]},"keep_original":false}));
        value["result"] = json!("reflected");
    }
    Recipe::parse(&value.to_string()).unwrap()
}
fn point(p: [&str; 3], mirror: bool) -> RPoint3 {
    let [mut x, y, z] = p.map(|v| Rational::parse(v).unwrap());
    if mirror {
        x = -x;
    }
    RPoint3::new(x, y, z)
}
#[test]
fn exact_material_controls_cover_contacts_cut_ends_and_excluded_arcs() {
    use Location::{Boundary, Inside, Outside};
    for (name, cases) in [
        (
            "convex",
            vec![
                (["21/5", "28/5", "23/2"], Boundary),
                (["39/10", "26/5", "23/2"], Inside),
                (["9/2", "6", "23/2"], Outside),
                (["-21/5", "28/5", "23/2"], Outside),
                (["6", "0", "12"], Boundary),
                (["8", "0", "11"], Boundary),
                (["0", "8", "12"], Outside),
            ],
        ),
        (
            "chord",
            vec![
                (["5", "0", "12"], Boundary),
                (["24/5", "18/5", "12"], Boundary),
                (["24/5", "32/5", "11"], Boundary),
                (["4", "3", "11"], Outside),
            ],
        ),
        (
            "seam",
            vec![
                (["7", "0", "23/2"], Boundary),
                (["7", "-1", "23/2"], Outside),
                (["6", "-1", "23/2"], Inside),
                (["-7", "0", "23/2"], Outside),
            ],
        ),
        (
            "major",
            vec![
                (["-7", "0", "23/2"], Boundary),
                (["0", "7", "23/2"], Boundary),
                (["21/5", "-28/5", "23/2"], Outside),
            ],
        ),
    ] {
        for mirror in [false, true] {
            let plan =
                Plan::compile_checked(&recipe(name, "2", "1", mirror), &mut decider()).unwrap();
            for (p, expected) in &cases {
                let got = plan
                    .proven_checked(
                        &point(*p, mirror),
                        &mut Workspace::default(),
                        &mut decider(),
                    )
                    .unwrap();
                assert_eq!(got.location, *expected, "{name} {mirror} {p:?}");
            }
        }
    }
}
#[test]
fn equal_unequal_and_edited_strips_resolve_the_same_named_geometry() {
    for name in ["convex", "chord", "seam", "major", "major-seam"] {
        for mirror in [false, true] {
            let mut old = None;
            for (distance, other) in [("2", "1"), ("2", "2"), ("1", "1/2")] {
                let plan =
                    Plan::compile_checked(&recipe(name, distance, other, mirror), &mut decider())
                        .unwrap();
                let mut origins: Vec<_> = plan.faces().iter().map(|f| f.origin.clone()).collect();
                origins.sort();
                origins.dedup();
                if let Some(old) = old {
                    assert_eq!(origins, old);
                }
                old = Some(origins);
                let bevels: Vec<_> = plan.faces().iter().filter(|f| {
                    if mirror { matches!(&f.origin.role,Role::Copy { of,.. } if matches!(of.role,Role::Bevel {..})) }
                    else { matches!(f.origin.role,Role::Bevel {..}) }
                }).collect();
                assert_eq!(bevels.len(), 1, "{name} {mirror}");
                let origin = &bevels[0].origin;
                let selected = plan.select(origin).unwrap();
                // Query an interior conical point; the transverse ends have
                // several incident faces and are checked by material controls.
                let z = if other == "2" {
                    "11"
                } else if other == "1/2" {
                    "47/4"
                } else {
                    "23/2"
                };
                let radius = if distance == "1" { "15/2" } else { "7" };
                let radius = Rational::parse(radius).unwrap();
                let (x, y) = if matches!(name, "chord" | "seam") {
                    (4, 3)
                } else {
                    (3, 4)
                };
                let mut x = radius.clone() * Rational::from_ratio(x, 5).unwrap();
                if mirror {
                    x = -x;
                }
                let sample = RPoint3::new(
                    x,
                    radius * Rational::from_ratio(y, 5).unwrap(),
                    Rational::parse(z).unwrap(),
                );
                let query = plan
                    .face_query_checked(
                        selected,
                        &sample,
                        &mut Workspace::default(),
                        &mut decider(),
                    )
                    .unwrap();
                assert_eq!(
                    query.state,
                    FaceState::OnFace,
                    "{name} {distance} {other} {mirror}"
                );
                let Direction::Unique(normal) = query.direction else {
                    panic!("{query:?}")
                };
                assert!(if mirror {
                    normal[0].hi() < 0.0
                } else {
                    normal[0].lo() > 0.0
                });
                assert!(normal[1].lo() > 0.0 && normal[2].lo() > 0.0);
                assert_ne!(
                    plan.face_query_checked(
                        selected,
                        &point(["-7", "-7", z], mirror),
                        &mut Workspace::default(),
                        &mut decider()
                    )
                    .unwrap()
                    .state,
                    FaceState::OnFace
                );
            }
        }
    }
}

#[test]
fn finite_concave_arcs_remain_typed_refusals_for_edits_and_reflections() {
    for mirror in [false, true] {
        for (plane, cylinder) in [("2", "1"), ("2", "2"), ("1", "1/2")] {
            let refusal =
                Plan::compile_checked(&recipe("concave", plane, cylinder, mirror), &mut decider())
                    .unwrap_err();
            assert_eq!(refusal.code, "chamfer.unsupported-concave-arc");
        }
    }
}

#[test]
fn circular_cap_charts_preserve_all_obj_groups_at_scoreboard_allowance() {
    use catallax_field::display::{Diagnostics, Options, mesh_plan_reported};
    use catallax_numbers::with_operation_budget;
    use std::collections::BTreeSet;

    // An independent point-set construction of the same quarter cylinder.
    // The box extends past both caps and the circle's tangent planes, so its
    // only active cuts are x >= 0 and y >= 0.
    let mut control = source("convex");
    control["steps"][0]["loops"] = json!([[{
        "kind":"circle", "id":"wall", "centre":["0","0"], "radius":"8"
    }]]);
    let mut bevel = control["steps"].as_array_mut().unwrap().pop().unwrap();
    let steps = control["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"clip", "op":"block", "size":["9","9","14"],
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","-1"]}}));
    steps.push(json!({"id":"quarter", "op":"boolean", "kind":"intersect",
        "target":"base", "tool":"clip"}));
    control["result"] = json!("quarter");
    let prefix = Recipe::parse(&control.to_string()).unwrap();
    bevel["body"] = json!("quarter");
    control["steps"].as_array_mut().unwrap().push(bevel);
    control["result"] = json!("bevel");
    let control = Recipe::parse(&control.to_string()).unwrap();

    for (name, recipe, expected_groups) in [
        ("full-circle quadrant control", prefix, 5),
        ("full-circle chamfer control", control, 6),
        ("original arc chamfer", recipe("convex", "2", "1", false), 6),
    ] {
        let options = Options::default();
        assert_eq!(options.tolerance, Rational::parse("0.01").unwrap());
        assert_eq!(options.max_cells, 2_000_000);
        assert_eq!(options.max_vertices, 1_000_000);
        let mut diagnostics = Diagnostics::default();
        let mut sources = Vec::new();
        let (answer, work) = with_operation_budget(1_000_000, || {
            let mut dec = Decider::default();
            let plan = Plan::compile_checked(&recipe, &mut dec)?;
            sources = plan
                .faces()
                .iter()
                .map(|face| face.origin.clone())
                .collect();
            let mesh = mesh_plan_reported(&plan, &options, &mut dec, &mut diagnostics)?;
            Ok::<_, catallax_numbers::Refusal>((plan, mesh))
        });
        if work.exhausted {
            for support in &diagnostics.supports {
                eprintln!(
                    "{}",
                    json!({"case":name,"origins":support.faces.iter().map(|face|&sources[face.index()]).collect::<Vec<_>>(),"cells":support.total.cells,"work":support.total.work,"queries":support.queries.work,"classification":support.classification.work})
                );
            }
        }
        assert!(!work.exhausted, "{name}: {work:?}");
        let (plan, mesh) = answer.unwrap();
        assert!(!mesh.positions.is_empty() && !mesh.triangles.is_empty());
        assert!(mesh.max_deviation <= options.tolerance);
        assert!(mesh.cells <= options.max_cells);
        assert!(mesh.positions.len() <= options.max_vertices);
        let origins: BTreeSet<_> = mesh
            .triangles
            .iter()
            .map(|t| plan.faces()[t.face.index()].origin.clone())
            .collect();
        assert_eq!(origins.len(), expected_groups, "{name}: {origins:?}");
        let mut obj = Vec::new();
        mesh.write_obj(&mut obj).unwrap();
        let obj = String::from_utf8(obj).unwrap();
        let groups: BTreeSet<_> = obj.lines().filter(|line| line.starts_with("g ")).collect();
        assert_eq!(groups.len(), expected_groups, "{name}");
        // A cap vertex may be within tolerance of the adjacent sloping cone.
        // Its radial error is then at most sqrt(1 + 2^2) * tolerance, bounded
        // above by 3*tolerance. Keep every cap triangle in this conservative
        // convex enclosure: vertex inclusion proves whole-triangle exclusion
        // from the removed annulus and the other three angular sectors.
        let margin = options.tolerance.clone();
        let cap_radius = Rational::from_i64(if expected_groups == 5 { 8 } else { 6 })
            + Rational::from_i64(3) * margin.clone();
        let mut cap_triangles = 0;
        for triangle in &mesh.triangles {
            let origin = &plan.faces()[triangle.face.index()].origin;
            if origin.step.as_str() != "base" || origin.role != Role::End {
                continue;
            }
            cap_triangles += 1;
            for &vertex in &triangle.vertices {
                let [x, y, z] =
                    mesh.positions[vertex].map(|v| Rational::from_f64_exact(v).unwrap());

                assert!(x >= -&margin && y >= -&margin, "{name}: excluded sector");
                assert!(
                    &x * &x + &y * &y <= &cap_radius * &cap_radius,
                    "{name}: removed cap annulus"
                );
                assert!((z - Rational::from_i64(12)).abs() <= margin);
            }
        }
        assert!(cap_triangles > 0, "{name}: missing end-cap geometry");
        eprintln!(
            "{name}: {} cells, {} OBJ groups, {cap_triangles} excluded-region-checked cap triangles, {} work",
            mesh.cells,
            groups.len(),
            work.used
        );
    }
}

#[test]
fn default_native_display_keeps_all_trimmed_faces_at_verifier_allowance() {
    use catallax_field::display::{Diagnostics, Options, mesh_plan_reported};
    use catallax_numbers::with_operation_budget;
    use std::collections::BTreeSet;
    for (name, count) in [
        ("convex", 6),
        ("chord", 6),
        ("seam", 5),
        ("major", 6),
        ("major-seam", 6),
    ] {
        for (distance, other) in [("2", "1"), ("1", "1/2")] {
            let options = Options::default();
            assert_eq!(options.tolerance, Rational::parse("0.01").unwrap());
            assert_eq!(options.max_cells, 2_000_000);
            assert_eq!(options.max_vertices, 1_000_000);
            let mut diagnostics = Diagnostics::default();
            let mut origins = Vec::new();
            let (result, work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::default();
                let plan = Plan::compile_checked(&recipe(name, distance, other, false), &mut dec)?;
                origins = plan
                    .faces()
                    .iter()
                    .map(|face| face.origin.clone())
                    .collect();
                let mesh = mesh_plan_reported(&plan, &options, &mut dec, &mut diagnostics)?;
                Ok::<_, catallax_numbers::Refusal>((plan, mesh))
            });
            if work.exhausted {
                for support in &diagnostics.supports {
                    eprintln!(
                        "{}",
                        json!({"case":name,"setbacks":[distance,other],
                        "origins":support.faces.iter().map(|face|&origins[face.index()]).collect::<Vec<_>>(),
                        "cells":support.total.cells,"depth":support.total.max_depth,
                        "work":support.total.work,"classification":support.classification.work,
                        "queries":support.queries.work})
                    );
                }
            }
            assert!(!work.exhausted, "{name} {distance}: {work:?}");
            let (plan, mesh) = result.unwrap();
            assert!(!mesh.positions.is_empty() && !mesh.triangles.is_empty());
            assert!(mesh.max_deviation <= options.tolerance);
            assert!(
                mesh.cells <= options.max_cells && mesh.positions.len() <= options.max_vertices
            );
            let origins: BTreeSet<_> = mesh
                .triangles
                .iter()
                .map(|t| plan.faces()[t.face.index()].origin.clone())
                .collect();
            assert_eq!(origins.len(), count, "{name} {distance}: {origins:?}");
            assert_eq!(
                origins
                    .iter()
                    .filter(|o| matches!(o.role, Role::Bevel { .. }))
                    .count(),
                1
            );
            eprintln!(
                "{name} {distance}/{other}: {} cells, {} vertices, {} triangles, {} work",
                mesh.cells,
                mesh.positions.len(),
                mesh.triangles.len(),
                work.used
            );
        }
    }
}

#[test]
fn chord_end_overflow_remains_a_typed_refusal() {
    let error =
        Plan::compile_checked(&recipe("chord", "4", "1", false), &mut decider()).unwrap_err();
    assert_eq!(error.code, "chamfer.unsupported-overflow");
}

#[test]
fn reversed_arc_reference_preserves_material_setback_order() {
    let mut value = source("convex");
    let bevel = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
    bevel["edges"][0]["between"]
        .as_array_mut()
        .unwrap()
        .reverse();
    bevel["distance"] = json!("1");
    bevel["other_distance"] = json!("2");
    let plan =
        Plan::compile_checked(&Recipe::parse(&value.to_string()).unwrap(), &mut decider()).unwrap();
    for (p, expected) in [
        (["21/5", "28/5", "23/2"], Location::Boundary),
        (["39/10", "26/5", "23/2"], Location::Inside),
        (["9/2", "6", "23/2"], Location::Outside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p, false), &mut Workspace::default(), &mut decider())
                .unwrap()
                .location,
            expected
        );
    }
}
