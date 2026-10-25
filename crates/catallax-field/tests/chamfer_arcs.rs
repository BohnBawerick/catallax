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
            ],
        ),
        (
            "concave",
            vec![
                (["27/5", "36/5", "9/2"], Boundary),
                (["51/10", "34/5", "9/2"], Inside),
                (["57/10", "38/5", "9/2"], Outside),
                (["-27/5", "36/5", "9/2"], Outside),
                (["10", "0", "4"], Boundary),
                (["8", "0", "5"], Boundary),
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
    for name in ["convex", "concave", "chord", "seam", "major", "major-seam"] {
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
                let z = if name == "concave" {
                    if other == "2" {
                        "5"
                    } else if other == "1/2" {
                        "17/4"
                    } else {
                        "9/2"
                    }
                } else if other == "2" {
                    "11"
                } else if other == "1/2" {
                    "47/4"
                } else {
                    "23/2"
                };
                let radius = if name == "concave" {
                    if distance == "1" { "17/2" } else { "9" }
                } else if distance == "1" {
                    "15/2"
                } else {
                    "7"
                };
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
fn default_native_display_keeps_all_trimmed_faces_at_verifier_allowance() {
    use catallax_field::display::{Options, mesh_plan};
    use catallax_numbers::with_operation_budget;
    use std::collections::BTreeSet;
    for (name, count) in [
        ("convex", 6),
        ("concave", 9),
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
            let (result, work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::default();
                let plan = Plan::compile_checked(&recipe(name, distance, other, false), &mut dec)?;
                let mesh = mesh_plan(&plan, &options, &mut dec)?;
                Ok::<_, catallax_numbers::Refusal>((plan, mesh))
            });
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
