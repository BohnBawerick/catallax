//! Public trimmed-arc controls: analytic material, topology, names and consumers.
use catallax_features::{Built, evaluate};
use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Curve, RPoint3, Surface};
use catallax_topology::{Loop, Sense};
use serde_json::{Value, json};
use std::collections::BTreeMap;

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn source(name: &str) -> Value {
    serde_json::from_str(match name {
        "convex" => include_str!("fixtures/chamfer-arcs/convex.json"),
        "concave" => include_str!("fixtures/chamfer-arcs/concave.json"),
        "chord" => include_str!("fixtures/chamfer-arcs/chord.json"),
        "seam" => include_str!("fixtures/chamfer-arcs/seam.json"),
        "major" => include_str!("fixtures/chamfer-arcs/major.json"),
        "major-seam" => include_str!("fixtures/chamfer-arcs/major-seam.json"),
        _ => unreachable!(),
    })
    .unwrap()
}
fn recipe(name: &str, plane: &str, cylinder: &str, mirror: bool) -> Recipe {
    let mut value = source(name);
    let bevel = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
    bevel["distance"] = json!(plane);
    bevel["other_distance"] = json!(cylinder);
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
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    assert!(built.model.body(built.result).marks.is_empty());
    built
}
fn point(p: [&str; 3], mirror: bool) -> RPoint3 {
    let [mut x, y, z] = p.map(|v| Rational::parse(v).unwrap());
    if mirror {
        x = -x;
    }
    RPoint3::new(x, y, z)
}

#[test]
fn finite_concave_arcs_remain_typed_refusals_for_edits_and_reflections() {
    for mirror in [false, true] {
        for (plane, cylinder) in [("2", "1"), ("2", "2"), ("1", "1/2")] {
            let error =
                evaluate(&recipe("concave", plane, cylinder, mirror), &mut decider()).unwrap_err();
            assert!(
                matches!(error,catallax_features::BuildError::Geometry{refusal,..} if refusal.code=="chamfer.unsupported-concave-arc")
            );
        }
    }
}

#[test]
fn finite_cones_share_exact_contact_arcs_and_transverse_end_boundaries() {
    for name in ["convex", "chord", "seam", "major", "major-seam"] {
        for cylinder in ["1", "2"] {
            for mirror in [false, true] {
                let built = build(&recipe(name, "2", cylinder, mirror));
                let faces = built.model.faces_of_body(built.result);
                let bevels: Vec<_> = faces.iter().copied().filter(|f| built.model.face_origins(*f).iter().any(|o| {
                    if mirror { matches!(&o.role, Role::Copy { of, .. } if matches!(of.role, Role::Bevel { .. })) }
                    else { matches!(o.role, Role::Bevel { .. }) }
                })).collect();
                assert_eq!(bevels.len(), 1, "{name} {cylinder} {mirror}");
                let face = built.model.face(bevels[0]);
                let surface = built.model.geom().surface(face.surface);
                assert!(matches!(surface, Surface::Cone(_)));
                let polynomial = surface.implicit().unwrap();
                assert_eq!(face.loops.len(), 1);
                let Loop::Coedges(cs) = built.model.loop_(face.loops[0]) else {
                    panic!("finite cone loop")
                };
                assert_eq!(cs.len(), 4);
                let mut circles = 0;
                for c in cs {
                    let edge = built.model.edge(built.model.coedge(*c).edge);
                    circles += usize::from(matches!(
                        built.model.geom().curve(edge.curve),
                        Curve::Circle(_)
                    ));
                    for v in [edge.start.unwrap(), edge.end.unwrap()] {
                        let p = built.model.geom().position(built.model.vertex(v).point);
                        assert!(matches!(
                            decider().sign(&polynomial.eval_real(&[p.x, p.y, p.z])),
                            Outcome::ProvenCoincident(_)
                        ));
                    }
                }
                assert_eq!(circles, 2);
                let mut uses: BTreeMap<_, Vec<_>> = BTreeMap::new();
                for f in &faces {
                    for lp in &built.model.face(*f).loops {
                        if let Loop::Coedges(cs) = built.model.loop_(*lp) {
                            for c in cs {
                                let c = built.model.coedge(*c);
                                uses.entry(c.edge).or_default().push(c.sense);
                            }
                        }
                    }
                }
                for senses in uses.values() {
                    assert_eq!(senses.len(), 2);
                    assert!(senses.contains(&Sense::Same) && senses.contains(&Sense::Reversed));
                }
            }
        }
    }
}

#[test]
fn quarter_strip_volumes_match_independent_revolution_integrals() {
    for (name, base, direction) in [("convex", 192, -1)] {
        for (plane, cylinder) in [("2", "1"), ("2", "2"), ("1", "1/2")] {
            let built = build(&recipe(name, plane, cylinder, false));
            let p = Rational::parse(plane).unwrap();
            let c = Rational::parse(cylinder).unwrap();
            // Integral of pi/4*((8 +/- p*t)^2 - 8^2), t=0..1,
            // with sign selected independently from the authored material.
            let delta = c
                * (Rational::from_i64(8) * p.clone()
                    + Rational::from_i64(direction) * p.clone() * p / Rational::from_i64(3))
                / Rational::from_i64(4);
            let factor = Rational::from_i64(base) + Rational::from_i64(direction) * delta;
            let expected = Ball::pi(256) * Ball::from_rational(&factor, 256);
            let measured = catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut decider(),
            )
            .unwrap();
            assert!(measured.marks.is_empty());
            assert!(
                measured.value.lo <= expected.hi() && measured.value.hi >= expected.lo(),
                "{name}: {:?}",
                measured.value
            );
            assert!(measured.value.width() <= Rational::parse("1e-9").unwrap());
            eprintln!(
                "{}",
                json!({"case":name,"plane_setback":plane,"cylinder_setback":cylinder,
                "volume":{"lo":measured.value.lo.to_string(),"hi":measured.value.hi.to_string()},
                "analytic_pi_factor":factor.to_string(),"max_work":1_000_000})
            );
        }
    }
}

#[test]
fn material_controls_include_selected_arc_opposite_arc_and_mirror() {
    use catallax_boolean::Containment::{Inside, Outside};
    for (name, cases) in [
        (
            "convex",
            vec![
                (["21/5", "28/5", "23/2"], None),
                (["39/10", "26/5", "23/2"], Some(Inside)),
                (["9/2", "6", "23/2"], Some(Outside)),
                (["-21/5", "28/5", "23/2"], Some(Outside)),
                (["3", "4", "1"], Some(Inside)),
            ],
        ),
        (
            "chord",
            vec![
                (["5", "0", "12"], None),
                (["24/5", "18/5", "12"], None),
                (["24/5", "32/5", "11"], None),
                (["4", "3", "11"], Some(Outside)),
            ],
        ),
        (
            "seam",
            vec![
                (["28/5", "-21/5", "23/2"], None),
                (["6", "-4", "23/2"], Some(Outside)),
                (["26/5", "-39/10", "23/2"], Some(Inside)),
                (["-28/5", "21/5", "23/2"], Some(Outside)),
            ],
        ),
        (
            "major",
            vec![
                (["-21/5", "-28/5", "23/2"], None),
                (["21/5", "28/5", "23/2"], None),
                (["21/5", "-28/5", "23/2"], Some(Outside)),
            ],
        ),
    ] {
        for mirror in [false, true] {
            let mut built = build(&recipe(name, "2", "1", mirror));
            for (p, expected) in &cases {
                let result = catallax_boolean::classify_point(
                    &mut built.model,
                    built.result,
                    &point(*p, mirror),
                    &mut decider(),
                );
                match expected {
                    Some(expected) => assert!(
                        matches!(result,Outcome::Proven(actual) if actual == *expected),
                        "{name} {mirror} {p:?}: {result:?}"
                    ),
                    None => assert!(
                        matches!(result, Outcome::ProvenCoincident(_)),
                        "{name} {mirror} {p:?}: {result:?}"
                    ),
                }
            }
        }
    }
}

#[test]
fn setback_edit_keeps_resolved_origins_and_default_step_supports() {
    for name in ["convex", "chord", "seam", "major", "major-seam"] {
        for mirror in [false, true] {
            let mut previous = None;
            for (plane, cylinder) in [("2", "1"), ("1", "1/2")] {
                let built = build(&recipe(name, plane, cylinder, mirror));
                let mut origins: Vec<_> = built
                    .model
                    .faces_of_body(built.result)
                    .into_iter()
                    .map(|f| built.model.face_origins(f).to_vec())
                    .collect();
                origins.sort();
                if let Some(old) = previous {
                    assert_eq!(origins, old);
                }
                previous = Some(origins);
                let mut step = Vec::new();
                let options = catallax_step::StepOptions::default();
                let report = catallax_step::write_step(
                    &built.model,
                    &[built.result],
                    &options,
                    &mut decider(),
                    &mut step,
                )
                .unwrap();
                assert!(report.marks.is_empty());
                assert!(report.value.max_deviation <= options.curve_tolerance);
                let step = String::from_utf8(step).unwrap();
                assert!(step.contains("CONICAL_SURFACE"));
                assert!(step.contains("CIRCLE"));
            }
        }
    }
}

#[test]
fn default_brep_display_keeps_the_conical_strip_and_shared_boundaries() {
    use catallax_mesh::{MeshOptions, mesh_body};
    use catallax_numbers::with_operation_budget;
    for name in ["convex", "chord", "seam", "major", "major-seam"] {
        let built = build(&recipe(name, "2", "1", false));
        let options = MeshOptions::default();
        assert_eq!(options.chord_tolerance, Rational::parse("0.01").unwrap());
        assert_eq!(options.max_vertices, 1_000_000);
        let (result, work) = with_operation_budget(100_000_000, || {
            mesh_body(
                &built.model,
                built.result,
                &options,
                &mut Decider::default(),
            )
        });
        assert!(!work.exhausted, "{name}: {work:?}");
        let mesh = result.unwrap().value;
        assert!(!mesh.triangles.is_empty() && mesh.positions.len() <= options.max_vertices);
        assert_eq!(
            mesh.faces.len(),
            built.model.faces_of_body(built.result).len()
        );
        for face in &mesh.faces {
            assert!(!face.triangles.is_empty(), "{name}: {face:?}");
            assert!(face.max_chord_error <= options.chord_tolerance);
        }
        eprintln!(
            "{name}: {} vertices, {} triangles, {} work",
            mesh.positions.len(),
            mesh.triangles.len(),
            work.used
        );
    }
}

#[test]
fn ordered_arc_references_assign_setbacks_to_the_same_supports() {
    let first = recipe("convex", "2", "1", false);
    let mut reversed = source("convex");
    let bevel = reversed["steps"]
        .as_array_mut()
        .unwrap()
        .last_mut()
        .unwrap();
    bevel["edges"][0]["between"]
        .as_array_mut()
        .unwrap()
        .reverse();
    bevel["distance"] = json!("1");
    bevel["other_distance"] = json!("2");
    let second = Recipe::parse(&reversed.to_string()).unwrap();
    let mut files = Vec::new();
    for recipe in [first, second] {
        let built = build(&recipe);
        let mut file = Vec::new();
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut decider(),
            &mut file,
        )
        .unwrap();
        files.push(file);
    }
    assert_eq!(files[0], files[1]);
}

#[test]
fn a_chord_end_consuming_the_contact_circle_refuses_overflow() {
    let error = evaluate(&recipe("chord", "4", "1", false), &mut decider()).unwrap_err();
    assert!(
        matches!(error,catallax_features::BuildError::Geometry { refusal,.. } if refusal.code == "chamfer.unsupported-overflow")
    );
}

#[test]
fn retained_generator_point_query_gap_remains_explicit() {
    let mut built = build(&recipe("seam", "2", "1", false));
    let result = catallax_boolean::classify_point(
        &mut built.model,
        built.result,
        &point(["7", "0", "23/2"], false),
        &mut decider(),
    );
    assert!(
        matches!(result,Outcome::Refused(r) if r.code == "intersect.partial-containment-required")
    );
}
