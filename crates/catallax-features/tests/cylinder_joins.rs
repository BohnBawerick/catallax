//! Public bounded cylinder-join offsets and finite selected normal sweeps.
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::Value;
fn source(name: &str) -> Value {
    serde_json::from_str(match name {
        "boss" => include_str!("fixtures/cylinder-joins/boss-plate.json"),
        "shaft" => include_str!("fixtures/cylinder-joins/stepped-shaft.json"),
        "boss-thicken" => include_str!("fixtures/cylinder-joins/boss-plate-thicken.json"),
        "shaft-thicken" => include_str!("fixtures/cylinder-joins/stepped-shaft-thicken.json"),
        _ => unreachable!(),
    })
    .unwrap()
}
fn recipe(name: &str) -> Recipe {
    Recipe::parse(&source(name).to_string()).unwrap()
}
#[test]
fn public_cylinder_join_sources_construct_complete_analytic_skins() {
    for name in ["boss", "shaft", "boss-thicken", "shaft-thicken"] {
        let r = recipe(name);
        let built = evaluate(&r, &mut Decider::default()).unwrap();
        assert!(
            built.model.check(built.result).is_empty(),
            "{name}: {:?}",
            built.model.check(built.result)
        );
        let _plan = Plan::compile(&r).unwrap();
    }
}
fn points(name: &str, points: &[([&str; 3], Location)]) {
    let r = recipe(name);
    let mut built = evaluate(&r, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&r).unwrap();
    for (p, want) in points {
        let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
        let point = RPoint3::new(x, y, z);
        let field = plan
            .proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
            .unwrap()
            .location;
        assert_eq!(field, *want, "field {name} {p:?}");
        let result = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point,
            &mut Decider::default(),
        );
        let brep = match result {
            Outcome::Proven(catallax_boolean::Containment::Inside) => Location::Inside,
            Outcome::Proven(catallax_boolean::Containment::Outside) => Location::Outside,
            Outcome::ProvenCoincident(_) => Location::Boundary,
            other => panic!("{name} {p:?}: {other:?}"),
        };
        assert_eq!(brep, *want, "B-rep {name} {p:?}");
    }
}
#[test]
fn outward_roots_caps_and_torus_sectors_have_exact_material_controls() {
    use Location::{Boundary, Inside, Outside};
    points(
        "boss",
        &[
            (["10", "8", "33/4"], Boundary),
            (["49/4", "8", "5"], Boundary),
            (["243/20", "8", "41/5"], Boundary),
            (["243/20", "8", "819/100"], Inside),
            (["243/20", "8", "821/100"], Outside),
            (["243/20", "8", "16/5"], Inside),
            (["49/4", "8", "13/4"], Boundary),
            (["123/10", "8", "33/10"], Outside),
            (["10", "8", "4"], Inside),
        ],
    );
    points(
        "shaft",
        &[
            (["0", "0", "-1/4"], Boundary),
            (["17/4", "0", "1"], Boundary),
            (["83/20", "0", "16/5"], Boundary),
            (["83/20", "0", "319/100"], Inside),
            (["83/20", "0", "321/100"], Outside),
            (["9/4", "0", "4"], Boundary),
            (["43/20", "0", "16/5"], Inside),
            (["23/10", "0", "33/10"], Outside),
            (["43/20", "0", "41/5"], Boundary),
            (["0", "0", "4"], Inside),
        ],
    );
}
#[test]
fn selected_wall_sweeps_stop_at_the_actual_source_rims() {
    use Location::{Boundary, Inside, Outside};
    for (name, x) in [("boss-thicken", "243/20"), ("shaft-thicken", "43/20")] {
        let x0 = if name == "boss-thicken" { "10" } else { "0" };
        let y = if name == "boss-thicken" { "8" } else { "0" };
        points(
            name,
            &[
                ([x, y, "4"], Inside),
                ([x, y, "3"], Boundary),
                ([x, y, "299/100"], Outside),
                ([x, y, "801/100"], Outside),
                ([x0, y, "4"], Outside),
            ],
        );
    }
}

fn oracle(name: &str) -> catallax_numbers::Ball {
    use catallax_numbers::Ball;
    let b = |n, d| Ball::from_rational(&Rational::from_ratio(n, d).unwrap(), 256);
    let pi = Ball::pi(256);
    match name {
        // Independent horizontal sections: rounded rectangle base plus the
        // cylinder above the displaced attachment plane and its top round.
        "boss" => b(1174, 1) + &pi * &b(1761, 64) + pi.square() * b(1, 16),
        "shaft" => &pi * &b(5585, 64) + pi.square() * b(5, 16),
        _ => &pi * &b(85, 16),
    }
}
#[test]
fn exact_volume_oracles_enclose_independent_meridian_integrals() {
    for name in ["boss", "shaft", "boss-thicken", "shaft-thicken"] {
        let built = evaluate(&recipe(name), &mut Decider::default()).unwrap();
        let answer = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(answer.marks.is_empty(), "{name}: {:?}", answer.marks);
        let expected = oracle(name);
        let volume = answer.value.volume;
        assert!(
            volume.lo <= expected.lo() && volume.hi >= expected.hi(),
            "{name}: {volume:?}, {expected:?}"
        );
    }
}
#[test]
fn field_volume_enclosures_contain_the_same_exact_integrals() {
    use catallax_field::mass;
    let mut failures = Vec::new();
    for name in ["boss", "shaft", "boss-thicken", "shaft-thicken"] {
        let started = std::time::Instant::now();
        let (answer, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            let plan = Plan::compile(&recipe(name))?;
            mass::properties(
                &plan,
                &mass::Accuracy::default(),
                &mass::Options::default(),
                &mut Decider::default(),
            )
        });
        eprintln!(
            "{name} mass: {} work, {:?} wall; {:?}",
            work.used,
            started.elapsed(),
            answer.as_ref().map(|v| v.cells)
        );
        match answer {
            Ok(answer) => {
                let expected = oracle(name);
                let volume = answer.volume;
                assert!(
                    volume.lo <= expected.lo() && volume.hi >= expected.hi(),
                    "{name}: {volume:?}, {expected:?}"
                );
            }
            Err(e) => failures.push((name, e)),
        }
    }
    assert!(failures.is_empty(), "{failures:?}");
}
#[test]
fn noncoaxial_overflow_and_multiface_selections_keep_typed_refusals() {
    let mut cases = Vec::new();
    let mut noncoaxial = source("shaft");
    noncoaxial["steps"][1]["pose"]["translation"][0] = "1".into();
    cases.push(noncoaxial);
    let mut large = source("shaft");
    large["steps"][3]["distance"] = "1".into();
    cases.push(large);
    let mut short = source("boss");
    short["steps"][1]["height"] = "1/2".into();
    cases.push(short);
    let mut close = source("boss");
    close["steps"][1]["pose"]["translation"][0] = "21/10".into();
    cases.push(close);
    let mut sharp = source("shaft-thicken");
    sharp["steps"][3]["faces"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({"face":{"step":"boss","role":"end"}}));
    cases.push(sharp);
    for value in cases {
        let r = Recipe::parse(&value.to_string()).unwrap();
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&r, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected typed refusal")
        };
        assert_eq!(refusal.code, "distance.unsupported-envelope", "{value}");
        assert_eq!(Plan::compile(&r).unwrap_err().code, refusal.code, "{value}");
    }
}
#[test]
fn selected_annular_shoulder_and_perforated_plate_sweep_only_the_trim() {
    use Location::{Inside, Outside};
    for (name, step) in [("boss-thicken", "plate"), ("shaft-thicken", "base")] {
        let mut value = source(name);
        value["steps"][3]["faces"] = serde_json::json!([{"face":{"step":step,"role":"end"}}]);
        let r = Recipe::parse(&value.to_string()).unwrap();
        let mut built = evaluate(&r, &mut Decider::default()).unwrap();
        assert!(
            built.model.check(built.result).is_empty(),
            "{:?}",
            built.model.check(built.result)
        );
        let plan = Plan::compile(&r).unwrap();
        let (x, y, core) = if step == "plate" {
            ("13", "8", "10")
        } else {
            ("3", "0", "0")
        };
        for (p, want) in [
            ([x, y, "31/10"], Inside),
            ([core, y, "31/10"], Outside),
            ([x, y, "33/10"], Outside),
        ] {
            let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
            let p = RPoint3::new(x, y, z);
            assert_eq!(
                plan.proven_checked(&p, &mut Workspace::default(), &mut Decider::default())
                    .unwrap()
                    .location,
                want
            );
            let answer = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &p,
                &mut Decider::default(),
            );
            assert_eq!(
                answer.proven(),
                Some(&if want == Inside {
                    catallax_boolean::Containment::Inside
                } else {
                    catallax_boolean::Containment::Outside
                })
            );
        }
    }
}

#[test]
fn posed_and_scaled_copies_preserve_exact_volume_and_complete_trims() {
    use catallax_shapes::Rotation;
    let quaternion = ["1", "2", "3", "4"];
    let [w, x, y, z] = quaternion.map(|s| Rational::parse(s).unwrap());
    let rotation = Rotation::from_quaternion(w, x, y, z).unwrap();
    let matrix = rotation.rational_matrix().unwrap();
    for name in ["boss", "shaft"] {
        for scale in ["1/100", "3/2"] {
            let scale = Rational::parse(scale).unwrap();
            let mut value = source(name);
            for stock in value["steps"].as_array_mut().unwrap().iter_mut().take(2) {
                for parameter in ["radius", "height"] {
                    if let Some(s) = stock[parameter].as_str() {
                        stock[parameter] =
                            (&Rational::parse(s).unwrap() * &scale).to_string().into();
                    }
                }
                if let Some(size) = stock.get_mut("size").and_then(Value::as_array_mut) {
                    for n in size {
                        *n = (&Rational::parse(n.as_str().unwrap()).unwrap() * &scale)
                            .to_string()
                            .into();
                    }
                }
                let t: Vec<_> = stock["pose"]["translation"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| Rational::parse(v.as_str().unwrap()).unwrap() * scale.clone())
                    .collect();
                let translated: Vec<_> = (0..3)
                    .map(|i| {
                        (0..3)
                            .fold(Rational::from_i64([3, -2, 5][i]), |a, j| {
                                a + matrix[i][j] * &t[j]
                            })
                            .to_string()
                    })
                    .collect();
                stock["pose"]["rotation"] = serde_json::json!(quaternion);
                stock["pose"]["translation"] = serde_json::json!(translated);
            }
            value["steps"][3]["distance"] = (&Rational::parse("1/4").unwrap() * &scale)
                .to_string()
                .into();
            let r = Recipe::parse(&value.to_string()).unwrap();
            let built = evaluate(&r, &mut Decider::default()).unwrap();
            assert!(built.model.check(built.result).is_empty());
            let _plan = Plan::compile(&r).unwrap();
            let answer = catallax_measure::mass_properties(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert!(
                answer.marks.is_empty(),
                "{name} {scale}: {:?}",
                answer.marks
            );
            let expected = oracle(name)
                * catallax_numbers::Ball::from_rational(&(&(&scale * &scale) * &scale), 256);
            assert!(
                answer.value.volume.lo <= expected.lo() && answer.value.volume.hi >= expected.hi()
            );
        }
    }
}
#[test]
fn brep_display_sews_analytic_contacts_and_spherical_plate_corners() {
    use std::collections::BTreeMap;
    for name in ["boss", "shaft", "boss-thicken", "shaft-thicken"] {
        let r = recipe(name);
        let built = evaluate(&r, &mut Decider::default()).unwrap();
        let options = catallax_mesh::MeshOptions {
            chord_tolerance: Rational::parse("1/200").unwrap(),
            ..Default::default()
        };
        let mesh = catallax_mesh::mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(mesh.marks.is_empty(), "{name}: {:?}", mesh.marks);
        let mesh = mesh.value;
        assert!(!mesh.triangles.is_empty());
        let mut edges = BTreeMap::new();
        for t in mesh.triangles {
            let [a, b, c] = t.vertices;
            for (a, b) in [(a, b), (b, c), (c, a)] {
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        assert!(
            edges.values().all(|uses| *uses == 2),
            "{name}: unsewn display edge"
        );
    }
}
