//! A concave round on the rings where a cylinder passes through a sphere.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate, evaluate_certified};
use catallax_measure::Accuracy;
use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::Value;

fn input() -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "sphere-cylinder")
        .unwrap()["recipe"]
        .clone()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 200_000_000,
        ..Default::default()
    })
}

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn on_thread(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn both_rings_round_with_certified_bands() {
    on_thread(|| {
        let recipe = Recipe::parse(&input().to_string()).unwrap();
        let certified = evaluate_certified(&recipe, &mut decider()).unwrap();
        let junction = certified.junction().unwrap();
        assert_eq!(junction.sphere_bands().len(), 2);
        assert!(junction.saddles().is_empty());
        for band in junction.sphere_bands() {
            assert!(band.bound().upper <= q("1/10000"));
            assert_ne!(band.rails()[0], band.rails()[1]);
        }
        let built = certified.built();
        let mut model = built.model.clone();
        assert!(model.check(built.result).is_empty());
        assert_eq!(model.faces_of_body(built.result).len(), 7);
        // The ring passes (3, 0, 4); the ball there has its centre at
        // (7/2, 0, sqrt(18)), so a point just off the sharp edge is new
        // material and the centre is outside.
        for (coordinates, expected) in [
            (["1", "0", "0"], Containment::Inside),
            (["1", "0", "9"], Containment::Inside),
            (["151/50", "0", "401/100"], Containment::Inside),
            (["151/50", "0", "-401/100"], Containment::Inside),
            (["7/2", "0", "17/4"], Containment::Outside),
            (["-26/5", "0", "0"], Containment::Outside),
        ] {
            let [x, y, z] = coordinates.map(q);
            let actual = classify_point(
                &mut model,
                built.result,
                &RPoint3::new(x, y, z),
                &mut decider(),
            );
            assert!(
                matches!(&actual, Outcome::Proven(found) if *found == expected),
                "at {coordinates:?}: {actual:?}"
            );
        }
        // No consumer measures or meshes the band yet, and each says so.
        let accuracy = Accuracy {
            relative: q("1/1000000"),
            absolute: q("1/1000000"),
        };
        let refusal =
            catallax_measure::volume(&model, built.result, &accuracy, &mut decider()).unwrap_err();
        assert_eq!(refusal.code, "measure.unsupported-surface", "{refusal:?}");
        let refusal =
            catallax_measure::bounding_box(&model, built.result, &accuracy, &mut decider())
                .unwrap_err();
        assert!(refusal.code.starts_with("measure."), "{refusal:?}");
        assert!(
            catallax_mesh::mesh_body(
                &model,
                built.result,
                &MeshOptions::default(),
                &mut decider()
            )
            .is_err()
        );
    });
}

#[test]
fn a_cylinder_that_does_not_pass_through_keeps_its_refusal() {
    on_thread(|| {
        let mut value = input();
        // At a gap of 16/5 the radius-2 cylinder leaves the radius-5 sphere
        // through its side, so the edge is one loop round neither.
        value["steps"][1]["pose"]["translation"] = serde_json::json!(["16/5", "0", "-10"]);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let error = evaluate(&recipe, &mut decider()).unwrap_err();
        let BuildError::Geometry { refusal, .. } = error else {
            panic!("expected a construction refusal, got {error:?}");
        };
        assert_eq!(refusal.code, "fillet.unsupported-contact", "{refusal:?}");
    });
}
