//! A plate-edge chamfer beside a boss whose rim is filleted.
//!
//! The fixture is the question SolidDraft sent for a 60 x 40 x 8 plate with a
//! 12 mm x 10 mm boss at (-15, 0), an R2 fillet on the boss's top rim and a
//! 1.5 mm chamfer on the plate's bottom edge at x = -30. The chamfer was refused
//! with `chamfer.unsupported-overflow` while the fillet existed and built
//! without it: the fillet's torus support, expanded about the world origin,
//! could not be separated from the chamfer's affected region box.

use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, Surface};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("fixtures/chamfer-after-fillet.json");

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

/// The fixture with the boss moved, the fillet optionally removed, and the
/// chamfer optionally left out of the result.
fn recipe(centre: [&str; 2], fillet: bool, chamfer: bool) -> Recipe {
    let mut value: Value = serde_json::from_str(FIXTURE).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    steps[2]["loops"][0][0]["centre"] = json!(centre);
    if !fillet {
        steps.retain(|s| s["op"] != "fillet");
        steps.last_mut().unwrap()["body"] = json!("op_40");
    }
    if !chamfer {
        steps.pop();
    }
    let result = steps.last().unwrap()["id"].clone();
    value["result"] = result;
    Recipe::parse(&value.to_string()).unwrap()
}

fn volume(recipe: &Recipe) -> (Rational, Rational) {
    let (built, work) =
        with_operation_budget(1_000_000, || evaluate(recipe, &mut Decider::default()));
    assert!(!work.exhausted, "{work:?}");
    let built = built.unwrap_or_else(|e| panic!("{e}: {work:?}"));
    assert!(built.model.check(built.result).is_empty());
    assert!(built.model.body(built.result).marks.is_empty());
    let v = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    (v.lo, v.hi)
}

#[test]
fn soliddraft_plate_chamfer_builds_beside_a_filleted_boss() {
    let recipe = Recipe::parse(FIXTURE).unwrap();
    let (built, work) =
        with_operation_budget(1_000_000, || evaluate(&recipe, &mut Decider::default()));
    assert!(!work.exhausted, "{work:?}");
    let mut built = built.unwrap_or_else(|e| panic!("{e}: {work:?}"));
    assert!(built.model.check(built.result).is_empty());
    assert!(built.model.body(built.result).marks.is_empty());
    assert_eq!(built.chamfers.len(), 1);
    assert!(matches!(built.chamfers[0].surface, Surface::Plane(_)));
    assert!(
        built
            .fillets
            .iter()
            .any(|f| matches!(f.surface, Surface::Torus(_)))
    );
    // The bevel plane is x + z = -57/2 between x = -30 and z = 0.
    let field = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, inside) in [
        (["-59/2", "0", "1/2"], false),
        (["-29", "0", "1"], true),
        (["-59/2", "19", "1/2"], false),
        (["-15", "0", "17"], true),
    ] {
        let point = RPoint3::new(q(p[0]), q(p[1]), q(p[2]));
        let (brep, location) = if inside {
            (Containment::Inside, Location::Inside)
        } else {
            (Containment::Outside, Location::Outside)
        };
        assert_eq!(
            classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default()
            )
            .proven(),
            Some(&brep),
            "{p:?}"
        );
        assert_eq!(
            field.proven(&point, &mut workspace).location,
            location,
            "{p:?}"
        );
    }
}

#[test]
fn the_chamfer_removes_the_same_prism_with_or_without_the_fillet() {
    // A boss near the chamfered edge separates only after recentring, with or
    // without the fillet's torus.
    for centre in [["-15", "0"], ["-22", "-12"]] {
        for fillet in [false, true] {
            let (lo, hi) = volume(&recipe(centre, fillet, true));
            let (before_lo, before_hi) = volume(&recipe(centre, fillet, false));
            // 1.5 * 1.5 / 2 across the 40 mm edge.
            let removed = q("45");
            assert!(&before_lo - &hi <= removed && removed <= &before_hi - &lo);
            assert!(&(&before_hi - &lo) - &(&before_lo - &hi) < q("1/1000000"));
        }
    }
}
