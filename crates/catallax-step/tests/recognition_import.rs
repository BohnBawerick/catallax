//! Feature recognition on a bare import: the P1 bracket written as STEP
//! with no manufacturing information, read back as a body, has the features
//! the built bracket has.

use catallax_features::evaluate;
use catallax_import::{import_step_bodies, sha256_hex};
use catallax_measure::recognition::{Recognised, recognise};
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, StepId};
use catallax_step::{StepOptions, write_step};
use catallax_topology::Model;

/// What recognition found, without the face ids, which differ between a
/// built and an imported body.
fn summary(found: &Recognised) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for hole in &found.holes {
        lines.push(format!(
            "hole {} through={} depth={:?} mouths={} bottoms={}",
            hole.diameter,
            hole.through,
            hole.depth,
            hole.opens_onto.len(),
            hole.bottoms.len()
        ));
    }
    for pattern in &found.hole_patterns {
        lines.push(format!(
            "pattern of {} holes {}",
            pattern.len(),
            found.holes[pattern[0]].diameter
        ));
    }
    for boss in &found.bosses {
        lines.push(format!("boss {} height={:?}", boss.diameter, boss.height));
    }
    for slot in &found.slots {
        lines.push(format!(
            "slot {} x {} through={}",
            slot.length, slot.width, slot.through
        ));
    }
    lines.sort();
    lines
}

#[test]
fn the_imported_bracket_has_the_features_of_the_built_one() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let native = recognise(&built.model, built.result, &mut Decider::default()).unwrap();

    let mut bytes = Vec::new();
    write_step(
        &built.model,
        &[built.result],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let text = String::from_utf8(bytes.clone()).unwrap();
    assert!(
        !text.contains("DATUM") && !text.contains("DIMENSIONAL_"),
        "the file carries no manufacturing information"
    );

    let mut model = Model::new();
    let bodies = import_step_bodies(
        &mut model,
        &bytes,
        &sha256_hex(&bytes),
        &StepId::new("imported").unwrap(),
    )
    .unwrap();
    assert_eq!(bodies.len(), 1);
    let imported = recognise(&model, bodies[0], &mut Decider::default()).unwrap();

    assert_eq!(
        summary(&imported),
        [
            "boss 16 height=Some(10.0)",
            "hole 4 through=true depth=Some(5.0) mouths=2 bottoms=0",
            "hole 4 through=true depth=Some(5.0) mouths=2 bottoms=0",
            "hole 4 through=true depth=Some(5.0) mouths=2 bottoms=0",
            "hole 4 through=true depth=Some(5.0) mouths=2 bottoms=0",
            "hole 8 through=true depth=Some(15.0) mouths=2 bottoms=0",
            "pattern of 4 holes 4",
            "slot 20 x 6 through=true",
        ]
    );
    assert_eq!(summary(&imported), summary(&native));
}
