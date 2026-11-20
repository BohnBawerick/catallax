//! Feature recognition on the P1 bracket: the bore, the boss, the four bolt
//! holes as one pattern and the flange slot, read off the built body.
//!
//! The expected records are those the bracket's drawing names, and the same
//! ones quiddity 0.5.4 reports on the bracket's STEP file.

use catallax_features::evaluate;
use catallax_measure::recognition::{Recognised, recognise};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

fn bracket() -> Recognised {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    recognise(&built.model, built.result, &mut Decider::default()).unwrap()
}

#[test]
fn the_bracket_has_a_bore_a_boss_four_patterned_holes_and_a_slot() {
    let found = bracket();

    let mut holes: Vec<(f64, bool, Option<f64>)> = found
        .holes
        .iter()
        .map(|hole| (hole.diameter, hole.through, hole.depth))
        .collect();
    holes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(
        holes,
        [
            (4.0, true, Some(5.0)),
            (4.0, true, Some(5.0)),
            (4.0, true, Some(5.0)),
            (4.0, true, Some(5.0)),
            (8.0, true, Some(15.0)),
        ],
        "four bolt holes and the bore, all through"
    );
    assert!(
        found
            .holes
            .iter()
            .all(|hole| hole.bottoms.is_empty() && hole.opens_onto.len() == 2)
    );

    assert_eq!(found.hole_patterns.len(), 1, "one pattern");
    let pattern = &found.hole_patterns[0];
    assert_eq!(pattern.len(), 4);
    assert!(
        pattern
            .iter()
            .all(|hole| found.holes[*hole].diameter == 4.0),
        "of the four bolt holes"
    );

    assert_eq!(found.bosses.len(), 1, "one boss");
    let boss = &found.bosses[0];
    assert_eq!((boss.diameter, boss.height), (16.0, Some(10.0)));
    let bore = found
        .holes
        .iter()
        .find(|hole| hole.diameter == 8.0)
        .unwrap();
    assert!(
        bore.opens_onto.contains(&boss.top),
        "the bore opens onto the boss top"
    );

    assert_eq!(found.slots.len(), 1, "one slot");
    let slot = &found.slots[0];
    assert_eq!((slot.length, slot.width, slot.through), (20.0, 6.0, true));
    assert!(slot.ends.iter().all(|end| !end.is_empty()));
    assert_ne!(slot.flats[0], slot.flats[1]);
}

#[test]
fn a_plain_block_has_no_feature() {
    let recipe = Recipe::parse(
        r#"{"format":"catallax-recipe/3","name":"block","steps":[
            {"id":"p","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
             "loops":[[{"kind":"line","id":"a","to":["10","0"]},{"kind":"line","id":"b","to":["10","10"]},
                       {"kind":"line","id":"c","to":["0","10"]},{"kind":"line","id":"d","to":["0","0"]}]]},
            {"id":"block","op":"extrude","profile":"p","distance":"5"}]}"#,
    )
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let found = recognise(&built.model, built.result, &mut Decider::default()).unwrap();
    assert_eq!(found, Recognised::default());
}
