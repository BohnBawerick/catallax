//! Feature recognition on the P1 bracket: the bore, the boss, the four bolt
//! holes as one pattern and the flange slot, read off the built body.
//!
//! The expected records are those the bracket's drawing names, and the same
//! ones quiddity 0.5.4 reports on the bracket's STEP file.

use catallax_features::evaluate;
use catallax_measure::recognition::{Recognised, recognise};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

/// True when a reported size is the expected one to the last few bits.
fn is(size: f64, expected: f64) -> bool {
    (size - expected).abs() < 1e-12
}

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
            .all(|hole| is(found.holes[*hole].diameter, 4.0)),
        "of the four bolt holes"
    );

    assert_eq!(found.bosses.len(), 1, "one boss");
    let boss = &found.bosses[0];
    assert_eq!((boss.diameter, boss.height), (16.0, Some(10.0)));
    let bore = found
        .holes
        .iter()
        .find(|hole| is(hole.diameter, 8.0))
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

/// A 20 x 10 x 6 block with the given cuts, each a `(tool, json)` step.
fn block_with(cuts: &[&str]) -> Recognised {
    let pose = |x: &str, y: &str, z: &str| {
        format!(r#"{{"rotation":["1","0","0","0"],"translation":["{x}","{y}","{z}"]}}"#)
    };
    let mut steps = vec![format!(
        r#"{{"id":"s0","op":"block","size":["20","10","6"],"pose":{}}}"#,
        pose("0", "0", "0")
    )];
    for (n, cut) in cuts.iter().enumerate() {
        let (kind, tool) = cut.split_once(' ').unwrap();
        steps.push(format!(r#"{{"id":"t{n}",{tool}}}"#));
        steps.push(format!(
            r#"{{"id":"s{}","op":"boolean","kind":"{kind}","target":"s{n}","tool":"t{n}"}}"#,
            n + 1
        ));
    }
    let text = format!(
        r#"{{"format":"catallax-recipe/3","name":"block","steps":[{}],"result":"s{}"}}"#,
        steps.join(","),
        cuts.len()
    );
    let built = evaluate(&Recipe::parse(&text).unwrap(), &mut Decider::default()).unwrap();
    recognise(&built.model, built.result, &mut Decider::default()).unwrap()
}

fn cylinder(radius: &str, height: &str, at: [&str; 3]) -> String {
    format!(
        r#""op":"cylinder","radius":"{radius}","height":"{height}","pose":{{"rotation":["1","0","0","0"],"translation":["{}","{}","{}"]}}"#,
        at[0], at[1], at[2]
    )
}

fn block(size: [&str; 3], at: [&str; 3]) -> String {
    format!(
        r#""op":"block","size":["{}","{}","{}"],"pose":{{"rotation":["1","0","0","0"],"translation":["{}","{}","{}"]}}"#,
        size[0], size[1], size[2], at[0], at[1], at[2]
    )
}

#[test]
fn a_plain_block_has_no_feature() {
    assert_eq!(block_with(&[]), Recognised::default());
}

#[test]
fn a_hole_with_a_floor_is_blind() {
    let found = block_with(&[&format!("subtract {}", cylinder("2", "5", ["5", "5", "2"]))]);
    assert_eq!(found.holes.len(), 1);
    let hole = &found.holes[0];
    assert_eq!(
        (hole.diameter, hole.through, hole.depth),
        (4.0, false, Some(4.0))
    );
    assert_eq!((hole.opens_onto.len(), hole.bottoms.len()), (1, 1));
    assert!(found.bosses.is_empty() && found.slots.is_empty() && found.hole_patterns.is_empty());
}

#[test]
fn holes_of_different_depth_are_not_one_pattern() {
    let found = block_with(&[
        &format!("subtract {}", cylinder("2", "5", ["5", "5", "2"])),
        &format!("subtract {}", cylinder("2", "5", ["15", "5", "3"])),
    ]);
    assert_eq!(found.holes.len(), 2);
    assert!(found.hole_patterns.is_empty());
}

#[test]
fn a_cylinder_on_a_face_is_a_boss_and_a_bare_cylinder_is_not() {
    let found = block_with(&[&format!("union {}", cylinder("3", "4", ["10", "5", "6"]))]);
    assert_eq!(found.bosses.len(), 1);
    assert_eq!(
        (found.bosses[0].diameter, found.bosses[0].height),
        (6.0, Some(4.0))
    );
    assert!(found.holes.is_empty() && found.slots.is_empty());

    let text = format!(
        r#"{{"format":"catallax-recipe/3","name":"rod","steps":[{{"id":"rod",{}}}],"result":"rod"}}"#,
        cylinder("3", "4", ["0", "0", "0"])
    );
    let built = evaluate(&Recipe::parse(&text).unwrap(), &mut Decider::default()).unwrap();
    let found = recognise(&built.model, built.result, &mut Decider::default()).unwrap();
    assert_eq!(
        found,
        Recognised::default(),
        "a rod is the part, not a boss on it"
    );
}

#[test]
fn a_slot_with_a_floor_is_blind_and_a_wider_pocket_is_not_a_slot() {
    // Two round ends 8 apart and the stretch between them, cut 3 deep from the top.
    let blind = [
        format!("subtract {}", cylinder("2", "4", ["6", "5", "3"])),
        format!("subtract {}", cylinder("2", "4", ["14", "5", "3"])),
        format!("subtract {}", block(["8", "4", "4"], ["6", "3", "3"])),
    ];
    let found = block_with(&blind.iter().map(String::as_str).collect::<Vec<_>>());
    assert_eq!(found.slots.len(), 1, "{found:?}");
    let slot = &found.slots[0];
    assert_eq!((slot.length, slot.width, slot.through), (12.0, 4.0, false));
    assert!(found.holes.is_empty(), "a slot's ends are not holes");

    // The same ends with a stretch wider than their diameter: the flats are
    // not tangent, so it is no slot.
    let wide = [
        format!("subtract {}", cylinder("2", "4", ["6", "5", "3"])),
        format!("subtract {}", cylinder("2", "4", ["14", "5", "3"])),
        format!("subtract {}", block(["8", "6", "4"], ["6", "2", "3"])),
    ];
    let found = block_with(&wide.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(found.slots.is_empty(), "{found:?}");
}
