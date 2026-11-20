//! Feature recognition on a bare import: the P1 bracket written as STEP
//! with no manufacturing information, read back as a body, has the features
//! the built bracket has, on the same faces.
//!
//! Counts and sizes alone would pass on a body whose features had moved or
//! swapped faces. So each feature is also compared by what it owns: the
//! surface of every face in it, and for a hole the very faces it opens onto
//! and the ones that close it, which is what makes it through or blind.

use catallax_features::evaluate;
use catallax_import::{import_step_bodies, sha256_hex};
use catallax_measure::recognition::{Recognised, recognise};
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::Surface;
use catallax_step::{StepOptions, write_step};
use catallax_topology::{FaceId, Model};

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

/// A face by the surface it lies on, in numbers that do not depend on the
/// body's ids: a plane by its unit normal and its distance from the origin, a
/// cylinder by the point of its axis nearest the origin, its unit direction
/// and its radius. Rounded to a nanometre; opposite normals and directions
/// read the same.
fn surface(model: &Model, face: FaceId) -> String {
    let r = |value: f64| format!("{:.6}", value + 0.0).replace("-0.000000", "0.000000");
    let unit = |v: [f64; 3]| {
        let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let first = v.iter().copied().find(|c| c.abs() > 1e-12).unwrap_or(1.0);
        let sign = if first < 0.0 { -1.0 } else { 1.0 };
        [
            sign * v[0] / length,
            sign * v[1] / length,
            sign * v[2] / length,
        ]
    };
    match model.geom().surface(model.face(face).surface) {
        Surface::Plane(plane) => {
            let n = unit([
                plane.normal.x.to_f64(),
                plane.normal.y.to_f64(),
                plane.normal.z.to_f64(),
            ]);
            let o = [
                plane.origin.x.to_f64(),
                plane.origin.y.to_f64(),
                plane.origin.z.to_f64(),
            ];
            let at = n[0] * o[0] + n[1] * o[1] + n[2] * o[2];
            format!("plane ({}, {}, {}) at {}", r(n[0]), r(n[1]), r(n[2]), r(at))
        }
        Surface::Cylinder(cylinder) => {
            let d = unit([
                cylinder.axis.x.to_f64(),
                cylinder.axis.y.to_f64(),
                cylinder.axis.z.to_f64(),
            ]);
            let o = [
                cylinder.origin.x.to_f64(),
                cylinder.origin.y.to_f64(),
                cylinder.origin.z.to_f64(),
            ];
            let along = d[0] * o[0] + d[1] * o[1] + d[2] * o[2];
            let p = [
                o[0] - along * d[0],
                o[1] - along * d[1],
                o[2] - along * d[2],
            ];
            format!(
                "cylinder through ({}, {}, {}) along ({}, {}, {}) radius {}",
                r(p[0]),
                r(p[1]),
                r(p[2]),
                r(d[0]),
                r(d[1]),
                r(d[2]),
                r(cylinder.radius_sq.to_f64().sqrt())
            )
        }
        _ => "another surface".to_owned(),
    }
}

/// Every feature with the faces it owns, each face by its surface.
fn ownership(model: &Model, found: &Recognised) -> Vec<String> {
    let faces = |faces: &[FaceId]| {
        let mut named: Vec<String> = faces.iter().map(|face| surface(model, *face)).collect();
        named.sort();
        named.dedup();
        named.join("; ")
    };
    let mut lines = Vec::new();
    for hole in &found.holes {
        lines.push(format!(
            "hole wall [{}] opens onto [{}] closed by [{}] {}",
            faces(&hole.faces),
            faces(&hole.opens_onto),
            faces(&hole.bottoms),
            if hole.through { "through" } else { "blind" }
        ));
    }
    for pattern in &found.hole_patterns {
        let mut walls: Vec<String> = pattern
            .iter()
            .map(|hole| faces(&found.holes[*hole].faces))
            .collect();
        walls.sort();
        lines.push(format!("pattern of walls [{}]", walls.join(" | ")));
    }
    for boss in &found.bosses {
        lines.push(format!(
            "boss wall [{}] stands on [{}] closed by [{}]",
            faces(&boss.faces),
            faces(&[boss.base]),
            faces(&[boss.top])
        ));
    }
    for slot in &found.slots {
        let mut ends = [faces(&slot.ends[0]), faces(&slot.ends[1])];
        ends.sort();
        lines.push(format!(
            "slot flats [{}] ends [{}] {}",
            faces(&slot.flats),
            ends.join(" | "),
            if slot.through { "through" } else { "blind" }
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

    // Ownership: every feature of the imported body holds the faces its
    // counterpart in the built body holds, and opens onto the same faces.
    let owned = ownership(&model, &imported);
    assert_eq!(owned, ownership(&built.model, &native));
    // The bore's wall is the cylinder of radius 4 about (30, 25), it opens onto
    // the plate's underside and the boss top, and nothing closes it.
    let bore = "hole wall [cylinder through (30.000000, 25.000000, 0.000000) along (0.000000, 0.000000, 1.000000) \
                radius 4.000000] opens onto [plane (0.000000, 0.000000, 1.000000) at 0.000000; plane (0.000000, \
                0.000000, 1.000000) at 15.000000] closed by [] through";
    assert!(owned.iter().any(|line| line == bore), "{owned:#?}");
    // The boss stands on the plate's top and is closed by the face the bore
    // opens onto.
    let boss = "boss wall [cylinder through (30.000000, 25.000000, 0.000000) along (0.000000, 0.000000, 1.000000) \
                radius 8.000000] stands on [plane (0.000000, 0.000000, 1.000000) at 5.000000] closed by [plane \
                (0.000000, 0.000000, 1.000000) at 15.000000]";
    assert!(owned.iter().any(|line| line == boss), "{owned:#?}");
    // The slot's flats are the planes at z 15 and z 21, and it goes through.
    let slot = owned
        .iter()
        .find(|line| line.starts_with("slot"))
        .expect("a slot");
    assert!(
        slot.contains(
            "flats [plane (0.000000, 0.000000, 1.000000) at 15.000000; plane (0.000000, 0.000000, \
                       1.000000) at 21.000000]"
        ) && slot.ends_with("through"),
        "{slot}"
    );
}

#[test]
fn a_blind_hole_stays_blind_through_a_bare_import() {
    // A 20 x 10 x 6 block with a hole of radius 2 from the top down to z 2.
    let text = r#"{"format":"catallax-recipe/3","name":"blind","steps":[
        {"id":"block","op":"block","size":["20","10","6"],
         "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"drill","op":"cylinder","radius":"2","height":"5",
         "pose":{"rotation":["1","0","0","0"],"translation":["5","5","2"]}},
        {"id":"drilled","op":"boolean","kind":"subtract","target":"block","tool":"drill"}],
        "result":"drilled"}"#;
    let built = evaluate(&Recipe::parse(text).unwrap(), &mut Decider::default()).unwrap();
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
    let mut model = Model::new();
    let bodies = import_step_bodies(
        &mut model,
        &bytes,
        &sha256_hex(&bytes),
        &StepId::new("imported").unwrap(),
    )
    .unwrap();
    let imported = recognise(&model, bodies[0], &mut Decider::default()).unwrap();
    let owned = ownership(&model, &imported);
    assert_eq!(owned, ownership(&built.model, &native));
    assert_eq!(
        owned,
        [
            "hole wall [cylinder through (5.000000, 5.000000, 0.000000) along (0.000000, 0.000000, 1.000000) radius \
          2.000000] opens onto [plane (0.000000, 0.000000, 1.000000) at 6.000000] closed by [plane (0.000000, \
          0.000000, 1.000000) at 2.000000] blind"
        ]
    );
}
