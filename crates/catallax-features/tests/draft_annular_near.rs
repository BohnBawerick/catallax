//! Nearest-face selection among the four trimmed faces of a hollow annulus.
use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_field::{Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::{Value as Json, json};

fn document() -> Json {
    serde_json::from_str(include_str!("fixtures/draft-annular/both.json")).unwrap()
}
fn recipe(source: &Json) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}
fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(|v| Rational::parse(v).unwrap());
    RPoint3::new(x, y, z)
}
fn with_faces(mut source: Json, faces: Json) -> Json {
    source["steps"][2]["faces"] = faces;
    source
}
fn near(xyz: [&str; 3]) -> Json {
    json!([{"face":{"step":"base"},"near":xyz}])
}
fn wall(segment: &str) -> Json {
    json!([{"face":{"step":"base","role":{"side":{"segment":segment}}}}])
}
fn brep_code(source: &Json) -> Option<String> {
    match evaluate(&recipe(source), &mut Decider::default()) {
        Ok(_) => None,
        Err(BuildError::Geometry { refusal, .. }) => Some(refusal.code.to_string()),
        Err(other) => panic!("{other:?}"),
    }
}
fn field_code(source: &Json) -> Option<String> {
    Plan::compile(&recipe(source))
        .err()
        .map(|e| e.code.to_string())
}

/// Both kernels refuse with `code`.
fn refuses(source: &Json, code: &str) {
    assert_eq!(brep_code(source).as_deref(), Some(code), "B-rep");
    assert_eq!(field_code(source).as_deref(), Some(code), "field");
}

/// Both kernels build the same solid as the explicit-role control.
fn selects_like(source: &Json, control: &Json) {
    let mut built = evaluate(&recipe(source), &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let mut expected = evaluate(&recipe(control), &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe(source)).unwrap();
    let reference = Plan::compile(&recipe(control)).unwrap();
    let mut storage = Workspace::default();
    let mut other = Workspace::default();
    let samples = [
        "-9", "-6", "-3", "-2", "0", "1", "2", "3", "4", "5", "7", "9",
    ];
    let heights = ["-1", "0", "1/2", "1", "2", "3"];
    let mut classes = std::collections::BTreeSet::new();
    for (i, x) in samples.iter().enumerate() {
        for (j, z) in heights.iter().enumerate() {
            let p = point([x, if (i + j) % 2 == 0 { "0" } else { "1/3" }, z]);
            let a = plan.proven(&p, &mut storage).location;
            assert_eq!(
                a,
                reference.proven(&p, &mut other).location,
                "field {x} {z}"
            );
            classes.insert(format!("{a:?}"));
            let b = classify_point(&mut built.model, built.result, &p, &mut Decider::default());
            let c = classify_point(
                &mut expected.model,
                expected.result,
                &p,
                &mut Decider::default(),
            );
            let tag = |o: &Outcome<Containment>| match o {
                Outcome::Proven(Containment::Inside) => 0,
                Outcome::Proven(Containment::Outside) => 1,
                Outcome::ProvenCoincident(_) => 2,
                _ => 3,
            };
            assert_eq!(tag(&b), tag(&c), "B-rep {x} {z}");
            assert!(tag(&b) < 3);
        }
    }
    assert!(classes.len() >= 2);
}

#[test]
fn hole_axis_and_outer_rim_points_select_the_inner_and_outer_walls() {
    // Hole axis: the inner wall is 2 away, the caps are no nearer than 5 and a
    // cap disc that filled the hole would have won at distance 1.
    selects_like(
        &with_faces(document(), near(["0", "0", "1"])),
        &with_faces(document(), wall("inner")),
    );
    selects_like(
        &with_faces(document(), near(["12", "0", "1"])),
        &with_faces(document(), wall("outer")),
    );
    // Wall height clamps to the finite interval: below the section the inner
    // wall is still reached through its rim, ahead of both caps.
    selects_like(
        &with_faces(document(), near(["2", "0", "1/2"])),
        &with_faces(document(), wall("inner")),
    );
}

#[test]
fn points_above_and_below_a_cap_reach_the_cap_and_its_refusal() {
    // The annular caps are not drafted faces; the existing construction owns
    // that refusal, reached only after the unique nearest face is a cap.
    for xyz in [["5", "0", "3"], ["5", "0", "-1"], ["5", "2", "5"]] {
        refuses(
            &with_faces(document(), near(xyz)),
            "draft.parallel-face-normal",
        );
    }
}

#[test]
fn exact_ties_remain_ambiguous_references() {
    // Mid-annulus the two caps are equally near; the hole centre ties the
    // inner wall with the start cap, and a point past the rim ties the outer
    // wall with the end cap.
    for xyz in [
        ["5", "0", "1"],
        ["0", "0", "0"],
        ["-5", "2", "1"],
        ["9", "0", "3"],
    ] {
        refuses(
            &with_faces(document(), near(xyz)),
            "recipe.ambiguous-reference",
        );
    }
}

#[test]
fn saved_near_point_is_reread_after_upstream_height_and_radius_edits() {
    let faces = near(["3", "0", "3/2"]);
    // Height 2: the end cap is nearest. Height 3: the point is inside the band
    // and the inner wall is nearest. The saved point is the same each time.
    let source = with_faces(document(), faces.clone());
    refuses(&source, "draft.parallel-face-normal");
    let mut tall = source.clone();
    tall["steps"][1]["distance"] = json!("3");
    selects_like(&tall, &with_faces(tall.clone(), wall("inner")));
    assert_eq!(tall["steps"][2]["faces"], faces);
    // Inner radius 4 puts the point in the hole, where the inner wall now
    // beats the end cap it lost to at height 2 before the edit.
    let mut wide = source.clone();
    wide["steps"][0]["loops"][1][0]["radius"] = json!("4");
    selects_like(&wide, &with_faces(wide.clone(), wall("inner")));
}
