//! Two rings of a two-face coaxial shell rounded at the inscribed radius.
//!
//! The plane z = 1 keeps the top of the torus R = 6, r = 2: a segment of
//! height 1 turned about the axis. Its inscribed circle has radius 1/2 about
//! (6, 3/2), so rounding both rings by 1/2 consumes the plane and the torus
//! wall and leaves the whole torus R = 6, r = 1/2 about z = 3/2, of volume
//! 3 pi^2. OpenCASCADE refuses this request.

use std::f64::consts::PI;

use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn real(value: &str) -> Real {
    Real::from(q(value))
}
fn same(a: &Real, b: &Real) -> bool {
    matches!(decider().sign(&(a - b)), Outcome::ProvenCoincident(_))
}

fn volume_of(built: &catallax_features::Built) -> (f64, f64) {
    let found = volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: q("1/1000000000"),
            absolute: q("1/1000000000"),
        },
        &mut decider(),
    )
    .unwrap();
    assert!(found.marks.is_empty());
    (
        found.value.lo.to_f64_nearest(),
        found.value.hi.to_f64_nearest(),
    )
}

fn assert_holds((lo, hi): (f64, f64), want: f64) {
    assert!(
        lo - 1e-7 <= want && want <= hi + 1e-7,
        "[{lo}, {hi}] does not hold {want}"
    );
}

// Every face lies on a torus of tube radius 1/2 about z = 3/2.
fn all_inscribed(built: &catallax_features::Built, faces: usize) {
    assert!(built.model.check(built.result).is_empty());
    let found = built.model.faces_of_body(built.result);
    assert_eq!(found.len(), faces);
    for face in found {
        let Surface::Torus(torus) = built.model.geom().surface(built.model.face(face).surface)
        else {
            continue;
        };
        assert!(same(&torus.minor, &real("1/2")));
        assert!(same(&torus.center.z, &real("3/2")));
    }
}

fn recipe(radius: &str, edges: &str) -> Recipe {
    let up = r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#;
    let block = r#"{"rotation":["1","0","0","0"],"translation":["-20","-20","1"]}"#;
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/4","name":"inscribed","result":"round","steps":[
        {{"id":"t","op":"torus","major":"6","minor":"2","pose":{up}}},
        {{"id":"b","op":"block","size":["40","40","10"],"pose":{block}}},
        {{"id":"j","op":"boolean","kind":"intersect","target":"t","tool":"b"}},
        {{"id":"round","op":"fillet","body":"j","radius":"{radius}","edges":{edges}}}]}}"#
    ))
    .unwrap()
}
const BOTH: &str = r#"[{"between":[{"step":"b","role":"start"},{"step":"t","role":{"side":{"segment":"wall"}}}]}]"#;

#[test]
fn both_rings_at_the_inscribed_radius_leave_the_whole_torus() {
    let mut built = evaluate(&recipe("1/2", BOTH), &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 2, "two half tubes");
    for face in &faces {
        let Surface::Torus(torus) = built.model.geom().surface(built.model.face(*face).surface)
        else {
            panic!("each half lies on the rolling torus");
        };
        assert!(same(&torus.major, &real("6")));
        assert!(same(&torus.minor, &real("1/2")));
        assert!(same(&torus.center.z, &real("3/2")));
        let [origin] = built.model.face_origins(*face) else {
            panic!("one blend origin per half");
        };
        assert!(matches!(origin.role, Role::Blend { .. }));
    }
    assert_eq!(built.fillets.len(), 2);
    // Both input faces are consumed, each by both blends.
    assert_eq!(built.consumed.len(), 2);
    for record in &built.consumed {
        assert_eq!(record.step.as_str(), "round");
        let mut by: Vec<_> = built.fillets.iter().map(|d| d.origin.clone()).collect();
        by.sort();
        assert_eq!(record.by, by);
    }
    let found = volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: q("1/1000000000"),
            absolute: q("1/1000000000"),
        },
        &mut decider(),
    )
    .unwrap();
    let want = 3.0 * PI * PI;
    let (lo, hi) = (
        found.value.lo.to_f64_nearest(),
        found.value.hi.to_f64_nearest(),
    );
    assert!(
        lo - 1e-7 <= want && want <= hi + 1e-7,
        "[{lo}, {hi}] does not hold {want}"
    );
    // The tube's centre and points just inside it are material; the old
    // plane contact and the old top beyond the tube are not.
    for (p, want) in [
        (["6", "0", "3/2"], Containment::Inside),
        (["129/20", "0", "3/2"], Containment::Inside),
        (["0", "-111/20", "3/2"], Containment::Inside),
        (["7", "0", "11/10"], Containment::Outside),
        (["5", "0", "11/10"], Containment::Outside),
        (["6", "0", "21/10"], Containment::Outside),
    ] {
        let [x, y, z] = p.map(q);
        let found = classify_point(
            &mut built.model,
            built.result,
            &RPoint3::new(x, y, z),
            &mut decider(),
        );
        assert!(
            matches!(found, Outcome::Proven(c) if c == want),
            "{p:?}: wanted {want:?}, got {found:?}"
        );
    }
}

#[test]
fn a_radius_off_the_limit_keeps_the_regular_answer() {
    // Below the limit the two strips stay apart and both faces survive.
    let built = evaluate(&recipe("2/5", BOTH), &mut decider()).unwrap();
    assert!(built.consumed.is_empty());
    assert_eq!(built.model.faces_of_body(built.result).len(), 4);
    // Past it no ball of the radius fits the section.
    let error = evaluate(&recipe("3/5", BOTH), &mut decider()).unwrap_err();
    let catallax_features::BuildError::Geometry { refusal, .. } = error else {
        panic!("a geometric refusal, got {error:?}");
    };
    assert_eq!(refusal.code, "fillet.radius-too-large");
}

/// A void of the same section inside a block: the concave rings fill all of
/// it but the torus, so the block keeps a toroidal void of volume 3 pi^2.
#[test]
fn concave_rings_at_the_inscribed_radius_leave_a_toroidal_void() {
    let up = r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#;
    let cap = r#"{"rotation":["1","0","0","0"],"translation":["-20","-20","1"]}"#;
    let stock = r#"{"rotation":["1","0","0","0"],"translation":["-20","-20","-5"]}"#;
    let recipe = Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/4","name":"void","result":"round","steps":[
        {{"id":"t","op":"torus","major":"6","minor":"2","pose":{up}}},
        {{"id":"b","op":"block","size":["40","40","10"],"pose":{cap}}},
        {{"id":"j","op":"boolean","kind":"intersect","target":"t","tool":"b"}},
        {{"id":"c","op":"block","size":["40","40","20"],"pose":{stock}}},
        {{"id":"v","op":"boolean","kind":"subtract","target":"c","tool":"j"}},
        {{"id":"round","op":"fillet","body":"v","radius":"1/2","edges":{BOTH}}}]}}"#
    ))
    .unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    all_inscribed(&built, 8);
    assert_eq!(built.consumed.len(), 2);
    assert_holds(volume_of(&built), 32000.0 - 3.0 * PI * PI);
    for (p, want) in [
        (["6", "0", "3/2"], Containment::Outside),
        (["7", "0", "11/10"], Containment::Inside),
        (["6", "0", "21/10"], Containment::Inside),
    ] {
        let [x, y, z] = p.map(q);
        let found = classify_point(
            &mut built.model,
            built.result,
            &RPoint3::new(x, y, z),
            &mut decider(),
        );
        assert!(
            matches!(found, Outcome::Proven(c) if c == want),
            "{p:?}: wanted {want:?}, got {found:?}"
        );
    }
}

/// Two such sections, as separate lumps, rounded together by `"all"`.
#[test]
fn every_inscribed_shell_of_a_selection_closes_its_own_torus() {
    let up = r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#;
    let aside = r#"{"rotation":["1","0","0","0"],"translation":["20","0","0"]}"#;
    let cap = r#"{"rotation":["1","0","0","0"],"translation":["-20","-20","1"]}"#;
    let recipe = Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/4","name":"pair","result":"round","steps":[
        {{"id":"t","op":"torus","major":"6","minor":"2","pose":{up}}},
        {{"id":"s","op":"torus","major":"6","minor":"2","pose":{aside}}},
        {{"id":"b","op":"block","size":["60","40","10"],"pose":{cap}}},
        {{"id":"j","op":"boolean","kind":"intersect","target":"t","tool":"b"}},
        {{"id":"k","op":"boolean","kind":"intersect","target":"s","tool":"b"}},
        {{"id":"u","op":"boolean","kind":"union","target":"j","tool":"k"}},
        {{"id":"round","op":"fillet","body":"u","radius":"1/2","edges":"all"}}]}}"#
    ))
    .unwrap();
    let built = evaluate(&recipe, &mut decider()).unwrap();
    all_inscribed(&built, 4);
    assert_eq!(built.fillets.len(), 4);
    assert_eq!(built.consumed.len(), 3, "the shared block origin coalesces");
    assert_holds(volume_of(&built), 6.0 * PI * PI);
}
