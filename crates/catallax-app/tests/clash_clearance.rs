//! Roadmap 3.3: distant, near and nested pairs of the public target parts get
//! clash or clear by a stated gap, the same at every thread count.
use catallax_app::queries::{BodyRelation, body_relation};
use catallax_features::evaluate;
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Answer, Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{Isometry, Pose, RPoint3, RVector3, Rotation, Vector3};
use catallax_topology::{BodyId, Model, block, cylinder, transformed_body};
use std::num::NonZeroUsize;

const P2: &str = include_str!("../../../parts/P2-collar.json");
const P3: &str = include_str!("../../../parts/P3-filleted-boss.json");
const P4: &str = include_str!("../../../parts/P4-shelled-tray.json");

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn at(x: &str, y: &str, z: &str) -> Pose {
    Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)))
}
fn shifted(model: &mut Model, body: BodyId, to: [&str; 3]) -> BodyId {
    let map = Isometry::translation(RVector3::new(q(to[0]), q(to[1]), q(to[2])));
    transformed_body(model, body, &map, Clone::clone)
}
/// A quarter or other turn about a coordinate axis through the origin, then a shift.
fn turned(
    model: &mut Model,
    body: BodyId,
    axis: [&str; 3],
    degrees: &str,
    to: [&str; 3],
) -> BodyId {
    let mut map = Isometry::around_axis(
        &RPoint3::new(q("0"), q("0"), q("0")),
        &RVector3::new(q(axis[0]), q(axis[1]), q(axis[2])),
        &q(degrees),
        &mut Decider::default(),
    )
    .unwrap();
    map.translation = map.translation.plus(&Vector3::new(
        q(to[0]).into(),
        q(to[1]).into(),
        q(to[2]).into(),
    ));
    transformed_body(model, body, &map, Clone::clone)
}
fn relation(model: &Model, a: BodyId, b: BodyId) -> Answer<BodyRelation> {
    body_relation(model, a, b, &Accuracy::default(), &mut Decider::default())
        .unwrap_or_else(|e| panic!("{e}"))
}
fn clear(model: &Model, a: BodyId, b: BodyId) -> Enclosure {
    match relation(model, a, b).value {
        BodyRelation::Clear { gap } => gap,
        other => panic!("{other:?}"),
    }
}
fn exact(model: &Model, a: BodyId, b: BodyId, value: &str) {
    let gap = clear(model, a, b);
    assert_eq!((gap.lo, gap.hi), (q(value), q(value)));
}
fn part(text: &str) -> catallax_features::Built {
    evaluate(&Recipe::parse(text).unwrap(), &mut Decider::default()).unwrap()
}

#[test]
fn distant_near_and_turned_copies_of_the_collar_are_clear() {
    // P2 is a ring about y between radii 6 and 12, 24 mm wide in x.
    let mut built = part(P2);
    let a = built.result;
    let model = &mut built.model;
    let far = shifted(model, a, ["500", "0", "0"]);
    let near = shifted(model, a, ["25", "0", "0"]);
    let quarter = turned(model, a, ["0", "0", "1"], "90", ["48", "0", "0"]);
    exact(model, a, far, "476");
    exact(model, a, near, "1");
    exact(model, a, quarter, "24");
}

#[test]
fn near_bosses_are_clear_where_the_boolean_refuses() {
    // P3 is 40 mm long in x. Its torus fillets make a non-coaxial torus pair
    // the Boolean cannot intersect; the proven boxes are 1 mm apart.
    let mut built = part(P3);
    let a = built.result;
    let model = &mut built.model;
    let near = shifted(model, a, ["41", "0", "0"]);
    exact(model, a, near, "1");
}

#[test]
fn a_pin_in_the_collar_bore_and_the_collar_in_the_tray_are_nested() {
    let mut built = part(P2);
    let collar = built.result;
    let model = &mut built.model;
    // Pins along y through the radius 6 bore: radius 5 is clear by 1 mm,
    // radius 7 overlaps the collar.
    for (radius, expected) in [("5", Some("1")), ("7", None)] {
        let pin = cylinder(model, q(radius), q("22"), &at("0", "0", "-5")).unwrap();
        let pin = turned(model, pin, ["1", "0", "0"], "-90", ["0", "0", "0"]);
        match (relation(model, collar, pin).value, expected) {
            (BodyRelation::Clear { gap }, Some(value)) => {
                assert_eq!((gap.lo, gap.hi), (q(value), q(value)));
            }
            (BodyRelation::Clash { overlap }, None) => assert!(overlap.lo > q("0")),
            (other, _) => panic!("radius {radius}: {other:?}"),
        }
    }
    // P2 standing on its end in P4's pocket, 1 mm above the floor and 1 mm
    // from the pocket's long walls. Both parts come from one recipe.
    let mut both: serde_json::Value = serde_json::from_str(P2).unwrap();
    let tray: serde_json::Value = serde_json::from_str(P4).unwrap();
    let mut steps = both["steps"].as_array().unwrap().clone();
    steps.extend(tray["steps"].as_array().unwrap().iter().cloned());
    both["steps"] = serde_json::Value::Array(steps);
    both["result"] = "rounded".into();
    let mut built = part(&both.to_string());
    let collar = built.bodies[&StepId::new("collar").unwrap()];
    let tray = built.bodies[&StepId::new("rounded").unwrap()];
    let model = &mut built.model;
    let standing = turned(model, collar, ["1", "0", "0"], "90", ["0", "0", "3"]);
    exact(model, tray, standing, "1");
    let inside = block(model, [q("10"), q("10"), q("5")], &at("-5", "-5", "4")).unwrap();
    exact(model, tray, inside, "2");
}

#[test]
fn answers_and_work_are_the_same_at_every_thread_count() {
    // The pin runs the overlap Boolean and the gap search; the copy is apart.
    let mut built = part(P2);
    let a = built.result;
    let model = &mut built.model;
    let pin = cylinder(model, q("5"), q("22"), &at("0", "0", "-5")).unwrap();
    let pin = turned(model, pin, ["1", "0", "0"], "-90", ["0", "0", "0"]);
    let quarter = turned(model, a, ["0", "0", "1"], "90", ["48", "0", "0"]);
    let model = &built.model;
    let run = || {
        [pin, quarter].map(|b| {
            let (answer, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                body_relation(model, a, b, &Accuracy::default(), &mut Decider::default())
            });
            let answer = answer.unwrap_or_else(|e| panic!("{e}"));
            let gap = match answer.value {
                BodyRelation::Clear { gap } => gap,
                other => panic!("{other:?}"),
            };
            (gap, answer.marks.len(), work.used)
        })
    };
    let records: Vec<_> = [1, 2, 4]
        .map(|n| {
            set_threads(NonZeroUsize::new(n));
            run()
        })
        .into();
    set_threads(None);
    assert!(records.windows(2).all(|w| w[0] == w[1]), "{records:?}");
    assert_eq!(records[0][0].0.lo, q("1"));
    assert_eq!(records[0][1].0.lo, q("24"));
}
