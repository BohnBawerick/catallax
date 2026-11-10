//! Fillets on the ring where a perpendicular cylinder meets a plane. The ball
//! rolls on the material side of both faces: a hole mouth is convex and its
//! torus lies outside the bore, a blind-hole floor is concave and its torus
//! lies inside, for bosses, holes, counterbores and both axis directions.

use std::f64::consts::PI;

use catallax_boolean::{Containment, classify_point};
use catallax_features::{
    BuildError, Built, evaluate,
    fillet::{RADIUS_TOO_LARGE, radius_limit},
};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface, Torus};
use catallax_topology::Sense;
use serde_json::{Value, json};

const UP: [&str; 4] = ["1", "0", "0", "0"];
// A half turn about x: the cylinder axis points down and its start is on top.
const DOWN: [&str; 4] = ["0", "1", "0", "0"];

fn block() -> Vec<Value> {
    vec![
        json!({"id": "block-section", "op": "sketch",
        "plane": {"rotation": UP, "translation": ["0", "0", "0"]},
        "loops": [[
            {"kind": "line", "id": "s", "to": ["20", "-20"]},
            {"kind": "line", "id": "e", "to": ["20", "20"]},
            {"kind": "line", "id": "n", "to": ["-20", "20"]},
            {"kind": "line", "id": "w", "to": ["-20", "-20"]}
        ]]}),
        json!({"id": "block", "op": "extrude", "profile": "block-section", "distance": "20"}),
    ]
}
fn cylinder(id: &str, radius: &str, height: &str, rotation: [&str; 4], z: &str) -> Value {
    json!({"id": id, "op": "cylinder", "radius": radius, "height": height,
        "pose": {"rotation": rotation, "translation": ["0", "0", z]}})
}
fn boolean(id: &str, kind: &str, target: &str, tool: &str) -> Value {
    json!({"id": id, "op": "boolean", "kind": kind, "target": target, "tool": tool})
}
fn wall(step: &str) -> Value {
    json!({"step": step, "role": {"side": {"segment": "wall"}}})
}
fn cap(step: &str, role: &str) -> Value {
    json!({"step": step, "role": role})
}
fn recipe(mut steps: Vec<Value>, body: &str, radius: &str, between: [Value; 2]) -> Recipe {
    steps.push(
        json!({"id": "round", "op": "fillet", "body": body, "radius": radius,
        "edges": [{"between": between}]}),
    );
    Recipe::parse(
        &json!({"format": "catallax-recipe/3", "name": "ring", "steps": steps, "result": "round"})
            .to_string(),
    )
    .unwrap()
}
fn build(recipe: &Recipe) -> Result<Built, catallax_numbers::Refusal> {
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    });
    evaluate(recipe, &mut dec).map_err(|error| match error {
        BuildError::Geometry { refusal, .. } => refusal,
        BuildError::Recipe(error) => panic!("invalid test recipe: {error:?}"),
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn f(value: &catallax_numbers::Real) -> f64 {
    value.to_f64()
}

/// The rolling ball expected at a ring: its centre circle's height and radius.
struct Ring {
    convex: bool,
    centre_z: f64,
    major: f64,
    minor: f64,
}

/// Checks the one ring strip against the exact rolling ball, the published
/// face's side, and the solid's volume against `volume` within its enclosure.
#[allow(
    clippy::float_cmp,
    reason = "fixture dimensions require exact equality"
)]
fn check(built: &mut Built, ring: &Ring, volume: f64) {
    assert_eq!(built.model.check(built.result), vec![]);
    let [definition] = built.fillets.as_slice() else {
        panic!("expected one strip, got {}", built.fillets.len());
    };
    let Surface::Torus(Torus {
        center,
        axis,
        major,
        minor,
        ..
    }) = &definition.surface
    else {
        panic!("expected a ring torus, got {:?}", definition.surface);
    };
    assert_eq!(definition.convex, ring.convex);
    assert_eq!((f(&center.x), f(&center.y)), (0.0, 0.0));
    assert_eq!(f(&center.z), ring.centre_z);
    assert_eq!((f(&axis.x), f(&axis.y)), (0.0, 0.0));
    assert_eq!(f(major), ring.major);
    assert_eq!(f(minor), ring.minor);
    // A convex round keeps the material in the tube, so the face follows the
    // torus's outward normal; a concave one adds material outside the tube.
    let face = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|face| {
            built
                .model
                .face_origins(*face)
                .iter()
                .any(|origin| matches!(origin.role, Role::Blend { .. }))
        })
        .unwrap();
    let sense = built.model.face(face).sense;
    assert_eq!(sense == Sense::Same, ring.convex, "torus face side");
    let measured = catallax_measure::volume_area_centroid(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .volume;
    let (lo, hi) = (measured.lo.to_f64_nearest(), measured.hi.to_f64_nearest());
    assert!(
        lo - 1e-9 <= volume && volume <= hi + 1e-9 && hi - lo < 1e-6,
        "volume [{lo}, {hi}] does not hold {volume}"
    );
}

/// Classifies the point at radius `x` on the positive x axis and height `z`.
fn at(built: &mut Built, x: &str, z: &str) -> Containment {
    let point = RPoint3::new(q(x), Rational::zero(), q(z));
    match classify_point(
        &mut built.model,
        built.result,
        &point,
        &mut Decider::default(),
    ) {
        Outcome::Proven(containment) => containment,
        other => panic!("unclassified {point:?}: {other:?}"),
    }
}

/// Points beside the round at 45 degrees from the ball centre `(s, z)`:
/// 0.99 radius from the centre and 1.018 radius from it, on the corner side.
fn beside(s: i64, z: i64, toward: [i64; 2]) -> [(String, String); 2] {
    let near = |offset: [i64; 2]| {
        (
            format!("{}/100", 100 * s + toward[0] * offset[0]),
            format!("{}/100", 100 * z + toward[1] * offset[1]),
        )
    };
    [near([71, 69]), near([72, 72])]
}

// Volume a radius-r round takes from or adds to a ring of radius `wall`: the
// corner region (square minus quarter disc) turned about the axis. Its
// centroid lies (10 - 3 pi) r / (3 (4 - pi)) from the corner, outward from
// the axis when `outward`.
fn ring_volume(wall: f64, r: f64, outward: bool) -> f64 {
    let area = r * r * (1.0 - PI / 4.0);
    let moment = r * r * r * (10.0 - 3.0 * PI) / 12.0;
    2.0 * PI * (wall * area + if outward { moment } else { -moment })
}
const BLOCK: f64 = 40.0 * 40.0 * 20.0;

#[test]
fn boss_top_is_convex_with_the_torus_inside_the_wall() {
    let mut built = build(&recipe(
        vec![cylinder("boss", "3", "10", UP, "0")],
        "boss",
        "1",
        [cap("boss", "end"), wall("boss")],
    ))
    .unwrap();
    let ring = Ring {
        convex: true,
        centre_z: 9.0,
        major: 2.0,
        minor: 1.0,
    };
    check(
        &mut built,
        &ring,
        PI * 9.0 * 10.0 - ring_volume(3.0, 1.0, false),
    );
    let [inside, outside] = beside(2, 9, [1, 1]);
    assert_eq!(at(&mut built, &inside.0, &inside.1), Containment::Inside);
    assert_eq!(at(&mut built, &outside.0, &outside.1), Containment::Outside);
}

fn through_hole(rotation: [&str; 4], z: &str) -> Vec<Value> {
    let mut steps = block();
    steps.push(cylinder("hole", "3", "30", rotation, z));
    steps.push(boolean("cut", "subtract", "block", "hole"));
    steps
}

/// The mouth of a through hole is convex: the ball sits in the material
/// outside the bore, so the torus's major radius is the hole radius plus the
/// fillet radius.
fn mouth(steps: Vec<Value>, block_face: &str, top: bool) {
    let mut built = build(&recipe(
        steps,
        "cut",
        "1",
        [cap("block", block_face), wall("hole")],
    ))
    .unwrap();
    let (z, toward) = if top { (19, [-1, 1]) } else { (1, [-1, -1]) };
    let ring = Ring {
        convex: true,
        centre_z: z as f64,
        major: 4.0,
        minor: 1.0,
    };
    check(
        &mut built,
        &ring,
        BLOCK - PI * 9.0 * 20.0 - ring_volume(3.0, 1.0, true),
    );
    let [inside, outside] = beside(4, z, toward);
    assert_eq!(at(&mut built, &inside.0, &inside.1), Containment::Inside);
    assert_eq!(at(&mut built, &outside.0, &outside.1), Containment::Outside);
    // The old torus of major radius 2 rounded into the bore; the bore stays empty.
    let rim = if top { "199/10" } else { "1/10" };
    assert_eq!(at(&mut built, "5/2", rim), Containment::Outside);
}

#[test]
fn hole_mouth_is_convex_with_the_torus_outside_the_bore() {
    // The reported recipe: the hole axis points up through the top face.
    mouth(through_hole(UP, "-5"), "end", true);
}

#[test]
fn hole_mouth_keeps_its_side_when_the_hole_axis_points_down() {
    mouth(through_hole(DOWN, "25"), "end", true);
}

#[test]
fn hole_mouth_on_the_bottom_face_is_convex() {
    mouth(through_hole(UP, "-5"), "start", false);
    mouth(through_hole(DOWN, "25"), "start", false);
}

/// The floor ring of a blind hole is concave: the ball sits in the bore, so the
/// torus's major radius is the hole radius minus the fillet radius.
fn floor(rotation: [&str; 4], z: &str, floor_role: &str, radius: &str) -> Built {
    let mut steps = block();
    steps.push(cylinder("hole", "3", "20", rotation, z));
    steps.push(boolean("cut", "subtract", "block", "hole"));
    let built = build(&recipe(
        steps,
        "cut",
        radius,
        [cap("hole", floor_role), wall("hole")],
    ));
    built.unwrap_or_else(|refusal| panic!("{refusal:?}"))
}

#[test]
fn blind_hole_floor_is_concave_with_the_torus_inside_the_bore() {
    // The reported recipe: drilled down from z = 30 to the floor at z = 10.
    for (rotation, z, role, centre_z, toward) in [
        (UP, "10", "start", 11, [1, -1]),
        // The same hole drilled up from below: floor at z = 10 facing down.
        (DOWN, "10", "start", 9, [1, 1]),
    ] {
        let mut built = floor(rotation, z, role, "1");
        let ring = Ring {
            convex: false,
            centre_z: centre_z as f64,
            major: 2.0,
            minor: 1.0,
        };
        check(
            &mut built,
            &ring,
            BLOCK - PI * 9.0 * 10.0 + ring_volume(3.0, 1.0, false),
        );
        let [air, material] = beside(2, centre_z, toward);
        assert_eq!(at(&mut built, &air.0, &air.1), Containment::Outside);
        assert_eq!(
            at(&mut built, &material.0, &material.1),
            Containment::Inside
        );
    }
}

#[test]
fn a_floor_round_wider_than_the_hole_is_refused() {
    let mut steps = block();
    steps.push(cylinder("hole", "3", "20", UP, "10"));
    steps.push(boolean("cut", "subtract", "block", "hole"));
    let Err(refusal) = build(&recipe(
        steps,
        "cut",
        "7/2",
        [cap("hole", "start"), wall("hole")],
    )) else {
        panic!("a radius 7/2 ball does not fit a radius 3 bore");
    };
    // The ball centre circle would have radius 3 - 7/2 < 0: no ball wider
    // than the bore rolls round its floor.
    assert_eq!(refusal.code, RADIUS_TOO_LARGE, "{refusal:?}");
    assert_eq!(radius_limit(&refusal), Some(q("3")));
}

#[test]
fn boss_base_is_concave_with_the_torus_outside_the_wall() {
    for (rotation, z) in [(UP, "10"), (DOWN, "30")] {
        let mut steps = block();
        steps.push(cylinder("boss", "3", "20", rotation, z));
        steps.push(boolean("joined", "union", "block", "boss"));
        let mut built = build(&recipe(
            steps,
            "joined",
            "1",
            [cap("block", "end"), wall("boss")],
        ))
        .unwrap();
        let ring = Ring {
            convex: false,
            centre_z: 21.0,
            major: 4.0,
            minor: 1.0,
        };
        check(
            &mut built,
            &ring,
            BLOCK + PI * 9.0 * 10.0 + ring_volume(3.0, 1.0, true),
        );
        let [air, material] = beside(4, 21, [-1, -1]);
        assert_eq!(at(&mut built, &air.0, &air.1), Containment::Outside);
        assert_eq!(
            at(&mut built, &material.0, &material.1),
            Containment::Inside
        );
    }
}

/// A radius 6 counterbore 5 deep over a radius 3 through hole. Its floor is
/// wider than two fillet radii, so a floor round stays clear of the hole.
fn counterbore(rotation: [&str; 4]) -> Vec<Value> {
    let mut steps = block();
    let (bore_z, hole_z) = if rotation == UP {
        ("15", "-5")
    } else {
        ("25", "25")
    };
    steps.push(cylinder("bore", "6", "10", rotation, bore_z));
    steps.push(boolean("bored", "subtract", "block", "bore"));
    steps.push(cylinder("hole", "3", "30", rotation, hole_z));
    steps.push(boolean("cut", "subtract", "bored", "hole"));
    steps
}
const COUNTERBORE: f64 = BLOCK - PI * 36.0 * 5.0 - PI * 9.0 * 15.0;

#[test]
fn counterbore_rings_take_their_own_sides() {
    for rotation in [UP, DOWN] {
        let floor_role = if rotation == UP { "start" } else { "end" };
        // The bore mouth on top is convex, outside the bore.
        let mut built = build(&recipe(
            counterbore(rotation),
            "cut",
            "1",
            [cap("block", "end"), wall("bore")],
        ))
        .unwrap();
        let ring = Ring {
            convex: true,
            centre_z: 19.0,
            major: 7.0,
            minor: 1.0,
        };
        check(&mut built, &ring, COUNTERBORE - ring_volume(6.0, 1.0, true));
        let [inside, outside] = beside(7, 19, [-1, 1]);
        assert_eq!(at(&mut built, &inside.0, &inside.1), Containment::Inside);
        assert_eq!(at(&mut built, &outside.0, &outside.1), Containment::Outside);
        // The bore's floor meets its wall in a concave ring, inside the bore.
        let mut built = build(&recipe(
            counterbore(rotation),
            "cut",
            "1",
            [cap("bore", floor_role), wall("bore")],
        ))
        .unwrap();
        let ring = Ring {
            convex: false,
            centre_z: 16.0,
            major: 5.0,
            minor: 1.0,
        };
        check(
            &mut built,
            &ring,
            COUNTERBORE + ring_volume(6.0, 1.0, false),
        );
        let [air, material] = beside(5, 16, [1, -1]);
        assert_eq!(at(&mut built, &air.0, &air.1), Containment::Outside);
        assert_eq!(
            at(&mut built, &material.0, &material.1),
            Containment::Inside
        );
        // The through hole's mouth in that floor is convex, outside the hole.
        let mut built = build(&recipe(
            counterbore(rotation),
            "cut",
            "1",
            [cap("bore", floor_role), wall("hole")],
        ))
        .unwrap();
        let ring = Ring {
            convex: true,
            centre_z: 14.0,
            major: 4.0,
            minor: 1.0,
        };
        check(&mut built, &ring, COUNTERBORE - ring_volume(3.0, 1.0, true));
        let [inside, outside] = beside(4, 14, [-1, 1]);
        assert_eq!(at(&mut built, &inside.0, &inside.1), Containment::Inside);
        assert_eq!(at(&mut built, &outside.0, &outside.1), Containment::Outside);
    }
}

/// The same rings made by the hole operation, drilled down from the top face.
fn drilled(depth: Value, head: Option<Value>) -> Vec<Value> {
    let mut steps = block();
    let mut hole = json!({"id": "drill", "op": "hole", "target": "block",
        "at": [{"id": "h", "point": ["0", "0", "20"], "direction": ["0", "0", "-1"]}],
        "diameter": "6", "depth": depth, "tip": "flat"});
    if let Some(head) = head {
        hole["head"] = head;
    }
    steps.push(hole);
    steps
}
fn part(name: &str) -> Value {
    json!({"step": "drill", "role": {"hole": {"at": "h", "part": name}}})
}

#[test]
fn hole_operation_rings_take_their_own_sides() {
    let counterbore = || Some(json!({"counterbore": {"diameter": "12", "depth": "5"}}));
    // Each case: the hole, the two faces, the ball, the volume and the
    // direction from the ball centre to the sharp corner.
    let cases = [
        (
            drilled(json!("through"), None),
            [cap("block", "end"), part("wall")],
            (true, 19, 4),
            BLOCK - PI * 9.0 * 20.0 - ring_volume(3.0, 1.0, true),
            [-1, 1],
        ),
        (
            drilled(json!({"blind": "10"}), None),
            [part("floor"), part("wall")],
            (false, 11, 2),
            BLOCK - PI * 9.0 * 10.0 + ring_volume(3.0, 1.0, false),
            [1, -1],
        ),
        (
            drilled(json!("through"), counterbore()),
            [cap("block", "end"), part("counterbore_wall")],
            (true, 19, 7),
            COUNTERBORE - ring_volume(6.0, 1.0, true),
            [-1, 1],
        ),
        (
            drilled(json!("through"), counterbore()),
            [part("counterbore_floor"), part("counterbore_wall")],
            (false, 16, 5),
            COUNTERBORE + ring_volume(6.0, 1.0, false),
            [1, -1],
        ),
        (
            drilled(json!("through"), counterbore()),
            [part("counterbore_floor"), part("wall")],
            (true, 14, 4),
            COUNTERBORE - ring_volume(3.0, 1.0, true),
            [-1, 1],
        ),
    ];
    for (steps, between, (convex, centre_z, major), volume, toward) in cases {
        let mut built = build(&recipe(steps, "drill", "1", between)).unwrap();
        let ring = Ring {
            convex,
            centre_z: centre_z as f64,
            major: major as f64,
            minor: 1.0,
        };
        check(&mut built, &ring, volume);
        // The point 0.99 radius from the centre is in the ball, the other is
        // in the corner: material kept by a convex round, added by a concave.
        let [ball, corner] = beside(major, centre_z, toward);
        let (ball_side, corner_side) = if convex {
            (Containment::Inside, Containment::Outside)
        } else {
            (Containment::Outside, Containment::Inside)
        };
        assert_eq!(at(&mut built, &ball.0, &ball.1), ball_side);
        assert_eq!(at(&mut built, &corner.0, &corner.1), corner_side);
    }
}
