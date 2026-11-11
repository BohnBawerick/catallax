//! A fillet radius proved too large for its edge refuses with its own code and
//! the bound past which every radius is too large. A radius just under the
//! bound never gets that code.

use catallax_features::{
    BuildError, evaluate,
    fillet::{RADIUS_TOO_LARGE, radius_limit},
};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

const LEVEL: [&str; 4] = ["1", "0", "0", "0"];
/// A quarter turn about x, which lays a cylinder along y.
const FLAT: [&str; 4] = ["1", "1", "0", "0"];
/// A half turn about y of `2 atan(7/24)`.
const TILT: [&str; 4] = ["24/25", "0", "7/25", "0"];

fn cylinder(id: &str, radius: &str, height: &str, rotation: [&str; 4], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"cylinder","radius":radius,"height":height,
        "pose":{"rotation":rotation,"translation":at}})
}
fn upright(id: &str, radius: &str, height: &str) -> Value {
    cylinder(id, radius, height, LEVEL, ["0", "0", "0"])
}
/// A rectangle extruded from a plane through `at`.
fn block(
    id: &str,
    [west, east, south, north]: [&str; 4],
    distance: &str,
    rotation: [&str; 4],
    at: [&str; 3],
) -> [Value; 2] {
    let section = format!("{id}-section");
    [
        json!({"id":section,"op":"sketch",
            "plane":{"rotation":rotation,"translation":at},
            "loops":[[{"kind":"line","id":"s","to":[east,south]},
                {"kind":"line","id":"e","to":[east,north]},
                {"kind":"line","id":"n","to":[west,north]},
                {"kind":"line","id":"w","to":[west,south]}]]}),
        json!({"id":id,"op":"extrude","profile":section,"distance":distance}),
    ]
}
fn boolean(id: &str, kind: &str, target: &str, tool: &str) -> Value {
    json!({"id":id,"op":"boolean","kind":kind,"target":target,"tool":tool})
}
fn wall(step: &str) -> Value {
    json!({"step":step,"role":{"side":{"segment":"wall"}}})
}
fn role(step: &str, role: Value) -> Value {
    json!({"step":step,"role":role})
}

/// A body and the one edge to round.
struct Part {
    steps: Vec<Value>,
    body: &'static str,
    edge: [Value; 2],
}
impl Part {
    fn rounded(&self, radius: &str) -> Recipe {
        let mut steps = self.steps.clone();
        steps.push(
            json!({"id":"round","op":"fillet","body":self.body,"radius":radius,
            "edges":[{"between":self.edge}]}),
        );
        Recipe::parse(
            &json!({"format":"catallax-recipe/3","name":"radius limit",
                "steps":steps,"result":"round"})
            .to_string(),
        )
        .unwrap()
    }
}

/// A plate 60 by 40 by 10, rounded along its top front edge.
fn plate() -> Part {
    let mut steps = Vec::new();
    steps.extend(block(
        "plate",
        ["-30", "30", "-20", "20"],
        "10",
        LEVEL,
        ["0", "0", "0"],
    ));
    Part {
        steps,
        body: "plate",
        edge: [
            role("plate", json!("end")),
            role("plate", json!({"side":{"segment":"s"}})),
        ],
    }
}
/// The same plate rounded along the inside corner of a 6 high wall on it.
fn ledge() -> Part {
    let mut steps = Vec::new();
    steps.extend(block(
        "plate",
        ["-30", "30", "-20", "20"],
        "10",
        LEVEL,
        ["0", "0", "0"],
    ));
    steps.extend(block(
        "wall",
        ["-30", "30", "10", "20"],
        "6",
        LEVEL,
        ["0", "0", "10"],
    ));
    steps.push(boolean("joined", "union", "plate", "wall"));
    Part {
        steps,
        body: "joined",
        edge: [
            role("plate", json!("end")),
            role("wall", json!({"side":{"segment":"s"}})),
        ],
    }
}
/// A cylinder rounded along its top rim.
fn rim(radius: &str, height: &str) -> Part {
    Part {
        steps: vec![upright("boss", radius, height)],
        body: "boss",
        edge: [role("boss", json!("end")), wall("boss")],
    }
}
/// A tube of radius 10 with a bore of 7, rounded along its top outer rim.
fn tube() -> Part {
    Part {
        steps: vec![
            upright("outer", "10", "10"),
            cylinder("bore", "7", "12", LEVEL, ["0", "0", "-1"]),
            boolean("tube", "subtract", "outer", "bore"),
        ],
        body: "tube",
        edge: [role("outer", json!("end")), wall("outer")],
    }
}
/// A block with a through hole of radius 3, rounded along the hole's mouth.
fn hole_mouth() -> Part {
    let mut steps = Vec::new();
    steps.extend(block(
        "block",
        ["-20", "20", "-20", "20"],
        "20",
        LEVEL,
        ["0", "0", "0"],
    ));
    steps.push(cylinder("hole", "3", "30", LEVEL, ["0", "0", "-5"]));
    steps.push(boolean("cut", "subtract", "block", "hole"));
    Part {
        steps,
        body: "cut",
        edge: [role("block", json!("end")), wall("hole")],
    }
}
/// A cylinder of radius 5 with a flat at y = 3, rounded along one flat edge.
fn d_shaft() -> Part {
    let mut steps = vec![upright("a", "5", "12")];
    steps.extend(block(
        "block",
        ["-12", "12", "-10", "3"],
        "12",
        LEVEL,
        ["0", "0", "0"],
    ));
    steps.push(boolean("joined", "intersect", "a", "block"));
    Part {
        steps,
        body: "joined",
        edge: [wall("a"), role("block", json!({"side":{"segment":"n"}}))],
    }
}
/// The lens of two radius 5 cylinders 6 apart.
fn lens() -> Part {
    Part {
        steps: vec![
            cylinder("a", "5", "12", LEVEL, ["-3", "0", "0"]),
            cylinder("b", "5", "12", LEVEL, ["3", "0", "0"]),
            boolean("joined", "intersect", "a", "b"),
        ],
        body: "joined",
        edge: [wall("a"), wall("b")],
    }
}
/// A radius 5 cylinder cut off by a tilted plane, rounded along the ellipse.
fn oblique() -> Part {
    let mut steps = vec![upright("wall", "5", "16")];
    steps.extend(block(
        "cut",
        ["-30", "30", "-30", "30"],
        "30",
        TILT,
        ["0", "0", "8"],
    ));
    steps.push(boolean("joined", "subtract", "wall", "cut"));
    Part {
        steps,
        body: "joined",
        edge: [wall("wall"), role("cut", json!("start"))],
    }
}
/// A bore of radius 4 through a radius 10 cylinder, their axes 1 apart,
/// rounded along the ring at one mouth.
fn ring() -> Part {
    Part {
        steps: vec![
            cylinder("wide", "10", "60", FLAT, ["0", "30", "0"]),
            cylinder("wall", "4", "30", LEVEL, ["1", "0", "0"]),
            boolean("joined", "subtract", "wide", "wall"),
        ],
        body: "joined",
        edge: [wall("wide"), wall("wall")],
    }
}
/// A radius 4 cylinder whose axis runs 8 from a radius 10 one's, so the two
/// cross in one loop, with the wide one cut out of it.
fn looped() -> Part {
    Part {
        steps: vec![
            cylinder("wide", "10", "60", FLAT, ["0", "30", "0"]),
            cylinder("wall", "4", "60", LEVEL, ["8", "0", "-30"]),
            boolean("joined", "subtract", "wall", "wide"),
        ],
        body: "joined",
        edge: [wall("wide"), wall("wall")],
    }
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
/// The refusal code and bound, or `None` when the fillet builds.
fn outcome(recipe: &Recipe) -> Option<(String, Option<Rational>)> {
    match evaluate(recipe, &mut decider()) {
        Ok(built) => {
            assert_eq!(built.model.check(built.result), vec![]);
            None
        }
        Err(BuildError::Geometry { refusal, .. }) => {
            Some((refusal.code.to_string(), radius_limit(&refusal)))
        }
        Err(other) => panic!("expected a geometric refusal: {other}"),
    }
}
#[track_caller]
fn too_large(part: &Part, radius: &str, limit: Option<&str>) {
    assert_eq!(
        outcome(&part.rounded(radius)),
        Some((RADIUS_TOO_LARGE.to_string(), limit.map(q))),
        "radius {radius}"
    );
}
/// A radius under the bound builds, or refuses for another reason.
#[track_caller]
fn not_too_large(part: &Part, radius: &str) {
    if let Some((code, _)) = outcome(&part.rounded(radius)) {
        assert_ne!(code, RADIUS_TOO_LARGE, "radius {radius}");
    }
}
#[track_caller]
fn builds(part: &Part, radius: &str) {
    assert_eq!(outcome(&part.rounded(radius)), None, "radius {radius}");
}

// The side face is 10 high, so a ball past 10 cannot touch it.
#[test]
fn a_plate_edge_takes_less_than_the_plate_thickness() {
    too_large(&plate(), "21/2", Some("10"));
    too_large(&plate(), "40", Some("10"));
    not_too_large(&plate(), "10");
    not_too_large(&plate(), "19/2");
    builds(&plate(), "2");
}

// A concave edge: the wall is 6 high.
#[test]
fn an_inside_corner_takes_less_than_its_wall() {
    too_large(&ledge(), "13/2", Some("6"));
    not_too_large(&ledge(), "11/2");
}

// A ball inside a cylinder must be narrower than it.
#[test]
fn a_rim_takes_less_than_its_cylinder_radius() {
    too_large(&rim("4", "10"), "9/2", Some("4"));
    too_large(&rim("4", "10"), "4", Some("4"));
    builds(&rim("4", "10"), "7/2");
}

// The wall below the rim is 3 high.
#[test]
fn a_rim_takes_less_than_its_wall_height() {
    too_large(&rim("6", "3"), "7/2", Some("3"));
    not_too_large(&rim("6", "3"), "5/2");
}

// The tube's wall is 3 thick, so a larger ball's contact lies over the bore.
#[test]
fn a_tube_rim_takes_less_than_its_wall_thickness() {
    too_large(&tube(), "7/2", Some("3"));
    not_too_large(&tube(), "5/2");
}

// A hole's mouth is convex, so the ball sits in the material and any radius
// can roll round it.
#[test]
fn a_hole_mouth_is_never_called_too_large() {
    not_too_large(&hole_mouth(), "4");
}

// The offset line meets the offset circle until r = 4.
#[test]
fn a_flat_on_a_shaft_takes_less_than_its_offset_reach() {
    too_large(&d_shaft(), "41/10", Some("4"));
    too_large(&d_shaft(), "6", Some("4"));
    too_large(&d_shaft(), "4", Some("4"));
    not_too_large(&d_shaft(), "39/10");
}

// The offset circles of radius 5 - r about centres 6 apart touch at r = 2.
#[test]
fn a_lens_takes_less_than_half_its_thickness() {
    too_large(&lens(), "21/10", Some("2"));
    too_large(&lens(), "2", Some("2"));
    not_too_large(&lens(), "19/10");
    builds(&lens(), "1");
}

// The ball sits inside the radius 5 cylinder.
#[test]
fn an_oblique_rim_takes_less_than_its_cylinder_radius() {
    too_large(&oblique(), "11/2", Some("5"));
    too_large(&oblique(), "5", Some("5"));
    builds(&oblique(), "9/2");
}

// The offsets 10 - r and 4 + r with axes 1 apart stop passing through each
// other at r = 5/2.
#[test]
fn a_ring_takes_less_than_where_its_offsets_stop_crossing() {
    too_large(&ring(), "13/5", Some("5/2"));
    too_large(&ring(), "10", Some("5/2"));
    too_large(&ring(), "5/2", Some("5/2"));
    builds(&ring(), "2");
}

// The offsets 10 + r and 4 - r with axes 8 apart cross in one loop until
// their difference 6 + 2r reaches 8.
#[test]
fn a_loop_takes_less_than_where_its_offsets_stop_crossing() {
    too_large(&looped(), "11/10", Some("1"));
    too_large(&looped(), "1", Some("1"));
    not_too_large(&looped(), "9/10");
}

#[test]
fn other_refusals_carry_no_bound() {
    let refusal =
        catallax_numbers::Refusal::new("fillet.unsupported-overflow", "radius limit 3 mm");
    assert_eq!(radius_limit(&refusal), None);
}

#[test]
fn a_rational_contact_bound_is_exact_even_when_the_legacy_hint_is_rounded() {
    let error = evaluate(&rim("1/3", "10").rounded("1/2"), &mut decider()).unwrap_err();
    let BuildError::Geometry { refusal, .. } = error else {
        panic!("expected geometry refusal: {error:?}")
    };
    assert_eq!(refusal.code, RADIUS_TOO_LARGE);
    assert_eq!(radius_limit(&refusal), Some(q("333/1000")));
    assert!(
        refusal
            .detail
            .contains("contact upper bound exactly 1/3 mm"),
        "{refusal:?}"
    );
    assert!(refusal.detail.contains("maximum unknown"), "{refusal:?}");
    assert!(
        refusal.detail.contains("limiting contact supports"),
        "{refusal:?}"
    );
    assert!(refusal.detail.contains(r#""step":"boss""#), "{refusal:?}");
}

#[test]
fn a_trim_bound_names_the_face_the_contact_leaves() {
    let error = evaluate(&plate().rounded("21/2"), &mut decider()).unwrap_err();
    let BuildError::Geometry { refusal, .. } = error else {
        panic!("expected geometry refusal: {error:?}")
    };
    assert_eq!(refusal.code, RADIUS_TOO_LARGE);
    assert_eq!(radius_limit(&refusal), Some(q("10")));
    assert!(
        refusal.detail.contains("contact leaves face"),
        "{refusal:?}"
    );
    assert!(refusal.detail.contains(r#""segment":"s""#), "{refusal:?}");
    assert!(
        refusal.detail.contains("contact upper bound exactly 10 mm"),
        "{refusal:?}"
    );
    assert!(refusal.detail.contains("maximum unknown"), "{refusal:?}");
}

#[test]
fn generator_bounds_name_both_supports_and_keep_equality_strict() {
    for radius in ["4", "41/10"] {
        let error = evaluate(&d_shaft().rounded(radius), &mut decider()).unwrap_err();
        let BuildError::Geometry { refusal, .. } = error else {
            panic!("expected geometry refusal: {error:?}")
        };
        assert_eq!(refusal.code, RADIUS_TOO_LARGE);
        assert_eq!(radius_limit(&refusal), Some(q("4")));
        assert!(
            refusal
                .detail
                .contains("strict affine offset-contact condition"),
            "{refusal:?}"
        );
        assert!(
            refusal.detail.contains("contact upper bound exactly 4 mm"),
            "{refusal:?}"
        );
        assert!(refusal.detail.contains(r#""step":"a""#), "{refusal:?}");
        assert!(refusal.detail.contains(r#""step":"block""#), "{refusal:?}");
        assert!(refusal.detail.contains("maximum unknown"), "{refusal:?}");
    }
}
