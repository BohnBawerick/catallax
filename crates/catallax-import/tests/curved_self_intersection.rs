//! Self-intersection of solids with cylinders, spheres and tori.

mod common;

use std::collections::BTreeMap;

use catallax_boolean::{BoolOp, boolean};
use catallax_import::{ImportRefusal, check_self_intersection, read_step};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::Pose;
use catallax_step::{StepOptions, write_step};
use catallax_topology::{Model, block, cylinder};
use common::*;

fn check(text: &str) -> Result<(), ImportRefusal> {
    check_self_intersection(&read_step(text.as_bytes()).unwrap())
}

fn crossing(r: &Result<(), ImportRefusal>) -> bool {
    matches!(r, Err(ImportRefusal::SelfIntersecting { .. }))
}

fn centre(x: &str, y: &str, z: &str) -> [Rational; 3] {
    [r(x), r(y), r(z)]
}

#[test]
fn a_cylinder_piercing_a_box_wall_intersects() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    // Along x from x = 3 to x = 6, through the wall at x = 4.
    s.cylinder(["3", "2", "2"], 0, "3", "0.5", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_cylinder_crossing_a_wall_edge_on_intersects() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    // Along z, centred on the wall plane x = 4: the wall cuts it along two lines.
    s.cylinder(["4", "2", "1"], 2, "2", "0.5", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn two_overlapping_cylinders_intersect() {
    let mut s = Soup::new();
    s.cylinder(["0", "0", "0"], 2, "2", "1", Caps::Both);
    s.cylinder(["0.5", "0", "0"], 2, "2", "1", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_sphere_poking_through_a_plane_face_intersects() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    s.sphere(["4", "2", "2"], "1");
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_torus_clipping_a_box_intersects() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    s.torus(["4", "2", "2"], 0, "1", "0.4");
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_box_with_a_through_hole_is_clear() {
    let mut holes: BTreeMap<usize, Vec<_>> = BTreeMap::new();
    holes.insert(4, vec![(centre("2", "2", "0"), 2, 0, r("1"), true)]);
    holes.insert(5, vec![(centre("2", "2", "1"), 2, 0, r("1"), false)]);
    let mut s = Soup::new();
    s.box_with_holes(["0", "0", "0"], ["4", "4", "1"], &holes);
    s.hole_wall(["2", "2", "0"], 2, "1", "1");
    let result = check(&s.finish());
    assert!(result.is_ok(), "{result:?}");
}

#[test]
fn a_boss_on_a_box_sharing_a_circular_edge_is_clear() {
    let mut holes: BTreeMap<usize, Vec<_>> = BTreeMap::new();
    holes.insert(5, vec![(centre("2", "2", "1"), 2, 0, r("1"), false)]);
    let mut s = Soup::new();
    s.box_with_holes(["0", "0", "0"], ["4", "4", "1"], &holes);
    s.cylinder(["2", "2", "1"], 2, "2", "1", Caps::TopOnly);
    let result = check(&s.finish());
    assert!(result.is_ok(), "{result:?}");
}

#[test]
fn a_boss_sunk_into_the_box_intersects() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "1"]);
    // Starts inside the box and passes through its top face.
    s.cylinder(["2", "2", "0.5"], 2, "2", "1", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_cylinder_standing_on_a_whole_face_has_its_edge_inside_that_face() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    // The wall's lower ring lies in the interior of the box top, which has no hole
    // for it. That is contact along a curve inside one face, as for planar faces.
    s.cylinder(["2", "2", "4"], 2, "1", "0.5", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_sphere_tangent_to_a_face_at_one_point_is_allowed() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    s.sphere(["2", "2", "5"], "1");
    let result = check(&s.finish());
    assert!(result.is_ok(), "{result:?}");
}

#[test]
fn separate_solids_are_not_compared_even_when_curved() {
    let mut a = Soup::new();
    a.cylinder(["0", "0", "0"], 2, "2", "1", Caps::Both);
    let mut b = Soup::new();
    b.cylinder(["0.5", "0", "0"], 2, "2", "1", Caps::Both);
    // Two files' solids joined by hand: the faces of each are checked alone.
    assert!(check(&a.finish()).is_ok());
    assert!(check(&b.finish()).is_ok());
}

#[test]
fn an_exporter_round_trip_of_a_bracket_is_clear() {
    // A plate with a boss and a through hole, made by Catallax's own booleans.
    let mut m = Model::new();
    let plate = block(
        &mut m,
        [20, 10, 4].map(Rational::from_i64),
        &Pose::identity(),
    )
    .unwrap();
    let boss = cylinder(
        &mut m,
        Rational::from_i64(3),
        Rational::from_i64(8),
        &Pose::identity(),
    )
    .unwrap();
    let mut dec = Decider::default();
    let body = boolean(&mut m, plate, boss, BoolOp::Unite, &mut dec)
        .expect("union")
        .body;
    let mut bytes = Vec::new();
    write_step(&m, &[body], &StepOptions::default(), &mut dec, &mut bytes).expect("export");
    let result = check(std::str::from_utf8(&bytes).unwrap());
    assert!(result.is_ok(), "{result:?}");
}

#[test]
fn an_elliptical_face_inside_a_box_is_clear_and_one_through_a_wall_intersects() {
    let mut inside = Soup::new();
    inside.box_(["0", "0", "0"], ["4", "4", "4"]);
    inside.elliptical_disc(["2", "2", "2"], 2, "1", "0.5");
    let result = check(&inside.finish());
    assert!(result.is_ok(), "{result:?}");
    let mut through = Soup::new();
    through.box_(["0", "0", "0"], ["4", "4", "4"]);
    through.elliptical_disc(["4", "2", "2"], 2, "1", "0.5");
    let result = check(&through.finish());
    assert!(crossing(&result), "{result:?}");
}

#[test]
fn a_pair_that_cannot_be_decided_is_not_verified_and_never_a_pass() {
    let mut s = Soup::new();
    s.box_(["0", "0", "0"], ["4", "4", "4"]);
    s.unbounded_plane(["2", "2", "2"], 0);
    let result = check(&s.finish());
    let Err(ImportRefusal::SelfIntersectionNotVerified {
        solid,
        faces,
        reason,
    }) = &result
    else {
        panic!("{result:?}");
    };
    assert!(solid.starts_with('#') && faces[0].starts_with('#') && faces[1].starts_with('#'));
    assert!(reason.contains("unbounded"), "{reason}");
    let json = serde_json::to_string(result.as_ref().unwrap_err()).unwrap();
    assert!(json.starts_with(r#"{"code":"self_intersection_not_verified""#));
}

#[test]
fn an_intersection_found_elsewhere_wins_over_an_undecided_pair() {
    let mut s = Soup::new();
    s.unbounded_plane(["10", "10", "10"], 0);
    s.cylinder(["0", "0", "0"], 2, "2", "1", Caps::Both);
    s.cylinder(["0.5", "0", "0"], 2, "2", "1", Caps::Both);
    let result = check(&s.finish());
    assert!(crossing(&result), "{result:?}");
}
