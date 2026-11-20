//! Minimum distance between two whole solids: face pairs far from the minimum
//! need no closed form, a line parallel to a circle's plane has one, turned
//! parallel faces answer from their edges, and a pair no closed form covers
//! widens the stated gap with a mark.
use catallax_measure::{Accuracy, solid_gap};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Isometry, Pose, RPoint3, RVector3, Rotation, Vector3};
use catallax_topology::{BodyId, Model, block, cone, cylinder, transformed_body};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn at(x: &str, y: &str, z: &str) -> Pose {
    Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)))
}
/// A turn by `degrees` about the axis through the origin, then a shift.
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
fn gap(
    model: &Model,
    a: BodyId,
    b: BodyId,
) -> catallax_numbers::Answer<catallax_measure::Enclosure> {
    solid_gap(model, a, b, &Accuracy::default(), &mut Decider::default()).unwrap()
}
fn within_accuracy(e: &catallax_measure::Enclosure) -> bool {
    e.width() <= &Accuracy::default().relative * &e.lo
}

#[test]
fn an_edge_parallel_to_a_rim_plane_meets_the_rim_in_closed_form() {
    // The block's lower front edge runs along x at y = 3, z = 2. The nearest
    // cylinder point is on its top rim at (0, 1, 1): the gap is sqrt(5).
    let mut model = Model::new();
    let pin = cylinder(&mut model, q("1"), q("1"), &at("0", "0", "0")).unwrap();
    let bar = block(&mut model, [q("10"), q("1"), q("1")], &at("-5", "3", "2")).unwrap();
    let answer = gap(&model, pin, bar);
    let e = answer.value;
    assert!(&e.lo * &e.lo <= q("5") && &e.hi * &e.hi >= q("5"), "{e:?}");
    assert!(within_accuracy(&e) && answer.marks.is_empty());
}

#[test]
fn faces_far_from_the_minimum_need_no_closed_form() {
    // Cone rims against block edges have no closed form, but every such pair
    // is farther than the base plane's 188 mm and is never visited.
    let mut model = Model::new();
    let base = block(&mut model, [q("10"), q("10"), q("10")], &at("0", "0", "0")).unwrap();
    let tip = cone(&mut model, q("2"), q("1"), q("4"), &at("200", "5", "3")).unwrap();
    let answer = gap(&model, base, tip);
    assert_eq!(answer.value.lo, q("188"));
    assert_eq!(answer.value.hi, q("188"));
    assert!(answer.marks.is_empty());
}

#[test]
fn turned_parallel_faces_answer_from_their_edges() {
    // A 4 mm cube turned 30 degrees about z: its top and bottom faces stay
    // parallel to the base block's, on planes with no rational frame. Its
    // nearest corner is 12 - sqrt(3) from the origin in x: the gap is 2 - sqrt(3).
    let mut model = Model::new();
    let base = block(&mut model, [q("10"), q("10"), q("10")], &at("0", "0", "0")).unwrap();
    let cube = block(&mut model, [q("4"), q("4"), q("4")], &at("-2", "-2", "-2")).unwrap();
    let cube = turned(&mut model, cube, ["0", "0", "1"], "30", ["13", "5", "5"]);
    let answer = gap(&model, base, cube);
    let e = answer.value;
    let (a, b) = (&q("2") - &e.lo, &q("2") - &e.hi);
    assert!(
        b > q("0") && &a * &a >= q("3") && &b * &b <= q("3"),
        "{e:?}"
    );
    assert!(within_accuracy(&e) && answer.marks.is_empty());
}

#[test]
fn a_pair_without_a_closed_form_states_a_wider_gap_with_a_mark() {
    // The second cylinder lies along y. Its near rim faces the first
    // cylinder's side at sqrt(26) - 2, a line against a circle in a parallel
    // plane. The two cylinders' rims meet at non-parallel normals, which no
    // closed form covers, so the stated gap reaches down to the 3 mm between
    // the rims' boxes.
    let mut model = Model::new();
    let a = cylinder(&mut model, q("2"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = cylinder(&mut model, q("2"), q("4"), &at("0", "0", "-2")).unwrap();
    let b = turned(&mut model, b, ["1", "0", "0"], "90", ["7", "3", "1"]);
    let answer = gap(&model, a, b);
    let e = answer.value;
    let top = &e.hi + &q("2");
    assert!(
        e.lo >= q("3") && e.lo < e.hi && &top * &top >= q("26"),
        "{e:?}"
    );
    assert!(
        answer
            .marks
            .iter()
            .any(|m| m.code == "measure.accuracy-not-reached"),
        "{:?}",
        answer.marks
    );
}
