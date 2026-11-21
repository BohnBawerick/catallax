//! Rule T5's pair search (section 4.3 of the Checkpoint 2 decisions) where
//! it has to read curved boundaries and crossings: a face bounded by
//! complete circles, an edge that passes through a face between two of its
//! samples, and two edges that cross. Every body here goes through the
//! production admission, and no threshold of the rule is changed.
#[path = "common/periodic.rs"]
mod periodic;
#[path = "common/rings.rs"]
mod rings;
use catallax_numbers::{Decider, Rational};
use catallax_topology::numerical::{
    self as n, NumericalDefinition, NumericalOptions, ToleranceWarning,
    tolerance::{self, EntityKind, FeaturePair},
};
use rings::{cap, hole, q, rim, wall};

/// The top and bottom faces of the square stock.
const BOTTOM: usize = 0;
const TOP: usize = 1;

fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
        vertex_tolerance_mm: None,
    }
}

/// The 4 by 4 by 1 mm square stock of the periodic fixture without its
/// bore: `f0` is the bottom at z = 0 and `f1` the top at z = 1.
fn block() -> NumericalDefinition {
    let mut d = periodic::bore();
    d.name = "authored-block".into();
    d.faces.pop();
    d.faces[BOTTOM].loops.truncate(1);
    d.faces[TOP].loops.truncate(1);
    d.edges.truncate(12);
    d.vertices.truncate(8);
    d
}

/// A disc of radius 2 and height 1 with a blind bore of radius 1 from below
/// whose ceiling is at `ceiling`. Every planar face is bounded by complete
/// circles only. The web between the ceiling and the top is `1 - ceiling`
/// thick, and the ceiling's rim lies under the middle of the top face, a
/// whole millimetre from the top face's own rim.
fn drilled_disc(ceiling: &str) -> NumericalDefinition {
    let mut d = rings::empty("authored-drilled-disc");
    let outer_0 = rim(&mut d, "outer-0", ["0", "0", "0"], "2");
    let outer_1 = rim(&mut d, "outer-1", ["0", "0", "1"], "2");
    let bore_0 = rim(&mut d, "bore-0", ["0", "0", "0"], "1");
    let bore_1 = rim(&mut d, "bore-1", ["0", "0", ceiling], "1");
    cap(&mut d, "bottom", false, &outer_0, &[&bore_0]);
    cap(&mut d, "top", true, &outer_1, &[]);
    cap(&mut d, "ceiling", false, &bore_1, &[]);
    wall(&mut d, "outer-wall", &outer_0, &outer_1, true);
    wall(&mut d, "bore-wall", &bore_0, &bore_1, false);
    d
}

/// The square stock with a pocket of radius 1 and depth 0.5 in its top, and
/// in the pocket a pin of radius `pin_radius` whose top is at `pin_top`. The
/// pin's top rim lies inside the hole of the stock's top face.
fn pinned_pocket(pin_radius: &str, pin_top: &str) -> NumericalDefinition {
    let mut d = block();
    d.name = "authored-pinned-pocket".into();
    let pocket_1 = rim(&mut d, "pocket-1", ["0", "0", "1"], "1");
    let pocket_0 = rim(&mut d, "pocket-0", ["0", "0", "0.5"], "1");
    let pin_0 = rim(&mut d, "pin-0", ["0", "0", "0.5"], pin_radius);
    let pin_1 = rim(&mut d, "pin-1", ["0", "0", pin_top], pin_radius);
    hole(&mut d, TOP, true, &pocket_1);
    wall(&mut d, "pocket-wall", &pocket_0, &pocket_1, false);
    cap(&mut d, "pocket-floor", true, &pocket_0, &[&pin_0]);
    wall(&mut d, "pin-wall", &pin_0, &pin_1, true);
    cap(&mut d, "pin-top", true, &pin_1, &[]);
    d
}

/// The square stock with a bore of radius 0.5 from below whose ceiling is at
/// `ceiling`. Above 1 the bore breaks through the top face, which its
/// topology does not know: the seam of the bore's wall passes through the
/// top face at (0.5, 0, 1).
fn bored_block(ceiling: &str) -> NumericalDefinition {
    let mut d = block();
    d.name = "authored-bored-block".into();
    let bore_0 = rim(&mut d, "bore-0", ["0", "0", "0"], "0.5");
    let bore_1 = rim(&mut d, "bore-1", ["0", "0", ceiling], "0.5");
    hole(&mut d, BOTTOM, false, &bore_0);
    wall(&mut d, "bore-wall", &bore_0, &bore_1, false);
    cap(&mut d, "ceiling", false, &bore_1, &[]);
    d
}

/// The square stock with two blind bores of radius 0.5 that meet at half
/// height: one from below about `lower`, one from above about `upper`, both
/// given as x and y. The rim of the first bore's ceiling and the rim of the
/// second bore's floor are two circles in the plane z = 0.5.
fn meeting_bores(lower: [&str; 2], upper: [&str; 2]) -> NumericalDefinition {
    let mut d = block();
    d.name = "authored-meeting-bores".into();
    let lower_0 = rim(&mut d, "lower-0", [lower[0], lower[1], "0"], "0.5");
    let lower_1 = rim(&mut d, "lower-1", [lower[0], lower[1], "0.5"], "0.5");
    hole(&mut d, BOTTOM, false, &lower_0);
    wall(&mut d, "lower-wall", &lower_0, &lower_1, false);
    cap(&mut d, "lower-ceiling", false, &lower_1, &[]);
    let upper_0 = rim(&mut d, "upper-0", [upper[0], upper[1], "0.5"], "0.5");
    let upper_1 = rim(&mut d, "upper-1", [upper[0], upper[1], "1"], "0.5");
    hole(&mut d, TOP, true, &upper_1);
    wall(&mut d, "upper-wall", &upper_0, &upper_1, false);
    cap(&mut d, "upper-floor", true, &upper_0, &[]);
    d
}

fn build(definition: NumericalDefinition) -> Result<n::NumericalBody, catallax_numbers::Refusal> {
    n::build(definition, options(), &mut Decider::default())
}

fn pair_warnings(body: &n::NumericalBody) -> Vec<FeaturePair> {
    body.tolerances()
        .warnings
        .iter()
        .filter_map(|warning| match warning {
            ToleranceWarning::Pair(pair) => Some(pair.clone()),
            _ => None,
        })
        .collect()
}

fn names(pair: &FeaturePair, a: (EntityKind, &str), b: (EntityKind, &str)) -> bool {
    let is = |entity: &(EntityKind, String), wanted: (EntityKind, &str)| {
        entity.0 == wanted.0 && entity.1 == wanted.1
    };
    (is(&pair.a, a) && is(&pair.b, b)) || (is(&pair.a, b) && is(&pair.b, a))
}

/// Whether a pair joins the stock's top face with a part of the pin.
fn joins_the_pin_and_the_top(pair: &FeaturePair) -> bool {
    let top = |entity: &(EntityKind, String)| entity.1 == "f1";
    let pin = |entity: &(EntityKind, String)| entity.1.starts_with("pin-");
    (top(&pair.a) && pin(&pair.b)) || (top(&pair.b) && pin(&pair.a))
}

fn near(value: &Rational, expected: f64) -> bool {
    (value.to_f64_nearest() - expected).abs() <= 1e-12
}

/// A face whose only loop is one complete circle has an interior. The
/// bore's ceiling rim sits under the middle of the disc's top face. Its
/// edge, its vertex and the seam that ends at it each store rule T1's floor,
/// 1e-7 mm, and the top face stores none, so the rule's numbers are: an
/// error under a web of 2e-7 mm, a warning under 5e-7 mm, nothing beyond.
#[test]
fn a_face_bounded_by_one_complete_circle_has_an_interior() {
    let top = (EntityKind::Face, "top");
    let under_it = [
        (EntityKind::Edge, "bore-1"),
        (EntityKind::Vertex, "bore-1-v"),
        (EntityKind::Edge, "bore-wall-seam"),
    ];
    // A web of 4e-7 mm: 1e-7 > 4e-7 / 5, and 1e-7 <= 4e-7 / 2.
    let warned = build(drilled_disc("0.9999996")).expect("a warning does not refuse");
    let pairs = pair_warnings(&warned);
    for entity in under_it {
        let pair = pairs
            .iter()
            .find(|pair| names(pair, top, entity))
            .unwrap_or_else(|| panic!("{entity:?} and the top face are a pair: {pairs:?}"));
        assert!(near(&pair.distance_mm, 4e-7), "{pair:?}");
        assert!(near(&pair.tolerance_sum_mm, 1e-7), "{pair:?}");
    }
    assert_eq!(pairs.len(), under_it.len(), "{pairs:?}");
    // A web of 1.5e-7 mm: 1e-7 > 1.5e-7 / 2.
    let refused = build(drilled_disc("0.99999985")).unwrap_err();
    assert_eq!(refused.code, "numerical.tolerance-pair");
    assert!(
        refused.detail.contains("face \"top\""),
        "{}",
        refused.detail
    );
    // A web of 6e-7 mm: 1e-7 <= 6e-7 / 5.
    let clear = build(drilled_disc("0.9999994")).expect("beyond five times the sum");
    assert!(
        pair_warnings(&clear).is_empty(),
        "{:?}",
        pair_warnings(&clear)
    );
    // The small features of rules T2 and T8 read the same faces: a web of
    // 0.01 mm shows as the same three pairs within 0.05 mm.
    let mut dec = Decider::default();
    let body = build(drilled_disc("0.99")).unwrap();
    let features = tolerance::features_within(&body, &q("0.05"), &mut dec).unwrap();
    for entity in under_it {
        let pair = features
            .iter()
            .find(|pair| names(pair, top, entity))
            .unwrap_or_else(|| panic!("{entity:?} is 0.01 mm under the top face: {features:?}"));
        assert!(near(&pair.distance_mm, 0.01), "{pair:?}");
    }
    assert_eq!(features.len(), under_it.len(), "{features:?}");
}

/// A hole that is one complete circle is not part of its face. The pin's
/// top is flush with the stock's top face and lies inside that face's hole,
/// half a millimetre from the hole's rim. Nothing of the pin is near the
/// face, so the body is admitted with no pair, flush or a hair under.
#[test]
fn a_hole_of_one_complete_circle_is_not_part_of_its_face() {
    for pin_top in ["1", "0.9999996"] {
        let body = build(pinned_pocket("0.5", pin_top))
            .unwrap_or_else(|refusal| panic!("pin top at {pin_top}: {refusal:?}"));
        assert!(
            pair_warnings(&body).is_empty(),
            "pin top at {pin_top}: {:?}",
            pair_warnings(&body)
        );
    }
    // Within 0.75 mm the pin's top rim is 0.5 mm from the hole's rim, and
    // no part of the pin is a pair with the face the hole is in.
    let mut dec = Decider::default();
    let body = build(pinned_pocket("0.5", "1")).unwrap();
    let features = tolerance::features_within(&body, &q("0.75"), &mut dec).unwrap();
    assert!(
        !features.iter().any(joins_the_pin_and_the_top),
        "{features:?}"
    );
    let rims = features
        .iter()
        .find(|pair| {
            names(
                pair,
                (EntityKind::Edge, "pin-1"),
                (EntityKind::Edge, "pocket-1"),
            )
        })
        .expect("the two rims are 0.5 mm apart");
    assert!(near(&rims.distance_mm, 0.5), "{rims:?}");
}

/// A hole is read on its circle, not on chords of it. The pin is 0.99 mm in
/// radius, so the groove round it is 0.01 mm wide, and its top is flush
/// with the stock's top face. The pin's top rim runs 0.01 mm inside the
/// hole's rim all the way round. A hole read as a polygon of fewer than 23
/// sides would have the pin's rim outside it in places, and so on the face,
/// at no distance from it. The body is admitted at a cap of 0.001 mm (a
/// complete circle keeps its neighbours a cap away), with no pair.
#[test]
fn a_narrow_groove_round_a_flush_pin_is_no_pair_with_the_face_about_it() {
    let narrow = NumericalOptions {
        tolerance_mm: q("0.001"),
        vertex_tolerance_mm: None,
    };
    let mut dec = Decider::default();
    let body = n::build(pinned_pocket("0.99", "1"), narrow, &mut dec)
        .unwrap_or_else(|refusal| panic!("{refusal:?}"));
    assert!(
        pair_warnings(&body).is_empty(),
        "{:?}",
        pair_warnings(&body)
    );
    let features = tolerance::features_within(&body, &q("0.05"), &mut dec).unwrap();
    assert!(
        !features.iter().any(joins_the_pin_and_the_top),
        "{features:?}"
    );
    let rims = features
        .iter()
        .find(|pair| {
            names(
                pair,
                (EntityKind::Edge, "pin-1"),
                (EntityKind::Edge, "pocket-1"),
            )
        })
        .expect("the two rims are 0.01 mm apart");
    assert!(near(&rims.distance_mm, 0.01), "{rims:?}");
}

/// An edge that passes through a face between two of its samples is at no
/// distance from it. The bore breaks through the top of the stock, and the
/// seam of its wall crosses the top face at (0.5, 0, 1), two thirds of the
/// way along. The seam's samples either side of the crossing are 0.06 mm
/// from the face, a hundred thousand times the reach of the rule, 5e-7 mm.
/// A bore that stops 0.25 mm short of the top is no pair.
#[test]
fn an_edge_that_crosses_a_face_between_two_samples_is_refused() {
    let refused = build(bored_block("1.5")).unwrap_err();
    assert_eq!(refused.code, "numerical.tolerance-pair");
    assert!(
        refused.detail.contains("edge \"bore-wall-seam\"")
            && refused.detail.contains("face \"f1\""),
        "{}",
        refused.detail
    );
    let short = build(bored_block("0.75")).expect("the bore stops short of the top");
    assert!(
        pair_warnings(&short).is_empty(),
        "{:?}",
        pair_warnings(&short)
    );
}

/// Two edges that cross are at no distance from each other. The two bores
/// meet at half height, and the rims of the lower one's ceiling and the
/// upper one's floor cross.
///
/// About (0, 0.49) and (0, -0.49) they cross at (0.0995, 0) and
/// (-0.0995, 0), at an angle of about 23 degrees, far from both rims'
/// vertices.
///
/// About (0, 0) and (0.5, -0.502) they cross at right angles at
/// (0.499996, -0.002), which on the lower rim is 0.004 radians before its
/// vertex at (0.5, 0): between the rim's last sample but one and its end,
/// where the nearest sample is the vertex itself.
///
/// In both, the first pair the rule refuses is the two rims: no vertex is
/// within the rule's reach of anything it does not touch.
#[test]
fn two_edges_that_cross_are_refused() {
    for (lower, upper) in [
        (["0", "0.49"], ["0", "-0.49"]),
        (["0", "0"], ["0.5", "-0.502"]),
    ] {
        let refused = build(meeting_bores(lower, upper)).unwrap_err();
        assert_eq!(refused.code, "numerical.tolerance-pair", "{lower:?}");
        assert!(
            refused.detail.contains("edge \"lower-1\"")
                && refused.detail.contains("edge \"upper-0\""),
            "{lower:?}: {}",
            refused.detail
        );
    }
}
