//! When a part is constructed, a tool plane within the 1 nm modelling
//! tolerance of a target plane is merged with it once, before the boolean
//! decides anything, and the answer is marked. The plain boolean stays exact.
use catallax_boolean::{BoolOp, SUB_TOLERANCE_MERGE, boolean, boolean_merging_gaps};
use catallax_numbers::{BoundKind, Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block};

fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

// A 12 by 9 by 1 base and a stand on its back edge, `gap` inside the base's
// left, right and back walls on each side. Their bottoms coincide exactly.
fn ledge(model: &mut Model, gap: &Rational) -> (BodyId, BodyId) {
    let at = |x: &Rational, y: &Rational| {
        Pose::from_rational(
            Rotation::identity(),
            RVector3::new(x.clone(), y.clone(), Rational::zero()),
        )
    };
    let zero = Rational::zero();
    let base = block(model, [r("12"), r("9"), r("1")], &at(&zero, &zero)).unwrap();
    let width = r("12") - gap.clone() - gap.clone();
    let stand = block(model, [width, r("1"), r("6")], &at(gap, gap)).unwrap();
    (base, stand)
}

fn faces(model: &Model, body: BodyId) -> usize {
    model.counts(body).faces
}

#[test]
fn a_sub_tolerance_ledge_is_one_flush_bracket_with_its_bound() {
    for gap in ["1e-30", "3.5164312091753774e-36", "4e-7"] {
        let gap = r(gap);
        for reverse in [false, true] {
            let mut model = Model::new();
            let (base, stand) = ledge(&mut model, &gap);
            let (a, b) = if reverse {
                (stand, base)
            } else {
                (base, stand)
            };
            let answer =
                boolean_merging_gaps(&mut model, a, b, BoolOp::Unite, &mut Decider::default())
                    .unwrap();
            // The L-shaped section's six sides and two ends; no sliver face.
            assert_eq!(faces(&model, answer.body), 8, "{gap} {reverse}");
            assert!(model.check(answer.body).is_empty());
            let marks: Vec<_> = answer.marks.iter().collect();
            let [mark] = marks.as_slice() else {
                panic!("{gap} {reverse}: {marks:?}");
            };
            assert_eq!(mark.code, SUB_TOLERANCE_MERGE);
            assert_eq!(mark.kind, BoundKind::Distance);
            // The corner moved across the side wall and the back wall.
            assert_eq!(mark.bound, &gap + &gap);
            assert_eq!(model.body(answer.body).marks, answer.marks);
            // The operands are unchanged.
            assert_eq!(faces(&model, a), 6);
            assert_eq!(faces(&model, b), 6);
            assert!(model.body(b).marks.is_empty());
        }
    }
}

#[test]
fn a_gap_the_tolerance_does_not_cover_stays_exact() {
    // 6e-7 mm walls are each within 1 nm, but the corner between two of them
    // would move 1.2e-6 mm (L1), so nothing merges; nor at 1e-3 mm.
    for gap in ["6e-7", "1e-6", "1e-3"] {
        let mut model = Model::new();
        let (base, stand) = ledge(&mut model, &r(gap));
        let answer = boolean_merging_gaps(
            &mut model,
            base,
            stand,
            BoolOp::Unite,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(answer.marks.is_empty(), "{gap}: {:?}", answer.marks);
        assert_eq!(faces(&model, answer.body), 11, "{gap}");
        assert!(model.check(answer.body).is_empty());
    }
}

#[test]
fn the_plain_boolean_keeps_the_exact_sliver() {
    let mut model = Model::new();
    let (base, stand) = ledge(&mut model, &r("1e-30"));
    let answer = boolean(
        &mut model,
        base,
        stand,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty());
    assert_eq!(faces(&model, answer.body), 11);
}

#[test]
fn a_merged_pocket_opens_through_the_top() {
    // The pocket tool's top lies 1e-30 below the block's top. Merged, the
    // pocket opens through the top instead of leaving a closed skin.
    let mut model = Model::new();
    let block_body = block(&mut model, [r("10"), r("10"), r("10")], &Pose::identity()).unwrap();
    let pocket = block(
        &mut model,
        [r("4"), r("4"), &r("5") - &r("1e-30")],
        &Pose::from_rational(Rotation::identity(), RVector3::new(r("3"), r("3"), r("5"))),
    )
    .unwrap();
    let mut dec = Decider::default();
    let answer =
        boolean_merging_gaps(&mut model, block_body, pocket, BoolOp::Subtract, &mut dec).unwrap();
    // An open pocket: the block's six faces, four walls and a floor.
    assert_eq!(faces(&model, answer.body), 11);
    assert!(model.check(answer.body).is_empty());
    let marks: Vec<_> = answer.marks.iter().collect();
    let [mark] = marks.as_slice() else {
        panic!("{marks:?}");
    };
    assert_eq!(mark.code, SUB_TOLERANCE_MERGE);
    assert_eq!(mark.bound, r("1e-30"));
    // The search's work is reported with the intersect stage.
    let reported: u64 = answer.report.stages.iter().map(|s| s.stats.work).sum();
    assert_eq!(reported, dec.stats().work);
    // Exact, the pocket stays closed under a skin 1e-30 thick.
    let exact = boolean(
        &mut model,
        block_body,
        pocket,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap();
    assert!(exact.marks.is_empty());
    assert_eq!(faces(&model, exact.body), 12);
}
