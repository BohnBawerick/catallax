#![allow(missing_docs)]

use catallax_features::{Profile, extrude};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{ElementId, Loop, Point2, Pose, Segment, StepId};
use catallax_topology::Model;

fn p(x: i64, y: i64) -> Point2 {
    Point2([x.into(), y.into()])
}
fn id(s: &str) -> ElementId {
    ElementId::new(s).unwrap()
}
fn arch(reverse: bool) -> Loop {
    let (controls, end, start) = if reverse {
        (vec![p(4, 4), p(0, 4)], p(0, 0), p(4, 0))
    } else {
        (vec![p(0, 4), p(4, 4)], p(4, 0), p(0, 0))
    };
    Loop(vec![
        Segment::Bezier {
            id: id("arch"),
            controls,
            to: end,
        },
        Segment::Line {
            id: id("base"),
            to: start,
        },
    ])
}
fn validate(loops: &[Loop]) -> Result<Profile, catallax_numbers::Refusal> {
    Profile::new(&Pose::identity(), loops, &mut Decider::default())
}

#[test]
fn quadratic_and_cubic_closed_profiles_validate_in_both_orders() {
    for reverse in [false, true] {
        let mut lp = arch(reverse);
        assert!(validate(&[lp.clone()]).is_ok());
        let Segment::Bezier { controls, .. } = &mut lp.0[0] else {
            unreachable!()
        };
        *controls = vec![p(2, 4)];
        assert!(validate(&[lp]).is_ok());
    }
}

#[test]
fn nested_and_disjoint_line_loops_use_the_curved_boundary() {
    let rectangle = |x: i64, y: i64| {
        Loop(
            [(x + 1, y), (x + 1, y + 1), (x, y + 1), (x, y)]
                .into_iter()
                .enumerate()
                .map(|(i, (x, y))| Segment::Line {
                    id: id(&format!("line{i}")),
                    to: p(x, y),
                })
                .collect(),
        )
    };
    assert!(validate(&[arch(false), rectangle(1, 1), rectangle(7, 7)]).is_ok());
    assert!(validate(&[arch(true), rectangle(1, 1)]).is_ok());
}

#[test]
fn crossing_or_touching_loops_never_pass_from_sampled_chords() {
    for y in [2, 3] {
        let crossing = Loop(vec![
            Segment::Line {
                id: id("a"),
                to: p(3, y),
            },
            Segment::Line {
                id: id("b"),
                to: p(3, y + 1),
            },
            Segment::Line {
                id: id("c"),
                to: p(1, y + 1),
            },
            Segment::Line {
                id: id("d"),
                to: p(1, y),
            },
        ]);
        let error = validate(&[arch(false), crossing]).unwrap_err();
        assert!(matches!(
            error.code,
            "sketch.loop-intersection" | "sketch.bezier-proof"
        ));
    }
}

#[test]
fn invalid_or_unproved_spans_keep_named_refusals() {
    for (controls, expected) in [
        (vec![], "sketch.bezier-controls"),
        (vec![p(-2, 4), p(4, 4)], "shapes.bezier-extrusion-embedding"),
        (vec![p(0, 0), p(4, 4)], "shapes.bezier-extrusion-singular"),
    ] {
        let mut lp = arch(false);
        let Segment::Bezier {
            controls: original, ..
        } = &mut lp.0[0]
        else {
            unreachable!()
        };
        *original = controls;
        assert_eq!(validate(&[lp]).unwrap_err().code, expected);
    }
    let circle = Loop(vec![Segment::Circle {
        id: id("circle"),
        centre: p(9, 9),
        radius: 1.into(),
    }]);
    assert_eq!(
        validate(&[arch(false), circle]).unwrap_err().code,
        "sketch.bezier-curve-mix"
    );
}

#[test]
fn unfinished_extrusion_refuses_before_writing_topology() {
    let profile = validate(&[arch(false)]).unwrap();
    let mut model = Model::new();
    assert_eq!(
        extrude(
            &mut model,
            &profile,
            &StepId::new("sweep").unwrap(),
            &Rational::from_i64(2),
            None
        )
        .unwrap_err()
        .code,
        "extrude.bezier-profile"
    );
}
