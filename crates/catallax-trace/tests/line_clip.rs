#![allow(missing_docs)]
//! Line clips against boxes and pseudo-boxes that the SMLib gate cannot
//! reach: a direction component below SMLib's zero along a huge finite box,
//! and a left-handed pseudo-box basis. A clip that misses prunes a branch of
//! the line-surface search, so these must not report a miss where the line
//! crosses the box.

use catallax_trace::extent::{Extent1, Extent3};
use catallax_trace::pseudo_box::PseudoBox;
use catallax_trace::vector::Vector3;

const BIG: f64 = 17_592_186_044_416.0; // 2^44

/// The box [-2^44, 2^44] x [0, 1/8] x [0, 1], crossed by the line from
/// (-2^44, -1/4, 1/2) along (1, 2^-44, 0) for t in [2^42, 3 * 2^41].
fn near_parallel() -> (Extent3, Vector3, Vector3) {
    let bbox = Extent3::new(Vector3::new(-BIG, 0.0, 0.0), Vector3::new(BIG, 0.125, 1.0));
    (
        bbox,
        Vector3::new(-BIG, -0.25, 0.5),
        Vector3::new(1.0, 1.0 / BIG, 0.0),
    )
}

fn covers(span: (f64, f64), low: f64, high: f64) -> bool {
    span.0 <= low && span.1 >= high
}

#[test]
fn a_near_parallel_line_still_crosses_a_long_box() {
    let (bbox, point, vector) = near_parallel();
    let span = bbox
        .intersect_line(point, vector)
        .expect("clip")
        .expect("the line crosses the box");
    assert!(
        covers(span, 4_398_046_511_104.0, 6_597_069_766_656.0),
        "{span:?}"
    );
}

#[test]
fn a_near_parallel_line_still_crosses_a_long_pseudo_box() {
    let (bbox, point, vector) = near_parallel();
    let mut pseudo = PseudoBox::new();
    pseudo.intervals = [
        Extent1::new(bbox.min.x, bbox.max.x),
        Extent1::new(bbox.min.y, bbox.max.y),
        Extent1::new(bbox.min.z, bbox.max.z),
    ];
    let span = pseudo
        .intersect_line(point, vector)
        .expect("clip")
        .expect("the line crosses the pseudo-box");
    assert!(
        covers(span, 4_398_046_511_104.0, 6_597_069_766_656.0),
        "{span:?}"
    );
}

#[test]
fn a_left_handed_pseudo_box_is_crossed_where_the_line_meets_it() {
    let mut pseudo = PseudoBox::new();
    pseudo
        .set_basis(
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        )
        .expect("basis");
    pseudo.intervals = [Extent1::new(0.0, 1.0); 3];
    let span = pseudo
        .intersect_line(Vector3::new(-2.0, 0.5, 0.5), Vector3::new(1.0, 0.0, 0.0))
        .expect("clip")
        .expect("the line crosses the pseudo-box");
    assert_eq!(span, (1.0, 2.0));
}

#[test]
fn a_line_beside_a_box_still_misses_it() {
    let bbox = Extent3::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
    let point = Vector3::new(-1.0, 2.0, 0.5);
    let vector = Vector3::new(1.0, 1.0e-14, 0.0);
    assert_eq!(bbox.intersect_line(point, vector).expect("clip"), None);
    let mut pseudo = PseudoBox::new();
    pseudo.intervals = [Extent1::new(0.0, 1.0); 3];
    assert_eq!(pseudo.intersect_line(point, vector).expect("clip"), None);
}
