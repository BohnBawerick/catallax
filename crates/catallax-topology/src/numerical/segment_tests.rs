use super::*;

#[test]
fn exact_segment_boxes_preserve_crossings_touches_overlap_and_orientation() {
    for (points, expected) in [
        ([[0, 0], [2, 0], [1, -1], [1, 1]], true),
        ([[0, 0], [2, 0], [2, 0], [3, 1]], true),
        ([[0, 0], [2, 0], [1, 0], [3, 0]], true),
        ([[0, 0], [2, 0], [3, 0], [4, 0]], false),
        ([[0, 0], [2, 0], [0, 1], [2, 1]], false),
        ([[0, 0], [4, 4], [0, 4], [1, 3]], false),
        ([[0, 0], [2, 0], [1, 0], [1, 0]], true),
        ([[0, 0], [2, 0], [3, 0], [3, 0]], false),
    ] {
        let points = points.map(|p| p.map(Rational::from_i64));
        for swap in [false, true] {
            for reverse_a in [false, true] {
                for reverse_b in [false, true] {
                    let mut a = [&points[0], &points[1]];
                    let mut b = [&points[2], &points[3]];
                    if reverse_a {
                        a.reverse();
                    }
                    if reverse_b {
                        b.reverse();
                    }
                    if swap {
                        std::mem::swap(&mut a, &mut b);
                    }
                    assert_eq!(intersects(a[0], a[1], b[0], b[1]), expected);
                }
            }
        }
    }
}

#[test]
fn sub_binary64_segment_gaps_are_not_rounded_into_contact() {
    let q = Rational::from_i64;
    let gap = Rational::parse(
        "1/100000000000000000000000000000000000000000000000000000000000000000000000000000000",
    )
    .unwrap();
    let a = [q(0), q(0)];
    let b = [q(1), q(0)];
    let c = [&q(1) + &gap, q(0)];
    let d = [q(2), q(0)];
    assert!(!intersects(&a, &b, &c, &d));
    assert!(intersects(&a, &b, &b, &d));
}
