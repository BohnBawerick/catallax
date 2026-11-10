use super::*;

fn check(endpoint: [i64; 2], count: usize) {
    for turn in [Turn::Ccw, Turn::Cw] {
        let endpoint = [
            endpoint[0],
            if turn == Turn::Ccw {
                endpoint[1]
            } else {
                -endpoint[1]
            },
        ];
        let start = [Rational::from_i64(5), Rational::zero()];
        let end = endpoint.map(Rational::from_i64);
        let mut dec = Decider::default();
        let (actual, c, s) = arc_step(&start, &end, turn, &mut dec).unwrap();
        assert_eq!(actual, count);
        assert_ne!(sign(&c, &mut dec).unwrap(), Sign::Negative);
        assert_eq!(
            sign(&s, &mut dec).unwrap(),
            if turn == Turn::Ccw {
                Sign::Positive
            } else {
                Sign::Negative
            }
        );
        let mut point = start.map(Real::from);
        for _ in 0..count {
            point = [
                &point[0] * &c - &point[1] * &s,
                &point[0] * &s + &point[1] * &c,
            ];
        }
        for (actual, expected) in point.iter().zip(end) {
            // Decider's separation certificate uses an annihilating
            // polynomial's proven nonzero-root bound, not a near-zero
            // tolerance. It is also exact authority for this identity.
            assert!(matches!(
                dec.sign(&(actual - &Real::from(expected))),
                Outcome::ProvenCoincident(_)
            ));
        }
        assert!(matches!(
            dec.sign(&(&c.square() + &s.square() - Real::one())),
            Outcome::ProvenCoincident(_)
        ));
    }
}

#[test]
fn one_span_uses_the_exact_directed_angle() {
    check([3, 4], 1);
}
#[test]
fn two_spans_use_equal_half_angles() {
    check([-3, 4], 2);
}
#[test]
fn three_spans_keep_the_exact_cubic_branch_selector() {
    check([-3, -4], 3);
}
#[test]
fn four_spans_keep_the_positive_quarter_angle_branch() {
    check([3, -4], 4);
}
#[test]
fn cardinal_endpoints_do_not_add_spans() {
    check([0, 5], 1);
    check([-5, 0], 2);
    check([0, -5], 3);
}
