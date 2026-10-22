#![allow(missing_docs)]
mod common;
use catallax_shapes::*;
use common::*;

#[test]
fn independent_fixture_math() {
    // These checks execute without any catallax-shapes behavior beyond public data fields.
    // They validate the answers fed to the ignored tests, not an implementation.
    let q1 = [q(1, 1), q(2, 1), q(3, 1), q(4, 1)];
    let q2 = [q(2, 1), q(-1, 1), q(0, 1), q(5, 1)];
    assert_eq!(
        hamilton(&q2, &q1),
        [q(-16, 1), q(-12, 1), q(20, 1), q(10, 1)]
    );
    let m = quaternion_matrix(&q1);
    assert_eq!(
        m,
        [
            [q(-2, 3), q(2, 15), q(11, 15)],
            [q(2, 3), q(-1, 3), q(2, 3)],
            [q(1, 3), q(14, 15), q(2, 15)]
        ]
    );
    let torus = Surface::Torus(Torus {
        center: p("0", "0", "0").into(),
        axis: v("0", "0", "1").into(),
        major: r("3").into(),
        minor: r("1").into(),
        ref_dir: v("1", "0", "0").into(),
    });
    for (point, answer) in [
        (p("54/25", "72/25", "4/5"), q(0, 1)),
        (p("4", "0", "0"), q(0, 1)),
        (p("2", "0", "0"), q(0, 1)),
        (p("3", "0", "1"), q(0, 1)),
        (p("0", "0", "0"), q(64, 1)),
        (p("3", "0", "0"), q(-35, 1)),
    ] {
        assert_eq!(implicit_oracle(&torus, &point), answer);
    }
    let c = oq(&r("1.414213562373095048801688724209698"));
    assert!(sq(&c) < q(2, 1));
    assert!(sq(&(c + oq(&r("1e-33")))) > q(2, 1));
    for n in 1..=7 {
        let c = q(64 - n * n, 64 + n * n);
        let s = q(16 * n, 64 + n * n);
        assert_eq!(sq(&c) + sq(&s), q(1, 1));
        assert!(c >= q(15, 113));
        assert!(s >= q(16, 65));
    }
}
