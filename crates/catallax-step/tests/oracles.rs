//! Implementation-independent checks for the test reader, maths and fixtures.
mod common;
use common::*;

#[test]
fn reader_respects_quotes_and_nested_references() {
    assert_eq!(
        split("'x,''#9'',y',(#21,#22),.T."),
        vec!["'x,''#9'',y'", "(#21,#22)", ".T."]
    );
    assert_eq!(
        refs("ENTITY('don''t #999',#21,(#22,#21),'#3')"),
        vec![21, 22, 21]
    );
    assert_eq!(string("a'b\\c"), "'a''b\\\\c'");
    assert_eq!(string("é"), "'\\X2\\00E9\\X0\\'");
    assert_eq!(string("😀"), "'\\X2\\D83DDE00\\X0\\'");
    assert_eq!(number("1.E+30"), r("1e30"));
    assert_eq!(number("-2.5E-7"), r("-2.5e-7"));
    assert_eq!(number("0."), n(0));
}
#[test]
fn de_boor_matches_quadratic_and_piecewise_linear_math() {
    // Quadratic Bezier P=(0,0),(1,2),(2,0) gives (2t,4t(1-t)).
    let quadratic = Spline {
        degree: 2,
        points: vec![[n(0), n(0), n(0)], [n(1), n(2), n(0)], [n(2), n(0), n(0)]],
        knots: vec![n(0), n(0), n(0), n(1), n(1), n(1)],
    };
    for i in 0..=100 {
        let t = Rational::from_ratio(i, 100).unwrap();
        assert_eq!(
            quadratic.at_span(2, &t),
            [&n(2) * &t, (&n(4) * &t) * (&n(1) - &t), n(0)]
        );
    }
    let piecewise = Spline {
        degree: 1,
        points: vec![[n(0), n(0), n(0)], [n(2), n(4), n(0)], [n(3), n(1), n(0)]],
        knots: vec![n(0), n(0), n(2), n(5), n(5)],
    };
    assert_eq!(piecewise.at_span(1, &n(1)), [n(1), n(2), n(0)]);
    assert_eq!(piecewise.at_span(2, &n(3)), [r("7/3"), n(3), n(0)]);
    let mut count = 0;
    piecewise.samples(|_| count += 1);
    assert_eq!(count, 2000);
}
#[test]
fn hand_built_fixtures_are_sound_before_export_exists() {
    for (m, b) in [void_block(), split_sphere(), touching_blocks()] {
        assert_eq!(m.check(b), vec![]);
    }
    for case in 1..=40 {
        let (m, a, b) = boolean_pair(case);
        assert_eq!(m.check(a), vec![], "case {case} first");
        assert_eq!(m.check(b), vec![], "case {case} second");
    }
    for case in 17..=26 {
        let (m, a, b) = branch_pair(case);
        assert_eq!(m.check(a), vec![]);
        assert_eq!(m.check(b), vec![]);
    }
}
#[test]
fn quaternion_oracle_and_rational_distance_boundaries() {
    assert_eq!(
        rotation_oracle([2, 1, 0, 0]),
        [
            [n(1), n(0), n(0)],
            [n(0), r("3/5"), r("-4/5")],
            [n(0), r("4/5"), r("3/5")]
        ]
    );
    radial_distance(n(25), &n(5), &n(0));
    radial_distance(n(36), &n(5), &n(1));
    radial_distance(n(16), &n(5), &n(1));
    assert!(std::panic::catch_unwind(|| radial_distance(n(37), &n(5), &n(1))).is_err());
    assert!(std::panic::catch_unwind(|| radial_distance(n(15), &n(5), &n(1))).is_err());
}

#[test]
fn reader_accepts_documented_file_and_rejects_broken_graphs() {
    let text = include_str!("fixtures/empty.step");
    let s = Step::parse(text.as_bytes(), "catallax");
    s.fixed("catallax", "1.E-7");
    assert_eq!(s.items(), vec![17]);
    let malformed = [
        text.replace("#20=DIRECTION", "#21=DIRECTION"),
        text.replace("#18,#19,#20", "#18,#19,#999"),
        text.replace(
            "ENDSEC;\nEND-ISO",
            "#21=CARTESIAN_POINT('',(0.,0.,0.));\nENDSEC;\nEND-ISO",
        ),
        text.replace("#18=CARTESIAN_POINT", "#18= CARTESIAN_POINT"),
        text.replace('\n', "\r\n"),
    ];
    for bad in malformed {
        assert!(std::panic::catch_unwind(|| Step::parse(bad.as_bytes(), "catallax")).is_err());
    }
    let out_of_walk_order = text
        .replace("(#17),#14", "(#17,#22,#21),#14")
        .replace("ENDSEC;\nEND-ISO", "#21=CARTESIAN_POINT('',(0.,0.,0.));\n#22=CARTESIAN_POINT('',(1.,0.,0.));\nENDSEC;\nEND-ISO");
    assert!(
        std::panic::catch_unwind(|| Step::parse(out_of_walk_order.as_bytes(), "catallax")).is_err()
    );
}
