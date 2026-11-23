const LOFT: &str = include_str!("../../../catallax-import/tests/fixtures/m27-loft-occt.step");

pub fn trimmed_loft() -> String {
    // Extend the first ruled support to v = -1. Its physical edges and
    // pcurves still bound the original [0, 1] rows.
    assert!(LOFT.contains("(0.,2.2360679775),(0.,1.)"));
    assert!(LOFT.contains("#33 = CARTESIAN_POINT('',(1.,-1.,0.));"));
    assert!(LOFT.contains("#35 = CARTESIAN_POINT('',(1.,1.,0.));"));
    LOFT.replacen("(0.,2.2360679775),(0.,1.)", "(0.,2.2360679775),(-1.,1.)", 1)
        .replace(
            "#33 = CARTESIAN_POINT('',(1.,-1.,0.));",
            "#33 = CARTESIAN_POINT('',(1.,-2.,-2.));",
        )
        .replace(
            "#35 = CARTESIAN_POINT('',(1.,1.,0.));",
            "#35 = CARTESIAN_POINT('',(2.,0.,-2.));",
        )
}
