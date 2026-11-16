//! Numerical source fits checked against independent polynomial identities.
use catallax_import::geometry_source::{AngleUnit, FitOptions, ParsedStep};
use catallax_numbers::{Decider, Rational};
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn source() -> ParsedStep {
    ParsedStep::read(
        b"DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=AXIS2_PLACEMENT_3D('',#1,$,$);
#3=PARABOLA('',#2,1.);
#4=DIRECTION('',(0.,0.,1.));
#5=VECTOR('',#4,7.);
#6=SURFACE_OF_LINEAR_EXTRUSION('',#3,#5);
#7=CIRCLE('',#2,1.);
#8=CARTESIAN_POINT('',(0.,0.));
#9=AXIS2_PLACEMENT_2D('',#8,$);
#10=CIRCLE('',#9,1.);
#11=DEFINITIONAL_REPRESENTATION('',(#10),#90);
#12=PCURVE('',#6,#11);
#13=PLANE('',#2);
#14=OFFSET_SURFACE('',#13,2.,.F.);
#15=OFFSET_CURVE_3D('',#7,0.25,.F.,#4);
#16=CARTESIAN_POINT('',(0.5,1.,0.));
#17=CARTESIAN_POINT('',(1.,0.,0.));
#18=B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#16,#17),.UNSPECIFIED.,.F.,.F.,(3,3),(0.1,0.9),.UNSPECIFIED.);
#90=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('param',''));
ENDSEC;",
        &mut Decider::default(),
    )
    .unwrap()
}
#[test]
fn secant_fits_bound_the_complete_parabola_interval() {
    let s = source();
    let mut dec = Decider::default();
    let charts = s
        .curve_charts(3, &q("1"), AngleUnit::Radians, &mut dec)
        .unwrap();
    let options = FitOptions {
        tolerance_mm: q("0.01"),
        maximum_pieces: 64,
    };
    let pieces = charts[0]
        .fit_curve([q("-1"), q("1")], &options, &mut dec)
        .unwrap();
    assert!(pieces.len() > 1);
    assert_eq!(pieces.first().unwrap().interval()[0], q("-1"));
    assert_eq!(pieces.last().unwrap().interval()[1], q("1"));
    for pair in pieces.windows(2) {
        assert_eq!(pair[0].interval()[1], pair[1].interval()[0]);
    }
    for p in pieces {
        assert_eq!(p.source(), 3);
        let [a, b] = p.interval();
        let middle = (a + b) / q("2");
        let truth = (&middle * &middle) / q("4");
        let controls = p.geometry().controls();
        let secant = (&controls[0][0] + &controls[1][0]) / q("2");
        let actual = (&secant - &truth).abs();
        assert!(actual > q("0"));
        assert!(actual <= *p.maximum_deviation_mm());
        assert!(*p.maximum_deviation_mm() <= options.tolerance_mm);
    }
}
#[test]
fn bilinear_fits_bound_the_complete_swept_rectangle() {
    let s = source();
    let mut dec = Decider::default();
    let charts = s
        .surface_charts(6, &q("1"), AngleUnit::Radians, &mut dec)
        .unwrap();
    let options = FitOptions {
        tolerance_mm: q("0.01"),
        maximum_pieces: 64,
    };
    let pieces = charts[0]
        .fit_surface([[q("-1"), q("1")], [q("0"), q("2")]], &options, &mut dec)
        .unwrap();
    assert!(pieces.len() > 1);
    for p in pieces {
        assert_eq!(p.source(), 6);
        let [a, b] = &p.rectangle()[0];
        let middle = (a + b) / q("2");
        let truth = (&middle * &middle) / q("4");
        let controls = p.geometry().controls();
        let secant = (&controls[0][0] + &controls[2][0]) / q("2");
        assert!((&secant - &truth).abs() <= *p.maximum_deviation_mm());
        assert!(*p.maximum_deviation_mm() <= options.tolerance_mm);
        assert_eq!(controls[0][2], q("0"));
        assert_eq!(controls[1][2], q("2"));
    }
}
#[test]
fn a_fit_never_returns_partial_geometry_after_exhausting_its_limit() {
    let s = source();
    let mut dec = Decider::default();
    let charts = s
        .curve_charts(7, &q("1"), AngleUnit::Degrees, &mut dec)
        .unwrap();
    let options = FitOptions {
        tolerance_mm: q("0.0001"),
        maximum_pieces: 1,
    };
    assert!(
        charts[0]
            .fit_curve([q("0"), q("360")], &options, &mut dec)
            .is_err()
    );
}

#[test]
fn raw_pcurve_units_cannot_receive_a_millimetre_fit_bound() {
    let s = source();
    let mut dec = Decider::default();
    let charts = s
        .curve_charts(12, &q("1"), AngleUnit::Radians, &mut dec)
        .unwrap();
    assert_eq!(
        charts[0].coordinate_space(),
        catallax_import::geometry_source::CoordinateSpace::SurfaceParameters
    );
    let options = FitOptions {
        tolerance_mm: q("0.01"),
        maximum_pieces: 64,
    };
    assert!(
        charts[0]
            .fit_curve([q("0"), q("1")], &options, &mut dec)
            .is_err()
    );
}

#[test]
fn finite_parabolas_and_simple_offsets_have_exact_native_carriers() {
    use catallax_numbers::{Outcome, Real};
    let s = source();
    let mut dec = Decider::default();
    let curve = s
        .parabola_segment(3, [q("-1"), q("1")], &q("1"), &mut dec)
        .unwrap();
    let point = curve.point_at(&Real::zero(), &mut dec).unwrap();
    for value in [&point.x, &point.y, &point.z] {
        assert!(matches!(dec.sign(value), Outcome::ProvenCoincident(_)));
    }
    let catallax_shapes::Surface::Plane(plane) = s.native_surface(14, &q("1"), &mut dec).unwrap()
    else {
        panic!("plane offset")
    };
    assert!(matches!(
        dec.sign(&(&plane.origin.z - &Real::from_i64(2))),
        Outcome::ProvenCoincident(_)
    ));
    let catallax_shapes::Curve::Circle(circle) = s.native_curve(15, &q("1"), &mut dec).unwrap()
    else {
        panic!("circle offset")
    };
    assert!(matches!(
        dec.sign(&(&circle.radius_sq - &Real::from_rational(q("0.5625")))),
        Outcome::ProvenCoincident(_)
    ));
}

#[test]
fn finite_spline_endpoints_keep_their_exact_decimal_domain() {
    let s = source();
    let mut dec = Decider::default();
    let charts = s
        .curve_charts(18, &q("1"), AngleUnit::Radians, &mut dec)
        .unwrap();
    let options = FitOptions {
        tolerance_mm: q("0.01"),
        maximum_pieces: 64,
    };
    let pieces = charts[0]
        .fit_curve([q("0.1"), q("0.9")], &options, &mut dec)
        .unwrap();
    assert_eq!(pieces.first().unwrap().interval()[0], q("0.1"));
    assert_eq!(pieces.last().unwrap().interval()[1], q("0.9"));
}
