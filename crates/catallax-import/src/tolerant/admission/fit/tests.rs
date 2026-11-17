use super::*;
use catallax_numbers::with_operation_budget;

#[test]
fn pcurve_chart_and_unusable_hints_preserve_native_continuation_and_budget() {
    let source = b"ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
#1=CARTESIAN_POINT('',(0.,0.));#2=CARTESIAN_POINT('',(1.,0.));#3=CARTESIAN_POINT('',(2.,0.));
#4=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.,(2,1,2),(0.,1.,2.),.UNSPECIFIED.);
#10=CARTESIAN_POINT('',(0.,0.,0.));#11=CARTESIAN_POINT('',(0.,1.,0.));
#12=CARTESIAN_POINT('',(2.,0.,0.));#13=CARTESIAN_POINT('',(2.,1.,0.));
#14=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#10,#11),(#12,#13)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,2.),(0.,1.),.UNSPECIFIED.);
#20=(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''));
#21=DEFINITIONAL_REPRESENTATION('',(#4),#20);#22=PCURVE('',#14,#21);
#31=CARTESIAN_POINT('',(0.,99.));#32=CARTESIAN_POINT('',(1.,99.));#33=CARTESIAN_POINT('',(2.,99.));
#34=B_SPLINE_CURVE_WITH_KNOTS('',1,(#31,#32,#33),.UNSPECIFIED.,.F.,.F.,(2,1,2),(0.,1.,2.),.UNSPECIFIED.);
#35=DEFINITIONAL_REPRESENTATION('',(#34),#20);#36=PCURVE('',#14,#35);
ENDSEC;END-ISO-10303-21;";
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(source, &mut dec).unwrap();
    let scale = Rational::one();
    let surface =
        NativeSurface::from_source(&parsed, 14, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let curve =
        NativeUvCurve::from_source(&parsed, 22, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let outside =
        NativeUvCurve::from_source(&parsed, 36, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let previous = Some(surface.location([0.75, 0.0]).unwrap());
    let later = pcurve_seed(&surface, Some(&curve), 1.5, previous, &mut dec)
        .unwrap()
        .unwrap();
    let point = surface.sample(later, &mut dec).unwrap().point_mm;
    assert!(
        norm(subtract(point, [1.5, 0.0, 0.0])) < 1.0e-12,
        "later spline span must not clamp to the first span endpoint"
    );
    for (curve, parameter) in [(&outside, 1.5), (&curve, 3.0)] {
        let seed = pcurve_seed(&surface, Some(curve), parameter, previous, &mut dec).unwrap();
        assert_eq!(
            seed.unwrap().parameters.map(f64::to_bits),
            previous.unwrap().parameters.map(f64::to_bits)
        );
        let fit = surface.project([1.5, 0.0, 0.0], seed, &mut dec).unwrap();
        assert!(
            fit.residual_mm < 1.0e-6,
            "unusable optional hints must preserve a usable continuation"
        );
    }
    let (result, work) = with_operation_budget(0, || {
        pcurve_seed(&surface, Some(&curve), 1.5, previous, &mut dec)
    });
    assert!(crate::tolerant::native::control(&result.unwrap_err()));
    assert!(
        work.exhausted,
        "work exhaustion must not become an absent hint"
    );
}
