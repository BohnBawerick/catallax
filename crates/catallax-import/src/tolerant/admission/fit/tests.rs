use super::*;
use catallax_numbers::with_operation_budget;

const SOURCE: &[u8] = b"ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
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
#41=OFFSET_CURVE_2D('',#4,0.25,.F.);
#42=DEFINITIONAL_REPRESENTATION('',(#41),#20);#43=PCURVE('',#14,#42);
ENDSEC;END-ISO-10303-21;";

#[test]
fn pcurve_chart_and_unusable_hints_preserve_native_continuation_and_budget() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SOURCE, &mut dec).unwrap();
    let scale = Rational::one();
    let surface =
        NativeSurface::from_source(&parsed, 14, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let curve = FaceCopy::File(Arc::new(
        NativeUvCurve::from_source(&parsed, 22, &scale, AngleUnit::Radians, &mut dec).unwrap(),
    ));
    let outside = FaceCopy::File(Arc::new(
        NativeUvCurve::from_source(&parsed, 36, &scale, AngleUnit::Radians, &mut dec).unwrap(),
    ));
    let continuation = surface.location([0.75, 0.0]).unwrap();
    let previous = Some(continuation);
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
            continuation.parameters.map(f64::to_bits)
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

#[test]
fn optional_pcurve_compilation_preserves_fitting_and_control_refusals() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SOURCE, &mut dec).unwrap();
    let scale = Rational::one();
    let mut context = Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), 0.01);
    let surface = context.surface(14, &mut dec).unwrap();
    let continuation = surface.location([0.75, 0.0]).unwrap();
    let hint = context.pcurve(43, &mut dec).unwrap().map(FaceCopy::File);
    let seed = pcurve_seed(&surface, hint.as_ref(), 1.5, Some(continuation), &mut dec).unwrap();
    assert_eq!(
        seed.unwrap().parameters.map(f64::to_bits),
        continuation.parameters.map(f64::to_bits)
    );
    let fit = surface.project([1.5, 0.0, 0.0], seed, &mut dec).unwrap();
    assert!(
        fit.residual_mm < 1.0e-6,
        "an unsupported optional 2D-offset pcurve must not prevent native 3D fitting"
    );
    let (result, work) = with_operation_budget(0, || context.pcurve(22, &mut dec));
    assert!(crate::tolerant::native::control(&result.unwrap_err()));
    assert!(
        work.exhausted,
        "compilation work refusal must not become an absent hint"
    );
}

/// A straight edge along `y = 0` of a flat spline patch, with three u/v
/// curves for it. `#22` traces the edge at the edge's own speed. `#26` traces
/// the same segment as `u = t/2 + t^2/4`: it is a quarter of a millimetre off
/// at the same parameter, at `t = 1`, and nowhere off as a curve. `#33` stops
/// at `t = 1.5`, short of the edge's end, and `#37` runs along `v = 99`, off
/// the patch.
const SPEEDS: &[u8] = b"ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
#1=CARTESIAN_POINT('',(0.,0.));#2=CARTESIAN_POINT('',(1.,0.));#3=CARTESIAN_POINT('',(2.,0.));
#4=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.,(2,1,2),(0.,1.,2.),.UNSPECIFIED.);
#5=CARTESIAN_POINT('',(0.5,0.));
#6=B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#5,#3),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,2.),.UNSPECIFIED.);
#10=CARTESIAN_POINT('',(0.,0.,0.));#11=CARTESIAN_POINT('',(0.,1.,0.));
#12=CARTESIAN_POINT('',(2.,0.,0.));#13=CARTESIAN_POINT('',(2.,1.,0.));
#14=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#10,#11),(#12,#13)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,2.),(0.,1.),.UNSPECIFIED.);
#20=(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''));
#21=DEFINITIONAL_REPRESENTATION('',(#4),#20);#22=PCURVE('',#14,#21);
#25=DEFINITIONAL_REPRESENTATION('',(#6),#20);#26=PCURVE('',#14,#25);
#31=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#3),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.5),.UNSPECIFIED.);
#32=DEFINITIONAL_REPRESENTATION('',(#31),#20);#33=PCURVE('',#14,#32);
#34=CARTESIAN_POINT('',(0.,99.));#35=CARTESIAN_POINT('',(2.,99.));
#36=B_SPLINE_CURVE_WITH_KNOTS('',1,(#34,#35),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,2.),.UNSPECIFIED.);
#38=DEFINITIONAL_REPRESENTATION('',(#36),#20);#37=PCURVE('',#14,#38);
#50=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#12),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,2.),.UNSPECIFIED.);
ENDSEC;END-ISO-10303-21;";

/// The gap rule for imports on one edge use (decision Q9): "measure the
/// geometric gap both ways, reparametrise, then continue with the
/// same-parameter gap". The checker reads the file's own carriers through
/// their binary64 evaluators and the u/v curve as the file wrote it.
#[test]
fn a_file_uv_curve_at_another_speed_is_given_the_parameter_of_its_edge() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SPEEDS, &mut dec).unwrap();
    let scale = Rational::one();
    let surface =
        NativeSurface::from_source(&parsed, 14, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let carrier = PhysicalCurve::Native(Arc::new(
        NativeCurve::from_source(&parsed, 50, &scale, AngleUnit::Radians, &mut dec).unwrap(),
    ));
    let mut points = Vec::new();
    for station in 0..=STATIONS {
        let t = 2.0 * station as f64 / STATIONS as f64;
        points.push(carrier.sample_at(t, 0, &mut dec).unwrap());
    }
    let edge = EdgeSamples {
        points,
        source_curve: Some(50),
        carrier,
    };
    let mut context = Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), 0.01);
    let mut read = |entity: u64, dec: &mut Decider| {
        let pcurve =
            NativeUvCurve::from_source(&parsed, entity, &scale, AngleUnit::Radians, dec).unwrap();
        import_gap::read_use(
            &mut context,
            &import_gap::FileUse {
                edge: &edge,
                surface: (14, &surface),
                pcurve: (entity, &pcurve),
                overshoot_mm: 1e-7,
            },
            dec,
        )
        .unwrap()
    };

    // At the edge's own speed nothing could be lowered: the same-parameter
    // gap stores rule T1's least floor or less, and the rule leaves the use.
    let same = read(22, &mut dec).unwrap_or_else(|why| panic!("{why}"));
    assert_eq!(same.rule.outcome, None);
    assert!(same.rule.geometric.is_none());
    assert!(same.copy.is_none());
    assert!(same.rule.stored_mm() <= 1e-7, "{}", same.rule.stored_mm());

    // At another speed the gap is all parameter.
    let (stretched, work) = with_operation_budget(20_000_000, || {
        read(26, &mut Decider::default()).unwrap_or_else(|why| panic!("{why}"))
    });
    let rule = &stretched.rule;
    assert!((rule.parametric.raw_mm - 0.25).abs() < 1e-9, "{rule:?}");
    let geometric = rule.geometric.expect("a geometric reading");
    assert!(geometric.reported_mm() < 1e-9, "{geometric:?}");
    assert_eq!(rule.outcome, Some(n::ImportGapOutcome::Reparametrised));
    assert!(rule.stored_mm() <= 1e-7, "{}", rule.stored_mm());
    let copy = stretched.copy.expect("the copy in the edge's parameter");
    for (location, _) in &edge.points {
        let t = location.parameters[0];
        let [u, v] = copy.uv(t);
        assert!((u - t).abs() < 1e-9 && v.abs() < 1e-12, "{t}: {u}, {v}");
    }
    // The checker's units and the chart work of the native evaluations are
    // both charged to the operation.
    assert!(work.used > rule.work, "{work:?} for {}", rule.work);
    println!(
        "stretched u/v curve: same-parameter {:.6e} mm, geometric {:.3e} mm, {:.3e} mm stored in {} span(s), {} units of the use, {} charged",
        rule.parametric.raw_mm,
        geometric.reported_mm(),
        rule.stored_mm(),
        copy.spans,
        rule.work,
        work.used
    );

    // A use the checker cannot read is said to be unread, with the reason.
    let short = read(33, &mut dec)
        .err()
        .expect("a u/v curve that stops short");
    assert!(
        short.contains("does not hold the edge's parameter range"),
        "{short}"
    );
    let off = read(37, &mut dec).err().expect("a u/v curve off the patch");
    assert!(off.contains("leaves the charts of its surface"), "{off}");
}

/// De Boor's recurrence on the trace a change of parameter leaves: a cubic
/// with triple interior knots, read at its knots, between them and past its
/// ends.
#[test]
fn a_reparametrised_trace_is_evaluated_as_its_spline() {
    // u = t^3 on [0, 1] and u = 1 + 3 (t - 1) on [1, 2], v = 2 t.
    let trace = crate::tolerant::admission::ReparametrisedTrace {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
        controls: vec![
            [0.0, 0.0],
            [0.0, 2.0 / 3.0],
            [0.0, 4.0 / 3.0],
            [1.0, 2.0],
            [2.0, 8.0 / 3.0],
            [3.0, 10.0 / 3.0],
            [4.0, 4.0],
        ],
        spans: 2,
    };
    for (t, u) in [
        (0.0, 0.0),
        (0.5, 0.125),
        (1.0, 1.0),
        (1.5, 2.5),
        (2.0, 4.0),
        (2.25, 4.75),
        (-0.5, -0.125),
    ] {
        let [found, v] = trace.uv(t);
        assert!((found - u).abs() < 1e-14, "{t}: {found} for {u}");
        assert!((v - 2.0 * t).abs() < 1e-14, "{t}: {v}");
    }
}
