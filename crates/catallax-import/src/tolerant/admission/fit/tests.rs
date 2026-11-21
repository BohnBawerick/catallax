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
/// at `t = 1.5`, short of the edge's end, `#37` runs along `v = 99`, off the
/// patch, and `#46` is the curve of `#22` on an interval of its own,
/// `[10, 12]`.
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
#44=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.,(2,1,2),(10.,11.,12.),.UNSPECIFIED.);
#45=DEFINITIONAL_REPRESENTATION('',(#44),#20);#46=PCURVE('',#14,#45);
#50=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#12),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,2.),.UNSPECIFIED.);
ENDSEC;END-ISO-10303-21;";

/// The patch of [`SPEEDS`] and the stations of its straight edge.
fn speeds_use(
    parsed: &ParsedStep,
    scale: &Rational,
    dec: &mut Decider,
) -> (NativeSurface, EdgeSamples) {
    let surface = NativeSurface::from_source(parsed, 14, scale, AngleUnit::Radians, dec).unwrap();
    let carrier = PhysicalCurve::Native(Arc::new(
        NativeCurve::from_source(parsed, 50, scale, AngleUnit::Radians, dec).unwrap(),
    ));
    let mut points = Vec::new();
    for station in 0..=STATIONS {
        let t = 2.0 * station as f64 / STATIONS as f64;
        points.push(carrier.sample_at(t, 0, dec).unwrap());
    }
    let edge = EdgeSamples {
        points,
        source_curve: Some(50),
        carrier,
    };
    (surface, edge)
}

/// The gap rule for imports on one edge use (decision Q9): "measure the
/// geometric gap both ways, reparametrise, then continue with the
/// same-parameter gap". The checker reads the file's own carriers through
/// their binary64 evaluators and the u/v curve as the file wrote it.
#[test]
fn a_file_uv_curve_at_another_speed_is_given_the_parameter_of_its_edge() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SPEEDS, &mut dec).unwrap();
    let scale = Rational::one();
    let (surface, edge) = speeds_use(&parsed, &scale, &mut dec);
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
                allowance: n::gap::Limits::default(),
            },
            dec,
        )
        .unwrap()
    };

    // At the edge's own speed nothing could be lowered: the same-parameter
    // gap stores rule T1's least floor or less, and the rule leaves the use.
    let same = read(22, &mut dec).unwrap_or_else(|unread| panic!("{}", unread.why));
    assert_eq!(same.rule.outcome, None);
    assert!(same.rule.geometric.is_none());
    assert!(same.copy.is_none());
    assert!(same.rule.stored_mm() <= 1e-7, "{}", same.rule.stored_mm());

    // At another speed the gap is all parameter.
    let (stretched, work) = with_operation_budget(20_000_000, || {
        read(26, &mut Decider::default()).unwrap_or_else(|unread| panic!("{}", unread.why))
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
    // The use's work holds the checker's units and the chart work of the
    // native evaluations, and all of it is charged to the operation.
    assert!(work.used >= rule.work, "{work:?} for {}", rule.work);
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
    // These three are found unread before the checker works on them, so
    // they carry no work.
    let short = read(33, &mut dec)
        .err()
        .expect("a u/v curve that stops short");
    assert!(
        short
            .why
            .contains("does not hold the edge's parameter range"),
        "{}",
        short.why
    );
    let off = read(37, &mut dec).err().expect("a u/v curve off the patch");
    assert!(
        off.why.contains("leaves the charts of its surface"),
        "{}",
        off.why
    );
    // The curve of `#22` on an interval of its own is the same segment, and
    // the rule does not read it: its interval is not its edge's.
    let shifted = read(46, &mut dec)
        .err()
        .expect("a u/v curve on a shifted interval");
    assert!(
        shifted
            .why
            .contains("does not hold the edge's parameter range"),
        "{}",
        shifted.why
    );
    assert_eq!([short.work, off.work, shifted.work], [0, 0, 0]);
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

/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." The work
/// a native chart charges for an evaluation is inside the use's one
/// allowance: the checker takes it, with its own unit, before it evaluates.
///
/// The stretched use of [`SPEEDS`] is read as it is, and again with its two
/// carriers stating no work, which leaves the checker's own units. Under an
/// allowance that fits those units exactly, the use as it is refuses by
/// name, while the request's own allowance is ample. So does a small cap,
/// and so does one unit short of the use's total.
#[test]
fn the_chart_work_of_a_native_use_is_inside_its_one_allowance() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SPEEDS, &mut dec).unwrap();
    let scale = Rational::one();
    let (surface, edge) = speeds_use(&parsed, &scale, &mut dec);
    let mut context = Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), 0.01);
    let pcurve =
        NativeUvCurve::from_source(&parsed, 26, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let used = import_gap::FileUse {
        edge: &edge,
        surface: (14, &surface),
        pcurve: (26, &pcurve),
        overshoot_mm: 1e-7,
        allowance: n::gap::Limits::default(),
    };
    let mut sampled = import_gap::sampled_use(&mut context, &used, &mut dec)
        .unwrap()
        .unwrap_or_else(|why| panic!("{why}"));
    let stated = [sampled.curve.work, sampled.surface.work];
    assert!(stated[0] > 0 && stated[1] > 0, "{stated:?}");
    let read = |sampled: &import_gap::SampledUse, allowance: u64| {
        with_operation_budget(20_000_000, || {
            sampled.read(n::gap::Limits { work: allowance }, &mut Decider::default())
        })
    };
    let total = |read: Result<Result<import_gap::UseReading, import_gap::Unread>, Error>| {
        let read = read
            .unwrap()
            .unwrap_or_else(|unread| panic!("{}", unread.why));
        (read.rule.work, read.rule.stored_mm().to_bits())
    };
    let refused =
        |read: Result<Result<import_gap::UseReading, import_gap::Unread>, Error>| match read {
            Err(Error::Numerical(refusal)) => refusal.code,
            Err(other) => panic!("{other:?}"),
            Ok(_) => panic!("the use was read"),
        };

    let (as_it_is, charged) = read(&sampled, n::gap::EDGE_WORK);
    let (whole, stored) = total(as_it_is);
    assert!(charged.used >= whole, "{charged:?} for {whole}");
    (sampled.curve.work, sampled.surface.work) = (0, 0);
    let (checker, same) = total(read(&sampled, n::gap::EDGE_WORK).0);
    assert_eq!(same, stored);
    assert!(whole > checker, "{whole} and {checker}");
    // The checker's own work fits an allowance of exactly that much.
    assert_eq!(total(read(&sampled, checker).0), (checker, stored));
    println!(
        "stretched u/v curve: {checker} units of the checker alone, {whole} with {} and {} units an evaluation of the curve and the surface",
        stated[0], stated[1]
    );

    // With the chart work the same allowance is not enough, and the use
    // refuses by name. The request has spent a small part of its own.
    (sampled.curve.work, sampled.surface.work) = (stated[0], stated[1]);
    let (over, request) = read(&sampled, checker);
    assert_eq!(refused(over), "numerical.gap-work-cap");
    assert!(
        !request.exhausted && request.used <= 2 * checker,
        "{request:?}"
    );
    // A small cap, and one unit short of the whole.
    assert_eq!(refused(read(&sampled, 100).0), "numerical.gap-work-cap");
    assert_eq!(
        refused(read(&sampled, whole - 1).0),
        "numerical.gap-work-cap"
    );
    assert_eq!(total(read(&sampled, whole).0), (whole, stored));
}

/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." The cap is
/// one allowance for the topological edge, shared by its uses.
///
/// The stretched use of [`SPEEDS`] is read twice, as the two uses of one
/// edge on two faces. An allowance that fits the use once, and not twice,
/// reads the first use and refuses the second by name. One unit short of
/// twice refuses the second as well. Exactly twice reads both. The uses of
/// two different edges each have an allowance of their own.
#[test]
fn the_uses_of_an_edge_share_its_one_allowance() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SPEEDS, &mut dec).unwrap();
    let scale = Rational::one();
    let (surface, edge) = speeds_use(&parsed, &scale, &mut dec);
    let mut context = Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), 0.01);
    let pcurve =
        NativeUvCurve::from_source(&parsed, 26, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let used = import_gap::FileUse {
        edge: &edge,
        surface: (14, &surface),
        pcurve: (26, &pcurve),
        overshoot_mm: 1e-7,
        allowance: n::gap::Limits::default(),
    };
    let sampled = import_gap::sampled_use(&mut context, &used, &mut dec)
        .unwrap()
        .unwrap_or_else(|why| panic!("{why}"));
    // The work of each use of the edges named, read in turn on `allowances`,
    // or the code the first refused use gives.
    let read_uses = |allowances: &mut import_gap::EdgeAllowances, edges: [usize; 2]| {
        let mut counted = Vec::new();
        for edge in edges {
            let read = allowances.read(edge, |allowance| {
                sampled.read(allowance, &mut Decider::default())
            });
            match read {
                Ok(read) => counted.push(
                    read.unwrap_or_else(|unread| panic!("{}", unread.why))
                        .rule
                        .work,
                ),
                Err(Error::Numerical(refusal)) => return Err((counted, refusal.code)),
                Err(other) => panic!("{other:?}"),
            }
        }
        Ok(counted)
    };
    let whole = read_uses(&mut import_gap::EdgeAllowances::new(), [7, 7]).unwrap();
    assert_eq!(whole[0], whole[1]);
    let once = whole[0];
    assert!(once > 0 && 2 * once < n::gap::EDGE_WORK);

    for each_edge in [once, 2 * once - 1] {
        assert_eq!(
            read_uses(&mut import_gap::EdgeAllowances::of(each_edge), [7, 7]),
            Err((vec![once], "numerical.gap-work-cap")),
            "{each_edge} units for an edge whose two uses count {}",
            2 * once
        );
    }
    assert_eq!(
        read_uses(&mut import_gap::EdgeAllowances::of(2 * once), [7, 7]),
        Ok(vec![once, once])
    );
    // Two edges, one use each.
    assert_eq!(
        read_uses(&mut import_gap::EdgeAllowances::of(once), [7, 8]),
        Ok(vec![once, once])
    );
}

/// A use that is found unread after the checker has worked on it has still
/// spent that work, and its edge's next use does not get it again (section
/// 4.2 step 8: "A hard per-edge work cap refuses by name").
///
/// The stretched use of [`SPEEDS`] is read as the two uses of one edge. On
/// the first, the rule's new copy is forced outside the charts of the
/// surface, so the rule finishes, counts its work, and the use is unread:
/// the caller rebuilds its trace and goes on. With an allowance one unit
/// short of the two uses together, the second use of the edge then refuses
/// by name. With exactly the two together it is read. An allowance that
/// took nothing off for the unread use would have read the second use under
/// both.
#[test]
fn an_unread_use_has_spent_its_work_on_the_edge() {
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(SPEEDS, &mut dec).unwrap();
    let scale = Rational::one();
    let (surface, edge) = speeds_use(&parsed, &scale, &mut dec);
    let mut context = Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), 0.01);
    let pcurve =
        NativeUvCurve::from_source(&parsed, 26, &scale, AngleUnit::Radians, &mut dec).unwrap();
    let used = import_gap::FileUse {
        edge: &edge,
        surface: (14, &surface),
        pcurve: (26, &pcurve),
        overshoot_mm: 1e-7,
        allowance: n::gap::Limits::default(),
    };
    let mut sampled = import_gap::sampled_use(&mut context, &used, &mut dec)
        .unwrap()
        .unwrap_or_else(|why| panic!("{why}"));
    let whole_charts = sampled.extent;
    let once = sampled
        .read(n::gap::Limits::default(), &mut dec)
        .unwrap()
        .unwrap_or_else(|unread| panic!("{}", unread.why))
        .rule
        .work;
    assert!(once > 0);

    // The same use with an extent in u that no point is inside, so the new
    // copy the rule makes is found outside the charts after the rule ran.
    sampled.extent[0][0] = Some(f64::INFINITY);
    let unread = sampled
        .read(n::gap::Limits::default(), &mut dec)
        .unwrap()
        .err()
        .expect("a new copy that leaves the charts");
    assert_eq!(
        unread,
        import_gap::Unread {
            why: "its copy in the edge's parameter leaves the charts of its surface",
            work: once,
        }
    );

    // The first use of edge 7 is unread and the second is read as it is.
    let as_it_is = import_gap::sampled_use(&mut context, &used, &mut dec)
        .unwrap()
        .unwrap_or_else(|why| panic!("{why}"));
    assert_eq!(as_it_is.extent, whole_charts);
    let second_use = |each_edge: u64| {
        let mut allowances = import_gap::EdgeAllowances::of(each_edge);
        let first = allowances
            .read(7, |allowance| {
                sampled.read(allowance, &mut Decider::default())
            })
            .unwrap();
        assert_eq!(first.err().map(|unread| unread.work), Some(once));
        assert_eq!(allowances.left(7).work, each_edge - once);
        match allowances.read(7, |allowance| {
            as_it_is.read(allowance, &mut Decider::default())
        }) {
            Ok(read) => Ok(read
                .unwrap_or_else(|unread| panic!("{}", unread.why))
                .rule
                .work),
            Err(Error::Numerical(refusal)) => Err(refusal.code),
            Err(other) => panic!("{other:?}"),
        }
    };
    assert_eq!(second_use(2 * once - 1), Err("numerical.gap-work-cap"));
    assert_eq!(second_use(once + 1), Err("numerical.gap-work-cap"));
    assert_eq!(second_use(2 * once), Ok(once));
}
