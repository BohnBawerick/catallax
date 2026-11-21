//! Check once (roadmap item 3.0A step 5; Checkpoint 2 decisions, section 4.1
//! item 5): "Reopening trusts the tolerance recorded at export."
//!
//! Each test writes a numerical body to STEP and opens the file again. The
//! same bytes with the stamp's first word changed are the control: a file
//! that records nothing is measured, as every import is (rule T3).
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;
#[path = "../../catallax-step/tests/common/n03_tube.rs"]
mod tube;
use catallax_import::{
    recorded::{RecordedField, recorded_tolerance},
    sha256_hex,
    tolerant::opened::{BodyView, OpenedStep, open_step},
};
use catallax_numbers::{Decider, NumericalCheck, OperationWork, Rational, with_operation_budget};
use catallax_recipe::StepId;
use catallax_step::numerical::{self, LengthScale, NumericalStepOptions};
use catallax_topology::numerical::{
    self as n, NumericalBody, NumericalOptions,
    tolerance::{RECORDED_USE, RECORDED_VERTEX, ReceiverProfile},
};
use fixture::{q, quarter};

fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
        vertex_tolerance_mm: None,
    }
}

fn written(body: &NumericalBody) -> (Vec<u8>, Rational) {
    let file = numerical::build(
        body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    (file.bytes().to_vec(), file.report().uncertainty_mm.clone())
}

/// The same file as another system's: the stamp no longer names a Catallax
/// numerical export, so it records nothing.
fn unstamped(bytes: &[u8]) -> Vec<u8> {
    let text = std::str::from_utf8(bytes).unwrap();
    assert!(text.contains("catallax numerical approximate;"));
    text.replacen(
        "catallax numerical approximate;",
        "another system approximate;",
        1,
    )
    .into_bytes()
}

/// Opens under an allowance large enough to count all the work.
fn open(bytes: &[u8]) -> (OpenedStep, OperationWork) {
    let (opened, work) = with_operation_budget(100_000_000_000, || {
        open_step(
            bytes,
            &sha256_hex(bytes),
            &StepId::new("import").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    (opened.unwrap(), work)
}

fn names(body: &NumericalBody) -> Vec<&str> {
    body.evidence()
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, .. } => Some(name.as_str()),
            NumericalCheck::Passed { .. } => None,
        })
        .collect()
}

/// The quarter cylinder with one of its side planes a micrometre off the
/// edges that bound it. With every face on its edges the written file is
/// exact and the certified import takes it; this one has a real gap, so the
/// numerical import takes it.
fn gapped_quarter() -> NumericalBody {
    gapped_quarter_named(&quarter().name)
}

/// [`gapped_quarter`] under a name of the caller's.
fn gapped_quarter_named(name: &str) -> NumericalBody {
    let mut definition = quarter();
    definition.name = name.into();
    let side = definition
        .faces
        .iter_mut()
        .find(|face| face.key == "f2")
        .expect("the quarter has a face f2");
    let catallax_shapes::Surface::Plane(plane) = &mut side.surface else {
        panic!("face f2 of the quarter is a plane");
    };
    plane.origin.y = &plane.origin.y + &catallax_numbers::Real::from(q("0.000001"));
    n::build(definition, options(), &mut Decider::default()).unwrap()
}

/// The numerical import of a file, under the default allowance.
fn import(bytes: &[u8]) -> (NumericalBody, OperationWork) {
    let (imported, work) = with_operation_budget(20_000_000, || {
        catallax_import::numerical::import_step(
            bytes,
            &sha256_hex(bytes),
            &StepId::new("import").unwrap(),
            &options(),
            &mut Decider::default(),
        )
    });
    let catallax_import::numerical::StepImport::Numerical(imported) = imported.unwrap() else {
        panic!("the gapped quarter is a numerical body");
    };
    assert!(!work.exhausted);
    (imported.bodies()[0].1.clone(), work)
}

#[test]
fn an_analytic_file_is_reopened_with_the_tolerance_it_records() {
    let body = gapped_quarter();
    let (bytes, recorded) = written(&body);
    assert!(recorded > q("0.000001"), "{recorded}");
    let read = recorded_tolerance(&bytes).unwrap();
    assert_eq!(read.tolerance_mm, recorded);
    assert_eq!(read.field, RecordedField::Recorded);

    let (reopened, work) = import(&bytes);
    // Every edge end and every edge use holds the recorded tolerance under
    // its own name, and nothing was measured.
    let names = names(&reopened);
    let uses: usize = reopened
        .definition()
        .faces
        .iter()
        .map(|face| face.loops.iter().map(Vec::len).sum::<usize>())
        .sum();
    assert_eq!(
        names.iter().filter(|name| **name == RECORDED_USE).count(),
        uses
    );
    assert_eq!(
        names
            .iter()
            .filter(|name| **name == RECORDED_VERTEX)
            .count(),
        reopened.definition().edges.len()
    );
    assert!(
        names
            .iter()
            .all(|name| *name == RECORDED_USE || *name == RECORDED_VERTEX),
        "{names:?}"
    );
    // The recorded tolerance is taken as it is: rule T4's margin is not
    // applied to it a second time.
    for entity in &reopened.tolerances().entities {
        assert_eq!(entity.tolerance_mm, recorded, "{entity:?}");
    }
    assert!(reopened.tolerances().warnings.is_empty());

    // The control is measured: its records are distances, with no recorded
    // name among them, and it costs more.
    let (control, measured) = import(&unstamped(&bytes));
    let control = self::names(&control);
    assert!(control.contains(&"trim.edge-surface-distance"));
    assert!(!control.contains(&RECORDED_USE) && !control.contains(&RECORDED_VERTEX));
    assert!(
        work.used < measured.used,
        "recorded {} and measured {}",
        work.used,
        measured.used
    );
    println!(
        "analytic quarter: imported with its record in {} work, measured in {} work",
        work.used, measured.used
    );

    // Opened as the application opens a file, the same body comes back.
    let opened = open_step(
        &bytes,
        &sha256_hex(&bytes),
        &StepId::new("import").unwrap(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let Some(BodyView::Numerical(reopened)) = opened.body(0) else {
        panic!(
            "the quarter reopens as a numerical body: {:?}, {:?}",
            opened.bodies()[0].tier(),
            opened.first_refusal()
        );
    };
    assert!(self::names(reopened).contains(&RECORDED_USE));
}

/// The square tube with its top face a flat spline patch a micrometre above
/// its edges. A spline face with a real gap is taken by neither the certified
/// import nor the analytic one, so the file opens as a native tolerant body,
/// the way the R2 rim does.
fn spline_topped_tube() -> NumericalBody {
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch, Surface};
    let mut definition = tube::tube();
    let mut dec = Decider::default();
    let lifted = catallax_numbers::Real::from(q("1.000001"));
    let controls = [(0, 0), (0, 4), (4, 0), (4, 4)]
        .map(|(x, y)| Point3::new(x.into(), y.into(), lifted.clone()))
        .to_vec();
    let knots = vec![q("0"), q("0"), q("1"), q("1")];
    let surface = BSplineSurface::new(
        [1, 1],
        [knots.clone(), knots],
        [2, 2],
        controls,
        vec![catallax_numbers::Real::one(); 4],
        &mut dec,
    )
    .unwrap();
    definition.faces[0].surface =
        Surface::BSpline(RationalPatch::new(surface, &mut dec).unwrap().into());
    for trim in definition.faces[0].loops.iter_mut().flatten() {
        for uv in trim.uv.iter_mut().flatten() {
            *uv = uv.clone() / q("4");
        }
    }
    n::build(definition, options(), &mut dec).unwrap()
}

#[test]
fn a_file_with_a_spline_face_is_reopened_with_the_tolerance_it_records() {
    let body = spline_topped_tube();
    let (bytes, recorded) = written(&body);
    assert!(recorded > q("0.000001"), "{recorded}");

    let (opened, work) = open(&bytes);
    let Some(BodyView::Tolerant {
        import,
        body: index,
    }) = opened.body(0)
    else {
        panic!(
            "the tube opens as a native tolerant body: {:?}, {:?}",
            opened.bodies()[0].tier(),
            opened.first_refusal()
        );
    };
    let solid = &import.bodies()[index];
    assert_eq!(solid.recorded_tolerance_mm(), Some(&recorded));
    let evidence = import.evidence();
    // Every vertex, edge and face holds the recorded tolerance, and no
    // residual was measured against it.
    let held = recorded.to_f64_nearest();
    assert!(!evidence.tolerances.is_empty());
    for entity in &evidence.tolerances {
        assert_eq!(entity.tolerance_mm.to_bits(), held.to_bits(), "{entity:?}");
        assert_eq!(entity.samples, 0, "{entity:?}");
    }
    let checks: Vec<_> = evidence.sampled_checks.iter().map(|c| c.check).collect();
    assert!(checks.contains(&"face.recorded-tolerance-written-uv"));
    assert!(checks.contains(&"edge.recorded-tolerance"));
    assert!(
        checks
            .iter()
            .all(|check| check.contains("recorded-tolerance")),
        "{checks:?}"
    );
    // The file's own u/v curves are the face copies: none was rebuilt.
    assert!(evidence.native_trims.iter().all(|trim| !trim.reconstructed));
    assert!(evidence.repairs.is_empty(), "{:?}", evidence.repairs);
    assert!(
        evidence
            .warnings
            .iter()
            .any(|warning| warning.code == "recorded-tolerance")
    );

    // The control is measured, and its native admission costs more.
    let (control, _) = open(&unstamped(&bytes));
    if let Some(BodyView::Tolerant {
        import: measured,
        body: index,
    }) = control.body(0)
    {
        assert_eq!(measured.bodies()[index].recorded_tolerance_mm(), None);
        assert!(
            measured
                .evidence()
                .sampled_checks
                .iter()
                .any(|check| check.check == "face.native-uv-sampled-fit")
        );
    }
    assert!(
        opened.native_admission_work_used() < control.native_admission_work_used(),
        "recorded {} and measured {}",
        opened.native_admission_work_used(),
        control.native_admission_work_used()
    );
    println!(
        "spline-topped tube: native admission {} work with its record, {} work measured ({:?}); whole reopen {} work",
        opened.native_admission_work_used(),
        control.native_admission_work_used(),
        control.bodies()[0].tier(),
        work.used
    );
}

/// A record is used only where the reader can take the file's entities as
/// they are written. Here the interval of the retained u/v curve is cut
/// short by hand, so its edge runs a tenth of its length past it: far more
/// than the recorded tolerance allows. The recorded pass refuses the body,
/// says why, and the measured passes admit it as they do any import.
#[test]
fn a_body_the_reader_cannot_take_as_written_is_measured() {
    let body = spline_topped_tube();
    let (bytes, _) = written(&body);
    let text = std::str::from_utf8(&bytes).unwrap();
    // One u/v curve of the spline face, with its interval cut from [0, 1]
    // to [0, 0.9].
    let curve = text
        .lines()
        .find(|line| {
            line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',1,") && line.contains("(2,2),(0.,1.),")
        })
        .expect("a straight u/v curve over the edge's interval");
    let cut = text
        .replace(curve, &curve.replace("(2,2),(0.,1.),", "(2,2),(0.,0.9),"))
        .into_bytes();
    assert_ne!(cut, bytes);
    assert!(recorded_tolerance(&cut).is_some());

    let (opened, _) = open(&cut);
    let import = opened
        .tolerant_admission()
        .expect("the native admission ran");
    let unused: Vec<_> = import
        .evidence()
        .warnings
        .iter()
        .filter(|warning| warning.code == "recorded-tolerance-unused")
        .collect();
    assert_eq!(unused.len(), 1, "{:?}", import.evidence().warnings);
    assert!(
        unused[0]
            .detail
            .contains("runs past the interval of its u/v curve"),
        "{}",
        unused[0].detail
    );
    // Whatever the measured passes make of the cut file, nothing in it
    // holds the recorded tolerance.
    assert!(
        import
            .bodies()
            .iter()
            .all(|solid| solid.recorded_tolerance_mm().is_none())
    );
    assert!(
        import
            .evidence()
            .sampled_checks
            .iter()
            .all(|check| !check.check.contains("recorded-tolerance"))
    );
}

/// One straight u/v curve of the spline face, written again as a quadratic
/// that traces the same segment at another speed: its middle control sits
/// four tenths of the way along, so the copy lags its edge by up to a
/// twentieth of the segment at the same parameter.
fn with_one_uv_curve_at_another_speed(text: &str) -> String {
    let curve = text
        .lines()
        .find(|line| {
            line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',1,") && line.contains("(2,2),(0.,1.),")
        })
        .expect("a straight u/v curve over the edge's interval");
    let controls = curve
        .split("'',1,(")
        .nth(1)
        .and_then(|rest| rest.split(')').next())
        .expect("two controls");
    let ends: Vec<&str> = controls.split(',').collect();
    assert_eq!(ends.len(), 2, "{curve}");
    let uv = |name: &str| -> [f64; 2] {
        let line = text
            .lines()
            .find(|line| line.starts_with(&format!("{name}=CARTESIAN_POINT(")))
            .expect("a control point");
        let inner = line
            .rsplit("(")
            .next()
            .and_then(|rest| rest.split(')').next())
            .expect("coordinates");
        let values: Vec<f64> = inner.split(',').map(|x| x.parse().unwrap()).collect();
        [values[0], values[1]]
    };
    let (first, last) = (uv(ends[0]), uv(ends[1]));
    let middle = [0, 1].map(|axis| first[axis] + 0.4 * (last[axis] - first[axis]));
    let next = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix('#')?
                .split('=')
                .next()?
                .parse::<u64>()
                .ok()
        })
        .max()
        .expect("an entity")
        + 1;
    let quadratic = curve
        .replace(
            &format!("'',1,({controls})"),
            &format!("'',2,({},#{next},{})", ends[0], ends[1]),
        )
        .replace("(2,2),(0.,1.),", "(3,3),(0.,1.),");
    assert_ne!(quadratic, curve);
    text.replace(
        curve,
        &format!(
            "#{next}=CARTESIAN_POINT('',({:.17},{:.17}));\n{quadratic}",
            middle[0], middle[1]
        ),
    )
}

/// A file that records nothing is an import, and its u/v curves are copies
/// of its own (decision Q9: "Imports: geometric gap first; if the parametric
/// gap is much larger, reparametrise instead of raising the tolerance").
/// The native tolerant admission reads every edge use that carries one
/// through the gap rule for imports, the code the numerical import runs.
///
/// The spline-topped tube as written has every u/v curve at its edge's
/// speed, and the rule keeps them all. With one of them at another speed the
/// same-parameter gap is tenths of a millimetre and the geometric gap is the
/// micrometre the top is off its edges. That copy is given its edge's
/// parameter. The edge's tolerance covers the micrometre and not the speed,
/// and the trace is not rebuilt from the 3D curve.
#[test]
fn a_file_that_records_nothing_is_read_by_the_gap_rule_for_imports() {
    let body = spline_topped_tube();
    let (bytes, _) = written(&body);
    let plain = unstamped(&bytes);

    let (opened, _) = open(&plain);
    let import = opened
        .tolerant_admission()
        .expect("the native admission ran");
    let evidence = import.evidence();
    assert!(!evidence.import_gaps.is_empty());
    for row in &evidence.import_gaps {
        assert!(row.kept, "{row:?}");
        assert_ne!(row.outcome, "reparametrised", "{row:?}");
        assert!(row.parametric_mm < 2e-6, "{row:?}");
    }
    assert!(
        evidence
            .native_trims
            .iter()
            .all(|trim| trim.reparametrised.is_none())
    );
    assert!(
        evidence
            .warnings
            .iter()
            .all(|warning| warning.code != "import-gap-rule-unread"),
        "{:?}",
        evidence.warnings
    );

    let stretched =
        with_one_uv_curve_at_another_speed(std::str::from_utf8(&plain).unwrap()).into_bytes();
    let (opened, work) = open(&stretched);
    let import = opened
        .tolerant_admission()
        .expect("the native admission ran");
    assert_eq!(import.bodies().len(), 1, "{:?}", opened.first_refusal());
    let evidence = import.evidence();
    let changed: Vec<_> = evidence
        .import_gaps
        .iter()
        .filter(|row| row.outcome == "reparametrised")
        .collect();
    assert_eq!(changed.len(), 1, "{:?}", evidence.import_gaps);
    let row = changed[0];
    let geometric = row.geometric_mm.expect("a geometric reading");
    let next = row
        .reparametrised_mm
        .expect("a copy in the edge's parameter");
    println!(
        "tube with one u/v curve at another speed: same-parameter {:.6e} mm stored, geometric {:.6e} mm, reparametrised {:.6e} mm stored, {} units; native admission {} work, whole reopen {} work",
        row.parametric_mm,
        geometric,
        next,
        row.work,
        opened.native_admission_work_used(),
        work.used
    );
    assert!(row.kept, "{row:?}");
    assert!(row.parametric_mm > 0.05, "{row:?}");
    // The top is a micrometre off its edges: that is the geometric gap, and
    // the new copy stores rule T4's tolerance of it.
    assert!((0.9e-6..1.1e-6).contains(&geometric), "{row:?}");
    assert!(next >= 1.1 * 0.9e-6 && next <= 1.21 * geometric, "{row:?}");

    // The use holds the file's curve at its edge's parameter: no source
    // PCURVE, and no trace rebuilt from the 3D curve.
    let trim = evidence
        .native_trims
        .iter()
        .find(|trim| {
            trim.face == row.face
                && trim.boundary == row.boundary
                && trim.use_index == row.use_index
        })
        .expect("the use");
    assert!(trim.reparametrised.is_some());
    assert!(trim.source_pcurve.is_none());
    assert!(!trim.reconstructed);
    assert!(evidence.repairs.is_empty(), "{:?}", evidence.repairs);
    // The edge's tolerance covers the new copy's gap, and nothing of the
    // tenths of a millimetre of speed.
    let edge = evidence
        .tolerances
        .iter()
        .find(|entity| entity.kind == "edge" && entity.entity == row.edge)
        .expect("the edge");
    assert!(edge.tolerance_mm >= next, "{edge:?}");
    assert!(edge.tolerance_mm < 1e-5, "{edge:?}");

    // The body's own STEP export writes the copy it holds.
    let exported = import.export_step(&mut Decider::default()).unwrap();
    assert!(
        std::str::from_utf8(exported.bytes())
            .unwrap()
            .contains("source trace in the parameter of its 3D curve")
    );
}

/// A file's u/v curve is kept as a use's face copy only where the gap rule
/// read it (section 4.2 step 2: the geometric gap is read both ways before
/// anything continues). One the rule does not read is never kept on the
/// strength of the stations: its trace is rebuilt from the 3D curve.
///
/// Two controls on the unstamped tube. In the first, one u/v curve of the
/// top is on an interval of its own, `[10, 11]`, off its edge's `[0, 1]`. In
/// the second it is the same segment written as a polyline, which agrees
/// with its edge at every station and which the rule has no reader for.
///
/// The pass that repairs nothing refuses the body for that one use. The
/// repair pass reads a straight edge in its line's own length, not in the
/// file's parameter, so there every file curve of the top is on another
/// interval than its edge, is unread, and is rebuilt.
#[test]
fn a_file_uv_curve_the_gap_rule_does_not_read_is_never_kept() {
    let body = spline_topped_tube();
    let (bytes, _) = written(&body);
    let plain = String::from_utf8(unstamped(&bytes)).unwrap();
    let curve = plain
        .lines()
        .find(|line| {
            line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',1,") && line.contains("(2,2),(0.,1.),")
        })
        .expect("a straight u/v curve over the edge's interval")
        .to_owned();
    let name = curve.split('=').next().expect("an entity name");
    let controls = curve
        .split("'',1,(")
        .nth(1)
        .and_then(|rest| rest.split(')').next())
        .expect("two controls");
    let shifted = plain.replace(&curve, &curve.replace("(2,2),(0.,1.),", "(2,2),(10.,11.),"));
    let polyline = plain.replace(&curve, &format!("{name}=POLYLINE('',({controls}));"));
    assert_ne!(shifted, plain);
    assert_ne!(polyline, plain);

    // As written, the first pass admits the tube and keeps every file curve.
    let (written_open, _) = open(plain.as_bytes());
    let as_written = written_open
        .tolerant_admission()
        .expect("the native admission ran")
        .evidence();
    let kept_as_written = as_written
        .native_trims
        .iter()
        .filter(|trim| trim.source_pcurve.is_some())
        .count();
    assert!(kept_as_written > 0);
    assert!(
        as_written
            .native_trims
            .iter()
            .all(|trim| !trim.reconstructed)
    );

    for (what, text, why) in [
        (
            "a shifted interval",
            shifted,
            "does not hold the edge's parameter range",
        ),
        ("a polyline", polyline, "neither a B-spline nor a line"),
    ] {
        let (opened, _) = open(text.as_bytes());
        let import = opened
            .tolerant_admission()
            .expect("the native admission ran");
        assert_eq!(
            import.bodies().len(),
            1,
            "{what}: {:?}",
            opened.first_refusal()
        );
        let evidence = import.evidence();
        let unread: Vec<_> = evidence
            .warnings
            .iter()
            .filter(|warning| warning.code == "import-gap-rule-unread")
            .collect();
        assert!(
            unread.iter().any(|warning| warning.detail.contains(why)),
            "{what}: {:?}",
            evidence.warnings
        );
        // Each unread use has a trace of this reader's own, and no file
        // curve.
        for warning in &unread {
            let trim = evidence
                .native_trims
                .iter()
                .find(|trim| warning.scope == format!("{}:{}", trim.face, trim.edge))
                .unwrap_or_else(|| panic!("{what}: no use for {}", warning.scope));
            assert!(trim.reconstructed, "{what}: {}", warning.scope);
            assert!(trim.source_pcurve.is_none(), "{what}: {}", warning.scope);
            assert!(trim.reparametrised.is_none(), "{what}: {}", warning.scope);
        }
        // Every file curve that was kept was read by the rule first, and
        // fewer are kept than in the file as written.
        let mut kept = 0;
        for trim in &evidence.native_trims {
            if trim.source_pcurve.is_none() {
                continue;
            }
            kept += 1;
            assert!(
                evidence.import_gaps.iter().any(|row| row.kept
                    && row.face == trim.face
                    && row.boundary == trim.boundary
                    && row.use_index == trim.use_index),
                "{what}: {} on {} kept its file curve unread",
                trim.edge,
                trim.face
            );
        }
        assert!(kept < kept_as_written, "{what}: {kept}");
        // The rebuilt traces are measured: no edge is looser than the
        // micrometre the top is off its edges allows.
        for entity in &evidence.tolerances {
            assert!(entity.tolerance_mm < 1e-5, "{what}: {entity:?}");
        }
    }
}

/// Nothing binds the record to the geometry under it (owner ruling of 20
/// November 2026: no binding and no fingerprint). A file whose wall was
/// moved by hand after export still opens with the tolerance it records, as
/// an XT file opens with the tolerances it stores. Without its stamp the
/// same file is measured, and a cap under the moved wall's gap refuses it.
#[test]
fn nothing_binds_the_record_to_the_geometry_under_it() {
    let body = n::build(quarter(), options(), &mut Decider::default()).unwrap();
    let (bytes, recorded) = written(&body);
    let text = std::str::from_utf8(&bytes).unwrap();
    let reference = |line: &str| -> String {
        let start = line.rfind('(').unwrap();
        let id: String = line[start..]
            .chars()
            .skip_while(|c| *c != '#')
            .skip(1)
            .take_while(char::is_ascii_digit)
            .collect();
        format!("#{id}=")
    };
    let plane = text.lines().find(|s| s.contains("=PLANE('f2',")).unwrap();
    let placement = text
        .lines()
        .find(|s| s.starts_with(&reference(plane)))
        .unwrap();
    let origin = placement.split(['(', ',']).nth(2).unwrap();
    let original = text
        .lines()
        .find(|s| s.starts_with(&format!("{origin}=")))
        .unwrap();
    let moved = text
        .replace(
            original,
            &format!("{origin}=CARTESIAN_POINT('',(0.,0.000001,0.));"),
        )
        .into_bytes();
    assert_ne!(moved, bytes);
    let import = |bytes: &[u8]| {
        catallax_import::numerical::import_step(
            bytes,
            &sha256_hex(bytes),
            &StepId::new("import").unwrap(),
            &NumericalOptions {
                tolerance_mm: q("0.0000002"),
                vertex_tolerance_mm: None,
            },
            &mut Decider::default(),
        )
    };
    let catallax_import::numerical::StepImport::Numerical(trusted) = import(&moved).unwrap() else {
        panic!("the moved quarter is a numerical body");
    };
    for entity in &trusted.bodies()[0].1.tolerances().entities {
        assert_eq!(entity.tolerance_mm, recorded);
    }
    let refusal = import(&unstamped(&moved)).unwrap_err();
    assert!(format!("{refusal:?}").contains("tolerance"), "{refusal:?}");
}

/// Every tolerance of a body reopened from a file, which must all be the one
/// the file records.
fn reopened_at(bytes: &[u8], recorded: &Rational, what: &str) {
    let read = recorded_tolerance(bytes).unwrap_or_else(|| panic!("{what}: no record"));
    assert_eq!(&read.tolerance_mm, recorded, "{what}");
    assert_eq!(read.field, RecordedField::Recorded, "{what}");
    let (reopened, _) = import(bytes);
    let entities = &reopened.tolerances().entities;
    assert!(!entities.is_empty(), "{what}");
    for entity in entities {
        assert_eq!(&entity.tolerance_mm, recorded, "{what}: {entity:?}");
    }
}

/// The stamp names the body and the receiver, and both names are a caller's
/// text. Here both look like the writer's own fields, a scale among them.
/// The file is not scaled, and it reopens with the tolerance its export
/// reports: the record is a string of its own, and the scale is the file's
/// length unit.
#[test]
fn names_that_look_like_fields_do_not_change_what_an_unscaled_file_records() {
    let name = "part;scaled_export=x1000;written_values=millimetres;rounding_error_mm=3;recorded_tolerance_mm=5;source_tolerance_mm=7";
    let body = gapped_quarter_named(name);
    let (_, plain) = written(&gapped_quarter());
    let receivers = [
        None,
        Some(ReceiverProfile {
            name: "NX;recorded_tolerance_mm=9;scaled_export=x10;written_values=millimetres".into(),
            tolerance_mm: q("0.01"),
        }),
    ];
    for receiver in receivers {
        let what = format!(
            "receiver {:?}",
            receiver.as_ref().map(|profile| &profile.name)
        );
        let options = NumericalStepOptions {
            receiver,
            ..Default::default()
        };
        let file = numerical::build(&body, &options, &mut Decider::default()).unwrap();
        let text = std::str::from_utf8(file.bytes()).unwrap();
        assert!(text.contains("SI_UNIT(.MILLI.,.METRE.)"), "{what}");
        assert!(text.contains(name), "{what}");
        // The name changes nothing the file records.
        let recorded = file.report().uncertainty_mm.clone();
        assert_eq!(recorded, plain, "{what}");
        reopened_at(file.bytes(), &recorded, &what);
    }
}

/// A scaled export writes the model's millimetre values under a larger
/// length unit, and its record is one of those values. The file reopens
/// with every length, and the record, that many times larger. The body's
/// name claims another scale, and it is not read.
#[test]
fn a_scaled_file_is_reopened_at_the_scale_of_its_length_unit() {
    let body = gapped_quarter_named("part;scaled_export=x10;written_values=millimetres");
    let (_, plain) = written(&body);
    for scale in LengthScale::ALL {
        let what = format!("{scale:?}");
        let options = NumericalStepOptions {
            scale: Some(scale),
            ..Default::default()
        };
        let file = numerical::build(&body, &options, &mut Decider::default()).unwrap();
        // The written record is the unscaled file's: the model's value.
        assert_eq!(file.report().uncertainty_mm, plain, "{what}");
        let recorded = &plain * &Rational::from_i64(scale.factor());
        reopened_at(file.bytes(), &recorded, &what);
    }

    // The native tolerant path, at a thousand times.
    let tube = spline_topped_tube();
    let (_, plain) = written(&tube);
    let options = NumericalStepOptions {
        scale: Some(LengthScale::Thousand),
        ..Default::default()
    };
    let file = numerical::build(&tube, &options, &mut Decider::default()).unwrap();
    let recorded = &plain * &Rational::from_i64(1000);
    let (opened, _) = open(file.bytes());
    let import = opened
        .tolerant_admission()
        .expect("the native admission ran");
    assert_eq!(import.bodies().len(), 1, "{:?}", opened.first_refusal());
    assert_eq!(import.bodies()[0].recorded_tolerance_mm(), Some(&recorded));
}

/// SolidDraft's STEP export of the owner's pocket rim rounded at R2, written
/// before check once: 19 faces, torus strips and fitted spline cuts. Its
/// stamp holds the limit its written values were admitted under.
#[test]
fn the_r2_rim_file_is_reopened_without_measuring_it_again() {
    let bytes =
        include_bytes!("../../catallax-mesh/tests/fixtures/reopened-rolling-ball-rim-r2.step");
    let read = recorded_tolerance(bytes).unwrap();
    assert_eq!(read.tolerance_mm, q("1/10000000"));
    assert_eq!(read.field, RecordedField::AdmittedLimit);

    let (opened, work) = open(bytes);
    let Some(BodyView::Tolerant {
        import,
        body: index,
    }) = opened.body(0)
    else {
        panic!("the rim opens as a native tolerant body");
    };
    let solid = &import.bodies()[index];
    assert_eq!(solid.recorded_tolerance_mm(), Some(&read.tolerance_mm));
    assert_eq!(solid.counts().faces, 19);
    let evidence = import.evidence();
    assert!(
        evidence
            .sampled_checks
            .iter()
            .all(|check| check.check.contains("recorded-tolerance")),
    );
    assert!(evidence.repairs.is_empty());
    assert!(
        evidence
            .tolerances
            .iter()
            .all(|entity| entity.tolerance_mm.to_bits() == 1e-7f64.to_bits())
    );
    // Measured, this file's native admission took 1,082,165,753 work units
    // and the whole reopening 1,085,655,229 (main at 945f2504a).
    println!(
        "R2 rim: reopened in {} work, native admission {}",
        work.used,
        opened.native_admission_work_used()
    );
    assert!(
        opened.native_admission_work_used() < 20_000_000,
        "{}",
        opened.native_admission_work_used()
    );
    assert!(work.used < 20_000_000, "{}", work.used);
}
