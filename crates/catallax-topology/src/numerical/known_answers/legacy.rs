//! Baseline: the sample locations of today's admission, read against the
//! known answers.
//!
//! Admission has two rules. Ordinary trims and the edges of declared numerical
//! carriers are sampled at nine equal locations. Retained native u/v curves
//! are sampled at the eighths of every clipped knot span of the two curves.
//! Both rules below take their locations from the admission code itself and
//! evaluate the gap as admission does, as a 128-bit upper bound. They store
//! no tolerance: admission compares each sample with a fixed limit.

use super::*;

fn edge(case: &EdgeCase) -> NumericalEdge {
    NumericalEdge {
        key: case.name.clone(),
        curve: case.curve.clone(),
        vertices: [0, 1],
        parameters: case.parameters.clone(),
    }
}

fn largest(
    case: &EdgeCase,
    parameters: Vec<Rational>,
    dec: &mut Decider,
) -> Result<Reading, Refusal> {
    let mut reading = Reading {
        raw_mm: Rational::zero(),
        allowance_mm: Rational::zero(),
        stored_mm: None,
        at: None,
        exact_zero: false,
        domain_exit: None,
        evaluations: parameters.len() as u64,
    };
    for t in parameters {
        let [_, upper] = exact::gap_enclosure(case, &t, dec)?;
        if reading.at.is_none() || upper > reading.raw_mm {
            reading.raw_mm = upper;
            reading.at = Some(t);
        }
    }
    Ok(reading)
}

/// Nine equal locations per edge use (`SUBDIVISIONS` in `numerical.rs`).
struct NinePoints;

impl GapChecker for NinePoints {
    fn name(&self) -> &'static str {
        "legacy-nine-points"
    }
    fn applies(&self, _: &EdgeCase) -> bool {
        true
    }
    fn measure(&self, case: &EdgeCase, dec: &mut Decider) -> Result<Reading, Refusal> {
        let parameters = native::numerical_parameters(&edge(case), Sense::Same, dec)?;
        largest(case, parameters, dec)
    }
}

/// Eighths of every clipped knot span of the 3D curve and the native u/v
/// curve (`native::parameters`).
struct KnotSpanEighths;

impl GapChecker for KnotSpanEighths {
    fn name(&self) -> &'static str {
        "legacy-knot-span-eighths"
    }
    fn applies(&self, case: &EdgeCase) -> bool {
        matches!(case.face_copy, FaceCopy::Native(_))
    }
    fn measure(&self, case: &EdgeCase, dec: &mut Decider) -> Result<Reading, Refusal> {
        let FaceCopy::Native(curve) = &case.face_copy else {
            unreachable!("applies to native face copies")
        };
        let parameters = native::parameters(&edge(case), curve, Sense::Same, dec)?;
        largest(case, parameters, dec)
    }
}

fn outcome(rows: &[Row]) -> Vec<(String, bool)> {
    rows.iter()
        .map(|row| (row.name.clone(), row.verdict.is_ok()))
        .collect()
}

/// The baseline reading is the admission code's own: on an analytic face the
/// values equal `native::samples`, sample for sample.
#[test]
fn the_baseline_reads_what_admission_reads() {
    let case = fixed::cases()
        .into_iter()
        .find(|case| case.name == "07a-seam-crossing")
        .expect("case");
    let (FaceCopy::Native(curve), CaseSurface::Kernel(surface)) = (&case.face_copy, &case.surface)
    else {
        unreachable!("an analytic face with a native curve")
    };
    let mut dec = Decider::default();
    let edge = edge(&case);
    let ends = case.parameters.each_ref().map(|t| {
        let at = exact::face_parameters(&case, t, &mut dec).expect("inside");
        at.map(|value| value.as_rational().expect("rational end"))
    });
    let face = NumericalFace {
        key: "face".into(),
        surface: surface.clone(),
        sense: Sense::Same,
        loops: Vec::new(),
    };
    let trim = NumericalTrim {
        edge: 0,
        sense: Sense::Same,
        uv: ends,
    };
    let (samples, _) =
        native::samples(&face, &trim, &edge, curve, false, &mut dec).expect("admission samples");
    let parameters = native::parameters(&edge, curve, Sense::Same, &mut dec).expect("parameters");
    assert_eq!(samples.len(), parameters.len());
    for (sample, t) in samples.iter().zip(&parameters) {
        let [_, upper] = exact::gap_enclosure(&case, t, &mut dec).expect("gap");
        assert_eq!(sample.parameters[0], *t);
        assert_eq!(sample.distance_upper_mm, upper);
    }
}

/// Today's sampling against the fourteen edges. The expected column records
/// which edges each rule already reads within the known-answer rule, so the
/// misses are stated and a change in the oracle or in admission's sample
/// locations shows up here.
#[test]
fn todays_sampling_misses_the_recorded_failures() {
    let cases = fixed::cases();
    let nine = run(&NinePoints, &cases);
    let eighths = run(&KnotSpanEighths, &cases);
    publish(
        "known-answers-legacy-nine-points.txt",
        &receipt(NinePoints.name(), None, &nine),
    );
    publish(
        "known-answers-legacy-knot-span-eighths.txt",
        &receipt(KnotSpanEighths.name(), None, &eighths),
    );
    let passing = |rows: &[Row]| -> Vec<String> {
        outcome(rows)
            .into_iter()
            .filter(|(_, pass)| *pass)
            .map(|(name, _)| name)
            .collect()
    };
    assert_eq!(passing(&nine), NINE_POINTS_PASS, "nine-point rule");
    assert_eq!(passing(&eighths), KNOT_SPAN_EIGHTHS_PASS, "knot-span rule");
}

/// Edges the nine-point rule reads within the known-answer rule.
const NINE_POINTS_PASS: &[&str] = &[
    "13a-zero-line-on-plane",
    "13b-zero-circle-on-cylinder",
    "13c-zero-boundary-rail",
    "13d-zero-stored-iso-curve",
];

/// Edges the knot-span rule reads within the known-answer rule.
const KNOT_SPAN_EIGHTHS_PASS: &[&str] = &[
    "01-narrow-peak-dense-knots",
    "06-peak-midway-between-equal-samples",
    "13c-zero-boundary-rail",
    "13d-zero-stored-iso-curve",
];

/// Today's sampling against new seeded edges. No outcome is expected of it;
/// the receipt records what it reads, with the seed.
#[test]
fn todays_sampling_on_fresh_edges_is_recorded() {
    let seed = fresh::seed();
    let cases = fresh::cases(seed, 3);
    for (checker, file) in [
        (
            &NinePoints as &dyn GapChecker,
            "known-answers-fresh-legacy-nine-points.txt",
        ),
        (
            &KnotSpanEighths,
            "known-answers-fresh-legacy-knot-span-eighths.txt",
        ),
    ] {
        let rows = run(checker, &cases);
        assert_eq!(rows.len(), cases.len());
        publish(file, &receipt(checker.name(), Some(seed), &rows));
    }
}
