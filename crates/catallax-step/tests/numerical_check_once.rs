//! Check once (roadmap item 3.0A step 5; Checkpoint 2 decisions, section 4.1
//! item 5): "Export writes the numbers and adds the rounding error to the
//! tolerance; it does not re-admit the body."
//!
//! Each test writes a body through one of the three writers and reads what
//! the file records. The audit is the admission of the written values that
//! an export ran before: it is the yardstick here, and no export runs it.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;
#[path = "common/n03_tube.rs"]
mod tube;
use catallax_numbers::{
    Decider, NumericalCheck, NumericalEvidence, Rational, with_operation_budget,
};
use catallax_step::numerical::{self, NumericalStepOptions};
use catallax_topology::numerical::{self as n, NumericalBody, NumericalOptions};
use fixture::{bowed_copy, q, quarter};

fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
        vertex_tolerance_mm: None,
    }
}

/// One body for each writer: lines on planes, an analytic quarter cylinder,
/// and the quarter with a retained u/v curve, which the curved writer takes.
fn bodies() -> Vec<(&'static str, NumericalBody)> {
    let mut dec = Decider::default();
    let source = quarter();
    let (bowed, _) = bowed_copy(&source, "0.006");
    vec![
        (
            "planar",
            n::build(tube::tube(), options(), &mut dec).unwrap(),
        ),
        (
            "analytic",
            n::build(quarter(), options(), &mut dec).unwrap(),
        ),
        (
            "curved",
            n::build_with_native_trims(source, vec![bowed], options(), &mut dec).unwrap(),
        ),
    ]
}

fn digits(significant_digits: u32) -> NumericalStepOptions {
    NumericalStepOptions {
        significant_digits,
        // A coarse file's lengths round by more than the default limit.
        max_scalar_rounding_error_mm: q("1"),
        ..Default::default()
    }
}

/// One field of the file's stamp, as it is written.
fn stamp(file: &numerical::NumericalStep, key: &str) -> String {
    let text = std::str::from_utf8(file.bytes()).unwrap();
    let start = text.find(&format!(";{key}=")).expect(key) + key.len() + 2;
    let end = text[start..].find([';', '\'']).unwrap();
    text[start..start + end].to_owned()
}

/// The stored records of the edge uses a set of checks measured.
fn measured_uses(evidence: &NumericalEvidence) -> Vec<Rational> {
    evidence
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, samples, .. }
                if name.starts_with("trim.") && name.ends_with("edge-surface-distance") =>
            {
                samples
                    .iter()
                    .map(|sample| sample.distance_upper_mm.clone())
                    .max()
            }
            _ => None,
        })
        .collect()
}

#[test]
fn an_export_writes_the_numbers_and_does_not_admit_them_again() {
    for (name, body) in bodies() {
        let options = NumericalStepOptions::default();
        let (file, export) = with_operation_budget(20_000_000, || {
            numerical::build(&body, &options, &mut Decider::default())
        });
        let file = file.unwrap();
        let (audit, audited) = with_operation_budget(20_000_000, || {
            numerical::audit_written(&body, &options, &mut Decider::default())
        });
        let audit = audit.unwrap();
        // The audit writes the same bytes and then measures every use of
        // them. The export stops at the bytes.
        assert_eq!(audit.file.bytes(), file.bytes(), "{name}");
        let uses: usize = body
            .definition()
            .faces
            .iter()
            .map(|face| face.loops.iter().map(Vec::len).sum::<usize>())
            .sum();
        assert_eq!(measured_uses(&audit.evidence).len(), uses, "{name}");
        assert!(
            export.used < audited.used,
            "{name}: export {} and audit {}",
            export.used,
            audited.used
        );
        println!(
            "{name}: export {} work, audit of the written values {} work",
            export.used, audited.used
        );
        let checks = &file.report().checks;
        assert!(
            checks.contains(&"file.rounding-error-added-to-tolerance".to_owned()),
            "{name}: {checks:?}"
        );
        assert!(
            checks.iter().all(|check| !check.contains("admission")),
            "{name}: {checks:?}"
        );
    }
}

#[test]
fn the_file_records_the_largest_entity_tolerance_plus_the_rounding_error() {
    for (name, body) in bodies() {
        let largest = body.tolerances().largest_mm();
        for significant_digits in [17, 9, 6] {
            let options = digits(significant_digits);
            let file = numerical::build(&body, &options, &mut Decider::default()).unwrap();
            let report = file.report();
            let exact = &largest + &report.rounding.bound_mm;
            // The recorded tolerance is that sum rounded up to the written
            // digits: no smaller, and larger by under one part in 10^(d-1).
            let slack = &exact * &q(&format!("1e-{}", significant_digits - 1));
            assert!(
                report.uncertainty_mm >= exact && report.uncertainty_mm <= &exact + &slack,
                "{name} at {significant_digits} digits: {} against {exact}",
                report.uncertainty_mm
            );
            // The stamp holds the same number, and the rounding error by it.
            assert_eq!(
                q(&stamp(&file, "recorded_tolerance_mm")),
                report.uncertainty_mm,
                "{name}"
            );
            assert!(
                q(&stamp(&file, "rounding_error_mm")) >= report.rounding.bound_mm,
                "{name}"
            );
            // Seventeen digits move nothing that matters; six move more
            // than the body's own tolerance.
            if significant_digits == 17 {
                assert!(report.rounding.bound_mm < q("1e-9"), "{name}");
            }
        }
    }
}

/// Every gap an audit measured on the written values: each edge use of the
/// written tables, with rule T4's margin taken off its stored record, and
/// each sample of a u/v line as an analytic face's file holds it.
fn measured_gaps(evidence: &NumericalEvidence) -> Vec<Rational> {
    let margin = q("1.1");
    let mut gaps: Vec<Rational> = measured_uses(evidence)
        .into_iter()
        .map(|stored| &stored / &margin)
        .collect();
    for check in evidence.checks() {
        if let NumericalCheck::Samples { name, samples, .. } = check
            && name == "file.actual-decimal-pcurve-edge-surface-distance"
        {
            gaps.extend(
                samples
                    .iter()
                    .map(|sample| sample.distance_upper_mm.clone()),
            );
        }
    }
    gaps
}

#[test]
fn the_values_a_coarse_file_holds_are_inside_what_it_records() {
    for (name, body) in bodies() {
        let own = body.tolerances().largest_mm();
        let mut opened = false;
        for significant_digits in [17, 9, 6, 4] {
            let options = digits(significant_digits);
            let audit = numerical::audit_written(&body, &options, &mut Decider::default())
                .unwrap_or_else(|refusal| panic!("{name} at {significant_digits}: {refusal:?}"));
            let recorded = &audit.file.report().uncertainty_mm;
            let largest = measured_gaps(&audit.evidence)
                .into_iter()
                .max()
                .expect("a measured gap");
            // No gap measured on the written values is over the record.
            assert!(
                largest <= *recorded,
                "{name} at {significant_digits} digits: measured {largest} over recorded {recorded}"
            );
            opened |= largest > own;
            println!(
                "{name} at {significant_digits} digits: largest measured gap {:.3e} mm, recorded {:.3e} mm, of which rounding {:.3e} mm",
                largest.to_f64_nearest(),
                recorded.to_f64_nearest(),
                audit.file.report().rounding.bound_mm.to_f64_nearest()
            );
        }
        // On the analytic body four digits open a gap wider than any the
        // body held, and the rounding error is what covers it. The planar
        // tube's whole numbers round to themselves at any digit count.
        if name == "analytic" {
            assert!(opened, "{name}");
        }
        if name == "planar" {
            assert!(!opened, "{name}");
        }
    }
}

#[test]
fn a_file_whose_numbers_are_exact_records_the_bodys_tolerance_alone() {
    // The square tube's coordinates, directions and magnitudes are whole
    // numbers, so nothing is rounded.
    let body = n::build(tube::tube(), options(), &mut Decider::default()).unwrap();
    let file = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let report = file.report();
    assert_eq!(report.rounding.bound_mm, Rational::zero());
    // Rule T4's 1.1 is in the entity tolerance already and is not applied to
    // the record again.
    assert_eq!(report.uncertainty_mm, body.tolerances().largest_mm());
    assert_eq!(stamp(&file, "rounding_error_mm"), "0");
}
