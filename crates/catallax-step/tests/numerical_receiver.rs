//! Rule T8 (Checkpoint 2 decisions, section 4.3): export with a receiver
//! profile. Each test writes a file and reads what the writer did.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;
use catallax_numbers::Decider;
use catallax_shapes::Surface;
use catallax_step::numerical::{self, LengthScale, NumericalStepOptions, ReceiverProfile};
use catallax_topology::numerical::{
    self as n, NumericalOptions,
    tolerance::{EntityKind, FeatureBand},
};
use fixture::{bowed_copy, q, quarter, scaled};

fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
    }
}

fn text(file: &numerical::NumericalStep) -> String {
    String::from_utf8(file.bytes().to_vec()).unwrap()
}

/// The file's global uncertainty, as written.
fn uncertainty(file: &str) -> String {
    let start = file
        .find("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(")
        .expect("an uncertainty")
        + "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(".len();
    file[start..].split(')').next().unwrap().to_owned()
}

/// A STEP decimal as a rational: STEP allows `1.E-7` and `2.`.
fn decimal(token: &str) -> catallax_numbers::Rational {
    let mut token = token.replace(".E", ".0E");
    if token.ends_with('.') {
        token.push('0');
    }
    q(&token)
}

/// A face copy over the receiver's tolerance is refitted before writing, and
/// the file's uncertainty is the largest entity tolerance of what is
/// written.
#[test]
fn a_receiver_profile_refits_copies_over_its_tolerance_and_writes_the_largest() {
    let mut dec = Decider::default();
    let source = quarter();
    let (bowed, edge) = bowed_copy(&source, "0.006");
    let body = n::build_with_native_trims(source, vec![bowed], options(), &mut dec).unwrap();
    let loose = body
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    assert!(loose > q("0.004"), "{loose}");

    // Without a profile the body is written as it is, and the header has no
    // uncertainty, as before.
    let plain = numerical::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
    assert!(plain.report().receiver.is_none());
    assert!(plain.report().uncertainty_mm >= loose);
    assert!(!text(&plain).contains("UNCERTAINTY_MEASURE_WITH_UNIT"));

    let catia = NumericalStepOptions {
        receiver: Some(ReceiverProfile::catia()),
        ..Default::default()
    };
    let file = numerical::build(&body, &catia, &mut dec).unwrap();
    let report = file.report().receiver.as_ref().unwrap();
    assert_eq!(report.profile, ReceiverProfile::catia());
    assert_eq!(report.refits.len(), 1, "{:?}", report.refits);
    let refit = &report.refits[0];
    assert_eq!(refit.edge, edge);
    assert!(refit.before_mm > q("0.001") && refit.after_mm <= q("0.001"));
    assert!(report.residuals.is_empty(), "{:?}", report.residuals);
    assert!(report.features.is_empty());
    assert_eq!(report.scaled_export_offer, None);
    assert!(file.report().uncertainty_mm <= q("0.001"));
    let written = text(&file);
    assert!(decimal(&uncertainty(&written)) >= file.report().uncertainty_mm);
    assert!(written.contains("receiver=CATIA;receiver_tolerance_mm=1/1000"));
    file.check_written(file.bytes(), &mut dec).unwrap();
}

/// A gap in the geometry survives the refit: each wall edge of a quarter
/// cylinder whose wall sits 0.005 mm off its edges stays over CATIA's
/// 0.001 mm and is listed as a residual, and the file declares it.
#[test]
fn a_gap_no_refit_can_lower_is_listed_and_declared() {
    let mut dec = Decider::default();
    let mut source = quarter();
    let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q("1.005") * &q("1.005")).into();
    let body = n::build(source, options(), &mut dec).unwrap();
    let catia = NumericalStepOptions {
        receiver: Some(ReceiverProfile::catia()),
        ..Default::default()
    };
    let file = numerical::build(&body, &catia, &mut dec).unwrap();
    let report = file.report().receiver.as_ref().unwrap();
    assert!(report.refits.is_empty(), "{:?}", report.refits);
    assert!(!report.residuals.is_empty());
    for residual in &report.residuals {
        assert!(residual.tolerance_mm > q("0.001"));
    }
    assert!(file.report().uncertainty_mm >= q("0.0055"));
    assert!(decimal(&uncertainty(&text(&file))) >= file.report().uncertainty_mm);
}

/// A micro part: the quarter cylinder at a hundredth of its size, 0.01 mm
/// across. For NX at 0.01 mm its features are under 4 times the receiver's
/// tolerance; the writer lists them and offers the smallest scaled export
/// that lifts them to 10 times it. Written at a thousand times,
/// the numbers stay and the length unit becomes the metre.
#[test]
fn a_micro_part_lists_its_small_features_and_is_offered_a_scaled_export() {
    let mut dec = Decider::default();
    let body = n::build(scaled(quarter(), &q("0.01")), options(), &mut dec).unwrap();
    let nx = NumericalStepOptions {
        receiver: Some(ReceiverProfile::nx()),
        ..Default::default()
    };
    let file = numerical::build(&body, &nx, &mut dec).unwrap();
    let report = file.report().receiver.as_ref().unwrap();
    assert!(
        report
            .features
            .iter()
            .any(|feature| feature.band == FeatureBand::UnderFour)
    );
    for feature in &report.features {
        assert!(feature.pair.distance_mm < q("0.1"));
    }
    // The offer is the smallest scale that lifts the smallest feature to
    // at least 10 times the receiver's tolerance.
    let smallest = report
        .features
        .iter()
        .map(|feature| feature.pair.distance_mm.clone())
        .min()
        .unwrap();
    let offer = report.scaled_export_offer.unwrap();
    let lifted =
        |scale: LengthScale| &smallest * &catallax_numbers::Rational::from_i64(scale.factor());
    assert!(lifted(offer) >= q("0.1"), "{smallest} x {offer:?}");
    for smaller in LengthScale::ALL.into_iter().filter(|scale| *scale < offer) {
        assert!(lifted(smaller) < q("0.1"));
    }
    assert!(text(&file).contains("SI_UNIT(.MILLI.,.METRE.)"));

    let metres = NumericalStepOptions {
        receiver: Some(ReceiverProfile::nx()),
        scale: Some(LengthScale::Thousand),
        ..Default::default()
    };
    let scaled_file = numerical::build(&body, &metres, &mut dec).unwrap();
    let written = text(&scaled_file);
    assert!(written.contains("SI_UNIT($,.METRE.)"));
    assert!(!written.contains(".MILLI."));
    assert!(written.contains("scaled_export=x1000;written_values=millimetres"));
    assert_eq!(scaled_file.report().scale, Some(LengthScale::Thousand));
    // The same numbers: the two files differ only in the unit and stamp.
    assert_eq!(
        uncertainty(&written),
        uncertainty(&text(&file)),
        "the uncertainty is written in the same numbers"
    );
}
