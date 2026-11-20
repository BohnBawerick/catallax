//! Rule T8 (Checkpoint 2 decisions, section 4.3): export with a receiver
//! profile. Each test writes a file and reads what the writer did.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;
#[path = "../../catallax-topology/tests/common/tent.rs"]
mod tent;
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
        vertex_tolerance_mm: None,
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

    let scaled = numerical::build(
        &body,
        &NumericalStepOptions {
            scale: Some(LengthScale::Ten),
            ..catia
        },
        &mut dec,
    )
    .unwrap();
    let report = scaled.report().receiver.as_ref().unwrap();
    assert_eq!(report.refits.len(), 1, "{:?}", report.refits);
    assert_eq!(report.refits[0].before_mm, &refit.before_mm * &q("10"));
    assert!(report.refits[0].after_mm <= report.profile.tolerance_mm);
    assert!(report.residuals.is_empty());
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
        assert!(residual.entity.tolerance_mm > q("0.001"));
        assert_eq!(residual.receiver_tolerance_mm, residual.entity.tolerance_mm);
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

/// A face's interior is a feature a receiver sees (reviewer hold on
/// 83858ed4): the apex of a tent 0.008 mm above the middle of its 100 mm
/// base, with every edge and vertex of the base more than 70 mm away. For
/// CATIA at 0.001 mm that is under 10 times and not under 4 times its
/// tolerance; for NX at 0.01 mm it is under 4 times; 0.0115 mm above, it is
/// no feature for CATIA.
#[test]
fn a_vertex_near_a_face_interior_is_a_receiver_feature() {
    let mut dec = Decider::default();
    let band_of = |height: &str, profile: ReceiverProfile, dec: &mut Decider| {
        let body = n::build(tent::tent(height, "0.002"), options(), dec).unwrap();
        let file = numerical::build(
            &body,
            &NumericalStepOptions {
                receiver: Some(profile),
                ..Default::default()
            },
            dec,
        )
        .unwrap();
        let report = file.report().receiver.clone().unwrap();
        report
            .features
            .iter()
            .find(|feature| {
                let ends = [&feature.pair.a, &feature.pair.b];
                ends.contains(&&(EntityKind::Vertex, tent::APEX.to_string()))
                    && ends.contains(&&(EntityKind::Face, tent::BASE.to_string()))
            })
            .map(|feature| feature.band)
    };
    assert_eq!(
        band_of("0.006", ReceiverProfile::catia(), &mut dec),
        Some(FeatureBand::UnderTen)
    );
    assert_eq!(
        band_of("0.006", ReceiverProfile::nx(), &mut dec),
        Some(FeatureBand::UnderFour)
    );
    assert_eq!(band_of("0.0095", ReceiverProfile::catia(), &mut dec), None);
}

/// A scaled export is read by its receiver at the scale (reviewer hold on
/// 83858ed4). Written a thousand times larger for CATIA at 0.001 mm, a
/// tolerance of 2e-6 mm reads as 0.002 mm and is a residual, though it is
/// none unscaled; and the tent's apex 0.008 mm over its base reads as 8 mm,
/// no small feature, though it is one unscaled.
#[test]
fn a_scaled_export_applies_the_receivers_limits_to_the_scaled_file() {
    let mut dec = Decider::default();
    let catia = |scale: Option<LengthScale>| NumericalStepOptions {
        receiver: Some(ReceiverProfile::catia()),
        scale,
        ..Default::default()
    };
    // The wall 1.8e-6 mm off its edges: about 2e-6 mm on each wall edge.
    let mut source = quarter();
    let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q("1.0000018") * &q("1.0000018")).into();
    let body = n::build(source, options(), &mut dec).unwrap();
    let plain = numerical::build(&body, &catia(None), &mut dec).unwrap();
    let plain = plain.report().receiver.as_ref().unwrap();
    assert_eq!(plain.scale, 1);
    assert!(plain.residuals.is_empty(), "{:?}", plain.residuals);
    let scaled = numerical::build(&body, &catia(Some(LengthScale::Thousand)), &mut dec).unwrap();
    let scaled = scaled.report().receiver.as_ref().unwrap();
    assert_eq!(scaled.scale, 1000);
    assert!(!scaled.residuals.is_empty());
    for residual in &scaled.residuals {
        assert!(residual.receiver_tolerance_mm > q("0.001"), "{residual:?}");
        assert_eq!(
            residual.receiver_tolerance_mm,
            &residual.entity.tolerance_mm * &q("1000")
        );
    }
    // The tent: a feature for CATIA unscaled, none at a thousand times.
    let tent = n::build(tent::tent("0.006", "0.002"), options(), &mut dec).unwrap();
    let features = |scale| {
        numerical::build(&tent, &catia(scale), &mut Decider::default())
            .unwrap()
            .report()
            .receiver
            .clone()
            .unwrap()
            .features
    };
    let unscaled = features(None);
    assert!(!unscaled.is_empty());
    for feature in &unscaled {
        assert_eq!(feature.receiver_distance_mm, feature.pair.distance_mm);
    }
    assert!(features(Some(LengthScale::Thousand)).is_empty());
    // At ten times, 0.08 mm reads as no small feature for CATIA either;
    // a band is only ever given on the receiver's distance.
    for feature in features(Some(LengthScale::Ten)) {
        assert_eq!(
            feature.receiver_distance_mm,
            &feature.pair.distance_mm * &q("10")
        );
        assert!(feature.receiver_distance_mm < q("0.01"));
    }
}
