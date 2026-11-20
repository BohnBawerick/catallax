//! Tolerance rules T0 to T8 (decision Q7, section 4.3 of the Checkpoint 2
//! decisions) on admitted bodies: the cap at its boundary, rule T1's size
//! warning, rule T3's refit before an opt-in, and rule T7 along a chain.
#[path = "common/analytic.rs"]
mod fixture;
use catallax_numbers::{Decider, NumericalCheck, Rational};
use catallax_shapes::Surface;
use catallax_topology::numerical::{
    self as n, NumericalOptions, ToleranceWarning,
    tolerance::{self, EntityKind},
};
use fixture::{bowed_copy, q, quarter, scaled};

/// The quarter cylinder's curved face.
const WALL: usize = 4;

fn options(cap: &str) -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q(cap),
    }
}

fn with_radius(value: &str) -> n::NumericalDefinition {
    let mut source = quarter();
    set_radius(&mut source, value);
    source
}

fn set_radius(source: &mut n::NumericalDefinition, value: &str) {
    let Surface::Cylinder(cylinder) = &mut source.faces[WALL].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q(value) * &q(value)).into();
}

fn wall(value: &str) -> Surface {
    with_radius(value).faces[WALL].surface.clone()
}

/// The largest measured edge-use record.
fn largest_record(body: &n::NumericalBody) -> Rational {
    body.evidence()
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, samples, .. }
                if name == "trim.edge-surface-distance" =>
            {
                Some(samples[0].distance_upper_mm.clone())
            }
            _ => None,
        })
        .fold(Rational::zero(), Rational::max)
}

/// The near-cap boundary (reviewer note relayed in inbox 003): a raw gap
/// under the cap is admitted only when rule T4's full allowance, 1.1 times
/// the gap plus its measurement error, fits under it. Off the known-answer
/// fixtures: the quarter cylinder's wall is widened, so each wall edge sits
/// a known distance off it, measured by the checker's closed forms.
#[test]
fn a_raw_gap_under_the_cap_is_admitted_only_when_its_full_allowance_fits() {
    let cap = q("0.01");
    // 0.009 mm: 1.1 x 0.009 = 0.0099 mm, under the cap with its error.
    let fits = n::build(
        with_radius("1.009"),
        options("0.01"),
        &mut Decider::default(),
    )
    .expect("the full allowance fits");
    let stored = largest_record(&fits);
    assert!(stored <= cap, "{stored}");
    assert!(stored >= q("0.0098999"), "{stored}");
    // 0.0095 mm: the raw gap is under the cap, but 1.1 x 0.0095 = 0.01045
    // mm is not. The existing named refusal, with no threshold changed.
    let raw = q("0.0095");
    assert!(raw < cap);
    let refused = n::build(
        with_radius("1.0095"),
        options("0.01"),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refused.code, "numbers.numerical-tolerance-exceeded");
}

/// Rule T1 on an admitted body: the floor of a body under 10 m is 1e-7 mm,
/// and a body over 1 km is admitted with a warning, its floor raised to
/// 1e-11 of its coordinate size.
#[test]
fn a_body_over_a_kilometre_is_admitted_with_a_warning_and_a_raised_floor() {
    let small = n::build(quarter(), options("0.01"), &mut Decider::default()).unwrap();
    assert_eq!(small.tolerances().floor_mm, q("0.0000001"));
    assert!(small.tolerances().warnings.is_empty());
    for entity in &small.tolerances().entities {
        assert!(entity.tolerance_mm >= small.tolerances().floor_mm);
    }
    let k = q("2000000");
    let big = n::build(
        scaled(quarter(), &k),
        options("0.01"),
        &mut Decider::default(),
    )
    .unwrap();
    let report = big.tolerances();
    // The frame's box holds each circle's centre widened by its radius, so
    // `R` is at least the 1 km of the quarter's own box.
    assert!(report.coordinate_size_mm >= q("1000000"));
    assert_eq!(
        report.floor_mm,
        tolerance::floor_mm(&report.coordinate_size_mm)
    );
    assert!(report.floor_mm >= q("0.00001"));
    assert!(report.warnings.iter().any(|warning| matches!(
        warning,
        ToleranceWarning::LargeBody { extent_mm } if *extent_mm > q("1000000")
    )));
}

/// Rule T7 along a chain of edits that starts at a near-zero input: the wall
/// of the quarter cylinder sits on its edges (each stores rule T1's floor),
/// then moves out in four steps. Each stage reports the absolute tolerance
/// before and after, the cap, and what is left under it. The warning fires
/// exactly where a tolerance rose more than ten times above its input, and
/// stays a warning: nothing is refused while the cap holds. The refit of a
/// grown edge cannot lower a gap that is in the geometry, so the warning of
/// the first step stands, small as its absolute value is.
#[test]
fn a_chain_reports_absolute_tolerances_and_headroom_at_every_stage() {
    let cap = q("0.01");
    let mut dec = Decider::default();
    let mut body = n::build(quarter(), options("0.01"), &mut dec).unwrap();
    let largest_edge = |body: &n::NumericalBody| {
        body.tolerances()
            .entities
            .iter()
            .filter(|entity| entity.kind == EntityKind::Edge)
            .map(|entity| entity.tolerance_mm.clone())
            .fold(Rational::zero(), Rational::max)
    };
    let start = largest_edge(&body);
    println!("stage 0: largest edge {start} mm (near zero, at the floor), cap {cap} mm");
    assert_eq!(start, q("0.0000001"));
    let steps = [
        ("1.000001", true),
        ("1.000003", false),
        ("1.0005", true),
        ("1.0009", false),
    ];
    for (stage, (radius, warns)) in steps.into_iter().enumerate() {
        let input = largest_edge(&body);
        body = body.replace_surface("f4", wall(radius), &mut dec).unwrap();
        let output = largest_edge(&body);
        let report = body.tolerances();
        let growth: Vec<_> = report
            .warnings
            .iter()
            .filter_map(|warning| match warning {
                ToleranceWarning::Growth {
                    kind: EntityKind::Edge,
                    input_mm,
                    output_mm,
                    cap_mm,
                    headroom_mm,
                    same_entity,
                    ..
                } => Some((
                    input_mm.clone(),
                    output_mm.clone(),
                    cap_mm.clone(),
                    headroom_mm.clone(),
                    *same_entity,
                )),
                _ => None,
            })
            .collect();
        println!(
            "stage {}: radius {radius}: largest edge {} mm -> {} mm (x{:.2}), cap {} mm, headroom {} mm, {} growth warnings",
            stage + 1,
            input.to_f64_nearest(),
            output.to_f64_nearest(),
            (&output / &input).to_f64_nearest(),
            report.cap_mm,
            report.headroom_mm().to_f64_nearest(),
            growth.len()
        );
        assert_eq!(report.cap_mm, cap);
        assert!(output <= cap);
        assert_eq!(report.headroom_mm(), &cap - &report.largest_mm());
        assert_eq!(!growth.is_empty(), warns, "stage {}", stage + 1);
        for (input_mm, output_mm, cap_mm, headroom_mm, same_entity) in growth {
            assert!(same_entity);
            assert!(tolerance::grew(&input_mm, &output_mm));
            assert_eq!(headroom_mm, &cap_mm - &output_mm);
            assert!(headroom_mm > Rational::zero());
        }
    }
}

/// Rule T7's refit: a continuation gives one circle use of the wall a face
/// copy bowed off its 3D curve (the curve itself still lies on the wall),
/// which raises that edge's tolerance more than ten times. The face copy is
/// refitted from the 3D curve, measured again, and the lower tolerance is
/// kept, so no growth is left to warn about.
#[test]
fn a_grown_face_copy_is_refitted_and_its_tolerance_lowered() {
    let mut dec = Decider::default();
    let body = n::build(quarter(), options("0.01"), &mut dec).unwrap();
    let (bowed, edge) = bowed_copy(body.definition(), "0.006");
    let before = body
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    // Without the refit the bowed copy stores about 1.1 x 0.0045 mm.
    let unrefitted = n::build_with_native_trims(
        body.definition().clone(),
        vec![bowed.clone()],
        options("0.01"),
        &mut dec,
    )
    .unwrap();
    let bowed_mm = unrefitted
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    assert!(bowed_mm > q("0.004"), "{bowed_mm}");
    let next = body
        .rebuild_with_native_trims(
            body.definition().clone(),
            vec![bowed],
            options("0.01"),
            &mut dec,
        )
        .unwrap();
    let after = next
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    println!("edge {edge}: {before} mm before, {bowed_mm} mm bowed, {after} mm after the refit");
    assert!(after < q("0.000001"), "{after}");
    assert!(!next.tolerances().warnings.iter().any(|warning| matches!(
        warning,
        ToleranceWarning::Growth { key, .. } if *key == edge
    )));
}

/// Rule T3: an import whose face copy puts a stored tolerance over the
/// automatic 0.01 mm cap is refitted first; where the refit brings it under,
/// the body is admitted at the automatic cap. A gap in the geometry itself
/// survives the refit and refuses, naming the opt-in; an opt-in up to 0.1
/// mm admits it, and one over 0.1 mm refuses.
#[test]
fn an_import_is_refitted_before_the_user_is_asked_to_opt_in() {
    let mut dec = Decider::default();
    let source = quarter();
    let (bowed, _) = bowed_copy(&source, "0.03");
    assert_eq!(
        n::build_with_native_trims(
            source.clone(),
            vec![bowed.clone()],
            options("0.01"),
            &mut dec
        )
        .unwrap_err()
        .code,
        "numbers.numerical-tolerance-exceeded"
    );
    let admitted = n::admit_import(source, vec![bowed], options("0.01"), &mut dec).unwrap();
    assert!(admitted.tolerances().largest_mm() <= q("0.01"));
    assert_eq!(admitted.options().tolerance_mm, q("0.01"));

    let refused =
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(refused.code, "numbers.numerical-tolerance-exceeded");
    assert!(refused.detail.contains("opt in"), "{}", refused.detail);
    let opted =
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.05"), &mut dec).unwrap();
    assert!(opted.tolerances().largest_mm() <= q("0.05"));
    assert_eq!(
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.2"), &mut dec)
            .unwrap_err()
            .code,
        "numerical.import-opt-in"
    );
}
