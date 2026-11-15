//! Actual Part 21 circular edges keep shared identity through numerical admission.
mod common;
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_import::{
    numerical::{self, StepImport},
    sha256_hex,
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::{q, quarter};

fn references(line: &str) -> Vec<u64> {
    line.split_once('=')
        .unwrap()
        .1
        .split('#')
        .skip(1)
        .map(|s| {
            s.chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>()
                .parse()
                .unwrap()
        })
        .collect()
}
fn displaced_wall(delta: &str) -> Vec<u8> {
    let body = n::build(
        quarter(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let file = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let text = std::str::from_utf8(file.bytes()).unwrap();
    let plane = text.lines().find(|s| s.contains("=PLANE('f2',")).unwrap();
    let placement = references(plane)[0];
    let placement = text
        .lines()
        .find(|s| s.starts_with(&format!("#{placement}=")))
        .unwrap();
    let origin = references(placement)[0];
    let original = text
        .lines()
        .find(|s| s.starts_with(&format!("#{origin}=")))
        .unwrap();
    text.replace(
        original,
        &format!("#{origin}=CARTESIAN_POINT('',(0.,{delta},0.));"),
    )
    .into_bytes()
}
fn checked(bytes: &[u8], tolerance: &str) -> Result<StepImport, numerical::ImportError> {
    numerical::import_step(
        bytes,
        &sha256_hex(bytes),
        &StepId::new("analytic-source").unwrap(),
        &n::NumericalOptions {
            tolerance_mm: q(tolerance),
        },
        &mut Decider::default(),
    )
}

fn combine_sources(first: &[u8], second: &str) -> Vec<u8> {
    let first = std::str::from_utf8(first).unwrap();
    let split = first.rfind("ENDSEC;").unwrap();
    let data = second
        .split_once("DATA;\n")
        .unwrap()
        .1
        .split_once("ENDSEC;")
        .unwrap()
        .0;
    let offset = first
        .lines()
        .filter_map(|line| line.strip_prefix('#'))
        .map(|line| line.split_once('=').unwrap().0.parse::<u64>().unwrap())
        .max()
        .unwrap();
    let mut shifted = String::new();
    let mut chars = data.chars().peekable();
    while let Some(ch) = chars.next() {
        shifted.push(ch);
        if ch == '#' {
            let mut digits = String::new();
            while chars.peek().is_some_and(char::is_ascii_digit) {
                digits.push(chars.next().unwrap());
            }
            shifted.push_str(&(digits.parse::<u64>().unwrap() + offset).to_string());
        }
    }
    format!("{}{}{}", &first[..split], shifted, &first[split..]).into_bytes()
}

#[test]
fn an_incidence_gap_cannot_hide_a_crossing_in_another_solid_in_either_order() {
    let gapped = displaced_wall("0.000001");
    let mut soup = common::Soup::new();
    soup.box_(["0", "0", "0"], ["2", "2", "2"]);
    soup.box_(["1", "1", "1"], ["3", "3", "3"]);
    let crossing = soup.finish();
    let combined = combine_sources(&gapped, &crossing);
    let text = std::str::from_utf8(&combined).unwrap();
    let shells = text
        .lines()
        .filter(|line| line.contains("=CLOSED_SHELL("))
        .collect::<Vec<_>>();
    let face_refs = shells
        .iter()
        .flat_map(|line| references(line))
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let name = shells[0].split_once('=').unwrap().0;
    let same_solid = text
        .replace(
            shells[0],
            &format!("{name}=CLOSED_SHELL('',({face_refs}));"),
        )
        .into_bytes();
    for source in [
        combined,
        same_solid,
        combine_sources(crossing.as_bytes(), std::str::from_utf8(&gapped).unwrap()),
    ] {
        let result = checked(&source, "0.01");
        assert!(
            matches!(
                result,
                Err(numerical::ImportError::Source(
                    catallax_import::ImportRefusal::SelfIntersecting { .. }
                ))
            ),
            "{result:?}"
        );
    }
}

#[test]
fn source_pair_coverage_includes_every_unaffected_solid() {
    let source = combine_sources(
        &displaced_wall("0.000001"),
        &common::box_step(common::BoxFixture::new(["2", "3", "4"])),
    );
    let recognized = catallax_import::read_step(&source).unwrap();
    let StepImport::Numerical(imported) = checked(&source, "0.01").unwrap() else {
        panic!("expected numerical bodies")
    };
    assert_eq!(imported.bodies().len(), 2);
    let checks = &imported.source_admission().unwrap().intersection_checks;
    assert_eq!(checks.len(), 10 + 15);
    let untouched = checks
        .iter()
        .filter(|check| {
            matches!(check,
        numerical::NumericalSourcePairCheck::CheckedNoCrossing { solid, .. }
        if solid == &recognized.solids[1].name)
        })
        .count();
    assert_eq!(untouched, 15);
}

#[test]
fn curved_dependency_rebuilds_finite_parameters_and_retains_shared_vertices() {
    let source = displaced_wall("0.000001");
    let before = source.clone();
    let StepImport::Numerical(imported) = checked(&source, "0.01").unwrap() else {
        panic!("expected numerical admission")
    };
    let admission = imported.source_admission().unwrap();
    assert_eq!(admission.method, "retained_shared_analytic_boundaries/3");
    assert_eq!(admission.intersection_checks.len(), 10);
    assert!(admission.intersection_checks.iter().any(|check| matches!(
        check,
        numerical::NumericalSourcePairCheck::CheckedNoCrossing { .. }
    )));
    let recognized = catallax_import::read_step(&source).unwrap();
    let affected = recognized.solids[0].faces.iter().find(|face|
        matches!(&face.surface, catallax_import::Surface::Plane(frame) if frame.origin[1] == q("0.000001")))
        .unwrap().name.clone();
    for check in &admission.intersection_checks {
        if let numerical::NumericalSourcePairCheck::Unresolved { faces, refusal, .. } = check
            && matches!(
                refusal.as_ref(),
                catallax_import::ImportRefusal::NotSolid { .. }
            )
        {
            assert!(faces.contains(&affected));
        }
    }
    assert_eq!(admission.tolerance_mm, q("0.01"));
    assert!(
        matches!(&admission.strict_refusal,catallax_import::ImportRefusal::PlanarRepairNotVerified { detail } if detail.stage=="nearest_feature" && detail.feature_size_squared_mm2.is_none())
    );
    let body = &imported.bodies()[0].1;
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (5, 9, 6));
    assert_eq!(
        body.definition()
            .vertices
            .iter()
            .map(|v| v.position.clone())
            .collect::<std::collections::BTreeSet<_>>(),
        quarter()
            .vertices
            .iter()
            .map(|v| v.position.clone())
            .collect::<std::collections::BTreeSet<_>>()
    );
    assert_eq!(source, before);
    assert_eq!(imported.record().source.sha256, sha256_hex(&source));
    assert!(imported.record().face_tolerances.is_empty());
    assert!(imported.record().derivative.is_none());
    let step = output::build(
        body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(
        std::str::from_utf8(step.bytes())
            .unwrap()
            .contains("numerical approximate")
    );
}

#[test]
fn over_cap_residual_and_a_smaller_requested_tolerance_refuse() {
    assert!(checked(&displaced_wall("0.0101"), "0.01").is_err());
    assert!(checked(&displaced_wall("0.000001"), "0.0000001").is_err());
    assert!(checked(&displaced_wall("0.000001"), "0.0101").is_err());
}

#[test]
fn exhausted_analytic_imports_return_no_body() {
    let source = displaced_wall("0.000001");
    let (result, work) = with_operation_budget(1, || checked(&source, "0.01"));
    assert!(result.is_err() && work.exhausted);
}

#[test]
fn source_bytes_are_charged_before_the_initial_verification_hash() {
    let source = vec![b' '; 100_000];
    let (result, work) = with_operation_budget(1, || {
        numerical::import_step(
            &source,
            &"0".repeat(64),
            &StepId::new("analytic-source").unwrap(),
            &n::NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(
        work.exhausted,
        "verification hashed unmetered bytes: {result:?}"
    );
    assert!(matches!(
        result,
        Err(numerical::ImportError::Source(
            catallax_import::ImportRefusal::ResourceLimit
        ))
    ));
}

#[test]
fn a_caller_deadline_keeps_its_actual_refusal_code() {
    let mut dec = Decider::default();
    dec.set_deadline(Some(
        std::time::Instant::now() - std::time::Duration::from_secs(1),
    ));
    let result = numerical::import_step(
        b"not STEP",
        &"0".repeat(64),
        &StepId::new("expired-source").unwrap(),
        &n::NumericalOptions::default(),
        &mut dec,
    );
    assert!(
        matches!(result,Err(numerical::ImportError::Numerical(ref refusal)) if refusal.code=="numbers.deadline")
    );
}
