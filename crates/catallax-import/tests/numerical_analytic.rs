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
    // A wall moved by hand is not what the export recorded. The stamp goes
    // with it, so the file is measured as any import is; with the stamp,
    // reopening would trust the tolerance it records (tests/check_once.rs).
    text.replace(
        original,
        &format!("#{origin}=CARTESIAN_POINT('',(0.,{delta},0.));"),
    )
    .replacen(
        "catallax numerical approximate;",
        "another system approximate;",
        1,
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
            vertex_tolerance_mm: None,
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

#[test]
fn a_gapped_cone_retains_its_fitted_carrier_and_records_both_tolerances() {
    let mut soup = common::Soup::new();
    soup.cone(["0", "0", "0"], "1", "30.", Ok(("3", "2.732050808")));
    let source = common::in_degrees(&soup.finish());
    let plane = source.lines().find(|s| s.contains("=PLANE(")).unwrap();
    let placement = references(plane)[0];
    let placement = source
        .lines()
        .find(|s| s.starts_with(&format!("#{placement}=")))
        .unwrap();
    let origin = references(placement)[0];
    let original = source
        .lines()
        .find(|s| s.starts_with(&format!("#{origin}=")))
        .unwrap();
    let gapped = source.replace(
        original,
        &format!("#{origin}=CARTESIAN_POINT('',(0.,0.,0.000001));"),
    );
    let (result, work) = with_operation_budget(20_000_000, || checked(gapped.as_bytes(), "0.01"));
    assert!(!work.exhausted, "{work:?}");
    let StepImport::Numerical(imported) = result.unwrap() else {
        panic!("expected marked numerical continuation")
    };
    assert_eq!(imported.bodies().len(), 1);
    let admission = imported.source_admission().unwrap();
    assert_eq!(admission.tolerance_mm, q("0.01"));
    assert_eq!(admission.topology_repairs.len(), 1);
    assert_eq!(
        admission.topology_repairs[0].method,
        "paired_cone_generatrix/1"
    );
    let tolerances = &imported.record().face_tolerances;
    assert_eq!(tolerances.len(), 1);
    assert!(q(&tolerances[0].deviation_upper_mm) < q("1e-9"));
    let body = &imported.bodies()[0].1;
    assert!(&body.options().tolerance_mm + &q(&tolerances[0].deviation_upper_mm) <= q("0.01"));
    let mut reaches_apex = body.definition().clone();
    for face in &mut reaches_apex.faces {
        if matches!(face.surface, catallax_shapes::Surface::Cone(_)) {
            for uv in face.loops.iter_mut().flatten().flat_map(|t| &mut t.uv) {
                uv[1] = q("0");
            }
        }
    }
    let refusal = n::build(
        reaches_apex,
        body.options().clone(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "numerical.invalid-trim");
    assert!(refusal.detail.contains("above its apex"));
}

fn cone_export_source(angle: &str, top: &str, mode: &str) -> Vec<u8> {
    let mut soup = common::Soup::new();
    soup.cone(["0", "0", "0"], "1", angle, Ok(("3", top)));
    let source = common::in_degrees(&soup.finish());
    let plane = source.lines().find(|s| s.contains("=PLANE(")).unwrap();
    let placement = references(plane)[0];
    let placement = source
        .lines()
        .find(|s| s.starts_with(&format!("#{placement}=")))
        .unwrap();
    let origin = references(placement)[0];
    let original = source
        .lines()
        .find(|s| s.starts_with(&format!("#{origin}=")))
        .unwrap();
    let source = source.replace(
        original,
        &format!("#{origin}=CARTESIAN_POINT('',(0.,0.,0.000000001));"),
    );
    let mut result = String::new();
    for line in source.lines() {
        if line.contains("=CARTESIAN_POINT(") || line.contains("=DIRECTION(") {
            let (prefix, xyz) = line.rsplit_once('(').unwrap();
            let p = xyz
                .trim_end_matches(&[')', ';'][..])
                .split(',')
                .map(|token| common::r(token.trim().trim_end_matches('.')))
                .collect::<Vec<_>>();
            let coordinates = match mode {
                "oblique" => {
                    assert!(p[1].is_zero());
                    [p[2].clone(), p[0].clone(), p[2].clone()]
                }
                "turned" => [p[0].clone(), -&p[1], -&p[2]],
                _ => [p[0].clone(), p[1].clone(), p[2].clone()],
            }
            .map(|v| v.to_decimal_exact().unwrap());
            result.push_str(&format!(
                "{prefix}({},{},{}));\n",
                coordinates[0], coordinates[1], coordinates[2]
            ));
        } else {
            result.push_str(line);
            result.push('\n');
        }
    }
    result.into_bytes()
}

#[test]
fn numerical_cones_export_step_axial_parameters_and_checked_geometry() {
    for (name, angle, top, mode) in [
        ("cone-30", "30.", "2.732050808", "identity"),
        ("cone-60", "60.", "6.196152423", "identity"),
        ("cone-turned", "30.", "2.732050808", "turned"),
        ("cone-oblique", "30.", "3.449489742783178", "oblique"),
    ] {
        let bytes = cone_export_source(angle, top, mode);
        let (result, work) = with_operation_budget(20_000_000, || {
            let StepImport::Numerical(imported) = checked(&bytes, "0.01")? else {
                panic!("expected numerical cone")
            };
            let body = &imported.bodies()[0].1;
            let audit = output::audit_written(
                body,
                &NumericalStepOptions::default(),
                &mut Decider::default(),
            )?;
            let file = &audit.file;
            let face = body
                .definition()
                .faces
                .iter()
                .find(|f| matches!(f.surface, catallax_shapes::Surface::Cone(_)))
                .unwrap();
            let native_min = face
                .loops
                .iter()
                .flatten()
                .flat_map(|t| &t.uv)
                .map(|uv| &uv[1])
                .min()
                .unwrap();
            let step_min = audit
                .evidence
                .checks()
                .iter()
                .filter_map(|c| {
                    let catallax_numbers::NumericalCheck::Samples {
                        name,
                        scope,
                        samples,
                    } = c
                    else {
                        return None;
                    };
                    (name == "file.actual-decimal-pcurve-edge-surface-distance"
                        && scope.ends_with(&format!("on {}", face.key)))
                    .then_some(samples)
                })
                .flatten()
                .map(|s| &s.parameters[2])
                .min()
                .unwrap();
            assert!(
                (step_min.clone() - native_min.clone()).abs() < q("1e-8"),
                "{name}"
            );
            assert!(
                file.report()
                    .checks
                    .iter()
                    .any(|s| s == "file.cone-step-v-axial-height")
            );
            let text = std::str::from_utf8(file.bytes()).unwrap();
            assert!(text.contains("CONICAL_SURFACE(") && text.contains("PCURVE("));
            assert!(
                text.contains("catallax numerical approximate")
                    && text.contains("complete_body_bound=unestablished")
            );
            file.check_written(file.bytes(), &mut Decider::default())
                .unwrap();
            if let Ok(directory) = std::env::var("CATALLAX_CONE_EVIDENCE") {
                std::fs::create_dir_all(&directory).unwrap();
                std::fs::write(
                    std::path::Path::new(&directory).join(format!("{name}.step")),
                    file.bytes(),
                )
                .unwrap();
            }
            Ok::<_, numerical::ImportError>(())
        });
        assert!(
            result.is_ok() && !work.exhausted,
            "{name}: {result:?}, {work:?}"
        );
    }
}
