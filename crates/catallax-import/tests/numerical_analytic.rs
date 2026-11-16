//! Actual Part 21 circular edges keep shared identity through numerical admission.
mod common;
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_import::{
    numerical::{self, StepImport},
    sha256_hex,
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
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

fn opened_checked(bytes: &[u8], tolerance: &str) -> catallax_import::tolerant::opened::OpenedStep {
    use catallax_import::tolerant::{ToleranceOptions, opened};
    opened::open_step(
        bytes,
        &sha256_hex(bytes),
        &StepId::new("retained-spline-image").unwrap(),
        ToleranceOptions {
            initial_mm: q("0.000000000001"),
            maximum_mm: q(tolerance),
        },
        &mut Decider::default(),
    )
    .unwrap()
}

fn spline_boundary(offset: &str, scales: [&str; 4]) -> Vec<u8> {
    let source = displaced_wall("0.000000000001");
    let recognized = catallax_import::read_step(&source).unwrap();
    let edge = recognized.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .find(|e| matches!(e.curve, catallax_import::Curve::Line(..)))
        .unwrap();
    let (a, b) = if edge.curve_same_sense {
        (&edge.start.point, &edge.end.point)
    } else {
        (&edge.end.point, &edge.start.point)
    };
    let direction: [Rational; 3] = std::array::from_fn(|i| &b[i] - &a[i]);
    let axis = (0..3).find(|&i| direction[i].is_zero()).unwrap();
    let offset = q(offset);
    let controls = scales.map(|scale| {
        let scale = q(scale);
        let mut p: [Rational; 3] = std::array::from_fn(|i| {
            let component = &scale * &direction[i];
            &a[i] + &component
        });
        p[axis] = &p[axis] + &offset;
        p
    });
    with_spline(&source, edge, &controls, &["1.", "2.", "3.", "1."])
}

fn circular_spline_boundary(weight: &str, height: &str, reverse: bool) -> Vec<u8> {
    let source = displaced_wall("0.000000000001");
    let recognized = catallax_import::read_step(&source).unwrap();
    let edge = recognized.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .find(|e| matches!(e.curve, catallax_import::Curve::Circle(..)))
        .unwrap();
    let catallax_import::Curve::Circle(frame, _) = &edge.curve else {
        unreachable!()
    };
    let (a, b) = if edge.curve_same_sense {
        (&edge.start.point, &edge.end.point)
    } else {
        (&edge.end.point, &edge.start.point)
    };
    let height = q(height);
    let middle = std::array::from_fn(|i| {
        let sum = &a[i] + &b[i];
        let flat = &sum - &frame.origin[i];
        let lift = &height * &frame.axis[i];
        &flat + &lift
    });
    let mut controls = [a.clone(), middle, b.clone()];
    let mut edge = edge.clone();
    if reverse {
        controls.reverse();
        edge.curve_same_sense = !edge.curve_same_sense;
    }
    with_spline(&source, &edge, &controls, &["1.", weight, "1."])
}

fn with_spline(
    source: &[u8],
    edge: &catallax_import::OrientedEdge,
    controls: &[[Rational; 3]],
    weights: &[&str],
) -> Vec<u8> {
    let degree = controls.len() - 1;
    let multiplicity = degree + 1;
    let text = std::str::from_utf8(source).unwrap();
    let edge_line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", edge.edge)))
        .unwrap();
    let curve = references(edge_line)[2];
    let original = text
        .lines()
        .find(|line| line.starts_with(&format!("#{curve}=")))
        .unwrap();
    let first = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix('#')?
                .split_once('=')?
                .0
                .parse::<u64>()
                .ok()
        })
        .max()
        .unwrap()
        + 1;
    let decimal = |value: &Rational| {
        let text = value.to_string();
        let value = if let Some((n, d)) = text.split_once('/') {
            n.parse::<f64>().unwrap() / d.parse::<f64>().unwrap()
        } else {
            text.parse::<f64>().unwrap()
        };
        format!("{value:.12}")
    };
    let mut points = String::new();
    let mut ids = Vec::new();
    for (index, p) in controls.iter().enumerate() {
        let id = first + index as u64;
        points.push_str(&format!(
            "#{id}=CARTESIAN_POINT('',({}));\n",
            p.each_ref().map(decimal).join(",")
        ));
        ids.push(format!("#{id}"));
    }
    let replacement = format!(
        "#{curve}=(BOUNDED_CURVE() B_SPLINE_CURVE({degree},({}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({multiplicity},{multiplicity}),(0.,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({})) REPRESENTATION_ITEM(''));",
        ids.join(","),
        weights.join(",")
    );
    let edge_prefix = edge_line.rsplit_once(',').unwrap().0;
    let edge_definition = format!(
        "{edge_prefix},{});",
        if edge.curve_same_sense { ".T." } else { ".F." }
    );
    let text = text
        .replace(original, &replacement)
        .replace(edge_line, &edge_definition);
    let end = text.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &text[..end], points, &text[end..]).into_bytes()
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
    assert_eq!(admission.method, "retained_shared_analytic_boundaries/4");
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
fn a_regular_weighted_spline_line_image_rebuilds_only_finite_parameters() {
    let source = spline_boundary("0.00000001", ["-1", "0", "1", "2"]);
    let opened = opened_checked(&source, "0.01");
    let Some(StepImport::Numerical(imported)) = opened.admission(0) else {
        panic!(
            "expected numerical admission: {:?}",
            opened.bodies()[0].refusal()
        );
    };
    let images = &imported.source_admission().unwrap().linear_spline_images;
    assert_eq!(images.len(), 1);
    let image = &images[0];
    assert_eq!(
        image.endpoint_distance_squared_mm2,
        [
            q("0.0000000000000001").to_string(),
            q("0.0000000000000001").to_string()
        ]
    );
    let recognized = catallax_import::read_step(&source).unwrap();
    assert!(
        recognized.solids[0]
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flat_map(|b| &b.edges)
            .any(|e| e.edge == image.edge && matches!(e.curve, catallax_import::Curve::BSpline(_)))
    );
    let body = &imported.bodies()[0].1;
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (5, 9, 6));
    assert!(
        body.definition()
            .edges
            .iter()
            .any(|e| e.key == image.edge && matches!(e.curve, catallax_shapes::Curve::Line(_)))
    );
    assert_eq!(imported.record().source.sha256, sha256_hex(&source));
    let refused = opened_checked(&source, "0.000000001");
    assert!(matches!(
        refused.body(0),
        Some(catallax_import::tolerant::opened::BodyView::Surfaces { .. })
    ));
    assert!(matches!(refused.bodies()[0].refusal(),
        Some(catallax_import::tolerant::Error::Numerical(r))
        if r.code == "import.numerical-spline-domain"));
}

#[test]
fn finite_spline_extent_and_backtracking_are_not_infinite_line_permissions() {
    use catallax_import::tolerant::{Error, opened::BodyView};
    for (offset, scales) in [
        ("0", ["0.00000002", "0.3", "0.6", "1"]),
        ("0.00000001", ["-1", "3", "-2", "2"]),
    ] {
        let source = spline_boundary(offset, scales);
        let tolerance = if offset == "0" { "0.000000001" } else { "0.01" };
        let opened = opened_checked(&source, tolerance);
        let expected = if offset == "0" {
            "import.numerical-spline-domain"
        } else {
            "import.numerical-unsupported-family"
        };
        let refusal = opened.bodies()[0].refusal();
        assert!(
            matches!(refusal, Some(Error::Numerical(r)) if r.code == expected),
            "{refusal:?}"
        );
        let Some(BodyView::Surfaces {
            source: retained, ..
        }) = opened.body(0)
        else {
            panic!("a declined line image must retain a surface body");
        };
        let original = catallax_import::read_step(&source).unwrap();
        let original = original.solids[0]
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flat_map(|b| &b.edges)
            .find(|e| matches!(e.curve, catallax_import::Curve::BSpline(_)))
            .unwrap();
        let retained = retained
            .edges()
            .iter()
            .find(|e| e.name == original.edge)
            .unwrap();
        assert_eq!(retained.curve, original.curve);
    }
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
            let file = output::build(
                body,
                &NumericalStepOptions::default(),
                &mut Decider::default(),
            )?;
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
            let step_min = file
                .report()
                .emitted_evidence
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

#[test]
fn cylindrical_spline_arcs_have_continuous_spatial_and_finite_image_bounds() {
    for (height, reverse) in [
        ("0", false),
        ("0", true),
        ("0.0001", false),
        ("0.0001", true),
    ] {
        let source = circular_spline_boundary("0.7071067811865475", height, reverse);
        let opened = opened_checked(&source, "0.01");
        let Some(StepImport::Numerical(imported)) = opened.admission(0) else {
            panic!(
                "expected circular image: {:?}",
                opened.bodies()[0].refusal()
            );
        };
        let images = &imported.source_admission().unwrap().circular_spline_images;
        assert_eq!(images.len(), 1);
        let image = &images[0];
        let mapped = q(&image.mapped_curve_deviation_upper_mm);
        let endpoint = q(&image.endpoint_distance_upper_mm);
        let finite = q(&image.finite_arc_deviation_upper_mm);
        assert!(finite <= q("0.01") && finite >= mapped && finite >= endpoint);
        if height != "0" {
            assert!(mapped >= q("0.00004"));
            let declined = opened_checked(&source, "0.000001");
            assert!(matches!(
                declined.body(0),
                Some(catallax_import::tolerant::opened::BodyView::Surfaces { .. })
            ));
        }
        let body = &imported.bodies()[0].1;
        assert!(
            body.definition()
                .edges
                .iter()
                .any(|edge| edge.key == image.edge
                    && matches!(edge.curve, catallax_shapes::Curve::Circle(_)))
        );
        let counts = body.counts();
        assert_eq!((counts.faces, counts.edges, counts.vertices), (5, 9, 6));
    }
}

#[test]
fn cylindrical_endpoints_do_not_approve_a_bowed_source_spline() {
    use catallax_import::tolerant::opened::BodyView;
    let source = circular_spline_boundary("0.5", "0", false);
    let opened = opened_checked(&source, "0.01");
    let Some(BodyView::Surfaces {
        source: retained, ..
    }) = opened.body(0)
    else {
        panic!("an unbounded circle fit must retain source surfaces");
    };
    let original = catallax_import::read_step(&source).unwrap();
    let original = original.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .find(|e| matches!(e.curve, catallax_import::Curve::BSpline(_)))
        .unwrap();
    let retained = retained
        .edges()
        .iter()
        .find(|e| e.name == original.edge)
        .unwrap();
    assert_eq!(retained.curve, original.curve);
}
