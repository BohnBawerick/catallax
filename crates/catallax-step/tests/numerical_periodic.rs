//! One shared cylinder seam carries two explicit periodic pcurves.
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_shapes::{Curve, Point3, Surface};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::{Sense, numerical as n};

#[test]
fn a_complete_bore_emits_one_seam_curve_with_two_uses() {
    let body = n::build(
        fixture::bore(),
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
    assert_eq!(text.matches("=SEAM_CURVE(").count(), 1);
    assert_eq!(file.report().counts, [7, 15, 10]);
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&path).join("periodic-bore.step"),
            file.bytes(),
        )
        .unwrap();
    }
}

fn shifted_center_bore() -> n::NumericalDefinition {
    let mut definition = fixture::bore();
    let shift = catallax_numbers::Real::from(fixture::q("-1.0004"));
    let shifted = |p: &Point3| Point3::new(p.x.clone(), &p.y + &shift, p.z.clone());
    for vertex in &mut definition.vertices {
        vertex.position[1] = &vertex.position[1] + &fixture::q("-1.0004");
    }
    for edge in &mut definition.edges {
        match &mut edge.curve {
            Curve::Line(line) => line.origin = shifted(&line.origin),
            Curve::Circle(circle) => circle.center = shifted(&circle.center),
            _ => unreachable!(),
        }
    }
    for face in &mut definition.faces {
        match &mut face.surface {
            Surface::Plane(plane) => plane.origin = shifted(&plane.origin),
            Surface::Cylinder(cylinder) => cylinder.origin = shifted(&cylinder.origin),
            _ => unreachable!(),
        }
    }
    definition
}

#[test]
fn a_shifted_center_uses_the_written_center_for_its_zero_phase() {
    let body = n::build(
        shifted_center_bore(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    for digits in [4, 17] {
        let options = NumericalStepOptions {
            significant_digits: digits,
            max_scalar_rounding_error_mm: fixture::q("0.001"),
            ..Default::default()
        };
        let audit = output::audit_written(&body, &options, &mut Decider::default()).unwrap();
        let file = audit.file;
        let checks = audit
            .evidence
            .checks()
            .iter()
            .filter_map(|check| match check {
                catallax_numbers::NumericalCheck::Samples { name, samples, .. }
                    if name == "file.closed-circle-zero-phase-tangent-bound" =>
                {
                    Some(samples)
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(checks.len(), 2);
        assert!(
            checks
                .iter()
                .flat_map(|samples| samples.iter())
                .all(|sample| sample.distance_upper_mm <= fixture::q("0.000000001"))
        );
        if digits == 17
            && let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE")
        {
            std::fs::write(
                std::path::Path::new(&path).join("periodic-bore-shifted-center.step"),
                file.bytes(),
            )
            .unwrap();
        }
    }
}

#[test]
fn a_direction_that_rounds_away_from_zero_phase_is_refused() {
    let mut definition = fixture::bore();
    let c = Rational::from_ratio(20, 29).unwrap();
    let s = Rational::from_ratio(21, 29).unwrap();
    let phase = Ball::atan2(&Ball::from_rational(&s, 128), &Ball::from_rational(&c, 128))
        .unwrap()
        .mid();
    let period = Ball::pi(128).mid() * Rational::from_i64(2);
    for vertex in &mut definition.vertices[8..10] {
        vertex.position[0] = c.clone();
        vertex.position[1] = s.clone();
    }
    for edge in &mut definition.edges {
        if matches!(edge.curve, Curve::Circle(_)) {
            edge.parameters = [phase.clone(), &phase + &period];
        } else if edge.key == "cylinder-u-seam"
            && let Curve::Line(line) = &mut edge.curve
        {
            line.origin = Point3::new(c.clone().into(), s.clone().into(), 0.into());
        }
    }
    for face in &mut definition.faces {
        if matches!(face.surface, Surface::Plane(_)) && face.loops.len() == 2 {
            face.loops[1][0].uv = [[c.clone(), s.clone()], [c.clone(), s.clone()]];
        } else if matches!(face.surface, Surface::Cylinder(_)) {
            for trim in face.loops.iter_mut().flatten() {
                for uv in &mut trim.uv {
                    uv[0] = &uv[0] + &phase;
                }
            }
        }
    }
    let body = n::build(
        definition,
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let coarse = NumericalStepOptions {
        significant_digits: 4,
        max_scalar_rounding_error_mm: fixture::q("0.001"),
        ..Default::default()
    };
    let error = output::build(&body, &coarse, &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "step.numerical-closed-circle-phase");
    let fine = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&path).join("periodic-bore-oblique-anchor.step"),
            fine.bytes(),
        )
        .unwrap();
    }
}

#[test]
fn equivalent_closed_rim_directions_and_periods_emit_the_same_step() {
    let period = Ball::pi(128).mid() * Rational::from_i64(2);
    let export = |definition| {
        let body = n::build(
            definition,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap()
    };
    // The negative direction is the real 00003390 failure. Period shifts also
    // exercise the natural range which a STEP reader reconstructs from vertices.
    for phased in [false, true] {
        let mut baseline = fixture::bore();
        let index = baseline
            .edges
            .iter()
            .position(|edge| edge.key == "rim-0")
            .unwrap();
        let start = if phased {
            let catallax_shapes::Curve::Circle(circle) = &mut baseline.edges[index].curve else {
                unreachable!()
            };
            circle.ref_dir = catallax_shapes::RVector3::new(0.into(), 1.into(), 0.into()).to_real();
            -(&period / &Rational::from_i64(4))
        } else {
            Rational::zero()
        };
        baseline.edges[index].parameters = [start.clone(), &start + &period];
        let expected = export(baseline.clone());
        if phased && let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&path).join("periodic-bore-phased.step"),
                expected.bytes(),
            )
            .unwrap();
        }
        for (direction, periods) in [(1, -1), (1, 1), (-1, -1), (-1, 0), (-1, 1)] {
            let mut definition = baseline.clone();
            let phase = &start + &(&period * &Rational::from_i64(periods));
            definition.edges[index].parameters = [
                phase.clone(),
                &phase + &(&period * &Rational::from_i64(direction)),
            ];
            if direction < 0 {
                for trim in definition
                    .faces
                    .iter_mut()
                    .flat_map(|face| face.loops.iter_mut().flatten())
                {
                    if trim.edge == index {
                        trim.sense = match trim.sense {
                            Sense::Same => Sense::Reversed,
                            Sense::Reversed => Sense::Same,
                        };
                    }
                }
            }
            let file = export(definition);
            assert!(
                file.bytes() == expected.bytes(),
                "phased={phased}, direction={direction}, periods={periods}"
            );
        }
    }
}

#[test]
fn multiple_bores_and_complete_annular_caps_write_checked_step() {
    for (name, definition, expected) in [
        (
            "multiple-bores",
            fixture::bores(&["-0.5", "0.5"], "0.4"),
            [8, 18, 12],
        ),
        ("annulus", fixture::annulus("2"), [4, 6, 4]),
    ] {
        let body = n::build(
            definition,
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
        assert_eq!(file.report().counts, expected);
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&path).join(format!("{name}.step")),
                file.bytes(),
            )
            .unwrap();
        }
    }
}

#[test]
fn cylinder_cells_keep_reader_compatible_charts_across_period_gauges() {
    let period = Ball::pi(128).mid() * Rational::from_i64(2);
    let export = |definition| {
        let body = n::build(
            definition,
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap()
    };
    for negative_height in [false, true] {
        let mut baseline = fixture::annulus("2");
        if negative_height {
            let height = fixture::q("2");
            let shift = catallax_numbers::Real::from(height.clone());
            for face in &mut baseline.faces {
                if let Surface::Cylinder(cylinder) = &mut face.surface {
                    cylinder.origin.z = &cylinder.origin.z + &shift;
                    for uv in face
                        .loops
                        .iter_mut()
                        .flatten()
                        .flat_map(|trim| &mut trim.uv)
                    {
                        uv[1] = &uv[1] - &height;
                    }
                }
            }
        }
        let expected = export(baseline.clone());
        for turns in [-2, -1, 0, 1, 2] {
            let mut definition = baseline.clone();
            let shift = &period * &Rational::from_i64(turns);
            for face in &mut definition.faces {
                if matches!(face.surface, Surface::Cylinder(_)) {
                    for uv in face
                        .loops
                        .iter_mut()
                        .flatten()
                        .flat_map(|trim| &mut trim.uv)
                    {
                        uv[0] = &uv[0] + &shift;
                    }
                }
            }
            let file = export(definition);
            assert_eq!(
                file.bytes(),
                expected.bytes(),
                "negative_height={negative_height}, turns={turns}"
            );
            if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                std::fs::write(
                    std::path::Path::new(&path)
                        .join(format!("period-gauge-{negative_height}-{turns}.step")),
                    file.bytes(),
                )
                .unwrap();
            }
        }
    }
}
