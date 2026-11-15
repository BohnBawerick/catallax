//! Approximate analytic STEP preserves circular boundaries and explicit pcurves.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, NumericalCheck, Real, with_operation_budget};
use catallax_shapes::{Curve, Point3, Surface, Vector3};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::{q, quarter};

fn save(name: &str, bytes: &[u8]) {
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(std::path::Path::new(&path).join(name), bytes).unwrap();
    }
}
#[test]
fn circular_boundaries_emit_approximate_step_with_checked_cylindrical_pcurves() {
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
    assert!(text.contains("catallax numerical approximate"));
    assert!(text.contains("CIRCLE("));
    assert!(text.contains("CYLINDRICAL_SURFACE("));
    assert!(text.contains("PCURVE("));
    assert!(!text.contains("UNCERTAINTY_MEASURE_WITH_UNIT"));
    assert!(file.report().emitted_evidence.checks().iter().any(|check| matches!(check,NumericalCheck::Samples { name,.. } if name=="file.actual-decimal-pcurve-edge-surface-distance")));
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    save("quarter-cylinder.step", file.bytes());
    let mut changed = file.bytes().to_vec();
    changed[0] ^= 1;
    assert_eq!(
        file.check_written(&changed, &mut Decider::default())
            .unwrap_err()
            .code,
        "step.numerical-file-changed"
    );
}

#[test]
fn analytic_pose_rebuild_retains_history_and_rechecks_the_written_file() {
    let body = n::build(
        quarter(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let mut d = body.definition().clone();
    let shift = [q("3"), q("5"), q("7")];
    let translated = |p: &Point3| {
        Point3::new(
            &p.x + &Real::from(shift[0].clone()),
            &p.y + &Real::from(shift[1].clone()),
            &p.z + &Real::from(shift[2].clone()),
        )
    };
    for vertex in &mut d.vertices {
        for (i, value) in vertex.position.iter_mut().enumerate() {
            *value = &*value + &shift[i];
        }
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            Curve::Line(line) => line.origin = translated(&line.origin),
            Curve::Circle(circle) => circle.center = translated(&circle.center),
            _ => unreachable!(),
        }
    }
    for face in &mut d.faces {
        match &mut face.surface {
            Surface::Plane(plane) => plane.origin = translated(&plane.origin),
            Surface::Cylinder(cylinder) => cylinder.origin = translated(&cylinder.origin),
            _ => unreachable!(),
        }
    }
    let moved = body
        .rebuild(d, body.options().clone(), &mut Decider::default())
        .unwrap();
    let file = output::build(
        &moved,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(file.report().source_history.len(), 1);
    assert_eq!(file.report().counts, [5, 9, 6]);
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    save("quarter-cylinder-moved.step", file.bytes());
}

#[test]
fn distant_periodic_coordinates_and_exhausted_writers_return_no_file() {
    let body = n::build(
        quarter(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let (result, work) = with_operation_budget(0, || {
        output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted && result.is_err());
    let mut d = quarter();
    let turns = Ball::pi(128).mid() * q("2000");
    for edge in &mut d.edges {
        if matches!(edge.curve, Curve::Circle(_)) {
            for parameter in &mut edge.parameters {
                *parameter = &*parameter + &turns;
            }
        }
    }
    let body = n::build(d, n::NumericalOptions::default(), &mut Decider::default()).unwrap();
    let refusal = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.numerical-unsupported");
    assert!(refusal.detail.contains("periodic chart"));
}

#[test]
fn partial_circle_exchange_uses_the_reader_parameter_chart() {
    let period = Ball::pi(128).mid() * q("2");
    for quarter_turns in [0, 1, 2, 3] {
        let mut canonical_bytes = None;
        for turns in [-1, 0, 1] {
            let mut definition = quarter();
            let reference = match quarter_turns {
                0 => [1, 0, 0],
                1 => [0, 1, 0],
                2 => [-1, 0, 0],
                _ => [0, -1, 0],
            };
            let phase = &(&period * &q(&turns.to_string()))
                - &(&period * &q(&format!("{quarter_turns}/4")));
            for edge in &mut definition.edges {
                if let Curve::Circle(circle) = &mut edge.curve {
                    circle.ref_dir = catallax_shapes::RVector3::new(
                        reference[0].into(),
                        reference[1].into(),
                        reference[2].into(),
                    )
                    .to_real();
                    for parameter in &mut edge.parameters {
                        *parameter = &*parameter + &phase;
                    }
                }
            }
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
            file.check_written(file.bytes(), &mut Decider::default())
                .unwrap();
            assert_eq!(file.report().counts, [5, 9, 6]);
            if let Some(bytes) = &canonical_bytes {
                assert_eq!(file.bytes(), bytes, "equivalent periodic arc chart");
            } else {
                canonical_bytes = Some(file.bytes().to_vec());
            }
            save(
                &format!("quarter-chart-{quarter_turns}-{turns}.step"),
                file.bytes(),
            );
        }
    }
}

#[test]
fn partial_circle_phase_is_checked_after_direction_rounding() {
    let mut definition = quarter();
    let c = q("20/29");
    let s = q("21/29");
    let rotate = |x: &Real, y: &Real| {
        [
            &(x * &Real::from(c.clone())) - &(y * &Real::from(s.clone())),
            &(x * &Real::from(s.clone())) + &(y * &Real::from(c.clone())),
        ]
    };
    let point = |p: &Point3| {
        let [x, y] = rotate(&p.x, &p.y);
        Point3::new(x, y, p.z.clone())
    };
    let vector = |v: &Vector3| {
        let [x, y] = rotate(&v.x, &v.y);
        Vector3::new(x, y, v.z.clone())
    };
    for vertex in &mut definition.vertices {
        let [x, y, _] = vertex.position.clone();
        vertex.position[0] = &(&x * &c) - &(&y * &s);
        vertex.position[1] = &(&x * &s) + &(&y * &c);
    }
    for edge in &mut definition.edges {
        match &mut edge.curve {
            Curve::Line(line) => {
                line.origin = point(&line.origin);
                line.dir = vector(&line.dir);
            }
            Curve::Circle(circle) => {
                circle.center = point(&circle.center);
                circle.normal = vector(&circle.normal);
                circle.ref_dir = vector(&circle.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    for face in &mut definition.faces {
        match &mut face.surface {
            Surface::Plane(plane) => {
                plane.origin = point(&plane.origin);
                plane.normal = vector(&plane.normal);
                plane.ref_dir = vector(&plane.ref_dir);
            }
            Surface::Cylinder(cylinder) => {
                cylinder.origin = point(&cylinder.origin);
                cylinder.axis = vector(&cylinder.axis);
                cylinder.ref_dir = vector(&cylinder.ref_dir);
            }
            _ => unreachable!(),
        }
    }
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
    save("quarter-oblique-anchor.step", file.bytes());
    let error = output::build(
        &body,
        &NumericalStepOptions {
            significant_digits: 4,
            max_scalar_rounding_error_mm: q("1/100"),
            ..NumericalStepOptions::default()
        },
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "step.numerical-arc-phase");
}
