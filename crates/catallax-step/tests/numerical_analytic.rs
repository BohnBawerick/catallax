//! Approximate analytic STEP preserves circular boundaries and explicit pcurves.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, NumericalCheck, Real, with_operation_budget};
use catallax_shapes::{Curve, Point3, Surface};
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
