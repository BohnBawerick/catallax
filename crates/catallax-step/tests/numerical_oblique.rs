//! Independently rotated N03 geometry exercises general line and plane frames.

#[path = "common/n03_tube.rs"]
mod tube;

use catallax_numbers::{Decider, NumericalCheck, Rational, Real, with_operation_budget};
use catallax_shapes::{Curve, Point3, Pose, Rotation, Surface, Vector3};
use catallax_step::numerical::{self, NumericalStepOptions};
use catallax_topology::numerical::{self as source, NumericalDefinition, NumericalOptions};

fn rotate(p: &[Rational; 3]) -> [Rational; 3] {
    let [a, b, c] = [7, 24, 25].map(Rational::from_i64);
    [
        (&a * &p[0] + &b * &p[2]) / c.clone(),
        p[1].clone(),
        (&a * &p[2] - &b * &p[0]) / c,
    ]
}
fn rotated_point(p: &Point3) -> Point3 {
    let value = rotate(&[&p.x, &p.y, &p.z].map(|q| q.as_rational().unwrap()));
    let [x, y, z] = value.map(Real::from);
    Point3::new(x, y, z)
}
fn rotated_vector(p: &Vector3) -> Vector3 {
    let p = rotated_point(&Point3::new(p.x.clone(), p.y.clone(), p.z.clone()));
    Vector3::new(p.x, p.y, p.z)
}
fn definition() -> NumericalDefinition {
    let mut d = tube::tube();
    d.name = "oblique two-loop tube".into();
    for vertex in &mut d.vertices {
        vertex.position = rotate(&vertex.position);
    }
    for edge in &mut d.edges {
        let Curve::Line(line) = &mut edge.curve else {
            unreachable!()
        };
        line.origin = rotated_point(&line.origin);
        line.dir = rotated_vector(&line.dir);
    }
    for face in &mut d.faces {
        let Surface::Plane(plane) = &mut face.surface else {
            unreachable!()
        };
        plane.origin = rotated_point(&plane.origin);
        plane.normal = rotated_vector(&plane.normal);
        plane.ref_dir = rotated_vector(&plane.ref_dir);
    }
    d
}

#[test]
fn oblique_round_trip_keeps_topology_and_separate_rounding_scopes() {
    let body = source::build(
        definition(),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let before = format!("{:?}", body.definition());
    let (artifact, work) = with_operation_budget(20_000_000, || {
        numerical::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
    });
    let artifact = artifact.unwrap();
    assert!(!work.exhausted);
    assert_eq!(artifact.report().counts, [10, 24, 16]);
    assert_eq!(artifact.report().source_evidence, *body.evidence());
    assert!(artifact.report().max_direction_ratio_rounding_error > Rational::zero());
    assert_eq!(
        artifact.report().max_scalar_rounding_error_mm,
        Rational::zero()
    );
    assert_eq!(before, format!("{:?}", body.definition()));
    let samples = artifact
        .report()
        .emitted_evidence
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { samples, .. } => Some(samples.as_slice()),
            _ => None,
        })
        .flatten()
        .collect::<Vec<_>>();
    assert!(!samples.is_empty());
    assert!(
        samples
            .iter()
            .all(|sample| sample.distance_upper_mm <= body.options().tolerance_mm)
    );
    let capture = std::env::var_os("CATALLAX_OBLIQUE_CAPTURE");
    let directory = capture
        .clone()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("cx-oblique-{}", std::process::id())));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("oblique-tube.step");
    std::fs::write(&path, artifact.bytes()).unwrap();
    artifact
        .check_written(&std::fs::read(&path).unwrap(), &mut Decider::default())
        .unwrap();
    if capture.is_none() {
        std::fs::remove_dir_all(directory).unwrap();
    }
}

#[test]
fn coarse_oblique_payload_must_pass_its_own_checks() {
    let body = source::build(
        definition(),
        NumericalOptions {
            tolerance_mm: tube::q("0.00000001"),
            vertex_tolerance_mm: None,
        },
        &mut Decider::default(),
    )
    .unwrap();
    let options = NumericalStepOptions {
        significant_digits: 1,
        max_scalar_rounding_error_mm: Rational::one(),
        ..Default::default()
    };
    let error = numerical::build(&body, &options, &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "numbers.numerical-tolerance-exceeded");
}

#[test]
fn oblique_export_keeps_tolerated_source_residuals_numerical() {
    let mut d = definition();
    d.vertices[0].position[0] = &d.vertices[0].position[0] + &tube::q("0.005");
    let body = source::build(d, NumericalOptions::default(), &mut Decider::default()).unwrap();
    let artifact = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let minimum_residual = tube::q("0.0049");
    assert!(
        artifact
            .report()
            .emitted_evidence
            .checks()
            .iter()
            .any(|check| matches!(check,
                NumericalCheck::Samples { samples, .. }
                    if samples.iter().any(|sample| sample.distance_upper_mm > minimum_residual)))
    );
}

#[test]
fn numerical_uv_roundoff_does_not_move_or_refuse_a_planar_boundary() {
    let mut d = tube::tube();
    for uv in d.faces[0]
        .loops
        .iter_mut()
        .flatten()
        .flat_map(|trim| &mut trim.uv)
    {
        uv[0] = &uv[0] + &tube::q("1e-38");
    }
    let body = source::build(d, NumericalOptions::default(), &mut Decider::default()).unwrap();
    let artifact = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(artifact.report().counts, [10, 24, 16]);
    assert!(artifact.report().emitted_evidence.repairs().is_empty());
}

#[test]
fn a_tolerated_plane_carrier_gap_remains_in_the_payload_checks() {
    let mut d = tube::tube();
    let Surface::Plane(plane) = &mut d.faces[0].surface else {
        unreachable!()
    };
    plane.origin.z = &plane.origin.z + &Real::from(tube::q("0.005"));
    let body = source::build(d, NumericalOptions::default(), &mut Decider::default()).unwrap();
    let artifact = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let minimum_residual = tube::q("0.0049");
    assert!(
        artifact
            .report()
            .emitted_evidence
            .checks()
            .iter()
            .any(|check| matches!(check,
                NumericalCheck::Samples { samples, .. }
                    if samples.iter().any(|sample| sample.distance_upper_mm > minimum_residual)))
    );
}

#[test]
fn successive_rotations_reconstruct_the_written_orthogonal_frame() {
    let initial = source::build(
        definition(),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let pose = Pose {
        rotation: Rotation::from_quaternion(
            Rational::from_i64(2),
            Rational::zero(),
            Rational::zero(),
            Rational::one(),
        )
        .unwrap(),
        translation: Vector3::zero(),
    };
    let mut d = initial.definition().clone();
    for vertex in &mut d.vertices {
        let [x, y, z] = vertex.position.clone().map(Real::from);
        let p = pose.apply_point(&Point3::new(x, y, z));
        vertex.position = [&p.x, &p.y, &p.z].map(|q| q.as_rational().unwrap());
    }
    for edge in &mut d.edges {
        let Curve::Line(line) = &mut edge.curve else {
            unreachable!()
        };
        line.origin = pose.apply_point(&line.origin);
        line.dir = pose.apply_vector(&line.dir);
    }
    for face in &mut d.faces {
        let Surface::Plane(plane) = &mut face.surface else {
            unreachable!()
        };
        plane.origin = pose.apply_point(&plane.origin);
        plane.normal = pose.apply_vector(&plane.normal);
        plane.ref_dir = pose.apply_vector(&plane.ref_dir);
    }
    let body = initial
        .rebuild(d, NumericalOptions::default(), &mut Decider::default())
        .unwrap();
    let artifact = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(artifact.report().counts, [10, 24, 16]);
    assert_eq!(artifact.report().source_history.len(), 1);
    assert!(artifact.report().max_direction_ratio_rounding_error > Rational::zero());
    if let Some(directory) = std::env::var_os("CATALLAX_OBLIQUE_CAPTURE") {
        let directory = std::path::PathBuf::from(directory);
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("oblique-3d-tube.step");
        std::fs::write(&path, artifact.bytes()).unwrap();
        artifact
            .check_written(&std::fs::read(path).unwrap(), &mut Decider::default())
            .unwrap();
    }
}
