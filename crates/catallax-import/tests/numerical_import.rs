//! Numerical source admission retains its category even for zero residuals.

mod common;

use catallax_import::{
    ImportRefusal, import_step_bodies,
    numerical::{self, ImportError, NumericalImport, StepImport},
    sha256_hex,
};
use catallax_numbers::{Decider, NumericalCheck, Real, with_operation_budget};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Point3, Pose, Surface};
use catallax_step::numerical::{NumericalStepOptions, build as export};
use catallax_topology::{Model, block, numerical as n};
use common::{BoxFixture, box_step, box_with_placement, boxes_step, r};

fn checked(source: &[u8]) -> Result<StepImport, ImportError> {
    numerical::import_step(
        source,
        &sha256_hex(source),
        &StepId::new("source").unwrap(),
        &n::NumericalOptions::default(),
        &mut Decider::default(),
    )
}
fn admitted(source: &[u8]) -> NumericalImport {
    let StepImport::Numerical(imported) = checked(source).unwrap() else {
        panic!("control must use the numerical path")
    };
    *imported
}
fn displaced_carrier(z: &str, width: &str, multi: bool) -> Vec<u8> {
    let fixture = BoxFixture::new([width, "1", "1"]);
    let text = if multi {
        let mut other = BoxFixture::new([width, "1", "1"]);
        other.origin[0] = r("5");
        boxes_step(&[fixture, other])
    } else {
        box_step(fixture)
    };
    let line = text
        .lines()
        .rfind(|line| line.contains("CARTESIAN_POINT("))
        .unwrap();
    let (prefix, _) = line.split_once(",(").unwrap();
    text.replacen(line, &format!("{prefix},(0.,0.,{z}));"), 1)
        .into_bytes()
}

#[test]
fn a_source_repair_has_numerical_checks_and_separate_source_history() {
    let source = box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
    let imported = admitted(source.as_bytes());
    assert_eq!(imported.bodies().len(), 1);
    assert_eq!(
        imported.record().source.sha256,
        sha256_hex(source.as_bytes())
    );
    assert!(!imported.record().source_incidence.is_empty());
    assert_eq!(
        imported.record().source_incidence[0].source_definition,
        "projected_vertex_boundary/1"
    );
    let body = &imported.bodies()[0].1;
    assert_eq!(body.options().tolerance_mm, r("0.01"));
    assert_eq!(body.counts().faces, 6);
    let samples = body
        .evidence()
        .checks()
        .iter()
        .flat_map(|check| match check {
            NumericalCheck::Samples { samples, .. } => samples.as_slice(),
            _ => &[],
        })
        .collect::<Vec<_>>();
    assert!(!samples.is_empty());
    assert!(
        samples
            .iter()
            .all(|sample| sample.distance_upper_mm.is_zero())
    );
    // Zero residual still came out of the numerical variant and numerical type.
    assert!(!body.exact_topology_faults().is_empty());
}

#[test]
fn a_previously_supported_source_keeps_the_certified_path() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let StepImport::Certified(imported) = checked(source.as_bytes()).unwrap() else {
        panic!("clean source changed category")
    };
    assert_eq!(imported.bodies().len(), 1);
    assert!(imported.model().check(imported.bodies()[0]).is_empty());
    assert_eq!(imported.source_sha256(), sha256_hex(source.as_bytes()));
}

#[test]
fn the_source_cap_is_inclusive_and_the_next_larger_gap_refuses() {
    let imported = admitted(&displaced_carrier("1.01", "1", false));
    assert_eq!(
        imported.record().face_tolerances[0].deviation_upper_mm,
        "0.01"
    );
    let error = checked(&displaced_carrier("1.010000000001", "1", false)).unwrap_err();
    let ImportError::Source(ImportRefusal::PlanarRepairNotVerified { detail }) = error else {
        panic!("{error:?}")
    };
    assert_eq!(detail.stage, "import_cap");
}

#[test]
fn the_neighbouring_feature_limit_and_requested_tolerance_still_apply() {
    let error = checked(&displaced_carrier("1.003", "0.02", false)).unwrap_err();
    let ImportError::Source(ImportRefusal::PlanarRepairNotVerified { detail }) = error else {
        panic!("{error:?}")
    };
    assert_eq!(detail.stage, "nearest_feature");
    let source = displaced_carrier("1.005", "1", false);
    let error = numerical::import_step(
        &source,
        &sha256_hex(&source),
        &StepId::new("source").unwrap(),
        &n::NumericalOptions {
            tolerance_mm: r("0.001"),
        },
        &mut Decider::default(),
    )
    .unwrap_err();
    assert!(
        matches!(error, ImportError::Numerical(ref e) if e.code == "import.numerical-source-limit")
    );
}

#[test]
fn the_exact_entry_rolls_back_all_source_solids_on_failure() {
    let mut model = Model::new();
    let prior = block(&mut model, [r("2"), r("3"), r("4")], &Pose::identity()).unwrap();
    let counts = model.counts(prior);
    let mut expected = model.clone();
    let expected_next = block(&mut expected, [r("7"), r("8"), r("9")], &Pose::identity()).unwrap();
    let source = displaced_carrier("1.02", "1", true);
    assert!(
        import_step_bodies(
            &mut model,
            &source,
            &sha256_hex(&source),
            &StepId::new("source").unwrap()
        )
        .is_err()
    );
    assert_eq!(model.counts(prior), counts);
    assert!(model.check(prior).is_empty());
    let next = block(&mut model, [r("7"), r("8"), r("9")], &Pose::identity()).unwrap();
    assert_eq!(next, expected_next);
    assert_eq!(
        model.faces_of_body(next),
        expected.faces_of_body(expected_next)
    );
}

#[test]
fn exhausted_source_queries_return_no_numerical_body() {
    let source = displaced_carrier("1.001", "1", false);
    let (result, work) = with_operation_budget(1, || checked(&source));
    assert!(work.exhausted);
    assert!(matches!(
        result,
        Err(ImportError::Source(ImportRefusal::ResourceLimit))
    ));
}

#[test]
fn true_crossings_are_not_repaired_numerically() {
    let mut soup = common::Soup::new();
    soup.box_(["0", "0", "0"], ["2", "2", "2"]);
    soup.box_(["1", "1", "1"], ["3", "3", "3"]);
    assert!(matches!(
        checked(soup.finish().as_bytes()),
        Err(ImportError::Source(ImportRefusal::SelfIntersecting { .. }))
    ));
}

fn translate(p: &Point3, delta: &[Real; 3]) -> Point3 {
    Point3::new(&p.x + &delta[0], &p.y + &delta[1], &p.z + &delta[2])
}

#[test]
fn a_numerical_source_takes_a_pose_edit_and_writes_approximate_step() {
    let source = box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
    let imported = admitted(source.as_bytes());
    let original = &imported.bodies()[0].1;
    let mut definition = original.definition().clone();
    let shift = [r("3"), r("5"), r("7")];
    let delta = shift.clone().map(Real::from);
    for vertex in &mut definition.vertices {
        for i in 0..3 {
            vertex.position[i] = &vertex.position[i] + &shift[i];
        }
    }
    for edge in &mut definition.edges {
        let Curve::Line(line) = &mut edge.curve else {
            unreachable!()
        };
        line.origin = translate(&line.origin, &delta);
    }
    for face in &mut definition.faces {
        let Surface::Plane(plane) = &mut face.surface else {
            unreachable!()
        };
        plane.origin = translate(&plane.origin, &delta);
    }
    let moved = n::build(
        definition,
        original.options().clone(),
        &mut Decider::default(),
    )
    .unwrap();
    let step = export(
        &moved,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(String::from_utf8_lossy(step.bytes()).contains("numerical approximate"));
    assert!(step.report().max_scalar_rounding_error_mm <= step.report().scalar_rounding_limit_mm);
    let dir = std::env::temp_dir().join(format!("cx-import-numerical-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("answer.step");
    std::fs::write(&path, step.bytes()).unwrap();
    step.check_written(&std::fs::read(&path).unwrap(), &mut Decider::default())
        .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}
