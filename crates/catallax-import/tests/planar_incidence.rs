//! Bounded, recorded repair of inconsistent planar STEP source data.

mod common;

use catallax_import::{
    Approval, Decision, Disposition, GapKind, ImportOptions, ImportRefusal, PlanarRepairMethod,
    Surface, apply, propose, read_step, rebuild, sha256_hex,
};
use catallax_numbers::{Rational, with_operation_budget};
use catallax_topology::Model;
use common::{BoxFixture, box_step, box_with_placement};

fn tilted_top(x: &str) -> Vec<u8> {
    let bytes = box_with_placement(5, Some([x, "0", "1"]), Some(["1", "0", "0"])).into_bytes();
    let file = read_step(&bytes).unwrap();
    let Surface::Plane(top) = &file.solids[0].faces[5].surface else {
        unreachable!();
    };
    assert_eq!(top.axis[0], Rational::parse(x).unwrap());
    bytes
}

fn displaced_carrier(z: &str, width: &str) -> Vec<u8> {
    let text = box_step(BoxFixture::new([width, "1", "1"]));
    let line = text
        .lines()
        .rfind(|line| line.contains("CARTESIAN_POINT("))
        .unwrap();
    let (prefix, _) = line.split_once(",(").unwrap();
    let changed = format!("{prefix},(0.,0.,{z}));");
    let source = text.replacen(line, &changed, 1).into_bytes();
    let file = read_step(&source).unwrap();
    let Surface::Plane(top) = &file.solids[0].faces[5].surface else {
        unreachable!()
    };
    assert_eq!(top.origin[2], Rational::parse(z).unwrap());
    source
}

fn warped_corner(z: &str) -> Vec<u8> {
    move_vertex(
        box_step(BoxFixture::new(["1", "1", "1"])),
        ["1", "1", "1"],
        z,
    )
}

fn move_vertex(text: String, point: [&str; 3], z: &str) -> Vec<u8> {
    let file = read_step(text.as_bytes()).unwrap();
    let point = point.map(|v| Rational::parse(v).unwrap());
    let name = &file.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start, &e.end])
        .find(|v| v.point == point)
        .unwrap()
        .name;
    let vertex = text
        .lines()
        .find(|line| line.starts_with(&format!("{name}=")))
        .unwrap();
    let point = vertex.rsplit('#').next().unwrap().trim_end_matches(");");
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{point}=")))
        .unwrap();
    let (prefix, _) = line.split_once(",(").unwrap();
    text.replacen(line, &format!("{prefix},(1.,1.,{z}));"), 1)
        .into_bytes()
}

fn repaired(source: &[u8]) -> catallax_import::Imported {
    apply(
        &mut Model::new(),
        source,
        "control.step",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap()
}

#[test]
fn a_small_carrier_discrepancy_is_repaired_with_a_face_record() {
    let source = tilted_top("0.000001");
    let mut model = Model::new();
    let options = ImportOptions::default();
    let imported = apply(&mut model, &source, "control.step", &options, &[], None)
        .expect("a bounded planar source discrepancy should be repairable");
    assert_eq!(imported.bodies.len(), 1);
    assert!(model.check(imported.bodies[0].1).is_empty());
    assert_eq!(imported.record.source.sha256, sha256_hex(&source));
    assert!(!imported.record.changes.is_empty());
    assert!(imported.derivative.is_some());
    assert_eq!(imported.record.format, "catallax-import-repair/3");
    assert_eq!(imported.record.source_incidence.len(), 1);
    let source_bound = &imported.record.source_incidence[0];
    assert_eq!(
        source_bound.source_definition,
        "projected_vertex_boundary/1"
    );
    assert_eq!(source_bound.method, PlanarRepairMethod::RefitPlanes);
    assert_eq!(source_bound.vertex_movement_upper_mm, "0");
    assert!(!imported.record.face_tolerances.is_empty());
    for face in &imported.record.face_tolerances {
        let upper = Rational::parse(&face.deviation_upper_mm).unwrap();
        assert!(upper > Rational::zero());
        assert!(upper <= Rational::parse("0.01").unwrap());
    }
    let reopened = rebuild(&mut Model::new(), &source, &options, &imported.record, None).unwrap();
    assert_eq!(reopened.record, imported.record);
    assert_eq!(reopened.derivative, imported.derivative);
}

#[test]
fn a_clean_source_keeps_its_existing_record_and_geometry() {
    let source = box_step(BoxFixture::new(["1", "1", "1"]));
    let imported = apply(
        &mut Model::new(),
        source.as_bytes(),
        "control.step",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap();
    assert_eq!(imported.record.format, "catallax-import-repair/2");
    assert!(imported.record.changes.is_empty());
    assert!(imported.record.face_tolerances.is_empty());
    assert!(imported.derivative.is_none());
    assert!(
        serde_json::to_value(&imported.record)
            .unwrap()
            .get("source_incidence")
            .is_none()
    );
}

#[test]
fn the_exact_ten_micrometre_boundary_is_allowed() {
    let imported = repaired(&displaced_carrier("1.01", "1"));
    assert_eq!(imported.record.face_tolerances.len(), 1);
    assert_eq!(
        imported.record.face_tolerances[0].deviation_upper_mm,
        "0.01"
    );
}

#[test]
fn a_bound_above_ten_micrometres_is_a_measured_refusal() {
    let source = displaced_carrier("1.010000000001", "1");
    let error = apply(
        &mut Model::new(),
        &source,
        "control.step",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap_err();
    let ImportRefusal::PlanarRepairNotVerified { detail } = error else {
        panic!("{error:?}")
    };
    assert_eq!(detail.stage, "import_cap");
    assert_eq!(detail.measurement, "whole_face_source_bound");
    assert_eq!(detail.limit_mm, "0.01");
    assert!(
        Rational::parse(detail.deviation_upper_mm.as_ref().unwrap()).unwrap()
            > Rational::parse("0.01").unwrap()
    );
}

#[test]
fn the_nearest_feature_limit_also_applies() {
    let source = displaced_carrier("1.003", "0.02");
    let error = apply(
        &mut Model::new(),
        &source,
        "control.step",
        &ImportOptions::default(),
        &[],
        None,
    )
    .unwrap_err();
    let ImportRefusal::PlanarRepairNotVerified { detail } = error else {
        panic!("{error:?}")
    };
    assert_eq!(detail.stage, "nearest_feature");
    assert_eq!(
        detail.feature_size_squared_mm2,
        Some(Rational::parse("0.0004").unwrap().to_string())
    );
}

#[test]
fn a_noncoplanar_shared_corner_is_snapped_once_and_recorded() {
    let source = warped_corner("1.000001");
    let before = read_step(&source).unwrap();
    let top = &before.solids[0].faces[5].name;
    let imported = repaired(&source);
    assert_eq!(imported.bodies.len(), 1);
    assert_eq!(imported.record.face_tolerances.len(), 3);
    assert!(
        imported
            .record
            .source_incidence
            .iter()
            .all(|r| r.method == PlanarRepairMethod::SnapSharedVertices)
    );
    let top = imported
        .record
        .source_incidence
        .iter()
        .find(|r| r.face == *top)
        .unwrap();
    assert_eq!(top.geometric_deviation_upper_mm, "0");
    assert!(Rational::parse(&top.source_to_result_upper_mm).unwrap() > Rational::zero());
    assert!(Rational::parse(&top.vertex_movement_upper_mm).unwrap() > Rational::zero());
    let rebuilt = rebuild(
        &mut Model::new(),
        &source,
        &ImportOptions::default(),
        &imported.record,
        None,
    )
    .unwrap();
    assert_eq!(rebuilt.record, imported.record);
}

#[test]
fn an_oblique_snap_retains_exact_fractional_bounds_and_rebuilds() {
    let mut rings = common::box_rings(["0", "0", "0"], ["1", "1", "1"]);
    for point in rings.iter_mut().flatten() {
        if point[0] == "1" && point[2] == "1" {
            point[2] = "1.5";
        }
    }
    let text = common::polyhedra_step(&[rings]);
    let placement = text
        .lines()
        .rfind(|line| line.contains("AXIS2_PLACEMENT_3D("))
        .unwrap();
    let args = placement.split_once("AXIS2_PLACEMENT_3D(").unwrap().1;
    let axis = args.split(',').nth(2).unwrap();
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{axis}=")))
        .unwrap();
    let text = text.replacen(line, &format!("{axis}=DIRECTION('',(-1.,0.,2.000001));"), 1);
    let source = move_vertex(text, ["1", "1", "1.5"], "1.500001");
    let imported = repaired(&source);
    assert!(
        imported
            .record
            .source_incidence
            .iter()
            .all(|r| r.method == PlanarRepairMethod::SnapSharedVertices)
    );
    assert!(
        imported.record.changes.iter().any(|c| c
            .bounds_after
            .iter()
            .flatten()
            .any(|v| v.contains('/')))
    );
    let rebuilt = rebuild(
        &mut Model::new(),
        &source,
        &ImportOptions::default(),
        &imported.record,
        None,
    )
    .unwrap();
    assert_eq!(rebuilt.record, imported.record);
    assert_eq!(rebuilt.derivative, imported.derivative);
}

#[test]
fn a_declined_incidence_repair_leaves_the_solid_unbuilt_and_unchanged() {
    let source = tilted_top("0.000001");
    let options = ImportOptions::default();
    let proposals = propose(&source, &options).unwrap();
    assert_eq!(proposals.format, "catallax-import-proposals/2");
    assert_eq!(proposals.proposals.len(), 1);
    let proposal = &proposals.proposals[0];
    assert_eq!(proposal.kind, GapKind::PlanarIncidence);
    assert_eq!(proposal.disposition, Disposition::Auto);
    let done = apply(
        &mut Model::new(),
        &source,
        "control.step",
        &options,
        &[Decision::Decline {
            id: proposal.id.clone(),
        }],
        None,
    )
    .unwrap();
    assert!(done.bodies.is_empty());
    assert_eq!(done.left_open.len(), 1);
    assert!(done.derivative.is_none());
    assert!(done.record.face_tolerances.is_empty());
    assert!(done.record.source_incidence.is_empty());
    assert!(
        done.record
            .changes
            .iter()
            .all(|c| c.approval == Approval::Declined)
    );
    assert_eq!(
        rebuild(&mut Model::new(), &source, &options, &done.record, None)
            .unwrap()
            .record,
        done.record
    );
}

#[test]
fn forged_source_bounds_are_rejected_without_mutating_the_model() {
    let source = tilted_top("0.000001");
    let mut record = repaired(&source).record;
    record.source_incidence[0].source_to_result_upper_mm = "0".into();
    let mut model = Model::new();
    let before = format!("{model:?}");
    let result = rebuild(
        &mut model,
        &source,
        &ImportOptions::default(),
        &record,
        None,
    );
    assert!(matches!(
        result,
        Err(ImportRefusal::ValidationFailed { .. })
    ));
    assert_eq!(format!("{model:?}"), before);
}

#[test]
fn the_last_unit_of_work_cannot_publish_a_partial_model() {
    let source = tilted_top("0.000001");
    let options = ImportOptions::default();
    let (result, full) = with_operation_budget(20_000_000, || {
        apply(
            &mut Model::new(),
            &source,
            "control.step",
            &options,
            &[],
            None,
        )
    });
    assert!(result.is_ok() && !full.exhausted && full.used > 0);
    let mut model = Model::new();
    let before = format!("{model:?}");
    let (result, work) = with_operation_budget(full.used - 1, || {
        apply(&mut model, &source, "control.step", &options, &[], None)
    });
    assert!(work.exhausted);
    assert!(matches!(result, Err(ImportRefusal::ResourceLimit)));
    assert_eq!(format!("{model:?}"), before);
}
