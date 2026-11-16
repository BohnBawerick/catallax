//! Geometric measurements and whole-candidate publication, not ledger-only repairs.
mod common;

use catallax_import::{
    sha256_hex,
    tolerant::{
        self, ToleranceOptions,
        healing::{self, HealingChange},
    },
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_topology::{Sense, numerical as n};
use common::{BoxFixture, box_step, r};

fn source(shift: &str) -> String {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0"), r("0"), r(shift)]));
    box_step(fixture)
}
fn read(text: &str) -> tolerant::TolerantFile {
    tolerant::read_step(
        text.as_bytes(),
        &sha256_hex(text.as_bytes()),
        ToleranceOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn a_shared_edge_source_is_sewn_with_measured_geometry_and_an_approximate_result() {
    let text = source("0.001");
    let input = read(&text);
    let before = input.clone();
    assert_eq!(input.bodies()[0].edges().len(), 16);
    assert_eq!(input.bodies()[0].vertices().len(), 12);
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert_eq!(healed.source(), &before);
    assert!(healed.report().approximate);
    assert_eq!(
        healed
            .report()
            .changes
            .iter()
            .filter(|c| matches!(c, HealingChange::EdgesSewn { .. }))
            .count(),
        4
    );
    assert_eq!(
        healed
            .report()
            .changes
            .iter()
            .filter(|c| matches!(c, HealingChange::VerticesJoined { .. }))
            .count(),
        4
    );
    assert_eq!(
        healed
            .report()
            .stages
            .iter()
            .map(|s| s.stage.as_str())
            .collect::<Vec<_>>(),
        ["sewing", "wire", "face", "shell", "solid"]
    );
    assert!(!healed.report().measurements.is_empty());
    let result = healed.definition(&mut Decider::default()).unwrap();
    let body = checked_planar_fixture(&result);
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
    let record = healed.tolerances();
    assert!(!record.bodies[0].changes.is_empty());
    assert!(
        record.bodies[0]
            .changes
            .iter()
            .all(|c| r(&c.after_mm) <= r("0.01"))
    );
    assert!(before.bodies()[0].tolerance_changes().is_empty());
}

// Independent consumer of the derived fixture through the public numerical
// B-rep builder. Retain the supplied 3D line carriers: an exact constructor
// would correctly refuse their nonzero, explicitly measured incidence gaps.
fn checked_planar_fixture(file: &catallax_import::StepFile) -> n::NumericalBody {
    use catallax_shapes::{Curve, Line, Plane, RPoint3, RVector3, Surface};
    use std::collections::BTreeMap;
    let dot = |a: &[catallax_numbers::Rational; 3], b: &[catallax_numbers::Rational; 3]| {
        (0..3).fold(r("0"), |sum, i| sum + &a[i] * &b[i])
    };
    let mut definition = n::NumericalDefinition {
        name: "measured-shared-edge-candidate".into(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let solid = &file.solids[0];
    for face in &solid.faces {
        let catallax_import::Surface::Plane(plane) = &face.surface else {
            panic!("authored planar fixture")
        };
        let y = [
            &plane.axis[1] * &plane.ref_dir[2] - &plane.axis[2] * &plane.ref_dir[1],
            &plane.axis[2] * &plane.ref_dir[0] - &plane.axis[0] * &plane.ref_dir[2],
            &plane.axis[0] * &plane.ref_dir[1] - &plane.axis[1] * &plane.ref_dir[0],
        ];
        assert_eq!(dot(&plane.ref_dir, &plane.ref_dir), r("1"));
        assert_eq!(dot(&y, &y), r("1"));
        let uv = |point: &[catallax_numbers::Rational; 3]| {
            let delta = std::array::from_fn(|i| &point[i] - &plane.origin[i]);
            [dot(&delta, &plane.ref_dir), dot(&delta, &y)]
        };
        let mut loops = Vec::new();
        for boundary in &face.loops {
            let mut uses = boundary.edges.clone();
            if !boundary.orientation {
                uses.reverse();
                for edge in &mut uses {
                    edge.forward = !edge.forward;
                }
            }
            let mut trims = Vec::new();
            for edge in &uses {
                let mut indices = Vec::new();
                for vertex in [&edge.start, &edge.end] {
                    let index = *vertices.entry(vertex.name.clone()).or_insert_with(|| {
                        let i = definition.vertices.len();
                        definition.vertices.push(n::NumericalVertex {
                            key: vertex.name.clone(),
                            position: vertex.point.clone(),
                        });
                        i
                    });
                    indices.push(index);
                }
                let index = *edges.entry(edge.edge.clone()).or_insert_with(|| {
                    let catallax_import::Curve::Line(origin, direction) = &edge.curve else {
                        panic!("authored line fixture")
                    };
                    let parameter = |p: &[catallax_numbers::Rational; 3]| {
                        let delta = std::array::from_fn(|i| &p[i] - &origin[i]);
                        dot(&delta, direction) / dot(direction, direction)
                    };
                    let i = definition.edges.len();
                    definition.edges.push(n::NumericalEdge {
                        key: edge.edge.clone(),
                        curve: Curve::Line(Line {
                            origin: RPoint3::new(
                                origin[0].clone(),
                                origin[1].clone(),
                                origin[2].clone(),
                            )
                            .to_real(),
                            dir: RVector3::new(
                                direction[0].clone(),
                                direction[1].clone(),
                                direction[2].clone(),
                            )
                            .to_real(),
                        }),
                        vertices: [indices[0], indices[1]],
                        parameters: [parameter(&edge.start.point), parameter(&edge.end.point)],
                    });
                    i
                });
                let (a, b) = if edge.forward {
                    (&edge.start.point, &edge.end.point)
                } else {
                    (&edge.end.point, &edge.start.point)
                };
                trims.push(n::NumericalTrim {
                    edge: index,
                    sense: if edge.forward {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                    uv: [uv(a), uv(b)],
                });
            }
            loops.push(trims);
        }
        definition.faces.push(n::NumericalFace {
            key: face.name.clone(),
            surface: Surface::Plane(Plane {
                origin: RPoint3::new(
                    plane.origin[0].clone(),
                    plane.origin[1].clone(),
                    plane.origin[2].clone(),
                )
                .to_real(),
                normal: RVector3::new(
                    plane.axis[0].clone(),
                    plane.axis[1].clone(),
                    plane.axis[2].clone(),
                )
                .to_real(),
                ref_dir: RVector3::new(
                    plane.ref_dir[0].clone(),
                    plane.ref_dir[1].clone(),
                    plane.ref_dir[2].clone(),
                )
                .to_real(),
            }),
            sense: if face.same_sense {
                Sense::Same
            } else {
                Sense::Reversed
            },
            loops,
        });
    }
    n::build(
        definition,
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn a_late_geometry_cap_failure_publishes_neither_earlier_sewing_nor_growth() {
    let mut text = source("0.001");
    let line = text
        .lines()
        .find(|l| l.contains("=LINE("))
        .unwrap()
        .to_owned();
    let origin = line.split(",#").nth(1).unwrap().split(',').next().unwrap();
    let changed = line.replacen(&format!(",#{origin}"), ",#900001", 1);
    text = text.replacen(&line, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#900001=CARTESIAN_POINT('',(100.,100.,100.));\nENDSEC;\nEND-ISO-10303-21;",
        1,
    );
    let input = read(&text);
    let before = input.clone();
    let error = healing::heal(&input, &mut Decider::default()).unwrap_err();
    assert!(matches!(error,tolerant::Error::Numerical(ref r) if r.code=="import.healing-cap"));
    assert_eq!(input, before);
    assert!(input.bodies()[0].tolerance_changes().is_empty());
}

#[test]
fn over_cap_open_geometry_is_retained_with_explicit_shell_issues() {
    let input = read(&source("0.02"));
    let before = input.clone();
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert!(
        !healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::EdgesSewn { .. }))
    );
    assert!(
        healed
            .report()
            .issues
            .iter()
            .any(|i| i.code == "shell.open-or-nonmanifold")
    );
}

#[test]
fn an_exhausted_candidate_preserves_the_original_source_and_ledger() {
    let input = read(&source("0.001"));
    let before = input.clone();
    let (complete, measured) = with_operation_budget(20_000_000, || {
        healing::heal(&input, &mut Decider::default())
    });
    complete.unwrap();
    for allowance in [0, 1, measured.used / 2, measured.used - 1] {
        let (result, work) =
            with_operation_budget(allowance, || healing::heal(&input, &mut Decider::default()));
        assert!(result.is_err() && work.exhausted);
        assert_eq!(input, before);
    }
}
