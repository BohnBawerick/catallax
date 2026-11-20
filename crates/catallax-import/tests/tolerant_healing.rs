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

fn native_open(text: &str) -> tolerant::opened::OpenedStep {
    let (result, work) = with_operation_budget(20_000_000, || {
        tolerant::opened::open_step(
            text.as_bytes(),
            &sha256_hex(text.as_bytes()),
            &catallax_recipe::StepId::new("wire-control").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    result.unwrap()
}

#[test]
fn coincident_circle_endpoints_sew_only_the_same_finite_arc() {
    for sense in [false, true] {
        for reverse in [false, true] {
            // The two faces share a chord, but own distinct circle edges.
            // Opposite carrier senses with the same endpoints select the
            // complementary quarter/three-quarter arcs, despite zero gap.
            let (start, end) = if reverse { (13, 12) } else { (12, 13) };
            let text = format!(
                "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
                 #1=CARTESIAN_POINT('',(0.,0.,0.));\n\
                 #2=DIRECTION('',(0.,0.,1.));\n\
                 #3=DIRECTION('',(1.,0.,0.));\n\
                 #4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
                 #5=CIRCLE('',#4,1.);\n#6=PLANE('',#4);\n\
                 #10=CARTESIAN_POINT('',(1.,0.,0.));\n\
                 #11=CARTESIAN_POINT('',(0.,1.,0.));\n\
                 #12=VERTEX_POINT('',#10);\n#13=VERTEX_POINT('',#11);\n\
                 #14=DIRECTION('',(-1.,1.,0.));\n\
                 #15=VECTOR('',#14,1.);\n#16=LINE('',#10,#15);\n\
                 #20=EDGE_CURVE('',#12,#13,#5,.T.);\n\
                 #21=EDGE_CURVE('',#{start},#{end},#5,.{}.);\n\
                 #22=EDGE_CURVE('',#12,#13,#16,.T.);\n\
                 #30=ORIENTED_EDGE('',*,*,#20,.T.);\n\
                 #31=ORIENTED_EDGE('',*,*,#22,.F.);\n\
                 #32=ORIENTED_EDGE('',*,*,#21,.{}.);\n\
                 #33=ORIENTED_EDGE('',*,*,#22,.F.);\n\
                 #40=EDGE_LOOP('',(#30,#31));\n\
                 #41=FACE_OUTER_BOUND('',#40,.T.);\n\
                 #42=ADVANCED_FACE('',(#41),#6,.T.);\n\
                 #43=EDGE_LOOP('',(#32,#33));\n\
                 #44=FACE_OUTER_BOUND('',#43,.T.);\n\
                 #45=ADVANCED_FACE('',(#44),#6,.T.);\n\
                 #50=CLOSED_SHELL('',(#42,#45));\n\
                 #51=MANIFOLD_SOLID_BREP('',#50);\n\
                 #60=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n\
                 #61=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#60)) REPRESENTATION_CONTEXT('',''));\n\
                 ENDSEC;\nEND-ISO-10303-21;\n",
                if sense { 'T' } else { 'F' },
                if reverse { 'F' } else { 'T' }
            );
            let input = read(&text);
            let before = input.clone();
            let healed = healing::heal(&input, &mut Decider::default()).unwrap();
            assert_eq!(input, before);
            let sewn = healed.report().changes.iter().any(|change| matches!(change,
                HealingChange::EdgesSewn { retained, replaced, .. }
                if [retained.as_str(), replaced.as_str()].iter().any(|s| ["#20", "#21"].contains(s))));
            assert_eq!(sewn, sense != reverse);
            let definition = healed.definition(&mut Decider::default()).unwrap();
            let arcs: Vec<_> = definition
                .named_faces()
                .flat_map(|(_, face)| &face.loops)
                .flat_map(|boundary| &boundary.edges)
                .filter(|edge| matches!(edge.curve, catallax_import::Curve::Circle(..)))
                .collect();
            assert_eq!(arcs.len(), 2);
            assert_eq!(arcs[0].edge, "#20");
            assert_eq!(arcs[1].edge, if sewn { "#20" } else { "#21" });
            assert!(arcs[0].curve_same_sense);
            assert_eq!(arcs[1].curve_same_sense, if sewn { true } else { sense });
        }
    }
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

#[test]
fn slightly_skew_line_carriers_sew_by_their_complete_finite_segments() {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let text = box_step(fixture);
    let baseline = healing::heal(&read(&text), &mut Decider::default()).unwrap();
    let replaced = baseline
        .report()
        .changes
        .iter()
        .find_map(|change| match change {
            HealingChange::EdgesSewn { replaced, .. } => Some(replaced.as_str()),
            _ => None,
        })
        .unwrap();
    let declaration = |name: &str| {
        text.lines()
            .find(|line| line.starts_with(&format!("{name}=")))
            .unwrap()
    };
    let last_reference = |line: &str| {
        format!(
            "#{}",
            line.rsplit('#')
                .next()
                .unwrap()
                .split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
        )
    };
    let carrier = last_reference(declaration(replaced));
    let vector = last_reference(declaration(&carrier));
    let direction = last_reference(declaration(&vector));
    let original = declaration(&direction);
    let text = text.replace(original, &original.replacen("0.", "0.0001", 1));
    let input = read(&text);
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(
        healed
            .report()
            .changes
            .iter()
            .filter(|change| matches!(change, HealingChange::EdgesSewn { .. }))
            .count(),
        4,
    );
    let body = checked_planar_fixture(&healed.definition(&mut Decider::default()).unwrap());
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
    let exported =
        catallax_step::numerical::build(&body, &Default::default(), &mut Decider::default())
            .unwrap();
    exported
        .check_written(exported.bytes(), &mut Decider::default())
        .unwrap();
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
                    tag: catallax_topology::numerical::CopyTag::Independent,
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
        let (result, work) = with_operation_budget(allowance, || {
            healing::attempt(&input, &mut Decider::default())
        });
        let failure = result.unwrap_err();
        assert!(work.exhausted);
        assert!(
            matches!(&failure.error, tolerant::Error::Numerical(r) if r.code == "numbers.operation-budget")
        );
        assert_eq!(
            failure.preparation_work_used
                + failure.stages.iter().map(|s| s.work_used).sum::<u64>()
                + failure.publication_work_used,
            work.used
        );
        assert_eq!(input, before);
    }
}

#[test]
fn an_existing_wire_is_reordered_without_losing_any_edge_use() {
    let mut text = box_step(BoxFixture::new(["2", "3", "4"]));
    let line = text
        .lines()
        .find(|l| l.contains("=EDGE_LOOP("))
        .unwrap()
        .to_owned();
    let (prefix, tail) = line.split_once(",(").unwrap();
    let mut refs = tail
        .split(')')
        .next()
        .unwrap()
        .split(',')
        .collect::<Vec<_>>();
    refs.swap(1, 2);
    text = text.replacen(&line, &format!("{prefix},({}));", refs.join(",")), 1);
    let input = read(&text);
    let before = input.clone();
    let result = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert_eq!(
        result
            .report()
            .changes
            .iter()
            .filter(|c| matches!(c, HealingChange::WireReordered { .. }))
            .count(),
        1
    );
    let body = checked_planar_fixture(&result.definition(&mut Decider::default()).unwrap());
    assert_eq!(
        (
            body.counts().faces,
            body.counts().edges,
            body.counts().vertices
        ),
        (6, 12, 8)
    );
    let opened = native_open(&text);
    assert_ne!(
        opened.bodies()[0].tier(),
        tolerant::opened::Tier::SurfaceBody
    );
    let repairs = match opened.body(0).unwrap() {
        tolerant::opened::BodyView::Tolerant { import, .. } => &import.evidence().repairs,
        _ => {
            &opened.bodies()[0]
                .healing()
                .expect("the checked admission must retain its wire-repair history")
                .report()
                .changes
        }
    };
    assert!(
        repairs
            .iter()
            .any(|change| matches!(change, HealingChange::WireReordered { .. }))
    );
    match opened.body(0).unwrap() {
        tolerant::opened::BodyView::Numerical(body) => {
            assert_eq!((body.counts().faces, body.counts().edges), (6, 12));
        }
        tolerant::opened::BodyView::Tolerant { import, body } => {
            let counts = import.bodies()[body].counts();
            assert_eq!((counts.faces, counts.edges), (6, 12));
        }
        tolerant::opened::BodyView::Certified { model, body } => {
            assert_eq!(model.body(body).kind, catallax_topology::BodyKind::Solid);
        }
        tolerant::opened::BodyView::Surfaces { .. } => panic!("wire repair must retain a solid"),
    }
}

fn reverse_last_flag(line: &str) -> String {
    if let Some(prefix) = line.strip_suffix(".T.);") {
        format!("{prefix}.F.);")
    } else if let Some(prefix) = line.strip_suffix(".F.);") {
        format!("{prefix}.T.);")
    } else {
        panic!("authored boolean flag")
    }
}

#[test]
fn a_reversed_face_boundary_gets_an_explicit_orientation_repair() {
    let mut text = box_step(BoxFixture::new(["2", "3", "4"]));
    let line = text
        .lines()
        .find(|l| l.contains("=FACE_OUTER_BOUND("))
        .unwrap()
        .to_owned();
    text = text.replacen(&line, &reverse_last_flag(&line), 1);
    let input = read(&text);
    let result = healing::heal(&input, &mut Decider::default()).unwrap();
    assert!(
        result
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::BoundaryReversed { .. }))
    );
    checked_planar_fixture(&result.definition(&mut Decider::default()).unwrap());
    assert!(input.bodies()[0].tolerance_changes().is_empty());
}

#[test]
fn an_inside_out_polyhedral_shell_is_oriented_outwards_as_one_candidate() {
    let text = box_step(BoxFixture::new(["2", "3", "4"]))
        .lines()
        .map(|line| {
            if line.contains("=FACE_OUTER_BOUND(") || line.contains("=ADVANCED_FACE(") {
                reverse_last_flag(line)
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let input = read(&text);
    let result = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(result.report().changes.iter().filter(|c|matches!(c,HealingChange::FaceReversed{cause,..} if cause=="solid-positive-oriented-polyhedral-volume")).count(),6);
    checked_planar_fixture(&result.definition(&mut Decider::default()).unwrap());
    assert!(result.report().approximate);
}

#[test]
fn a_true_crossing_remains_an_unresolved_surface_condition() {
    let mut soup = common::Soup::new();
    let ring = [
        ["0", "0", "0"],
        ["2", "2", "0"],
        ["0", "2", "0"],
        ["2", "0", "0"],
    ]
    .map(|p| p.map(r));
    soup.plane_face(&ring, &[], 2, true, &[r("0"), r("0"), r("0")]);
    let input = read(&soup.finish());
    let before = input.clone();
    let result = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert!(
        result
            .report()
            .issues
            .iter()
            .any(|i| i.code == "face.crossing-or-collinear-nonadjacent-boundaries")
    );
    assert!(
        !result
            .report()
            .checks
            .iter()
            .any(|c| c.check == "solid.positive-oriented-polyhedral-volume")
    );
}

#[test]
fn sewing_keeps_independent_source_bodies_and_their_bindings_separate() {
    let mut first = BoxFixture::new(["2", "3", "4"]);
    first.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let text = common::boxes_step(&[first, BoxFixture::new(["2", "3", "4"])]);
    let input = read(&text);
    let original = input.definition(&mut Decider::default()).unwrap();
    let result = healing::heal(&input, &mut Decider::default()).unwrap();
    let definition = result.definition(&mut Decider::default()).unwrap();
    assert_eq!(definition.solids.len(), 2);
    assert_eq!(definition.solids[1], original.solids[1]);
    let unchanged = input.bodies()[1].name();
    assert!(
        !result
            .report()
            .changes
            .iter()
            .any(|change| change.body() == unchanged)
    );
}

#[test]
fn normal_shared_edge_sewing_exports_a_checked_approximate_body() {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let text = box_step(fixture);
    let input = read(&text);
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    let body = checked_planar_fixture(&healed.definition(&mut Decider::default()).unwrap());
    let step = catallax_step::numerical::build(&body, &Default::default(), &mut Decider::default())
        .unwrap();
    step.check_written(step.bytes(), &mut Decider::default())
        .unwrap();
    assert!(healed.report().approximate);
    assert_eq!(step.report().counts, [6, 12, 8]);
    if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
        let folder = std::path::Path::new(&folder);
        std::fs::write(folder.join("normal-shared-edge.step"), step.bytes()).unwrap();
        let record = serde_json::json!({"format":"catallax-healing-geometry-control/1","source_sha256":input.source_sha256(),"healing":healed.report(),"tolerances":healed.tolerances(),"step_sha256":sha256_hex(step.bytes()),"numerical_step":true,"complete_source_body_bound_established":false});
        std::fs::write(
            folder.join("normal-shared-edge.json"),
            serde_json::to_vec_pretty(&record).unwrap(),
        )
        .unwrap();
    }
}

#[test]
fn a_replaced_source_carrier_cannot_escape_its_original_geometry_cap() {
    fn refs(value: &str) -> Vec<usize> {
        value
            .split('#')
            .skip(1)
            .map(|s| {
                s[..s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len())]
                    .parse()
                    .unwrap()
            })
            .collect()
    }
    let mut text = source("0.001");
    let input = read(&text);
    let baseline = healing::heal(&input, &mut Decider::default()).unwrap();
    let replaced = baseline
        .report()
        .changes
        .iter()
        .find_map(|c| match c {
            HealingChange::EdgesSewn { replaced, .. } => Some(replaced),
            _ => None,
        })
        .unwrap();
    let edge = text
        .lines()
        .find(|l| l.starts_with(&format!("{replaced}=")))
        .unwrap();
    let carrier = refs(edge.split_once('=').unwrap().1)[2];
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{carrier}=")))
        .unwrap()
        .to_owned();
    let origin = refs(line.split_once('=').unwrap().1)[0];
    let replacement = line.replacen(&format!(",#{origin}"), ",#900002", 1);
    text = text.replacen(&line, &replacement, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#900002=CARTESIAN_POINT('',(100.,100.,100.));\nENDSEC;\nEND-ISO-10303-21;",
        1,
    );
    let input = read(&text);
    let before = input.clone();
    let error = healing::heal(&input, &mut Decider::default()).unwrap_err();
    assert!(matches!(error,tolerant::Error::Numerical(ref r) if r.code=="import.healing-cap"));
    assert_eq!(input, before);
}

fn split_wire_endpoint(gap: &str) -> String {
    let mut soup = common::Soup::new();
    let ring = [
        ["0", "0", "0"],
        ["2", "0", "0"],
        ["2", "2", "0"],
        ["0", "2", "0"],
    ]
    .map(|p| p.map(r));
    soup.plane_face(&ring, &[], 2, true, &[r("0"), r("0"), r("0")]);
    let mut text = soup.finish();
    let file = catallax_import::read_step(text.as_bytes()).unwrap();
    let edge = &file.solids[0].faces[0].loops[0].edges[0];
    assert!(edge.forward);
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("{}=", edge.edge)))
        .unwrap()
        .to_owned();
    let mut fields: Vec<_> = line.split(',').map(str::to_owned).collect();
    fields[2] = "#900004".into();
    let x = (&r("2") + &r(gap)).to_decimal_exact().unwrap();
    let points = format!(
        "#900003=CARTESIAN_POINT('',({x},0.,0.));\n#900004=VERTEX_POINT('',#900003);\nENDSEC;\nEND-ISO-10303-21;"
    );
    text = text.replacen(&line, &fields.join(","), 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &points,
        1,
    );
    text
}

#[test]
fn an_ordered_wire_gap_is_joined_from_measured_endpoints_without_exact_incidence() {
    let input = read(&split_wire_endpoint("0.000001"));
    let before = input.clone();
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert!(healed.report().approximate);
    assert!(
        healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::VerticesJoined { .. }))
    );
    assert!(
        healed
            .report()
            .measurements
            .iter()
            .any(|m| m.check == "ordered-adjacent-source-endpoint-gap")
    );
    assert!(
        !healed
            .report()
            .issues
            .iter()
            .any(|i| i.code.starts_with("wire.open"))
    );
    assert!(
        healed
            .report()
            .issues
            .iter()
            .any(|i| i.code == "shell.open-or-nonmanifold")
    );
    assert!(!healed.tolerances().bodies[0].changes.is_empty());
    let opened = native_open(&split_wire_endpoint("0.000001"));
    assert_eq!(
        opened.bodies()[0].first_blocker(),
        Some(tolerant::admission::NativeStage::ShellIncidence)
    );
}

#[test]
fn an_over_cap_wire_gap_stays_unresolved_and_does_not_grow_tolerances() {
    let input = read(&split_wire_endpoint("0.02"));
    let before = input.clone();
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert!(
        healed
            .report()
            .issues
            .iter()
            .any(|i| i.code.starts_with("wire.open"))
    );
    assert!(healed.tolerances().bodies[0].changes.is_empty());
    assert!(
        !healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::VerticesJoined { .. }))
    );
    let opened = native_open(&split_wire_endpoint("0.02"));
    assert_eq!(
        opened.bodies()[0].first_blocker(),
        Some(tolerant::admission::NativeStage::WireConnection)
    );
    assert!(
        opened
            .tolerant_admission()
            .unwrap()
            .evidence()
            .repairs
            .is_empty()
    );
}

fn split_box_edge(fraction: &str) -> String {
    let text = box_step(BoxFixture::new(["2", "3", "4"]));
    let file = catallax_import::read_step(text.as_bytes()).unwrap();
    let edge = &file.solids[0].faces[0].loops[0].edges[0];
    let point: [catallax_numbers::Rational; 3] = std::array::from_fn(|i| {
        &edge.start.point[i] + &((&edge.end.point[i] - &edge.start.point[i]) * r(fraction))
    });
    let edge_line = text
        .lines()
        .find(|l| l.starts_with(&format!("{}=", edge.edge)))
        .unwrap();
    let mut fields: Vec<_> = edge_line.split(',').map(str::to_owned).collect();
    let end = fields[2].clone();
    fields[2] = "#900004".into();
    let mut additions = vec![
        format!(
            "#900003=CARTESIAN_POINT('',({},{},{}));",
            common::dec(&point[0]),
            common::dec(&point[1]),
            common::dec(&point[2])
        ),
        "#900004=VERTEX_POINT('',#900003);".into(),
        format!(
            "#900006=EDGE_CURVE('',#900004,{end},{},{}",
            fields[3], fields[4]
        ),
    ];
    let mut uses = std::collections::BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| l.contains("=ORIENTED_EDGE(") && l.contains(&format!(",{},", edge.edge)))
    {
        let (id, _) = line.split_once('=').unwrap();
        let forward = line.ends_with(".T.);");
        let next = format!("#{}", 900010 + uses.len());
        additions.push(format!(
            "{next}=ORIENTED_EDGE('',*,*,#900006,{});",
            if forward { ".T." } else { ".F." }
        ));
        uses.insert(id.to_owned(), (next, forward));
    }
    let text = text
        .lines()
        .map(|line| {
            if line == edge_line {
                return fields.join(",");
            }
            if line.contains("=EDGE_LOOP(") {
                let (prefix, tail) = line.split_once(",(").unwrap();
                let mut expanded = Vec::new();
                for reference in tail.split(')').next().unwrap().split(',') {
                    if let Some((next, forward)) = uses.get(reference) {
                        if *forward {
                            expanded.extend([reference.to_owned(), next.clone()]);
                        } else {
                            expanded.extend([next.clone(), reference.to_owned()]);
                        }
                    } else {
                        expanded.push(reference.to_owned());
                    }
                }
                return format!("{prefix},({}));", expanded.join(","));
            }
            line.to_owned()
        })
        .collect::<Vec<_>>()
        .join("\n");
    text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{}\nENDSEC;\nEND-ISO-10303-21;", additions.join("\n")),
        1,
    )
}

#[test]
fn only_finite_lines_below_the_initial_uncertainty_are_removed() {
    for (index, fraction) in ["0", "0.000000001"].iter().enumerate() {
        let input = read(&split_box_edge(fraction));
        assert_eq!(
            (
                input.bodies()[0].edges().len(),
                input.bodies()[0].vertices().len()
            ),
            (13, 9)
        );
        let healed = healing::heal(&input, &mut Decider::default()).unwrap();
        assert_eq!(
            healed
                .report()
                .changes
                .iter()
                .filter(|c| matches!(c, HealingChange::SmallLineRemoved { .. }))
                .count(),
            1
        );
        assert!(healed.report().approximate);
        let body = checked_planar_fixture(&healed.definition(&mut Decider::default()).unwrap());
        assert_eq!((body.counts().edges, body.counts().vertices), (12, 8));
        let step =
            catallax_step::numerical::build(&body, &Default::default(), &mut Decider::default())
                .unwrap();
        step.check_written(step.bytes(), &mut Decider::default())
            .unwrap();
        if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&folder).join(format!("small-line-{index}.step")),
                step.bytes(),
            )
            .unwrap();
        }
    }
    let input = read(&split_box_edge("0.000001"));
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert!(
        !healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::SmallLineRemoved { .. }))
    );
    let body = checked_planar_fixture(&healed.definition(&mut Decider::default()).unwrap());
    assert_eq!((body.counts().edges, body.counts().vertices), (13, 9));
    let step = catallax_step::numerical::build(&body, &Default::default(), &mut Decider::default())
        .unwrap();
    if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&folder).join("preserved-line.step"),
            step.bytes(),
        )
        .unwrap();
    }
}

#[test]
fn a_closed_circle_is_not_a_small_line_because_its_endpoint_chord_is_zero() {
    let mut soup = common::Soup::new();
    soup.cylinder(["0", "0", "0"], 2, "2", "1", common::Caps::Both);
    let input = read(&soup.finish());
    let original = input.definition(&mut Decider::default()).unwrap();
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert!(
        !healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::SmallLineRemoved { .. }))
    );
    assert_eq!(
        healed.definition(&mut Decider::default()).unwrap(),
        original
    );
}
