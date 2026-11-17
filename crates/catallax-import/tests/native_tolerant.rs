//! Native admission must not inherit proof gates or invent material closure.
mod common;
#[path = "common/shared_curve.rs"]
mod shared_curve;

use catallax_import::{
    geometry_source::{ParsedStep, Topology},
    sha256_hex,
    tolerant::{
        self,
        opened::{self, BodyView, Tier},
    },
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use common::{BoxFixture, box_step, r};

fn open(source: &str) -> opened::OpenedStep {
    let (result, work) = with_operation_budget(20_000_000, || {
        opened::open_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("native-case").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(
        !work.exhausted,
        "the regression must exercise geometry, not a resource cutoff"
    );
    result.unwrap()
}

#[test]
fn native_entity_tolerances_absorb_shared_vertex_error_beyond_the_checked_cap() {
    let original = box_step(BoxFixture::new(["2", "3", "4"]));
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(original.as_bytes(), &mut dec).unwrap();
    let vertex = parsed
        .entities()
        .find(|&entity| {
            parsed
                .entity_types(entity)
                .unwrap()
                .contains(&"VERTEX_POINT")
        })
        .unwrap();
    let Topology::Vertex(point) = parsed.topology(vertex, &mut dec).unwrap() else {
        unreachable!()
    };
    let prefix = format!("#{point}=");
    let source = original
        .lines()
        .map(|line| {
            if line.starts_with(&prefix) {
                format!("#{point}=CARTESIAN_POINT('',(0.1,0.,0.));")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let opened = open(&source);
    assert_eq!(opened.bodies()[0].tier(), Tier::TolerantSolid);
    let BodyView::Tolerant { import, body } = opened.body(0).unwrap() else {
        panic!("measured native geometry required")
    };
    let solid = &import.bodies()[body];
    let counts = solid.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
    assert!(
        solid
            .definition()
            .vertices()
            .iter()
            .any(|vertex| vertex.tolerance_mm() > &r("0.01"))
    );
    assert!(
        import
            .evidence()
            .sampled_checks
            .iter()
            .any(|check| check.maximum_residual_mm > 0.09)
    );
    let controls = opened.saved_controls();
    let replayed = opened::reopen(
        source.as_bytes(),
        &controls,
        &StepId::new("native-case").unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(replayed.bodies()[0].tier(), Tier::TolerantSolid);
    let BodyView::Tolerant { import: replay, .. } = replayed.body(0).unwrap() else {
        unreachable!()
    };
    assert_eq!(
        serde_json::to_value(import.evidence()).unwrap(),
        serde_json::to_value(replay.evidence()).unwrap()
    );
}

const TUBE: &str = "ISO-10303-21;\nHEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,1.));\n\
#3=CARTESIAN_POINT('',(1.,0.,0.));#4=CARTESIAN_POINT('',(1.,0.,1.));\n\
#5=CARTESIAN_POINT('',(0.,1.,0.));#6=CARTESIAN_POINT('',(0.,1.,1.));\n\
#7=CARTESIAN_POINT('',(0.,0.,0.));#8=CARTESIAN_POINT('',(0.,0.,1.));\n\
#9=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4),(#5,#6),(#7,#8)),.UNSPECIFIED.,.T.,.F.,.F.,(2,1,1,2),(2,2),(0.,1.,2.,3.),(0.,1.),.UNSPECIFIED.);\n\
#10=ADVANCED_FACE('',(),#9,.T.);#11=CLOSED_SHELL('',(#10));#12=MANIFOLD_SOLID_BREP('',#11);\n\
#90=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
#91=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#90))REPRESENTATION_CONTEXT('',''));\n\
ENDSEC;END-ISO-10303-21;";

#[test]
fn declared_closed_spline_directions_do_not_cap_an_edgeless_open_tube() {
    for source in [
        TUBE.to_owned(),
        TUBE.replace(".UNSPECIFIED.,.T.,.F.,.F.,", ".UNSPECIFIED.,.T.,.T.,.F.,"),
    ] {
        let opened = open(&source);
        assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
        let BodyView::Surfaces { source, .. } = opened.body(0).unwrap() else {
            panic!("open tube must remain retained geometry")
        };
        let source = source.materialize(&mut Decider::default()).unwrap();
        assert_eq!(source.shells()[0].faces.len(), 1);
        assert!(
            matches!(opened.bodies()[0].refusal(), Some(tolerant::Error::Numerical(error)) if error.code == "import.tolerant-broken-topology")
        );
    }
}

#[test]
fn sewing_does_not_match_bowed_spline_edges_by_endpoints_alone() {
    let source = "ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(1.,0.,0.));\n\
#3=CARTESIAN_POINT('',(0.5,1.,0.));#4=CARTESIAN_POINT('',(0.5,-1.,0.));\n\
#5=DIRECTION('',(1.,0.,0.));#6=DIRECTION('',(-1.,0.,0.));#7=DIRECTION('',(0.,0.,1.));\n\
#8=VECTOR('',#5,1.);#9=VECTOR('',#6,1.);#10=AXIS2_PLACEMENT_3D('',#1,#7,#5);#11=PLANE('',#10);\n\
#12=B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#3,#2),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);\n\
#13=B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#4,#2),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);\n\
#14=LINE('',#2,#9);#15=LINE('',#1,#8);#16=VERTEX_POINT('',#1);#17=VERTEX_POINT('',#2);\n\
#18=EDGE_CURVE('',#16,#17,#12,.T.);#19=EDGE_CURVE('',#17,#16,#14,.T.);\n\
#20=EDGE_CURVE('',#17,#16,#13,.F.);#21=EDGE_CURVE('',#16,#17,#15,.T.);\n\
#22=ORIENTED_EDGE('',*,*,#18,.T.);#23=ORIENTED_EDGE('',*,*,#19,.T.);\n\
#24=ORIENTED_EDGE('',*,*,#20,.T.);#25=ORIENTED_EDGE('',*,*,#21,.T.);\n\
#26=EDGE_LOOP('',(#22,#23));#27=EDGE_LOOP('',(#24,#25));#28=FACE_OUTER_BOUND('',#26,.T.);#29=FACE_OUTER_BOUND('',#27,.T.);\n\
#30=ADVANCED_FACE('',(#28),#11,.T.);#31=ADVANCED_FACE('',(#29),#11,.F.);\n\
#32=CLOSED_SHELL('',(#30,#31));#33=MANIFOLD_SOLID_BREP('',#32);\n\
#90=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
#91=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#90))REPRESENTATION_CONTEXT('',''));ENDSEC;END-ISO-10303-21;";
    let opened = open(source);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    let BodyView::Surfaces { source, .. } = opened.body(0).unwrap() else {
        panic!("different finite spline images must not become a closed solid")
    };
    let source = source.materialize(&mut Decider::default()).unwrap();
    assert_eq!(source.shells()[0].faces.len(), 2);
    assert_eq!(source.edges().len(), 4);
}

#[test]
fn exported_ellipse_extrusion_has_a_continuous_side_uv_wire() {
    use catallax_import::geometry_source::{AngleUnit, Curve, Geometry, Surface};
    use catallax_import::tolerant::native::{NativeLocation, NativeUvCurve};
    use std::collections::BTreeMap;

    let opened = open(common::elliptical_extrusion::SOURCE);
    let BodyView::Tolerant { import, body } = opened.body(0).unwrap() else {
        panic!(
            "elliptical extrusion must reach native admission: {:?}",
            opened.bodies()[0].refusal()
        );
    };
    let mut dec = Decider::default();
    let bytes = import.export_body(body, &mut dec).unwrap().into_bytes();
    if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&folder).join("native-ellipse.step"),
            &bytes,
        )
        .unwrap();
    }
    let parsed = ParsedStep::read(&bytes, &mut dec).unwrap();
    let side = parsed
        .entities()
        .find(|&id| {
            matches!(
                parsed.geometry(id, &mut dec),
                Ok(Geometry::Surface(Surface::Extrusion(..)))
            )
        })
        .unwrap();
    let bounds = parsed
        .entities()
        .find_map(|id| match parsed.topology(id, &mut dec) {
            Ok(Topology::Face {
                surface: Some(surface),
                bounds,
                ..
            }) if surface == side => Some(bounds),
            _ => None,
        })
        .unwrap();
    let Topology::FaceBound { loop_id, .. } = parsed.topology(bounds[0], &mut dec).unwrap() else {
        unreachable!()
    };
    let Topology::EdgeLoop(uses) = parsed.topology(loop_id, &mut dec).unwrap() else {
        unreachable!()
    };
    let mut occurrences = BTreeMap::<u64, usize>::new();
    let mut segments = Vec::new();
    for edge_use in uses {
        let Topology::OrientedEdge(edge, forward) = parsed.topology(edge_use, &mut dec).unwrap()
        else {
            unreachable!()
        };
        let Topology::Edge {
            curve: Some(curve),
            same_sense: Some(sense),
            ..
        } = parsed.topology(edge, &mut dec).unwrap()
        else {
            unreachable!()
        };
        let Geometry::Curve(Curve::SurfaceCurve(_, associated, _, _)) =
            parsed.geometry(curve, &mut dec).unwrap()
        else {
            panic!("export must supply the side's curve-on-surface geometry");
        };
        let pcurves = associated
            .into_iter()
            .filter(|&id| {
                matches!(
                    parsed.geometry(id, &mut dec),
                    Ok(Geometry::Curve(Curve::Pcurve(surface, _))) if surface == side
                )
            })
            .collect::<Vec<_>>();
        let occurrence = occurrences.entry(edge).or_default();
        let pcurve = pcurves[*occurrence];
        *occurrence += 1;
        let charts = parsed
            .curve_charts(pcurve, &r("1"), AngleUnit::Radians, &mut dec)
            .unwrap();
        let limits = &charts[0].domain()[0];
        let carrier =
            NativeUvCurve::from_source(&parsed, pcurve, &r("1"), AngleUnit::Radians, &mut dec)
                .unwrap();
        let mut endpoints =
            [limits[0].as_ref().unwrap(), limits[1].as_ref().unwrap()].map(|limit| {
                carrier
                    .sample(
                        NativeLocation {
                            chart: 0,
                            parameters: [limit.to_f64_nearest(), 0.0],
                        },
                        &mut dec,
                    )
                    .unwrap()
                    .0
            });
        if forward != sense {
            endpoints.reverse();
        }
        segments.push(endpoints);
    }
    for index in 0..segments.len() {
        let end = segments[index][1];
        let next = segments[(index + 1) % segments.len()][0];
        assert!(
            (end[0] - next[0]).hypot(end[1] - next[1]) < 1e-5,
            "UV wire jumps sheets: {end:?} -> {next:?}"
        );
    }
}

#[test]
fn selected_body_step_report_does_not_disclose_sibling_geometry_or_repairs() {
    let mut first = BoxFixture::new(["2", "3", "4"]);
    first.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let mut second = BoxFixture::new(["2", "3", "4"]);
    second.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    second.origin = [r("20"), r("0"), r("0")];
    let source = common::boxes_step(&[first, second]);
    let mut dec = Decider::default();
    let native = tolerant::admission::import_step(
        source.as_bytes(),
        &sha256_hex(source.as_bytes()),
        Default::default(),
        &mut dec,
    )
    .unwrap();
    let selected = native.bodies()[0].name();
    let sibling = native.bodies()[1].name();
    assert!(
        native
            .evidence()
            .native_trims
            .iter()
            .any(|trim| trim.body == sibling)
    );
    assert!(
        serde_json::to_value(&native.evidence().repairs)
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .any(|repair| repair["body"] == sibling)
    );
    let bytes = native.export_body(0, &mut dec).unwrap().into_bytes();
    let parsed = ParsedStep::read(&bytes, &mut dec).unwrap();
    let roots = parsed
        .entities()
        .filter(|&entity| {
            matches!(
                parsed.topology(entity, &mut dec),
                Ok(Topology::Solid { .. } | Topology::SurfaceModel(_) | Topology::FaceModel(_))
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        roots,
        vec![selected.strip_prefix('#').unwrap().parse::<u64>().unwrap()],
        "body-local STEP exposed an unselected body root"
    );
    let text = std::str::from_utf8(&bytes).unwrap();
    let json = text
        .split("/* CATALLAX_REPORT_JSON ")
        .skip(1)
        .map(|part| part.split_once(" */").unwrap().0)
        .collect::<String>();
    let report: serde_json::Value = serde_json::from_str(&json).unwrap();
    let evidence = &report["source_document_evidence"];
    for field in [
        "tolerances",
        "sampled_checks",
        "repairs",
        "topology_checks",
        "topology_measurements",
        "warnings",
        "native_trims",
        "native_domains",
    ] {
        for record in evidence[field].as_array().unwrap() {
            if !record["body"].is_null() {
                assert_eq!(
                    record["body"], selected,
                    "sibling record escaped through {field}"
                );
            }
        }
    }
    assert!(
        evidence["native_trims"]
            .as_array()
            .unwrap()
            .iter()
            .any(|trim| trim["body"] == selected)
    );
    assert!(
        evidence["repairs"]
            .as_array()
            .unwrap()
            .iter()
            .any(|repair| repair["body"] == selected)
    );
    for check in report["export_checks"].as_array().unwrap() {
        assert_eq!(check["body"], selected);
    }
}

#[test]
fn a_finite_open_spline_cannot_become_closed_by_repeating_its_vertex() {
    let source = common::elliptical_extrusion::SOURCE.replace(
        "#12=ELLIPSE('opposite-axis lower ellipse',#9,2.,1.);",
        "#12=B_SPLINE_CURVE_WITH_KNOTS('',1,(#3,#100),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);",
    ).replace("ENDSEC;\nEND-ISO-10303-21;",
        "#100=CARTESIAN_POINT('',(8.,0.,0.));\nENDSEC;\nEND-ISO-10303-21;");
    let result = tolerant::admission::import_step(
        source.as_bytes(),
        &sha256_hex(source.as_bytes()),
        Default::default(),
        &mut Decider::default(),
    );
    assert!(matches!(&result, Err(tolerant::Error::Numerical(error))
        if error.code == "import.tolerant-broken-topology"));
}

#[test]
fn mixed_angular_contexts_refuse_before_the_default_checked_upgrade() {
    let source = common::in_degrees(common::elliptical_extrusion::SOURCE).replace(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#99994=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#42,#99990,#44))REPRESENTATION_CONTEXT('',''));\nENDSEC;\nEND-ISO-10303-21;",
    );
    let (result, work) = with_operation_budget(20_000_000, || {
        opened::open_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("mixed-units").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    assert!(matches!(&result, Err(tolerant::Error::Numerical(error))
        if error.code == "import.tolerant-ambiguous-units"));
}

#[test]
fn shared_curve_pcurves_stay_with_each_exported_body() {
    use catallax_import::geometry_source::{Curve as SourceCurve, Geometry};
    use std::collections::BTreeSet;

    let fixture = shared_curve::boxes();
    let mut dec = Decider::default();
    let native = tolerant::admission::import_step(
        fixture.source.as_bytes(),
        &sha256_hex(fixture.source.as_bytes()),
        Default::default(),
        &mut dec,
    )
    .unwrap();
    let selected = native.export_body(1, &mut dec).unwrap().into_bytes();
    let parsed = ParsedStep::read(&selected, &mut dec).unwrap();
    for surface in fixture.first_surfaces {
        assert!(
            parsed.entity_types(surface).is_err(),
            "a sibling-only carrier remained reachable from body-local STEP"
        );
    }
    let document = native.export_step(&mut dec).unwrap();
    let parsed = ParsedStep::read(document.bytes(), &mut dec).unwrap();
    for root in fixture.roots {
        for (edge, surfaces) in shared_curve::uses(&parsed, root, &mut dec) {
            let Topology::Edge {
                curve: Some(curve), ..
            } = parsed.topology(edge, &mut dec).unwrap()
            else {
                unreachable!()
            };
            let Geometry::Curve(SourceCurve::SurfaceCurve(_, associated, _, _)) =
                parsed.geometry(curve, &mut dec).unwrap()
            else {
                panic!("each body's edge needs its own surface associations");
            };
            let actual = associated
                .into_iter()
                .map(|id| match parsed.geometry(id, &mut dec).unwrap() {
                    Geometry::Curve(SourceCurve::Pcurve(surface, _)) => surface,
                    Geometry::Surface(_) => id,
                    _ => panic!("invalid curve-on-surface association"),
                })
                .collect::<BTreeSet<_>>();
            assert_eq!(
                actual,
                surfaces.into_iter().collect::<BTreeSet<_>>(),
                "the writer reused another body's pcurve associations"
            );
        }
    }
}

#[test]
fn a_bad_first_body_cannot_refuse_a_usable_later_sibling() {
    let mut later = BoxFixture::new(["2", "3", "4"]);
    later.origin[0] = r("20");
    let source = common::boxes_step(&[BoxFixture::new(["2", "3", "4"]), later]).replacen(
        "MANIFOLD_SOLID_BREP",
        "CSG_SOLID",
        1,
    );
    let opened = open(&source);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    assert_eq!(
        opened.bodies()[0].first_blocker(),
        Some(tolerant::admission::NativeStage::BodyBuild)
    );
    assert_eq!(opened.bodies()[1].tier(), Tier::ValidSolid);
    let native = opened.tolerant_admission().unwrap();
    assert_eq!(native.bodies().len(), 1);
    assert_eq!(native.bodies()[0].source_index(), 1);
    let name = native.bodies()[0].name();
    assert!(
        native
            .evidence()
            .tolerances
            .iter()
            .all(|row| row.body == name)
    );
    assert!(native.export_step(&mut Decider::default()).is_err());
    let file = native.export_body(0, &mut Decider::default()).unwrap();
    let parsed = ParsedStep::read(file.bytes(), &mut Decider::default()).unwrap();
    let declarations = parsed.model_declarations(&mut Decider::default()).unwrap();
    assert_eq!(declarations.len(), 1);
    assert_eq!(format!("#{}", declarations[0].entity), name);
}

#[test]
fn checked_upgrade_exhaustion_keeps_the_already_built_native_body() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let baseline = open(&source);
    let allowance = baseline.read_work_used() + baseline.native_admission_work_used() + 10;
    let (opened, work) = with_operation_budget(allowance, || {
        opened::open_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("upgrade").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let opened = opened.unwrap();
    assert!(work.exhausted);
    assert_eq!(opened.bodies()[0].tier(), Tier::TolerantSolid);
    assert!(opened.bodies()[0].upgrade_refusal().is_some());
    assert!(opened.bodies()[0].first_blocker().is_none());
    let BodyView::Tolerant { import, body } = opened.body(0).unwrap() else {
        panic!("optional proof must not discard the native body")
    };
    assert_eq!(import.bodies()[body].counts().faces, 6);
    import.export_body(body, &mut Decider::default()).unwrap();
}

#[test]
fn missing_sphere_pole_edges_are_reconstructed_from_measured_collapse() {
    let source = include_str!("fixtures/half-sphere-missing-poles.step");
    let (native, work) = with_operation_budget(20_000_000, || {
        tolerant::admission::import_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let native = native.unwrap();
    assert!(!work.exhausted);
    let inserted = native
        .evidence()
        .repairs
        .iter()
        .filter(|change| {
            matches!(
                change,
                tolerant::healing::HealingChange::DegenerateEdgeInserted { .. }
            )
        })
        .count();
    assert_eq!(inserted, 2);
    assert_eq!(native.bodies()[0].counts().edges, 5);
    let file = native.export_body(0, &mut Decider::default()).unwrap();
    if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&folder).join("native-poles.step"),
            file.bytes(),
        )
        .unwrap();
    }
}
