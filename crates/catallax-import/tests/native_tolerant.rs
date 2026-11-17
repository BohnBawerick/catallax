//! Native admission must not inherit proof gates or invent material closure.
mod common;

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
    let bytes = import.export_body(body, &mut dec).unwrap();
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
    let mut second = first;
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
    let bytes = native.export_body(0, &mut dec).unwrap();
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
