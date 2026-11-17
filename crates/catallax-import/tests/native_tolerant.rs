//! Native admission must not inherit proof gates or invent material closure.
mod common;

use catallax_import::{geometry_source::{ParsedStep, Topology}, sha256_hex, tolerant::{self, opened::{self, BodyView, Tier}}};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use common::{BoxFixture, box_step, r};

fn open(source: &str) -> opened::OpenedStep {
    let (result, work) = with_operation_budget(20_000_000, || opened::open_step(
        source.as_bytes(), &sha256_hex(source.as_bytes()), &StepId::new("native-case").unwrap(),
        Default::default(), &mut Decider::default(),
    ));
    assert!(!work.exhausted, "the regression must exercise geometry, not a resource cutoff");
    result.unwrap()
}

#[test]
fn native_entity_tolerances_absorb_shared_vertex_error_beyond_the_checked_cap() {
    let original = box_step(BoxFixture::new(["2", "3", "4"]));
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(original.as_bytes(), &mut dec).unwrap();
    let vertex = parsed.entities().find(|&entity| parsed.entity_types(entity).unwrap().contains(&"VERTEX_POINT")).unwrap();
    let Topology::Vertex(point) = parsed.topology(vertex, &mut dec).unwrap() else { unreachable!() };
    let prefix = format!("#{point}=");
    let source = original.lines().map(|line| {
        if line.starts_with(&prefix) { format!("#{point}=CARTESIAN_POINT('',(0.1,0.,0.));") } else { line.to_owned() }
    }).collect::<Vec<_>>().join("\n");
    let opened = open(&source);
    assert_eq!(opened.bodies()[0].tier(), Tier::TolerantSolid);
    let BodyView::Tolerant { import, body } = opened.body(0).unwrap() else { panic!("measured native geometry required") };
    let solid = &import.bodies()[body];
    let counts = solid.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
    assert!(solid.definition().vertices().iter().any(|vertex| vertex.tolerance_mm() > &r("0.01")));
    assert!(import.evidence().sampled_checks.iter().any(|check| check.maximum_residual_mm > 0.09));
    let controls = opened.saved_controls();
    let replayed = opened::reopen(source.as_bytes(), &controls, &StepId::new("native-case").unwrap(), &mut Decider::default()).unwrap();
    assert_eq!(replayed.bodies()[0].tier(), Tier::TolerantSolid);
    let BodyView::Tolerant { import: replay, .. } = replayed.body(0).unwrap() else { unreachable!() };
    assert_eq!(serde_json::to_value(import.evidence()).unwrap(), serde_json::to_value(replay.evidence()).unwrap());
}

const TUBE: &str = "ISO-10303-21;\nHEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,1.));\n\
#3=CARTESIAN_POINT('',(1.,0.,0.));#4=CARTESIAN_POINT('',(1.,0.,1.));\n\
#5=CARTESIAN_POINT('',(0.,1.,0.));#6=CARTESIAN_POINT('',(0.,1.,1.));\n\
#7=CARTESIAN_POINT('',(0.,0.,0.));#8=CARTESIAN_POINT('',(0.,0.,1.));\n\
#9=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4),(#5,#6),(#7,#8)),.UNSPECIFIED.,.T.,.F.,.F.,(2,1,1,2),(2,2),(0.,1.,2.,3.),(0.,1.),.UNSPECIFIED.);\n\
#10=ADVANCED_FACE('',(),#9,.T.);#11=CLOSED_SHELL('',(#10));#12=MANIFOLD_SOLID_BREP('',#11);\n\
ENDSEC;END-ISO-10303-21;";

#[test]
fn one_closed_spline_direction_does_not_cap_an_edgeless_tube() {
    let opened = open(TUBE);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    let BodyView::Surfaces { source, .. } = opened.body(0).unwrap() else { panic!("open tube must remain retained geometry") };
    assert_eq!(source.shells()[0].faces.len(), 1);
    assert!(matches!(opened.bodies()[0].refusal(), Some(tolerant::Error::Numerical(error)) if error.code == "import.tolerant-broken-topology"));
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
#32=CLOSED_SHELL('',(#30,#31));#33=MANIFOLD_SOLID_BREP('',#32);ENDSEC;END-ISO-10303-21;";
    let opened = open(source);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    let BodyView::Surfaces { source, .. } = opened.body(0).unwrap() else { panic!("different finite spline images must not become a closed solid") };
    assert_eq!(source.shells()[0].faces.len(), 2);
    assert_eq!(source.edges().len(), 4);
}
