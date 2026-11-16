//! New source definitions enter the existing per-body selector once.
mod common;
use catallax_import::{
    ImportRefusal, read_step, sha256_hex,
    tolerant::{
        self,
        opened::{self, BodyView, Tier},
    },
};
use catallax_numbers::Decider;
use catallax_recipe::StepId;
use catallax_topology::BodyKind;

const SQUARE: &str = include_str!("fixtures/surface-square.step");
fn open(source: &str) -> opened::OpenedStep {
    opened::open_step(
        source.as_bytes(),
        &sha256_hex(source.as_bytes()),
        &StepId::new("adapter").unwrap(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
}
#[test]
fn face_model_polygon_reaches_the_existing_selector_and_exports_a_surface() {
    let source = SQUARE
        .replace(
            "EDGE_LOOP('',(#24,#25,#26,#27))",
            "POLY_LOOP('',(#1,#2,#3,#4))",
        )
        .replace("ADVANCED_FACE('',(#29),#9,.T.)", "FACE('',(#29))")
        .replace("OPEN_SHELL('',(#30))", "CONNECTED_FACE_SET('',(#30))")
        .replace("SHELL_BASED_SURFACE_MODEL", "FACE_BASED_SURFACE_MODEL");
    assert!(matches!(
        read_step(source.as_bytes()),
        Err(ImportRefusal::Unsupported { .. })
    ));
    let opened = open(&source);
    assert_eq!(opened.bodies().len(), 1);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    let Some(BodyView::Certified { model, body }) = opened.body(0) else {
        panic!("finite polygon must build as a surface: {opened:?}");
    };
    assert_eq!(model.body(body).kind, BodyKind::Surface);
    let mut bytes = Vec::new();
    catallax_step::write_step(
        model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert_eq!(read_step(&bytes).unwrap().surface_models.len(), 1);
    let replayed = opened::reopen(
        source.as_bytes(),
        &opened.saved_controls(),
        &StepId::new("adapter").unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(replayed.source(), opened.source());
}
#[test]
fn spindle_torus_is_retained_complete_and_never_replaced_by_its_outer_sheet() {
    let source = SQUARE.replace("PLANE('',#8)", "TOROIDAL_SURFACE('',#8,1.,2.)");
    assert!(matches!(
        read_step(source.as_bytes()),
        Err(ImportRefusal::Unsupported { .. })
    ));
    let opened = open(&source);
    assert_eq!(opened.bodies().len(), 1);
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    let Some(BodyView::Surfaces { source, .. }) = opened.body(0) else {
        panic!("invalid trims must remain source surfaces")
    };
    let catallax_import::Surface::Retained(carrier) = &source.shells()[0].faces[0].surface else {
        panic!("whole source torus required")
    };
    let catallax_import::geometry_source::Geometry::Surface(
        catallax_import::geometry_source::Surface::Torus(_, a, b),
    ) = carrier
        .source()
        .geometry(carrier.entity(), &mut Decider::default())
        .unwrap()
    else {
        panic!("torus source declaration required")
    };
    assert_eq!((a.to_string(), b.to_string()), ("1".into(), "2".into()));
    assert_eq!(source.edges().len(), 4);
}
#[test]
fn oriented_inner_shells_are_retained_and_not_silently_discarded() {
    let source = SQUARE.replace(
        "#32=SHELL_BASED_SURFACE_MODEL('',(#31));",
        "#32=BREP_WITH_VOIDS('',#31,(#36));\n#36=ORIENTED_CLOSED_SHELL('',*,#31,.F.);",
    );
    let opened = open(&source);
    assert_eq!(opened.bodies().len(), 1);
    assert_eq!(opened.source().bodies()[0].shells().len(), 2);
    let shells = opened.source().bodies()[0].shells();
    assert_ne!(shells[0].faces[0].same_sense, shells[1].faces[0].same_sense);
    assert!(matches!(
        opened.source().definition(&mut Decider::default()),
        Err(tolerant::Error::Source(ImportRefusal::NotSolid { .. }))
    ));
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
}
#[test]
fn complex_face_and_edge_inheritance_preserve_the_source_inventory() {
    let source = SQUARE
        .replace(
            "ADVANCED_FACE('',(#29),#9,.T.)",
            "(ADVANCED_FACE() FACE((#29)) FACE_SURFACE(#9,.T.) REPRESENTATION_ITEM(''))",
        )
        .replace(
            "EDGE_CURVE('',#10,#11,#16,.T.)",
            "(EDGE(#10,#11) EDGE_CURVE(#16,.T.) REPRESENTATION_ITEM(''))",
        );
    let opened = open(&source);
    assert_eq!(opened.bodies().len(), 1);
    assert_eq!(opened.source().bodies()[0].edges().len(), 4);
    assert!(matches!(opened.body(0), Some(BodyView::Certified { .. })));
}
