//! Preserve shared prototype identity and distinct placements through STEP export.

mod common;

use catallax_import::{
    geometry_source::{Geometry, ParsedStep, PlacementLink},
    sha256_hex,
    tolerant::opened,
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use common::{BoxFixture, r};

struct RepeatedFixture {
    source: String,
    prototype: u64,
    base: u64,
}

fn repeated_fixture(unit: common::Unit, with_orphan: bool) -> RepeatedFixture {
    let mut prototype_box = BoxFixture::new(["2", "3", "4"]);
    prototype_box.unit = unit;
    let mut boxes = vec![prototype_box];
    if with_orphan {
        let mut orphan = BoxFixture::new(["2", "3", "4"]);
        orphan.unit = unit;
        orphan.origin = [r("20"), r("0"), r("0")];
        boxes.push(orphan);
    }
    let mut source = common::boxes_step(&boxes);
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(source.as_bytes(), &mut dec).unwrap();
    let models = parsed.model_declarations(&mut dec).unwrap();
    let prototype = models[0].entity;
    let context = parsed
        .entities()
        .find(|&id| {
            parsed
                .entity_types(id)
                .unwrap()
                .contains(&"GEOMETRIC_REPRESENTATION_CONTEXT")
        })
        .unwrap();
    let base = parsed.entities().last().unwrap() + 1;
    let id = |offset| base + offset;
    let extra = format!(
        r#"
#{origin}=CARTESIAN_POINT('map origin',(0.,0.,0.));
#{x}=DIRECTION('x',(1.,0.,0.));
#{y}=DIRECTION('y',(0.,1.,0.));
#{z}=DIRECTION('z',(0.,0.,1.));
#{frame}=AXIS2_PLACEMENT_3D('map frame',#{origin},#{z},#{x});
#{part}=ADVANCED_BREP_SHAPE_REPRESENTATION('unit part',(#{prototype},#{frame}),#{context});
#{map}=REPRESENTATION_MAP(#{frame},#{part});
#{translated}=CARTESIAN_POINT('translated',(10.,0.,0.));
#{turned}=CARTESIAN_POINT('turned',(0.,20.,0.));
#{negative_x}=DIRECTION('negative x',(-1.,0.,0.));
#{identity}=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','identity',#{x},#{y},#{origin},1.,#{z});
#{translation}=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','translation',#{x},#{y},#{translated},1.,#{z});
#{rotation}=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','quarter turn',#{y},#{negative_x},#{turned},1.,#{z});
#{first}=MAPPED_ITEM('at origin',#{map},#{identity});
#{second}=MAPPED_ITEM('translated',#{map},#{translation});
#{third}=MAPPED_ITEM('turned',#{map},#{rotation});
#{assembly}=SHAPE_REPRESENTATION('three instances',(#{first},#{second},#{third}),#{context});
#{application}=APPLICATION_CONTEXT('automotive_design');
#{product_context}=PRODUCT_CONTEXT('',#{application},'mechanical');
#{product}=PRODUCT('repeated-control','repeated-control','',(#{product_context}));
#{formation}=PRODUCT_DEFINITION_FORMATION('','',#{product});
#{definition_context}=PRODUCT_DEFINITION_CONTEXT('part definition',#{application},'design');
#{definition}=PRODUCT_DEFINITION('design','',#{formation},#{definition_context});
#{shape}=PRODUCT_DEFINITION_SHAPE('','',#{definition});
#{binding}=SHAPE_DEFINITION_REPRESENTATION(#{shape},#{assembly});
"#,
        origin = id(0),
        x = id(1),
        y = id(2),
        z = id(3),
        frame = id(4),
        part = id(5),
        map = id(6),
        translated = id(7),
        turned = id(8),
        negative_x = id(9),
        identity = id(10),
        translation = id(11),
        rotation = id(12),
        first = id(13),
        second = id(14),
        third = id(15),
        assembly = id(16),
        application = id(17),
        product_context = id(18),
        product = id(19),
        formation = id(20),
        definition_context = id(21),
        definition = id(22),
        shape = id(23),
        binding = id(24)
    );
    source.insert_str(source.rfind("ENDSEC;").unwrap(), &extra);
    if with_orphan {
        let orphan = format!(
            "#{}=ADVANCED_BREP_SHAPE_REPRESENTATION('detached',(#{}),#{context});\n",
            id(25),
            models[1].entity
        );
        source.insert_str(source.rfind("ENDSEC;").unwrap(), &orphan);
    }
    RepeatedFixture {
        source,
        prototype,
        base,
    }
}

fn assert_inventory(parsed: &ParsedStep, fixture: &RepeatedFixture, geometry_root: u64) {
    let mut dec = Decider::default();
    let id = |offset| fixture.base + offset;
    let prototypes = parsed.model_declarations(&mut dec).unwrap();
    assert_eq!(
        prototypes
            .iter()
            .map(|model| model.entity)
            .collect::<Vec<_>>(),
        [geometry_root]
    );
    if geometry_root != fixture.prototype {
        let PlacementLink::MappedItem(map, transform) =
            parsed.placement_link(fixture.prototype, &mut dec).unwrap()
        else {
            panic!("selected prototype mapping missing")
        };
        let PlacementLink::RepresentationMap(_, representation) =
            parsed.placement_link(map, &mut dec).unwrap()
        else {
            panic!("selected prototype representation missing")
        };
        assert!(
            parsed
                .representation(representation)
                .unwrap()
                .items
                .contains(&geometry_root)
        );
        let Geometry::Transformation(mapping) = parsed.geometry(transform, &mut dec).unwrap()
        else {
            panic!("selected prototype unit mapping missing")
        };
        assert_eq!(mapping.scales[0], Some(r("1")));
        assert_eq!(
            parsed.geometry(mapping.origin, &mut dec).unwrap(),
            Geometry::Point(vec![r("0"); 3])
        );
    }
    assert_eq!(
        parsed.representation(id(5)).unwrap().items,
        [fixture.prototype, id(4)]
    );
    assert_eq!(
        parsed.representation(id(16)).unwrap().items,
        [id(13), id(14), id(15)]
    );
    assert_eq!(
        parsed.placement_link(id(6), &mut dec).unwrap(),
        PlacementLink::RepresentationMap(id(4), id(5))
    );
    assert_eq!(
        parsed.geometry(id(4), &mut dec).unwrap(),
        Geometry::Axis2Placement3d(id(0), Some(id(3)), Some(id(1)))
    );
    for (offset, values) in [
        (0, ["0", "0", "0"]),
        (7, ["10", "0", "0"]),
        (8, ["0", "20", "0"]),
    ] {
        assert_eq!(
            parsed.geometry(id(offset), &mut dec).unwrap(),
            Geometry::Point(values.map(r).to_vec())
        );
    }
    for (offset, values) in [
        (1, ["1", "0", "0"]),
        (2, ["0", "1", "0"]),
        (3, ["0", "0", "1"]),
        (9, ["-1", "0", "0"]),
    ] {
        assert_eq!(
            parsed.geometry(id(offset), &mut dec).unwrap(),
            Geometry::Direction(values.map(r).to_vec())
        );
    }
    for index in 0..3 {
        assert_eq!(
            parsed.placement_link(id(13 + index), &mut dec).unwrap(),
            PlacementLink::MappedItem(id(6), id(10 + index))
        );
        let Geometry::Transformation(transform) =
            parsed.geometry(id(10 + index), &mut dec).unwrap()
        else {
            panic!("placement transform missing")
        };
        assert_eq!(transform.dimension, Some(3));
        assert_eq!(transform.origin, id([0, 7, 8][index as usize]));
        let axes = if index == 2 {
            [id(2), id(9), id(3)]
        } else {
            [id(1), id(2), id(3)]
        };
        assert_eq!(transform.axes, axes.map(Some));
        assert_eq!(transform.scales[0], Some(r("1")));
        assert!(
            transform
                .scales
                .iter()
                .all(|scale| scale.as_ref().is_none_or(|value| value == &r("1")))
        );
    }
}

#[test]
fn whole_document_roundtrip_preserves_shared_prototype_and_distinct_instances() {
    for unit in [common::Unit::Millimetre, common::Unit::Metre] {
        let fixture = repeated_fixture(unit, false);
        let bytes = fixture.source.as_bytes();
        let hash = sha256_hex(bytes);
        let (result, work) = with_operation_budget(20_000_000, || {
            opened::open_step(
                bytes,
                &hash,
                &StepId::new("repeated-placement").unwrap(),
                Default::default(),
                &mut Decider::default(),
            )
        });
        assert!(!work.exhausted);
        let opened = result.unwrap();
        assert_eq!(opened.source().bodies().len(), 1);
        assert!(!matches!(
            opened.body(0),
            Some(opened::BodyView::Surfaces { .. }) | None
        ));
        assert_eq!(
            opened.source().bodies()[0].name(),
            format!("#{}", fixture.prototype)
        );
        assert_inventory(opened.parsed_source().unwrap(), &fixture, fixture.prototype);
        let (result, work) = with_operation_budget(20_000_000, || {
            opened.export_step_document(&mut Decider::default())
        });
        assert!(!work.exhausted);
        let exported = result.unwrap();
        assert_eq!(exported.bodies().len(), 1);
        let binding = &exported.bodies()[0];
        assert_eq!(binding.source_entity, fixture.prototype);
        assert_eq!(binding.source_index, 0);
        assert_eq!(binding.roots.len(), 1);
        let parsed = ParsedStep::read(exported.bytes(), &mut Decider::default()).unwrap();
        assert_inventory(&parsed, &fixture, binding.roots[0]);
        let (result, work) = with_operation_budget(20_000_000, || {
            opened::open_step(
                exported.bytes(),
                &sha256_hex(exported.bytes()),
                &StepId::new("reopened-placement").unwrap(),
                Default::default(),
                &mut Decider::default(),
            )
        });
        assert!(!work.exhausted);
        let reopened = result.unwrap();
        assert_eq!(reopened.bodies().len(), 1);
        assert!(!matches!(
            reopened.body(0),
            Some(opened::BodyView::Surfaces { .. }) | None
        ));
        assert_inventory(
            reopened.parsed_source().unwrap(),
            &fixture,
            binding.roots[0],
        );
        let (result, work) = with_operation_budget(20_000_000, || {
            reopened.export_step_document(&mut Decider::default())
        });
        assert!(!work.exhausted);
        let roundtrip = result.unwrap();
        assert_eq!(roundtrip.bodies().len(), 1);
        assert_eq!(roundtrip.bodies()[0].counts, [6, 12, 8]);
    }
}

#[test]
fn incomplete_source_document_refuses_without_exporting_a_sibling() {
    let complete = BoxFixture::new(["2", "3", "4"]);
    let mut incomplete = BoxFixture::new(["2", "3", "4"]);
    incomplete.omit = Some(5);
    let source = common::boxes_step(&[complete, incomplete]);
    let (result, work) = with_operation_budget(20_000_000, || {
        opened::open_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("incomplete-document").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    let opened = result.unwrap();
    assert!(!matches!(
        opened.body(0),
        Some(opened::BodyView::Surfaces { .. }) | None
    ));
    assert!(matches!(
        opened.body(1),
        Some(opened::BodyView::Surfaces { .. })
    ));
    let refusal = opened
        .export_step_document(&mut Decider::default())
        .unwrap_err();
    assert_eq!(refusal.code, "import.selected-document");
}

#[test]
fn an_orphan_representation_does_not_establish_source_document_coverage() {
    let fixture = repeated_fixture(common::Unit::Millimetre, true);
    let (result, work) = with_operation_budget(20_000_000, || {
        opened::open_step(
            fixture.source.as_bytes(),
            &sha256_hex(fixture.source.as_bytes()),
            &StepId::new("orphan-document").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    let opened = result.unwrap();
    assert_eq!(opened.bodies().len(), 2);
    for index in 0..2 {
        assert!(!matches!(
            opened.body(index),
            Some(opened::BodyView::Surfaces { .. }) | None
        ));
    }
    let refusal = opened
        .export_step_document(&mut Decider::default())
        .unwrap_err();
    assert_eq!(refusal.code, "import.selected-document");
}

#[test]
fn a_prototype_bound_to_conflicting_length_contexts_still_refuses() {
    let fixture = repeated_fixture(common::Unit::Millimetre, false);
    let id = |offset| fixture.base + offset;
    let extra = format!(
        "#{}=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\n\
         #{}=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#{}))REPRESENTATION_CONTEXT('',''));\n\
         #{}=ADVANCED_BREP_SHAPE_REPRESENTATION('conflicting units',(#{}),#{});\n\
         #{}=SHAPE_DEFINITION_REPRESENTATION(#{},#{});\n",
        id(25),
        id(26),
        id(25),
        id(27),
        fixture.prototype,
        id(26),
        id(28),
        id(23),
        id(27),
    );
    let mut source = fixture.source;
    source.insert_str(source.rfind("ENDSEC;").unwrap(), &extra);
    let (result, work) = with_operation_budget(20_000_000, || {
        opened::open_step(
            source.as_bytes(),
            &sha256_hex(source.as_bytes()),
            &StepId::new("conflicting-length-contexts").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    assert!(matches!(
        result,
        Err(catallax_import::tolerant::Error::Source(
            catallax_import::ImportRefusal::UnknownUnit { .. }
        ))
    ));
}
