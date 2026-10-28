//! The Part 21 reader and analytic recognition through the public entry.

mod common;

use catallax_import::{Curve, ImportRefusal, Surface, read_step};
use common::*;

#[test]
fn reads_a_box_as_exact_planes_and_lines() {
    let file = read_step(box_step(BoxFixture::new(["10.5", "20", "0.001"])).as_bytes()).unwrap();
    assert_eq!(file.mm_per_unit, r("1"));
    let [solid] = file.solids.as_slice() else {
        panic!("one solid")
    };
    assert_eq!(solid.faces.len(), 6);
    for face in &solid.faces {
        assert!(matches!(face.surface, Surface::Plane(_)));
        assert_eq!(face.loops.len(), 1);
        assert!(face.loops[0].outer && face.loops[0].edges.len() == 4);
        assert!(
            face.loops[0]
                .edges
                .iter()
                .all(|e| matches!(e.curve, Curve::Line(..)))
        );
    }
    let xs: Vec<_> = solid
        .faces
        .iter()
        .flat_map(|f| f.loops[0].edges.iter().map(|e| e.start.point[2].clone()))
        .collect();
    assert!(xs.contains(&r("0.001")) && xs.contains(&r("0")));
    // Shared edges are the same EDGE_CURVE entity from both faces.
    let mut counts = std::collections::BTreeMap::new();
    for f in &solid.faces {
        for e in &f.loops[0].edges {
            *counts.entry(e.edge.clone()).or_insert(0) += 1;
        }
    }
    assert_eq!(counts.len(), 12);
    assert!(counts.values().all(|&n| n == 2));
}

#[test]
fn converts_metres_and_inches_exactly() {
    let mut m = BoxFixture::new(["1", "2", "3"]);
    m.unit = Unit::Metre;
    let file = read_step(box_step(m).as_bytes()).unwrap();
    assert_eq!(file.mm_per_unit, r("1000"));
    let max = file.solids[0]
        .faces
        .iter()
        .flat_map(|f| f.loops[0].edges.iter().flat_map(|e| e.start.point.clone()))
        .max()
        .unwrap();
    assert_eq!(max, r("3000"));
    let mut i = BoxFixture::new(["1", "2", "3"]);
    i.unit = Unit::Inch;
    let file = read_step(box_step(i).as_bytes()).unwrap();
    assert_eq!(file.mm_per_unit, r("25.4"));
    let max = file.solids[0]
        .faces
        .iter()
        .flat_map(|f| f.loops[0].edges.iter().flat_map(|e| e.start.point.clone()))
        .max()
        .unwrap();
    assert_eq!(max, r("76.2"));
}

#[test]
fn reads_two_solids() {
    let mut b = BoxFixture::new(["1", "1", "1"]);
    b.origin = [r("1.005"), r("0"), r("0")];
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), b]);
    assert_eq!(read_step(text.as_bytes()).unwrap().solids.len(), 2);
}

#[test]
fn refuses_with_the_entity_name() {
    let text = box_step(BoxFixture::new(["1", "1", "1"]));
    let swapped = text.replace("PLANE('',", "B_SPLINE_SURFACE('',");
    let refusal = read_step(swapped.as_bytes()).unwrap_err();
    let ImportRefusal::Unsupported { entity, type_name } = refusal else {
        panic!("expected unsupported, got {refusal:?}");
    };
    assert!(entity.starts_with('#'));
    assert_eq!(type_name, "B_SPLINE_SURFACE");
    assert!(matches!(
        read_step(b"not a step file").unwrap_err(),
        ImportRefusal::Malformed { .. }
    ));
}

#[test]
fn refuses_unknown_units_and_missing_units() {
    let text = box_step(BoxFixture::new(["1", "1", "1"]));
    let bad = text.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.KILO.,.METRE.)");
    assert!(matches!(
        read_step(bad.as_bytes()).unwrap_err(),
        ImportRefusal::UnknownUnit { .. }
    ));
    let none = text.replace("LENGTH_UNIT()", "AREA_UNIT()");
    assert!(matches!(
        read_step(none.as_bytes()).unwrap_err(),
        ImportRefusal::UnknownUnit { .. }
    ));
}

#[test]
fn quoted_text_and_comments_do_not_confuse_the_reader() {
    let text = box_step(BoxFixture::new(["1", "1", "1"]))
        .replace("'box'", "'a;b #1 ''q'' (x)'")
        .replace("DATA;", "DATA; /* a ; comment */");
    assert_eq!(read_step(text.as_bytes()).unwrap().solids[0].faces.len(), 6);
}
