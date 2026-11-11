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
fn placement_projects_a_non_perpendicular_reference_exactly() {
    for (axis, reference, projected) in [
        (["2", "0", "0"], ["1", "3", "0"], ["0", "3", "0"]),
        (
            ["2", "3", "6"],
            ["1", "2", "3"],
            ["-3/49", "20/49", "-9/49"],
        ),
    ] {
        let text = box_with_placement(0, Some(axis), Some(reference));
        let file = read_step(text.as_bytes()).unwrap();
        let Surface::Plane(p) = &file.solids[0].faces[0].surface else {
            panic!("plane")
        };
        assert_eq!(p.axis, axis.map(r));
        assert_eq!(p.ref_dir, projected.map(r));
    }
}

#[test]
fn placement_defaults_each_optional_direction_independently() {
    for (axis, reference, expected_axis, expected_reference) in [
        (None, None, ["0", "0", "1"], ["1", "0", "0"]),
        (
            Some(["0", "0", "2"]),
            None,
            ["0", "0", "2"],
            ["1", "0", "0"],
        ),
        (
            Some(["2", "0", "0"]),
            None,
            ["2", "0", "0"],
            ["0", "1", "0"],
        ),
        (
            Some(["-2", "0", "0"]),
            None,
            ["-2", "0", "0"],
            ["0", "1", "0"],
        ),
        (
            Some(["1", "1", "0"]),
            None,
            ["1", "1", "0"],
            ["1/2", "-1/2", "0"],
        ),
        (
            None,
            Some(["2", "3", "4"]),
            ["0", "0", "1"],
            ["2", "3", "0"],
        ),
    ] {
        let text = box_with_placement(0, axis, reference);
        let file = read_step(text.as_bytes()).unwrap();
        let Surface::Plane(p) = &file.solids[0].faces[0].surface else {
            panic!("plane")
        };
        assert_eq!(p.axis, expected_axis.map(r));
        assert_eq!(p.ref_dir, expected_reference.map(r));
    }
}

#[test]
fn placement_refuses_parallel_directions_with_the_entity_name() {
    for (axis, reference) in [
        (Some(["2", "3", "6"]), ["4", "6", "12"]),
        (Some(["2", "3", "6"]), ["-4", "-6", "-12"]),
        (None, ["0", "0", "2"]),
    ] {
        let text = box_with_placement(0, axis, Some(reference));
        let line = text
            .lines()
            .find(|l| l.contains("AXIS2_PLACEMENT_3D("))
            .unwrap();
        let id = line.split_once('=').unwrap().0;
        assert_eq!(
            read_step(text.as_bytes()).unwrap_err(),
            ImportRefusal::Malformed {
                detail: format!("placement {id} axis and reference direction are parallel"),
            }
        );
    }
}

#[test]
fn placement_keeps_arbitrarily_small_nonzero_projections() {
    let text = box_with_placement(0, Some(["1", "0", "0"]), Some(["1", "1e-100", "0"]));
    let file = read_step(text.as_bytes()).unwrap();
    let Surface::Plane(p) = &file.solids[0].faces[0].surface else {
        panic!("plane")
    };
    assert_eq!(p.ref_dir, ["0", "1e-100", "0"].map(r));
}

#[test]
fn placement_keeps_perpendicular_direction_ratios_unchanged() {
    let axis = ["2", "3", "6"];
    let reference = ["3", "-2", "0"];
    let text = box_with_placement(0, Some(axis), Some(reference));
    let file = read_step(text.as_bytes()).unwrap();
    let Surface::Plane(p) = &file.solids[0].faces[0].surface else {
        panic!("plane")
    };
    assert_eq!(p.axis, axis.map(r));
    assert_eq!(p.ref_dir, reference.map(r));
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
    let bad = text.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.LIGHTYEAR.,.METRE.)");
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

/// The millimetre box with its length unit replaced. `bodies[0]` is the unit entity
/// and the rest are extra entities, each written after it; `@n` in any body is the
/// id of `bodies[n]`. The representation context assigns the entities in `assigned`.
fn box_with_length_units(bodies: &[&str], assigned: &[usize]) -> String {
    let text = box_step(BoxFixture::new(["10", "20", "30"]));
    let unit_line = text
        .lines()
        .find(|l| l.contains("SI_UNIT(.MILLI.,.METRE.)"))
        .expect("the box declares its unit");
    let first: u64 = unit_line[1..unit_line.find('=').unwrap()].parse().unwrap();
    let resolve = |body: &str| {
        let mut out = body.to_string();
        for n in (0..bodies.len()).rev() {
            out = out.replace(&format!("@{n}"), &(first + n as u64).to_string());
        }
        out
    };
    let lines: String = bodies
        .iter()
        .enumerate()
        .map(|(n, b)| format!("#{}={};\n", first + n as u64, resolve(b)))
        .collect();
    let list: Vec<String> = assigned
        .iter()
        .map(|n| format!("#{}", first + *n as u64))
        .collect();
    text.replace(&format!("{unit_line}\n"), &lines).replace(
        &format!("GLOBAL_UNIT_ASSIGNED_CONTEXT((#{first}))"),
        &format!("GLOBAL_UNIT_ASSIGNED_CONTEXT(({}))", list.join(",")),
    )
}

fn box_with_length_unit(bodies: &[&str]) -> String {
    box_with_length_units(bodies, &[0])
}

const CONVERTED: &str = "( CONVERSION_BASED_UNIT('METRE',#@1) LENGTH_UNIT() NAMED_UNIT(*) )";
const MILLIMETRE: &str = "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )";
const BARE_METRE: &str = "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )";

fn mm_per_unit(text: &str) -> catallax_numbers::Rational {
    read_step(text.as_bytes()).unwrap().mm_per_unit
}

#[test]
fn reads_a_conversion_unit_named_metre_over_millimetres() {
    let text = box_with_length_unit(&[
        CONVERTED,
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1000.),#@2)",
        MILLIMETRE,
    ]);
    assert_eq!(mm_per_unit(&text), r("1000"));
}

#[test]
fn reads_a_conversion_unit_over_a_bare_metre() {
    let text = box_with_length_unit(&[
        CONVERTED,
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.),#@2)",
        BARE_METRE,
    ]);
    assert_eq!(mm_per_unit(&text), r("1000"));
    let inch = box_with_length_unit(&[
        "( CONVERSION_BASED_UNIT('inch',#@1) LENGTH_UNIT() NAMED_UNIT(*) )",
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#@2)",
        BARE_METRE,
    ]);
    assert_eq!(mm_per_unit(&inch), r("25.4"));
}

#[test]
fn reads_a_chain_of_conversion_units_exactly() {
    let text = box_with_length_unit(&[
        "( CONVERSION_BASED_UNIT('FOOT',#@1) LENGTH_UNIT() NAMED_UNIT(*) )",
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(12.),#@2)",
        "( CONVERSION_BASED_UNIT('INCH',#@3) LENGTH_UNIT() NAMED_UNIT(*) )",
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#@4)",
        MILLIMETRE,
    ]);
    assert_eq!(mm_per_unit(&text), r("304.8"));
}

#[test]
fn reads_every_si_prefix_on_a_metre() {
    let cases = [
        ("EXA", "1000000000000000000000"),
        ("PETA", "1000000000000000000"),
        ("TERA", "1000000000000000"),
        ("GIGA", "1000000000000"),
        ("MEGA", "1000000000"),
        ("KILO", "1000000"),
        ("HECTO", "100000"),
        ("DECA", "10000"),
        ("DECI", "100"),
        ("CENTI", "10"),
        ("MILLI", "1"),
        ("MICRO", "0.001"),
        ("NANO", "0.000001"),
        ("PICO", "0.000000001"),
        ("FEMTO", "0.000000000001"),
        ("ATTO", "0.000000000000001"),
    ];
    for (prefix, mm) in cases {
        let text = box_step(BoxFixture::new(["1", "1", "1"])).replace(
            "SI_UNIT(.MILLI.,.METRE.)",
            &format!("SI_UNIT(.{prefix}.,.METRE.)"),
        );
        assert_eq!(mm_per_unit(&text), r(mm), "{prefix}");
    }
}

#[test]
fn a_conversion_prefix_scales_its_value() {
    let text = box_with_length_unit(&[
        CONVERTED,
        "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(5.),#@2)",
        "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.KILO.,.METRE.) )",
    ]);
    assert_eq!(mm_per_unit(&text), r("5000000"));
}

fn assert_unit_refusal(text: &str, expected: &str) {
    match read_step(text.as_bytes()).unwrap_err() {
        ImportRefusal::UnknownUnit { found } => assert_eq!(found, expected),
        other => panic!("expected an unknown unit, got {other:?}"),
    }
}

#[test]
fn refuses_conflicting_length_units() {
    let two = box_with_length_units(&[MILLIMETRE, BARE_METRE], &[0, 1]);
    assert_unit_refusal(&two, "conflicting length units");
    // Two assigned units that agree are one unit.
    let agree = box_with_length_units(
        &[
            MILLIMETRE,
            "( CONVERSION_BASED_UNIT('METRE',#@2) LENGTH_UNIT() NAMED_UNIT(*) )",
            "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.),#@3)",
            MILLIMETRE,
        ],
        &[0, 1],
    );
    assert_eq!(mm_per_unit(&agree), r("1"));
}

#[test]
fn refuses_a_conversion_unit_without_a_usable_length() {
    let measure = "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.),#@2)";
    // The conversion names a measure that is absent.
    assert_unit_refusal(&box_with_length_unit(&[CONVERTED, MILLIMETRE]), "METRE");
    // Its base unit is not a length.
    let seconds = "( NAMED_UNIT(*) SI_UNIT($,.SECOND.) TIME_UNIT() )";
    assert_unit_refusal(
        &box_with_length_unit(&[CONVERTED, measure, seconds]),
        "SECOND",
    );
    // A zero or negative factor is no unit.
    for value in ["0.", "-1."] {
        let m = measure.replace("1.", value);
        assert_unit_refusal(&box_with_length_unit(&[CONVERTED, &m, MILLIMETRE]), "METRE");
    }
    // No unit is declared at all.
    let none = box_step(BoxFixture::new(["1", "1", "1"])).replace("LENGTH_UNIT()", "AREA_UNIT()");
    assert_unit_refusal(&none, "");
}

/// `n` conversion units of factor one, each over the next, ending on a millimetre.
fn conversion_chain(n: usize) -> String {
    let mut bodies = Vec::new();
    for i in 0..n {
        bodies.push(format!(
            "( CONVERSION_BASED_UNIT('U{i}',#@{}) LENGTH_UNIT() NAMED_UNIT(*) )",
            2 * i + 1
        ));
        bodies.push(format!(
            "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.),#@{})",
            2 * i + 2
        ));
    }
    bodies.push(MILLIMETRE.to_string());
    let refs: Vec<&str> = bodies.iter().map(String::as_str).collect();
    box_with_length_unit(&refs)
}

#[test]
fn refuses_conversion_cycles_and_chains_beyond_the_depth_limit() {
    // The unit converts through a measure that is expressed in the unit itself.
    let cycle = "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.),#@0)";
    assert_unit_refusal(&box_with_length_unit(&[CONVERTED, cycle]), "METRE");
    assert_eq!(mm_per_unit(&conversion_chain(8)), r("1"));
    assert_unit_refusal(&conversion_chain(9), "U8");
}
