//! Exact STEP spline source recognition, units and resource boundaries.

use catallax_import::{
    Curve, ImportRefusal, Surface, build_file, canonical_bytes, read_step, read_step_with_decider,
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_topology::Model;

const SOURCE: &[u8] = include_bytes!("fixtures/spline-roof.step");
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

#[test]
fn reads_weighted_step_geometry_without_rounding_or_certifying_a_solid() {
    let file = read_step(SOURCE).unwrap();
    let roof = file.solids[0]
        .faces
        .iter()
        .find_map(|f| match &f.surface {
            Surface::BSpline(s) => Some(s),
            _ => None,
        })
        .unwrap();
    assert_eq!(roof.geometry.control_count(), [3, 3]);
    assert_eq!(roof.geometry.controls()[4], [q("0.5"), q("0.5"), q("1.25")]);
    assert_eq!(roof.geometry.weights()[4], q("2"));
    assert_eq!(
        roof.geometry.domain(),
        [[&q("0.1"), &q("0.9")], [&q("0.2"), &q("0.8")]]
    );
    let mut dec = Decider::default();
    let surface = roof.geometry.to_surface(&mut dec).unwrap();
    let t = q("0.5");
    let jet = surface.jet_box([&t, &t], [&t, &t], 128, &dec).unwrap();
    for (actual, expected) in jet.value.iter().zip([q("0.5"), q("0.5"), q("1.1")]) {
        assert!(actual.lo() <= expected && actual.hi() >= expected);
    }
    let curve = file.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|lp| &lp.edges)
        .find_map(|e| match &e.curve {
            Curve::BSpline(c) => Some(c),
            _ => None,
        })
        .unwrap();
    assert_eq!(curve.geometry.weights(), &[q("1"), q("2"), q("1")]);
    let evaluator = curve.geometry.to_curve(&mut dec).unwrap();
    let value = evaluator.evaluate_box([&t, &t], 128, &dec).unwrap();
    for (actual, expected) in value.iter().zip([q("0.5"), q("0"), q("1")]) {
        assert!(actual.lo() <= expected && actual.hi() >= expected);
    }
    let mut model = Model::new();
    assert!(matches!(
        build_file(&mut model, &file, None),
        Err(ImportRefusal::Unsupported { .. })
    ));
}

#[test]
fn scales_only_euclidean_lengths_and_retains_the_decimal_definition() {
    let metres = String::from_utf8(SOURCE.to_vec())
        .unwrap()
        .replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)");
    let file = read_step(metres.as_bytes()).unwrap();
    let roof = file.solids[0]
        .faces
        .iter()
        .find_map(|f| match &f.surface {
            Surface::BSpline(s) => Some(s),
            _ => None,
        })
        .unwrap();
    assert_eq!(roof.geometry.controls()[4], [q("500"), q("500"), q("1250")]);
    assert_eq!(roof.geometry.weights()[4], q("2"));
    assert_eq!(roof.geometry.knots()[0][0], q("0.1"));
    let mut changed = file.clone();
    let Surface::BSpline(roof) = &mut changed.solids[0].faces.last_mut().unwrap().surface else {
        panic!("roof");
    };
    roof.closed[0] = Some(true);
    assert_ne!(canonical_bytes(&file), canonical_bytes(&changed));
}

#[test]
fn clamps_uniform_and_simple_knotted_step_encodings_exactly() {
    let original = std::str::from_utf8(SOURCE).unwrap();
    // Use the actual curve's control references, while changing only its STEP
    // subtype and knot storage. Uniform active-domain clamping is exact.
    let line = original
        .lines()
        .find(|s| s.contains("B_SPLINE_CURVE(2,"))
        .unwrap();
    let start = line.find("B_SPLINE_CURVE(2,").unwrap();
    let rest = &line[start + "B_SPLINE_CURVE(2,".len()..];
    let controls = &rest[..rest.find(')').unwrap() + 1];
    let replacement = format!("UNIFORM_CURVE('uniform',2,{controls},.UNSPECIFIED.,.F.,.F.)");
    let id = line.split('=').next().unwrap();
    let simple = original.replace(line, &format!("{id}={replacement};"));
    let file = read_step(simple.as_bytes()).unwrap();
    let c = file.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|lp| &lp.edges)
        .find_map(|e| match &e.curve {
            Curve::BSpline(c) => Some(c),
            _ => None,
        })
        .unwrap();
    assert_eq!(c.geometry.domain(), [&q("2"), &q("3")]);
    assert_eq!(c.geometry.controls()[0], [q("0.25"), q("0"), q("1")]);
    assert_eq!(c.geometry.controls()[2], [q("0.75"), q("0"), q("1")]);
}

#[test]
fn refuses_malformed_nets_and_exhausted_work_without_partial_source_data() {
    let original = std::str::from_utf8(SOURCE).unwrap();
    for (from, to) in [
        (
            "((1.,1.,1.),(1.,2.,1.),(1.,1.,1.))",
            "((1.,1.,1.),(1.,0.,1.),(1.,1.,1.))",
        ),
        ("(3,3),(3,3),(0.1,0.9)", "(3,2),(3,3),(0.1,0.9)"),
        ("(3,3),(3,3),(0.1,0.9)", "(3,3),(3,3),(0.9,0.1)"),
    ] {
        assert!(matches!(
            read_step(original.replace(from, to).as_bytes()),
            Err(ImportRefusal::Malformed { .. })
        ));
    }
    let (result, work) =
        with_operation_budget(1, || read_step_with_decider(SOURCE, &Decider::default()));
    assert!(matches!(result, Err(ImportRefusal::ResourceLimit)));
    assert!(work.exhausted);
}
