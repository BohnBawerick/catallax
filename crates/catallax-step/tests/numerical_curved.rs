//! The actual loft crossing retains its native curves and physical patch hole.

use catallax_numbers::{
    Decider, NumericalCheck, Rational, current_operation_work, with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::numerical::{self, NumericalOptions};
use serde_json::{Value, json};

fn cut_body(
    name: &str,
    size: [&str; 3],
    at: [&str; 3],
    dec: &mut Decider,
) -> numerical::NumericalBody {
    let mut input: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/loft-spline-prism.json"
    ))
    .unwrap();
    let steps = input["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"box","op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":at}}));
    steps.push(json!({"id":"cut","op":"boolean","kind":"subtract","target":"loft","tool":"box"}));
    input["result"] = json!("cut");
    let recipe = Recipe::parse(&input.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, dec).unwrap();
    let definition =
        numerical::definition_from_model(&built.model, built.result, name.into(), dec).unwrap();
    numerical::build(
        definition,
        NumericalOptions {
            tolerance_mm: Rational::parse("1/10000000").unwrap(),
            vertex_tolerance_mm: None,
        },
        dec,
    )
    .unwrap()
}

#[test]
fn the_loft_crossing_reaches_approximate_step_with_its_hole() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let body = cut_body(
            "loft crossing",
            ["10", "10", "1"],
            ["-3", "1999999/1000000", "1"],
            &mut dec,
        );
        let before = current_operation_work().used;
        let file = step::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
        let text = std::str::from_utf8(file.bytes()).unwrap();
        assert!(text.contains("numerical approximate"));
        assert!(text.contains("complete_body_bound=unestablished"));
        assert_eq!(text.matches("FACE_OUTER_BOUND(").count(), 7);
        assert_eq!(text.matches("=FACE_BOUND(").count(), 1);
        assert_eq!(text.matches("=EDGE_CURVE(").count(), 12);
        assert!(text.contains("RATIONAL_B_SPLINE_SURFACE("));
        assert!(text.contains("RATIONAL_B_SPLINE_CURVE("));
        assert!(text.contains("=PCURVE("));
        assert_eq!(file.report().counts, [7, 12, 8]);
        assert_eq!(file.report().source_definition.faces[0].loops.len(), 2);
        assert_eq!(
            file.report().source_evidence.tolerance_mm(),
            body.evidence().tolerance_mm()
        );
        assert!(file.report().spline_rounding.is_some());
        assert!(file.report().emitted_evidence.checks().iter().any(|check| matches!(check,NumericalCheck::Passed {name,scope} if name=="trim.hole-contained-disjoint" && scope==&body.definition().faces[0].key)));
        assert!(file.report().emitted_evidence.checks().iter().any(|check| matches!(check,NumericalCheck::Passed {name,..} if name=="trim.sampled-simple-oriented-parameter-loop")));
        let output = std::env::temp_dir().join(format!(
            "catallax-numerical-curved-{}.step",
            std::process::id()
        ));
        std::fs::write(&output, file.bytes()).unwrap();
        let reread = std::fs::read(&output).unwrap();
        std::fs::remove_file(&output).unwrap();
        file.check_written(&reread, &mut dec).unwrap();
        eprintln!(
            "loft admission={before} approximate STEP={} work units",
            current_operation_work().used - before
        );
        let coarse = NumericalStepOptions {
            significant_digits: 1,
            max_scalar_rounding_error_mm: Rational::one(),
            ..Default::default()
        };
        assert!(step::build(&body, &coarse, &mut dec).is_err());
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[test]
fn recorded_m10_008_keeps_its_hole_and_rejects_broken_loop_controls() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let body = cut_body("m7-cut", ["2", "4", "1"], ["1", "1", "1"], &mut dec);
        let hole = body
            .definition()
            .faces
            .iter()
            .position(|face| face.loops.len() == 2)
            .unwrap();
        let file = step::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
        let text = std::str::from_utf8(file.bytes()).unwrap();
        assert!(text.contains("RATIONAL_B_SPLINE_SURFACE("));
        assert!(text.contains("=PCURVE("));
        assert_eq!(text.matches("=FACE_BOUND(").count(), 1);
        assert_eq!(file.report().source_definition.faces[hole].loops.len(), 2);
        assert_eq!(
            file.report().counts,
            [
                body.counts().faces,
                body.counts().edges,
                body.counts().vertices
            ]
        );
        let mut reversed = body.definition().clone();
        let inner = &mut reversed.faces[hole].loops[1];
        inner.reverse();
        for trim in inner {
            trim.uv.reverse();
            trim.sense = if trim.sense == catallax_topology::Sense::Same {
                catallax_topology::Sense::Reversed
            } else {
                catallax_topology::Sense::Same
            };
        }
        assert!(numerical::build(reversed, body.options().clone(), &mut dec).is_err());
        let mut missing = body.definition().clone();
        missing.faces[hole].loops.pop();
        assert!(numerical::build(missing, body.options().clone(), &mut dec).is_err());
        file.check_written(file.bytes(), &mut dec).unwrap();
        eprintln!(
            "M10-008 numerical STEP counts={:?} total_work={}",
            file.report().counts,
            current_operation_work().used
        );
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[path = "common/n03_tube.rs"]
mod tube;

#[test]
fn repeated_endpoint_carrier_keeps_the_admitted_finite_traversal() {
    use catallax_numbers::Real;
    use catallax_shapes::{BSplineCurve, Curve, Point3};
    let mut dec = Decider::default();
    let mut definition = tube::tube();
    let edge = &definition.edges[0];
    let point = |id: usize| {
        let p = &definition.vertices[id].position;
        Point3::new(
            p[0].clone().into(),
            p[1].clone().into(),
            p[2].clone().into(),
        )
    };
    let a = point(edge.vertices[0]);
    let b = point(edge.vertices[1]);
    let c = Point3::new(9.into(), 9.into(), 9.into());
    definition.edges[0].curve = Curve::BSpline(
        BSplineCurve::new(
            1,
            [0, 0, 1, 2, 3, 4, 4].map(Rational::from).to_vec(),
            vec![a.clone(), b.clone(), c, a, b],
            vec![Real::one(); 5],
            &mut dec,
        )
        .unwrap(),
    );
    let body = numerical::build(
        definition,
        NumericalOptions {
            tolerance_mm: tube::q("1/10000000"),
            vertex_tolerance_mm: None,
        },
        &mut dec,
    )
    .unwrap();
    let file = step::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
    let text = std::str::from_utf8(file.bytes()).unwrap();
    assert!(
        text.contains("=TRIMMED_CURVE("),
        "finite admitted interval missing from STEP"
    );
    assert!(text.contains("PARAMETER_VALUE(0.)") && text.contains("PARAMETER_VALUE(1.)"));
    if let Ok(directory) = std::env::var("CATALLAX_N03_CAPTURE") {
        std::fs::write(
            std::path::Path::new(&directory).join("repeated-endpoint.step"),
            file.bytes(),
        )
        .unwrap();
    }
    file.check_written(file.bytes(), &mut dec).unwrap();
}

#[path = "common/copy_tags.rs"]
mod copy_tags;

/// The writer of spline faces and retained u/v curves.
#[test]
fn the_written_spline_body_measures_every_use_whatever_its_tag() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let body = cut_body(
            "loft crossing",
            ["10", "10", "1"],
            ["-3", "1999999/1000000", "1"],
            &mut Decider::default(),
        );
        copy_tags::assert_the_written_body_is_measured(&body, "1/10000000");
    });
    assert!(!work.exhausted, "{} work units", work.used);
}
