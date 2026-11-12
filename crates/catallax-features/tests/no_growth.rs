//! M12 no-growth evidence for later features on P1-P5 and retained-law offsets.

use catallax_features::{Built, evaluate, evaluate_certified};
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{SplineSupport, Surface};
use catallax_topology::{BodyId, FIT_MARK, fitted_face_bound};
use serde_json::{Value, json};
use std::num::NonZeroUsize;

const ALLOWANCE: u64 = 20_000_000;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn later_features(source: &str) -> (Recipe, Vec<StepId>) {
    let mut value: Value = serde_json::from_str(source).unwrap();
    let original = value["result"].as_str().unwrap().to_owned();
    value["steps"].as_array_mut().unwrap().extend([
        json!({"id":"far_tool","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["1000","1000","1000"]}}),
        json!({"id":"cut_once","op":"boolean","kind":"subtract","target":original,"tool":"far_tool"}),
        json!({"id":"reflected","op":"mirror","body":"cut_once","plane":{"point":["7","11","13"],"normal":["1","2","3"]},"keep_original":false}),
        json!({"id":"returned","op":"mirror","body":"reflected","plane":{"point":["7","11","13"],"normal":["1","2","3"]},"keep_original":false}),
        json!({"id":"cut_twice","op":"boolean","kind":"subtract","target":"returned","tool":"far_tool"}),
    ]);
    value["result"] = json!("cut_twice");
    let stages = [
        original.as_str(),
        "cut_once",
        "reflected",
        "returned",
        "cut_twice",
    ]
    .map(|id| StepId::new(id).unwrap())
    .to_vec();
    (Recipe::parse(&value.to_string()).unwrap(), stages)
}

fn measured_faces(built: &Built, body: BodyId) -> Vec<Rational> {
    let model = &built.model;
    assert!(model.check(body).is_empty());
    let mut bounds = Vec::new();
    for face in model.faces_of_body(body) {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RollingEnvelope(fit)) => {
                let measured = fit.remeasure(&mut Decider::default()).unwrap();
                assert_eq!(measured, fit.bound().upper);
                assert!(measured <= q("0.0001"));
                bounds.push(measured);
            }
            Surface::BSpline(_) => panic!("unmeasured fitted family in target part"),
            _ => {}
        }
    }
    bounds.sort();
    let worst = bounds.last().cloned();
    assert_eq!(fitted_face_bound(model, body), worst);
    let marks: Vec<_> = model
        .body(body)
        .marks
        .iter()
        .filter(|mark| mark.code == FIT_MARK)
        .map(|mark| mark.bound.clone())
        .collect();
    assert_eq!(marks, worst.into_iter().collect::<Vec<_>>());
    bounds
}

fn signature(built: &Built, body: BodyId) -> Vec<String> {
    let model = &built.model;
    let mut out = vec![format!("{:?}", model.counts(body))];
    let mut dec = Decider::default();
    for face in model.faces_of_body(body) {
        let f = model.face(face);
        out.push(format!("{:?} {:?}", model.face_origins(face), f.sense));
        if let Surface::BSpline(support) = model.geom().surface(f.surface) {
            let carrier = support.carrier();
            out.push(format!(
                "{:?} {:?} {:?}",
                carrier.degrees(),
                carrier.knots(),
                carrier.weights()
            ));
            for p in carrier.controls() {
                out.push(format!(
                    "{:?}",
                    [&p.x, &p.y, &p.z].map(|x| dec.enclose_at(x, 128).unwrap())
                ));
            }
        }
    }
    for edge in model.edges_of_body(body) {
        let e = model.edge(edge);
        out.push(format!("{:?} {:?}", model.edge_name(body, edge), e.sense));
        for v in [e.start, e.end].into_iter().flatten() {
            let p = model.geom().position(model.vertex(v).point);
            out.push(format!(
                "{:?}",
                [&p.x, &p.y, &p.z].map(|x| dec.enclose_at(x, 128).unwrap())
            ));
        }
    }
    out
}

fn offset_recipe(chained: bool) -> Recipe {
    let mut value = json!({
        "format":"catallax-recipe/3", "name":"M12 retained-law offset chain",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["8","8","4"],"pose":{"rotation":["3","0","1","0"],"translation":["-2","-4","4"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"}
        ], "result":"offset_0"
    });
    let mut parent = "source".to_owned();
    for i in 0..if chained { 4 } else { 1 } {
        let id = format!("offset_{i}");
        value["steps"].as_array_mut().unwrap().push(json!({
            "id":id,"op":"offset","body":parent,
            "distance":if chained { "-1/16" } else { "-1/4" }
        }));
        parent = id;
    }
    value["result"] = json!(parent);
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn recorded_error_does_not_accumulate_at_one_two_or_four_workers() {
    let parts = [
        ("P1", include_str!("../../../parts/P1-bracket.json"), 0),
        ("P2", include_str!("../../../parts/P2-collar.json"), 0),
        (
            "P3",
            include_str!("../../../parts/P3-filleted-boss.json"),
            0,
        ),
        ("P4", include_str!("../../../parts/P4-shelled-tray.json"), 0),
        (
            "P5",
            include_str!("../../../parts/P5-flanged-housing.json"),
            4,
        ),
    ];
    let mut baseline = None;
    for workers in [1, 2, 4] {
        set_threads(NonZeroUsize::new(workers));
        let mut run = Vec::new();
        for (label, source, fitted_count) in parts {
            let base = Recipe::parse(source).unwrap();
            let (base_result, base_work) =
                with_operation_budget(ALLOWANCE, || evaluate(&base, &mut Decider::default()));
            let (recipe, stages) = later_features(source);
            let (result, work) =
                with_operation_budget(ALLOWANCE, || evaluate(&recipe, &mut Decider::default()));
            if let Err(base_error) = base_result {
                let error = result
                    .err()
                    .expect("the same base recipe must still refuse");
                assert_eq!(error.to_string(), base_error.to_string());
                assert_eq!(work.used, base_work.used);
                assert_eq!(work.exhausted, base_work.exhausted);
                eprintln!("{label}: workers={workers}, no-growth unproved: {error}; {work:?}");
                run.push((vec![error.to_string()], vec![], work.used));
                continue;
            }
            let built = result.unwrap_or_else(|err| panic!("{label}: {err:?}; {work:?}"));
            assert!(!work.exhausted);
            let original = built.bodies[&stages[0]];
            let before = measured_faces(&built, original);
            assert_eq!(before.len(), fitted_count, "{label}");
            for stage in stages {
                let body = built.bodies[&stage];
                let measured = measured_faces(&built, body);
                assert_eq!(measured, before, "{label} {stage}");
                assert_eq!(built.model.counts(body), built.model.counts(original));
                run.push((signature(&built, body), measured, work.used));
            }
            eprintln!(
                "{label}: workers={workers}, stages=5, fitted={fitted_count}, bound_mm={before:?}, work={}",
                work.used
            );
        }
        let mut final_bounds = Vec::new();
        for chained in [false, true] {
            let recipe = offset_recipe(chained);
            let (result, work) = with_operation_budget(ALLOWANCE, || {
                evaluate_certified(&recipe, &mut Decider::default())
            });
            let certified =
                result.unwrap_or_else(|err| panic!("offset chained={chained}: {err:?}; {work:?}"));
            assert!(!work.exhausted);
            certified.verify(&mut Decider::default()).unwrap();
            let mut bounds = Vec::new();
            for i in 0..if chained { 4 } else { 1 } {
                let certificate = certified
                    .certificate(&StepId::new(format!("offset_{i}")).unwrap())
                    .unwrap();
                let target = certificate.target();
                let measured = catallax_field::measurement::measure_elliptic_offset_sweep(
                    certificate.spline().unwrap().carrier(),
                    target.definition(),
                    target.axial(),
                    &catallax_field::measurement::MeasurementOptions {
                        width: q("0.0001"),
                        max_cells: 1,
                        path_segments: 1,
                    },
                    &mut Decider::default(),
                )
                .unwrap();
                assert!(measured.converged);
                for region in certificate.record().regions.values() {
                    assert_eq!(
                        region.exception.as_ref().unwrap().bound.upper,
                        measured.upper
                    );
                }
                bounds.push(measured.upper);
            }
            assert_eq!(
                certified.result_certificate().target().distance_mm(),
                &q("1/4")
            );
            final_bounds.push(bounds.last().unwrap().clone());
            eprintln!(
                "offset: workers={workers}, chained={chained}, bound_mm={bounds:?}, work={}",
                work.used
            );
            run.push((
                signature(certified.built(), certified.built().result),
                bounds,
                work.used,
            ));
        }
        assert_eq!(
            final_bounds[0], final_bounds[1],
            "refit the same final law without inherited error"
        );
        if let Some(before) = &baseline {
            assert_eq!(&run, before, "workers={workers}");
        } else {
            baseline = Some(run);
        }
    }
    set_threads(None);
}
