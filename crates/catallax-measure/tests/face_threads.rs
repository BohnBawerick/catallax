use std::num::NonZeroUsize;

use catallax_measure::{Accuracy, face_areas, mass_properties};
use catallax_numbers::{
    Decider, OperationWork, Rational, set_run_ahead_wait, set_threads, with_operation_budget,
};
use catallax_recipe::Recipe;

const LOFT: &str = r#"{
  "format": "catallax-recipe/5", "name": "Asymmetric ruled loft", "result": "loft",
  "steps": [
    {"id":"a", "op":"sketch", "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
     "loops":[[{"kind":"line","id":"a1","to":["1","-1"]},
               {"kind":"line","id":"a2","to":["1","1"]},
               {"kind":"line","id":"a3","to":["-1","1"]},
               {"kind":"line","id":"a4","to":["-1","-1"]}]]},
    {"id":"b", "op":"sketch", "plane":{"rotation":["1","0","0","0"],"translation":["0","0","2"]},
     "loops":[[{"kind":"line","id":"b1","to":["1","0"]},
               {"kind":"line","id":"b2","to":["0","2"]},
               {"kind":"line","id":"b3","to":["-2","0"]},
               {"kind":"line","id":"b4","to":["0","-1"]}]]},
    {"id":"loft", "op":"loft", "shape":"ruled",
     "sections":[{"profile":"a","start":["a1"]},{"profile":"b","start":["b1"]}]}
  ]
}"#;

fn work_fields(work: OperationWork) -> (u64, bool) {
    (work.used, work.exhausted)
}

#[test]
fn loft_face_areas_and_mass_keep_values_and_work_in_request_order() {
    set_run_ahead_wait(true);
    let run = |threads| {
        set_threads(NonZeroUsize::new(threads));
        // Start with a fresh contribution cache on each run; retain it for the
        // second mass call, as an application request would.
        std::thread::spawn(|| {
            let recipe: Recipe = serde_json::from_str(LOFT).unwrap();
            let app = catallax_app::build_incremental_without_display(
                &recipe,
                None,
                &mut Decider::default(),
            )
            .unwrap();
            let body = app.checkpoint.built();
            let faces = body.model.faces_of_body(body.result);
            let accuracy = Accuracy {
                relative: Rational::zero(),
                absolute: Rational::parse("1e-2").unwrap(),
            };
            let areas = |limit| {
                with_operation_budget(limit, || {
                    face_areas(&body.model, &faces, &accuracy, &mut Decider::default())
                })
            };
            let mass = |limit| {
                with_operation_budget(limit, || {
                    mass_properties(&body.model, body.result, &accuracy, &mut Decider::default())
                })
            };
            let before_areas = catallax_numbers::run_ahead_adopted();
            let (full, face_work) = areas(20_000_000);
            assert!(full.iter().all(Result::is_ok), "{full:?}");
            if catallax_numbers::threads() > 1 {
                // There are only two planar caps. At least one worker must
                // have supplied a spline integral, not just a cap.
                assert!(catallax_numbers::run_ahead_adopted() - before_areas > 2);
            }
            let (short, short_work) = areas(face_work.used / 2);
            assert!(short.iter().any(Result::is_err));
            let before_mass = catallax_numbers::run_ahead_adopted();
            let (cold, cold_work) = mass(20_000_000);
            assert!(cold.is_ok(), "{cold:?}");
            if catallax_numbers::threads() > 1 {
                assert!(catallax_numbers::run_ahead_adopted() > before_mass);
            }
            let (warm, warm_work) = mass(20_000_000);
            assert!(warm.is_ok(), "{warm:?}");
            let (refusal, refusal_work) = mass(cold_work.used / 2);
            assert!(refusal.is_err(), "{refusal:?}");
            (
                format!("{full:?}\n{short:?}\n{cold:?}\n{warm:?}\n{refusal:?}"),
                [face_work, short_work, cold_work, warm_work, refusal_work].map(work_fields),
            )
        })
        .join()
        .unwrap()
    };
    let serial = run(1);
    for threads in [2, 4, 24] {
        assert_eq!(run(threads), serial, "{threads} threads");
    }
}
