//! Actual skew collar through the numerical recipe and parametric replay path.

use super::*;
use catallax_features::{NumericalHoleBuilt, evaluate_numerical_hole};
use catallax_recipe::Recipe;
use catallax_topology::numerical::NumericalBody;

fn request<T>(label: &str, run: impl FnOnce() -> Result<T, String>) -> (T, u64) {
    let (result, work) = with_operation_budget(20_000_000, run);
    assert!(!work.exhausted, "{label}: {work:?}");
    (
        result.unwrap_or_else(|error| panic!("{label}: {error}")),
        work.used,
    )
}
fn names(body: &NumericalBody) -> Vec<String> {
    body.definition()
        .faces
        .iter()
        .map(|face| face.key.clone())
        .collect()
}
fn recipe(diameter: &Rational) -> Recipe {
    let mut value = serde_json::to_value(super::skew::source()).unwrap();
    value["steps"].as_array_mut().unwrap().push(json!({"id":"hole","op":"hole","target":"collar","diameter":diameter.to_string(),"depth":"through","tip":"flat","at":[{"id":"radial-0","point":["21/2","12","0"],"direction":["1","-1","1"]}]}));
    value["result"] = json!("hole");
    Recipe::parse(&value.to_string()).unwrap()
}
#[test]
fn skew_collar_builds_edits_saves_reopens_and_exports_with_native_names() {
    let mut previous: Option<NumericalHoleBuilt> = None;
    let mut records = Vec::new();
    for (label, diameter) in [("initial", q("2")), ("edited", q("3/2"))] {
        let mut dec = Decider::default();
        let recipe = recipe(&diameter);
        let options = NumericalOptions {
            tolerance_mm: q("1e-7"),
            vertex_tolerance_mm: None,
        };
        let (next, build_work) = request(label, || {
            let next = if let Some(previous) = &previous {
                previous.edit(&recipe, options.clone(), &mut dec)
            } else {
                evaluate_numerical_hole(&recipe, options.clone(), &mut dec)
            }
            .map_err(|error| error.to_string())?
            .expect("actual two-bevel collar input family");
            let hole = next.hole();
            let body = hole.body();
            assert_eq!(body.kind(), catallax_topology::BodyKind::Solid);
            assert!(
                body.exact_topology_faults()
                    .iter()
                    .all(|fault| fault.code == "topology.missing-incidence")
            );
            assert_eq!(body.options().tolerance_mm, q("1e-7"));
            assert_eq!(body.history().len(), usize::from(label == "edited"));
            assert_eq!(
                hole.source_faces().len() + hole.tool_faces().len(),
                body.definition().faces.len()
            );
            let mut source_counts = std::collections::BTreeMap::new();
            for source in hole.source_faces().values() {
                *source_counts.entry(*source).or_insert(0) += 1;
            }
            assert!(
                source_counts.values().any(|count| *count == 5),
                "all five regions of the intersected cone must survive"
            );
            assert_eq!(next.face_origins().len(), body.definition().faces.len());
            assert!(
                next.face_origins()
                    .values()
                    .all(|origins| !origins.is_empty())
            );
            if let Some(previous) = &previous {
                assert_eq!(names(body), names(previous.hole().body()));
                assert_eq!(next.face_origins(), previous.face_origins());
            }
            Ok(next)
        });
        let (live, export_work) = request("live STEP", || {
            catallax_step::numerical::build(
                next.hole().body(),
                &catallax_step::numerical::NumericalStepOptions::default(),
                &mut dec,
            )
            .map_err(|error| format!("{error:?}"))
        });
        // The native part file stores authored features. Reopen rebuilds from
        // those definitions rather than serializing a numerical B-rep cache.
        let (saved, save_work) = request("save", || {
            serde_json::to_vec(&json!({"format":"catallax-numerical-step/1","recipe":next.recipe(),"tolerance_mm":options.tolerance_mm.to_string()})).map_err(|error|error.to_string())
        });
        let (restored, reopen_work) = request("reopen", || {
            let value: serde_json::Value =
                serde_json::from_slice(&saved).map_err(|error| error.to_string())?;
            let recipe =
                Recipe::parse(&value["recipe"].to_string()).map_err(|error| error.to_string())?;
            let options = NumericalOptions {
                tolerance_mm: Rational::parse(value["tolerance_mm"].as_str().unwrap()).unwrap(),
                vertex_tolerance_mm: None,
            };
            evaluate_numerical_hole(&recipe, options, &mut dec)
                .map_err(|error| error.to_string())?
                .ok_or_else(|| "saved numerical Hole lost its input family".into())
        });
        assert_eq!(names(restored.hole().body()), names(next.hole().body()));
        assert_eq!(restored.face_origins(), next.face_origins());
        assert_eq!(
            restored.hole().body().options(),
            next.hole().body().options()
        );
        let (reopened, cold_export_work) = request("reopened STEP", || {
            catallax_step::numerical::build(
                restored.hole().body(),
                &catallax_step::numerical::NumericalStepOptions::default(),
                &mut dec,
            )
            .map_err(|error| format!("{error:?}"))
        });
        assert_eq!(live.bytes(), reopened.bytes());
        if let Some(path) = std::env::var_os("CATALLAX_SKEW_WALL_PROOF") {
            let path = std::path::PathBuf::from(path);
            std::fs::create_dir_all(&path).unwrap();
            std::fs::write(path.join(format!("collar-{label}.step")), live.bytes()).unwrap();
            std::fs::write(
                path.join(format!("collar-{label}-reopened.step")),
                reopened.bytes(),
            )
            .unwrap();
            std::fs::write(path.join(format!("collar-{label}.recipe.json")), &saved).unwrap();
        }
        records.push(json!({"state":label,"diameter_mm":diameter.to_string(),"build_work":build_work,"build_stages":next.work_stages(),"live_export_work":export_work,"save_work":save_work,"reopen_work":reopen_work,"cold_export_work":cold_export_work,"faces":names(next.hole().body()),"origins":next.face_origins(),"native_trims":next.hole().body().native_trims().len(),"status":"numerical","tolerance_mm":"1e-7","target_work_per_request":20_000_000}));
        previous = Some(next);
    }
    if let Some(path) = std::env::var_os("CATALLAX_SKEW_WALL_PROOF") {
        std::fs::write(
            std::path::PathBuf::from(path).join("summary.json"),
            serde_json::to_vec_pretty(&records).unwrap(),
        )
        .unwrap();
    }
    println!("{}", serde_json::to_string(&records).unwrap());
}
