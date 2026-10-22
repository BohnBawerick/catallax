//! Complete-check measurements for fixed representative recipes.
//! Run `cargo run --release -p catallax-compare --example measure -- [runs] [case]`.
use catallax_compare::{Completion, Controller, Options, Reason};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::time::Instant;
fn pose(x: &str, y: &str, z: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,y,z]})
}
fn recipe(name: &str, steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(
        &json!({"format":catallax_recipe::FORMAT,"name":name,"steps":steps,"result":result})
            .to_string(),
    )
    .unwrap()
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let runs: usize = args.first().map(|s| s.parse()).transpose()?.unwrap_or(3);
    let block =
        || json!({"id":"stock","op":"block","size":["20","10","5"],"pose":pose("0","0","0")});
    let recipes = vec![
        recipe("block", vec![block()]),
        recipe(
            "cylinder",
            vec![
                json!({"id":"stock","op":"cylinder","radius":"2","height":"5","pose":pose("0","0","0")}),
            ],
        ),
        recipe(
            "sphere",
            vec![json!({"id":"stock","op":"sphere","radius":"2","pose":pose("0","0","0")})],
        ),
        recipe(
            "frustum",
            vec![
                json!({"id":"stock","op":"cone","base_radius":"2","top_radius":"1","height":"4","pose":pose("0","0","0")}),
            ],
        ),
        recipe(
            "torus",
            vec![
                json!({"id":"stock","op":"torus","major":"3","minor":"1","pose":pose("0","0","0")}),
            ],
        ),
        recipe(
            "joined_blocks",
            vec![
                block(),
                json!({"id":"tool","op":"block","size":["5","10","5"],"pose":pose("20","0","0")}),
                json!({"id":"result","op":"boolean","kind":"union","target":"stock","tool":"tool"}),
            ],
        ),
        recipe(
            "bored_block",
            vec![
                block(),
                json!({"id":"tool","op":"cylinder","radius":"3/2","height":"7","pose":pose("5","5","-1")}),
                json!({"id":"result","op":"boolean","kind":"subtract","target":"stock","tool":"tool"}),
            ],
        ),
        recipe(
            "intersecting_spheres",
            vec![
                json!({"id":"stock","op":"sphere","radius":"2","pose":pose("0","0","0")}),
                json!({"id":"tool","op":"sphere","radius":"2","pose":pose("2","0","0")}),
                json!({"id":"result","op":"boolean","kind":"intersect","target":"stock","tool":"tool"}),
            ],
        ),
    ];
    for recipe in recipes {
        if args.get(1).is_some_and(|name| name != &recipe.name) {
            continue;
        }
        let controller = Controller::default();
        for run in 1..=runs {
            let revision = controller.begin();
            let start = Instant::now();
            assert_eq!(
                revision.check(&recipe, &Options::default()),
                Completion::Published
            );
            let elapsed = start.elapsed();
            let report = controller.latest().unwrap();
            let disagreements = report
                .findings
                .iter()
                .filter(|f| !matches!(f.reason, Reason::Unresolved(_)))
                .count();
            println!(
                "case={} run={} elapsed_ms={:.3} verdict={:?} forward={} reverse={} disagreements={} unresolved={}",
                recipe.name,
                run,
                elapsed.as_secs_f64() * 1000.,
                report.verdict,
                report.forward_vertices,
                report.reverse_samples,
                disagreements,
                report.findings.len() - disagreements
            );
            for finding in report
                .findings
                .iter()
                .filter(|f| !matches!(f.reason, Reason::Unresolved(_)))
                .take(1)
            {
                println!("first_disagreement={finding:?}");
            }
        }
    }
    Ok(())
}
