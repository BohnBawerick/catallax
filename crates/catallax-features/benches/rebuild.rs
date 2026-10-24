//! Public history-growth benchmark. Timings are observations, not kernel policy.
use catallax_features::{Built, evaluate, evaluate_incremental};
use catallax_numbers::{Decider, FactKind, with_operation_log};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::GeomKey;
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    hint::black_box,
    time::Instant,
};

fn rectangle(id: &str, x: usize, width: usize, height: usize) -> Value {
    json!({"id":id,"op":"sketch",
    "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
    "loops":[[
        {"kind":"line","id":"bottom","to":[(x+width).to_string(),"8"]},
        {"kind":"line","id":"right","to":[(x+width).to_string(),(8+height).to_string()]},
        {"kind":"line","id":"top","to":[x.to_string(),(8+height).to_string()]},
        {"kind":"line","id":"left","to":[x.to_string(),"8"]}
    ]]})
}

fn recipe(distinct: bool) -> Recipe {
    let mut base = rectangle("base_profile", 0, 320, 20);
    for line in base["loops"][0].as_array_mut().unwrap() {
        let y = line["to"][1].as_str().unwrap().parse::<usize>().unwrap() - 8;
        line["to"][1] = json!(y.to_string());
    }
    let mut steps = vec![
        base,
        json!({"id":"base","op":"extrude","profile":"base_profile","distance":"10"}),
    ];
    let mut target = "base".to_string();
    for cut in 1..=30 {
        let profile = format!("profile{cut}");
        let tool = format!("tool{cut}");
        let result = format!("cut{cut}");
        steps.push(rectangle(
            &profile,
            if distinct { 9 + 10 * (cut - 1) } else { 9 },
            2,
            4,
        ));
        steps.push(json!({"id":tool,"op":"extrude","profile":profile,"distance":"12"}));
        steps.push(
            json!({"id":result,"op":"boolean","kind":"subtract","target":target,"tool":tool}),
        );
        target = result;
    }
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"thirty-cuts","steps":steps,"result":target})
            .to_string(),
    )
    .unwrap()
}

fn export(built: &Built) -> Vec<u8> {
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
}

// Previous closure algorithm, over exactly the same retained facts. This isolates
// query indexing from geometry pruning, boolean shortcuts and numeric decisions.
fn graph_holds(
    edges: &BTreeMap<(FactKind, GeomKey), BTreeSet<GeomKey>>,
    kind: FactKind,
    a: GeomKey,
    b: GeomKey,
) -> bool {
    if kind != FactKind::Incident && a == b {
        return true;
    }
    if let Some(neighbors) = edges.get(&(kind, a)) {
        if neighbors.contains(&b) {
            return true;
        }
    } else if kind != FactKind::Incident || !edges.contains_key(&(FactKind::Coincident, a)) {
        return false;
    }
    let mut seen = BTreeSet::from([(a, false)]);
    let mut queue = VecDeque::from([(a, false)]);
    while let Some((key, used)) = queue.pop_front() {
        for edge_kind in [kind, FactKind::Coincident] {
            if edge_kind != kind && kind != FactKind::Incident {
                continue;
            }
            if let Some(neighbors) = edges.get(&(edge_kind, key)) {
                for &next in neighbors {
                    let used = used || edge_kind == FactKind::Incident;
                    if next == b && (kind != FactKind::Incident || used) {
                        return true;
                    }
                    if seen.insert((next, used)) {
                        queue.push_back((next, used));
                    }
                }
            }
        }
    }
    false
}

fn counterfactual(built: &Built) {
    let facts = built.model.geom().facts();
    let mut edges: BTreeMap<_, BTreeSet<_>> = BTreeMap::new();
    for fact in facts.iter() {
        edges.entry((fact.kind, fact.a)).or_default().insert(fact.b);
        if fact.kind != FactKind::Incident {
            edges.entry((fact.kind, fact.b)).or_default().insert(fact.a);
        }
    }
    let points: Vec<_> = built
        .model
        .vertices_of_body(built.result)
        .into_iter()
        .map(|v| GeomKey::Point(built.model.vertex(v).point))
        .collect();
    let surfaces: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|f| GeomKey::Surface(built.model.face(f).surface))
        .collect();
    let queries: Vec<_> = points
        .iter()
        .flat_map(|&a| {
            points
                .iter()
                .map(move |&b| (FactKind::Coincident, a, b))
                .chain(surfaces.iter().map(move |&b| (FactKind::Incident, a, b)))
        })
        .collect();
    for &(kind, a, b) in &queries {
        assert_eq!(facts.holds(kind, a, b), graph_holds(&edges, kind, a, b));
    }
    for indexed in [false, true] {
        let start = Instant::now();
        for _ in 0..20 {
            for &(kind, a, b) in &queries {
                black_box(if indexed {
                    facts.holds(kind, a, b)
                } else {
                    graph_holds(&edges, kind, a, b)
                });
            }
        }
        println!(
            "{}",
            json!({"counterfactual":"same-facts","indexed":indexed,
            "queries":20*queries.len(),"wall_ms":start.elapsed().as_secs_f64()*1000.0})
        );
    }
}

fn main() {
    let flag = |name| std::env::args().any(|arg| arg == name);
    let distinct = flag("--distinct");
    let recipe = recipe(distinct);
    let fork = flag("--fork");
    if flag("--recipe") {
        print!("{}", recipe.to_text());
        return;
    }
    let expected_faces = if distinct { 126 } else { 10 };
    let built = if flag("--incremental") || fork {
        let mut prefix = recipe.clone();
        prefix.steps.truncate(2);
        prefix.result = StepId::new("base").unwrap();
        let mut checkpoint = evaluate_incremental(&prefix, None, &mut Decider::default()).unwrap();
        let mut retained = Vec::new();
        let mut undo_bytes = Vec::new();
        for cut in 1..=30 {
            let mut edited = recipe.clone();
            edited.steps.truncate(2 + 3 * cut);
            edited.result = StepId::new(format!("cut{cut}")).unwrap();
            let mut fork_ms = Vec::new();
            let next = if fork {
                for _ in 0..5 {
                    let start = Instant::now();
                    let snapshot = black_box(checkpoint.clone());
                    fork_ms.push(start.elapsed().as_secs_f64() * 1000.0);
                    drop(snapshot);
                }
                fork_ms.sort_by(f64::total_cmp);
                let next = checkpoint.clone();
                if flag("--verify") {
                    undo_bytes.push(export(checkpoint.built()));
                }
                retained.push(checkpoint);
                next
            } else {
                checkpoint
            };
            let start = Instant::now();
            checkpoint =
                evaluate_incremental(&edited, Some(next), &mut Decider::default()).unwrap();
            let built = checkpoint.built();
            assert!(built.model.check(built.result).is_empty());
            if flag("--display") {
                catallax_mesh::mesh_body(
                    &built.model,
                    built.result,
                    &Default::default(),
                    &mut Decider::default(),
                )
                .unwrap();
            }
            let elapsed = start.elapsed().as_secs_f64() * 1000.0;
            println!(
                "{}",
                json!({"mode":if fork {"fork"} else {"incremental"},"distinct":distinct,"cut":cut,
                "fork_ms":fork_ms.get(2),"display":flag("--display"),"wall_ms":elapsed})
            );
            if flag("--assert-budget") {
                assert!(elapsed < 1000.0, "cut{cut}: {elapsed} ms");
            }
            if flag("--verify") {
                let full = evaluate(&edited, &mut Decider::default()).unwrap();
                assert_eq!(export(built), export(&full), "cut{cut}");
                assert_eq!(built.bodies, full.bodies);
            }
        }
        if flag("--verify") {
            for (saved, expected) in retained.iter().zip(undo_bytes) {
                assert_eq!(export(saved.built()), expected, "retained undo model");
            }
        }
        checkpoint.into_built()
    } else {
        let started = Instant::now();
        let (built, log) = with_operation_log(|| evaluate(&recipe, &mut Decider::default()));
        let built = built.unwrap();
        let elapsed = started.elapsed().as_secs_f64() * 1000.0;
        let cuts: Vec<_> = log
            .iter()
            .filter(|entry| entry.operation == "boolean")
            .collect();
        for entry in &cuts {
            println!(
                "{}",
                json!({"mode":"full","distinct":distinct,"cut":entry.id,
                "wall_ms":entry.wall_ms,"work":entry.work_after-entry.work_before,
                "decisions":entry.after.as_ref().unwrap().decisions-entry.before.decisions})
            );
        }
        println!(
            "{}",
            json!({"total_ms":elapsed,"faces":expected_faces,
            "surfaces":built.model.geom().surface_count(),"points":built.model.geom().point_count(),
            "facts":built.model.geom().facts().len()})
        );
        if flag("--assert-budget") && !distinct {
            assert!(
                cuts[29].wall_ms <= 1.5 * cuts[1].wall_ms,
                "cut30/cut2 > 1.5"
            );
        }
        built
    };
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(built.model.counts(built.result).faces, expected_faces);
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let volume = if distinct { 61600 } else { 63920 }.into();
    assert!(mass.value.volume.lo <= volume && mass.value.volume.hi >= volume);
    if flag("--fact-counterfactual") {
        counterfactual(&built);
    }
    if flag("--operand-counterfactual") {
        for input in ["base", "cut29"] {
            let mut model = built.model.clone();
            let a = built.bodies[&StepId::new(input).unwrap()];
            let b = built.bodies[&StepId::new("tool30").unwrap()];
            let start = Instant::now();
            let answer = catallax_boolean::boolean(
                &mut model,
                a,
                b,
                catallax_boolean::BoolOp::Subtract,
                &mut Decider::default(),
            )
            .unwrap();
            println!(
                "{}",
                json!({"operand":input,"same_store":true,
                "wall_ms":start.elapsed().as_secs_f64()*1000.0})
            );
            for stage in answer.report.stages {
                println!(
                    "{}",
                    json!({"operand":input,"stage":format!("{:?}",stage.stage),
                    "decisions":stage.stats.decisions,"work":stage.stats.work})
                );
            }
        }
    }
}
