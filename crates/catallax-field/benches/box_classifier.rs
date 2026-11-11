//! Shared-machine coverage and cost of the prepared box classifier, without speed gates.
//!
//! Tiles the plan's certified bounds with `CELLS` boxes per axis, classifies every box
//! with 1, 2 and 4 workers, and fails unless the three answer lists are identical.

use catallax_field::{
    Plan,
    printing::{BoxClassifier, BoxState, BoxStorage, bounds},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use serde_json::json;
use std::time::Instant;

type Corners = ([Rational; 3], [Rational; 3]);

fn grid(bounds: &[[Rational; 2]; 3], cells: i64) -> Vec<Corners> {
    let at = |axis: usize, k: i64| {
        let [lo, hi] = &bounds[axis];
        lo + &((hi - lo) * Rational::from_ratio(k, cells).expect("positive cells"))
    };
    let mut boxes = Vec::new();
    for i in 0..cells {
        for j in 0..cells {
            for k in 0..cells {
                let index = [i, j, k];
                boxes.push((
                    std::array::from_fn(|axis| at(axis, index[axis])),
                    std::array::from_fn(|axis| at(axis, index[axis] + 1)),
                ));
            }
        }
    }
    boxes
}

fn classify(
    classifier: &BoxClassifier<'_>,
    boxes: &[Corners],
    workers: usize,
) -> Result<(Vec<BoxState>, f64), Box<dyn std::error::Error + Send + Sync>> {
    let chunk = boxes.len().div_ceil(workers);
    let start = Instant::now();
    let parts = std::thread::scope(|scope| {
        let handles: Vec<_> = boxes
            .chunks(chunk)
            .map(|part| {
                scope.spawn(move || {
                    let mut storage = BoxStorage::default();
                    let mut dec = Decider::default();
                    part.iter()
                        .map(|(lo, hi)| classifier.classify(lo, hi, &mut storage, &mut dec))
                        .collect::<Result<Vec<_>, _>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("worker"))
            .collect::<Result<Vec<_>, _>>()
    })?;
    let seconds = start.elapsed().as_secs_f64();
    Ok((parts.into_iter().flatten().collect(), seconds))
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: box_classifier RECIPE CELLS".into());
    }
    let recipe = Recipe::parse(&std::fs::read_to_string(&args[1])?)?;
    let cells: i64 = args[2].parse()?;
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec)?;
    let bounds = bounds(&plan, &mut dec)?;
    let start = Instant::now();
    let classifier = BoxClassifier::new(&plan, &mut dec)?;
    let prepare_ms = start.elapsed().as_secs_f64() * 1000.0;
    let boxes = grid(&bounds, cells);
    let (states, one) = classify(&classifier, &boxes, 1)?;
    let mut timings = vec![json!({"workers": 1, "ms": one * 1000.0})];
    for workers in [2, 4] {
        let (other, seconds) = classify(&classifier, &boxes, workers)?;
        if other != states {
            return Err(format!("{workers} workers changed an answer").into());
        }
        timings.push(json!({"workers": workers, "ms": seconds * 1000.0}));
    }
    let count = |state| states.iter().filter(|s| **s == state).count();
    println!(
        "{}",
        json!({
            "recipe": recipe.name,
            "cells_per_axis": cells,
            "boxes": boxes.len(),
            "prepare_ms": prepare_ms,
            "inside": count(BoxState::Inside),
            "outside": count(BoxState::Outside),
            "unknown": count(BoxState::Unknown),
            "ns_per_box_one_worker": one * 1e9 / boxes.len() as f64,
            "timings": timings,
            "same_answers_at_1_2_4_workers": true,
        })
    );
    Ok(())
}
