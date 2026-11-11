//! The prepared box classifier of PLAN M13 step F6: no coverage refusal on the target
//! parts P1 to P4, every decided box confirmed by exact point queries, and the same
//! answers with 1, 2 and 4 workers.

use std::collections::BTreeMap;

use catallax_field::{
    Location, Plan, Workspace,
    printing::{BoxClassifier, BoxState, BoxStorage, bounds, classify_box},
};
use catallax_import::sha256_hex;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::json;

type Corners = ([Rational; 3], [Rational; 3]);

fn part(name: &str) -> Plan {
    let path = format!("{}/../../parts/{name}.json", env!("CARGO_MANIFEST_DIR"));
    let recipe = Recipe::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
    Plan::compile_checked(&recipe, &mut Decider::default()).unwrap()
}

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

/// Splits every box the classifier leaves unresolved into eight, `depth` times,
/// starting from the plan's certified bounds. Returns every visited box and its state.
fn octree(plan: &Plan, depth: usize) -> Vec<(Corners, BoxState)> {
    let mut dec = Decider::default();
    let [x, y, z] = bounds(plan, &mut dec).unwrap();
    let classifier = BoxClassifier::new(plan, &mut dec).unwrap();
    let mut storage = BoxStorage::default();
    let half = q("1/2");
    let mut level = vec![(
        [x[0].clone(), y[0].clone(), z[0].clone()],
        [x[1].clone(), y[1].clone(), z[1].clone()],
    )];
    let mut visited = Vec::new();
    for step in 0..=depth {
        let mut next = Vec::new();
        for (lo, hi) in level {
            let state = classifier
                .classify(&lo, &hi, &mut storage, &mut dec)
                .unwrap();
            if state == BoxState::Unknown && step < depth {
                let mid: [Rational; 3] = std::array::from_fn(|i| &(&lo[i] + &hi[i]) * &half);
                for corner in 0..8 {
                    let pick = |i: usize, upper: bool| {
                        if (corner >> i & 1 == 1) == upper {
                            mid[i].clone()
                        } else if upper {
                            hi[i].clone()
                        } else {
                            lo[i].clone()
                        }
                    };
                    next.push((
                        std::array::from_fn(|i| pick(i, false)),
                        std::array::from_fn(|i| pick(i, true)),
                    ));
                }
            }
            visited.push(((lo, hi), state));
        }
        level = next;
    }
    visited
}

/// A decided box is a proof about every point in it: each corner and the centre
/// must have that exact location.
fn confirm(plan: &Plan, boxes: &[(Corners, BoxState)]) -> [usize; 3] {
    let mut workspace = Workspace::default();
    let mut dec = Decider::default();
    let half = q("1/2");
    let mut counts = [0; 3];
    for ((lo, hi), state) in boxes {
        let expected = match state {
            BoxState::Inside => Location::Inside,
            BoxState::Outside => Location::Outside,
            BoxState::Unknown => {
                counts[2] += 1;
                continue;
            }
        };
        counts[usize::from(*state == BoxState::Outside)] += 1;
        let mut points: Vec<[Rational; 3]> = (0..8)
            .map(|corner| {
                std::array::from_fn(|i| {
                    if corner >> i & 1 == 0 {
                        lo[i].clone()
                    } else {
                        hi[i].clone()
                    }
                })
            })
            .collect();
        points.push(std::array::from_fn(|i| &(&lo[i] + &hi[i]) * &half));
        for [x, y, z] in points {
            let point = RPoint3::new(x, y, z);
            let query = plan
                .proven_checked(&point, &mut workspace, &mut dec)
                .unwrap();
            assert_eq!(
                query.location, expected,
                "{state:?} box {lo:?} {hi:?} at {point:?}"
            );
        }
    }
    counts
}

fn grid(plan: &Plan, cells: i64) -> Vec<Corners> {
    let bounds = bounds(plan, &mut Decider::default()).unwrap();
    let at = |axis: usize, k: i64| {
        let [lo, hi] = &bounds[axis];
        lo + &((hi - lo) * Rational::from_ratio(k, cells).unwrap())
    };
    let mut boxes = Vec::new();
    for i in 0..cells {
        for j in 0..cells {
            for k in 0..cells {
                let index = [i, j, k];
                boxes.push((
                    std::array::from_fn(|a| at(a, index[a])),
                    std::array::from_fn(|a| at(a, index[a] + 1)),
                ));
            }
        }
    }
    boxes
}

fn with_workers(plan: &Plan, boxes: &[Corners], workers: usize) -> Vec<BoxState> {
    let classifier = BoxClassifier::new(plan, &mut Decider::default()).unwrap();
    let classifier = &classifier;
    std::thread::scope(|scope| {
        let handles: Vec<_> = boxes
            .chunks(boxes.len().div_ceil(workers))
            .map(|chunk| {
                scope.spawn(move || {
                    let mut storage = BoxStorage::default();
                    let mut dec = Decider::default();
                    chunk
                        .iter()
                        .map(|(lo, hi)| {
                            classifier.classify(lo, hi, &mut storage, &mut dec).unwrap()
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|handle| handle.join().unwrap())
            .collect()
    })
}

fn covered(name: &str, depth: usize) {
    let plan = part(name);
    let boxes = octree(&plan, depth);
    let [inside, outside, unknown] = confirm(&plan, &boxes);
    println!("{name}: {inside} inside, {outside} outside, {unknown} unresolved boxes");
    assert!(
        inside > 0 && outside > 0,
        "{name} decides boxes on both sides"
    );
}

#[test]
fn p1_bracket_is_covered_and_every_decided_box_is_confirmed() {
    covered("P1-bracket", 4);
}

#[test]
fn p2_collar_is_covered_and_every_decided_box_is_confirmed() {
    covered("P2-collar", 4);
}

#[test]
fn p3_filleted_boss_is_covered_and_every_decided_box_is_confirmed() {
    covered("P3-filleted-boss", 5);
}

#[test]
fn p4_shelled_tray_is_covered_and_every_decided_box_is_confirmed() {
    covered("P4-shelled-tray", 4);
}

#[test]
fn answers_are_the_same_with_one_two_and_four_workers() {
    for name in [
        "P1-bracket",
        "P2-collar",
        "P3-filleted-boss",
        "P4-shelled-tray",
    ] {
        let plan = part(name);
        let boxes = grid(&plan, 9);
        let one = with_workers(&plan, &boxes, 1);
        assert_eq!(with_workers(&plan, &boxes, 2), one, "{name}");
        assert_eq!(with_workers(&plan, &boxes, 4), one, "{name}");
        // The one-box entry point prepares the same classifier.
        for ((lo, hi), state) in boxes.iter().zip(&one).step_by(37) {
            assert_eq!(
                classify_box(&plan, lo, hi, &mut Decider::default()).unwrap(),
                *state
            );
        }
    }
}

#[test]
fn an_inverted_box_still_refuses() {
    let plan = part("P2-collar");
    let classifier = BoxClassifier::new(&plan, &mut Decider::default()).unwrap();
    let refused = classifier
        .classify(
            &[q("1"), q("0"), q("0")],
            &[q("0"), q("1"), q("1")],
            &mut BoxStorage::default(),
            &mut Decider::default(),
        )
        .unwrap_err();
    assert_eq!(refused.code, "field.print-box");
}

fn cuboid(lo: [i64; 3], hi: [i64; 3]) -> Vec<u8> {
    let corner = |i: usize| -> [i64; 3] {
        std::array::from_fn(|a| if i >> a & 1 == 0 { lo[a] } else { hi[a] })
    };
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut s = String::from("solid part\n");
    for f in faces {
        for t in [[f[0], f[1], f[2]], [f[0], f[2], f[3]]] {
            s.push_str("facet normal 0 0 0\nouter loop\n");
            for v in t {
                let [x, y, z] = corner(v);
                s.push_str(&format!("vertex {x} {y} {z}\n"));
            }
            s.push_str("endloop\nendfacet\n");
        }
    }
    s.push_str("endsolid part\n");
    s.into_bytes()
}

#[test]
fn a_mesh_in_a_boolean_is_covered() {
    let bytes = cuboid([0, 0, 0], [2, 2, 2]);
    let hash = sha256_hex(&bytes);
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let recipe = Recipe::parse(
        &json!({
            "format": "catallax-recipe/4",
            "name": "mesh boolean",
            "steps": [
                {"id": "scan", "op": "import_mesh", "file": {"name": "part.stl", "sha256": hash},
                 "pose": identity, "scale": "1"},
                {"id": "block", "op": "block", "size": ["1", "1", "4"], "pose": identity},
                {"id": "cut", "op": "boolean", "kind": "subtract", "target": "scan", "tool": "block"}
            ],
            "result": "cut"
        })
        .to_string(),
    )
    .unwrap();
    let files = BTreeMap::from([(hash, &bytes[..])]);
    let plan = Plan::compile_with_files_checked(&recipe, &files, &mut Decider::default()).unwrap();
    let state = |lo: [&str; 3], hi: [&str; 3]| {
        classify_box(&plan, &lo.map(q), &hi.map(q), &mut Decider::default()).unwrap()
    };
    assert_eq!(
        state(["5/4", "5/4", "1/2"], ["7/4", "7/4", "3/2"]),
        BoxState::Inside
    );
    assert_eq!(
        state(["1/4", "1/4", "1/2"], ["3/4", "3/4", "3/2"]),
        BoxState::Outside
    );
    assert_eq!(state(["3", "0", "0"], ["4", "1", "1"]), BoxState::Outside);
    assert_eq!(
        state(["3/2", "1/2", "1/2"], ["5/2", "3/2", "3/2"]),
        BoxState::Unknown
    );
    let boxes = octree(&plan, 3);
    let [inside, outside, _] = confirm(&plan, &boxes);
    assert!(inside > 0 && outside > 0);
}
