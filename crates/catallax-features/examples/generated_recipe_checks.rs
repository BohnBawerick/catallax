//! Bounded R3 recipe checks; see tools/edge-hunt/GENERATED-RECIPES.md.
use catallax_features::{BuildError, evaluate};
use catallax_measure::{Accuracy, VolumeAreaCentroid, volume_area_centroid};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{collections::BTreeMap, fs, path::PathBuf};

const WORK: u64 = 20_000_000;

#[derive(Clone, Debug)]
struct Case {
    lo: [[i64; 3]; 2],
    size: [[i64; 3]; 2],
    denominator: i64,
    kind: usize,
}

// SplitMix64, per-index stream: prefixes do not depend on worker count or count.
fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn generate(seed: u64, index: usize) -> Case {
    let mut state = seed ^ (index as u64).wrapping_mul(0xd1342543de82ef95);
    let mut c = Case {
        lo: [[0; 3]; 2],
        size: [[0; 3]; 2],
        denominator: [1, 2, 8][(index / 24) % 3],
        kind: index % 4,
    };
    for axis in 0..3 {
        c.lo[0][axis] = (next(&mut state) % 9) as i64 - 4;
        c.size[0][axis] = 2 + (next(&mut state) % 5) as i64;
        c.size[1][axis] = 1 + (next(&mut state) % 5) as i64;
        c.lo[1][axis] = c.lo[0][axis] + (next(&mut state) % 7) as i64 - 3;
    }
    match (index / 4) % 6 {
        0 => {
            c.lo[1] = c.lo[0];
            c.size[1] = c.size[0];
        } // identity / empty
        1 => {
            c.lo[1] = c.lo[0];
            c.lo[1][0] += c.size[0][0] + 1;
        } // gap
        2 => {
            c.lo[1] = c.lo[0];
            c.lo[1][0] += c.size[0][0];
        } // face touch
        3 => {
            c.lo[1] = c.lo[0];
            c.size[1] = [1; 3];
        } // containment
        _ => {} // general overlapping or separated boxes
    }
    c
}
fn rat(n: i64, d: i64) -> Rational {
    Rational::parse(&format!("{n}/{d}")).unwrap()
}
impl Case {
    fn recipe(&self, shift: i64, swap: bool) -> String {
        let block = |i: usize, id: &str| {
            json!({"id":id,"op":"block",
                "size":self.size[i].map(|v| rat(v,self.denominator).to_string()),
                "pose":{"rotation":["1","0","0","0"],
                    "translation":self.lo[i].map(|v| rat(v+shift,self.denominator).to_string())}})
        };
        let mut steps = vec![block(0, "a")];
        if self.kind != 0 {
            steps.push(block(1, "b"));
            steps.push(json!({"id":"result","op":"boolean",
                "kind":(["", "union","intersect","subtract"][self.kind]),
                "target":if swap {"b"} else {"a"},"tool":if swap {"a"} else {"b"}}));
        }
        json!({"format":"catallax-recipe/3","name":"R3 generated boxes",
            "steps":steps,"result":if self.kind==0 {"a"} else {"result"}})
        .to_string()
    }
    fn metadata(&self) -> Value {
        json!({"lo":self.lo,"size":self.size,"denominator":self.denominator,"kind":self.kind})
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Expected {
    volume: Rational,
    area: Rational,
    centroid: Option<[Rational; 3]>,
}

// Independent integer arrangement: up to 27 open cells, with exact membership
// truth tables. Sum volumes and first moments; count only exposed cell faces.
// No kernel construction, topology, integration or enclosure routines are used.
fn oracle(c: &Case, shift: i64) -> Expected {
    let cuts: [Vec<i64>; 3] = std::array::from_fn(|a| {
        let mut v = vec![
            c.lo[0][a],
            c.lo[0][a] + c.size[0][a],
            c.lo[1][a],
            c.lo[1][a] + c.size[1][a],
        ];
        v.sort();
        v.dedup();
        v
    });
    let occupied = |p: [usize; 3]| {
        let inside = |b: usize| {
            (0..3).all(|a| {
                cuts[a][p[a]] >= c.lo[b][a] && cuts[a][p[a] + 1] <= c.lo[b][a] + c.size[b][a]
            })
        };
        match c.kind {
            0 => inside(0),
            1 => inside(0) || inside(1),
            2 => inside(0) && inside(1),
            3 => inside(0) && !inside(1),
            _ => unreachable!(),
        }
    };
    let mut volume = 0;
    let mut area = 0;
    let mut twice_moment = [0; 3];
    for x in 0..cuts[0].len() - 1 {
        for y in 0..cuts[1].len() - 1 {
            for z in 0..cuts[2].len() - 1 {
                let p = [x, y, z];
                if !occupied(p) {
                    continue;
                }
                let widths: [i64; 3] = std::array::from_fn(|a| cuts[a][p[a] + 1] - cuts[a][p[a]]);
                let v: i64 = widths.iter().product();
                volume += v;
                for a in 0..3 {
                    twice_moment[a] += v * (cuts[a][p[a]] + cuts[a][p[a] + 1] + 2 * shift);
                    for upper in [false, true] {
                        let mut q = p;
                        let outside = if upper {
                            q[a] += 1;
                            q[a] + 1 >= cuts[a].len()
                        } else if q[a] == 0 {
                            true
                        } else {
                            q[a] -= 1;
                            false
                        };
                        if outside || !occupied(q) {
                            area += widths[(a + 1) % 3] * widths[(a + 2) % 3];
                        }
                    }
                }
            }
        }
    }
    Expected {
        volume: rat(volume, c.denominator.pow(3)),
        area: rat(area, c.denominator.pow(2)),
        centroid: (volume != 0).then(|| twice_moment.map(|m| rat(m, 2 * volume * c.denominator))),
    }
}
fn check(p: &VolumeAreaCentroid, e: &Expected) -> Result<(), String> {
    if !p.volume.contains(&e.volume) {
        return Err("wrong-proven-volume".into());
    }
    if !p.area.contains(&e.area) {
        return Err("wrong-proven-area".into());
    }
    match (&p.centroid, &e.centroid) {
        (None, None) => {}
        (Some(p), Some(e)) if (0..3).all(|a| p[a].contains(&e[a])) => {}
        _ => return Err("wrong-proven-centroid".into()),
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Outcome {
    status: String,
    detail: String,
    work: u64,
}
fn run(recipe: &str, expected: &Expected) -> Outcome {
    run_with_limit(recipe, expected, WORK)
}
fn run_with_limit(recipe: &str, expected: &Expected, limit: u64) -> Outcome {
    let caught = std::panic::catch_unwind(|| {
        with_operation_budget(limit, || {
            let recipe =
                Recipe::parse(recipe).map_err(|e| ("invalid-recipe".into(), format!("{e:?}")))?;
            if Recipe::parse(&recipe.to_text()).unwrap() != recipe {
                return Err(("serialization".into(), "round trip changed recipe".into()));
            }
            let mut dec = Decider::default();
            let built = evaluate(&recipe, &mut dec).map_err(|e| match e {
                BuildError::Geometry { step, refusal } => (
                    format!("refused:build:{}", refusal.code),
                    format!("{step}:{}", refusal.detail),
                ),
                BuildError::Recipe(e) => ("invalid-recipe".into(), format!("{e:?}")),
            })?;
            let faults = built.model.check(built.result);
            if !faults.is_empty() {
                return Err(("invalid-topology".into(), format!("{faults:?}")));
            }
            let p =
                volume_area_centroid(&built.model, built.result, &Accuracy::default(), &mut dec)
                    .map_err(|e| (format!("refused:measure:{}", e.code), e.detail.to_string()))?;
            check(&p.value, expected)
                .map_err(|e| (e, format!("expected={expected:?}; actual={p:?}")))?;
            Ok(format!("{p:?}"))
        })
    });
    match caught {
        Err(_) => Outcome {
            status: "panic".into(),
            detail: "kernel or harness panic".into(),
            work: 0,
        },
        Ok((result, work)) => {
            let (status, detail) = match result {
                Ok(d) => ("checked".into(), d),
                Err(e) => e,
            };
            // A result returned after exhaustion must not be counted as checked.
            let status = if work.exhausted && status == "checked" {
                "budget-contract".into()
            } else {
                status
            };
            Outcome {
                status,
                detail,
                work: work.used,
            }
        }
    }
}
fn failing(o: &Outcome) -> bool {
    o.status != "checked" && !o.status.starts_with("refused:")
}
fn variants(c: &Case) -> Vec<(String, Expected)> {
    let mut v = vec![
        (c.recipe(0, false), oracle(c, 0)),
        (c.recipe(7, false), oracle(c, 7)),
    ];
    if c.kind == 1 || c.kind == 2 {
        v.push((c.recipe(0, true), oracle(c, 0)));
    }
    v
}
fn batch(cases: &[Case], workers: usize) -> Vec<Vec<Outcome>> {
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..workers)
            .map(|w| {
                s.spawn(move || {
                    (w..cases.len())
                        .step_by(workers)
                        .map(|i| {
                            (
                                i,
                                variants(&cases[i]).iter().map(|(r, e)| run(r, e)).collect(),
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let mut rows = vec![Vec::new(); cases.len()];
        for h in handles {
            for (i, row) in h.join().expect("worker panicked") {
                rows[i] = row;
            }
        }
        rows
    })
}
// Bounded, class-preserving greedy reduction. Not a claim of global minimality.
fn shrink(mut c: Case, variant: usize, class: &str) -> Case {
    for field in 0..12 {
        let mut smaller = c.clone();
        let b = (field % 6) / 3;
        let a = field % 3;
        if field < 6 {
            smaller.lo[b][a] = 0;
        } else {
            smaller.size[b][a] = 1;
        }
        let v = variants(&smaller);
        if run(&v[variant].0, &v[variant].1).status == class {
            c = smaller;
        }
    }
    c
}
fn write_new(path: PathBuf, value: &str) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    f.write_all(value.as_bytes())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 4 {
        return Err(
            "usage: generated_recipe_checks SEED START COUNT NEW_OUTPUT_DIRECTORY (COUNT 1..10000)"
                .into(),
        );
    }
    let seed = args[0].parse::<u64>()?;
    let start = args[1].parse::<usize>()?;
    let count = args[2].parse::<usize>()?;
    if !(1..=10000).contains(&count) || start.checked_add(count).is_none() {
        return Err("invalid range".into());
    }
    let out = PathBuf::from(&args[3]);
    fs::create_dir(&out)?;
    let cases: Vec<_> = (start..start + count).map(|i| generate(seed, i)).collect();
    // Save all portable inputs BEFORE execution, including those behind a crash.
    for (i, c) in cases.iter().enumerate() {
        for (v, (r, _)) in variants(c).iter().enumerate() {
            write_new(out.join(format!("case-{}-{v}.json", start + i)), r)?;
        }
    }
    let baseline = batch(&cases, 1);
    let mut failures = 0;
    let mut counts = BTreeMap::<String, usize>::new();
    let mut rows = Vec::new();
    for (i, result) in baseline.iter().enumerate() {
        for (v, o) in result.iter().enumerate() {
            *counts.entry(o.status.clone()).or_default() += 1;
            rows.push(json!({"index":start+i,"variant":v,"parameters":cases[i].metadata(),"status":o.status,"detail":o.detail,"work":o.work}));
            if failing(o) {
                failures += 1;
                let small = shrink(cases[i].clone(), v, &o.status);
                write_new(
                    out.join(format!("reduced-{}-{v}.json", start + i)),
                    &variants(&small)[v].0,
                )?;
                rows.push(json!({"reduction_of":[start+i,v],"class":o.status,"parameters":small.metadata()}));
            }
        }
    }
    for workers in [2, 4] {
        let other = batch(&cases, workers);
        for (i, (a, b)) in baseline.iter().zip(&other).enumerate() {
            if a != b {
                failures += 1;
                rows.push(json!({"index":start+i,"status":"worker-mismatch","workers":workers,"one":format!("{a:?}"),"other":format!("{b:?}")}));
            }
        }
    }
    // Non-vacuity: every requested operation must have a successful checked case.
    for kind in 0..4 {
        if cases.iter().any(|c| c.kind == kind)
            && !cases
                .iter()
                .zip(&baseline)
                .any(|(c, r)| c.kind == kind && r[0].status == "checked")
        {
            failures += 1;
            rows.push(json!({"status":"coverage-failure","kind":kind}));
        }
    }
    let summary = json!({"format":"catallax-generated-recipe-checks/1","seed":seed,"start":start,"count":count,
        "workers":[1,2,4],"work_limit":WORK,"outcomes":counts,"failures":failures,
        "evidence":"finite checks, not a proof of the kernel","rows":rows});
    write_new(
        out.join("results.json"),
        &serde_json::to_string_pretty(&summary)?,
    )?;
    println!("seed={seed} start={start} cases={count} outcomes={counts:?} failures={failures}");
    if failures != 0 {
        return Err("R3 checks failed; see results.json and portable recipes".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_measure::Enclosure;
    fn exact(v: Rational) -> Enclosure {
        Enclosure {
            lo: v.clone(),
            hi: v,
        }
    }
    #[test]
    fn analytic_controls_and_fault_injection() {
        let mut c = Case {
            lo: [[0; 3]; 2],
            size: [[2, 3, 4]; 2],
            denominator: 1,
            kind: 0,
        };
        let e = oracle(&c, 0);
        assert_eq!(
            e,
            Expected {
                volume: 24.into(),
                area: 52.into(),
                centroid: Some([1.into(), rat(3, 2), 2.into()])
            }
        );
        c.kind = 3;
        assert_eq!(
            oracle(&c, 0),
            Expected {
                volume: 0.into(),
                area: 0.into(),
                centroid: None
            }
        );
        c.kind = 1;
        c.lo[1][0] = 2;
        assert_eq!(oracle(&c, 0).area, Rational::from(80)); // 4 x 3 x 4 block
        c.lo[1][0] = 3;
        assert_eq!(oracle(&c, 0).area, Rational::from(104)); // disjoint
        let mut p = VolumeAreaCentroid {
            volume: exact(e.volume.clone()),
            area: exact(e.area.clone()),
            centroid: e.centroid.clone().map(|p| p.map(exact)),
        };
        assert!(check(&p, &e).is_ok());
        p.volume = exact(23.into());
        assert_eq!(check(&p, &e).unwrap_err(), "wrong-proven-volume");
        p.volume = exact(e.volume.clone());
        p.area = exact(51.into());
        assert_eq!(check(&p, &e).unwrap_err(), "wrong-proven-area");
        p.area = exact(e.area.clone());
        p.centroid = None;
        assert_eq!(check(&p, &e).unwrap_err(), "wrong-proven-centroid");
    }
    #[test]
    fn cavity_has_inner_area_and_subtracted_first_moment() {
        let c = Case {
            lo: [[0; 3], [1; 3]],
            size: [[4; 3], [1; 3]],
            denominator: 1,
            kind: 3,
        };
        assert_eq!(
            oracle(&c, 0),
            Expected {
                volume: 63.into(),
                area: 102.into(),
                centroid: Some(std::array::from_fn(|_| rat(253, 126))),
            }
        );
    }

    #[test]
    fn oracle_translation_and_inclusion_exclusion() {
        for index in 0..96 {
            let mut c = generate(1414213562, index);
            c.kind = 1;
            let union = oracle(&c, 0);
            let translated = oracle(&c, 7);
            assert_eq!(union.volume, translated.volume);
            assert_eq!(union.area, translated.area);
            if let (Some(a), Some(b)) = (union.centroid, translated.centroid) {
                for axis in 0..3 {
                    assert_eq!(&a[axis] + &rat(7, c.denominator), b[axis]);
                }
            }
            c.kind = 2;
            let intersection = oracle(&c, 0);
            let volumes: [i64; 2] = c.size.map(|s| s.iter().product());
            assert_eq!(
                union.volume + intersection.volume,
                rat(volumes[0] + volumes[1], c.denominator.pow(3))
            );
        }
    }

    #[test]
    fn exhausted_request_is_a_typed_refusal_not_an_answer() {
        let c = generate(20261111, 17);
        let result = run_with_limit(&c.recipe(0, false), &oracle(&c, 0), 0);
        assert_eq!(result.status, "refused:build:numbers.operation-budget");
        assert!(!failing(&result));
    }

    #[test]
    fn generated_smoke_workers_and_non_vacuity() {
        let cases: Vec<_> = (0..24).map(|i| generate(20261111, i)).collect();
        let one = batch(&cases, 1);
        for workers in [2, 4] {
            assert_eq!(one, batch(&cases, workers));
        }
        for (c, row) in cases.iter().zip(one) {
            assert!(row.iter().all(|r| !failing(r)), "{c:?}: {row:?}");
            assert_eq!(row[0].status, "checked", "{c:?}: {row:?}");
        }
    }
}
