//! F4/F5 evidence: scalar signs against deterministic `Plan::fast`, threaded
//! bits and signs against the scalar preview, and informational per-point times.
//!
//! `cargo bench -p catallax-field --bench fast_preview [-- POINTS]` samples
//! POINTS (default 2,000,000) seeded uniform points in each part's box grown by
//! a tenth on every side, plus a twentieth as many points next to the
//! deterministic zero set, found by bisection. Every point whose two signs
//! differ is printed and classified by the proven evaluator. Times are from
//! whatever machine runs this and are never a timing pass.

use catallax_field::{Location, Plan, UnprovedBatch, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::json;
use std::time::Instant;

const PARTS: [(&str, &str, [f64; 3]); 2] = [
    ("P1", "P1-bracket.json", [60.0, 40.0, 30.0]),
    ("P3", "P3-filleted-boss.json", [40.0, 30.0, 22.0]),
];
const ROUNDS: usize = 5;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = std::env::args()
        .skip(1)
        .find(|a| !a.starts_with('-'))
        .map_or(Ok(2_000_000), |a| a.parse())?;
    let parts = concat!(env!("CARGO_MANIFEST_DIR"), "/../../parts/");
    for (name, file, size) in PARTS {
        let recipe = Recipe::parse(&std::fs::read_to_string(format!("{parts}{file}"))?)?;
        let plan = Plan::compile(&recipe)?;
        let points = sample(n, size);
        let result = run(name, &plan, &points)?;
        println!("{}", serde_json::to_string(&result)?);
        for set in ["threaded_uniform", "threaded_near_surface"] {
            for row in result[set].as_array().expect("threaded records") {
                if row["bit_mismatches"] != 0 || row["sign_mismatches"] != 0 {
                    return Err("threaded preview differs from scalar preview".into());
                }
            }
        }
    }
    Ok(())
}

fn sample(n: usize, size: [f64; 3]) -> Vec<[f64; 3]> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| std::array::from_fn(|i| size[i] * (1.2 * next() - 0.1)))
        .collect()
}

fn sign(v: f64) -> i8 {
    if v < 0.0 {
        -1
    } else if v > 0.0 {
        1
    } else {
        0
    }
}

fn run(
    name: &str,
    plan: &Plan,
    points: &[[f64; 3]],
) -> Result<serde_json::Value, Box<dyn std::error::Error>> {
    let uniform = compare(plan, points)?;
    let near_points = near_surface(plan, points, points.len() / 20)?;
    let near = compare(plan, &near_points)?;
    let threaded_uniform = threaded(plan, points)?;
    let threaded_near_surface = threaded(plan, &near_points)?;
    let fast = time(points, |p, w| plan.fast(p, w).map(|v| v.value))?;
    let processor = time(points, |p, w| plan.unproved_preview(p, w).map(|v| v.value))?;
    let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
    Ok(json!({
        "part": name,
        "proof": "unproved",
        "timing": "shared machine, not a timing pass",
        "load_average": load.split_whitespace().take(3).collect::<Vec<_>>(),
        "uniform": uniform,
        "near_surface": near,
        "threaded_uniform": threaded_uniform,
        "threaded_near_surface": threaded_near_surface,
        "fast_ns_per_point": fast,
        "processor_ns_per_point": processor,
    }))
}

fn threaded(
    plan: &Plan,
    points: &[[f64; 3]],
) -> Result<Vec<serde_json::Value>, Box<dyn std::error::Error>> {
    let mut w = Workspace::default();
    let scalar: Vec<_> = points
        .iter()
        .map(|p| plan.unproved_preview(*p, &mut w))
        .collect::<Result<_, _>>()?;
    let mut batch = UnprovedBatch::default();
    let mut rows = Vec::new();
    for workers in [1, 2, 4, 8] {
        let values = plan.unproved_preview_batch(points, workers, &mut batch)?;
        let mut bit_mismatches = 0;
        let mut sign_mismatches = 0;
        for (a, b) in values.iter().zip(&scalar) {
            bit_mismatches += usize::from(a.value.to_bits() != b.value.to_bits());
            sign_mismatches += usize::from(sign(a.value) != sign(b.value));
        }
        let mut rounds = Vec::new();
        for round in 0..=ROUNDS {
            let started = Instant::now();
            std::hint::black_box(plan.unproved_preview_batch(points, workers, &mut batch)?);
            if round > 0 {
                rounds.push(started.elapsed().as_nanos() as f64 / points.len() as f64);
            }
        }
        rounds.sort_by(f64::total_cmp);
        rows.push(json!({
            "workers": workers,
            "points": points.len(),
            "scalar_zero": scalar.iter().filter(|v| v.value == 0.0).count(),
            "bit_mismatches": bit_mismatches,
            "sign_mismatches": sign_mismatches,
            "ns_per_point": rounds[ROUNDS / 2],
        }));
    }
    Ok(rows)
}

// Sign agreement on one point set. Every mismatch is listed with the proven
// location at that exact binary64 point.
fn compare(
    plan: &Plan,
    points: &[[f64; 3]],
) -> Result<serde_json::Value, Box<dyn std::error::Error>> {
    let mut w = Workspace::default();
    let (mut inside, mut outside, mut zero) = (0usize, 0usize, 0usize);
    let mut max_difference = 0.0f64;
    let mut max_value = 0.0f64;
    let mut deterministic = Vec::with_capacity(points.len());
    let mut mismatches = Vec::new();
    let mut max_mismatch_value = 0.0f64;
    // Mismatches whose proven location is inside or outside, by which
    // evaluator's sign agrees with it (neither, when one of them is zero);
    // and those proven on the boundary or left indeterminate.
    let (mut processor_right, mut deterministic_right, mut neither, mut unsettled) = (0, 0, 0, 0);
    for p in points {
        let d = plan.fast(*p, &mut w)?.value;
        let h = plan.unproved_preview(*p, &mut w)?.value;
        deterministic.push(d);
        max_difference = max_difference.max((h - d).abs());
        max_value = max_value.max(d.abs());
        match sign(d) {
            -1 => inside += 1,
            1 => outside += 1,
            _ => zero += 1,
        }
        if sign(d) != sign(h) {
            let exact = RPoint3::new(
                Rational::from_f64_exact(p[0]).expect("finite"),
                Rational::from_f64_exact(p[1]).expect("finite"),
                Rational::from_f64_exact(p[2]).expect("finite"),
            );
            let proven: Location = plan.proven(&exact, &mut w).location;
            max_mismatch_value = max_mismatch_value.max(d.abs());
            let truth = match proven {
                Location::Inside => -1,
                Location::Outside => 1,
                _ => 0,
            };
            if truth == 0 {
                unsettled += 1;
            } else if truth == sign(h) {
                processor_right += 1;
            } else if truth == sign(d) {
                deterministic_right += 1;
            } else {
                neither += 1;
            }
            mismatches.push(json!({"point": p, "deterministic": d, "processor": h,
                "proven": format!("{proven:?}")}));
        }
    }
    // Only points this close to the deterministic zero set can change sign
    // under the largest difference observed on this set.
    let band = deterministic
        .iter()
        .filter(|d| d.abs() <= max_difference)
        .count();
    Ok(json!({
        "points": points.len(),
        "deterministic_signs": {"inside": inside, "outside": outside, "zero": zero},
        "max_abs_deterministic_mm": max_value,
        "max_abs_difference_mm": max_difference,
        "points_within_max_difference_of_zero": band,
        "sign_mismatches": mismatches.len(),
        "max_abs_deterministic_at_mismatch_mm": max_mismatch_value,
        "mismatch_proven_sign_matches_processor": processor_right,
        "mismatch_proven_sign_matches_deterministic": deterministic_right,
        "mismatch_proven_sign_matches_neither": neither,
        "mismatch_proven_boundary_or_indeterminate": unsettled,
        "mismatches": mismatches,
    }))
}

// Up to `n` points next to the deterministic zero set: consecutive sample
// points of opposite sign are bisected with `Plan::fast` until the two ends
// are adjacent binary64 points along the segment, and both ends are kept.
fn near_surface(
    plan: &Plan,
    points: &[[f64; 3]],
    n: usize,
) -> Result<Vec<[f64; 3]>, Box<dyn std::error::Error>> {
    let mut w = Workspace::default();
    let mut found = Vec::new();
    let bits = |p: [f64; 3]| p.map(f64::to_bits);
    for &[mut a, mut b] in points.as_chunks::<2>().0 {
        if found.len() >= n {
            break;
        }
        let sa = sign(plan.fast(a, &mut w)?.value);
        if sa == 0 || sa == sign(plan.fast(b, &mut w)?.value) {
            continue;
        }
        loop {
            let m: [f64; 3] = std::array::from_fn(|i| 0.5 * (a[i] + b[i]));
            if bits(m) == bits(a) || bits(m) == bits(b) {
                break;
            }
            if sign(plan.fast(m, &mut w)?.value) == sa {
                a = m;
            } else {
                b = m;
            }
        }
        found.push(a);
        found.push(b);
    }
    Ok(found)
}

// Median of ROUNDS passes after one warm-up pass, in nanoseconds per point.
fn time(
    points: &[[f64; 3]],
    mut f: impl FnMut([f64; 3], &mut Workspace) -> Result<f64, catallax_numbers::Refusal>,
) -> Result<f64, Box<dyn std::error::Error>> {
    let mut w = Workspace::default();
    let mut rounds = Vec::new();
    for round in 0..=ROUNDS {
        let started = Instant::now();
        let mut checksum = 0.0;
        for &p in points {
            checksum += f(p, &mut w)?;
        }
        std::hint::black_box(checksum);
        if round > 0 {
            rounds.push(started.elapsed().as_nanos() as f64 / points.len() as f64);
        }
    }
    rounds.sort_by(f64::total_cmp);
    Ok(rounds[ROUNDS / 2])
}
