//! Checks the graphics-card first paint against the processor preview on the
//! target parts P1 to P4. Everything it reports is unproved display evidence:
//! agreement, mismatches and timings are information, not a pass.
use catallax_field::Plan;
use catallax_field::display::{PreviewMesh, PreviewOptions, preview_plan};
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use catallax_preview::batch::{Evaluator, Grid, Hardware, MAX_DISPATCH, Options, Samples, Side};
use catallax_preview::first_paint::{MAX_POINTS, SampledMesh, check_size, grid, surface_nets};
use catallax_preview::gpu::Gpu;
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{collections::HashMap, error::Error, fs, path::Path, time::Instant};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const PARTS: [(&str, &str); 4] = [
    ("P1", include_str!("../../../../parts/P1-bracket.json")),
    ("P2", include_str!("../../../../parts/P2-collar.json")),
    (
        "P3",
        include_str!("../../../../parts/P3-filleted-boss.json"),
    ),
    ("P4", include_str!("../../../../parts/P4-shelled-tray.json")),
];
/// The application allowance for the processor chart preview.
const PREVIEW_ALLOWANCE: u64 = 20_000_000;
/// Mismatching grid points listed in full per pitch.
const LISTED: usize = 64;
/// Seeded uniform points for the arbitrary-point route.
const POINTS: usize = 1_000_000;
/// A dispatch size that splits every grid unevenly, for the partition check.
const ODD_DISPATCH: usize = 100_003;
/// Nearest-triangle searches stop beyond this distance in millimetres.
const SEARCH_CAP: f64 = 5.0;
/// Farthest points listed per distance direction.
const FARTHEST: usize = 16;

fn ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn bits(samples: &Samples) -> Vec<u64> {
    samples.values.iter().map(|v| v.to_bits()).collect()
}

fn hardware(h: &Hardware) -> Value {
    match h {
        Hardware::GraphicsCard { adapter } => json!({"route":"graphics card","adapter":adapter}),
        Hardware::Processor { workers, reason } => {
            json!({"route":"processor","workers":workers,"reason":reason})
        }
    }
}

/// Sign and value agreement between two sample sets at the same points.
fn agreement(points: impl Fn(usize) -> [f64; 3], gpu: &Samples, cpu: &Samples) -> Value {
    let mut mismatches = 0usize;
    let mut beyond = 0usize;
    let mut max_difference = 0.0f64;
    let mut listed = Vec::new();
    let mut zero = [0usize; 2];
    for i in 0..cpu.values.len() {
        let (g, c) = (gpu.values[i], cpu.values[i]);
        max_difference = max_difference.max((g - c).abs());
        zero[0] += usize::from(gpu.side(i) == Side::Zero);
        zero[1] += usize::from(cpu.side(i) == Side::Zero);
        if gpu.side(i) != cpu.side(i) {
            mismatches += 1;
            beyond += usize::from(c.abs() > 1e-4);
            if listed.len() < LISTED {
                listed.push(json!({"index":i,"point":points(i),"gpu":g,"processor":c}));
            }
        }
    }
    json!({
        "points": cpu.values.len(),
        "side_mismatches": mismatches,
        "side_mismatches_beyond_0_0001": beyond,
        "gpu_zero": zero[0],
        "processor_zero": zero[1],
        "max_absolute_difference": max_difference,
        "listed_mismatches": listed,
    })
}

/// Nearest-triangle distances through a uniform bucket grid.
struct Nearest<'a> {
    positions: &'a [[f64; 3]],
    triangles: Vec<[usize; 3]>,
    cell: f64,
    buckets: HashMap<[i64; 3], Vec<u32>>,
}

impl<'a> Nearest<'a> {
    fn new(positions: &'a [[f64; 3]], triangles: Vec<[usize; 3]>, cell: f64) -> Self {
        let mut buckets: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        for (t, tri) in triangles.iter().enumerate() {
            let p = tri.map(|v| positions[v]);
            let lo: [i64; 3] = std::array::from_fn(|a| {
                (p.iter().map(|q| q[a]).fold(f64::INFINITY, f64::min) / cell).floor() as i64
            });
            let hi: [i64; 3] = std::array::from_fn(|a| {
                (p.iter().map(|q| q[a]).fold(f64::NEG_INFINITY, f64::max) / cell).floor() as i64
            });
            for x in lo[0]..=hi[0] {
                for y in lo[1]..=hi[1] {
                    for z in lo[2]..=hi[2] {
                        buckets.entry([x, y, z]).or_default().push(t as u32);
                    }
                }
            }
        }
        Self {
            positions,
            triangles,
            cell,
            buckets,
        }
    }

    /// Distance to the nearest triangle, or `None` beyond [`SEARCH_CAP`].
    fn distance(&self, p: [f64; 3]) -> Option<f64> {
        let home = p.map(|c| (c / self.cell).floor() as i64);
        let rings = (SEARCH_CAP / self.cell).ceil() as i64 + 1;
        let mut best = f64::INFINITY;
        for r in 0..=rings {
            for x in -r..=r {
                for y in -r..=r {
                    for z in -r..=r {
                        if x.abs().max(y.abs()).max(z.abs()) != r {
                            continue;
                        }
                        let Some(list) = self.buckets.get(&[home[0] + x, home[1] + y, home[2] + z])
                        else {
                            continue;
                        };
                        for &t in list {
                            let tri = self.triangles[t as usize].map(|v| self.positions[v]);
                            best = best.min(point_triangle(p, tri));
                        }
                    }
                }
            }
            // Anything in a later ring is at least r cells away.
            if best <= r as f64 * self.cell {
                break;
            }
        }
        (best <= SEARCH_CAP).then_some(best)
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean distance from `p` to a triangle (Ericson's closest-point regions).
fn point_triangle(p: [f64; 3], [a, b, c]: [[f64; 3]; 3]) -> f64 {
    let at = |q: [f64; 3]| dot(sub(p, q), sub(p, q)).sqrt();
    let along =
        |o: [f64; 3], d: [f64; 3], t: f64| [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return at(a);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return at(b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return at(along(a, ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return at(c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return at(along(a, ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return at(along(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))));
    }
    let denom = va + vb + vc;
    if denom.abs() < f64::MIN_POSITIVE {
        // A degenerate triangle: its edges are the whole set.
        return at(a).min(at(b)).min(at(c));
    }
    let v = vb / denom;
    let w = vc / denom;
    at([
        a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w,
    ])
}

/// Distances from every point to a mesh, in input order, on `workers` threads.
fn distances(points: &[[f64; 3]], to: &Nearest<'_>, workers: usize) -> Vec<Option<f64>> {
    let chunk = points.len().div_ceil(workers).max(1);
    std::thread::scope(|scope| {
        let handles: Vec<_> = points
            .chunks(chunk)
            .map(|part| {
                scope.spawn(move || part.iter().map(|p| to.distance(*p)).collect::<Vec<_>>())
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("distance worker"))
            .collect()
    })
}

fn statistics(points: &[[f64; 3]], values: &[Option<f64>], pitch: f64) -> Value {
    let mut order: Vec<usize> = (0..values.len()).collect();
    // Beyond the search cap sorts first, then by distance; ties by index.
    order.sort_by(|&a, &b| {
        let key = |i: usize| values[i].unwrap_or(f64::INFINITY);
        key(b).total_cmp(&key(a)).then(a.cmp(&b))
    });
    let farthest: Vec<_> = order
        .iter()
        .take(FARTHEST)
        .map(|&i| json!({"point": points[i], "distance_mm": values[i]}))
        .collect();
    let mut found: Vec<f64> = values.iter().flatten().copied().collect();
    found.sort_by(f64::total_cmp);
    let quantile = |q: f64| {
        found
            .get(((found.len() as f64 - 1.0) * q).round() as usize)
            .copied()
    };
    json!({
        "points": values.len(),
        "beyond_search_cap": values.len() - found.len(),
        "search_cap_mm": SEARCH_CAP,
        "max_mm": found.last(),
        "median_mm": quantile(0.5),
        "p99_mm": quantile(0.99),
        "mean_mm": (!found.is_empty()).then(|| found.iter().sum::<f64>() / found.len() as f64),
        "beyond_one_pitch": found.iter().filter(|d| **d > pitch).count(),
        "beyond_half_pitch": found.iter().filter(|d| **d > 0.5 * pitch).count(),
        "farthest": farthest,
    })
}

fn mesh_summary(mesh: &SampledMesh) -> Value {
    json!({
        "status": SampledMesh::STATUS,
        "vertices": mesh.positions.len(),
        "triangles": mesh.triangles.len(),
        "open_edges": mesh.open_edges,
        "volume_mm3": mesh.volume(),
    })
}

fn chart_preview(plan: &Plan) -> (std::result::Result<PreviewMesh, Refusal>, f64, Value) {
    let start = Instant::now();
    let (answer, work) = with_operation_budget(PREVIEW_ALLOWANCE, || {
        preview_plan(plan, &PreviewOptions::default(), &mut Decider::default())
    });
    (
        answer,
        ms(start),
        json!({"used": work.used, "exhausted": work.exhausted, "allowance": PREVIEW_ALLOWANCE}),
    )
}

fn seeded_points(lo: [f64; 3], hi: [f64; 3]) -> Vec<[f64; 3]> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    (0..POINTS)
        .map(|_| {
            std::array::from_fn(|axis| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let u = (state >> 11) as f64 / (1u64 << 53) as f64;
                let pad = 0.1 * (hi[axis] - lo[axis]);
                lo[axis] - pad + u * (hi[axis] - lo[axis] + 2.0 * pad)
            })
        })
        .collect()
}

struct Run {
    workers: usize,
    max_points: usize,
    pitches: Vec<i32>,
    failures: Vec<String>,
}

fn part(name: &str, text: &str, gpu: Option<&Gpu>, run: &mut Run) -> Result<Value> {
    let plan = Plan::compile(&Recipe::parse(text)?)?;
    let options = |workers, dispatch| Options { workers, dispatch };
    let device = Evaluator::new(&plan, gpu, options(run.workers, MAX_DISPATCH))?;
    let split = Evaluator::new(&plan, gpu, options(run.workers, ODD_DISPATCH))?;
    let processor = Evaluator::new(&plan, None, options(run.workers, MAX_DISPATCH))?;
    let mut row = json!({
        "part": name,
        "status": "unproved",
        "route": hardware(device.hardware()),
        "pitches": [],
    });
    eprintln!("{name}: {:?}", device.hardware());

    let (chart, chart_ms, chart_work) = chart_preview(&plan);
    let chart_index = chart.as_ref().ok().map(|mesh| {
        let triangles = mesh.triangles.iter().map(|t| t.vertices).collect();
        Nearest::new(&mesh.positions, triangles, 0.5)
    });
    row["processor_chart_preview"] = match &chart {
        Ok(mesh) => json!({
            "status": "unproved",
            "ms": chart_ms,
            "work": chart_work,
            "vertices": mesh.positions.len(),
            "triangles": mesh.triangles.len(),
            "cells": mesh.cells,
            "unresolved_cells": mesh.unresolved_cells,
            "skipped_checks": mesh.skipped_checks,
        }),
        Err(refusal) => {
            json!({"refusal":{"code":refusal.code,"detail":refusal.detail},"ms":chart_ms,"work":chart_work})
        }
    };

    let mut checked_workers = false;
    for &exponent in &run.pitches.clone() {
        let grid: Grid = grid(&plan, exponent, &mut Decider::default())?;
        let pitch = grid.pitch();
        if let Err(refusal) = check_size(&grid, run.max_points) {
            eprintln!("{name} pitch {pitch} mm: {}", refusal.detail);
            row["pitches"].as_array_mut().unwrap().push(json!({
                "pitch_mm": pitch,
                "grid_counts": grid.counts,
                "grid_points": grid.len(),
                "refusal": {"code": refusal.code, "detail": refusal.detail},
            }));
            continue;
        }
        let start = Instant::now();
        let first = device.grid(&grid)?;
        let gpu_ms = ms(start);
        let start = Instant::now();
        let second = device.grid(&grid)?;
        let gpu_repeat_ms = ms(start);
        let split_values = split.grid(&grid)?;
        let repeat_identical = bits(&first) == bits(&second);
        let split_identical = bits(&first) == bits(&split_values);
        if !repeat_identical || !split_identical {
            run.failures.push(format!(
                "{name} pitch {pitch}: graphics card values differ between runs"
            ));
        }
        let start = Instant::now();
        let cpu = processor.grid(&grid)?;
        let cpu_ms = ms(start);
        let mut workers = json!({});
        // Worker-count identity on the first sampled pitch only.
        if !checked_workers {
            checked_workers = true;
            for w in [1, 2, 4] {
                let start = Instant::now();
                let values = Evaluator::new(&plan, None, options(w, MAX_DISPATCH))?.grid(&grid)?;
                let identical = bits(&values) == bits(&cpu);
                if !identical {
                    run.failures.push(format!(
                        "{name} pitch {pitch}: processor values differ at {w} workers"
                    ));
                }
                workers[w.to_string()] = json!({"ms": ms(start), "identical": identical});
            }
        }
        let start = Instant::now();
        let gpu_mesh = surface_nets(&grid, &first.values)?;
        let mesh_ms = ms(start);
        let cpu_mesh = surface_nets(&grid, &cpu.values)?;
        let same_triangles = gpu_mesh.triangles == cpu_mesh.triangles;
        let vertex_shift = (gpu_mesh.positions.len() == cpu_mesh.positions.len()).then(|| {
            gpu_mesh
                .positions
                .iter()
                .zip(&cpu_mesh.positions)
                .map(|(a, b)| dot(sub(*a, *b), sub(*a, *b)).sqrt())
                .fold(0.0, f64::max)
        });
        let mut pitch_row = json!({
            "pitch_mm": pitch,
            "grid_counts": grid.counts,
            "grid_points": grid.len(),
            "gpu_ms": gpu_ms,
            "gpu_repeat_ms": gpu_repeat_ms,
            "gpu_repeat_bit_identical": repeat_identical,
            "gpu_split_dispatch_bit_identical": split_identical,
            "split_dispatch_points": ODD_DISPATCH,
            "processor_ms": cpu_ms,
            "processor_workers": run.workers,
            "processor_worker_identity": workers,
            "surface_nets_ms": mesh_ms,
            "agreement": agreement(|i| grid.point(i), &first, &cpu),
            "gpu_mesh": mesh_summary(&gpu_mesh),
            "processor_grid_mesh": mesh_summary(&cpu_mesh),
            "meshes_same_triangles": same_triangles,
            "max_vertex_shift_mm": vertex_shift,
        });
        if let (Ok(chart), Some(index)) = (&chart, &chart_index) {
            let start = Instant::now();
            let triangles = gpu_mesh
                .triangles
                .iter()
                .map(|t| t.map(|v| v as usize))
                .collect();
            let gpu_index = Nearest::new(&gpu_mesh.positions, triangles, 0.5);
            pitch_row["gpu_mesh_to_chart_preview"] = statistics(
                &gpu_mesh.positions,
                &distances(&gpu_mesh.positions, index, run.workers),
                pitch,
            );
            pitch_row["chart_preview_to_gpu_mesh"] = statistics(
                &chart.positions,
                &distances(&chart.positions, &gpu_index, run.workers),
                pitch,
            );
            pitch_row["distance_ms"] = json!(ms(start));
        }
        eprintln!(
            "{name} pitch {pitch} mm: {} points, gpu {gpu_ms:.1} ms, processor {cpu_ms:.1} ms, side mismatches {}, triangles {}",
            grid.len(),
            pitch_row["agreement"]["side_mismatches"],
            gpu_mesh.triangles.len()
        );
        row["pitches"].as_array_mut().unwrap().push(pitch_row);
    }

    // Arbitrary points: the route the display's batch queries use.
    let bounds = catallax_field::printing::bounds(&plan, &mut Decider::default())?;
    let lo = std::array::from_fn(|i| bounds[i][0].to_f64_nearest());
    let hi = std::array::from_fn(|i| bounds[i][1].to_f64_nearest());
    let points = seeded_points(lo, hi);
    let start = Instant::now();
    let gpu_points = device.points(&points)?;
    let gpu_ms = ms(start);
    let repeat = device.points(&points)?;
    let split_points = split.points(&points)?;
    let identical = bits(&gpu_points) == bits(&repeat) && bits(&gpu_points) == bits(&split_points);
    if !identical {
        run.failures.push(format!(
            "{name} points: graphics card values differ between runs"
        ));
    }
    // The graphics card sees each point rounded to binary32; compare there.
    let rounded: Vec<[f64; 3]> = points
        .iter()
        .map(|p| p.map(|c| f64::from(c as f32)))
        .collect();
    let start = Instant::now();
    let cpu_points = processor.points(&rounded)?;
    let cpu_ms = ms(start);
    row["points"] = json!({
        "count": points.len(),
        "box": "material bounds grown by a tenth per side, seeded uniform",
        "gpu_ms": gpu_ms,
        "processor_ms": cpu_ms,
        "gpu_bit_identical_across_runs_and_dispatch_sizes": identical,
        "agreement_at_binary32_points": agreement(|i| rounded[i], &gpu_points, &cpu_points),
    });
    Ok(row)
}

fn run() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let usage = "usage: catallax-first-paint OUTPUT_DIRECTORY [--cpu] [--workers N] [--pitch EXPONENT]... [--part P1|P2|P3|P4]... [--max-points N]";
    let out = Path::new(args.first().ok_or(usage)?);
    let mut workers = std::thread::available_parallelism().map_or(1, usize::from);
    let mut pitches = Vec::new();
    let mut parts = Vec::new();
    let mut cpu_only = false;
    let mut max_points = MAX_POINTS;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--cpu" => cpu_only = true,
            "--workers" => {
                i += 1;
                workers = args.get(i).ok_or(usage)?.parse()?;
            }
            "--pitch" => {
                i += 1;
                pitches.push(args.get(i).ok_or(usage)?.parse::<i32>()?);
            }
            "--max-points" => {
                i += 1;
                max_points = args.get(i).ok_or(usage)?.parse()?;
            }
            "--part" => {
                i += 1;
                parts.push(args.get(i).ok_or(usage)?.clone());
            }
            _ => return Err(usage.into()),
        }
        i += 1;
    }
    if pitches.is_empty() {
        pitches = vec![-1, -2, -3];
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()?;
    fs::create_dir_all(out)?;
    if out.canonicalize()?.starts_with(&root) {
        return Err(Refusal::new(
            "preview.artifacts",
            "review artifacts must be outside the repository",
        )
        .into());
    }
    let gpu = Gpu::new(cpu_only);
    let mut report = json!({
        "status": "unproved preview",
        "export": false,
        "proven_measurement": false,
        "timing_context": "shared machine, information only",
        "device": match &gpu { Ok(g) => json!(g.name()), Err(e) => json!(format!("none: {e}")) },
        "workers": workers,
        "max_points": max_points,
        "parts": [],
    });
    let mut state = Run {
        workers,
        max_points,
        pitches,
        failures: Vec::new(),
    };
    for (name, text) in PARTS {
        if !parts.is_empty() && !parts.iter().any(|p| p == name) {
            continue;
        }
        let row = part(name, text, gpu.as_ref().ok(), &mut state)?;
        report["parts"].as_array_mut().unwrap().push(row);
        fs::write(out.join("report.json"), serde_json::to_vec_pretty(&report)?)?;
    }
    report["determinism_failures"] = json!(state.failures);
    fs::write(out.join("report.json"), serde_json::to_vec_pretty(&report)?)?;
    if !state.failures.is_empty() {
        return Err(Refusal::new("preview.nondeterministic", state.failures.join("; ")).into());
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
