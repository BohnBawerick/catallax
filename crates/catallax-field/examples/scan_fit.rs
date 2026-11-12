//! Fits seeded noisy point clouds sampled from analytic parts and prints the
//! candidate surfaces, refusals, recipe and both deviation numbers.
//!
//! `cargo run --release -p catallax-field --example scan_fit -- <part> [count] [sigma] [seed] [--map <dir>]`
//! where `<part>` is `bore`, `cone`, `boss` or `p3`. With `--map`, the
//! deviation map is written to `<dir>` as coloured PLY points, a PNG with
//! three views per direction and a text summary, with the scan itself as `scan.obj`.

#[path = "scan_fit/map.rs"]
mod map;
#[path = "scan_fit/sampler.rs"]
mod sampler;

use catallax_field::reconstruction::scan::{ScanFitOptions, fit_scan};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{FileRef, Number, Operation, Pose};

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let map_dir = args.iter().position(|a| a == "--map").map(|i| {
        let dir = std::path::PathBuf::from(args.remove(i + 1));
        args.remove(i);
        dir
    });
    let part = args.first().map_or("boss", String::as_str);
    let count: usize = args.get(1).map_or(30_000, |s| s.parse().unwrap());
    let sigma: f64 = args.get(2).map_or(0.01, |s| s.parse().unwrap());
    let seed: u64 = args.get(3).map_or(1, |s| s.parse().unwrap());
    let patches = match part {
        "bore" => sampler::plate_with_bore(false),
        "cone" => sampler::plate_with_bore(true),
        "boss" => sampler::boss_block(false),
        "p3" => sampler::boss_block(true),
        other => panic!("unknown part {other}"),
    };
    let samples = sampler::sample(&patches, count, sigma, seed);
    let bytes = sampler::obj(&samples).into_bytes();
    let mut per_patch = std::collections::BTreeMap::new();
    for (_, k) in &samples {
        *per_patch.entry(patches[*k].name).or_insert(0usize) += 1;
    }
    println!("samples per patch {per_patch:?}");
    let step = Operation::ScanImport {
        file: FileRef {
            name: format!("{part}.obj"),
            sha256: catallax_import::sha256_hex(&bytes),
        },
        pose: Pose::identity(),
        scale: Number::new(Rational::one()),
    };
    let start = std::time::Instant::now();
    let fit = fit_scan(
        &step,
        &bytes,
        &ScanFitOptions::default(),
        &mut Decider::default(),
    )
    .expect("fit");
    println!(
        "fit in {:.2?}; sigma {:.5} mm",
        start.elapsed(),
        fit.sigma_mm
    );
    for (i, s) in fit.surfaces.iter().enumerate() {
        let mut truth = std::collections::BTreeMap::new();
        for (k, a) in fit.assignment.iter().enumerate() {
            if *a == Some(i as u32) {
                *truth.entry(patches[samples[k].1].name).or_insert(0usize) += 1;
            }
        }
        println!(
            "surface {i}: {:?} n={} rms={:.4} max={:.4} {:?}\n    truth {:?}",
            s.role, s.samples, s.rms_mm, s.max_mm, s.surface, truth
        );
    }
    for r in &fit.refusals {
        println!("refusal {}: {}", r.code, r.detail);
    }
    for (i, region) in fit
        .unexplained
        .iter()
        .enumerate()
        .filter(|(_, r)| r.samples.len() > 5)
    {
        let mut truth = std::collections::BTreeMap::new();
        for &k in &region.samples {
            *truth
                .entry(patches[samples[k as usize].1].name)
                .or_insert(0usize) += 1;
        }
        println!(
            "unexplained {i}: {} at {:?} truth {:?}",
            region.samples.len(),
            region.centroid,
            truth
        );
    }
    match &fit.part {
        Ok(p) => {
            println!("{}", p.recipe.to_text());
            let d = &p.deviation;
            println!(
                "scan->part n={} max={:.4} rms={:.4} p95={:.4}; part->scan n={} max={:.4} rms={:.4} p95={:.4}; chord {:.5}",
                d.scan_to_part.samples,
                d.scan_to_part.max_mm,
                d.scan_to_part.rms_mm,
                d.scan_to_part.p95_mm,
                d.part_to_scan.samples,
                d.part_to_scan.max_mm,
                d.part_to_scan.rms_mm,
                d.part_to_scan.p95_mm,
                d.chord_mm
            );
        }
        Err(r) => {
            println!("part refused {}: {}", r.code, r.detail);
            if let Some(recipe) = &fit.recipe {
                println!("{}", recipe.to_text());
            }
        }
    }
    if let Ok(p) = &fit.part {
        let mut worst: Vec<_> = p.deviation.part_map.iter().collect();
        worst.sort_by(|a, b| b.1.total_cmp(&a.1));
        for w in worst.iter().take(5) {
            println!("worst part sample {:?} d={:.4} face {}", w.0, w.1, w.2);
        }
    }
    if let (Some(dir), Ok(part)) = (&map_dir, &fit.part) {
        map::write(dir, &fit, part).expect("write the map");
        std::fs::write(dir.join("scan.obj"), &bytes).expect("write the scan");
        println!("map written to {}", dir.display());
    }
    println!("total {:.2?}", start.elapsed());
}
