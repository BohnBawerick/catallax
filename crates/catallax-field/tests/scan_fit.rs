//! Classic fitting of noisy point clouds to a clean part (PLAN M23): region
//! growing, then plane, cylinder and cone fits, an ordinary recipe built in
//! both kernels, and both deviation numbers every time. The clouds are seeded
//! samples of analytic parts written independently of the fitter.

#[path = "../examples/scan_fit/sampler.rs"]
mod sampler;

use catallax_field::reconstruction::scan::{
    ScanFit, ScanFitOptions, Surface, SurfaceKind, SurfaceUse, fit_scan,
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{FileRef, Number, Operation, Pose};

const SIGMA: f64 = 0.01;

fn step(bytes: &[u8]) -> Operation {
    Operation::ScanImport {
        file: FileRef {
            name: "cloud.obj".into(),
            sha256: catallax_import::sha256_hex(bytes),
        },
        pose: Pose::identity(),
        scale: Number::new(Rational::one()),
    }
}

fn fit(bytes: &[u8], workers: usize) -> ScanFit {
    let options = ScanFitOptions {
        workers,
        ..ScanFitOptions::default()
    };
    fit_scan(&step(bytes), bytes, &options, &mut Decider::default()).expect("fit")
}

fn cloud(patches: &[sampler::Patch], count: usize, seed: u64) -> Vec<u8> {
    sampler::obj(&sampler::sample(patches, count, SIGMA, seed)).into_bytes()
}

fn kinds(fit: &ScanFit, role: SurfaceUse) -> Vec<SurfaceKind> {
    fit.surfaces
        .iter()
        .filter(|s| s.role == role)
        .map(|s| s.surface.kind())
        .collect()
}

fn number(n: &Number) -> f64 {
    n.value().to_f64_nearest()
}

/// The stock block's size and the first bore's operation.
fn stock_and_tool(fit: &ScanFit) -> ([f64; 3], Operation) {
    let recipe = &fit.part.as_ref().expect("a part").recipe;
    let Operation::Block { size, .. } = &recipe.steps[1].op else {
        panic!("the stock is a block");
    };
    (size.each_ref().map(number), recipe.steps[2].op.clone())
}

#[test]
fn a_noisy_plate_with_a_bore_comes_back_as_six_planes_and_a_cylinder() {
    let bytes = cloud(&sampler::plate_with_bore(false), 20_000, 1);
    let fit = fit(&bytes, 4);
    assert!(
        (fit.sigma_mm - SIGMA).abs() < 0.2 * SIGMA,
        "{}",
        fit.sigma_mm
    );
    assert_eq!(
        kinds(&fit, SurfaceUse::StockFace),
        vec![SurfaceKind::Plane; 6]
    );
    assert_eq!(kinds(&fit, SurfaceUse::Bore), vec![SurfaceKind::Cylinder]);
    assert!(kinds(&fit, SurfaceUse::Unplaced).is_empty());
    let unexplained: usize = fit.unexplained.iter().map(|r| r.samples.len()).sum();
    assert!(unexplained <= 20, "{unexplained} unexplained samples");
    let (size, tool) = stock_and_tool(&fit);
    for (got, want) in size.iter().zip([30.0, 20.0, 8.0]) {
        assert!((got - want).abs() <= 0.005, "{size:?}");
    }
    let Operation::Cylinder { radius, .. } = tool else {
        panic!("a cylindrical bore");
    };
    assert!((number(&radius) - 4.0).abs() <= 0.005);
    let part = fit.part.as_ref().unwrap();
    assert!(part.brep.model.check(part.brep.result).is_empty());
    let d = &part.deviation;
    assert_eq!(d.scan_to_part.samples, 20_000);
    assert!(d.scan_to_part.max_mm <= 6.0 * SIGMA + 0.01, "{d:?}");
    assert!(d.scan_to_part.rms_mm <= 1.2 * SIGMA, "{d:?}");
    assert!(
        d.part_to_scan.samples > 0 && d.part_to_scan.max_mm < 1.0,
        "{d:?}"
    );
    assert!(d.chord_mm <= 0.001);
}

#[test]
fn a_conical_bore_is_fitted_as_a_cone() {
    let bytes = cloud(&sampler::plate_with_bore(true), 20_000, 2);
    let fit = fit(&bytes, 4);
    assert_eq!(kinds(&fit, SurfaceUse::Bore), vec![SurfaceKind::Cone]);
    let bore = fit
        .surfaces
        .iter()
        .find(|s| s.role == SurfaceUse::Bore)
        .unwrap();
    let Surface::Cone { slope, .. } = bore.surface else {
        unreachable!()
    };
    // Radius 3 at the bottom to 5 at the top over 8 mm.
    assert!((slope.abs() - 0.25).abs() <= 0.005, "{slope}");
    let (_, tool) = stock_and_tool(&fit);
    let Operation::Cone {
        base_radius,
        top_radius,
        height,
        ..
    } = tool
    else {
        panic!("a conical bore");
    };
    // The tool runs 1 mm past both faces, so its radii follow the same law.
    assert!((number(&base_radius) - 2.75).abs() <= 0.01);
    assert!((number(&top_radius) - 5.25).abs() <= 0.01);
    assert!((number(&height) - 10.0).abs() <= 0.005);
    let d = &fit.part.as_ref().unwrap().deviation;
    assert!(d.scan_to_part.max_mm <= 6.0 * SIGMA + 0.01, "{d:?}");
}

#[test]
fn a_block_with_a_boss_and_rounded_top_edges_is_rebuilt_with_its_rounds() {
    let bytes = cloud(&sampler::boss_block(false), 30_000, 3);
    let fit = fit(&bytes, 4);
    assert_eq!(
        kinds(&fit, SurfaceUse::StockFace),
        vec![SurfaceKind::Plane; 6]
    );
    assert_eq!(kinds(&fit, SurfaceUse::Boss), vec![SurfaceKind::Cylinder]);
    assert_eq!(kinds(&fit, SurfaceUse::Cap), vec![SurfaceKind::Plane]);
    assert!(kinds(&fit, SurfaceUse::Fillet).len() >= 4);
    let part = fit.part.as_ref().expect("a part");
    let fillets: Vec<_> = part
        .recipe
        .steps
        .iter()
        .filter_map(|s| match &s.op {
            Operation::Fillet { radius, edges, .. } => Some((number(radius), edges.clone())),
            _ => None,
        })
        .collect();
    assert_eq!(fillets.len(), 1, "one round for the four edges");
    assert!((fillets[0].0 - 1.0).abs() <= 0.02, "{}", fillets[0].0);
    let d = &part.deviation;
    assert!(d.scan_to_part.max_mm <= 6.0 * SIGMA + 0.02, "{d:?}");
}

#[test]
fn a_hole_in_the_scan_shows_in_part_to_scan_only() {
    let patches = sampler::plate_with_bore(false);
    // Remove every sample within 3 mm of (6, 6) on the top face: an occlusion.
    let samples: Vec<_> = sampler::sample(&patches, 20_000, SIGMA, 4)
        .into_iter()
        .filter(|(p, _)| !((p[0] - 6.0).hypot(p[1] - 6.0) < 3.0 && p[2] > 7.9))
        .collect();
    let bytes = sampler::obj(&samples).into_bytes();
    let fit = fit(&bytes, 4);
    let d = &fit.part.as_ref().expect("a part").deviation;
    assert!(d.scan_to_part.max_mm <= 6.0 * SIGMA + 0.01, "{d:?}");
    assert!(d.part_to_scan.max_mm >= 2.5, "{d:?}");
}

#[test]
fn samples_no_plane_cylinder_or_cone_explains_are_counted_refusals() {
    let mut patches = sampler::plate_with_bore(false);
    patches.push(sampler::dome("dome", [6.0, 14.0], 3.0, 8.0));
    let samples = sampler::sample(&patches, 20_000, SIGMA, 5);
    let dome = samples.iter().filter(|s| s.1 == patches.len() - 1).count();
    let bytes = sampler::obj(&samples).into_bytes();
    let fit = fit(&bytes, 4);
    let refusal = fit
        .refusals
        .iter()
        .find(|r| r.code == "reconstruction.unexplained-region")
        .expect("the dome is refused, not dropped");
    let unexplained: usize = fit.unexplained.iter().map(|r| r.samples.len()).sum();
    assert!(
        refusal
            .detail
            .starts_with(&format!("{unexplained} of 20000"))
    );
    assert!(
        unexplained * 2 >= dome,
        "{unexplained} of {dome} dome samples"
    );
    let explained: usize = fit.surfaces.iter().map(|s| s.samples).sum();
    assert_eq!(
        explained + unexplained,
        20_000,
        "every sample is accounted for"
    );
}

#[test]
fn one_two_and_four_workers_give_identical_fits() {
    let bytes = cloud(&sampler::boss_block(false), 12_000, 6);
    let runs: Vec<ScanFit> = [1, 2, 4].into_iter().map(|w| fit(&bytes, w)).collect();
    for run in &runs[1..] {
        assert_eq!(run.assignment, runs[0].assignment);
        assert_eq!(run.sigma_mm.to_bits(), runs[0].sigma_mm.to_bits());
        let (a, b) = (run.part.as_ref().unwrap(), runs[0].part.as_ref().unwrap());
        assert_eq!(a.recipe.to_text(), b.recipe.to_text());
        let bits = |d: &catallax_field::reconstruction::scan::ScanDeviation| {
            (
                d.scan_to_part.max_mm.to_bits(),
                d.part_to_scan.max_mm.to_bits(),
                d.scan_map.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                d.part_map.len(),
            )
        };
        assert_eq!(bits(&a.deviation), bits(&b.deviation));
    }
}

#[test]
fn too_few_samples_and_a_wrong_hash_refuse() {
    let bytes = cloud(&sampler::plate_with_bore(false), 50, 7);
    let refusal = fit_scan(
        &step(&bytes),
        &bytes,
        &ScanFitOptions::default(),
        &mut Decider::default(),
    )
    .err()
    .unwrap();
    assert_eq!(refusal.code, "reconstruction.too-few-samples");
    let other = step(b"other");
    let refusal = fit_scan(
        &other,
        &bytes,
        &ScanFitOptions::default(),
        &mut Decider::default(),
    )
    .err()
    .unwrap();
    assert_eq!(refusal.code, "import.scan-hash");
}
