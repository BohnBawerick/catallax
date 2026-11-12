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
    let runs: Vec<_> = [1, 2, 4]
        .into_iter()
        .map(|workers| {
            let mut dec = Decider::default();
            let limit = dec.policy().max_work;
            let (result, work) = catallax_numbers::with_operation_budget(limit, || {
                catallax_field::workers::with_workers(workers, || {
                    fit_scan(
                        &step(&bytes),
                        &bytes,
                        &ScanFitOptions {
                            workers,
                            ..ScanFitOptions::default()
                        },
                        &mut dec,
                    )
                })
                .expect("a positive field worker limit")
            });
            assert!(!work.exhausted);
            (result.expect("fit"), dec.take_stats(), work.used)
        })
        .collect();
    let base = &runs[0].0;
    let surface_bits = |surface: &Surface| match surface {
        Surface::Plane { normal, offset } => normal
            .iter()
            .chain([offset])
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        Surface::Cylinder {
            point,
            axis,
            radius,
        } => point
            .iter()
            .chain(axis)
            .chain([radius])
            .map(|v| v.to_bits())
            .collect(),
        Surface::Cone {
            point,
            axis,
            radius,
            slope,
        } => point
            .iter()
            .chain(axis)
            .chain([radius, slope])
            .map(|v| v.to_bits())
            .collect(),
    };
    for (run, stats, work) in &runs[1..] {
        assert_eq!(*work, runs[0].2);
        assert_eq!(stats, &runs[0].1);
        assert_eq!(run.sha256, base.sha256);
        assert_eq!(
            run.samples
                .iter()
                .map(|p| p.map(f64::to_bits))
                .collect::<Vec<_>>(),
            base.samples
                .iter()
                .map(|p| p.map(f64::to_bits))
                .collect::<Vec<_>>()
        );
        assert_eq!(run.assignment, base.assignment);
        assert_eq!(run.sigma_mm.to_bits(), base.sigma_mm.to_bits());
        let surfaces = |fit: &ScanFit| {
            fit.surfaces
                .iter()
                .map(|s| {
                    (
                        s.surface.kind(),
                        surface_bits(&s.surface),
                        s.samples,
                        s.rms_mm.to_bits(),
                        s.max_mm.to_bits(),
                        s.role,
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(surfaces(run), surfaces(base));
        let regions = |fit: &ScanFit| {
            fit.unexplained
                .iter()
                .map(|r| (r.samples.clone(), r.centroid.map(f64::to_bits)))
                .collect::<Vec<_>>()
        };
        assert_eq!(regions(run), regions(base));
        let refusals = |fit: &ScanFit| {
            fit.refusals
                .iter()
                .map(|r| (r.code.to_string(), r.detail.clone()))
                .collect::<Vec<_>>()
        };
        assert_eq!(refusals(run), refusals(base));
        let (a, b) = (run.part.as_ref().unwrap(), base.part.as_ref().unwrap());
        assert_eq!(a.recipe.to_text(), b.recipe.to_text());
        let side = |s: &catallax_field::reconstruction::scan::DeviationSide| {
            (
                s.samples,
                s.max_mm.to_bits(),
                s.rms_mm.to_bits(),
                s.p95_mm.to_bits(),
            )
        };
        let bits = |d: &catallax_field::reconstruction::scan::ScanDeviation| {
            (
                side(&d.scan_to_part),
                side(&d.part_to_scan),
                d.chord_mm.to_bits(),
                d.scan_map.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                d.part_map
                    .iter()
                    .map(|(p, d, f)| (p.map(f64::to_bits), d.to_bits(), *f))
                    .collect::<Vec<_>>(),
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

#[test]
fn six_near_orthogonal_plane_classes_refuse_without_a_panic() {
    let s = 0.9_f64.to_radians().sin();
    let z = (1.0 - 2.0 * s * s).sqrt();
    let mut patches = Vec::new();
    for axis in 0..3 {
        for high in [false, true] {
            patches.push(sampler::Patch {
                name: "near-orthogonal face",
                area: 400.0,
                sample: Box::new(move |u, v| {
                    let mut q = [0.0; 3];
                    q[axis] = if high { 20.0 } else { 0.0 };
                    q[(axis + 1) % 3] = 20.0 * u;
                    q[(axis + 2) % 3] = 20.0 * v;
                    let p = [q[0], q[1], (q[2] - s * (q[0] + q[1])) / z];
                    let mut normal = [0.0; 3];
                    normal[axis] = 1.0;
                    if axis == 2 {
                        normal = [s, s, z];
                    }
                    let sign = if high { 1.0 } else { -1.0 };
                    Some((p, normal.map(|n| sign * n)))
                }),
            });
        }
    }
    let bytes = sampler::obj(&sampler::sample(&patches, 12_000, 0.0, 19)).into_bytes();
    let fit = fit_scan(
        &step(&bytes),
        &bytes,
        &ScanFitOptions {
            workers: 1,
            noise_mm: Some(0.001),
            ..ScanFitOptions::default()
        },
        &mut Decider::default(),
    )
    .expect("read and fit the six planes");
    assert_eq!(fit.surfaces.len(), 6);
    assert_eq!(
        fit.part.err().expect("the stock frame does not fit").code,
        "reconstruction.no-stock"
    );
}

#[test]
fn a_nonuniform_partial_cylinder_is_unplaced() {
    let mut patches = sampler::plate_with_bore(false);
    patches.last_mut().unwrap().sample = Box::new(|u, v| {
        let angle = (2.0 * std::f64::consts::PI / 3.0) * u.powi(4);
        let (s, c) = angle.sin_cos();
        Some(([15.0 + 4.0 * c, 10.0 + 4.0 * s, 8.0 * v], [-c, -s, 0.0]))
    });
    let fit = fit(&cloud(&patches, 20_000, 23), 1);
    assert!(
        fit.surfaces
            .iter()
            .any(|s| s.surface.kind() == SurfaceKind::Cylinder),
        "the partial cylinder must reach assembly"
    );
    assert!(kinds(&fit, SurfaceUse::Bore).is_empty());
    assert!(fit.refusals.iter().any(|r| {
        r.code == "reconstruction.unplaced-surface" && r.detail.contains("more than half way round")
    }));
}

#[test]
fn collapsed_binary64_neighbours_exhaust_the_actual_work_allowance() {
    use catallax_numbers::with_operation_budget;

    let bytes: Vec<u8> = (1..=1024)
        .map(|i| format!("v 1.{i:030} 0 0\n"))
        .collect::<String>()
        .into_bytes();
    let scan_step = step(&bytes);
    let (scan, read_work) = with_operation_budget(1_000_000, || {
        catallax_import::scan::ScanReference::from_step(&scan_step, &bytes, &mut Decider::default())
    });
    let scan = scan.expect("the exact vertices are valid reference data");
    assert!(!read_work.exhausted);
    assert_eq!(scan.report().isolated_points, 1024);
    assert!(scan.approx_points().iter().all(|p| {
        p.map(f64::to_bits) == [1.0_f64.to_bits(), 0.0_f64.to_bits(), 0.0_f64.to_bits()]
    }));
    let limit = read_work.used + 300_000;
    let (result, work) = with_operation_budget(limit, || {
        fit_scan(
            &scan_step,
            &bytes,
            &ScanFitOptions {
                workers: 1,
                noise_mm: Some(0.01),
                ..ScanFitOptions::default()
            },
            &mut Decider::default(),
        )
    });
    assert!(
        work.exhausted,
        "nearest visits must exhaust the allowance; used {} of {limit}",
        work.used
    );
    assert_eq!(
        result
            .err()
            .expect("fitting must stop on exhausted work")
            .code,
        "numbers.operation-budget"
    );
}
