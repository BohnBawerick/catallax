//! PLAN M23's round trip: 30,000 noisy points taken from P3 come back as P3's
//! planes and cylinders, as a clean part in both kernels, with a deviation map.
//! P3's concave junction round is a torus, which classic plane, cylinder and
//! cone fitting cannot explain: its samples are a counted refusal, and the part
//! keeps a sharp junction there, which both deviation numbers show.

#[path = "../examples/scan_fit/sampler.rs"]
mod sampler;

use catallax_field::reconstruction::scan::{
    ScanFit, ScanFitOptions, Surface, SurfaceKind, SurfaceUse, fit_scan,
};
use catallax_import::scan::{FastSurface, ScanItem};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{FileRef, Number, Operation, Pose, Recipe};

const COUNT: usize = 30_000;
const SIGMA: f64 = 0.01;
const SEED: u64 = 3;

fn p3() -> Recipe {
    Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap()
}

fn fit(bytes: &[u8], workers: usize) -> ScanFit {
    let step = Operation::ScanImport {
        file: FileRef {
            name: "p3-scan.obj".into(),
            sha256: catallax_import::sha256_hex(bytes),
        },
        pose: Pose::identity(),
        scale: Number::new(Rational::one()),
    };
    let options = ScanFitOptions {
        workers,
        ..ScanFitOptions::default()
    };
    fit_scan(&step, bytes, &options, &mut Decider::default()).expect("fit")
}

#[test]
fn the_sampler_lies_on_p3s_own_boundary() {
    let built = catallax_features::evaluate(&p3(), &mut Decider::default()).unwrap();
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &catallax_mesh::MeshOptions {
            chord_tolerance: Rational::from_ratio(1, 10_000).unwrap(),
            ..catallax_mesh::MeshOptions::default()
        },
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let surface = FastSurface::new(
        mesh.positions.clone(),
        mesh.triangles
            .iter()
            .map(|t| ScanItem::Triangle(t.vertices.map(|v| v as u32)))
            .collect(),
    );
    let worst = sampler::sample(&sampler::boss_block(true), 5_000, 0.0, 11)
        .iter()
        .map(|(p, _)| surface.distance(*p))
        .fold(0.0, f64::max);
    assert!(
        worst <= 1e-4 + 1e-6,
        "noise-free samples are {worst} mm off P3"
    );
}

#[test]
fn thirty_thousand_noisy_points_from_p3_come_back_as_its_planes_and_cylinders() {
    let patches = sampler::boss_block(true);
    let samples = sampler::sample(&patches, COUNT, SIGMA, SEED);
    let junction = patches.iter().position(|p| p.name == "junction").unwrap();
    let bytes = sampler::obj(&samples).into_bytes();
    let fit = fit(&bytes, 4);
    let roles = |role| {
        fit.surfaces
            .iter()
            .filter(move |s| s.role == role)
            .map(|s| s.surface.clone())
            .collect::<Vec<_>>()
    };
    // P3's seven planes: six block faces and the boss top.
    let faces = roles(SurfaceUse::StockFace);
    assert_eq!(faces.len(), 6);
    let mut offsets: Vec<f64> = faces
        .iter()
        .map(|s| match *s {
            Surface::Plane { normal, offset } => {
                let k = (0..3).max_by(|&a, &b| normal[a].abs().total_cmp(&normal[b].abs()));
                assert!(normal[k.unwrap()].abs() >= 1.0 - 1e-6, "{normal:?}");
                offset
            }
            _ => panic!("a stock face is a plane"),
        })
        .collect();
    offsets.sort_by(f64::total_cmp);
    for (got, want) in offsets.iter().zip([0.0, 0.0, 0.0, 10.0, 30.0, 40.0]) {
        assert!((got - want).abs() <= 0.005, "{offsets:?}");
    }
    let [Surface::Plane { offset: cap, .. }] = roles(SurfaceUse::Cap)[..] else {
        panic!("one cap plane");
    };
    assert!((cap - 22.0).abs() <= 0.005, "{cap}");
    // P3's five cylinders: the boss and the four top-edge rounds.
    let [
        Surface::Cylinder {
            point,
            axis,
            radius,
        },
    ] = roles(SurfaceUse::Boss)[..]
    else {
        panic!("one cylindrical boss");
    };
    assert!(axis[2] >= 1.0 - 1e-6, "{axis:?}");
    assert!((point[0] - 20.0).abs() <= 0.005 && (point[1] - 15.0).abs() <= 0.005);
    assert!((radius - 8.0).abs() <= 0.005, "{radius}");
    let rounds = roles(SurfaceUse::Fillet);
    assert_eq!(rounds.len(), 4);
    for round in &rounds {
        let Surface::Cylinder { radius, .. } = round else {
            panic!("a round is a cylinder");
        };
        assert!((radius - 1.0).abs() <= 0.02, "{radius}");
    }
    assert!(fit.surfaces.iter().all(|s| s.role != SurfaceUse::Unplaced));
    assert!(
        fit.surfaces
            .iter()
            .all(|s| s.surface.kind() != SurfaceKind::Cone)
    );
    // The torus is refused with its count, never dropped or forced.
    let unexplained: Vec<u32> = fit
        .unexplained
        .iter()
        .flat_map(|r| r.samples.iter().copied())
        .collect();
    let refusal = fit
        .refusals
        .iter()
        .find(|r| r.code == "reconstruction.unexplained-region")
        .expect("the junction round is refused");
    assert!(
        refusal
            .detail
            .starts_with(&format!("{} of {COUNT}", unexplained.len()))
    );
    let on_torus = samples.iter().filter(|s| s.1 == junction).count();
    let refused_on_torus = unexplained
        .iter()
        .filter(|&&i| samples[i as usize].1 == junction)
        .count();
    assert!(refused_on_torus * 100 >= unexplained.len() * 99);
    assert!(
        refused_on_torus * 2 >= on_torus,
        "{refused_on_torus} of {on_torus}"
    );
    let explained: usize = fit.surfaces.iter().map(|s| s.samples).sum();
    assert_eq!(explained + unexplained.len(), COUNT);
    // The clean part, in both kernels, and P3's dimensions.
    let part = fit.part.as_ref().expect("a part");
    assert!(part.brep.model.check(part.brep.result).is_empty());
    let steps: Vec<&Operation> = part.recipe.steps.iter().map(|s| &s.op).collect();
    let Operation::Block { size, .. } = steps[1] else {
        panic!("a block stock");
    };
    for (got, want) in size.iter().zip([40.0, 30.0, 10.0]) {
        assert!(
            (got.value().to_f64_nearest() - want).abs() <= 0.005,
            "{size:?}"
        );
    }
    let Operation::Cylinder { radius, height, .. } = steps[2] else {
        panic!("a cylindrical boss");
    };
    assert!((radius.value().to_f64_nearest() - 8.0).abs() <= 0.005);
    assert!((height.value().to_f64_nearest() - 12.0).abs() <= 0.005);
    let Operation::Fillet { radius, .. } = steps[4] else {
        panic!("one fillet step for the four top edges");
    };
    assert!((radius.value().to_f64_nearest() - 1.0).abs() <= 0.02);
    // Both deviation numbers, and the map that locates the difference.
    let d = &part.deviation;
    assert_eq!(d.scan_map.len(), COUNT);
    let (mut torus, mut elsewhere) = (0.0f64, 0.0f64);
    for (i, &v) in d.scan_map.iter().enumerate() {
        if samples[i].1 == junction {
            torus = torus.max(v);
        } else {
            elsewhere = elsewhere.max(v);
        }
    }
    // Away from the torus the part is within noise; on it the sharp
    // junction misses the round by up to 2 (1 - 1/sqrt 2) = 0.586 mm.
    assert!(elsewhere <= 6.0 * SIGMA, "{elsewhere}");
    assert!((torus - 0.586).abs() <= 6.0 * SIGMA, "{torus}");
    assert_eq!(d.scan_to_part.max_mm, torus);
    assert!(
        d.part_to_scan.max_mm >= 0.8 && d.part_to_scan.max_mm <= 1.2,
        "{d:?}"
    );
    let corner = d
        .part_map
        .iter()
        .filter(|(p, _, _)| {
            ((p[0] - 20.0).hypot(p[1] - 15.0) - 8.0).abs() < 0.05 && (p[2] - 10.0).abs() < 0.05
        })
        .map(|s| s.1)
        .fold(0.0, f64::max);
    // The sharp junction lies 2 (sqrt 2 - 1) = 0.828 mm from the round.
    assert!(corner >= 0.75, "{corner}");
}

#[test]
fn the_p3_round_trip_is_identical_at_one_two_and_four_workers() {
    let bytes = sampler::obj(&sampler::sample(
        &sampler::boss_block(true),
        COUNT,
        SIGMA,
        SEED,
    ))
    .into_bytes();
    let runs: Vec<ScanFit> = [1, 2, 4].into_iter().map(|w| fit(&bytes, w)).collect();
    for run in &runs[1..] {
        assert_eq!(run.assignment, runs[0].assignment);
        let (a, b) = (run.part.as_ref().unwrap(), runs[0].part.as_ref().unwrap());
        assert_eq!(a.recipe.to_text(), b.recipe.to_text());
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.deviation.scan_map), bits(&b.deviation.scan_map));
        let part = |m: &[([f64; 3], f64, usize)]| {
            m.iter().map(|s| (s.1.to_bits(), s.2)).collect::<Vec<_>>()
        };
        assert_eq!(part(&a.deviation.part_map), part(&b.deviation.part_map));
    }
}
