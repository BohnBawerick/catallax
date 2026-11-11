//! Classic fitting of a scan to a clean part (PLAN M23): region growing, then
//! plane, cylinder and cone fits, then an ordinary recipe built in both
//! kernels, with both deviation numbers every time.
//!
//! Every fitted surface is a candidate. The part is defined by its recipe and
//! is exact under the existing B-rep and field definitions; nothing here claims
//! it equals the scanned object. Both deviations are unproved binary64 samples
//! over stated domains. Samples no plane, cylinder or cone explains are typed
//! refusals with counts, and fitted surfaces the block family cannot place
//! are refusals too, never silent drops. Free-form regions wait for M10 and
//! M11.

mod assemble;
mod fit;
mod knn;
mod linalg;
mod segment;

use catallax_import::scan::{FastSurface, ScanItem, ScanReference};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Operation, Recipe};

use crate::Plan;
pub use fit::{Surface, SurfaceKind};
use linalg::V3;

/// Most samples a scan may give the fitter. Above it the fit refuses with
/// `reconstruction.too-many-samples`.
pub const MAX_SAMPLES: usize = 1_000_000;
/// Fewest samples the fitter accepts.
pub const MIN_SAMPLES: usize = 100;
/// The part's display mesh for the deviation is certified within this many
/// millimetres of the B-rep.
pub const CHORD_MM: f64 = 0.001;

/// Fitting and sampling controls. None of them is a geometric tolerance of
/// the part; they choose how the scan is read and how densely deviation is sampled.
#[derive(Clone, Debug)]
pub struct ScanFitOptions {
    /// Threads for per-sample work. Answers are identical at any count.
    pub workers: usize,
    /// Noise deviation normal to the surface, in millimetres; estimated from
    /// the samples' local planes when absent.
    pub noise_mm: Option<f64>,
    /// Spacing of part-to-scan samples on the part, in millimetres.
    pub part_spacing_mm: f64,
    /// Barycentric subdivisions of each scan triangle for scan-to-part samples.
    pub triangle_steps: usize,
}

impl Default for ScanFitOptions {
    fn default() -> Self {
        Self {
            workers: std::thread::available_parallelism().map_or(1, usize::from),
            noise_mm: None,
            part_spacing_mm: 0.25,
            triangle_steps: 4,
        }
    }
}

/// How the part uses a fitted surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceUse {
    /// An outer face of the block stock.
    StockFace,
    /// The flat end of a boss or of a blind bore.
    Cap,
    /// A cylinder or cone added to the stock.
    Boss,
    /// A cylinder or cone removed from the stock.
    Bore,
    /// A convex round on a stock edge.
    Fillet,
    /// Not placed by the block family; reported as a refusal.
    Unplaced,
}

/// One fitted candidate surface.
#[derive(Clone, Debug)]
pub struct FittedSurface {
    /// The least-squares surface.
    pub surface: Surface,
    /// Samples it explains.
    pub samples: usize,
    /// Root-mean-square residual of those samples, in millimetres.
    pub rms_mm: f64,
    /// Largest absolute residual of those samples, in millimetres.
    pub max_mm: f64,
    /// How the part uses it.
    pub role: SurfaceUse,
}

pub(crate) struct Unplaced {
    surface: usize,
    samples: usize,
    reason: &'static str,
}

/// A connected group of samples no fitted surface explains.
#[derive(Clone, Debug)]
pub struct UnexplainedRegion {
    /// Its samples, by index into [`ScanFit::samples`].
    pub samples: Vec<u32>,
    /// Their mean, in millimetres.
    pub centroid: [f64; 3],
}

/// One direction of the sampled deviation. Unproved: binary64 distances at
/// finite samples, not a Hausdorff bound.
#[derive(Clone, Debug, PartialEq)]
pub struct DeviationSide {
    /// Samples measured.
    pub samples: usize,
    /// Largest sampled distance, in millimetres.
    pub max_mm: f64,
    /// Root-mean-square sampled distance, in millimetres.
    pub rms_mm: f64,
    /// The 95th percentile of sampled distances, in millimetres.
    pub p95_mm: f64,
}

/// Both deviation numbers and the maps they come from.
#[derive(Clone, Debug)]
pub struct ScanDeviation {
    /// From every scan sample to the part's boundary.
    pub scan_to_part: DeviationSide,
    /// From samples on the part's boundary to the scan.
    pub part_to_scan: DeviationSide,
    /// Proven bound on the display mesh's distance from the B-rep, in
    /// millimetres; distances to the part are measured on that mesh.
    pub chord_mm: f64,
    /// Scan-to-part distance of each fitting sample, in [`ScanFit::samples`] order.
    pub scan_map: Vec<f64>,
    /// Each part sample, its distance to the scan and its B-rep face index.
    pub part_map: Vec<([f64; 3], f64, usize)>,
}

/// The clean part: an ordinary recipe and its two kernels' values.
pub struct FittedPart {
    /// The recipe: the scan step, then the reconstruction's ordinary steps.
    pub recipe: Recipe,
    /// The validated B-rep.
    pub brep: catallax_features::Built,
    /// The field of the same recipe.
    pub field: Plan,
    /// Both deviation numbers, always present.
    pub deviation: ScanDeviation,
}

/// The fit of one scan.
pub struct ScanFit {
    /// The scan's SHA-256.
    pub sha256: String,
    /// The samples fitted: the scan's points in its order.
    pub samples: Vec<[f64; 3]>,
    /// Noise deviation used, in millimetres.
    pub sigma_mm: f64,
    /// Candidate surfaces, in a deterministic order.
    pub surfaces: Vec<FittedSurface>,
    /// Per sample: index into `surfaces`, or `None` when unexplained.
    pub assignment: Vec<Option<u32>>,
    /// Connected groups of unexplained samples.
    pub unexplained: Vec<UnexplainedRegion>,
    /// Typed refusals: unexplained regions and unplaced surfaces, with counts.
    pub refusals: Vec<Refusal>,
    /// The assembled recipe, also when a kernel then refuses to build it.
    pub recipe: Option<Recipe>,
    /// The part, or why the block family could not make one.
    pub part: Result<FittedPart, Refusal>,
}

/// Maps `f` over `0..n` on up to `workers` threads, in index order.
pub(crate) fn par_map<T: Send>(workers: usize, n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let workers = workers.clamp(1, n.max(1));
    if workers == 1 {
        return (0..n).map(f).collect();
    }
    let chunk = n.div_ceil(workers);
    std::thread::scope(|scope| {
        let f = &f;
        let handles: Vec<_> = (0..workers)
            .map(|w| {
                scope.spawn(move || {
                    (w * chunk..((w + 1) * chunk).min(n))
                        .map(f)
                        .collect::<Vec<T>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("fitting worker"))
            .collect()
    })
}

fn stats(mut values: Vec<f64>) -> DeviationSide {
    let samples = values.len();
    let rms = (values.iter().map(|v| v * v).sum::<f64>() / samples.max(1) as f64).sqrt();
    values.sort_by(f64::total_cmp);
    DeviationSide {
        samples,
        max_mm: values.last().copied().unwrap_or(0.0),
        rms_mm: rms,
        p95_mm: values
            .get((samples * 95).div_ceil(100).saturating_sub(1))
            .copied()
            .unwrap_or(0.0),
    }
}

/// Deterministic samples on a triangle: its corners, its edges at `spacing`
/// and its interior at about one per `spacing` squared, on an additive
/// recurrence with irrational steps.
fn triangle_samples(t: [V3; 3], spacing: f64, out: &mut Vec<V3>) {
    let at = |a: f64, b: f64| -> V3 {
        std::array::from_fn(|k| t[0][k] + a * (t[1][k] - t[0][k]) + b * (t[2][k] - t[0][k]))
    };
    for e in 0..3 {
        let (p, q) = (t[e], t[(e + 1) % 3]);
        let n = (linalg::norm(linalg::sub(q, p)) / spacing).ceil().max(1.0) as usize;
        for i in 0..n {
            let s = i as f64 / n as f64;
            out.push(std::array::from_fn(|k| p[k] + s * (q[k] - p[k])));
        }
    }
    let area = linalg::norm(linalg::cross(
        linalg::sub(t[1], t[0]),
        linalg::sub(t[2], t[0]),
    )) / 2.0;
    let m = (area / (spacing * spacing)).floor() as usize;
    // The plastic-number recurrence (R2) folded into the triangle.
    const A1: f64 = 0.754_877_666_246_692_7;
    const A2: f64 = 0.569_840_290_998_053_3;
    for i in 1..=m {
        let (mut a, mut b) = ((i as f64 * A1).fract(), (i as f64 * A2).fract());
        if a + b > 1.0 {
            (a, b) = (1.0 - a, 1.0 - b);
        }
        out.push(at(a, b));
    }
}

fn deviation(
    scan: &ScanReference,
    samples: &[V3],
    brep: &catallax_features::Built,
    options: &ScanFitOptions,
    dec: &mut Decider,
) -> Result<ScanDeviation, Refusal> {
    let mesh = catallax_mesh::mesh_body(
        &brep.model,
        brep.result,
        &catallax_mesh::MeshOptions {
            chord_tolerance: Rational::from_ratio(1, 1000).expect("nonzero"),
            ..catallax_mesh::MeshOptions::default()
        },
        dec,
    )?
    .value;
    let chord_mm = mesh
        .faces
        .iter()
        .map(|f| f.max_chord_error.to_f64_nearest())
        .fold(0.0, f64::max);
    let faces: Vec<usize> = mesh
        .triangles
        .iter()
        .map(|t| t.face.index() as usize)
        .collect();
    let part = FastSurface::new(
        mesh.positions.clone(),
        mesh.triangles
            .iter()
            .map(|t| ScanItem::Triangle(t.vertices.map(|v| v as u32)))
            .collect(),
    );
    let workers = options.workers;
    // Scan to part: every sample, plus interior grid points of scan triangles.
    let scan_map = par_map(workers, samples.len(), |i| part.distance(samples[i]));
    let steps = options.triangle_steps.clamp(1, 16);
    let mut extra = Vec::new();
    let points = scan.approx_points();
    for item in scan.items() {
        if let ScanItem::Triangle(t) = item {
            let p = t.map(|v| points[v as usize]);
            for i in 1..steps {
                for j in 1..steps - i {
                    let (a, b) = (i as f64 / steps as f64, j as f64 / steps as f64);
                    extra.push(std::array::from_fn(|k| {
                        p[0][k] + a * (p[1][k] - p[0][k]) + b * (p[2][k] - p[0][k])
                    }));
                }
            }
        }
    }
    let extra_map = par_map(workers, extra.len(), |i| part.distance(extra[i]));
    dec.charge_work((samples.len() + extra.len()) as u64)?;
    let mut part_samples = Vec::new();
    let mut part_faces = Vec::new();
    let spacing = options.part_spacing_mm.max(1e-3);
    for (t, face) in mesh.triangles.iter().zip(&faces) {
        let before = part_samples.len();
        triangle_samples(
            t.vertices.map(|v| mesh.positions[v]),
            spacing,
            &mut part_samples,
        );
        part_faces.resize(part_samples.len(), *face);
        dec.charge_work((part_samples.len() - before) as u64)?;
    }
    let part_distances = par_map(workers, part_samples.len(), |i| {
        scan.surface().distance(part_samples[i])
    });
    if scan_map
        .iter()
        .chain(&extra_map)
        .chain(&part_distances)
        .any(|d| !d.is_finite())
    {
        return Err(Refusal::new(
            "reconstruction.nonfinite-deviation",
            "a sampled deviation is not a finite number",
        ));
    }
    let scan_to_part = stats(scan_map.iter().chain(&extra_map).copied().collect());
    let part_to_scan = stats(part_distances.clone());
    Ok(ScanDeviation {
        scan_to_part,
        part_to_scan,
        chord_mm,
        scan_map,
        part_map: part_samples
            .into_iter()
            .zip(part_distances)
            .zip(part_faces)
            .map(|((p, d), f)| (p, d, f))
            .collect(),
    })
}

/// Fits the scan named by `scan_step`, a `scan_import` step, whose file bytes
/// are `bytes`. The result always carries the candidate surfaces and every
/// refusal; the part, when the block family can make one, carries both
/// deviation numbers.
pub fn fit_scan(
    scan_step: &Operation,
    bytes: &[u8],
    options: &ScanFitOptions,
    dec: &mut Decider,
) -> Result<ScanFit, Refusal> {
    let scan = ScanReference::from_step(scan_step, bytes, dec)?;
    let samples: Vec<V3> = scan.approx_points().to_vec();
    if samples.len() < MIN_SAMPLES {
        return Err(Refusal::new(
            "reconstruction.too-few-samples",
            format!(
                "the scan has {} points; fitting needs at least {MIN_SAMPLES}",
                samples.len()
            ),
        ));
    }
    if samples.len() > MAX_SAMPLES {
        return Err(Refusal::new(
            "reconstruction.too-many-samples",
            format!(
                "the scan has {} points; the stated limit is {MAX_SAMPLES}",
                samples.len()
            ),
        ));
    }
    let workers = options.workers.max(1);
    let cloud = segment::Cloud::new(&samples, options.noise_mm, workers);
    dec.charge_work((samples.len() * segment::K_NORMAL * 8) as u64)?;
    let found = segment::segment(&cloud, workers);
    dec.charge_work(found.work)?;
    // Fillets refitted tangent to their stock faces, and small surfaces the
    // block family cannot place, change what explains each sample, so the
    // samples are assigned again before the final assembly.
    let mut found = found;
    if let Ok(first) = assemble::assemble(scan_step, &found.surfaces, &samples, cloud.sigma) {
        let small = (samples.len() / 100).max(segment::MIN_SURFACE);
        let mut surfaces = found.surfaces.clone();
        for (s, surface) in first.refined {
            surfaces[s].0 = surface;
        }
        let mut s = 0;
        surfaces.retain(|(_, members)| {
            let keep = first.uses[s] != SurfaceUse::Unplaced || members.len() >= small;
            s += 1;
            keep
        });
        let mut work = 0;
        found = segment::assign(&cloud, &mut surfaces, workers, &mut work);
        dec.charge_work(work)?;
    }
    let assembly = assemble::assemble(scan_step, &found.surfaces, &samples, cloud.sigma);
    let unexplained: Vec<UnexplainedRegion> = segment::unexplained(&cloud, &found.assignment)
        .into_iter()
        .map(|group| {
            let (centroid, _, _) = linalg::moments(group.iter().map(|&i| samples[i as usize]));
            UnexplainedRegion {
                samples: group,
                centroid,
            }
        })
        .collect();
    let uses = match &assembly {
        Ok(a) => a.uses.clone(),
        Err(_) => vec![SurfaceUse::Unplaced; found.surfaces.len()],
    };
    // A round is reported as the part uses it: tangent to its stock faces.
    let mut reported = found.surfaces.clone();
    if let Ok(a) = &assembly {
        for (s, surface) in &a.refined {
            reported[*s].0 = surface.clone();
        }
    }
    let surfaces = reported
        .iter()
        .zip(&uses)
        .map(|((surface, members), &role)| {
            let residuals: Vec<f64> = members
                .iter()
                .map(|&i| surface.residual(samples[i as usize]).abs())
                .collect();
            FittedSurface {
                surface: surface.clone(),
                samples: members.len(),
                rms_mm: fit::rms(surface, &samples, members),
                max_mm: residuals.iter().copied().fold(0.0, f64::max),
                role,
            }
        })
        .collect();
    let mut refusals = Vec::new();
    if !unexplained.is_empty() {
        let total: usize = unexplained.iter().map(|r| r.samples.len()).sum();
        let largest = unexplained
            .iter()
            .map(|r| r.samples.len())
            .max()
            .unwrap_or(0);
        refusals.push(Refusal::new(
            "reconstruction.unexplained-region",
            format!(
                "{total} of {} samples in {} connected regions (largest {largest}) are explained by no plane, cylinder or cone",
                samples.len(),
                unexplained.len()
            ),
        ));
    }
    let recipe = assembly.as_ref().ok().map(|a| a.recipe.clone());
    let part = assembly.and_then(|assembly| {
        for u in &assembly.unplaced {
            refusals.push(Refusal::new(
                "reconstruction.unplaced-surface",
                format!(
                    "surface {} ({} samples): {}",
                    u.surface, u.samples, u.reason
                ),
            ));
        }
        let recipe = assembly.recipe;
        let brep = catallax_features::evaluate(&recipe, dec).map_err(|e| match e {
            catallax_features::BuildError::Geometry { refusal, .. } => refusal,
            catallax_features::BuildError::Recipe(e) => {
                Refusal::new("reconstruction.recipe", e.to_string())
            }
        })?;
        if !brep.model.check(brep.result).is_empty() {
            return Err(Refusal::new(
                "reconstruction.invalid-part",
                "the reconstructed B-rep failed structural validation",
            ));
        }
        let field = Plan::compile_checked(&recipe, dec)?;
        let deviation = deviation(&scan, &samples, &brep, options, dec)?;
        Ok(FittedPart {
            recipe,
            brep,
            field,
            deviation,
        })
    });
    let sigma_mm = cloud.sigma;
    drop(cloud);
    Ok(ScanFit {
        sha256: scan.sha256().to_string(),
        samples,
        sigma_mm,
        surfaces,
        assignment: found.assignment,
        recipe,
        unexplained,
        refusals,
        part,
    })
}
