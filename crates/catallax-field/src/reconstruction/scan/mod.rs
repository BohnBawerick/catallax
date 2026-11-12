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
use catallax_numbers::{
    Decider, Rational, Refusal, operation_budget_active, with_operation_budget,
};
use catallax_recipe::{Operation, Recipe};

use crate::Plan;
pub use fit::{Surface, SurfaceKind};
use linalg::V3;

/// Most samples a scan may give the fitter. Above it the fit refuses with
/// `reconstruction.too-many-samples`.
pub const MAX_SAMPLES: usize = 1_000_000;
/// Fewest samples the fitter accepts.
pub const MIN_SAMPLES: usize = 100;
/// Maximum samples retained in either directed deviation report.
pub const MAX_DEVIATION_SAMPLES: usize = MAX_SAMPLES;
/// The part's display mesh for the deviation is certified within this many
/// millimetres of the B-rep.
pub const CHORD_MM: f64 = 0.001;

/// Fitting and sampling controls. None of them is a geometric tolerance of
/// the part; they choose how the scan is read and how densely deviation is sampled.
#[derive(Clone, Debug)]
pub struct ScanFitOptions {
    /// Threads for normal estimation. Answers are identical at any count.
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
    let workers = workers
        .clamp(1, n.max(1))
        .min(crate::workers::worker_count());
    if workers == 1 {
        return (0..n).map(f).collect();
    }
    let indices: Vec<usize> = (0..n).collect();
    crate::workers::with_workers(workers, || crate::workers::map(&indices, |&i| f(i)))
        .expect("a positive worker count")
}

pub(super) fn charge_sort(count: usize, dec: &mut Decider) -> Result<(), Refusal> {
    let levels = u64::from(count.checked_ilog2().unwrap_or(0)) + 1;
    dec.charge_work((count as u64).saturating_mul(levels))
}

fn stats(mut values: Vec<f64>, dec: &mut Decider) -> Result<DeviationSide, Refusal> {
    let samples = values.len();
    let mut sum = 0.0;
    for chunk in values.chunks(256) {
        dec.charge_work(chunk.len() as u64)?;
        for v in chunk {
            sum += v * v;
        }
    }
    let rms = (sum / samples.max(1) as f64).sqrt();
    charge_sort(samples, dec)?;
    values.sort_by(f64::total_cmp);
    dec.check_deadline()?;
    Ok(DeviationSide {
        samples,
        max_mm: values.last().copied().unwrap_or(0.0),
        rms_mm: rms,
        p95_mm: values
            .get((samples * 95).div_ceil(100).saturating_sub(1))
            .copied()
            .unwrap_or(0.0),
    })
}

fn sampling_limit() -> Refusal {
    Refusal::new(
        "reconstruction.deviation-too-large",
        format!(
            "each deviation direction is limited to {MAX_DEVIATION_SAMPLES} samples; sampling counts and spacing must be finite"
        ),
    )
}

fn sample_count(count: f64) -> Result<usize, Refusal> {
    if !count.is_finite() || count < 0.0 || count > MAX_DEVIATION_SAMPLES as f64 {
        return Err(sampling_limit());
    }
    Ok(count as usize)
}

/// Reserves the complete sample count before allocating, with checkpoints
/// during generation. The recurrence and sample order are unchanged.
fn triangle_samples(
    t: [V3; 3],
    spacing: f64,
    out: &mut Vec<V3>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(sampling_limit());
    }
    let at = |a: f64, b: f64| -> V3 {
        std::array::from_fn(|k| t[0][k] + a * (t[1][k] - t[0][k]) + b * (t[2][k] - t[0][k]))
    };
    let mut edges = [0; 3];
    for e in 0..3 {
        let length = linalg::norm(linalg::sub(t[(e + 1) % 3], t[e]));
        if !length.is_finite() {
            return Err(sampling_limit());
        }
        edges[e] = sample_count((length / spacing).ceil().max(1.0))?;
    }
    let area = linalg::norm(linalg::cross(
        linalg::sub(t[1], t[0]),
        linalg::sub(t[2], t[0]),
    )) / 2.0;
    let m = sample_count((area / (spacing * spacing)).floor())?;
    let count = edges
        .iter()
        .sum::<usize>()
        .checked_add(m)
        .ok_or_else(sampling_limit)?;
    if out
        .len()
        .checked_add(count)
        .is_none_or(|n| n > MAX_DEVIATION_SAMPLES)
    {
        return Err(sampling_limit());
    }
    dec.charge_work(count as u64)?;
    out.try_reserve(count).map_err(|_| sampling_limit())?;
    for e in 0..3 {
        let (p, q) = (t[e], t[(e + 1) % 3]);
        let n = edges[e];
        for i in 0..n {
            if i % 256 == 0 {
                dec.check_deadline()?;
            }
            let s = i as f64 / n as f64;
            out.push(std::array::from_fn(|k| p[k] + s * (q[k] - p[k])));
        }
    }
    const A1: f64 = 0.754_877_666_246_692_7;
    const A2: f64 = 0.569_840_290_998_053_3;
    for i in 1..=m {
        if i % 256 == 1 {
            dec.check_deadline()?;
        }
        let (mut a, mut b) = ((i as f64 * A1).fract(), (i as f64 * A2).fract());
        if a + b > 1.0 {
            (a, b) = (1.0 - a, 1.0 - b);
        }
        out.push(at(a, b));
    }
    dec.check_deadline()
}

fn deviation(
    scan: &ScanReference,
    samples: &[V3],
    brep: &catallax_features::Built,
    options: &ScanFitOptions,
    dec: &mut Decider,
) -> Result<ScanDeviation, Refusal> {
    if !options.part_spacing_mm.is_finite() {
        return Err(sampling_limit());
    }
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
    if mesh.positions.len() > MAX_DEVIATION_SAMPLES || mesh.triangles.len() > MAX_DEVIATION_SAMPLES
    {
        return Err(sampling_limit());
    }
    let chord_mm = mesh
        .faces
        .iter()
        .map(|f| f.max_chord_error.to_f64_nearest())
        .fold(0.0, f64::max);
    dec.charge_work(mesh.triangles.len() as u64)?;
    let faces: Vec<usize> = mesh
        .triangles
        .iter()
        .map(|t| t.face.index() as usize)
        .collect();
    dec.charge_work((mesh.positions.len() + mesh.triangles.len()) as u64)?;
    let part = FastSurface::new_checked(
        mesh.positions.clone(),
        mesh.triangles
            .iter()
            .map(|t| ScanItem::Triangle(t.vertices.map(|v| v as u32)))
            .collect(),
        dec,
    )?;
    let steps = options.triangle_steps.clamp(1, 16);
    let interiors = (steps - 1) * steps.saturating_sub(2) / 2;
    dec.charge_work(scan.items().len() as u64)?;
    let triangles = scan
        .items()
        .iter()
        .filter(|i| matches!(i, ScanItem::Triangle(_)))
        .count();
    let extra_count = triangles
        .checked_mul(interiors)
        .ok_or_else(sampling_limit)?;
    if samples
        .len()
        .checked_add(extra_count)
        .is_none_or(|n| n > MAX_DEVIATION_SAMPLES)
    {
        return Err(sampling_limit());
    }
    dec.charge_work((samples.len() + extra_count) as u64)?;
    let mut scan_map = Vec::with_capacity(samples.len());
    for &p in samples {
        scan_map.push(part.distance_checked(p, dec)?);
    }
    let mut extra_map = Vec::with_capacity(extra_count);
    let points = scan.approx_points();
    for item in scan.items() {
        dec.charge_work(1)?;
        if let ScanItem::Triangle(t) = item {
            let p = t.map(|v| points[v as usize]);
            for i in 1..steps {
                for j in 1..steps - i {
                    let (a, b) = (i as f64 / steps as f64, j as f64 / steps as f64);
                    let point = std::array::from_fn(|k| {
                        p[0][k] + a * (p[1][k] - p[0][k]) + b * (p[2][k] - p[0][k])
                    });
                    extra_map.push(part.distance_checked(point, dec)?);
                }
            }
        }
    }
    let mut part_samples = Vec::new();
    let mut part_faces = Vec::new();
    let spacing = options.part_spacing_mm.max(1e-3);
    for (t, face) in mesh.triangles.iter().zip(&faces) {
        let before = part_samples.len();
        triangle_samples(
            t.vertices.map(|v| mesh.positions[v]),
            spacing,
            &mut part_samples,
            dec,
        )?;
        dec.charge_work((part_samples.len() - before) as u64)?;
        part_faces.resize(part_samples.len(), *face);
    }
    dec.charge_work(part_samples.len() as u64)?;
    let mut part_distances = Vec::with_capacity(part_samples.len());
    for &p in &part_samples {
        part_distances.push(scan.surface().distance_checked(p, dec)?);
    }
    for d in scan_map.iter().chain(&extra_map).chain(&part_distances) {
        dec.charge_work(1)?;
        if !d.is_finite() {
            return Err(Refusal::new(
                "reconstruction.nonfinite-deviation",
                "a sampled deviation is not a finite number",
            ));
        }
    }
    dec.charge_work((scan_map.len() + extra_map.len() + part_distances.len()) as u64)?;
    let scan_to_part = stats(scan_map.iter().chain(&extra_map).copied().collect(), dec)?;
    let part_to_scan = stats(part_distances.clone(), dec)?;
    dec.charge_work(part_samples.len() as u64)?;
    let part_map = part_samples
        .into_iter()
        .zip(part_distances)
        .zip(part_faces)
        .map(|((p, d), f)| (p, d, f))
        .collect();
    dec.check_deadline()?;
    Ok(ScanDeviation {
        scan_to_part,
        part_to_scan,
        chord_mm,
        scan_map,
        part_map,
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
    if operation_budget_active() {
        let result = fit_scan_inner(scan_step, bytes, options, dec);
        dec.check_deadline()?;
        return result;
    }
    let (result, work) = with_operation_budget(dec.policy().max_work, || {
        fit_scan_inner(scan_step, bytes, options, dec)
    });
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "scan fitting exhausted the caller allowance",
        ));
    }
    dec.check_deadline()?;
    result
}

fn fit_scan_inner(
    scan_step: &Operation,
    bytes: &[u8],
    options: &ScanFitOptions,
    dec: &mut Decider,
) -> Result<ScanFit, Refusal> {
    let scan = ScanReference::from_step(scan_step, bytes, dec)?;
    dec.charge_work(scan.approx_points().len() as u64)?;
    let samples: Vec<V3> = scan.approx_points().to_vec();
    for chunk in samples.chunks(256) {
        dec.charge_work(chunk.len() as u64)?;
        if chunk.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Refusal::new(
                "reconstruction.nonfinite-sample",
                "a scan point is outside the finite binary64 fitting range",
            ));
        }
    }
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
    let cloud = segment::Cloud::new(&samples, options.noise_mm, workers, dec)?;
    let found = segment::segment(&cloud, dec)?;
    // Fillets refitted tangent to their stock faces, and small surfaces the
    // block family cannot place, change what explains each sample, so the
    // samples are assigned again before the final assembly.
    let mut found = found;
    if let Ok(first) = assemble::assemble(scan_step, &found.surfaces, &samples, cloud.sigma, dec) {
        let small = (samples.len() / 100).max(segment::MIN_SURFACE);
        dec.charge_work(samples.len() as u64)?;
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
        found = segment::assign(&cloud, &mut surfaces, dec)?;
    }
    let assembly = assemble::assemble(scan_step, &found.surfaces, &samples, cloud.sigma, dec);
    dec.check_deadline()?;
    let unexplained: Vec<UnexplainedRegion> = segment::unexplained(&cloud, &found.assignment, dec)?
        .into_iter()
        .map(|group| {
            let (centroid, _, _) =
                linalg::moments_checked(group.iter().map(|&i| samples[i as usize]), dec)?;
            Ok(UnexplainedRegion {
                samples: group,
                centroid,
            })
        })
        .collect::<Result<_, Refusal>>()?;
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
    let mut surfaces = Vec::with_capacity(reported.len());
    for ((surface, members), &role) in reported.iter().zip(&uses) {
        let mut max_mm = 0.0_f64;
        for &i in members {
            dec.charge_work(1)?;
            max_mm = max_mm.max(surface.residual(samples[i as usize]).abs());
        }
        surfaces.push(FittedSurface {
            surface: surface.clone(),
            samples: members.len(),
            rms_mm: fit::rms(surface, &samples, members, dec)?,
            max_mm,
            role,
        });
    }
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
    dec.check_deadline()?;
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

#[cfg(test)]
mod resource_tests {
    use super::*;

    #[test]
    fn scan_batches_respect_an_enclosing_single_worker_limit() {
        let caller = std::thread::current().id();
        let threads =
            crate::workers::with_workers(1, || par_map(4, 16, |_| std::thread::current().id()))
                .unwrap();
        assert!(threads.iter().all(|&thread| thread == caller));
    }

    #[test]
    fn huge_triangle_sample_counts_refuse_before_allocation() {
        let mut points = Vec::new();
        let refusal = triangle_samples(
            [[0.0; 3], [1000.0, 0.0, 0.0], [0.0, 1000.0, 0.0]],
            0.001,
            &mut points,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(refusal.code, "reconstruction.deviation-too-large");
        assert_eq!(points.capacity(), 0);
    }

    #[test]
    fn sample_generation_reserves_work_before_allocation() {
        let mut points = Vec::new();
        let (result, work) = with_operation_budget(16, || {
            triangle_samples(
                [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                0.001,
                &mut points,
                &mut Decider::default(),
            )
        });
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
        assert_eq!(points.capacity(), 0);
    }

    #[test]
    fn expired_sampling_deadlines_refuse_before_allocation() {
        let mut dec = Decider::default();
        dec.set_deadline(Some(
            std::time::Instant::now() - std::time::Duration::from_secs(1),
        ));
        let mut points = Vec::new();
        let refusal = triangle_samples(
            [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            0.25,
            &mut points,
            &mut dec,
        )
        .unwrap_err();
        assert_eq!(refusal.code, "numbers.deadline");
        assert_eq!(points.capacity(), 0);
    }
}
