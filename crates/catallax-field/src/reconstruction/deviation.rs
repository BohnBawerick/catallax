//! Explicitly unproved, finite-sample fidelity diagnostics in both directions.
use super::{Reconstruction, fail};
use catallax_import::mesh::ClosedMesh;
use catallax_numbers::{Decider, Refusal};

/// Deterministic sampling controls; these are not accuracy tolerances.
#[derive(Clone, Copy)]
pub struct Sampling {
    /// Barycentric subdivisions on every source triangle (1..=16).
    pub triangle_steps: usize,
    /// Uniform subdivisions on each analytic surface parameter (2..=64).
    pub surface_steps: usize,
    /// Angular subdivisions of the bore (8..=512).
    pub angles: usize,
}
/// Finite binary64 sample maxima, NEVER certified global Hausdorff bounds.
#[derive(Debug, Clone, PartialEq)]
pub struct SampledDeviation {
    /// Maximum sampled source-to-analytic-boundary distance in mm, unproved.
    pub scan_to_part_mm: f64,
    /// Maximum sampled analytic-boundary-to-source distance in mm, unproved.
    pub part_to_scan_mm: f64,
    /// Source samples, including duplicates at triangle edges.
    pub scan_samples: usize,
    /// Analytic boundary samples, including duplicates at surface joins.
    pub part_samples: usize,
}
struct Shape {
    lo: [f64; 3],
    hi: [f64; 3],
    c: [f64; 2],
    r: [f64; 2],
}
impl Shape {
    fn distance(&self, p: [f64; 3], bore: bool) -> f64 {
        let mut best = f64::INFINITY;
        for face in 0..6 {
            let a = face / 2;
            let mut q: [f64; 3] = std::array::from_fn(|i| p[i].clamp(self.lo[i], self.hi[i]));
            q[a] = if face % 2 == 0 {
                self.lo[a]
            } else {
                self.hi[a]
            };
            if a == 2 && bore {
                let r = self.r[face % 2];
                let dx = q[0] - self.c[0];
                let dy = q[1] - self.c[1];
                let rho = dx.hypot(dy);
                if rho < r {
                    if rho > 0. {
                        q[0] = self.c[0] + r * dx / rho;
                        q[1] = self.c[1] + r * dy / rho;
                    } else {
                        q[0] = self.c[0] + r;
                        q[1] = self.c[1];
                    }
                }
            }
            best = best.min((0..3).map(|i| (p[i] - q[i]).powi(2)).sum::<f64>().sqrt());
        }
        if bore {
            let rho = (p[0] - self.c[0]).hypot(p[1] - self.c[1]);
            let h = self.hi[2] - self.lo[2];
            let dr = self.r[1] - self.r[0];
            let denominator = h * h + dr * dr;
            if !denominator.is_finite() || denominator <= 0. {
                return f64::NAN;
            }
            let t =
                (((rho - self.r[0]) * dr + (p[2] - self.lo[2]) * h) / denominator).clamp(0., 1.);
            best = best.min((rho - self.r[0] - t * dr).hypot(p[2] - self.lo[2] - t * h));
        }
        best
    }
}
/// Compare the entire source triangle set and all seven fitted surface domains
/// using finite samples. Source bytes/scale must match the reconstruction.
pub fn sampled_deviation(
    bytes: &[u8],
    part: &Reconstruction,
    sampling: Sampling,
    dec: &mut Decider,
) -> Result<SampledDeviation, Refusal> {
    compare(bytes, part, sampling, true, dec)
}
/// Adverse diagnostic only: compare the source against the filled outer box.
/// This deliberately removes the inferred bore; it does not return a reconstructed
/// solid and must never replace the fitted part or its normal fidelity report.
pub fn filled_envelope_control(
    bytes: &[u8],
    part: &Reconstruction,
    sampling: Sampling,
    dec: &mut Decider,
) -> Result<SampledDeviation, Refusal> {
    compare(bytes, part, sampling, false, dec)
}
fn compare(
    bytes: &[u8],
    part: &Reconstruction,
    s: Sampling,
    bore: bool,
    dec: &mut Decider,
) -> Result<SampledDeviation, Refusal> {
    if !(1..=16).contains(&s.triangle_steps)
        || !(2..=64).contains(&s.surface_steps)
        || !(8..=512).contains(&s.angles)
    {
        return Err(fail("sampling controls exceed the explicit bounded ranges"));
    }
    if catallax_import::sha256_hex(bytes) != part.source_sha256 {
        return Err(fail("deviation source identity differs from fitted source"));
    }
    let mesh = ClosedMesh::read_with_limit(bytes, &part.source_scale, super::MAX_TRIANGLES, dec)?;
    let shape = Shape {
        lo: part.bounds.each_ref().map(|b| b[0].to_f64_nearest()),
        hi: part.bounds.each_ref().map(|b| b[1].to_f64_nearest()),
        c: part.centre.each_ref().map(|v| v.to_f64_nearest()),
        r: part.radii.each_ref().map(|v| v.to_f64_nearest()),
    };
    if shape
        .lo
        .iter()
        .chain(&shape.hi)
        .chain(&shape.c)
        .chain(&shape.r)
        .any(|v| !v.is_finite())
        || (0..3).any(|i| shape.hi[i] - shape.lo[i] <= 0.)
        || shape.r.iter().any(|r| *r <= 0.)
    {
        return Err(fail(
            "fidelity sampling requires finite resolved binary64 dimensions",
        ));
    }
    let mut result = SampledDeviation {
        scan_to_part_mm: 0.,
        part_to_scan_mm: 0.,
        scan_samples: 0,
        part_samples: 0,
    };
    for t in mesh.triangles() {
        let p = t.map(|v| mesh.vertex(v).map(|r| r.to_f64_nearest()));
        for i in 0..=s.triangle_steps {
            for j in 0..=s.triangle_steps - i {
                dec.charge_work(1)?;
                let a = i as f64 / s.triangle_steps as f64;
                let b = j as f64 / s.triangle_steps as f64;
                let q = std::array::from_fn(|k| a * p[0][k] + b * p[1][k] + (1. - a - b) * p[2][k]);
                let d = shape.distance(q, bore);
                if !d.is_finite() {
                    return Err(fail("nonfinite sample distance"));
                }
                result.scan_to_part_mm = result.scan_to_part_mm.max(d);
                result.scan_samples += 1;
            }
        }
    }
    let mut sample = |p: [f64; 3]| -> Result<(), Refusal> {
        // The unproved fast distance has no Decider; conservatively charge every
        // triangle per sample rather than hiding this diagnostic's work.
        dec.charge_work(mesh.triangles().len() as u64)?;
        let d = mesh.fast_squared_distance(p).sqrt();
        if !d.is_finite() {
            return Err(fail("nonfinite source sample distance"));
        }
        result.part_to_scan_mm = result.part_to_scan_mm.max(d);
        result.part_samples += 1;
        Ok(())
    };
    for face in 0..6 {
        let axis = face / 2;
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for i in 0..=s.surface_steps {
            for j in 0..=s.surface_steps {
                let mut p = shape.lo;
                p[axis] = if face % 2 == 0 {
                    shape.lo[axis]
                } else {
                    shape.hi[axis]
                };
                p[u] += i as f64 / s.surface_steps as f64 * (shape.hi[u] - shape.lo[u]);
                p[v] += j as f64 / s.surface_steps as f64 * (shape.hi[v] - shape.lo[v]);
                if axis != 2
                    || !bore
                    || (p[0] - shape.c[0]).hypot(p[1] - shape.c[1]) >= shape.r[face % 2]
                {
                    sample(p)?;
                }
            }
        }
    }
    if bore {
        for k in 0..=s.surface_steps {
            for i in 0..s.angles {
                let t = k as f64 / s.surface_steps as f64;
                let angle = std::f64::consts::TAU * i as f64 / s.angles as f64;
                let r = shape.r[0] + t * (shape.r[1] - shape.r[0]);
                sample([
                    shape.c[0] + r * angle.cos(),
                    shape.c[1] + r * angle.sin(),
                    shape.lo[2] + t * (shape.hi[2] - shape.lo[2]),
                ])?;
            }
        }
    }
    Ok(result)
}
