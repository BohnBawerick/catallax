//! Unproved named-face selection from field chart preview triangles.
//!
//! This measures the nearest triangle in the supplied picture, including its
//! gaps. It does not certify the nearest defined face or a border distance.
//! Use `Prepared::nearest_face` for that proof and its typed ambiguities.
use crate::{
    batch::{Hardware, MAX_DISPATCH},
    gpu::Gpu,
};
use catallax_field::{FaceId, display::PreviewMesh};
use catallax_numbers::Refusal;

/// Maximum triangles in one selection buffer, within wgpu's default 128 MiB.
pub const MAX_TRIANGLES: usize = (128 << 20) / 48;

/// One unproved nearest-triangle selection. Ties use the first triangle.
#[derive(Clone, Debug, PartialEq)]
pub struct Pick {
    /// Field registry handle from the selected preview triangle.
    pub face: FaceId,
    /// Distance to that triangle in millimetres, without an error bound.
    pub distance: f64,
}

/// A batch of display selections, in input-point order.
#[derive(Clone, Debug, PartialEq)]
pub struct Picks {
    /// Unproved picks, never an export or measurement input.
    pub points: Vec<Pick>,
    /// Backend that produced these display estimates.
    pub hardware: Hardware,
    /// Cells missing from the supplied picture. Selections cannot account for
    /// those cells, even when a nearby triangle has a name.
    pub unresolved_cells: usize,
    /// Source triangles that collapsed in binary32 and were omitted.
    pub omitted_triangles: usize,
}
impl Picks {
    /// Status of every selection: always unproved.
    pub const STATUS: &'static str = "unproved";
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Triangle {
    vertices: [[f32; 4]; 3],
}

/// A snapshot of named field-preview triangles, rounded to binary32 for both
/// backends. It uses no P3 dimensions or hardcoded recipe names.
pub struct Picker {
    triangles: Vec<Triangle>,
    faces: Vec<FaceId>,
    unresolved_cells: usize,
    omitted_triangles: usize,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let edge = sub(b, a);
    let t = (dot(sub(p, a), edge) / dot(edge, edge)).clamp(0.0, 1.0);
    let d = sub(p, std::array::from_fn(|i| a[i] + t * edge[i]));
    dot(d, d)
}
fn distance(p: [f64; 3], triangle: &Triangle) -> f64 {
    let [a, b, c] = triangle
        .vertices
        .map(|v| [v[0] as f64, v[1] as f64, v[2] as f64]);
    let normal = cross(sub(b, a), sub(c, a));
    let norm = dot(normal, normal).sqrt();
    let n = normal.map(|v| v / norm);
    let value = dot(sub(p, a), n);
    let foot = std::array::from_fn(|i| p[i] - n[i] * value);
    if dot(cross(sub(b, a), sub(foot, a)), n) >= 0.0
        && dot(cross(sub(c, b), sub(foot, b)), n) >= 0.0
        && dot(cross(sub(a, c), sub(foot, c)), n) >= 0.0
    {
        value * value
    } else {
        segment_distance(p, a, b)
            .min(segment_distance(p, b, c))
            .min(segment_distance(p, c, a))
    }
}

fn rounded(point: [f64; 3]) -> Result<[f32; 4], Refusal> {
    let q = point.map(|v| v as f32);
    // The coordinate bound keeps squared cross products finite in binary32.
    if q.iter().any(|v| !v.is_finite() || v.abs() > 1e8) {
        return Err(Refusal::new(
            "preview.ownership-range",
            format!("point {point:?} exceeds finite binary32 selection range +/-1e8 mm"),
        ));
    }
    Ok([q[0], q[1], q[2], 0.0])
}

impl Picker {
    /// Copies the named triangles of one revision's field chart preview.
    /// An empty or fully collapsed mesh refuses with its counts. Gaps remain
    /// visible in every batch's metadata; no missing patch is filled.
    pub fn new(mesh: &PreviewMesh) -> Result<Self, Refusal> {
        if mesh.triangles.len() > MAX_TRIANGLES {
            return Err(Refusal::new(
                "preview.ownership-size",
                format!("{} triangles exceed {MAX_TRIANGLES}", mesh.triangles.len()),
            ));
        }
        let mut triangles = Vec::new();
        let mut faces = Vec::new();
        let mut omitted_triangles = 0;
        for t in &mesh.triangles {
            let mut vertices = [[0.0; 4]; 3];
            for (i, &index) in t.vertices.iter().enumerate() {
                vertices[i] = rounded(*mesh.positions.get(index).ok_or_else(|| {
                    Refusal::new(
                        "preview.ownership-mesh",
                        format!("vertex {index} exceeds {} positions", mesh.positions.len()),
                    )
                })?)?;
            }
            let [a, b, c] = vertices.map(|v| [v[0] as f64, v[1] as f64, v[2] as f64]);
            let edge =
                |a: [f64; 3], b: [f64; 3]| std::array::from_fn::<_, 3, _>(|i| (a[i] - b[i]) as f32);
            let ab = edge(b, a);
            let ac = edge(c, a);
            let n = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            let norm = (n[0] * n[0] + n[1] * n[1]) + n[2] * n[2];
            // Tiny triangles can collapse a squared normal in binary32.
            if !norm.is_finite() || norm < 1e-30 {
                omitted_triangles += 1;
                continue;
            }
            triangles.push(Triangle { vertices });
            faces.push(t.face);
        }
        if triangles.is_empty() {
            return Err(Refusal::new(
                "preview.ownership-empty",
                format!(
                    "{} source triangles, {omitted_triangles} collapsed, {} unresolved cells",
                    mesh.triangles.len(),
                    mesh.unresolved_cells
                ),
            ));
        }
        Ok(Self {
            triangles,
            faces,
            unresolved_cells: mesh.unresolved_cells,
            omitted_triangles,
        })
    }

    /// Selects the nearest named triangle to each point. The graphics-card
    /// route has one invocation per point. The processor route visits triangles
    /// in the same order at every worker count. Neither route is proved.
    /// `max_tests` bounds point-triangle comparisons before any dispatch.
    pub fn points(
        &self,
        points: &[[f64; 3]],
        gpu: Option<&Gpu>,
        workers: usize,
        max_tests: usize,
    ) -> Result<Picks, Refusal> {
        if workers == 0 {
            return Err(Refusal::new(
                "preview.ownership-workers",
                "selection requires at least one worker",
            ));
        }
        let count = points.len().checked_mul(self.triangles.len());
        if count.is_none_or(|n| n > max_tests) || points.len() > MAX_DISPATCH {
            return Err(Refusal::new(
                "preview.ownership-size",
                format!(
                    "{} points x {} triangles exceed {max_tests} comparisons or {MAX_DISPATCH} points",
                    points.len(),
                    self.triangles.len()
                ),
            ));
        }
        let points = points
            .iter()
            .copied()
            .map(rounded)
            .collect::<Result<Vec<_>, _>>()?;
        let selected = if points.is_empty() {
            Vec::new()
        } else if let Some(gpu) = gpu {
            let output = vec![[0u32; 4]; points.len()];
            let values: Vec<[u32; 4]> = gpu.dispatch(
                include_str!("picking.wgsl"),
                "select_triangles",
                &[
                    bytemuck::cast_slice(&points),
                    bytemuck::cast_slice(&self.triangles),
                    bytemuck::cast_slice(&output),
                ],
                2,
                [points.len().div_ceil(64) as u32, 1, 1],
            )?;
            values
                .into_iter()
                .map(|v| (v[0] as usize, f64::from(f32::from_bits(v[1]))))
                .collect()
        } else {
            let mut output = vec![(0, 0.0); points.len()];
            let chunk = points.len().div_ceil(workers);
            std::thread::scope(|scope| {
                let mut handles = Vec::new();
                for (points, output) in points.chunks(chunk).zip(output.chunks_mut(chunk)) {
                    handles.push(
                        std::thread::Builder::new()
                            .spawn_scoped(scope, move || {
                                for (p, result) in points.iter().zip(output) {
                                    let p = [p[0] as f64, p[1] as f64, p[2] as f64];
                                    let mut best = f64::INFINITY;
                                    for (i, t) in self.triangles.iter().enumerate() {
                                        let d = distance(p, t);
                                        if d < best {
                                            best = d;
                                            *result = (i + 1, d);
                                        }
                                    }
                                }
                            })
                            .map_err(|e| Refusal::new("preview.worker-start", e.to_string()))?,
                    );
                }
                for handle in handles {
                    handle
                        .join()
                        .map_err(|_| Refusal::new("preview.worker", "selection worker panicked"))?;
                }
                Ok::<_, Refusal>(())
            })?;
            output
        };
        let points = selected
            .into_iter()
            .map(|(i, d)| {
                if i == 0 || i > self.faces.len() || !d.is_finite() || d < 0.0 {
                    return Err(Refusal::new(
                        "preview.ownership-non-finite",
                        format!("triangle={i}, squared distance={d}"),
                    ));
                }
                Ok(Pick {
                    face: self.faces[i - 1],
                    distance: d.sqrt(),
                })
            })
            .collect::<Result<_, _>>()?;
        Ok(Picks {
            points,
            hardware: match gpu {
                Some(gpu) => Hardware::GraphicsCard {
                    adapter: gpu.name().to_string(),
                },
                None => Hardware::Processor {
                    workers,
                    reason: "processor selection requested".into(),
                },
            },
            unresolved_cells: self.unresolved_cells,
            omitted_triangles: self.omitted_triangles,
        })
    }
}
