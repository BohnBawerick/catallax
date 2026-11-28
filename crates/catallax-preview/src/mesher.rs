//! Unproved B-rep preview meshing on an application-owned wgpu device.

use std::{fmt, ops::Range};

use bytemuck::{Pod, Zeroable, cast_slice};
use catallax_mesh::float::{DISPLAY_DESCRIPTION_VERSION, DisplayDescription, DisplayFace, trim};
use catallax_topology::FaceId;

use crate::gpu::Gpu;

mod surface;

/// Lattice divisions in U and V. Trim chords are independent of this setting.
#[derive(Clone, Copy, Debug)]
pub struct FaceDetail(pub [u32; 2]);

impl Default for FaceDetail {
    fn default() -> Self {
        Self([48, 24])
    }
}

/// Why no GPU mesh is available for this body.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GpuRefusal {
    /// The app explicitly opted out of GPU preview.
    Disabled,
    /// The app could not select a GPU adapter for its renderer.
    NoAdapter,
    /// The supplied device cannot run the preview compute pipeline.
    UnsupportedDevice,
    /// The display description uses an unknown format version.
    DescriptionVersion(u32),
    /// The body's world-space origin is not finite.
    InvalidOrigin,
    /// Pipeline validation failed.
    Shader(String),
    /// The GPU could not allocate or submit the mesh.
    Device(String),
}

impl GpuRefusal {
    /// Stable code for logging a CPU fallback.
    pub fn code(&self) -> &'static str {
        match self {
            Self::Disabled => "preview.gpu-disabled",
            Self::NoAdapter => "preview.no-adapter",
            Self::UnsupportedDevice => "preview.no-device",
            Self::DescriptionVersion(_) => "preview.description-version",
            Self::InvalidOrigin => "preview.gpu-non-finite",
            Self::Shader(_) => "preview.shader",
            Self::Device(_) => "preview.gpu-device",
        }
    }
}

impl fmt::Display for GpuRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.code())?;
        match self {
            Self::DescriptionVersion(version) => write!(f, ": {version}"),
            Self::Shader(error) | Self::Device(error) => write!(f, ": {error}"),
            _ => Ok(()),
        }
    }
}

impl std::error::Error for GpuRefusal {}

/// Why this face needs today's CPU preview instead.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FaceRefusal {
    /// The chart bounds, period or requested detail are invalid.
    InvalidChart,
    /// The chart requires a trim shape not yet handled (e.g. a pole boundary).
    UnsupportedTrim,
    /// The trim did not close or could not be triangulated without changing a chord.
    InvalidTrim,
    /// The surface parameters or sample layout are invalid.
    InvalidSurface,
    /// A body-relative coordinate or a surface parameter is outside binary32 range.
    GpuRange,
    /// The face would exceed the available device buffer or dispatch limits.
    MeshSize,
}

impl FaceRefusal {
    /// Stable code for logging a CPU fallback.
    pub fn code(self) -> &'static str {
        match self {
            Self::InvalidChart => "preview.grid-range",
            Self::UnsupportedTrim => "preview.unsupported-trim",
            Self::InvalidTrim => "preview.invalid-trim",
            Self::InvalidSurface => "preview.invalid-surface",
            Self::GpuRange => "preview.gpu-range",
            Self::MeshSize => "preview.mesh-size",
        }
    }
}

impl From<trim::TrimRefusal> for FaceRefusal {
    fn from(error: trim::TrimRefusal) -> Self {
        use trim::TrimRefusal::*;
        match error {
            InvalidChart => Self::InvalidChart,
            UnsupportedTrim => Self::UnsupportedTrim,
            InvalidTrim(_) => Self::InvalidTrim,
            MeshSize => Self::MeshSize,
            InvalidPosition => Self::GpuRange,
        }
    }
}

impl fmt::Display for FaceRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.code())
    }
}

impl std::error::Error for FaceRefusal {}

/// Renderer vertex, with body-relative position and outward normal.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct PreviewVertex {
    /// XYZ in millimetres; W is 1.
    pub position: [f32; 4],
    /// Normal XYZ; W is 0.
    pub normal: [f32; 4],
}

impl PreviewVertex {
    /// Use with a vertex-step buffer and no base-vertex offset on indexed draws.
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        const ATTRIBUTES: [wgpu::VertexAttribute; 2] =
            wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3];
        wgpu::VertexBufferLayout {
            array_stride: size_of::<Self>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                ATTRIBUTES[0],
                wgpu::VertexAttribute {
                    offset: 16,
                    ..ATTRIBUTES[1]
                },
            ],
        }
    }
}

/// Buffers on the device passed to [`Mesher::new`]. No host copy is produced.
pub struct MeshBuffers {
    /// [`PreviewVertex`] records.
    pub vertex: wgpu::Buffer,
    /// `Uint32` triangle indices, relative to the start of `vertex`.
    pub index: wgpu::Buffer,
    /// Number of vertex records written by the compute pass.
    pub vertex_count: u32,
}

/// Successful indexed draw for one face.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FaceDraw {
    /// Index range in [`MeshBuffers::index`], in indices rather than bytes.
    pub indices: Range<u32>,
    /// True when the GPU interpolated a sampled surface instead of evaluating a carrier.
    pub sampled_grid: bool,
}

/// Outcome for a face, in display-description order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FacePreview<F = FaceId> {
    /// Face handle or stable definition key to use for names and the pick ID.
    pub face: F,
    /// Draw this range, or draw this face from the CPU preview.
    pub draw: Result<FaceDraw, FaceRefusal>,
}

/// GPU preview and the status of every requested face.
pub struct GpuMesh<F = FaceId> {
    /// World-space body origin in millimetres; add in the renderer's local frame.
    pub origin: [f64; 3],
    /// None when the description is empty or every face refused.
    pub buffers: Option<MeshBuffers>,
    /// One entry for every input face, even if a face refused.
    pub faces: Vec<FacePreview<F>>,
}

/// Result of [`preview_or_cpu`]. A preview is never a checked mesh.
pub enum Preview<T, F = FaceId> {
    /// All faces came from the GPU.
    Gpu(GpuMesh<F>),
    /// Draw successful GPU ranges and use the CPU mesh only for refused faces.
    Mixed {
        /// Includes the typed fallback reason for each refused face.
        gpu: GpuMesh<F>,
        /// Today's CPU preview for looking up each refused face.
        cpu: T,
    },
    /// Use today's CPU preview for the whole body.
    Cpu {
        /// Reason the GPU route could not run.
        reason: GpuRefusal,
        /// Today's CPU preview.
        mesh: T,
    },
}

/// Uses the app's current CPU preview on body-wide or per-face GPU refusal.
///
/// The app passes a refusal here when it has disabled GPU preview or could not
/// acquire a device. `cpu_preview` runs at most once, and never on GPU success.
pub async fn preview_or_cpu<T, E, F: Clone, K>(
    mesher: Result<&Mesher, GpuRefusal>,
    description: &DisplayDescription<F, K>,
    detail: impl FnMut(&F) -> FaceDetail,
    cpu_preview: impl FnOnce() -> Result<T, E>,
) -> Result<Preview<T, F>, E> {
    let result = match mesher {
        Ok(mesher) => mesher.mesh(description, detail).await,
        Err(reason) => Err(reason),
    };
    match result {
        Ok(gpu) if gpu.faces.iter().any(|f| f.draw.is_err()) => Ok(Preview::Mixed {
            gpu,
            cpu: cpu_preview()?,
        }),
        Ok(gpu) => Ok(Preview::Gpu(gpu)),
        Err(reason) => Ok(Preview::Cpu {
            reason,
            mesh: cpu_preview()?,
        }),
    }
}

/// Reusable compute pipeline on one renderer device.
pub struct Mesher {
    gpu: Gpu,
}

impl Mesher {
    /// Compile the preview pipeline using the app's existing device and queue.
    /// `queue` must be the queue obtained with `device`.
    pub async fn new(device: &wgpu::Device, queue: &wgpu::Queue) -> Result<Self, GpuRefusal> {
        Ok(Self {
            gpu: Gpu::new(device, queue).await?,
        })
    }

    /// Triangulate the trims in binary64, then write binary32 vertices on the GPU.
    /// Use the returned buffers on the same device after this call; it does not
    /// wait for GPU completion or map either buffer.
    pub async fn mesh<F: Clone, K>(
        &self,
        description: &DisplayDescription<F, K>,
        mut detail: impl FnMut(&F) -> FaceDetail,
    ) -> Result<GpuMesh<F>, GpuRefusal> {
        if description.version != DISPLAY_DESCRIPTION_VERSION {
            return Err(GpuRefusal::DescriptionVersion(description.version));
        }
        if description.origin.iter().any(|x| !x.is_finite()) {
            return Err(GpuRefusal::InvalidOrigin);
        }
        let mut packed = PackedMesh::new(&self.gpu, description);
        let faces = description
            .faces
            .iter()
            .map(|face| FacePreview {
                face: face.face.clone(),
                draw: packed.append(face, detail(&face.face).0),
            })
            .collect();
        let buffers = if packed.inputs.is_empty() {
            None
        } else {
            Some(
                self.gpu
                    .dispatch(
                        cast_slice(&packed.inputs),
                        cast_slice(&packed.positions),
                        cast_slice(&packed.surfaces),
                        cast_slice(&packed.indices),
                        cast_slice(&packed.knots),
                        packed.inputs.len() as u32,
                    )
                    .await?,
            )
        };
        Ok(GpuMesh {
            origin: description.origin,
            buffers,
            faces,
        })
    }
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct InputVertex {
    uv: [f32; 2],
    boundary: u32,
    surface: u32,
}

impl InputVertex {
    fn new(uv: [f64; 2], boundary: u32) -> Self {
        Self {
            uv: uv.map(|x| x as f32),
            boundary,
            surface: 0,
        }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default, Pod, Zeroable)]
struct SurfaceParams {
    origin: [f32; 4],
    x: [f32; 4],
    y: [f32; 4],
    z: [f32; 4],
    mode: [u32; 4],
    // SplineNet in eval/spline.wgsl, followed by padding to a vec4 boundary.
    spline: [u32; 8],
}

struct PackedMesh {
    edges: Vec<Result<Range<u32>, FaceRefusal>>,
    positions: Vec<[f32; 4]>,
    inputs: Vec<InputVertex>,
    indices: Vec<u32>,
    surfaces: Vec<SurfaceParams>,
    knots: Vec<f32>,
    max_positions: usize,
    max_indices: usize,
    max_vertices: usize,
    max_surfaces: usize,
    max_knots: usize,
}

impl trim::TrimBoundaries for PackedMesh {
    fn edge_range(&self, edge: usize) -> Option<Range<u32>> {
        self.edges.get(edge)?.as_ref().ok().cloned()
    }

    fn same_position(&self, a: u32, b: u32) -> bool {
        self.positions[a as usize].map(f32::to_bits) == self.positions[b as usize].map(f32::to_bits)
    }

    fn position(&self, sample: u32) -> [f64; 3] {
        let [x, y, z, _] = self.positions[sample as usize];
        [x, y, z].map(f64::from)
    }

    fn add_pole(&mut self, position: [f64; 3]) -> Result<u32, trim::TrimRefusal> {
        if self.positions.len() == self.max_positions {
            return Err(trim::TrimRefusal::MeshSize);
        }
        let index = self.positions.len() as u32;
        self.positions
            .push(position32(position).map_err(|_| trim::TrimRefusal::InvalidPosition)?);
        Ok(index)
    }

    fn grid_uses_f32(&self) -> bool {
        true
    }
}

impl PackedMesh {
    fn new<F, K>(gpu: &Gpu, desc: &DisplayDescription<F, K>) -> Self {
        let max_positions = gpu.storage_bytes / size_of::<[f32; 4]>();
        let mut positions = Vec::new();
        let edges = desc
            .edges
            .iter()
            .map(|edge| {
                let first = positions.len();
                if edge.points.len() < 2 {
                    return Err(FaceRefusal::InvalidTrim);
                }
                if edge.points.len() > max_positions - first {
                    return Err(FaceRefusal::MeshSize);
                }
                positions.reserve(edge.points.len());
                for &point in &edge.points {
                    match position32(point) {
                        Ok(point) => positions.push(point),
                        Err(error) => {
                            positions.truncate(first);
                            return Err(error);
                        }
                    }
                }
                Ok(first as u32..positions.len() as u32)
            })
            .collect();
        Self {
            edges,
            positions,
            inputs: Vec::new(),
            indices: Vec::new(),
            surfaces: Vec::new(),
            knots: Vec::new(),
            max_positions,
            max_indices: gpu.index_bytes / size_of::<u32>(),
            max_vertices: gpu.max_vertices(),
            max_surfaces: gpu.storage_bytes / size_of::<SurfaceParams>(),
            max_knots: gpu.storage_bytes / size_of::<f32>(),
        }
    }

    fn append<F>(
        &mut self,
        face: &DisplayFace<F>,
        detail: [u32; 2],
    ) -> Result<FaceDraw, FaceRefusal> {
        if self.surfaces.len() == self.max_surfaces {
            return Err(FaceRefusal::MeshSize);
        }
        for edge in face.loops.iter().flat_map(|lp| &lp.uses) {
            self.edges
                .get(edge.edge)
                .ok_or(FaceRefusal::InvalidTrim)?
                .as_ref()
                .map_err(|&error| error)?;
        }
        let first = self.positions.len();
        let max_vertices = self.max_vertices - self.inputs.len();
        let max_indices = self.max_indices - self.indices.len();
        let trim = match trim::triangulate(face, self, detail, max_vertices, max_indices) {
            Ok(trim) => trim,
            Err(error) => {
                self.positions.truncate(first);
                return Err(error.into());
            }
        };
        let first_surface = self.positions.len();
        let first_knot = self.knots.len();
        let params = match surface::pack(
            face,
            trim.bounds,
            trim.sphere_chart,
            &mut self.positions,
            self.max_positions,
            &mut self.knots,
            self.max_knots,
        ) {
            Ok(params) => params,
            Err(error) => {
                self.positions.truncate(first_surface);
                self.knots.truncate(first_knot);
                if trim.sphere_chart.is_some()
                    || !matches!(error, FaceRefusal::GpuRange | FaceRefusal::UnsupportedTrim)
                {
                    self.positions.truncate(first);
                    return Err(error);
                }
                match surface::sampled(
                    face,
                    trim.bounds,
                    detail,
                    &mut self.positions,
                    self.max_positions,
                ) {
                    Ok(params) => params,
                    Err(_) => {
                        self.positions.truncate(first);
                        return Err(error);
                    }
                }
            }
        };
        let base = self.inputs.len() as u32;
        let slot = self.surfaces.len() as u32;
        let start = self.indices.len() as u32;
        self.inputs.extend(trim.vertices.into_iter().map(|sample| {
            let mut vertex = InputVertex::new(sample.uv, sample.boundary);
            vertex.surface = slot;
            vertex
        }));
        self.indices
            .extend(trim.indices.into_iter().map(|index| base + index));
        self.surfaces.push(params);
        Ok(FaceDraw {
            indices: start..self.indices.len() as u32,
            sampled_grid: params.mode[0] == surface::GRID,
        })
    }
}

fn scalar32(x: f64) -> Result<f32, FaceRefusal> {
    let x = x as f32;
    if x.is_finite() {
        // Endpoints may arrive through either of their topological edges.
        Ok(if x == 0.0 { 0.0 } else { x })
    } else {
        Err(FaceRefusal::GpuRange)
    }
}

fn position32(p: [f64; 3]) -> Result<[f32; 4], FaceRefusal> {
    Ok([scalar32(p[0])?, scalar32(p[1])?, scalar32(p[2])?, 1.0])
}
