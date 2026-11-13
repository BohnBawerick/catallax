//! Unproved batch point evaluation on the graphics card, with the processor
//! preview as the fallback route.
//!
//! Every value here is a display estimate. It carries no proof or error bound
//! and never feeds an export, a proof or a proven measurement.
use crate::gpu::Gpu;
use catallax_field::{Plan, UnprovedBatch, Workspace};
use catallax_numbers::Refusal;

/// Largest point count one dispatch covers: 65,535 workgroups of 64.
pub const MAX_DISPATCH: usize = 65_535 * 64;
/// Processor grid points generated per batch call.
const PROCESSOR_CHUNK: usize = 1 << 20;
/// Grid indices stay below this magnitude so binary32 holds `index + 1/2` exactly.
const INDEX_LIMIT: i64 = 1 << 23;
/// Pitch exponents keep every grid coordinate a normal binary32 number.
const PITCH_EXPONENTS: std::ops::RangeInclusive<i32> = -40..=40;

/// Where a batch of unproved preview values was computed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Hardware {
    /// Binary32 field arithmetic in a wgpu compute shader.
    GraphicsCard {
        /// Adapter name, backend and device type.
        adapter: String,
    },
    /// Processor display arithmetic in binary64. Field batches use
    /// [`Plan::unproved_preview_batch`]; triangle picks measure the picture.
    Processor {
        /// Worker threads used.
        workers: usize,
        /// Why the graphics card was not used.
        reason: String,
    },
}

/// Side of the surface reported by one unproved value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    /// The value is negative.
    Inside,
    /// The value is positive.
    Outside,
    /// The value is zero. A Boolean residual can be zero away from the boundary.
    Zero,
}

/// A cell-centred regular grid whose spacing is a power of two. Point
/// `(i, j, k)` is `(origin + (i, j, k) + 1/2) * 2^pitch_exponent`, exact in
/// binary32 and binary64, so both routes evaluate identical points. The half
/// step keeps samples off faces at whole multiples of the pitch, where a value
/// of zero would make the side depend on the last bit. Index `i + nx * (j + ny * k)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Grid {
    /// Integer coordinates of the first point's cell, in pitches.
    pub origin: [i32; 3],
    /// Points along each axis.
    pub counts: [u32; 3],
    /// The spacing is `2^pitch_exponent` millimetres.
    pub pitch_exponent: i32,
}

impl Grid {
    /// The smallest grid with this pitch that covers `[lo, hi]` plus `margin`
    /// whole pitches on every side: the outermost points lie at least
    /// `margin - 1/2` pitches outside the bounds. Refuses with
    /// `preview.grid-range` when an index or coordinate would leave the exact
    /// binary32 range.
    pub fn covering(
        lo: [f64; 3],
        hi: [f64; 3],
        pitch_exponent: i32,
        margin: u32,
    ) -> Result<Self, Refusal> {
        let range = |detail: String| Refusal::new("preview.grid-range", detail);
        if !PITCH_EXPONENTS.contains(&pitch_exponent) {
            return Err(range(format!(
                "pitch exponent {pitch_exponent} outside {PITCH_EXPONENTS:?}"
            )));
        }
        let pitch = 2f64.powi(pitch_exponent);
        let mut origin = [0; 3];
        let mut counts = [0; 3];
        for axis in 0..3 {
            if !(lo[axis].is_finite() && hi[axis].is_finite() && lo[axis] <= hi[axis]) {
                return Err(range(format!(
                    "axis {axis} bounds are not a finite interval"
                )));
            }
            let first = (lo[axis] / pitch).floor() - f64::from(margin);
            let last = (hi[axis] / pitch).ceil() + f64::from(margin) - 1.0;
            if first <= -(INDEX_LIMIT as f64) || last >= INDEX_LIMIT as f64 || last < first {
                return Err(range(format!(
                    "axis {axis} needs grid indices beyond 2^23 at this pitch"
                )));
            }
            origin[axis] = first as i32;
            counts[axis] = (last - first) as u32 + 1;
        }
        let grid = Self {
            origin,
            counts,
            pitch_exponent,
        };
        grid.check_count()?;
        Ok(grid)
    }

    /// Spacing in millimetres.
    pub fn pitch(&self) -> f64 {
        2f64.powi(self.pitch_exponent)
    }

    /// Number of points, or `None` when the product of the axis counts
    /// overflows `usize`.
    pub fn checked_len(&self) -> Option<usize> {
        self.counts
            .iter()
            .try_fold(1usize, |n, &c| n.checked_mul(c as usize))
    }

    /// Number of points, saturating at `usize::MAX` when the product of the
    /// axis counts overflows. Every evaluation and size check refuses such a
    /// grid with a typed code.
    pub fn len(&self) -> usize {
        self.checked_len().unwrap_or(usize::MAX)
    }

    /// Refuses with `preview.grid-range` unless the point count fits u32
    /// indices, using checked multiplication.
    fn check_count(&self) -> Result<(), Refusal> {
        match self.checked_len() {
            Some(n) if n <= u32::MAX as usize => Ok(()),
            n => Err(Refusal::new(
                "preview.grid-range",
                format!(
                    "{} grid points ({} x {} x {}) exceed u32 indices",
                    n.map_or_else(|| "more than usize::MAX".to_string(), |n| n.to_string()),
                    self.counts[0],
                    self.counts[1],
                    self.counts[2]
                ),
            )),
        }
    }

    /// Whether the grid has no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The point at integer offsets `(i, j, k)` from the origin.
    pub fn at(&self, i: u32, j: u32, k: u32) -> [f64; 3] {
        let pitch = self.pitch();
        [
            (f64::from(self.origin[0] + i as i32) + 0.5) * pitch,
            (f64::from(self.origin[1] + j as i32) + 0.5) * pitch,
            (f64::from(self.origin[2] + k as i32) + 0.5) * pitch,
        ]
    }

    /// The point at a linear index.
    pub fn point(&self, index: usize) -> [f64; 3] {
        let [nx, ny, _] = self.counts.map(|n| n as usize);
        self.at(
            (index % nx) as u32,
            ((index / nx) % ny) as u32,
            (index / (nx * ny)) as u32,
        )
    }

    fn check(&self) -> Result<(), Refusal> {
        let range = |detail: &str| Refusal::new("preview.grid-range", detail.to_string());
        if !PITCH_EXPONENTS.contains(&self.pitch_exponent) {
            return Err(range("pitch exponent out of range"));
        }
        for axis in 0..3 {
            let first = i64::from(self.origin[axis]);
            let last = first + i64::from(self.counts[axis]) - 1;
            if first <= -INDEX_LIMIT || last >= INDEX_LIMIT {
                return Err(range("grid indices beyond 2^23"));
            }
        }
        self.check_count()
    }
}

/// Unproved values in input order, labelled with the hardware that made them.
#[derive(Clone, Debug, PartialEq)]
pub struct Samples {
    /// Negative inside, positive outside, canonical `+0.0` for zero. Graphics
    /// card values are binary32 widened exactly to binary64.
    pub values: Vec<f64>,
    /// Whether the real-arithmetic formula is a true signed distance. When
    /// false, the magnitude is a Boolean residual, not a distance.
    pub true_distance_formula: bool,
    /// Where the values were computed.
    pub hardware: Hardware,
}

impl Samples {
    /// Proof status of every value: always unproved.
    pub const STATUS: &'static str = "unproved";

    /// Inside, outside or zero for value `index`.
    pub fn side(&self, index: usize) -> Side {
        let v = self.values[index];
        if v < 0.0 {
            Side::Inside
        } else if v > 0.0 {
            Side::Outside
        } else {
            Side::Zero
        }
    }

    /// Unsigned estimate for value `index`: a distance estimate only when
    /// [`Self::true_distance_formula`] holds, otherwise a residual.
    pub fn distance_estimate(&self, index: usize) -> f64 {
        self.values[index].abs()
    }
}

/// Tuning for an [`Evaluator`]. These settings never change a value.
#[derive(Clone, Debug)]
pub struct Options {
    /// Processor threads for the fallback route. Zero refuses.
    pub workers: usize,
    /// Points per graphics-card dispatch, at most [`MAX_DISPATCH`].
    pub dispatch: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            workers: std::thread::available_parallelism().map_or(1, usize::from),
            dispatch: MAX_DISPATCH,
        }
    }
}

struct Kernels<'g> {
    gpu: &'g Gpu,
    points: wgpu::ComputePipeline,
    grid: wgpu::ComputePipeline,
}

/// Batch evaluator for one compiled plan. It uses the graphics card when one
/// is given and the plan translates to WGSL, and the processor preview
/// otherwise. The chosen route is fixed for the evaluator's lifetime.
pub struct Evaluator<'p, 'g> {
    plan: &'p Plan,
    kernels: Option<Kernels<'g>>,
    hardware: Hardware,
    options: Options,
    true_distance_formula: bool,
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    origin: [i32; 4],
    counts: [u32; 4],
    start: u32,
    len: u32,
    pitch: f32,
    pad: u32,
}

impl<'p, 'g> Evaluator<'p, 'g> {
    /// Chooses the route. Without `gpu`, or when the plan has no WGSL
    /// translation (`field.preview-wgsl-*`), the processor preview is used and
    /// the reason is recorded in [`Self::hardware`]. A shader that fails
    /// validation is the typed refusal `preview.shader`, not a fallback.
    pub fn new(plan: &'p Plan, gpu: Option<&'g Gpu>, options: Options) -> Result<Self, Refusal> {
        if options.workers == 0 {
            return Err(Refusal::new(
                "field.preview-workers",
                "unproved preview requires at least one worker",
            ));
        }
        if options.dispatch == 0 || options.dispatch > MAX_DISPATCH {
            return Err(Refusal::new(
                "preview.dispatch",
                format!("dispatch size must be between 1 and {MAX_DISPATCH}"),
            ));
        }
        // The processor preview reports the plan's formula flag with each value.
        let true_distance_formula = plan
            .unproved_preview([0.0; 3], &mut Workspace::default())
            .is_ok_and(|v| v.true_distance_formula);
        let processor = |reason: String| Hardware::Processor {
            workers: options.workers,
            reason,
        };
        let Some(gpu) = gpu else {
            return Ok(Self {
                plan,
                kernels: None,
                hardware: processor("no graphics card device".into()),
                options,
                true_distance_formula,
            });
        };
        let field = match plan.unproved_preview_wgsl() {
            Ok(source) => source,
            Err(refusal) if refusal.code.starts_with("field.preview-wgsl") => {
                return Ok(Self {
                    plan,
                    kernels: None,
                    hardware: processor(format!("{}: {}", refusal.code, refusal.detail)),
                    options,
                    true_distance_formula,
                });
            }
            Err(refusal) => return Err(refusal),
        };
        let source = format!("{field}\n{}", include_str!("batch.wgsl"));
        let kernels = Kernels {
            gpu,
            points: gpu.pipeline(&source, "evaluate_points")?,
            grid: gpu.pipeline(&source, "evaluate_grid")?,
        };
        Ok(Self {
            plan,
            kernels: Some(kernels),
            hardware: Hardware::GraphicsCard {
                adapter: gpu.name().to_string(),
            },
            options,
            true_distance_formula,
        })
    }

    /// The route every batch from this evaluator takes.
    pub fn hardware(&self) -> &Hardware {
        &self.hardware
    }

    fn samples(&self, values: Vec<f64>) -> Samples {
        Samples {
            values,
            true_distance_formula: self.true_distance_formula,
            hardware: self.hardware.clone(),
        }
    }

    /// Evaluates arbitrary points in input order. Non-finite coordinates
    /// refuse with `field.non-finite-point` on both routes. On the graphics
    /// card each point is rounded to binary32 first; a coordinate outside the
    /// finite binary32 range refuses with `preview.gpu-range`. Any non-finite
    /// result refuses the whole batch with `preview.gpu-non-finite`.
    pub fn points(&self, points: &[[f64; 3]]) -> Result<Samples, Refusal> {
        if let Some(i) = points.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
            return Err(Refusal::new(
                "field.non-finite-point",
                format!("point {i} has a non-finite coordinate"),
            ));
        }
        let Some(kernels) = &self.kernels else {
            let mut storage = UnprovedBatch::default();
            let values = self
                .plan
                .unproved_preview_batch(points, self.options.workers, &mut storage)?
                .iter()
                .map(|v| v.value)
                .collect();
            return Ok(self.samples(values));
        };
        let mut packed = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            let q = p.map(|c| c as f32);
            if q.iter().any(|c| !c.is_finite()) {
                return Err(Refusal::new(
                    "preview.gpu-range",
                    format!("point {i} is outside the finite binary32 range"),
                ));
            }
            packed.push([q[0], q[1], q[2], 0.0]);
        }
        let values = self.dispatch(kernels, &kernels.points, points.len(), |start, len| {
            (
                Params {
                    origin: [0; 4],
                    counts: [0; 4],
                    start: start as u32,
                    len: len as u32,
                    pitch: 0.0,
                    pad: 0,
                },
                Some(bytemuck::cast_slice(&packed[start..start + len])),
            )
        })?;
        Ok(self.samples(values))
    }

    /// Evaluates every point of `grid` in index order.
    pub fn grid(&self, grid: &Grid) -> Result<Samples, Refusal> {
        grid.check()?;
        let Some(kernels) = &self.kernels else {
            let mut storage = UnprovedBatch::default();
            let mut values = Vec::with_capacity(grid.len());
            let mut points = Vec::with_capacity(PROCESSOR_CHUNK.min(grid.len()));
            let mut start = 0;
            while start < grid.len() {
                let end = (start + PROCESSOR_CHUNK).min(grid.len());
                points.clear();
                points.extend((start..end).map(|i| grid.point(i)));
                values.extend(
                    self.plan
                        .unproved_preview_batch(&points, self.options.workers, &mut storage)?
                        .iter()
                        .map(|v| v.value),
                );
                start = end;
            }
            return Ok(self.samples(values));
        };
        let origin = [grid.origin[0], grid.origin[1], grid.origin[2], 0];
        let counts = [grid.counts[0], grid.counts[1], grid.counts[2], 0];
        let pitch = grid.pitch() as f32;
        let values = self.dispatch(kernels, &kernels.grid, grid.len(), |start, len| {
            (
                Params {
                    origin,
                    counts,
                    start: start as u32,
                    len: len as u32,
                    pitch,
                    pad: 0,
                },
                None,
            )
        })?;
        Ok(self.samples(values))
    }

    /// Runs `pipeline` over `total` values in dispatches of at most
    /// `options.dispatch`, in order, and checks every result is finite.
    fn dispatch<'a>(
        &self,
        kernels: &Kernels<'_>,
        pipeline: &wgpu::ComputePipeline,
        total: usize,
        chunk: impl Fn(usize, usize) -> (Params, Option<&'a [u8]>),
    ) -> Result<Vec<f64>, Refusal> {
        let gpu = kernels.gpu;
        let mut values = Vec::with_capacity(total);
        if total == 0 {
            return Ok(values);
        }
        let size = self.options.dispatch.min(total);
        let output_bytes = (size * 4) as u64;
        let params = gpu.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unproved batch parameters"),
            size: std::mem::size_of::<Params>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let output = gpu.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unproved batch values"),
            size: output_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let readback = gpu.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unproved batch readback"),
            size: output_bytes,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let uses_points = chunk(0, size).1.is_some();
        let input = uses_points.then(|| {
            gpu.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("unproved batch points"),
                size: (size * 16) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        });
        let mut entries = vec![
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: output.as_entire_binding(),
            },
        ];
        if let Some(input) = &input {
            entries.push(wgpu::BindGroupEntry {
                binding: 2,
                resource: input.as_entire_binding(),
            });
        }
        let group = gpu.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("unproved batch"),
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        });
        let mut start = 0;
        while start < total {
            let len = size.min(total - start);
            let (header, points) = chunk(start, len);
            gpu.queue
                .write_buffer(&params, 0, bytemuck::bytes_of(&header));
            if let (Some(input), Some(points)) = (&input, points) {
                gpu.queue.write_buffer(input, 0, points);
            }
            let mut encoder = gpu.device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(pipeline);
                pass.set_bind_group(0, &group, &[]);
                pass.dispatch_workgroups(len.div_ceil(64) as u32, 1, 1);
            }
            let bytes = (len * 4) as u64;
            encoder.copy_buffer_to_buffer(&output, 0, &readback, 0, bytes);
            gpu.queue.submit([encoder.finish()]);
            let chunk_values: Vec<f32> = gpu.read(&readback, bytes)?;
            for (offset, v) in chunk_values.into_iter().enumerate() {
                if !v.is_finite() {
                    return Err(Refusal::new(
                        "preview.gpu-non-finite",
                        format!(
                            "graphics card value {} is non-finite; no partial batch is kept",
                            start + offset
                        ),
                    ));
                }
                // Canonical zero, as on the processor route.
                values.push(if v == 0.0 { 0.0 } else { f64::from(v) });
            }
            start += len;
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_recipe::Recipe;

    const P1: &str = include_str!("../../../parts/P1-bracket.json");
    const P3: &str = include_str!("../../../parts/P3-filleted-boss.json");

    fn plan(text: &str) -> Plan {
        Plan::compile(&Recipe::parse(text).unwrap()).unwrap()
    }

    fn options(workers: usize, dispatch: usize) -> Options {
        Options { workers, dispatch }
    }

    fn small_grid() -> Grid {
        Grid::covering([-1.0; 3], [41.0, 31.0, 23.0], 0, 1).unwrap()
    }

    #[test]
    fn batch_kernels_validate_without_a_gpu() {
        for text in [P1, P3] {
            let source = format!(
                "{}\n{}",
                plan(text).unproved_preview_wgsl().unwrap(),
                include_str!("batch.wgsl")
            );
            let module = naga::front::wgsl::parse_str(&source)
                .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
            naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::empty(),
            )
            .validate(&module)
            .unwrap();
        }
    }

    #[test]
    fn oversized_axis_products_refuse_instead_of_wrapping() {
        // 4,194,304 points per axis: 2^66 in all, beyond usize.
        let refusal = Grid::covering([-2_097_150.0; 3], [2_097_150.0; 3], 0, 2).unwrap_err();
        assert_eq!(refusal.code, "preview.grid-range");
        let grid = Grid {
            origin: [-2_097_152; 3],
            counts: [4_194_304; 3],
            pitch_exponent: 0,
        };
        assert_eq!(grid.checked_len(), None);
        assert_eq!(grid.len(), usize::MAX);
        assert_eq!(
            crate::first_paint::check_size(&grid, usize::MAX)
                .unwrap_err()
                .code,
            "preview.grid-size"
        );
        let plan = plan(P3);
        let evaluator = Evaluator::new(&plan, None, options(1, MAX_DISPATCH)).unwrap();
        assert_eq!(
            evaluator.grid(&grid).unwrap_err().code,
            "preview.grid-range"
        );
    }

    #[test]
    fn grid_points_are_exact_in_binary32() {
        let grid = Grid::covering([-60.3, 0.0, 1e-3], [60.0, 40.0, 30.0], -3, 2).unwrap();
        assert_eq!(grid.origin, [-485, -2, -2]);
        for index in [0, 1, grid.len() / 2, grid.len() - 1] {
            let p = grid.point(index);
            assert!(
                p.iter()
                    .all(|c| f64::from(*c as f32).to_bits() == c.to_bits())
            );
        }
        assert_eq!(
            grid.point(grid.len() - 1),
            grid.at(grid.counts[0] - 1, grid.counts[1] - 1, grid.counts[2] - 1)
        );
        for refused in [
            Grid::covering([0.0; 3], [1.0; 3], 41, 0),
            Grid::covering([0.0; 3], [1e9; 3], -1, 0),
            Grid::covering([f64::NAN; 3], [1.0; 3], 0, 0),
        ] {
            assert_eq!(refused.unwrap_err().code, "preview.grid-range");
        }
    }

    #[test]
    fn processor_route_matches_the_scalar_preview_at_any_worker_count() {
        let plan = plan(P3);
        let grid = small_grid();
        let reference = Evaluator::new(&plan, None, options(1, MAX_DISPATCH))
            .unwrap()
            .grid(&grid)
            .unwrap();
        let mut w = Workspace::default();
        for (i, v) in reference.values.iter().enumerate().step_by(97) {
            let scalar = plan.unproved_preview(grid.point(i), &mut w).unwrap().value;
            assert_eq!(v.to_bits(), scalar.to_bits());
        }
        assert!(matches!(reference.hardware, Hardware::Processor { .. }));
        assert_eq!(Samples::STATUS, "unproved");
        assert!(reference.values.iter().any(|v| *v < 0.0));
        for workers in [2, 4] {
            let values = Evaluator::new(&plan, None, options(workers, MAX_DISPATCH))
                .unwrap()
                .grid(&grid)
                .unwrap()
                .values;
            assert_eq!(
                values.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                reference
                    .values
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>()
            );
        }
        let points: Vec<_> = (0..grid.len()).map(|i| grid.point(i)).collect();
        assert_eq!(
            Evaluator::new(&plan, None, options(3, MAX_DISPATCH))
                .unwrap()
                .points(&points)
                .unwrap()
                .values,
            reference.values
        );
    }

    #[test]
    fn invalid_requests_refuse_with_their_codes() {
        let plan = plan(P3);
        let refused = |o| Evaluator::new(&plan, None, o).err().unwrap().code;
        assert_eq!(refused(options(0, MAX_DISPATCH)), "field.preview-workers");
        assert_eq!(refused(options(1, 0)), "preview.dispatch");
        assert_eq!(refused(options(1, MAX_DISPATCH + 1)), "preview.dispatch");
        let evaluator = Evaluator::new(&plan, None, options(2, MAX_DISPATCH)).unwrap();
        assert_eq!(
            evaluator
                .points(&[[0.0, f64::INFINITY, 0.0]])
                .unwrap_err()
                .code,
            "field.non-finite-point"
        );
        assert!(evaluator.points(&[]).unwrap().values.is_empty());
    }

    #[test]
    fn gpu_route_repeats_bit_for_bit_and_agrees_in_sign_when_available() {
        let Ok(gpu) = Gpu::new(false) else {
            eprintln!("GPU check skipped: no adapter; kernel validation still runs");
            return;
        };
        let plan = plan(P3);
        let grid = small_grid();
        let device = Evaluator::new(&plan, Some(&gpu), options(4, MAX_DISPATCH)).unwrap();
        assert!(matches!(device.hardware(), Hardware::GraphicsCard { .. }));
        let first = device.grid(&grid).unwrap();
        let bits = |s: &Samples| s.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&first), bits(&device.grid(&grid).unwrap()));
        let split = Evaluator::new(&plan, Some(&gpu), options(4, 1_001)).unwrap();
        assert_eq!(bits(&first), bits(&split.grid(&grid).unwrap()));
        let points: Vec<_> = (0..grid.len()).map(|i| grid.point(i)).collect();
        assert_eq!(bits(&first), bits(&split.points(&points).unwrap()));
        let cpu = Evaluator::new(&plan, None, options(4, MAX_DISPATCH))
            .unwrap()
            .grid(&grid)
            .unwrap();
        for (i, (g, c)) in first.values.iter().zip(&cpu.values).enumerate() {
            assert!((g - c).abs() < 1e-3, "point {i}: {g} vs {c}");
            if c.abs() > 1e-4 {
                assert_eq!(first.side(i), cpu.side(i), "point {i}");
            }
        }
        assert_eq!(
            device.points(&[[1e39, 0.0, 0.0]]).unwrap_err().code,
            "preview.gpu-range"
        );
    }
}
