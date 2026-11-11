//! Throwaway unproved field viewer. No export or proven measurement API.
use catallax_field::{Plan, Workspace};
mod ownership;
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{error::Error, fs, path::Path, time::Instant};
use wgpu::util::DeviceExt;

const WIDTH: usize = 960;
const HEIGHT: usize = 720;
const SCALE: f64 = 0.065;
const FORWARD: [f64; 3] = [-0.577350269; 3];
const RIGHT: [f64; 3] = [
    std::f64::consts::FRAC_1_SQRT_2,
    -std::f64::consts::FRAC_1_SQRT_2,
    0.0,
];
const UP: [f64; 3] = [-0.408248290, -0.408248290, 0.816496581];
const TARGET: [f64; 3] = [20.0, 15.0, 10.0];
const ALLOWANCE: u64 = 1_000_000;
// The existing F8 nearest-face regression allowance, per exact query.
const OWNERSHIP_ALLOWANCE: u64 = 100_000_000;
type Result<T> = std::result::Result<T, Box<dyn Error>>;

struct Gpu {
    device: wgpu::Device,
    queue: wgpu::Queue,
    name: String,
}
impl Gpu {
    fn new(disabled: bool) -> std::result::Result<Self, String> {
        if disabled {
            return Err("preview.gpu-disabled".into());
        }
        pollster::block_on(async {
            let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    ..Default::default()
                })
                .await
                .map_err(|e| format!("preview.no-adapter: {e}"))?;
            let info = adapter.get_info();
            let (device, queue) = adapter
                .request_device(&wgpu::DeviceDescriptor::default())
                .await
                .map_err(|e| format!("preview.no-device: {e}"))?;
            Ok(Self {
                device,
                queue,
                name: format!("{} {:?} {:?}", info.name, info.backend, info.device_type),
            })
        })
    }
    fn dispatch<T: bytemuck::Pod>(
        &self,
        source: &str,
        entry: &str,
        data: &[&[u8]],
        output: usize,
        groups: [u32; 3],
    ) -> Result<Vec<T>> {
        self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("unproved field preview"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let pipeline = self
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("unproved preview"),
                layout: None,
                module: &shader,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            });
        if let Some(error) = pollster::block_on(self.device.pop_error_scope()) {
            return Err(Refusal::new("preview.shader", error.to_string()).into());
        }
        let buffers: Vec<_> = data
            .iter()
            .map(|bytes| {
                self.device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some("unproved preview buffer"),
                        contents: bytes,
                        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                    })
            })
            .collect();
        let buffer = &buffers[output];
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unproved preview readback"),
            size: buffer.size(),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &buffers
                .iter()
                .enumerate()
                .map(|(i, buffer)| wgpu::BindGroupEntry {
                    binding: i as u32,
                    resource: buffer.as_entire_binding(),
                })
                .collect::<Vec<_>>(),
        });
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &group, &[]);
            pass.dispatch_workgroups(groups[0], groups[1], groups[2]);
        }
        encoder.copy_buffer_to_buffer(buffer, 0, &readback, 0, buffer.size());
        self.queue.submit([encoder.finish()]);
        let (tx, rx) = std::sync::mpsc::channel();
        readback.slice(..).map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: Some(std::time::Duration::from_secs(120)),
        })?;
        rx.recv()??;
        let mapped = readback.slice(..).get_mapped_range();
        let result = bytemuck::cast_slice(&mapped).to_vec();
        drop(mapped);
        readback.unmap();
        Ok(result)
    }
    fn compute(&self, source: &str, entry: &str, data: &[[f32; 4]]) -> Result<Vec<[f32; 4]>> {
        let groups = if entry == "render" {
            [(WIDTH as u32).div_ceil(8), (HEIGHT as u32).div_ceil(8), 1]
        } else {
            [(data.len() as u32).div_ceil(64), 1, 1]
        };
        self.dispatch(source, entry, &[bytemuck::cast_slice(data)], 0, groups)
    }
    fn owners(&self, source: &str, pixels: &[[f32; 4]]) -> Result<Vec<u32>> {
        let ids = vec![0u32; pixels.len()];
        self.dispatch(
            source,
            "select_faces",
            &[bytemuck::cast_slice(pixels), bytemuck::cast_slice(&ids)],
            1,
            [(pixels.len() as u32).div_ceil(64), 1, 1],
        )
    }
}

fn origin(x: usize, y: usize) -> [f64; 3] {
    std::array::from_fn(|i| {
        TARGET[i] - FORWARD[i] * 80.0
            + RIGHT[i] * (x as f64 + 0.5 - WIDTH as f64 * 0.5) * SCALE
            + UP[i] * (HEIGHT as f64 * 0.5 - y as f64 - 0.5) * SCALE
    })
}
fn cpu_pixel(
    plan: &Plan,
    x: usize,
    y: usize,
    workspace: &mut Workspace,
) -> std::result::Result<[f32; 4], Refusal> {
    let o = origin(x, y);
    let at = |t: f64| std::array::from_fn(|i| o[i] + t * FORWARD[i]);
    let bounds = [40.1, 30.1, 22.1];
    let mut t = 0.0f64;
    let mut end = 140.0f64;
    for i in 0..3 {
        let a = (-0.1 - o[i]) / FORWARD[i];
        let b = (bounds[i] - o[i]) / FORWARD[i];
        t = t.max(a.min(b));
        end = end.min(a.max(b));
    }
    let mut previous = t;
    for _ in 0..16384 {
        if t > end {
            return Ok([0.0; 4]);
        }
        let d = plan.unproved_preview(at(t), workspace)?.value;
        if d < 0.0 {
            let (mut l, mut h) = (previous, t);
            for _ in 0..16 {
                let m = (l + h) * 0.5;
                if plan.unproved_preview(at(m), workspace)?.value < 0.0 {
                    h = m;
                } else {
                    l = m;
                }
            }
            let p = at((l + h) * 0.5);
            return Ok([p[0] as f32, p[1] as f32, p[2] as f32, 1.0]);
        }
        previous = t;
        t += (d * 0.8).max(0.002);
    }
    Err(Refusal::new(
        "preview.ray-budget",
        "ray step allowance exhausted",
    ))
}
fn cpu_render(plan: &Plan, workers: usize) -> Result<Vec<[f32; 4]>> {
    let mut pixels = vec![[0.0; 4]; WIDTH * HEIGHT];
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (chunk, data) in pixels
            .chunks_mut((WIDTH * HEIGHT).div_ceil(workers))
            .enumerate()
        {
            let start = chunk * (WIDTH * HEIGHT).div_ceil(workers);
            handles.push(
                std::thread::Builder::new()
                    .spawn_scoped(scope, move || {
                        let mut w = Workspace::default();
                        for (j, p) in data.iter_mut().enumerate() {
                            let k = start + j;
                            *p = cpu_pixel(plan, k % WIDTH, k / WIDTH, &mut w)?;
                        }
                        Ok::<_, Refusal>(())
                    })
                    .map_err(|e| Refusal::new("preview.worker-start", e.to_string()))?,
            );
        }
        for h in handles {
            h.join()
                .map_err(|_| Refusal::new("preview.worker", "preview worker panicked"))??;
        }
        Ok::<_, Refusal>(())
    })?;
    Ok(pixels)
}

fn reference(recipe: &Recipe) -> Result<Vec<bool>> {
    let (answer, work) = with_operation_budget(ALLOWANCE, || {
        catallax_features::evaluate(recipe, &mut Decider::default())
    });
    let built = answer?;
    if work.exhausted || !built.model.check(built.result).is_empty() {
        return Err(Refusal::new("preview.reference-brep", "B-rep reference refused").into());
    }
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &catallax_mesh::MeshOptions::default(),
        &mut Decider::default(),
    )?;
    if !mesh.marks.is_empty() {
        return Err(Refusal::new("preview.reference-mesh", "marked B-rep reference mesh").into());
    }
    let mut mask = vec![false; WIDTH * HEIGHT];
    let project = |p: [f64; 3]| {
        let d = std::array::from_fn::<_, 3, _>(|i| p[i] - TARGET[i]);
        [
            WIDTH as f64 * 0.5 + dot(d, RIGHT) / SCALE,
            HEIGHT as f64 * 0.5 - dot(d, UP) / SCALE,
        ]
    };
    for triangle in &mesh.value.triangles {
        let p = triangle.vertices.map(|v| project(mesh.value.positions[v]));
        let minx = p
            .iter()
            .map(|v| v[0])
            .fold(f64::INFINITY, f64::min)
            .floor()
            .max(0.0) as usize;
        let maxx = p
            .iter()
            .map(|v| v[0])
            .fold(0.0, f64::max)
            .ceil()
            .min(WIDTH as f64 - 1.0) as usize;
        let miny = p
            .iter()
            .map(|v| v[1])
            .fold(f64::INFINITY, f64::min)
            .floor()
            .max(0.0) as usize;
        let maxy = p
            .iter()
            .map(|v| v[1])
            .fold(0.0, f64::max)
            .ceil()
            .min(HEIGHT as f64 - 1.0) as usize;
        for y in miny..=maxy {
            for x in minx..=maxx {
                let mut signs = [0.0; 3];
                for i in 0..3 {
                    let a = p[i];
                    let b = p[(i + 1) % 3];
                    signs[i] = (b[0] - a[0]) * (y as f64 + 0.5 - a[1])
                        - (b[1] - a[1]) * (x as f64 + 0.5 - a[0]);
                }
                if signs.iter().all(|v| *v >= 0.0) || signs.iter().all(|v| *v <= 0.0) {
                    mask[y * WIDTH + x] = true;
                }
            }
        }
    }
    Ok(mask)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn boundary(mask: &[bool]) -> Vec<(i32, i32)> {
    (1..HEIGHT - 1)
        .flat_map(|y| {
            (1..WIDTH - 1).filter_map(move |x| {
                let i = y * WIDTH + x;
                (mask[i]
                    && [i - 1, i + 1, i - WIDTH, i + WIDTH]
                        .iter()
                        .any(|j| !mask[*j]))
                .then_some((x as i32, y as i32))
            })
        })
        .collect()
}
fn outline_error(a: &[bool], b: &[bool]) -> Result<f64> {
    let a = boundary(a);
    let b = boundary(b);
    if a.is_empty() || b.is_empty() {
        return Err(Refusal::new(
            "preview.empty-outline",
            "empty reference or preview outline",
        )
        .into());
    }
    let directed = |a: &[(i32, i32)], b: &[(i32, i32)]| {
        a.iter()
            .map(|(x, y)| {
                b.iter()
                    .map(|(u, v)| ((x - u).pow(2) + (y - v).pow(2)) as f64)
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
            .sqrt()
    };
    Ok(directed(&a, &b).max(directed(&b, &a)))
}
fn ppm(path: &Path, ids: &[u32]) -> Result<()> {
    let mut data =
        format!("P6\n# unproved preview, never an export\n{WIDTH} {HEIGHT}\n255\n").into_bytes();
    for &id in ids {
        data.extend_from_slice(&match id {
            0 => [24, 28, 34],
            u32::MAX => [255, 48, 80],
            _ => [
                70 + (id.wrapping_mul(67) % 160) as u8,
                70 + (id.wrapping_mul(97) % 160) as u8,
                70 + (id.wrapping_mul(43) % 160) as u8,
            ],
        });
    }
    let header = data.len() - WIDTH * HEIGHT * 3;
    for (letter, ch) in "UNPROVED PREVIEW".chars().enumerate() {
        let rows = match ch {
            'U' => [17, 17, 17, 17, 17, 17, 14],
            'N' => [17, 25, 25, 21, 19, 19, 17],
            'P' => [30, 17, 17, 30, 16, 16, 16],
            'R' => [30, 17, 17, 30, 20, 18, 17],
            'O' => [14, 17, 17, 17, 17, 17, 14],
            'V' => [17, 17, 17, 17, 17, 10, 4],
            'E' => [31, 16, 16, 30, 16, 16, 31],
            'D' => [30, 17, 17, 17, 17, 17, 30],
            'I' => [31, 4, 4, 4, 4, 4, 31],
            'W' => [17, 17, 17, 21, 21, 21, 10],
            _ => [0; 7],
        };
        for (y, row) in rows.iter().enumerate() {
            for x in 0..5 {
                if row & (1 << (4 - x)) == 0 {
                    continue;
                }
                for dy in 0..2 {
                    for dx in 0..2 {
                        let pixel = (16 + y * 2 + dy) * WIDTH + 16 + letter * 12 + x * 2 + dx;
                        data[header + pixel * 3..header + pixel * 3 + 3]
                            .copy_from_slice(&[245, 205, 110]);
                    }
                }
            }
        }
    }
    fs::write(path, data)?;
    Ok(())
}
fn probes(plan: &Plan, gpu: &Gpu, source: &str) -> Result<Value> {
    let mut state = 7u64;
    let points: Vec<[f32; 4]> = (0..20000)
        .map(|_| {
            let mut p = [0.0; 4];
            for i in 0..3 {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                p[i] = ([40.0, 30.0, 22.0][i]
                    * (1.2 * ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.1))
                    as f32;
            }
            p
        })
        .collect();
    let values = gpu.compute(source, "probe", &points)?;
    let mut max_error = 0.0f64;
    let mut sign_errors = 0;
    let mut w = Workspace::default();
    for p in values {
        let reference = plan
            .unproved_preview([p[0] as f64, p[1] as f64, p[2] as f64], &mut w)?
            .value;
        if !p[3].is_finite() {
            return Err(Refusal::new("preview.non-finite", "GPU probe is non-finite").into());
        }
        max_error = max_error.max((p[3] as f64 - reference).abs());
        if reference.abs() > 1e-4 && (p[3] < 0.0) != (reference < 0.0) {
            sign_errors += 1;
        }
    }
    if max_error > 0.001 || sign_errors != 0 {
        return Err(Refusal::new(
            "preview.evaluator-mismatch",
            format!("max error {max_error}, sign mismatches {sign_errors}"),
        )
        .into());
    }
    Ok(
        json!({"points":points.len(),"max_absolute_difference_mm":max_error,"sign_mismatches_beyond_0_0001_mm":sign_errors}),
    )
}
fn run() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() {
        return Err("usage: catallax-preview OUTPUT_DIRECTORY [--cpu] [--workers N]".into());
    }
    let workers = if let Some(i) = args.iter().position(|s| s == "--workers") {
        args.get(i + 1)
            .ok_or("missing worker count")?
            .parse::<usize>()?
    } else {
        std::thread::available_parallelism().map_or(1, usize::from)
    };
    if workers == 0 || workers > 1024 {
        return Err(Refusal::new("preview.workers", "workers must be between 1 and 1024").into());
    }
    let out = Path::new(&args[0]);
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()?;
    fs::create_dir_all(out)?;
    if out.canonicalize()?.starts_with(&root) {
        return Err(Refusal::new(
            "preview.artifacts",
            "review artifacts must be outside the repository",
        )
        .into());
    }
    let gpu = Gpu::new(args.iter().any(|s| s == "--cpu"));
    let mut report = json!({"status":"unproved preview","timing_context":"shared machine, information only","export":false,"proven_measurement":false,"width":WIDTH,"height":HEIGHT,"workers":workers,"backend":gpu.as_ref().map(|g|g.name.clone()).unwrap_or_else(|e|format!("processor fallback: {e}")),"frames":[]});
    let mut frames = Vec::new();
    for (name, radius) in [("p3", "2"), ("p3-radius-edit", "3")] {
        let mut definition: Value =
            serde_json::from_str(include_str!("../../../parts/P3-filleted-boss.json"))?;
        let junction = definition["steps"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|step| step["id"] == "junction")
            .unwrap();
        junction["radius"] = json!(radius);
        let recipe = Recipe::parse(&definition.to_string())?;
        let plan = Plan::compile(&recipe)?;
        let source = format!(
            "{}\n{}",
            plan.unproved_preview_wgsl()?,
            include_str!("ray.wgsl")
        );
        let comparison = if let Ok(g) = &gpu {
            probes(&plan, g, &source)?
        } else {
            Value::Null
        };
        let start = Instant::now();
        let pixels = if let Ok(g) = &gpu {
            g.compute(&source, "render", &vec![[0.0; 4]; WIDTH * HEIGHT])?
        } else {
            cpu_render(&plan, workers)?
        };
        let frame_ms = start.elapsed().as_secs_f64() * 1000.0;
        let unfinished = pixels.iter().filter(|p| p[3] < 0.0).count();
        if unfinished > 0 {
            return Err(Refusal::new(
                "preview.ray-budget",
                format!("{unfinished} GPU rays exhausted 16384 steps"),
            )
            .into());
        }
        let reference = reference(&recipe)?;
        let mask: Vec<_> = pixels
            .iter()
            .map(|p| p[3].to_bits() == 1.0f32.to_bits())
            .collect();
        let error = outline_error(&mask, &reference)?;
        let named = ownership::named_faces(&plan)?;
        let radius_number = radius.parse::<f64>()?;
        let start = Instant::now();
        let ids = if let Ok(g) = &gpu {
            g.owners(&ownership::source(&named, radius_number), &pixels)?
        } else {
            ownership::cpu(&pixels, &named, radius_number, workers)?
        };
        let ownership_ms = start.elapsed().as_secs_f64() * 1000.0;
        if pixels
            .iter()
            .zip(&ids)
            .any(|(p, id)| p[3] > 0.0 && *id == 0)
        {
            return Err(Refusal::new(
                "preview.missing-owner",
                "surface pixel has no named preview owner",
            )
            .into());
        }
        eprintln!(
            "{name}: ray frame {frame_ms:.1} ms, ownership {ownership_ms:.1} ms, outline {error:.3} px"
        );
        ppm(&out.join(format!("{name}.ppm")), &ids)?;
        ppm(
            &out.join(format!("{name}-brep.ppm")),
            &reference.iter().map(|b| u32::from(*b)).collect::<Vec<_>>(),
        )?;
        fs::write(
            out.join(format!("{name}.face-ids.u32le")),
            ids.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>(),
        )?;
        fs::write(
            out.join(format!("{name}.hits.f32le")),
            pixels
                .iter()
                .flatten()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        )?;
        let names: Vec<_> = plan
            .faces()
            .iter()
            .enumerate()
            .filter(|(i, _)| plan.result_face(*i).is_some())
            .map(|(i, f)| json!({"id":i+1,"origin":&f.origin}))
            .collect();
        fs::write(
            out.join(format!("{name}.names.json")),
            serde_json::to_vec_pretty(&json!({"background":0,"refused":u32::MAX,"faces":names}))?,
        )?;
        let row = json!({"name":name,"junction_radius_mm":radius,"outline_error_px":error,"frame_and_readback_ms":frame_ms,"ownership_ms":ownership_ms,"ownership":{"status":"unproved","method":"nearest closed P3 patch","refused_pixels":0},"f8_check":{"status":"pending","sample_count":ownership::SAMPLE_COUNT},"gpu_processor_probe":comparison,"hit_pixels":mask.iter().filter(|b|**b).count()});
        frames.push((name.to_string(), plan, pixels, ids));
        report["frames"].as_array_mut().unwrap().push(row);
        fs::write(out.join("report.json"), serde_json::to_vec_pretty(&report)?)?;
        if error > 1.0 {
            return Err(Refusal::new(
                "preview.outline",
                format!("{name} outline error {error} pixels exceeds one pixel"),
            )
            .into());
        }
    }
    let mut wrong = 0;
    for (i, (name, plan, pixels, ids)) in frames.iter().enumerate() {
        let samples = ownership::sample_pixels(pixels, ownership::SAMPLE_COUNT)?;
        let prepared = plan.prepare_distances(&mut Decider::default())?;
        prepared.prepare_nearest(&mut Decider::default())?;
        let start = Instant::now();
        let mut check = ownership::check(plan, &prepared, pixels, ids, &samples, workers)?;
        fs::write(
            out.join(format!("{name}.f8-samples.json")),
            serde_json::to_vec_pretty(&check)?,
        )?;
        let sample_ms = start.elapsed().as_secs_f64() * 1000.0;
        check.as_object_mut().unwrap().remove("samples");
        check["status"] = json!("checked");
        wrong += check["proved"].as_u64().unwrap() - check["matching_proved"].as_u64().unwrap();
        eprintln!(
            "{name}: F8 sample agreement {}/{}; refused {}",
            check["matching_proved"], check["sample_count"], check["refused"]
        );
        report["frames"][i]["f8_check"] = check;
        report["frames"][i]["f8_sample_ms"] = json!(sample_ms);
        fs::write(out.join("report.json"), serde_json::to_vec_pretty(&report)?)?;
    }
    if wrong > 0 {
        return Err(Refusal::new("preview.ownership-mismatch",format!("{wrong} GPU or processor preview IDs disagree with proved F8 answers; all pixels and both answers retained in sample reports")).into());
    }
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn outline_metric_rejects_two_pixel_displacement() {
        let mut a = vec![false; WIDTH * HEIGHT];
        let mut b = a.clone();
        for y in 100..200 {
            for x in 100..200 {
                a[y * WIDTH + x] = true;
                b[y * WIDTH + x + 2] = true;
            }
        }
        assert_eq!(outline_error(&a, &a).unwrap(), 0.0);
        assert_eq!(outline_error(&a, &b).unwrap(), 2.0);
    }
}

#[cfg(test)]
mod gpu_tests {
    use super::*;
    const P1: &str = include_str!("../../../parts/P1-bracket.json");
    const P3: &str = include_str!("../../../parts/P3-filleted-boss.json");
    fn shader(plan: &Plan) -> String {
        format!(
            "{}\n{}",
            plan.unproved_preview_wgsl().unwrap(),
            include_str!("ray.wgsl")
        )
    }
    #[test]
    fn generated_wgsl_validates_without_a_gpu() {
        for text in [
            P1.to_string(),
            P3.to_string(),
            P3.replace("\"radius\": \"2\"", "\"radius\": \"3\""),
        ] {
            let plan = Plan::compile(&Recipe::parse(&text).unwrap()).unwrap();
            let source = shader(&plan);
            let module = naga::front::wgsl::parse_str(&source)
                .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
            naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::empty(),
            )
            .validate(&module)
            .unwrap();
            assert_eq!(source, shader(&plan));
        }
    }
    #[test]
    fn gpu_values_match_processor_when_adapter_is_available() {
        let Ok(gpu) = Gpu::new(false) else {
            eprintln!("GPU check skipped: no adapter; WGSL validation still runs");
            return;
        };
        for text in [P1, P3] {
            let plan = Plan::compile(&Recipe::parse(text).unwrap()).unwrap();
            probes(&plan, &gpu, &shader(&plan)).unwrap();
        }
    }
    #[test]
    fn cpu_ray_matches_with_one_two_and_four_workers() {
        let plan = Plan::compile(&Recipe::parse(P3).unwrap()).unwrap();
        let points: Vec<_> = (0..48).map(|i| (100 + i * 15, 100 + i * 9)).collect();
        let render = |workers: usize| {
            std::thread::scope(|scope| {
                let handles: Vec<_> = points
                    .chunks(points.len().div_ceil(workers))
                    .map(|chunk| {
                        let plan = &plan;
                        scope.spawn(move || {
                            let mut w = Workspace::default();
                            chunk
                                .iter()
                                .map(|&(x, y)| cpu_pixel(plan, x, y, &mut w).unwrap())
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().unwrap())
                    .collect::<Vec<_>>()
            })
        };
        let expected = render(1);
        assert!(expected.iter().any(|p| p[3] > 0.0));
        assert!(expected.iter().any(|p| p[3] == 0.0));
        for workers in [2, 4] {
            assert_eq!(render(workers), expected);
        }
        assert!(matches!(Gpu::new(true),Err(e) if e=="preview.gpu-disabled"));
    }
}
