//! P5 flanged housing: the same compute and reference checks on native and WebGPU.

#[cfg(any(target_os = "linux", target_os = "macos"))]
mod corpus;

use std::collections::{BTreeSet, HashMap};

use catallax_app::{build_certified, preview_body};
use catallax_mesh::float::{DisplayDescription, DisplayFace, DisplayOptions, describe_body};
use catallax_mesh::{Mesh, MeshOptions, PreviewOptions};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_preview::mesher::{
    FaceDetail, FaceRefusal, GpuMesh, GpuRefusal, MeshBuffers, Mesher, Preview, PreviewVertex,
    preview_or_cpu,
};
use catallax_recipe::Recipe;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::{wasm_bindgen_test, wasm_bindgen_test_configure};
#[cfg(target_arch = "wasm32")]
wasm_bindgen_test_configure!(run_in_browser);

const TOLERANCE: f64 = 0.01;
const P5_DETAIL: FaceDetail = FaceDetail([128, 32]);

#[cfg(target_os = "linux")]
#[test]
fn corpus_software_compute() {
    pollster::block_on(async {
        let Some(adapter) = software_adapter().await else {
            return;
        };
        check_corpus(adapter).await;
    });
}

#[cfg(target_os = "macos")]
#[test]
#[ignore = "Metal acceptance with the held-out STEP sources"]
fn corpus_metal() {
    assert!(std::env::var_os("P5_SOURCES").is_some());
    assert!(std::env::var_os("P5_SELECTION").is_some());
    pollster::block_on(async { check_corpus(native_adapter().await).await });
}
#[cfg(target_os = "linux")]
#[test]
fn p5_software_compute() {
    pollster::block_on(async {
        let Some(adapter) = software_adapter().await else {
            return;
        };
        check_p5(adapter).await;
    });
}

#[cfg(target_os = "linux")]
async fn software_adapter() -> Option<wgpu::Adapter> {
    let adapter = adapter(wgpu::Backends::VULKAN, true).await;
    if let Some(adapter) = &adapter {
        assert_eq!(adapter.get_info().device_type, wgpu::DeviceType::Cpu);
    } else {
        assert!(
            std::env::var_os("CI").is_none(),
            "no Vulkan software adapter in CI"
        );
        eprintln!("SKIP: no Vulkan software adapter (install Mesa lavapipe)");
    }
    adapter
}

#[cfg(any(target_os = "windows", target_os = "macos"))]
#[test]
#[ignore = "native hardware acceptance; run explicitly on the GPU host"]
fn p5_native_hardware() {
    pollster::block_on(async { check_p5(native_adapter().await).await });
}

#[cfg(any(target_os = "windows", target_os = "macos"))]
async fn native_adapter() -> wgpu::Adapter {
    let (backends, backend) = if cfg!(target_os = "windows") {
        (wgpu::Backends::DX12, wgpu::Backend::Dx12)
    } else {
        (wgpu::Backends::METAL, wgpu::Backend::Metal)
    };
    let adapter = adapter(backends, false).await.expect("hardware adapter");
    let info = adapter.get_info();
    assert_eq!(info.backend, backend);
    assert!(
        matches!(
            info.device_type,
            wgpu::DeviceType::IntegratedGpu | wgpu::DeviceType::DiscreteGpu
        ),
        "software or unknown adapter: {info:?}",
    );
    adapter
}

#[cfg(target_os = "windows")]
#[test]
#[ignore = "short native DX12 hardware dispatch preflight"]
fn dx12_dispatch_smoke() {
    use catallax_mesh::float::{
        DISPLAY_DESCRIPTION_VERSION, DisplayEdge, DisplayEdgeUse, DisplayFace, DisplayFrame,
        DisplayLoop, DisplaySurface, DisplayTrimSample,
    };

    pollster::block_on(async {
        let adapter = native_adapter().await;
        eprintln!("DX12 preview adapter: {:?}", adapter.get_info());
        let (device, queue) = adapter.request_device(&Default::default()).await.unwrap();
        let mesher = Mesher::new(&device, &queue).await.unwrap();
        let corners = [[0., 0.], [2., 0.], [2., 3.], [0., 3.], [0., 0.]];
        let description = DisplayDescription {
            version: DISPLAY_DESCRIPTION_VERSION,
            origin: [1e12, -1e12, 7.],
            edges: vec![DisplayEdge {
                edge: 7,
                points: corners.map(|[x, y]| [x, y, 0.]).to_vec(),
            }],
            faces: vec![DisplayFace {
                face: 42,
                reversed: false,
                surface: DisplaySurface::Plane(DisplayFrame {
                    origin: [0.; 3],
                    x: [1., 0., 0.],
                    y: [0., 1., 0.],
                    z: [0., 0., 1.],
                }),
                bounds: [[0., 2.], [0., 3.]],
                periods: [None; 2],
                loops: vec![DisplayLoop {
                    uses: vec![DisplayEdgeUse {
                        edge: 0,
                        samples: corners
                            .into_iter()
                            .enumerate()
                            .map(|(point, uv)| DisplayTrimSample { point, uv })
                            .collect(),
                    }],
                    vertex: None,
                }],
            }],
        };
        let mesh = mesher
            .mesh(&description, |_| FaceDetail([2, 2]))
            .await
            .unwrap();
        assert_eq!(mesh.faces[0].face, 42);
        assert!(!mesh.faces[0].draw.as_ref().unwrap().sampled_grid);
        let (vertices, indices) = readback(&device, &queue, &mesh.buffers.unwrap()).await;
        assert!(
            vertices
                .iter()
                .any(|v| v.position.map(f32::to_bits) == [1., 1.5, 0., 1.].map(f32::to_bits))
        );
        assert!(
            vertices
                .iter()
                .all(|v| v.normal.map(f32::to_bits) == [0., 0., 1., 0.].map(f32::to_bits))
        );
        let (triangles, remainder) = indices.as_chunks::<3>();
        assert!(remainder.is_empty());
        let area: f32 = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| vertices[i as usize].position);
                ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) * 0.5
            })
            .sum();
        assert!((area - 6.).abs() < f32::EPSILON);
    });
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen_test]
async fn p5_browser_webgpu() {
    let adapter = adapter(wgpu::Backends::BROWSER_WEBGPU, false)
        .await
        .expect("WebGPU adapter");
    let info = adapter.get_info();
    assert_eq!(info.backend, wgpu::Backend::BrowserWebGpu);
    assert_ne!(info.device_type, wgpu::DeviceType::Cpu);
    check_p5(adapter).await;
}

async fn adapter(backends: wgpu::Backends, fallback: bool) -> Option<wgpu::Adapter> {
    let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle();
    descriptor.backends = backends;
    let instance = wgpu::Instance::new(descriptor);
    instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: fallback,
            ..Default::default()
        })
        .await
        .ok()
}

async fn check_p5(adapter: wgpu::Adapter) {
    eprintln!("P5 preview adapter: {:?}", adapter.get_info());
    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor::default())
        .await
        .unwrap();
    let gpu = Mesher::new(&device, &queue).await.unwrap();
    check_rational_patch(&gpu, &device, &queue).await;
    check_ellipse_patch(&gpu, &device, &queue).await;
    check_analytic_patches(&gpu, &device, &queue).await;
    let recipe = Recipe::parse(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    });
    let options = MeshOptions {
        chord_tolerance: Rational::from_ratio(1, 200).unwrap(),
        ..Default::default()
    };
    let app = build_certified(&recipe, Some(&options), &mut dec).unwrap();
    let checked = &app.display().unwrap().as_ref().unwrap().value;
    let built = app.built().built();
    let description =
        describe_body(&built.model, built.result, &DisplayOptions::default()).unwrap();
    let preview = preview_or_cpu(Ok(&gpu), &description, |_| P5_DETAIL, || Ok::<_, ()>(()))
        .await
        .unwrap();
    let mesh = match preview {
        Preview::Gpu(mesh) => mesh,
        Preview::Mixed { gpu, .. } => panic!(
            "P5 faces refused: {:?}",
            gpu.faces
                .iter()
                .filter(|f| f.draw.is_err())
                .collect::<Vec<_>>()
        ),
        Preview::Cpu { reason, .. } => panic!("P5 GPU refused: {reason}"),
    };
    let buffers = mesh.buffers.as_ref().unwrap();
    let (vertices, indices) = readback(&device, &queue, buffers).await;
    assert_geometry(&description, &mesh, &vertices, &indices, checked);

    let first = description.faces[0].face;
    let cpu = || preview_body(built, &PreviewOptions::default(), &mut dec);
    let fallback = preview_or_cpu(
        Ok(&gpu),
        &description,
        |face| {
            if *face == first {
                FaceDetail([0, 1])
            } else {
                FaceDetail::default()
            }
        },
        cpu,
    )
    .await
    .unwrap();
    let Preview::Mixed { gpu, cpu } = fallback else {
        panic!("one refused face must retain other GPU draws");
    };
    assert_eq!(gpu.faces[0].draw, Err(FaceRefusal::InvalidChart));
    assert!(gpu.faces[1..].iter().all(|face| face.draw.is_ok()));
    assert!(cpu.mesh.faces.iter().any(|face| face.face == first));
    let cpu_only = preview_or_cpu(
        Err(GpuRefusal::NoAdapter),
        &description,
        |_| unreachable!(),
        || Ok::<_, ()>(cpu),
    )
    .await
    .unwrap();
    assert!(
        matches!(cpu_only, Preview::Cpu { reason: GpuRefusal::NoAdapter, mesh } if mesh.mesh.faces.len() == checked.faces.len())
    );
}

async fn check_rational_patch(gpu: &Mesher, device: &wgpu::Device, queue: &wgpu::Queue) {
    use catallax_mesh::float::{DisplaySpline, DisplaySurface};

    let description = rectangular_patch(
        DisplaySurface::BSpline(DisplaySpline {
            degree: [2, 1],
            knots: [vec![8., 8., 8., 14., 14., 14.], vec![-3., -3., 1., 1.]],
            size: [3, 2],
            points: vec![
                [0., 0., 0.],
                [2., 0., 2.],
                [4., 0., 0.],
                [0., 2., 0.],
                [2., 2., 2.],
                [4., 2., 0.],
            ],
            weights: Some(vec![1., 2., 1., 1., 2., 1.]),
        }),
        [[8., 14.], [-3., 1.]],
        [[0., 0., 0.], [4., 0., 0.], [4., 2., 0.], [0., 2., 0.]],
    );
    // At (1/3, 1/2), the weighted Bernstein numerators are [4, 8, 1].
    let p = [20. / 13., 1., 16. / 13., 1.];
    let n = [-3., 0., 7., 0.].map(|x| x / 58_f32.sqrt());
    check_patch(gpu, device, queue, &description, FaceDetail([3, 2]), p, n).await;
}

async fn check_analytic_patches(gpu: &Mesher, device: &wgpu::Device, queue: &wgpu::Queue) {
    use catallax_mesh::float::{DisplayFrame, DisplaySurface};

    let frame = DisplayFrame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    };
    let d = std::f32::consts::FRAC_1_SQRT_2;
    let cases = [
        (
            DisplaySurface::Cone { frame, slope: 0.5 },
            [2., 4.],
            [[1., 0., 2.], [0., 1., 2.], [0., 2., 4.], [2., 0., 4.]],
            [1.5 * d, 1.5 * d, 3., 1.],
            [d, d, -0.5, 0.].map(|x| x / 1.25_f32.sqrt()),
        ),
        (
            DisplaySurface::Sphere {
                frame,
                radius: std::f64::consts::SQRT_2,
            },
            [-0.125, 0.125],
            [[1., 0., -1.], [0., 1., -1.], [0., 1., 1.], [1., 0., 1.]],
            [1., 1., 0., 1.],
            [d, d, 0., 0.],
        ),
        (
            DisplaySurface::Torus {
                frame,
                major: 3.,
                minor: 1.,
            },
            [-0.25, 0.25],
            [[3., 0., -1.], [0., 3., -1.], [0., 3., 1.], [3., 0., 1.]],
            [4. * d, 4. * d, 0., 1.],
            [d, d, 0., 0.],
        ),
    ];
    for (surface, v, points, p, n) in cases {
        let description = rectangular_patch(surface, [[0., 0.25], v], points);
        check_patch(gpu, device, queue, &description, FaceDetail([2, 2]), p, n).await;
    }
}

async fn check_ellipse_patch(gpu: &Mesher, device: &wgpu::Device, queue: &wgpu::Queue) {
    use catallax_mesh::float::{DisplayFrame, DisplaySurface};

    let description = rectangular_patch(
        DisplaySurface::EllipticCylinder {
            frame: DisplayFrame {
                origin: [0.; 3],
                x: [1., 0., 0.],
                y: [0., 1., 0.],
                z: [0., 0., 1.],
            },
            radii: [3., 2.],
        },
        [[0., 0.25], [0., 2.]],
        [[3., 0., 0.], [0., 2., 0.], [0., 2., 2.], [3., 0., 2.]],
    );
    let p = [3. / 2_f32.sqrt(), 2. / 2_f32.sqrt(), 1., 1.];
    let n = [2., 3., 0., 0.].map(|x| x / 13_f32.sqrt());
    check_patch(gpu, device, queue, &description, FaceDetail([2, 2]), p, n).await;
}

fn rectangular_patch(
    surface: catallax_mesh::float::DisplaySurface,
    bounds: [[f64; 2]; 2],
    points: [[f64; 3]; 4],
) -> DisplayDescription<u32, u32> {
    use catallax_mesh::float::{
        DISPLAY_DESCRIPTION_VERSION, DisplayEdge, DisplayEdgeUse, DisplayLoop, DisplayTrimSample,
    };

    let [[u0, u1], [v0, v1]] = bounds;
    let corners = [[u0, v0], [u1, v0], [u1, v1], [u0, v1], [u0, v0]];
    DisplayDescription {
        version: DISPLAY_DESCRIPTION_VERSION,
        origin: [0.; 3],
        edges: vec![DisplayEdge {
            edge: 0,
            points: [points[0], points[1], points[2], points[3], points[0]].to_vec(),
        }],
        faces: vec![DisplayFace {
            face: 0,
            reversed: false,
            surface,
            bounds,
            periods: [None; 2],
            loops: vec![DisplayLoop {
                uses: vec![DisplayEdgeUse {
                    edge: 0,
                    samples: corners
                        .into_iter()
                        .enumerate()
                        .map(|(point, uv)| DisplayTrimSample { point, uv })
                        .collect(),
                }],
                vertex: None,
            }],
        }],
    }
}

async fn check_patch(
    gpu: &Mesher,
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    description: &DisplayDescription<u32, u32>,
    detail: FaceDetail,
    p: [f32; 4],
    n: [f32; 4],
) {
    let mesh = gpu.mesh(description, |_| detail).await.unwrap();
    assert!(!mesh.faces[0].draw.as_ref().unwrap().sampled_grid);
    let (vertices, _) = readback(device, queue, &mesh.buffers.unwrap()).await;
    let vertex = vertices
        .iter()
        .find(|v| v.position.iter().zip(p).all(|(&a, b)| (a - b).abs() < 1e-5))
        .expect("missing patch interior vertex");
    assert!(
        vertex
            .normal
            .iter()
            .zip(n)
            .all(|(&a, b)| (a - b).abs() < 1e-5)
    );
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
async fn check_corpus(adapter: wgpu::Adapter) {
    use catallax_features::EvaluationOutcome;
    use catallax_mesh::{float::describe_numerical_body, mesh_body};
    use catallax_numbers::with_operation_budget;

    eprintln!("Corpus adapter: {:?}", adapter.get_info());
    let (device, queue) = adapter.request_device(&Default::default()).await.unwrap();
    let mesher = Mesher::new(&device, &queue).await.unwrap();
    let mut failures = Vec::new();
    let mut rows = 0;
    let mut grid_faces = 0;
    let selected = std::env::var("CATALLAX_DISPLAY_CASES").ok();
    let digest = std::env::var_os("CATALLAX_DISPLAY_DIGEST").is_some();
    for (id, outcome) in corpus::bodies(selected.as_deref()) {
        rows += 1;
        let started = std::time::Instant::now();
        let result: Result<(), String> = async {
            match outcome {
                EvaluationOutcome::Exact(checkpoint) => {
                    let built = checkpoint.built();
                    let description =
                        describe_body(&built.model, built.result, &corpus_display_options())
                            .map_err(|error| error.to_string())?;
                    let gpu = corpus_gpu(id, &mesher, &device, &description, started, &mut grid_faces)
                        .await;
                    let (reference, _) = with_operation_budget(u64::MAX, || {
                        mesh_body(
                            &built.model,
                            built.result,
                            &MeshOptions {
                                chord_tolerance: Rational::from_ratio(1, 200).unwrap(),
                                ..Default::default()
                            },
                            &mut Decider::default(),
                        )
                    });
                    let reference = match reference {
                        Ok(reference) => reference,
                        Err(error)
                            if matches!(id, "fillet-vertex-4-two" | "m10-top-cut")
                                && error.code == "mesh.precision-or-budget" =>
                        {
                            // These rows are reported separately from checked geometry
                            // acceptance. They must still produce a draw for every face.
                            gpu?;
                            eprintln!("{id}: known checked-mesh refusal: {error}; GPU drew all faces");
                            return Ok(());
                        }
                        Err(error) => return Err(error.to_string()),
                    };
                    let gpu = gpu?;
                    let (vertices, indices) =
                        readback(&device, &queue, gpu.buffers.as_ref().unwrap()).await;
                    if digest {
                        report_digest(id, &vertices, &indices, Some(&reference.value));
                    }
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        assert_geometry(&description, &gpu, &vertices, &indices, &reference.value);
                    }))
                    .map_err(|_| "geometry check failed".to_owned())
                }
                EvaluationOutcome::Numerical(built) => {
                    let description = describe_numerical_body(built.body(), &corpus_display_options())
                        .map_err(|error| error.to_string())?;
                    let gpu = corpus_gpu(id, &mesher, &device, &description, started, &mut grid_faces)
                        .await?;
                    let (vertices, indices) =
                        readback(&device, &queue, gpu.buffers.as_ref().unwrap()).await;
                    if digest {
                        report_digest(id, &vertices, &indices, None);
                    }
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        assert_closed(
                            &description,
                            &gpu,
                            &vertices,
                            &indices,
                            |face, face_vertices| {
                                let mut distance = face.carrier_distance().unwrap();
                                for i in face_vertices {
                                    let [x, y, z, _] = vertices[i as usize].position;
                                    let p = [x, y, z].map(f64::from);
                                    let gap = distance(p);
                                    assert!(
                                        gap <= TOLERANCE,
                                        "{id}: {:?}: vertex {p:?} is {gap:.6} mm from the admitted carrier",
                                        face.face
                                    );
                                }
                            },
                        );
                    }))
                    .map_err(|_| "geometry check failed".to_owned())
                }
            }
        }.await;
        if let Err(error) = result {
            eprintln!("{id}: {error}");
            failures.push(format!("{id}: {error}"));
        }
    }
    assert!(rows > 0, "no selected corpus rows");
    eprintln!("GPU corpus: {rows} rows, {grid_faces} sampled-grid faces");
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
fn corpus_display_options() -> DisplayOptions {
    DisplayOptions {
        chord_tolerance: TOLERANCE / 4.0,
        ..Default::default()
    }
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
fn report_digest(id: &str, vertices: &[PreviewVertex], indices: &[u32], reference: Option<&Mesh>) {
    use std::hash::{DefaultHasher, Hasher};

    let mut gpu = DefaultHasher::new();
    gpu.write(bytemuck::cast_slice(vertices));
    gpu.write(bytemuck::cast_slice(indices));
    let checked = reference.map(|reference| {
        let mut bytes = Vec::new();
        reference.write_obj(&mut bytes).unwrap();
        let mut hasher = DefaultHasher::new();
        hasher.write(&bytes);
        hasher.finish()
    });
    eprintln!(
        "{id}: gpu_digest={:016x}, checked_digest={checked:016x?}",
        gpu.finish()
    );
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
async fn corpus_gpu<F: Clone + std::fmt::Debug, E>(
    id: &str,
    mesher: &Mesher,
    device: &wgpu::Device,
    description: &DisplayDescription<F, E>,
    started: std::time::Instant,
    grid_faces: &mut usize,
) -> Result<GpuMesh<F>, String> {
    use catallax_mesh::float::DisplaySurface::*;

    let gpu = mesher
        .mesh(description, |_| P5_DETAIL)
        .await
        .map_err(|error| error.to_string())?;
    device
        .poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: Some(std::time::Duration::from_secs(30)),
        })
        .map_err(|error| error.to_string())?;
    let mut grids = std::collections::BTreeMap::<_, usize>::new();
    let mut refusals = Vec::new();
    for (face, result) in description.faces.iter().zip(&gpu.faces) {
        let kind = match face.surface {
            Plane(_) => "plane",
            Cylinder { .. } => "cylinder",
            Cone { .. } => "cone",
            Sphere { .. } => "sphere",
            Torus { .. } => "torus",
            EllipticCylinder { .. } => "elliptic-cylinder",
            BSpline(_) => "spline",
            Grid(_) => "grid",
        };
        match &result.draw {
            Ok(draw) if draw.sampled_grid => {
                *grids.entry(kind).or_default() += 1;
            }
            Ok(_) => (),
            Err(error) => {
                chart_diagnostic(face);
                refusals.push(format!("{:?} {kind}: {error}", face.face));
            }
        }
    }
    *grid_faces += grids.values().sum::<usize>();
    eprintln!(
        "{id}: preview_ms={:.3}, faces={}, sampled_grid={grids:?}",
        started.elapsed().as_secs_f64() * 1_000.0,
        gpu.faces.len()
    );
    if refusals.is_empty() {
        Ok(gpu)
    } else {
        Err(refusals.join(", "))
    }
}

fn chart_diagnostic<F: std::fmt::Debug>(face: &DisplayFace<F>) {
    if std::env::var_os("CATALLAX_TRIM_DIAGNOSTICS").is_none() {
        return;
    }
    use catallax_mesh::float::DisplaySurface::*;
    let kind = match &face.surface {
        Plane(_) => "plane",
        Cylinder { .. } => "cylinder",
        Cone { .. } => "cone",
        Sphere { .. } => "sphere",
        Torus { .. } => "torus",
        EllipticCylinder { .. } => "elliptic-cylinder",
        BSpline(_) => "spline",
        Grid(_) => "grid",
    };
    eprintln!(
        "chart {:?} {kind}: {:?}, periods {:?}, reversed {}",
        face.face, face.bounds, face.periods, face.reversed
    );
    for ring in &face.loops {
        if let Some((_, uv)) = ring.vertex {
            eprintln!("  pole {uv:?}");
        }
        for edge in &ring.uses {
            let mut lo = [f64::INFINITY; 2];
            let mut hi = [f64::NEG_INFINITY; 2];
            for sample in &edge.samples {
                for i in 0..2 {
                    lo[i] = lo[i].min(sample.uv[i]);
                    hi[i] = hi[i].max(sample.uv[i]);
                }
            }
            eprintln!(
                "  edge {}: {} samples {:?} to {:?}, extent {lo:?} through {hi:?}",
                edge.edge,
                edge.samples.len(),
                edge.samples.first().map(|p| p.uv),
                edge.samples.last().map(|p| p.uv),
            );
        }
    }
}

async fn readback(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    buffers: &MeshBuffers,
) -> (Vec<PreviewVertex>, Vec<u32>) {
    let vertex_bytes = buffers.vertex.size();
    let index_bytes = buffers.index.size();
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("P5 test readback"),
        size: vertex_bytes + index_bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder = device.create_command_encoder(&Default::default());
    encoder.copy_buffer_to_buffer(&buffers.vertex, 0, &staging, 0, vertex_bytes);
    encoder.copy_buffer_to_buffer(&buffers.index, 0, &staging, vertex_bytes, index_bytes);
    queue.submit([encoder.finish()]);
    let slice = staging.slice(..);
    let (send, recv) = futures_channel::oneshot::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = send.send(result);
    });
    #[cfg(not(target_arch = "wasm32"))]
    device
        .poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: Some(std::time::Duration::from_secs(30)),
        })
        .unwrap();
    recv.await.unwrap().unwrap();
    let mapped = slice.get_mapped_range().unwrap();
    let (vertices, indices) = mapped.split_at(vertex_bytes as usize);
    let vertices = bytemuck::cast_slice::<_, PreviewVertex>(vertices).to_vec();
    let indices = bytemuck::cast_slice::<_, u32>(indices).to_vec();
    assert_eq!(vertices.len(), buffers.vertex_count as usize);
    (vertices, indices)
}

type VertexKey = [u32; 3];
type EdgeKey = (VertexKey, VertexKey);

fn assert_geometry(
    description: &DisplayDescription,
    gpu: &GpuMesh,
    vertices: &[PreviewVertex],
    indices: &[u32],
    checked: &Mesh,
) {
    assert_eq!(gpu.faces.len(), checked.faces.len());
    assert_closed(
        description,
        gpu,
        vertices,
        indices,
        |face, face_vertices| {
            let checked_face = checked.faces.iter().find(|f| f.face == face.face).unwrap();
            let reference = DistanceGrid::new(
                &checked.triangles[checked_face.triangles.clone()],
                &checked.positions,
            );
            for index in face_vertices {
                let vertex = vertices[index as usize];
                let position =
                    std::array::from_fn(|i| f64::from(vertex.position[i]) + gpu.origin[i]);
                if !reference.within_tolerance(position) {
                    let [x, y, z, _] = vertex.position;
                    let local = [x, y, z].map(f64::from);
                    let carrier = face.carrier_distance().map(|mut distance| distance(local));
                    chart_diagnostic(face);
                    let closest = reference.distance(position);
                    panic!(
                        "face {:?}: vertex {position:?} is {closest:.6} mm from checked face; carrier {carrier:?}",
                        face.face
                    );
                }
            }
        },
    );
}

fn assert_closed<F: Eq + std::fmt::Debug, E: std::fmt::Debug>(
    description: &DisplayDescription<F, E>,
    gpu: &GpuMesh<F>,
    vertices: &[PreviewVertex],
    indices: &[u32],
    mut check_face: impl FnMut(&DisplayFace<F>, BTreeSet<u32>),
) {
    assert_eq!(
        gpu.origin.map(f64::to_bits),
        description.origin.map(f64::to_bits)
    );
    assert_eq!(gpu.faces.len(), description.faces.len());
    let key = |i: u32| {
        let p = vertices[i as usize].position;
        [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
    };
    let mut edges = HashMap::<EdgeKey, (u32, i32)>::new();
    let mut next_index = 0;
    let mut used = vec![false; vertices.len()];
    for (face, source) in gpu.faces.iter().zip(&description.faces) {
        assert_eq!(face.face, source.face);
        let draw = face.draw.as_ref().unwrap();
        assert_eq!(draw.indices.start, next_index);
        next_index = draw.indices.end;
        let draw_indices = &indices[draw.indices.start as usize..draw.indices.end as usize];
        assert!(!draw_indices.is_empty() && draw_indices.len().is_multiple_of(3));
        let face_vertices: BTreeSet<_> = draw_indices.iter().copied().collect();
        for &index in &face_vertices {
            used[index as usize] = true;
            let vertex = vertices[index as usize];
            let position: [f64; 3] =
                std::array::from_fn(|i| f64::from(vertex.position[i]) + gpu.origin[i]);
            let length_sq = vertex.normal[..3].iter().map(|x| x * x).sum::<f32>();
            assert!(
                (length_sq - 1.0).abs() < 1e-4,
                "face {:?}: bad normal {:?} at {position:?}",
                face.face,
                vertex.normal
            );
        }
        check_face(source, face_vertices);
        for t in draw_indices.as_chunks::<3>().0 {
            for i in 0..3 {
                let (a, b) = (key(t[i]), key(t[(i + 1) % 3]));
                assert_ne!(a, b, "face {:?}: collapsed triangle edge", face.face);
                let use_ = edges.entry((a.min(b), a.max(b))).or_default();
                use_.0 += 1;
                use_.1 += if a < b { 1 } else { -1 };
            }
        }
    }
    assert_eq!(next_index as usize, indices.len());
    assert!(used.into_iter().all(|used| used));
    if let Some((edge, usage)) = edges.iter().find(|&(_, &usage)| usage != (2, 0)) {
        panic!("GPU preview is not a closed, oriented mesh: {edge:?} used {usage:?}");
    }
    for edge in &description.edges {
        for pair in edge.points.windows(2) {
            let key = |p: [f64; 3]| {
                p.map(|x| {
                    let x = x as f32;
                    (if x == 0.0 { 0.0 } else { x }).to_bits()
                })
            };
            let (a, b) = (key(pair[0]), key(pair[1]));
            assert_eq!(
                edges.get(&(a.min(b), a.max(b))),
                Some(&(2, 0)),
                "missing shared chord on {:?}",
                edge.edge
            );
        }
    }
}

struct DistanceGrid {
    triangles: Vec<[[f64; 3]; 3]>,
    cells: HashMap<[i32; 3], Vec<usize>>,
}

impl DistanceGrid {
    fn new(triangles: &[catallax_mesh::Triangle], positions: &[[f64; 3]]) -> Self {
        let triangles: Vec<_> = triangles
            .iter()
            .map(|t| t.vertices.map(|v| positions[v]))
            .collect();
        let mut cells = HashMap::<_, Vec<_>>::new();
        for (i, [a, b, c]) in triangles.iter().copied().enumerate() {
            let low: [i32; 3] =
                std::array::from_fn(|j| (a[j].min(b[j]).min(c[j]) - TOLERANCE).floor() as i32);
            let high: [i32; 3] =
                std::array::from_fn(|j| (a[j].max(b[j]).max(c[j]) + TOLERANCE).floor() as i32);
            for z in low[2]..=high[2] {
                for y in low[1]..=high[1] {
                    for x in low[0]..=high[0] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        Self { triangles, cells }
    }

    fn within_tolerance(&self, p: [f64; 3]) -> bool {
        self.cells
            .get(&p.map(|x| x.floor() as i32))
            .is_some_and(|indices| {
                indices
                    .iter()
                    .any(|&i| distance_sq(p, self.triangles[i]) <= TOLERANCE * TOLERANCE)
            })
    }

    fn distance(&self, p: [f64; 3]) -> f64 {
        self.triangles
            .iter()
            .map(|&t| distance_sq(p, t))
            .fold(f64::INFINITY, f64::min)
            .sqrt()
    }
}

fn distance_sq(p: [f64; 3], [a, b, c]: [[f64; 3]; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| a[i] - b[i]);
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let (u, v, d) = (sub(b, a), sub(c, a), sub(p, a));
    let (uu, uv, vv, du, dv) = (dot(u, u), dot(u, v), dot(v, v), dot(d, u), dot(d, v));
    let det = uu * vv - uv * uv;
    if det > 0.0 {
        let s = (vv * du - uv * dv) / det;
        let t = (uu * dv - uv * du) / det;
        if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
            let offset = std::array::from_fn(|i| d[i] - s * u[i] - t * v[i]);
            return dot(offset, offset);
        }
    }
    [(a, b), (b, c), (c, a)]
        .into_iter()
        .map(|(a, b)| {
            let (u, d) = (sub(b, a), sub(p, a));
            let len = dot(u, u);
            let t = if len > 0.0 {
                (dot(u, d) / len).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let offset = std::array::from_fn(|i| d[i] - t * u[i]);
            dot(offset, offset)
        })
        .fold(f64::INFINITY, f64::min)
}
