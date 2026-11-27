//! P5 flanged housing: the same compute and reference checks on native and WebGPU.

use std::collections::{BTreeSet, HashMap};

use catallax_app::{build_certified, preview_body};
use catallax_mesh::float::{DisplayDescription, DisplayOptions, describe_body};
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
fn p5_software_compute() {
    pollster::block_on(async {
        let Some(adapter) = adapter(wgpu::Backends::VULKAN, true).await else {
            assert!(
                std::env::var_os("CI").is_none(),
                "no Vulkan software adapter in CI"
            );
            eprintln!("SKIP: no Vulkan software adapter (install Mesa lavapipe)");
            return;
        };
        let info = adapter.get_info();
        assert_eq!(info.device_type, wgpu::DeviceType::Cpu);
        check_p5(adapter).await;
    });
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
    assert_eq!(
        gpu.origin.map(f64::to_bits),
        description.origin.map(f64::to_bits)
    );
    assert_eq!(gpu.faces.len(), checked.faces.len());
    assert!(
        gpu.faces
            .iter()
            .any(|face| face.draw.as_ref().unwrap().sampled_grid)
    );
    assert!(
        gpu.faces
            .iter()
            .any(|face| !face.draw.as_ref().unwrap().sampled_grid)
    );
    let key = |i: u32| {
        let p = vertices[i as usize].position;
        [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
    };
    let mut edges = HashMap::<EdgeKey, (u32, i32)>::new();
    let mut next_index = 0;
    let mut used = vec![false; vertices.len()];
    for face in &gpu.faces {
        let draw = face.draw.as_ref().unwrap();
        assert_eq!(draw.indices.start, next_index);
        next_index = draw.indices.end;
        let draw_indices = &indices[draw.indices.start as usize..draw.indices.end as usize];
        assert!(!draw_indices.is_empty() && draw_indices.len().is_multiple_of(3));
        let checked_face = checked.faces.iter().find(|f| f.face == face.face).unwrap();
        let reference = DistanceGrid::new(
            &checked.triangles[checked_face.triangles.clone()],
            &checked.positions,
        );
        let face_vertices: BTreeSet<_> = draw_indices.iter().copied().collect();
        for index in face_vertices {
            used[index as usize] = true;
            let vertex = vertices[index as usize];
            let position = std::array::from_fn(|i| f64::from(vertex.position[i]) + gpu.origin[i]);
            let length_sq = vertex.normal[..3].iter().map(|x| x * x).sum::<f32>();
            assert!(
                (length_sq - 1.0).abs() < 1e-4,
                "face {:?}: bad normal {:?} at {position:?}",
                face.face,
                vertex.normal
            );
            assert!(
                reference.within_tolerance(position),
                "face {:?}: vertex {position:?} outside {TOLERANCE} mm of its checked face",
                face.face
            );
        }
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
    assert!(
        edges.values().all(|&usage| usage == (2, 0)),
        "GPU preview is not a closed, oriented mesh"
    );
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
