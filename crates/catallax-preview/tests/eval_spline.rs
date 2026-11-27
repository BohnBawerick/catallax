//! Read back the same WGSL fragment the mesher will use. Set
//! `CATALLAX_SPLINE_ADAPTER=cpu` or `hardware` to require that adapter class.
//! Set `P5_SOURCES` and `P5_SELECTION` to include the held-out STEP rows.

#![cfg(not(target_arch = "wasm32"))]

mod corpus;

#[path = "../src/eval/spline.rs"]
mod spline;

use bytemuck::{Pod, Zeroable, cast_slice};
use catallax_features::EvaluationOutcome;
use catallax_numbers::Decider;
use catallax_shapes::{BSplineSurface, SplineSupport, Surface};
use wgpu::util::DeviceExt;

use spline::{MAX_DEGREE, Net};

const HARNESS: &str = r#"
struct Query {
    net: SplineNet,
    u: f32,
    v: f32,
}

struct Answer {
    position: vec4f,
    normal: vec4f,
}

@group(0) @binding(0) var<storage, read> spline_controls: array<vec4f>;
@group(0) @binding(1) var<storage, read> spline_knots: array<f32>;
@group(0) @binding(2) var<storage, read> queries: array<Query>;
@group(0) @binding(3) var<storage, read_write> answers: array<Answer>;

@compute @workgroup_size(64)
fn evaluate(@builtin(global_invocation_id) id: vec3u) {
    if id.x >= arrayLength(&queries) { return; }
    let q = queries[id.x];
    let sample = spline_sample(q.net, vec2f(q.u, q.v));
    answers[id.x] = Answer(vec4f(sample.position, 0.0), vec4f(sample.normal, 0.0));
}
"#;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Query {
    net: Net,
    uv: [f32; 2],
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Answer {
    position: [f32; 4],
    normal: [f32; 4],
}

#[test]
fn spline_shader_validates_without_a_device() {
    let source = [spline::WGSL, HARNESS].concat();
    let module = naga::front::wgsl::parse_str(&source)
        .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
    naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::empty(),
    )
    .validate(&module)
    .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
}

struct Carrier {
    degree: [usize; 2],
    count: [usize; 2],
    knots: [Vec<f64>; 2],
    controls: Vec<[f64; 4]>,
    size: f64,
}

impl Carrier {
    fn stored(surface: &BSplineSurface) -> Self {
        let points: Vec<_> = surface
            .controls()
            .iter()
            .map(|p| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()])
            .collect();
        let first = points[0];
        let size = norm(std::array::from_fn(|axis| {
            let (lo, hi) = points
                .iter()
                .fold((first[axis], first[axis]), |(lo, hi), p| {
                    (lo.min(p[axis]), hi.max(p[axis]))
                });
            hi - lo
        }));
        let controls = points
            .iter()
            .zip(surface.weights())
            .map(|(p, w)| {
                let w = w.to_f64();
                [
                    (p[0] - first[0]) * w,
                    (p[1] - first[1]) * w,
                    (p[2] - first[2]) * w,
                    w,
                ]
            })
            .collect();
        Self {
            degree: surface.degrees(),
            count: surface.control_count(),
            knots: surface
                .knots()
                .map(|knots| knots.iter().map(|k| k.to_f64_nearest()).collect()),
            controls,
            size,
        }
    }

    fn reference(&self, uv: [f64; 2]) -> Reference {
        let [p, q] = self.degree;
        let [ku, kv] = self.knots.each_ref().map(Vec::as_slice);
        let i = span(ku, p, self.count[0], uv[0]);
        let j = span(kv, q, self.count[1], uv[1]);
        let mut rows = [[0.0; 4]; MAX_DEGREE + 1];
        let mut v_slopes = rows;
        for a in 0..=p {
            let offset = (i - p + a) * self.count[1] + j - q;
            (rows[a], v_slopes[a]) = de_boor(kv, q, j, uv[1], &self.controls[offset..=offset + q]);
        }
        let (h, hu) = de_boor(ku, p, i, uv[0], &rows[..=p]);
        let (hv, _) = de_boor(ku, p, i, uv[0], &v_slopes[..=p]);
        let position = std::array::from_fn(|a| h[a] / h[3]);
        let tangent = |d: [f64; 4]| -> [f64; 3] {
            std::array::from_fn(|a| (d[a] - position[a] * d[3]) / h[3])
        };
        let u = tangent(hu);
        let v = tangent(hv);
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let length = norm(n);
        let normal =
            (length > 64.0 * f64::EPSILON * self.size.powi(2)).then(|| n.map(|x| x / length));
        Reference { position, normal }
    }
}

fn span(knots: &[f64], degree: usize, count: usize, t: f64) -> usize {
    (degree..count)
        .rev()
        .find(|&i| knots[i] <= t && knots[i] < knots[i + 1])
        .unwrap()
}

// Independent binary64 reference: differentiate the de Boor interpolation
// steps, rather than evaluating basis functions and their differences.
fn de_boor(
    knots: &[f64],
    p: usize,
    span: usize,
    t: f64,
    controls: &[[f64; 4]],
) -> ([f64; 4], [f64; 4]) {
    let mut value = [[0.0; 4]; MAX_DEGREE + 1];
    let mut slope = value;
    value[..=p].copy_from_slice(controls);
    for level in 1..=p {
        for j in (level..=p).rev() {
            let i = span - p + j;
            let width = knots[i + p + 1 - level] - knots[i];
            let a = (t - knots[i]) / width;
            for axis in 0..4 {
                let delta = value[j][axis] - value[j - 1][axis];
                slope[j][axis] =
                    (1.0 - a) * slope[j - 1][axis] + a * slope[j][axis] + delta / width;
                value[j][axis] = (1.0 - a) * value[j - 1][axis] + a * value[j][axis];
            }
        }
    }
    (value[p], slope[p])
}

fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn error(actual: [f32; 4], expected: [f64; 3]) -> f64 {
    norm(std::array::from_fn(|i| f64::from(actual[i]) - expected[i]))
}

struct Reference {
    position: [f64; 3],
    normal: Option<[f64; 3]>,
}

#[derive(Default)]
struct Batch {
    controls: Vec<[f32; 4]>,
    knots: Vec<f32>,
    queries: Vec<Query>,
    expected: Vec<Reference>,
    faces: Vec<(String, f64, std::ops::Range<usize>)>,
}

impl Batch {
    fn add(&mut self, name: String, carrier: Carrier) {
        let [p, q] = carrier.degree;
        let [u, v] = carrier.count;
        assert!(
            p <= MAX_DEGREE && q <= MAX_DEGREE && carrier.size > 0.0,
            "{name}: {p}x{q}, {u}x{v}"
        );
        let net = Net {
            control_offset: self.controls.len().try_into().unwrap(),
            knot_u_offset: self.knots.len().try_into().unwrap(),
            knot_v_offset: (self.knots.len() + carrier.knots[0].len())
                .try_into()
                .unwrap(),
            degree_u: p as u32,
            degree_v: q as u32,
            count_u: u as u32,
            count_v: v as u32,
        };
        self.controls
            .extend(carrier.controls.iter().map(|c| c.map(|x| x as f32)));
        for knots in &carrier.knots {
            for w in knots.windows(2) {
                assert!(
                    w[0].total_cmp(&w[1]).is_eq() || (w[0] as f32) < (w[1] as f32),
                    "{name}: binary32 collapsed a knot span: {w:?}"
                );
            }
            self.knots.extend(knots.iter().map(|&k| k as f32));
        }
        let samples = |axis: usize| {
            let k = &carrier.knots[axis][carrier.degree[axis]..=carrier.count[axis]];
            let mut samples = vec![k[0]];
            for w in k.windows(2).filter(|w| w[0] < w[1]) {
                for f in [0.125, 0.5, 0.875] {
                    samples.push(w[0] + f * (w[1] - w[0]));
                }
                // Test the right-hand convention at exactly representable knots.
                if f64::from(w[1] as f32).total_cmp(&w[1]).is_eq() {
                    samples.push(w[1]);
                }
            }
            if samples.last() != k.last() {
                samples.push(*k.last().unwrap());
            }
            samples
        };
        let vs = samples(1);
        let start = self.queries.len();
        for u in samples(0) {
            for &v in &vs {
                self.queries.push(Query {
                    net,
                    uv: [u as f32, v as f32],
                });
                self.expected.push(carrier.reference([u, v]));
            }
        }
        self.faces
            .push((name, carrier.size, start..self.queries.len()));
    }

    async fn compare(&self, adapter: &wgpu::Adapter) {
        let (device, queue) = adapter.request_device(&Default::default()).await.unwrap();
        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Spline evaluation test"),
            source: wgpu::ShaderSource::Wgsl([spline::WGSL, HARNESS].concat().into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Spline evaluation test"),
            layout: None,
            module: &shader,
            entry_point: Some("evaluate"),
            compilation_options: Default::default(),
            cache: None,
        });
        assert!(
            scope.pop().await.is_none(),
            "spline pipeline validation failed"
        );
        let input = |contents| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents,
                usage: wgpu::BufferUsages::STORAGE,
            })
        };
        let output = |usage| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: (self.queries.len() * std::mem::size_of::<Answer>()) as u64,
                usage,
                mapped_at_creation: false,
            })
        };
        let readback = output(wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ);
        let answers = output(wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC);
        let buffers = [
            input(cast_slice(&self.controls)),
            input(cast_slice(&self.knots)),
            input(cast_slice(&self.queries)),
            answers,
        ];
        let entries = std::array::from_fn::<_, 4, _>(|i| wgpu::BindGroupEntry {
            binding: i as u32,
            resource: buffers[i].as_entire_binding(),
        });
        let group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        });
        let mut encoder = device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &group, &[]);
            pass.dispatch_workgroups((self.queries.len() as u32).div_ceil(64), 1, 1);
        }
        encoder.copy_buffer_to_buffer(&buffers[3], 0, &readback, 0, readback.size());
        let submission = queue.submit([encoder.finish()]);
        let slice = readback.slice(..);
        let (send, recv) = futures_channel::oneshot::channel();
        slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = send.send(result);
        });
        device
            .poll(wgpu::PollType::Wait {
                submission_index: Some(submission),
                timeout: Some(std::time::Duration::from_secs(30)),
            })
            .unwrap();
        recv.await.unwrap().unwrap();
        let mapped = slice.get_mapped_range().unwrap();
        let result: &[Answer] = cast_slice(&mapped);
        assert_eq!(result.len(), self.expected.len());
        let mut worst = [0.0_f64; 2];
        for (name, size, range) in &self.faces {
            let mut face = [0.0_f64; 2];
            for i in range.clone() {
                let uv = self.queries[i].uv;
                let got = &result[i];
                let want = &self.expected[i];
                let position = error(got.position, want.position) / size;
                assert!(
                    position <= 5e-6,
                    "{name} sample {i} {uv:?}: relative position error {position:e}, got {:?}, expected {:?}",
                    got.position,
                    want.position
                );
                face[0] = face[0].max(position);
                if let Some(n) = want.normal {
                    let normal = error(got.normal, n);
                    assert!(
                        normal <= 3e-3,
                        "{name} {uv:?}: normal chord error {normal:e}, got {:?}, expected {n:?}",
                        got.normal
                    );
                    face[1] = face[1].max(normal);
                } else {
                    assert!(
                        got.normal.into_iter().all(f32::is_finite),
                        "{name} {uv:?}: singular normal is not finite"
                    );
                }
            }
            eprintln!(
                "{name}: position/size {:.3e}, unit normal chord {:.3e}",
                face[0], face[1]
            );
            for i in 0..2 {
                worst[i] = worst[i].max(face[i]);
            }
        }
        let info = adapter.get_info();
        eprintln!(
            "{} {:?} {:?}: {} nets, {} samples, worst position/size {:.3e}, normal chord {:.3e}",
            info.name,
            info.backend,
            info.device_type,
            self.faces.len(),
            result.len(),
            worst[0],
            worst[1]
        );
    }
}

fn add_face(batch: &mut Batch, row: &str, index: usize, surface: &Surface) {
    let Surface::BSpline(support) = surface else {
        return;
    };
    if let SplineSupport::RationalAtlas(atlas) = support {
        for (chart, entry) in atlas.charts().iter().enumerate() {
            batch.add(
                format!("{row}/{index}/{chart}"),
                Carrier::stored(entry.patch().carrier()),
            );
        }
    } else {
        batch.add(format!("{row}/{index}"), Carrier::stored(support.carrier()));
    }
}

fn add_corpus(batch: &mut Batch) {
    let mut exact = 0;
    let mut fitted = 0;
    for (id, built) in corpus::bodies() {
        match built {
            EvaluationOutcome::Exact(checkpoint) => {
                exact += 1;
                let built = checkpoint.built();
                for (index, face) in built
                    .model
                    .faces_of_body(built.result)
                    .into_iter()
                    .enumerate()
                {
                    add_face(
                        batch,
                        id,
                        index,
                        built.model.geom().surface(built.model.face(face).surface),
                    );
                }
            }
            EvaluationOutcome::Numerical(built) => {
                fitted += 1;
                for (index, face) in built.body().definition().faces.iter().enumerate() {
                    add_face(batch, id, index, &face.surface);
                }
            }
        }
    }
    eprintln!("Display corpus: {exact} exact, {fitted} fitted");
}

#[test]
fn stored_splines_match_binary64() {
    pollster::block_on(check_stored_splines());
}

async fn check_stored_splines() {
    let required = std::env::var("CATALLAX_SPLINE_ADAPTER").ok();
    let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle();
    let software = required.as_deref() == Some("cpu");
    if software {
        descriptor.backends = wgpu::Backends::VULKAN;
    }
    let instance = wgpu::Instance::new(descriptor);
    let adapter = match instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: software,
            ..Default::default()
        })
        .await
    {
        Ok(adapter) => adapter,
        Err(reason) if required.is_none() => {
            eprintln!("SKIP spline readback: {reason}");
            return;
        }
        Err(reason) => panic!("{reason}"),
    };
    let info = adapter.get_info();
    match required.as_deref() {
        Some("cpu") => assert_eq!(info.device_type, wgpu::DeviceType::Cpu, "{info:?}"),
        Some("hardware") => assert!(
            matches!(
                info.device_type,
                wgpu::DeviceType::IntegratedGpu | wgpu::DeviceType::DiscreteGpu
            ),
            "{info:?}"
        ),
        None => (),
        Some(value) => panic!("unknown CATALLAX_SPLINE_ADAPTER={value}"),
    }
    let mut batch = Batch::default();
    add_corpus(&mut batch);
    assert!(
        batch
            .faces
            .iter()
            .any(|(name, _, _)| name.starts_with("fillet-oblique-cone/"))
    );
    assert!(
        batch
            .faces
            .iter()
            .any(|(name, _, _)| name.starts_with("fillet-cross-unequal/"))
    );
    // Exercise the size limit even if no individual fixture reaches both limits.
    let large = Carrier {
        degree: [8, 2],
        count: [257, 3],
        knots: [
            std::iter::repeat_n(0.0, 9)
                .chain((1..249).map(|i| f64::from(i) / 256.0))
                .chain(std::iter::repeat_n(1.0, 9))
                .collect(),
            vec![0., 0., 0., 1., 1., 1.],
        ],
        controls: (0..257)
            .flat_map(|i| {
                (0..3).map(move |j| {
                    let u = f64::from(i) / 256.0;
                    let v = f64::from(j) / 2.0;
                    let w = 1.0 + u * v;
                    [u * w, v * w, u * v * w, w]
                })
            })
            .collect(),
        size: 3.0_f64.sqrt(),
    };
    batch.add("degree-8x2-net-257x3".to_owned(), large);
    batch.compare(&adapter).await;
}
