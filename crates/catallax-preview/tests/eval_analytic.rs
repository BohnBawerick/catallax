//! Compares the actual compute evaluator with an independent binary64 chart evaluation.
// Set CATALLAX_ANALYTIC_CORPUS=1 to add the predictor corpus. Imported
// rows need CATALLAX_ANALYTIC_SOURCES, a JSON object from model IDs to STEP paths.
#[path = "../src/eval/analytic.rs"]
mod analytic;

use analytic::{AnalyticSurface, AnalyticVertex, Definition, Kind, WGSL};
use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational, Real};
use catallax_preview::gpu::Gpu;
use catallax_recipe::Recipe;
use catallax_shapes::{
    BezierExtrusion, Cone, Cylinder, Plane, Point3, RPoint3, RVector3, RollingPipe, Sphere,
    Surface, Torus, Vector3,
};

const DRIVER: &str = r"
@group(0) @binding(0) var<storage, read> faces: array<AnalyticSurface>;
@group(0) @binding(1) var<storage, read> samples: array<vec3<f32>>;
@group(0) @binding(2) var<storage, read_write> vertices: array<AnalyticVertex>;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if id.x < arrayLength(&samples) {
        let sample = samples[id.x];
        vertices[id.x] = eval_analytic(faces[u32(sample.z)], sample.xy);
    }
}
";

fn shader() -> String {
    [WGSL, DRIVER].concat()
}

#[test]
fn analytic_shader_validates_without_an_adapter() {
    let source = shader();
    let module = naga::front::wgsl::parse_str(&source)
        .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
    naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::empty(),
    )
    .validate(&module)
    .unwrap();
}

type V = [f64; 3];

fn xyz([x, y, z, _]: [f64; 4]) -> V {
    [x, y, z]
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mul(k: f64, a: V) -> V {
    a.map(|x| k * x)
}
fn norm(a: V) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn unit(a: V) -> V {
    mul(norm(a).recip(), a)
}
fn cross(a: V, b: V) -> V {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn circle(u: f64) -> (f64, f64) {
    (std::f64::consts::TAU * u).sin_cos()
}
fn quarters(u: f64) -> (f64, f64) {
    (std::f64::consts::FRAC_PI_2 * u.floor() + 2.0 * u.fract().rem_euclid(1.0).atan()).sin_cos()
}

fn point(face: &Definition, [u, v]: [f64; 2]) -> V {
    let [p, axis, x, y, extra, parameters, section] = face.data;
    let [p, axis, x, y, e] = [p, axis, x, y, extra].map(xyz);
    let [a, b, c, d] = parameters;
    let (sin, cos) = match face.kind {
        Kind::Pipe | Kind::Ruled => quarters(u),
        _ => circle(u),
    };
    let radial = add(mul(cos, x), mul(sin, y));
    match face.kind {
        Kind::Cone => add(p, add(mul(v, axis), mul(a * v, radial))),
        Kind::Sphere => {
            let (sin, cos) = circle(v);
            add(p, mul(a, add(mul(cos, radial), mul(sin, axis))))
        }
        Kind::Torus => {
            let (sin, cos) = circle(v);
            add(p, add(mul(a + b * cos, radial), mul(b * sin, axis)))
        }
        Kind::Sweep => {
            let t = 1.0 - u;
            let profile = if a == 2.0 {
                add(mul(t * t, p), add(mul(2.0 * t * u, axis), mul(u * u, x)))
            } else {
                add(
                    add(mul(t.powi(3), p), mul(3.0 * t * t * u, axis)),
                    add(mul(3.0 * t * u * u, x), mul(u.powi(3), y)),
                )
            };
            add(profile, mul(v, e))
        }
        Kind::Pipe => {
            let centre = add(
                p,
                mul(
                    (1.0 + d * cos).recip(),
                    add(mul(a, radial), mul(b * cos, axis)),
                ),
            );
            let last = mul(section[0], sub(radial, mul(c, axis)));
            let chord = unit(add(mul(2.0 - v, e), mul(v, last)));
            // Reflection across the normalized half-angle chord.
            let projection: f64 = (0..3).map(|i| chord[i] * e[i]).sum();
            add(
                centre,
                mul(section[1], sub(mul(2.0 * projection, chord), e)),
            )
        }
        Kind::Ruled => {
            let z = a + (b - a) * v;
            add(
                p,
                add(
                    mul(z, axis),
                    mul(c - d * z + extra[0] * cos + extra[1] * sin, radial),
                ),
            )
        }
        Kind::LinearFillet => {
            let t = 1.0 - u;
            let scale = (a + (b - a) * v) / (t * t + 2.0 * c * u * t + u * u);
            add(
                add(p, mul(v, axis)),
                add(mul(scale * t * t, x), mul(scale * u * u, y)),
            )
        }
    }
}

fn reference(face: &Definition, same: bool, [u, v]: [f64; 2]) -> (V, V) {
    let p = point(face, [u, v]);
    let [origin, axis, x, y, _, radii, _] = face.data;
    let (sin, cos) = circle(if face.kind == Kind::Cone && v == 0.0 {
        0.0
    } else {
        u
    });
    let radial = add(mul(cos, xyz(x)), mul(sin, xyz(y)));
    let n = match face.kind {
        Kind::Cone => sub(radial, mul(radii[0], xyz(axis))),
        Kind::Sphere => sub(p, xyz(origin)),
        Kind::Torus => sub(p, add(xyz(origin), mul(radii[0], radial))),
        _ => {
            // At a collapsed fillet row the normal is the interior limit.
            // Elsewhere use position differences instead of the shader's formulas.
            let h = 1.0e-6;
            let normal_v = if face.kind == Kind::LinearFillet {
                0.5
            } else {
                v
            };
            cross(
                sub(
                    point(face, [u + h, normal_v]),
                    point(face, [u - h, normal_v]),
                ),
                sub(point(face, [u, v + h]), point(face, [u, v - h])),
            )
        }
    };
    (p, mul(if same { 1.0 } else { -1.0 }, unit(n)))
}

struct Case {
    name: String,
    definition: Definition,
    same: bool,
}

fn case_recipe(text: &str, id: &str) -> serde_json::Value {
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap_or_else(|| panic!("No recipe for case {id}"))["recipe"]
        .take()
}

fn corpus(mut visit: impl FnMut(&str, &catallax_features::Built)) {
    use serde_json::{Value, json};
    use std::{fs, path::Path};

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let manifest: Value =
        serde_json::from_str(include_str!("../../../tools/predictor/corpus.json")).unwrap();
    let sources: Value = std::env::var_os("CATALLAX_ANALYTIC_SOURCES")
        .map_or(Value::Null, |path| {
            serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
        });
    let mut exact = 0;
    let mut refused = 0;
    let mut missing = 0;
    for row in manifest["rows"].as_array().unwrap() {
        let id = row["id"].as_str().unwrap();
        let mut files = catallax_features::ImportFiles::new();
        let bytes;
        let mut recipe = if let Some(path) = row["recipe"].as_str() {
            serde_json::from_slice(&fs::read(root.join(path)).unwrap()).unwrap()
        } else if let Some(path) = row["cases"].as_str() {
            case_recipe(
                &fs::read_to_string(root.join(path)).unwrap(),
                row["case"].as_str().unwrap(),
            )
        } else if let Some(model) = row["import"]["model"].as_str() {
            let Some(path) = sources[model].as_str() else {
                assert!(sources.is_null(), "No supplied STEP path for {model}");
                missing += 1;
                continue;
            };
            let hash = row["import"]["sha256"].as_str().unwrap();
            bytes = fs::read(path).unwrap();
            files.insert(hash.into(), bytes.as_slice());
            json!({"format": "catallax-recipe/3", "name": id, "steps": [{
                "id": "imported", "op": "import_step",
                "file": {"name": "part.step", "sha256": hash},
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}
            }], "result": "imported"})
        } else {
            assert!(row["inline"].is_object(), "{id}: missing recipe");
            row["inline"].clone()
        };
        if let Some(steps) = row["append"].as_array() {
            recipe["steps"]
                .as_array_mut()
                .unwrap()
                .extend(steps.iter().cloned());
            recipe["result"] = row["result"].clone();
        }
        if row["format"].is_string() {
            recipe["format"] = row["format"].clone();
        }
        let recipe = Recipe::parse(&recipe.to_string()).unwrap_or_else(|e| panic!("{id}: {e:?}"));
        match catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default()) {
            Ok(built) => {
                exact += 1;
                visit(id, &built);
            }
            Err(reason) => {
                refused += 1;
                eprintln!("{id}: no exact faces ({reason})");
            }
        }
    }
    eprintln!(
        "Corpus: {exact} exact rows, {refused} without an exact body, {missing} missing STEP sources"
    );
}

fn fixtures() -> Vec<Case> {
    let mut cases = Vec::new();
    let mut add = |name: &str, surface: &Surface| {
        for same in [true, false] {
            if let Some(definition) = Definition::from_surface(surface, [13.0, -4.0, 9.0], same) {
                cases.push(Case {
                    name: name.into(),
                    definition,
                    same,
                });
            }
        }
    };
    let mut add_body = |name: &str, built: &catallax_features::Built| {
        let model = &built.model;
        for &lump in &model.body(built.result).lumps {
            for &shell in &model.lump(lump).shells {
                for &face in &model.shell(shell).faces {
                    let surface = model.geom().surface(model.face(face).surface);
                    add(&format!("{name} {face:?}"), surface);
                }
            }
        }
    };
    for (name, input) in [
        (
            "P0 sphere",
            include_str!("../../../parts/P0/equal-cylinder-sphere.json"),
        ),
        ("P2", include_str!("../../../parts/P2-collar.json")),
        ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
        (
            "oblique draft",
            include_str!("../../catallax-features/tests/fixtures/draft-oblique.json"),
        ),
        (
            "zero-end fillet",
            &case_recipe(
                include_str!("../../../parts/fillets/cases.jsonl"),
                "variable-zero-end",
            )
            .to_string(),
        ),
    ] {
        let built = evaluate(&Recipe::parse(input).unwrap(), &mut Decider::default())
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        add_body(name, &built);
    }
    if std::env::var_os("CATALLAX_ANALYTIC_CORPUS").is_some() {
        corpus(add_body);
    }

    let r = Real::from_i64;
    let q = Rational::from_i64;
    let vector = |xyz: [i64; 3]| {
        let [x, y, z] = xyz.map(r);
        Vector3::new(x, y, z)
    };
    let mut dec = Decider::default();
    add(
        "rotated sphere poles",
        &Surface::Sphere(Sphere {
            center: Point3::new(r(100), r(-31), r(62)),
            radius: r(7),
            axis: vector([3, 0, 4]),
            ref_dir: vector([0, 1, 0]),
        }),
    );
    add(
        "inclined ring torus",
        &Surface::Torus(Torus {
            center: Point3::origin(),
            major: r(8),
            minor: r(3),
            axis: vector([1, 2, 2]),
            ref_dir: vector([0, 1, -1]),
        }),
    );
    let plane = Plane {
        origin: Point3::new(r(0), r(0), r(3)),
        normal: vector([3, 0, 4]),
        ref_dir: vector([0, 1, 0]),
    };
    let cylinder = Cylinder {
        origin: Point3::origin(),
        axis: vector([0, 0, 1]),
        radius_sq: r(36),
        ref_dir: vector([1, 0, 0]),
    };
    let cone = Cone {
        apex: Point3::new(r(0), r(0), r(8)),
        axis: vector([0, 0, -1]),
        slope: Real::from(Rational::from_ratio(3, 4).unwrap()),
        ref_dir: vector([1, 0, 0]),
    };
    let radius = Real::from(Rational::from_ratio(1, 2).unwrap());
    for (name, pipe) in [
        (
            "cylindrical pipe",
            RollingPipe::new(&plane, &cylinder, [true; 2], radius.clone(), true, &mut dec),
        ),
        (
            "conical pipe",
            RollingPipe::conical(&plane, &cone, [true; 2], radius, true, &mut dec),
        ),
    ] {
        add(name, &Surface::BSpline(pipe.unwrap().into()));
    }
    add("cone apex", &Surface::Cone(cone));
    for count in [3, 4] {
        let controls = [[0, 0], [1, 3], [4, 3], [6, 1]][..count]
            .iter()
            .map(|&[x, y]| RPoint3::new(q(x), q(y), q(0)))
            .collect();
        let sweep =
            BezierExtrusion::new(controls, RVector3::new(q(0), q(0), q(5)), &mut dec).unwrap();
        add("Bezier sweep", &Surface::BSpline(sweep.into()));
    }
    cases
}

fn grid(kind: Kind) -> (Vec<f32>, Vec<f32>) {
    // Exact seam aliases, one sample either side of each wrap, and asymmetric
    // interior points. Sphere latitude endpoints are the exact chart poles.
    let turn = || {
        vec![
            -1.0, -0.015625, 0.0, 0.015625, 0.1875, 0.25, 0.5, 0.75, 0.984375, 1.0,
        ]
    };
    let segment = || vec![0.0, 0.125, 0.625, 0.9375, 1.0];
    match kind {
        Kind::Cone => (turn(), vec![0.0, 0.03125, 0.5, 2.0]),
        Kind::Sphere => (
            turn(),
            vec![-0.25, -0.234375, -0.125, 0.0, 0.1875, 0.234375, 0.25],
        ),
        Kind::Torus => (turn(), turn()),
        Kind::Sweep | Kind::LinearFillet => (segment(), segment()),
        Kind::Pipe | Kind::Ruled => (
            vec![
                -0.015625, 0.0, 0.015625, 0.375, 1.0, 1.625, 2.0, 3.0, 3.984375, 4.0,
            ],
            segment(),
        ),
    }
}

#[test]
fn analytic_compute_matches_binary64_charts() {
    let requirement = std::env::var("CATALLAX_ANALYTIC_ADAPTER").ok();
    let gpu = match Gpu::new(false) {
        Ok(gpu) => gpu,
        Err(reason) => {
            assert!(
                requirement.is_none() && std::env::var_os("CATALLAX_REQUIRE_GPU").is_none(),
                "{reason}"
            );
            eprintln!("Analytic compute skipped: {reason}");
            return;
        }
    };
    eprintln!("Analytic compute adapter: {}", gpu.name());
    let device_type = gpu.name().split_whitespace().last().unwrap();
    match requirement.as_deref() {
        Some("hardware") => assert!(
            matches!(device_type, "IntegratedGpu" | "DiscreteGpu"),
            "{}",
            gpu.name()
        ),
        Some("software") => assert_eq!(device_type, "Cpu", "{}", gpu.name()),
        None => (),
        Some(other) => panic!("Unknown CATALLAX_ANALYTIC_ADAPTER: {other}"),
    }

    let cases = fixtures();
    let faces: Vec<AnalyticSurface> = cases.iter().map(|case| case.definition.upload()).collect();
    let mut samples = Vec::new();
    let mut expected = Vec::new();
    for (i, case) in cases.iter().enumerate() {
        let (us, vs) = grid(case.definition.kind);
        let domain = us.iter().flat_map(|&u| vs.iter().map(move |&v| [u, v]));
        let start = samples.len();
        for [u, v] in domain {
            let (p, n) = reference(&case.definition, case.same, [f64::from(u), f64::from(v)]);
            samples.push([u, v, i as f32, 0.0]);
            expected.push((p, n));
        }
        let first = expected[start].0;
        let size = expected[start..]
            .iter()
            .map(|(p, _)| norm(sub(*p, first)))
            .fold(0.0, f64::max);
        assert!(size > 0.0, "{}", case.name);
    }
    let zero = vec![[0.0_f32; 8]; samples.len()];
    let actual: Vec<AnalyticVertex> = gpu
        .dispatch(
            &shader(),
            "main",
            &[
                bytemuck::cast_slice(&faces),
                bytemuck::cast_slice(&samples),
                bytemuck::cast_slice(&zero),
            ],
            2,
            [(samples.len() as u32).div_ceil(64), 1, 1],
        )
        .unwrap();

    let mut start = 0;
    assert_eq!(actual.len(), samples.len());
    for case in &cases {
        let (us, vs) = grid(case.definition.kind);
        let end = start + us.len() * vs.len();
        let first = expected[start].0;
        let size = expected[start..end]
            .iter()
            .map(|(p, _)| norm(sub(*p, first)))
            .fold(0.0, f64::max);
        // Position error is relative to the sampled face's extent, plus one
        // body-local upload rounding allowance. Normal error is unit-vector chord distance.
        let position_bound = 4.0e-6 * size + 4.0 * f64::from(f32::EPSILON) * norm(first);
        let normal_bound = 5.0e-5;
        for j in start..end {
            let [u, v, ..] = samples[j];
            let (p, n) = expected[j];
            let position_error = norm(sub(xyz(actual[j].position.map(f64::from)), p));
            let normal_error = norm(sub(xyz(actual[j].normal.map(f64::from)), n));
            assert!(
                position_error <= position_bound && normal_error <= normal_bound,
                "{} {:?} same={} [{u}, {v}]: position {position_error:e} > {position_bound:e} or normal {normal_error:e} > {normal_bound:e}",
                case.name,
                case.definition.kind,
                case.same
            );
        }
        start = end;
    }
    for kind in [
        Kind::Cone,
        Kind::Sphere,
        Kind::Torus,
        Kind::Sweep,
        Kind::Pipe,
        Kind::Ruled,
        Kind::LinearFillet,
    ] {
        assert!(
            cases.iter().any(|case| case.definition.kind == kind),
            "missing {kind:?}"
        );
    }
    eprintln!(
        "Compared {} samples across {} oriented faces",
        samples.len(),
        cases.len()
    );
}
