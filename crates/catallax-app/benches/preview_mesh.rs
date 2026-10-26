//! Development-only paired timings and sampled diagnostics, never proof evidence.
use std::{env, fs, time::Instant};

use catallax_mesh::{MeshOptions, PreviewMesh, PreviewOptions, mesh_body};
use catallax_numbers::{Decider, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use catallax_topology::{Model, Sense};
use serde_json::{Value, json};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).filter(|arg| arg != "--bench").collect();
    if args.is_empty() || args.len() > 3 {
        return Err("usage: preview_mesh OUTPUT.json [REPETITIONS [FIXTURE]]".into());
    }
    let repetitions: usize = args.get(1).map_or(Ok(3), |v| v.parse())?;
    if repetitions == 0 {
        return Err("repetitions must be positive".into());
    }
    let fixtures = [
        (
            "near-miss-flange",
            include_str!("../tests/fixtures/preview-near-miss-flange.json"),
            "0.01",
        ),
        (
            "fillet",
            include_str!("fine_mesh/fillet.json"),
            "0.0005462600113499065",
        ),
        (
            "chamfer-P2",
            include_str!("../../../parts/P2-collar.json"),
            "0.00036",
        ),
        ("chamfer", include_str!("fine_mesh/chamfer.json"), "0.00036"),
        (
            "countersink",
            include_str!("fine_mesh/countersink.json"),
            "0.000722841614740048",
        ),
        (
            "P1",
            include_str!("../../../parts/P1-bracket.json"),
            "0.0016638509548634457",
        ),
        (
            "P3",
            include_str!("../../../parts/P3-filleted-boss.json"),
            "0.0005462600113499065",
        ),
    ];
    let mut rows = Vec::new();
    for (name, input, chord) in fixtures {
        if args.get(2).is_some_and(|filter| filter != name) {
            continue;
        }
        let recipe = Recipe::parse(input)?;
        let (built, build_work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            catallax_features::evaluate(&recipe, &mut Decider::default())
        });
        let built = match built {
            Ok(built) if !build_work.exhausted => built,
            other => {
                rows.push(json!({"fixture":name,"build_refusal":format!("{other:?}"),"build_work":build_work.used}));
                continue;
            }
        };
        let options = MeshOptions {
            chord_tolerance: chord.parse()?,
            angular_deviation: Some("0.03490658503988659".parse()?),
            ..MeshOptions::default()
        };
        for iteration in 0..=repetitions {
            // One warmup, then alternating order on the same immutable body.
            let mut row = json!({"fixture":name,"iteration":iteration,"warmup":iteration == 0,
                "chord_target_mm":chord,"angle_target_rad":"0.03490658503988659",
                "allowance_each":catallax_app::DEFAULT_BUILD_WORK,"build_work":build_work.used});
            for preview in if iteration % 2 == 0 {
                [true, false]
            } else {
                [false, true]
            } {
                let start = Instant::now();
                let mut decider = Decider::default();
                decider.set_deadline(Some(start + std::time::Duration::from_secs(10)));
                if preview {
                    let (result, work) =
                        with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                            catallax_app::preview_body(
                                &built,
                                &PreviewOptions::from(&options),
                                &mut decider,
                            )
                        });
                    let ms = start.elapsed().as_secs_f64() * 1000.0;
                    row["preview_ms"] = json!(ms);
                    row["preview_work"] = json!(work.used);
                    row["preview_exhausted"] = json!(work.exhausted);
                    match result {
                        Ok(display) if !work.exhausted => {
                            row["preview_vertices"] = json!(display.mesh.positions.len());
                            row["preview_triangles"] = json!(display.mesh.triangles.len());
                            row["diagnostic_not_a_proof"] =
                                diagnostics(&built.model, &display.mesh)?;
                        }
                        other => row["preview_refusal"] = json!(format!("{other:?}")),
                    }
                } else {
                    let (result, work) =
                        with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                            mesh_body(&built.model, built.result, &options, &mut decider)
                        });
                    row["certified_ms"] = json!(start.elapsed().as_secs_f64() * 1000.0);
                    row["certified_work"] = json!(work.used);
                    row["certified_exhausted"] = json!(work.exhausted);
                    match result {
                        Ok(answer) if !work.exhausted => {
                            row["certified_vertices"] = json!(answer.value.positions.len());
                            row["certified_triangles"] = json!(answer.value.triangles.len());
                        }
                        other => row["certified_refusal"] = json!(format!("{other:?}")),
                    }
                }
            }
            println!("{row}");
            rows.push(row);
        }
    }
    fs::write(&args[0], serde_json::to_vec_pretty(&rows)?)?;
    Ok(())
}

fn diagnostics(model: &Model, mesh: &PreviewMesh) -> Result<Value, Box<dyn std::error::Error>> {
    let mut max_distance = 0.0_f64;
    let mut max_angle = 0.0_f64;
    let mut samples = 0;
    let mut degenerate = 0;
    let mut worst = Value::Null;
    let dec = Decider::default();
    for face in &mesh.faces {
        let definition = model.face(face.face);
        let surface = model.geom().surface(definition.surface);
        let scalar = |value: &Real| -> Result<f64, catallax_numbers::Refusal> {
            Ok(dec.enclose_at(value, 128)?.mid().to_f64_nearest())
        };
        let (origin, axis, first, second) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal, 0.0, 0.0),
            Surface::Cylinder(s) => (&s.origin, &s.axis, scalar(&s.radius())?, 0.0),
            Surface::Cone(s) => (&s.apex, &s.axis, scalar(&s.slope)?, 0.0),
            Surface::Sphere(s) => (&s.center, &s.axis, scalar(&s.radius)?, 0.0),
            Surface::Torus(s) => (&s.center, &s.axis, scalar(&s.major)?, scalar(&s.minor)?),
            _ => return Err("diagnostic does not cover this retained surface".into()),
        };
        let origin = [scalar(&origin.x)?, scalar(&origin.y)?, scalar(&origin.z)?];
        let mut axis = [scalar(&axis.x)?, scalar(&axis.y)?, scalar(&axis.z)?];
        let length = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
        axis = axis.map(|x| x / length);
        for triangle in &mesh.triangles[face.triangles.clone()] {
            let p = triangle.vertices.map(|i| mesh.positions[i]);
            let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
            let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
            let cross: [f64; 3] = std::array::from_fn(|i| {
                a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
            });
            let norm = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm == 0.0 {
                degenerate += 1;
                continue;
            }
            for w in [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [1.0 / 3.0; 3],
                [0.2, 0.3, 0.5],
            ] {
                let d: [f64; 3] = std::array::from_fn(|i| {
                    (0..3).map(|j| w[j] * p[j][i]).sum::<f64>() - origin[i]
                });
                let z = (0..3).map(|i| d[i] * axis[i]).sum::<f64>();
                let radial: [f64; 3] = std::array::from_fn(|i| d[i] - z * axis[i]);
                let r = radial.iter().map(|x| x * x).sum::<f64>().sqrt();
                let (distance, normal) = match surface {
                    Surface::Plane(_) => (z.abs(), axis),
                    Surface::Cylinder(_) => ((r - first).abs(), radial.map(|x| x / r)),
                    Surface::Cone(_) => (
                        (r - first * z).abs() / (1.0 + first * first).sqrt(),
                        std::array::from_fn(|i| radial[i] / r - first * axis[i]),
                    ),
                    Surface::Sphere(_) => {
                        let radius = d.iter().map(|x| x * x).sum::<f64>().sqrt();
                        ((radius - first).abs(), d.map(|x| x / radius))
                    }
                    Surface::Torus(_) => {
                        let dr = r - first;
                        (
                            (dr.hypot(z) - second).abs(),
                            std::array::from_fn(|i| radial[i] * dr / r + axis[i] * z),
                        )
                    }
                    _ => unreachable!(),
                };
                let n = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
                let sign = if definition.sense == Sense::Same {
                    1.0
                } else {
                    -1.0
                };
                let cosine = sign * (0..3).map(|i| cross[i] * normal[i]).sum::<f64>() / (norm * n);
                if cosine.is_finite() {
                    max_angle = max_angle.max(cosine.clamp(-1.0, 1.0).acos());
                }
                if distance > max_distance {
                    max_distance = distance;
                    worst = json!({"face":face.face.index(),
                        "triangle":p,"weights":w});
                }
                samples += 1;
            }
        }
    }
    Ok(
        json!({"sampled_max_distance_mm":max_distance,"sampled_max_normal_angle_rad":max_angle,
        "samples":samples,"degenerate_triangles":degenerate,
        "worst_sample":worst,
        "scope":"eight points per triangle against retained analytic carrier; not trim, Hausdorff, watertightness or deviation proof"}),
    )
}
