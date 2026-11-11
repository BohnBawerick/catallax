//! Unproved closed-patch ownership for P3, and its independent exact F8 checks.
use crate::{OWNERSHIP_ALLOWANCE, Result, WIDTH};
use catallax_field::{
    Plan,
    distance::{NearestOptions, Ownership, Prepared},
};
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
use catallax_recipe::Role;
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

pub const SAMPLE_COUNT: usize = 1024;
pub const SEED: u64 = 0xf9_2026_1111;

pub fn named_faces(plan: &Plan) -> Result<[u32; 13]> {
    let mut ids = [0; 13];
    for (i, face) in plan.faces().iter().enumerate() {
        if plan.result_face(i).is_none() {
            continue;
        }
        let side = |s: &str| match s {
            "x_min" => Some(0),
            "x_max" => Some(1),
            "y_min" => Some(2),
            "y_max" => Some(3),
            _ => None,
        };
        let slot = match (face.origin.step.as_str(), &face.origin.role) {
            ("block", Role::Start) => Some(0),
            ("block", Role::End) => Some(1),
            ("block", Role::Side { segment, .. }) => side(segment.as_str()).map(|s| s + 2),
            ("boss", Role::End) => Some(6),
            ("boss", Role::Side { .. }) => Some(7),
            ("junction", Role::Blend { .. }) => Some(8),
            ("rounded", Role::Blend { edge }) => edge.between().iter().find_map(|o| {
                if let Role::Side { segment, .. } = &o.role {
                    side(segment.as_str()).map(|s| s + 9)
                } else {
                    None
                }
            }),
            _ => None,
        };
        if let Some(slot) = slot {
            ids[slot] = plan.select(&face.origin).unwrap().index() as u32 + 1;
        }
    }
    if ids.contains(&0) {
        return Err(Refusal::new(
            "preview.ownership-shape",
            "P3 requires all thirteen named closed patches",
        )
        .into());
    }
    Ok(ids)
}
pub fn source(ids: &[u32; 13], radius: f64) -> String {
    format!(
        "const JUNCTION_RADIUS:f32={radius:.9};\nconst NAMED_FACES:array<u32,13> = array<u32,13>({});\n{}",
        ids.iter()
            .map(|id| format!("{id}u"))
            .collect::<Vec<_>>()
            .join(","),
        include_str!("ownership.wgsl")
    )
}
fn rim_at(p: [f64; 3], extent: f64, a: f64) -> f64 {
    let d = 1.0 - a.cos();
    let q = [d, p[1].clamp(d, extent - d), 9.0 + a.sin()];
    (0..3).map(|i| (p[i] - q[i]).powi(2)).sum()
}
fn rim(p: [f64; 3], extent: f64) -> f64 {
    let mut best = f64::INFINITY;
    let mut index = 0;
    let angle = std::f64::consts::FRAC_PI_2;
    for i in 0..=32 {
        let d = rim_at(p, extent, angle * i as f64 / 32.0);
        if d < best {
            best = d;
            index = i;
        }
    }
    let mut lo = angle * (index - 1).max(0) as f64 / 32.0;
    let mut hi = angle * (index + 1).min(32) as f64 / 32.0;
    for _ in 0..28 {
        let a = (2.0 * lo + hi) / 3.0;
        let b = (lo + 2.0 * hi) / 3.0;
        if rim_at(p, extent, a) < rim_at(p, extent, b) {
            hi = b;
        } else {
            lo = a;
        }
    }
    best.min(rim_at(p, extent, (lo + hi) * 0.5))
}
fn rect(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| (p[i] - p[i].clamp(lo[i], hi[i])).powi(2))
        .sum()
}
pub fn owner(p: [f64; 3], ids: &[u32; 13], radius: f64) -> u32 {
    let rho = ((p[0] - 20.0).powi(2) + (p[1] - 15.0).powi(2)).sqrt();
    let blend = 8.0 + radius;
    let mut top = [p[0].clamp(1.0, 39.0), p[1].clamp(1.0, 29.0)];
    if ((top[0] - 20.0).powi(2) + (top[1] - 15.0).powi(2)).sqrt() < blend {
        let direction = if rho > 0.0 {
            [(p[0] - 20.0) / rho, (p[1] - 15.0) / rho]
        } else {
            [1.0, 0.0]
        };
        top = [20.0 + blend * direction[0], 15.0 + blend * direction[1]];
    }
    let torus = [rho - blend, p[2] - (10.0 + radius)];
    let junction = if torus[0] <= 0.0 && torus[1] <= 0.0 {
        ((torus[0].powi(2) + torus[1].powi(2)).sqrt() - radius).powi(2)
    } else {
        ((torus[0] + radius).powi(2) + torus[1].powi(2))
            .min(torus[0].powi(2) + (torus[1] + radius).powi(2))
    };
    let distances = [
        rect(p, [0.0; 3], [40.0, 30.0, 0.0]),
        (p[0] - top[0]).powi(2) + (p[1] - top[1]).powi(2) + (p[2] - 10.0).powi(2),
        rect(p, [0.0; 3], [0.0, 30.0, 9.0]),
        rect(p, [40.0, 0.0, 0.0], [40.0, 30.0, 9.0]),
        rect(p, [0.0; 3], [40.0, 0.0, 9.0]),
        rect(p, [0.0, 30.0, 0.0], [40.0, 30.0, 9.0]),
        (rho - 8.0).max(0.0).powi(2) + (p[2] - 22.0).powi(2),
        (rho - 8.0).powi(2) + (p[2] - p[2].clamp(10.0 + radius, 22.0)).powi(2),
        junction,
        rim(p, 30.0),
        rim([40.0 - p[0], p[1], p[2]], 30.0),
        rim([p[1], p[0], p[2]], 40.0),
        rim([30.0 - p[1], p[0], p[2]], 40.0),
    ];
    let mut best = f64::INFINITY;
    let mut result = 0;
    for (i, d) in distances.into_iter().enumerate() {
        if d < best {
            best = d;
            result = ids[i];
        }
    }
    result
}
pub fn cpu(pixels: &[[f32; 4]], ids: &[u32; 13], radius: f64, workers: usize) -> Result<Vec<u32>> {
    let mut output = vec![0; pixels.len()];
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (points, values) in pixels
            .chunks(pixels.len().div_ceil(workers))
            .zip(output.chunks_mut(pixels.len().div_ceil(workers)))
        {
            handles.push(
                std::thread::Builder::new()
                    .spawn_scoped(scope, move || {
                        for (p, id) in points.iter().zip(values) {
                            if p[3] > 0.0 {
                                *id = owner([p[0] as f64, p[1] as f64, p[2] as f64], ids, radius);
                            }
                        }
                    })
                    .map_err(|e| Refusal::new("preview.worker-start", e.to_string()))?,
            );
        }
        for h in handles {
            h.join()
                .map_err(|_| Refusal::new("preview.worker", "selection worker panicked"))?;
        }
        Ok::<_, Refusal>(())
    })?;
    Ok(output)
}
pub fn sample_pixels(pixels: &[[f32; 4]], count: usize) -> Result<Vec<usize>> {
    let mut indices: Vec<_> = pixels
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (p[3] > 0.0).then_some(i))
        .collect();
    if indices.len() < count {
        return Err(Refusal::new(
            "preview.sample-count",
            format!(
                "only {} surface pixels for {count} required samples",
                indices.len()
            ),
        )
        .into());
    }
    let mut state = SEED;
    for i in 0..count {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = i + ((state % ((indices.len() - i) as u64)) as usize);
        indices.swap(i, j);
    }
    indices.truncate(count);
    indices.sort_unstable();
    Ok(indices)
}
pub fn check(
    plan: &Plan,
    prepared: &Prepared<'_>,
    pixels: &[[f32; 4]],
    ids: &[u32],
    samples: &[usize],
    workers: usize,
) -> Result<Value> {
    let mut rows = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for chunk in samples.chunks(samples.len().div_ceil(workers)) {
            handles.push(std::thread::Builder::new().spawn_scoped(scope,move||{
                chunk.iter().map(|&index|{
                    let p=pixels[index];let point=[p[0] as f64,p[1] as f64,p[2] as f64];
                    let [x,y,z]=point.map(|v|Rational::from_f64_exact(v).unwrap());
                    let (answer,work)=with_operation_budget(OWNERSHIP_ALLOWANCE,||prepared.nearest_face(&RPoint3::new(x,y,z).to_real(),&NearestOptions::default(),&mut Decider::default()));
                    let mut matches=false;
                    let reference=match answer {
                        Ok(Ownership::Proved(o))=>{let expected=plan.select(&o.origin).unwrap().index() as u32+1;matches=expected==ids[index];json!({"status":"proved","face_id":expected,"origin":o.origin})},
                        Ok(Ownership::Refused(a))=>json!({"status":"refused","code":a.refusal().code,"detail":a.refusal().detail,"cells":a.cells,"contenders":a.contenders.iter().map(|c|&c.origin).collect::<Vec<_>>()}),
                        Err(e)=>json!({"status":"refused","code":e.code,"detail":e.detail}),
                    };
                    let origin=ids[index].checked_sub(1).and_then(|i|plan.faces().get(i as usize)).map(|f|&f.origin);
                    json!({"pixel":[index%WIDTH,index/WIDTH],"index":index,"point_mm":point,"gpu_or_processor_preview":{"face_id":ids[index],"origin":origin},"f8":reference,"matches":matches,"work":work.used})
                }).collect::<Vec<_>>()
            }).map_err(|e|Refusal::new("preview.worker-start",e.to_string()))?);
        }
        let mut rows = Vec::with_capacity(samples.len());
        for h in handles {
            rows.extend(
                h.join()
                    .map_err(|_| Refusal::new("preview.worker", "F8 sample worker panicked"))?,
            );
        }
        Ok::<_, Refusal>(rows)
    })?;
    rows.sort_by_key(|r| r["index"].as_u64().unwrap());
    let matches = rows.iter().filter(|r| r["matches"] == true).count();
    let proved = rows
        .iter()
        .filter(|r| r["f8"]["status"] == "proved")
        .count();
    let mismatches: Vec<_> = rows
        .iter()
        .filter(|r| r["matches"] != true)
        .cloned()
        .collect();
    Ok(
        json!({"sample_count":samples.len(),"seed":SEED,"agreement_rate":matches as f64/samples.len() as f64,"matching_proved":matches,"proved":proved,"refused":samples.len()-proved,"mismatches":mismatches,"samples":rows}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Viewer;
    fn fixture() -> Plan {
        Plan::compile(
            &catallax_recipe::Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json"))
                .unwrap(),
        )
        .unwrap()
    }
    pub(super) fn patch_points(radius: f64) -> (Vec<[f32; 4]>, Vec<usize>) {
        let mut points = Vec::new();
        let mut slots = Vec::new();
        for i in 0..1040 {
            let slot = i % 13;
            let u = 0.05 + 0.9 * ((i * 37 % 997) as f64 / 997.0);
            let v = 0.05 + 0.9 * ((i * 73 % 991) as f64 / 991.0);
            let theta = u * std::f64::consts::FRAC_PI_2;
            let phi = v * std::f64::consts::TAU;
            let d = 1.0 - theta.cos();
            let z = 9.0 + theta.sin();
            let p = match slot {
                0 => [2.0 + 36.0 * u, 2.0 + 26.0 * v, 0.0],
                1 => [1.5, 2.0 + 26.0 * v, 10.0],
                2 => [0.0, 2.0 + 26.0 * u, 1.0 + 7.0 * v],
                3 => [40.0, 2.0 + 26.0 * u, 1.0 + 7.0 * v],
                4 => [2.0 + 36.0 * u, 0.0, 1.0 + 7.0 * v],
                5 => [2.0 + 36.0 * u, 30.0, 1.0 + 7.0 * v],
                6 => [20.0 + 7.0 * u * phi.cos(), 15.0 + 7.0 * u * phi.sin(), 22.0],
                7 => [
                    20.0 + 8.0 * phi.cos(),
                    15.0 + 8.0 * phi.sin(),
                    10.5 + radius + (11.0 - radius) * u,
                ],
                8 => {
                    let r = 8.0 + radius - radius * theta.cos();
                    [
                        20.0 + r * phi.cos(),
                        15.0 + r * phi.sin(),
                        10.0 + radius - radius * theta.sin(),
                    ]
                }
                9 => [d, 1.1 + 27.8 * v, z],
                10 => [40.0 - d, 1.1 + 27.8 * v, z],
                11 => [1.1 + 37.8 * v, d, z],
                _ => [1.1 + 37.8 * v, 30.0 - d, z],
            };
            points.push([p[0] as f32, p[1] as f32, p[2] as f32, 1.0]);
            slots.push(slot);
        }
        (points, slots)
    }
    #[test]
    fn patch_names_and_processor_workers_agree() {
        let plan = fixture();
        let ids = named_faces(&plan).unwrap();
        for radius in [2.0, 3.0] {
            let (points, slots) = patch_points(radius);
            let expected: Vec<_> = slots.iter().map(|s| ids[*s]).collect();
            for workers in [1, 2, 4] {
                assert!(cpu(&points, &ids, radius, workers).unwrap() == expected);
            }
        }
    }
    #[test]
    fn seeded_surface_sample_is_fixed_unique_and_large_enough() {
        let (points, _) = patch_points(2.0);
        let a = sample_pixels(&points, SAMPLE_COUNT).unwrap();
        assert_eq!(a, sample_pixels(&points, SAMPLE_COUNT).unwrap());
        assert!(a.windows(2).all(|w| w[0] < w[1]));
        assert!(a.len() >= 1000);
        assert!(sample_pixels(&points[..12], SAMPLE_COUNT).is_err());
    }
    #[test]
    fn exact_f8_sample_answers_and_work_match_all_worker_counts() {
        let plan = fixture();
        let prepared = plan.prepare_distances(&mut Decider::default()).unwrap();
        prepared.prepare_nearest(&mut Decider::default()).unwrap();
        let pixels = vec![[5.0, 6.0, -3.0, 1.0]; 4];
        let names = named_faces(&plan).unwrap();
        let ids = vec![names[0]; 4];
        let samples = vec![0, 1, 2, 3];
        let expected = check(&plan, &prepared, &pixels, &ids, &samples, 1).unwrap();
        assert_eq!(expected["agreement_rate"], 1.0);
        for workers in [2, 4] {
            assert_eq!(
                check(&plan, &prepared, &pixels, &ids, &samples, workers).unwrap(),
                expected
            );
        }
    }
    #[test]
    fn closed_patch_wgsl_validates_without_gpu() {
        let ids = named_faces(&fixture()).unwrap();
        for radius in [2.0, 3.0] {
            let source = source(&ids, radius);
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
    fn gpu_closed_patch_ids_match_independent_patch_controls() {
        let plan = fixture();
        let ids = named_faces(&plan).unwrap();
        let Ok(gpu) = crate::Gpu::new(false) else {
            eprintln!("GPU ownership check skipped: no adapter");
            return;
        };
        for radius in [2.0, 3.0] {
            let source = source(&ids, radius);
            let module = naga::front::wgsl::parse_str(&source)
                .unwrap_or_else(|e| panic!("{}", e.emit_to_string(&source)));
            naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::empty(),
            )
            .validate(&module)
            .unwrap();
            let (points, slots) = patch_points(radius);
            let actual = gpu.owners(&source, &points).unwrap();
            for (i, (a, s)) in actual.iter().zip(slots).enumerate() {
                assert_eq!(*a, ids[s], "patch sample {i}");
            }
        }
    }
}
