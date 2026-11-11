//! A concave round on the rings where a cylinder passes through a sphere.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate, evaluate_certified};
use catallax_measure::{Accuracy, Enclosure, VolumeAreaCentroid, volume_area_centroid};
use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::Value;

fn input() -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "sphere-cylinder")
        .unwrap()["recipe"]
        .clone()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 200_000_000,
        ..Default::default()
    })
}

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn on_thread(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn both_rings_round_with_certified_bands() {
    on_thread(|| {
        let recipe = Recipe::parse(&input().to_string()).unwrap();
        let certified = evaluate_certified(&recipe, &mut decider()).unwrap();
        let junction = certified.junction().unwrap();
        assert_eq!(junction.sphere_bands().len(), 2);
        assert!(junction.saddles().is_empty());
        for band in junction.sphere_bands() {
            assert!(band.bound().upper <= q("1/10000"));
            assert_ne!(band.rails()[0], band.rails()[1]);
        }
        let built = certified.built();
        let mut model = built.model.clone();
        assert!(model.check(built.result).is_empty());
        assert_eq!(model.faces_of_body(built.result).len(), 7);
        // The ring passes (3, 0, 4); the ball there has its centre at
        // (7/2, 0, sqrt(18)), so a point just off the sharp edge is new
        // material and the centre is outside.
        for (coordinates, expected) in [
            (["1", "0", "0"], Containment::Inside),
            (["1", "0", "9"], Containment::Inside),
            (["151/50", "0", "401/100"], Containment::Inside),
            (["151/50", "0", "-401/100"], Containment::Inside),
            (["7/2", "0", "17/4"], Containment::Outside),
            (["-26/5", "0", "0"], Containment::Outside),
        ] {
            let [x, y, z] = coordinates.map(q);
            let actual = classify_point(
                &mut model,
                built.result,
                &RPoint3::new(x, y, z),
                &mut decider(),
            );
            assert!(
                matches!(&actual, Outcome::Proven(found) if *found == expected),
                "at {coordinates:?}: {actual:?}"
            );
        }
        // Volume, area and centroid are measured; see the test below. A box and
        // a mesh still need the carrier itself, and each says so.
        let accuracy = Accuracy {
            relative: q("1/1000000"),
            absolute: q("1/1000000"),
        };
        let refusal =
            catallax_measure::bounding_box(&model, built.result, &accuracy, &mut decider())
                .unwrap_err();
        assert!(refusal.code.starts_with("measure."), "{refusal:?}");
        assert!(
            catallax_mesh::mesh_body(
                &model,
                built.result,
                &MeshOptions::default(),
                &mut decider()
            )
            .is_err()
        );
    });
}

// The fixture is symmetric in y = 0, so its mirror image there is the same
// solid with a left-handed band frame. Each retained face keeps its points.
fn mirrored() -> Value {
    let mut value = input();
    value["steps"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({
            "id": "flipped", "op": "mirror", "body": "round",
            "plane": {"point": ["0", "0", "0"], "normal": ["0", "1", "0"]},
            "keep_original": false
        }));
    value["result"] = "flipped".into();
    value
}

#[test]
fn a_mirrored_copy_keeps_its_caps() {
    on_thread(|| {
        let recipe = Recipe::parse(&mirrored().to_string()).unwrap();
        let built = evaluate(&recipe, &mut decider()).unwrap();
        let mut model = built.model.clone();
        assert!(model.check(built.result).is_empty());
        assert_eq!(model.faces_of_body(built.result).len(), 7);
        // (5, 0, 0) is on the retained sphere face, between the two rings,
        // so material ends there.
        for (coordinates, expected) in [
            (["49/10", "0", "0"], Containment::Inside),
            (["51/10", "0", "0"], Containment::Outside),
            (["-49/10", "0", "0"], Containment::Inside),
            (["0", "49/10", "0"], Containment::Inside),
            (["0", "-51/10", "0"], Containment::Outside),
            (["1", "0", "0"], Containment::Inside),
            (["1", "0", "9"], Containment::Inside),
            (["151/50", "0", "401/100"], Containment::Inside),
            (["151/50", "0", "-401/100"], Containment::Inside),
            (["7/2", "0", "17/4"], Containment::Outside),
            (["-26/5", "0", "0"], Containment::Outside),
        ] {
            let [x, y, z] = coordinates.map(q);
            let actual = classify_point(
                &mut model,
                built.result,
                &RPoint3::new(x, y, z),
                &mut decider(),
            );
            assert!(
                matches!(&actual, Outcome::Proven(found) if *found == expected),
                "at {coordinates:?}: {actual:?}"
            );
        }
    });
}

#[test]
fn a_cylinder_that_does_not_pass_through_keeps_its_refusal() {
    on_thread(|| {
        let mut value = input();
        // At a gap of 16/5 the radius-2 cylinder leaves the radius-5 sphere
        // through its side, so the edge is one loop round neither.
        value["steps"][1]["pose"]["translation"] = serde_json::json!(["16/5", "0", "-10"]);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let error = evaluate(&recipe, &mut decider()).unwrap_err();
        let BuildError::Geometry { refusal, .. } = error else {
            panic!("expected a construction refusal, got {error:?}");
        };
        assert_eq!(refusal.code, "fillet.unsupported-contact", "{refusal:?}");
    });
}

// The fixture's numbers: a radius-5 sphere at the origin, a radius-2 cylinder
// whose axis runs through (1, 0, z), and a ball of radius 1/2 outside both.
const SPHERE: f64 = 5.0;
const CYLINDER: f64 = 2.0;
const GAP: f64 = 1.0;
const BALL: f64 = 0.5;

// The ball centre at angle `t` round the cylinder axis, on the upper ring.
fn centre(t: f64) -> [f64; 3] {
    let (big, small) = (SPHERE + BALL, CYLINDER + BALL);
    let (x, y) = (GAP + small * t.cos(), small * t.sin());
    [x, y, (big * big - x * x - y * y).sqrt()]
}

// How far the balls' lowest point stands above the sphere, up the column at
// `(radial, angle)` round the cylinder axis, or `None` where no ball reaches
// down past the sphere. Read off the balls themselves, not off the band.
fn column(radial: f64, angle: f64) -> Option<f64> {
    let (x, y) = (GAP + radial * angle.cos(), radial * angle.sin());
    let floor = (SPHERE * SPHERE - x * x - y * y).sqrt();
    let reach = (BALL / (CYLINDER + BALL)).asin() * 1.5;
    let under = |t: f64| {
        let c = centre(t);
        let left = BALL * BALL - (x - c[0]).powi(2) - (y - c[1]).powi(2);
        (left > 0.0).then(|| c[2] - left.sqrt())
    };
    const SAMPLES: usize = 64;
    let mut best: Option<(f64, f64)> = None;
    for i in 0..=SAMPLES {
        let t = angle - reach + 2.0 * reach * i as f64 / SAMPLES as f64;
        if let Some(hit) = under(t)
            && best.is_none_or(|(_, low)| hit < low)
        {
            best = Some((t, hit));
        }
    }
    let (mut t, mut low) = best?;
    let mut step = 2.0 * reach / SAMPLES as f64;
    for _ in 0..30 {
        step /= 2.0;
        for next in [t - step, t + step] {
            if let Some(hit) = under(next)
                && hit < low
            {
                (t, low) = (next, hit);
            }
        }
    }
    (low > floor).then_some(low - floor)
}

// The upper pocket's volume and its first moment along x, by columns.
fn pocket() -> [f64; 2] {
    const ANGLES: usize = 720;
    const RADII: usize = 600;
    let width = BALL / RADII as f64;
    let turn = 2.0 * std::f64::consts::PI / ANGLES as f64;
    let (mut volume, mut along) = (0.0, 0.0);
    for a in 0..ANGLES {
        let angle = (a as f64 + 0.5) * turn;
        let mut last = f64::INFINITY;
        for r in 0..RADII {
            let radial = CYLINDER + (r as f64 + 0.5) * width;
            // The balls touch the sphere at the rail and stand off it again
            // beyond, where the gap is outside the pocket.
            let Some(height) = column(radial, angle).filter(|height| *height <= last + 1e-9) else {
                break;
            };
            last = height;
            let size = height * radial * width * turn;
            volume += size;
            along += size * (GAP + radial * angle.cos());
        }
    }
    [volume, along]
}

// The upper band's area less the strips of sphere and cylinder it covers, from
// a fine triangulation of the band and the strips' own integrals.
fn area_change() -> f64 {
    let (big, small) = (SPHERE + BALL, CYLINDER + BALL);
    let point = |t: f64, s: f64| {
        let c = centre(t);
        let first = c.map(|v| -v / big);
        let toward = [-t.cos(), -t.sin(), 0.0];
        let cosine: f64 = (0..3).map(|i| first[i] * toward[i]).sum();
        let angle = cosine.acos();
        let sine = angle.sin();
        // From the sphere ray at `s = 0` to the cylinder ray at `s = 1`.
        let (a, b) = (((1.0 - s) * angle).sin() / sine, (s * angle).sin() / sine);
        std::array::from_fn::<f64, 3, _>(|i| c[i] + BALL * (a * first[i] + b * toward[i]))
    };
    const TURNS: usize = 4000;
    const ACROSS: usize = 200;
    let mut band = 0.0;
    for i in 0..TURNS {
        let [t0, t1] = [i, i + 1].map(|k| 2.0 * std::f64::consts::PI * k as f64 / TURNS as f64);
        for j in 0..ACROSS {
            let [s0, s1] = [j, j + 1].map(|k| k as f64 / ACROSS as f64);
            let corners = [point(t0, s0), point(t1, s0), point(t1, s1), point(t0, s1)];
            for [p, q, r] in [[0, 1, 2], [0, 2, 3]].map(|k| k.map(|i| corners[i])) {
                let u: [f64; 3] = std::array::from_fn(|i| q[i] - p[i]);
                let v: [f64; 3] = std::array::from_fn(|i| r[i] - p[i]);
                let n = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                band += (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt() / 2.0;
            }
        }
    }
    // The sphere's area over the shadow between the cylinder's circle and the
    // rail's, by polar cells round the sphere's centre.
    let scale = SPHERE / big;
    let (rail_x, rail_r) = (GAP * scale, small * scale);
    let mut sphere = 0.0;
    const RINGS: usize = 3000;
    const SPOKES: usize = 3000;
    let outer = rail_x + rail_r;
    for i in 0..RINGS {
        let rho = (i as f64 + 0.5) * outer / RINGS as f64;
        for j in 0..SPOKES {
            let phi = (j as f64 + 0.5) * 2.0 * std::f64::consts::PI / SPOKES as f64;
            let (x, y) = (rho * phi.cos(), rho * phi.sin());
            let off_cylinder = (x - GAP).powi(2) + y * y > CYLINDER * CYLINDER;
            let in_rail = (x - rail_x).powi(2) + y * y < rail_r * rail_r;
            if off_cylinder && in_rail {
                let cell =
                    rho * (outer / RINGS as f64) * (2.0 * std::f64::consts::PI / SPOKES as f64);
                sphere += cell * SPHERE / (SPHERE * SPHERE - rho * rho).sqrt();
            }
        }
    }
    // The cylinder between the ring and the rail, at the ball centre's height.
    let mut cylinder = 0.0;
    for i in 0..TURNS {
        let t = (i as f64 + 0.5) * 2.0 * std::f64::consts::PI / TURNS as f64;
        let (x, y) = (GAP + CYLINDER * t.cos(), CYLINDER * t.sin());
        let ring = (SPHERE * SPHERE - x * x - y * y).sqrt();
        cylinder += CYLINDER * (centre(t)[2] - ring) * 2.0 * std::f64::consts::PI / TURNS as f64;
    }
    band - sphere - cylinder
}

#[test]
fn both_rings_measure_against_their_columns() {
    on_thread(|| {
        let recipe = Recipe::parse(&input().to_string()).unwrap();
        let built = evaluate(&recipe, &mut decider()).unwrap();
        let base = built
            .model
            .blend_base(built.result)
            .expect("a recorded base");
        let accuracy = Accuracy {
            relative: q("1/1000000"),
            absolute: q("1/1000"),
        };
        let measure = |body| {
            volume_area_centroid(&built.model, body, &accuracy, &mut decider())
                .unwrap()
                .value
        };
        let [rounded, plain] = [measure(built.result), measure(base)];
        let mid = |value: &Enclosure| ((&value.lo + &value.hi) / q("2")).to_f64_nearest();
        let width = |value: &Enclosure| (&value.hi - &value.lo).to_f64_nearest();
        // Each ring's ball fills a concave corner, so both add material. The
        // lower ring mirrors the upper in z = 0.
        let [volume, along] = pocket();
        let change = mid(&rounded.volume) - mid(&plain.volume);
        assert!(volume > 0.1, "{volume}");
        assert!(
            (change - 2.0 * volume).abs() < 2.0 * volume * 1e-3,
            "{change} {volume}"
        );
        // The body's volume holds the carrier as well as the definition, so its
        // width is the stated bound the fitted bands add, not quadrature: here
        // under a thousandth of a cubic millimetre for bands within 1e-4 mm. The
        // definition alone measures far tighter.
        let stated = width(&rounded.volume);
        assert!(stated < 2e-3, "{stated}");
        assert!(
            (change - 2.0 * volume).abs() < stated / 2.0 + 1e-4,
            "{change} {volume} {stated}"
        );
        let definition = catallax_measure::blend_definition_properties(
            &built.model,
            built.result,
            &accuracy,
            &mut decider(),
        )
        .unwrap()
        .value;
        assert!(width(&definition.volume) < volume * 1e-3);
        assert!((mid(&definition.volume) - mid(&rounded.volume)).abs() < stated);
        let moment = |value: &VolumeAreaCentroid| -> [f64; 3] {
            let centroid = value.centroid.as_ref().unwrap();
            std::array::from_fn(|i| mid(&centroid[i]) * mid(&value.volume))
        };
        let [a, b] = [moment(&rounded), moment(&plain)];
        let target = [2.0 * along, 0.0, 0.0];
        for i in 0..3 {
            assert!(
                (a[i] - b[i] - target[i]).abs() < 2e-3 * (volume + along.abs()),
                "{i}: {} {}",
                a[i] - b[i],
                target[i]
            );
        }
        let area = mid(&rounded.area) - mid(&plain.area);
        let target = 2.0 * area_change();
        assert!((area - target).abs() < 2e-3, "{area} {target}");
        // A mirrored copy reads its left-handed frame to the same measures.
        let recipe = Recipe::parse(&mirrored().to_string()).unwrap();
        let flipped = evaluate(&recipe, &mut decider()).unwrap();
        let copy = volume_area_centroid(&flipped.model, flipped.result, &accuracy, &mut decider())
            .unwrap()
            .value;
        assert!((mid(&copy.volume) - mid(&rounded.volume)).abs() < 1e-6);
        assert!((mid(&copy.area) - mid(&rounded.area)).abs() < 1e-6);
        let [a, b] = [&copy, &rounded].map(|value| value.centroid.as_ref().unwrap());
        for i in 0..3 {
            let sign = if i == 1 { -1.0 } else { 1.0 };
            assert!((mid(&a[i]) - sign * mid(&b[i])).abs() < 1e-6, "{i}");
        }
        eprintln!(
            "volume change {change:.6} columns {:.6} stated {stated:.2e}; area change {area:.6} grid {target:.6}",
            2.0 * volume
        );
    });
}
