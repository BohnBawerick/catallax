//! Rolling-ball fillets on the ring where two nonparallel cylinders cross.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, certify_junction, evaluate, evaluate_incremental};
use catallax_field::Plan;
use catallax_measure::{
    Accuracy, Enclosure, bounding_box, edge_length, face_area, face_bounding_box, mass_properties,
    volume, volume_area_centroid,
};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, with_operation_log};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{EnvelopeGrid, RPoint3, SplineSupport, Surface};
use catallax_step::write_step;
use catallax_topology::{BodyId, FaceId};
use serde_json::{Value, json};

const LEVEL: [&str; 4] = ["1", "0", "0", "0"];
/// A quarter turn about x, which lays a cylinder along y.
const FLAT: [&str; 4] = ["1", "1", "0", "0"];
/// Leans a cylinder's axis to `(0, 3/5, 4/5)`.
const LEAN: [&str; 4] = ["3", "-1", "0", "0"];

fn cylinder(id: &str, radius: &str, height: &str, rotation: [&str; 4], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"cylinder","radius":radius,"height":height,
        "pose":{"rotation":rotation,"translation":at}})
}
fn wide() -> Value {
    cylinder("wide", "10", "60", FLAT, ["0", "30", "0"])
}
fn wall() -> Value {
    cylinder("wall", "4", "30", LEVEL, ["1", "0", "0"])
}
fn ring() -> Value {
    json!([{"between":[{"step":"wide","role":{"side":{"segment":"wall"}}},
        {"step":"wall","role":{"side":{"segment":"wall"}}}]}])
}
fn parse(steps: Vec<Value>, result: &str) -> Recipe {
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"ring","steps":steps,"result":result})
            .to_string(),
    )
    .unwrap()
}
/// The two cylinders joined, and rounded when a radius is given.
fn recipe(parts: [Value; 2], kind: &str, order: [&str; 2], radius: Option<&str>) -> Recipe {
    let [first, second] = parts;
    let mut steps = vec![
        first,
        second,
        json!({"id":"joined","op":"boolean","kind":kind,"target":order[0],"tool":order[1]}),
    ];
    if let Some(radius) = radius {
        steps.push(
            json!({"id":"round","op":"fillet","body":"joined","radius":radius,
            "edges":ring()}),
        );
    }
    parse(steps, if radius.is_some() { "round" } else { "joined" })
}
fn boss(radius: &str) -> Recipe {
    recipe([wide(), wall()], "union", ["wide", "wall"], Some(radius))
}
fn bore(radius: &str) -> Recipe {
    recipe([wide(), wall()], "subtract", ["wide", "wall"], Some(radius))
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut decider()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    built
}
fn classify(built: &mut Built, points: &[[&str; 3]], want: Containment) {
    for p in points {
        let [x, y, z] = p.map(q);
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &RPoint3::new(x, y, z),
            &mut decider(),
        );
        assert!(
            matches!(actual, Outcome::Proven(found) if found == want),
            "at {p:?}, wanted {want:?}: {actual:?}"
        );
    }
}
fn bands(built: &Built) -> Vec<FaceId> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|face| {
            matches!(
                built.model.geom().surface(built.model.face(*face).surface),
                Surface::BSpline(SplineSupport::RollingSaddle(_))
            )
        })
        .collect()
}
fn fine(built: &Built) -> Mesh {
    let options = MeshOptions {
        chord_tolerance: q("1/4"),
        max_vertices: 4_000_000,
        ..Default::default()
    };
    mesh_body(&built.model, built.result, &options, &mut decider())
        .unwrap()
        .value
}

// Every mesh edge is used once each way, so the mesh is closed and consistently
// wound, and the positive volume says the winding faces out.
fn assert_closed(mesh: &Mesh) {
    let mut uses: BTreeMap<(usize, usize), i32> = BTreeMap::new();
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices;
        for (from, to) in [(a, b), (b, c), (c, a)] {
            let sign = if from < to { 1 } else { -1 };
            *uses.entry((from.min(to), from.max(to))).or_default() += sign;
        }
    }
    let open = uses.values().filter(|count| **count != 0).count();
    assert_eq!(open, 0, "{open} of {} mesh edges are unmatched", uses.len());
    let size: f64 = mesh
        .triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum();
    assert!(size > 0.0, "the mesh is wound inward: {size}");
}

/// Two cylinders in the frame of their common perpendicular, written out here
/// from the rolling ball's two distance conditions and sharing no code with the
/// kernel. `across` runs from the wide axis to the narrow one, `up` is the
/// narrow axis, pointing out of the wide cylinder at this ring, and `round` with
/// `across` and `up` makes a right-handed frame about the narrow axis.
struct Setting {
    wide: f64,
    narrow: f64,
    ball: f64,
    gap: f64,
    cos: f64,
    sin: f64,
    across: [f64; 3],
    round: [f64; 3],
    up: [f64; 3],
    foot: [f64; 3],
}
fn dot(a: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * c[0] + a[1] * c[1] + a[2] * c[2]
}
impl Setting {
    // Heights along the narrow axis of the wide cylinder.
    fn floor(&self, radial: f64, angle: f64) -> f64 {
        let x = self.gap + radial * angle.cos();
        let y = (self.wide * self.wide - x * x).sqrt();
        (y - radial * self.cos * angle.sin()) / self.sin
    }
    // A ball centre at angle `t` round the narrow axis.
    fn centre(&self, t: f64, side: f64) -> (f64, f64) {
        let [big, small] = [self.wide + side * self.ball, self.narrow + self.ball];
        let x = self.gap + small * t.cos();
        let y = (big * big - x * x).sqrt();
        (small, (y - small * self.cos * t.sin()) / self.sin)
    }
    // How far the ball set stands off the wide cylinder, up the line at
    // `(radial, angle)`. Above the wide cylinder a ball of an added fillet
    // leaves material under it. Below it, a ball of a removed one takes the
    // material above it.
    fn column(&self, radial: f64, angle: f64, side: f64) -> f64 {
        let meets = |t: f64| {
            let (small, height) = self.centre(t, side);
            let (dx, dy) = (
                radial * angle.cos() - small * t.cos(),
                radial * angle.sin() - small * t.sin(),
            );
            let left = self.ball * self.ball - dx * dx - dy * dy;
            (left > 0.0).then(|| height - side * left.sqrt())
        };
        let reach = (self.ball / (self.narrow + self.ball)).asin();
        let beats = |a: f64, c: f64| if side > 0.0 { a < c } else { a > c };
        let mut best: Option<(f64, f64)> = None;
        const SAMPLES: usize = 64;
        for i in 0..=SAMPLES {
            let t = angle - reach + 2.0 * reach * i as f64 / SAMPLES as f64;
            if let Some(hit) = meets(t)
                && best.is_none_or(|(_, value)| beats(hit, value))
            {
                best = Some((t, hit));
            }
        }
        let Some((mut t, mut value)) = best else {
            return 0.0;
        };
        // Close in on the extreme by a shrinking local search.
        let mut step = 2.0 * reach / SAMPLES as f64;
        for _ in 0..24 {
            step /= 2.0;
            for candidate in [t - step, t + step] {
                if let Some(hit) = meets(candidate)
                    && beats(hit, value)
                {
                    (t, value) = (candidate, hit);
                }
            }
        }
        (side * (value - self.floor(radial, angle))).max(0.0)
    }
    // The volume between the wide cylinder and the rolling ball, by midpoint
    // sections round the narrow axis. Each section is a set of columns.
    // The ball touches the wide cylinder where a column falls to zero. Past that
    // point a column stands in open space, so each section stops there.
    fn section_integral(&self, side: f64) -> f64 {
        self.pocket(side)[0]
    }
    // The pocket's volume, then its first moments about the world origin.
    fn pocket(&self, side: f64) -> [f64; 4] {
        const ANGLES: usize = 720;
        const RADII: usize = 480;
        let width = 2.0 * self.ball / RADII as f64;
        let turn = 2.0 * PI / ANGLES as f64;
        let (mut sum, mut along, mut round, mut up) = (0.0, 0.0, 0.0, 0.0);
        for a in 0..ANGLES {
            let angle = (a as f64 + 0.5) * turn;
            let mut last = f64::INFINITY;
            for r in 0..RADII {
                let radial = self.narrow + (r as f64 + 0.5) * width;
                let height = self.column(radial, angle, side);
                if height > last + 1e-9 {
                    break;
                }
                last = height;
                let size = height * radial * width * turn;
                sum += size;
                along += size * radial * angle.cos();
                round += size * radial * angle.sin();
                // A column runs from the floor up by `side * height`.
                up += size * (self.floor(radial, angle) + side * height / 2.0);
            }
        }
        let moment: [f64; 3] = std::array::from_fn(|i| {
            self.foot[i] * sum + self.across[i] * along + self.round[i] * round + self.up[i] * up
        });
        [sum, moment[0], moment[1], moment[2]]
    }
    // The area the band adds less the two strips it covers, read off its grid: the
    // facets of the band, and each rail against the ring on its cylinder's flat
    // chart. The ring comes from the two cylinders alone.
    fn area_change(&self, grid: &EnvelopeGrid) -> f64 {
        let at = |column: usize, line: usize| {
            let p = &grid.points[column * grid.rows + line];
            let d: [f64; 3] = std::array::from_fn(|k| p[k].mid().to_f64_nearest() - self.foot[k]);
            [dot(d, self.across), dot(d, self.round), dot(d, self.up)]
        };
        let mut band = 0.0;
        for column in 0..grid.columns - 1 {
            for line in 0..grid.rows - 1 {
                let corners = [
                    at(column, line),
                    at(column + 1, line),
                    at(column + 1, line + 1),
                    at(column, line + 1),
                ];
                for [a, b, c] in [
                    [corners[0], corners[1], corners[2]],
                    [corners[0], corners[2], corners[3]],
                ] {
                    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                    let n = [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ];
                    band += dot(n, n).sqrt() / 2.0;
                }
            }
        }
        // Angle round the wide axis, and height along it.
        let chart = |[u1, u2, h]: [f64; 3]| {
            let x = self.gap + u1;
            let y = self.cos * u2 + self.sin * h;
            let z = -self.sin * u2 + self.cos * h;
            [x.atan2(y), z]
        };
        let shoelace = |points: &[[f64; 2]]| {
            let count = points.len();
            (0..count)
                .map(|i| {
                    let (a, b) = (points[i], points[(i + 1) % count]);
                    a[0] * b[1] - a[1] * b[0]
                })
                .sum::<f64>()
                / 2.0
        };
        let on_rail: Vec<_> = (0..grid.columns - 1).map(|c| chart(at(c, 0))).collect();
        const STEPS: usize = 20_000;
        let ring: Vec<[f64; 3]> = (0..STEPS)
            .map(|i| {
                let angle = 2.0 * PI * i as f64 / STEPS as f64;
                let radial = self.narrow;
                [
                    radial * angle.cos(),
                    radial * angle.sin(),
                    self.floor(radial, angle),
                ]
            })
            .collect();
        let on_ring: Vec<_> = ring.iter().map(|p| chart(*p)).collect();
        let wide = self.wide * (shoelace(&on_rail).abs() - shoelace(&on_ring).abs()).abs();
        let mut narrow = 0.0;
        let count = grid.columns - 1;
        for c in 0..count {
            let [a, b] = [at(c, grid.rows - 1), at((c + 1) % count, grid.rows - 1)];
            let [from, to] = [a[1].atan2(a[0]), b[1].atan2(b[0])];
            let turn = (to - from + 3.0 * PI).rem_euclid(2.0 * PI) - PI;
            let middle = from + turn / 2.0;
            let rail = (a[2] + b[2]) / 2.0;
            narrow += (rail - self.floor(self.narrow, middle)).abs() * turn.abs() * self.narrow;
        }
        band - wide - narrow
    }
    // The same volume read off the fitted band alone: each triangle of its grid
    // casts a shadow on the narrow axis's plane and stands over the wide cylinder.
    fn band_integral(&self, grid: &EnvelopeGrid) -> f64 {
        let at = |column: usize, line: usize| {
            let p = &grid.points[column * grid.rows + line];
            let d: [f64; 3] = std::array::from_fn(|k| p[k].mid().to_f64_nearest() - self.foot[k]);
            [dot(d, self.across), dot(d, self.round), dot(d, self.up)]
        };
        let mut sum = 0.0;
        for column in 0..grid.columns - 1 {
            for line in 0..grid.rows - 1 {
                let corners = [
                    at(column, line),
                    at(column + 1, line),
                    at(column + 1, line + 1),
                    at(column, line + 1),
                ];
                for [a, b, c] in [
                    [corners[0], corners[1], corners[2]],
                    [corners[0], corners[2], corners[3]],
                ] {
                    let area =
                        ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.0;
                    let [x, y, h] = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.0);
                    sum += area * (h - self.floor(x.hypot(y), y.atan2(x))).abs();
                }
            }
        }
        sum
    }
}

// The square tee: wide along y, narrow along z through (1, 0, 0).
fn square(ball: f64) -> Setting {
    Setting {
        wide: 10.0,
        narrow: 4.0,
        ball,
        gap: 1.0,
        cos: 0.0,
        sin: 1.0,
        across: [1.0, 0.0, 0.0],
        round: [0.0, 1.0, 0.0],
        up: [0.0, 0.0, 1.0],
        foot: [1.0, 0.0, 0.0],
    }
}
// The narrow axis leans to (0, 3/5, 4/5) through (1, 0, 0). With the wide axis
// taken as (0, 1, 0), the common perpendicular is along -x.
fn leaning(ball: f64) -> Setting {
    Setting {
        wide: 10.0,
        narrow: 4.0,
        ball,
        gap: -1.0,
        cos: 0.6,
        sin: 0.8,
        across: [-1.0, 0.0, 0.0],
        round: [0.0, -0.8, 0.6],
        up: [0.0, 0.6, 0.8],
        foot: [1.0, 0.0, 0.0],
    }
}

// The band's chord is 1/1000, so the volume over a band of about 40 mm2 cannot
// be out by more than 0.04 mm3 from the mesh alone, and it is far closer.
// The solid's volume and centroid against its base and the pocket's sections.
fn assert_measures(built: &Built, base: BodyId, setting: &Setting, side: f64) {
    let [want, moment @ ..] = setting.pocket(side);
    let accuracy = Accuracy {
        relative: q("1/1000000"),
        absolute: q("1/1000"),
    };
    let measure = |body| {
        volume_area_centroid(&built.model, body, &accuracy, &mut decider())
            .unwrap()
            .value
    };
    let start = std::time::Instant::now();
    let plain = measure(base);
    let base_time = start.elapsed();
    let rounded = measure(built.result);
    eprintln!(
        "measurement: base {base_time:?}, rounded {:?} with its base",
        start.elapsed() - base_time
    );
    let mid = |value: &Enclosure| ((&value.lo + &value.hi) / q("2")).to_f64_nearest();
    // The ball removes material inside a cylinder it rolls in, and adds outside.
    let change = mid(&rounded.volume) - mid(&plain.volume);
    assert!(
        (change - side * want).abs() < want * 1e-3,
        "{change} {want}"
    );
    let of_body = |value: &catallax_measure::VolumeAreaCentroid| -> [f64; 3] {
        let centroid = value.centroid.as_ref().unwrap();
        std::array::from_fn(|i| mid(&centroid[i]) * mid(&value.volume))
    };
    let [a, b] = [of_body(&rounded), of_body(&plain)];
    for i in 0..3 {
        let found = a[i] - b[i];
        let target = side * moment[i];
        assert!(
            (found - target).abs() < 2e-3 * (want + moment[i].abs()),
            "{i}: {found} {target}"
        );
    }
    eprintln!("change {change:.5}, sections {:.5}", side * want);
    let [band] = bands(built)[..] else {
        panic!("one band");
    };
    let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) =
        built.model.geom().surface(built.model.face(band).surface)
    else {
        panic!("a saddle");
    };
    let grid = saddle
        .display_grid(&q("1/1000"), 4_000_000, &decider())
        .unwrap();
    let area = mid(&rounded.area) - mid(&plain.area);
    let target = setting.area_change(&grid);
    eprintln!("area change {area:.5}, grid {target:.5}");
    assert!((area - target).abs() < 2e-2, "{area} {target}");
}

fn assert_volume(built: &Built, setting: &Setting, side: f64) {
    let [band] = bands(built)[..] else {
        panic!("one band");
    };
    assert_closed(&fine(built));
    let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) =
        built.model.geom().surface(built.model.face(band).surface)
    else {
        panic!("a saddle");
    };
    let grid = saddle
        .display_grid(&q("1/1000"), 4_000_000, &decider())
        .unwrap();
    let want = setting.section_integral(side);
    let found = setting.band_integral(&grid);
    eprintln!("band {found:.5}, sections {want:.5}");
    assert!(
        (found - want).abs() < want * 1e-3,
        "band {found}, sections {want}"
    );
}

#[test]
fn a_boss_gains_exactly_its_band() {
    let mut built = build(&boss("1"));
    assert_eq!(built.model.faces_of_body(built.result).len(), 6);
    assert_eq!(bands(&built).len(), 1);
    assert_eq!(built.fillets.len(), 1);
    classify(
        &mut built,
        // Round the ring, in the corner: added.
        &[
            ["1", "81/20", "10"],
            ["1", "-81/20", "10"],
            ["101/20", "0", "87/10"],
            ["-61/20", "0", "48/5"],
        ],
        Containment::Inside,
    );
    classify(
        &mut built,
        // The ball's own centre at those stations, and well out.
        &[
            ["1", "5", "219/20"],
            ["1", "-5", "219/20"],
            ["6", "0", "46/5"],
            ["-4", "0", "41/4"],
            ["20", "0", "0"],
        ],
        Containment::Outside,
    );
    classify(
        &mut built,
        &[["1", "0", "20"], ["0", "20", "0"]],
        Containment::Inside,
    );
    assert_volume(&built, &square(1.0), 1.0);
}

#[test]
fn a_bore_loses_exactly_its_lip() {
    let mut built = build(&bore("1"));
    assert_eq!(built.model.faces_of_body(built.result).len(), 6);
    classify(
        &mut built,
        // At the mouth, just inside the wide cylinder and clear of the bore:
        // removed.
        &[
            ["1", "81/20", "99/10"],
            ["1", "-81/20", "99/10"],
            ["101/20", "0", "43/5"],
        ],
        Containment::Outside,
    );
    classify(
        &mut built,
        // The ball's centres stay in the material.
        &[["1", "5", "89/10"], ["6", "0", "71/10"], ["-6", "0", "0"]],
        Containment::Inside,
    );
    assert_volume(&built, &square(1.0), -1.0);
}

#[test]
fn a_larger_ball_and_either_face_order_round_the_same_ring() {
    let built = build(&recipe(
        [wide(), wall()],
        "union",
        ["wall", "wide"],
        Some("2"),
    ));
    assert_volume(&built, &square(2.0), 1.0);
}

#[test]
fn a_leaning_boss_matches_its_sections() {
    let lean = cylinder("wall", "4", "30", LEAN, ["1", "0", "0"]);
    let built = build(&recipe(
        [wide(), lean],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_volume(&built, &leaning(1.0), 1.0);
}

#[test]
fn a_through_bore_rounds_both_mouths() {
    let hole = cylinder("wall", "4", "60", LEVEL, ["1", "0", "-30"]);
    let mut built = build(&recipe(
        [wide(), hole],
        "subtract",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(bands(&built).len(), 2);
    assert_eq!(built.model.faces_of_body(built.result).len(), 6);
    assert_closed(&fine(&built));
    // Both mouths lose their lip, and the ball's centres stay in the material.
    classify(
        &mut built,
        &[["1", "81/20", "99/10"], ["1", "81/20", "-99/10"]],
        Containment::Outside,
    );
    classify(
        &mut built,
        &[["1", "5", "89/10"], ["1", "5", "-89/10"]],
        Containment::Inside,
    );
}

#[test]
fn a_cached_ring_and_loop_record_their_phases() {
    for (name, parts, kind) in [
        ("ring", [wide(), wall()], "union"),
        ("loop boss", [wide(), past("8")], "union"),
        ("loop notch", [wide(), past("8")], "subtract"),
    ] {
        let order = ["wide", "wall"];
        let joined = recipe(parts.clone(), kind, order, None);
        let prefix = evaluate_incremental(&joined, None, &mut decider()).unwrap();
        let mut dec = decider();
        let rounded = recipe(parts, kind, order, Some("1"));
        let (result, phases) =
            with_operation_log(|| evaluate_incremental(&rounded, Some(prefix), &mut dec));
        let built = result.unwrap().into_built();
        for phase in ["fillet.contact", "fillet.fitting", "fillet.trimming"] {
            let measured: Vec<_> = phases.iter().filter(|p| p.operation == phase).collect();
            assert!(!measured.is_empty(), "{phase} was not recorded");
            let total: f64 = measured.iter().map(|p| p.wall_ms).sum();
            eprintln!("{name}: {phase}={total:.3}ms over {} spans", measured.len());
        }
        let start = std::time::Instant::now();
        assert_eq!(built.model.check(built.result), vec![]);
        eprintln!("{name}: certification={:?}", start.elapsed());
        let start = std::time::Instant::now();
        fine(&built);
        eprintln!("{name}: display={:?}", start.elapsed());
    }
}

#[test]
fn a_later_step_reuses_the_fit() {
    let rounded = evaluate_incremental(&boss("1"), None, &mut decider()).unwrap();
    let next = parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
            cylinder("far", "2", "4", LEVEL, ["50", "0", "0"]),
        ],
        "round",
    );
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&next, Some(rounded), &mut dec));
    let built = result.unwrap().into_built();
    assert_eq!(bands(&built).len(), 1);
    assert!(phases.iter().all(|p| p.operation != "fillet.fitting"));
}

#[test]
fn a_radius_edit_keeps_the_blend_name() {
    let names = |radius: &str| {
        let rounded = build(&boss(radius));
        let mut found: Vec<_> = rounded
            .model
            .faces_of_body(rounded.result)
            .into_iter()
            .flat_map(|f| rounded.model.face_origins(f).to_vec())
            .filter(|o| matches!(o.role, Role::Blend { .. }))
            .collect();
        found.sort();
        found
    };
    let first = names("1");
    assert_eq!(first.len(), 1);
    assert_eq!(first, names("3/2"));
}

#[test]
fn a_mirrored_copy_keeps_its_band() {
    let rounded = |last: Value, result: &str| {
        let mut steps = vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
        ];
        steps.push(last);
        build(&parse(steps, result))
    };
    let mut flipped = rounded(
        json!({"id":"flipped","op":"mirror","body":"round",
            "plane":{"point":["20","0","0"],"normal":["3","4","0"]},"keep_original":false}),
        "flipped",
    );
    assert_eq!(bands(&flipped).len(), 1);
    assert_closed(&fine(&flipped));
    // The boss's witnesses, reflected in that plane.
    classify(
        &mut flipped,
        &[["5396/500", "8553/500", "10"]],
        Containment::Inside,
    );
    classify(
        &mut flipped,
        &[["247/25", "421/25", "219/20"]],
        Containment::Outside,
    );
}

fn refusal(recipe: &Recipe) -> (String, String) {
    match evaluate(recipe, &mut decider()) {
        Err(BuildError::Geometry { refusal, .. }) => {
            (refusal.code.to_string(), refusal.detail.to_string())
        }
        Err(other) => panic!("expected a geometric refusal: {other}"),
        Ok(_) => panic!("expected a refusal, and it built"),
    }
}

#[test]
fn a_ball_that_does_not_fit_refuses() {
    // The bore's offset reaches the wide offset.
    assert_eq!(refusal(&bore("3")).0, "fillet.unsupported-contact");
    // The wide offset closes.
    assert_eq!(refusal(&bore("10")).0, "fillet.unsupported-contact");
}

#[test]
fn a_band_past_the_end_of_the_narrow_cylinder_refuses() {
    let short = cylinder("wall", "4", "11", LEVEL, ["1", "0", "0"]);
    let (code, _) = refusal(&recipe(
        [wide(), short],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(code, "fillet.unsupported-overflow");
}

#[test]
fn a_band_past_the_end_of_the_wide_cylinder_refuses() {
    let short = cylinder("wide", "10", "12", FLAT, ["0", "6", "0"]);
    let (code, _) = refusal(&recipe(
        [short, wall()],
        "union",
        ["wide", "wall"],
        Some("3"),
    ));
    assert_eq!(code, "fillet.unsupported-overflow");
}

// Equal radii cross in two ellipses, whose ball folds at their meeting points.
#[test]
fn equal_radii_refuse() {
    let same = cylinder("wall", "10", "30", LEVEL, ["0", "0", "0"]);
    let (code, _) = refusal(&recipe(
        [wide(), same],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(code, "fillet.unsupported-contact");
}

// A narrow cylinder that stops inside the wide one leaves open edges, not a ring.
#[test]
fn an_open_overlap_refuses() {
    let side = cylinder("wall", "4", "30", LEVEL, ["8", "0", "0"]);
    let (code, _) = refusal(&recipe(
        [wide(), side],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(code, "fillet.unsupported-topology");
}

fn with_neighbour(x: &str, radius: &str) -> Recipe {
    parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            cylinder("next", radius, "30", LEVEL, [x, "0", "0"]),
            json!({"id":"both","op":"boolean","kind":"union","target":"joined","tool":"next"}),
            json!({"id":"round","op":"fillet","body":"both","radius":"1","edges":ring()}),
        ],
        "round",
    )
}

// A second boss in the way is a curved face nothing clears the band of.
#[test]
fn a_curved_neighbour_in_the_way_refuses() {
    assert_eq!(
        refusal(&with_neighbour("7", "7/4")).0,
        "fillet.unsupported-overflow"
    );
}

// The band's contact on the wide cylinder stops short of x = 11/2 here, and its
// hull proves that.
#[test]
fn a_curved_neighbour_clear_of_the_band_stays() {
    let mut built = build(&with_neighbour("15/2", "2"));
    assert_eq!(bands(&built).len(), 1);
    assert_closed(&fine(&built));
    classify(&mut built, &[["101/20", "0", "87/10"]], Containment::Inside);
}

fn measured(recipe: &Recipe, setting: &Setting, side: f64) {
    let built = build(recipe);
    let base = built.model.blend_base(built.result).unwrap();
    assert_measures(&built, base, setting, side);
}

#[test]
fn a_boss_measures_as_its_base_and_sections() {
    measured(&boss("1"), &square(1.0), 1.0);
}

#[test]
fn a_bore_measures_as_its_base_and_sections() {
    measured(&bore("1"), &square(1.0), -1.0);
}

#[test]
fn a_leaning_boss_measures_as_its_base_and_sections() {
    let lean = cylinder("wall", "4", "30", LEAN, ["1", "0", "0"]);
    measured(
        &recipe([wide(), lean], "union", ["wide", "wall"], Some("1")),
        &leaning(1.0),
        1.0,
    );
}

#[test]
fn a_mirrored_copy_measures_its_volume() {
    let original = build(&boss("1"));
    let flipped = build(&parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
            json!({"id":"flipped","op":"mirror","body":"round",
                "plane":{"point":["20","0","0"],"normal":["3","4","0"]},"keep_original":false}),
        ],
        "flipped",
    ));
    let accuracy = Accuracy {
        relative: q("1/1000000"),
        absolute: q("1/1000"),
    };
    let measure = |built: &Built| {
        volume(&built.model, built.result, &accuracy, &mut decider())
            .unwrap()
            .value
    };
    let [a, c] = [measure(&original), measure(&flipped)];
    assert!(a.lo <= c.hi && c.lo <= a.hi, "{a:?} {c:?}");
}

#[test]
fn a_through_bore_changes_by_both_mouths() {
    let hole = cylinder("wall", "4", "60", LEVEL, ["1", "0", "-30"]);
    let built = build(&recipe(
        [wide(), hole],
        "subtract",
        ["wide", "wall"],
        Some("1"),
    ));
    let base = built.model.blend_base(built.result).unwrap();
    let accuracy = Accuracy {
        relative: q("1/1000000"),
        absolute: q("1/1000"),
    };
    let measure = |body| {
        let value = volume(&built.model, body, &accuracy, &mut decider())
            .unwrap()
            .value;
        ((&value.lo + &value.hi) / q("2")).to_f64_nearest()
    };
    let change = measure(built.result) - measure(base);
    let want = -2.0 * square(1.0).section_integral(-1.0);
    assert!((change - want).abs() < -want * 1e-3, "{change} {want}");
}

// The base reads a step as rings alone, so anything else in the step takes it away.
#[test]
fn a_ring_with_another_strip_has_no_base() {
    let top = json!({"between":[{"step":"wall","role":{"side":{"segment":"wall"}}},
        {"step":"wall","role":"end"}]});
    let mut edges = ring();
    edges.as_array_mut().unwrap().push(top);
    let built = build(&parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":edges}),
        ],
        "round",
    ));
    assert_eq!(built.fillets.len(), 2);
    assert_eq!(built.model.blend_base(built.result), None);
    let accuracy = Accuracy {
        relative: q("1/1000"),
        absolute: q("1"),
    };
    let refused = volume(&built.model, built.result, &accuracy, &mut decider());
    assert_eq!(refused.unwrap_err().code, "measure.unsupported-surface");
}

#[test]
fn a_pattern_of_rounded_bodies_builds_with_a_box() {
    let copies = build(&parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
            json!({"id":"copies","op":"pattern","body":"round",
                "layout":{"linear":{"vector":["100","0","0"],"count":2}}}),
        ],
        "copies",
    ));
    assert_eq!(bands(&copies).len(), 2);
    let found = bounding_box(
        &copies.model,
        copies.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value
    .unwrap();
    // The first copy's wall sets the low end of x, and the second copy's the high.
    assert!((&found.min[0] + &q("10")).abs() < q("1/1000000"));
    assert!((&found.max[0] - &q("110")).abs() < q("1/1000000"));
}

#[test]
fn face_measures_and_the_field_stay_refused() {
    let built = build(&boss("1"));
    let accuracy = Accuracy {
        relative: q("1/1000"),
        absolute: q("1"),
    };
    let mass = mass_properties(&built.model, built.result, &accuracy, &mut decider());
    assert_eq!(mass.unwrap_err().code, "measure.unsupported-surface");
    let [band] = bands(&built)[..] else {
        panic!("one band");
    };
    let area = face_area(&built.model, band, &accuracy, &mut decider());
    assert_eq!(area.unwrap_err().code, "measure.unsupported-surface");
    // The walls its rails trim, and a rail itself.
    for face in built.model.faces_of_body(built.result) {
        if matches!(
            built.model.geom().surface(built.model.face(face).surface),
            Surface::Cylinder(_)
        ) {
            let area = face_area(&built.model, face, &accuracy, &mut decider());
            assert_eq!(area.unwrap_err().code, "measure.unsupported-surface");
            let found = face_bounding_box(&built.model, face, &accuracy, &mut decider());
            assert_eq!(found.unwrap_err().code, "measure.unsupported-surface");
        }
    }
    let rail = catallax_topology::saddle_rows(&built.model, band).unwrap()[0].0;
    let length = edge_length(&built.model, rail, &accuracy, &mut decider());
    assert_eq!(length.unwrap_err().code, "measure.unsupported-surface");
    let field = Plan::compile_checked(&boss("1"), &mut decider()).unwrap_err();
    assert_eq!(field.code, "field.fillet-contact");
}

// A file's whole face, edge and vertex count, and its bound against the body.
fn written(built: &Built) -> catallax_step::StepReport {
    let mut bytes = Vec::new();
    let answer = write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut decider(),
        &mut bytes,
    )
    .unwrap();
    let file = String::from_utf8(bytes).unwrap();
    assert!(file.contains("SEAM_CURVE"));
    assert!(answer.value.max_deviation <= q("1/10000000"));
    answer.value
}

#[test]
fn a_ring_certifies_and_writes_each_cylinder_chart() {
    let built = build(&boss("1"));
    let report = written(&built);
    // The band and each wall add a seam, and each closed edge gains a vertex.
    assert_eq!(
        (report.faces, report.edges, report.vertices),
        (6, 8, 5),
        "{report:?}"
    );
    assert_eq!((report.added_edges, report.added_vertices), (3, 5));
    let certified = certify_junction(built, &mut decider()).unwrap();
    let junction = certified.junction().unwrap();
    let [saddle] = junction.saddles() else {
        panic!("one fitted ring");
    };
    assert_eq!(junction.faces().len(), 6);
    assert!(junction.envelopes().is_empty());
    let fitted = junction
        .record()
        .regions
        .values()
        .filter(|region| region.exception.is_some())
        .count();
    assert_eq!(fitted, 1);
    // The band and both walls it trims carry the fit, and the end caps do not.
    let trimmed = junction
        .faces()
        .iter()
        .filter(|face| !junction.face_bound(**face).is_zero())
        .count();
    assert_eq!(trimmed, 3);
    assert_eq!(junction.face_bound(saddle.face()), saddle.bound().upper);
    assert!(saddle.bound().upper <= *junction.fit_target_mm());
}

#[test]
fn a_through_bore_a_mirrored_copy_and_a_lean_certify() {
    let hole = cylinder("wall", "4", "60", LEVEL, ["1", "0", "-30"]);
    let through = build(&recipe(
        [wide(), hole],
        "subtract",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(written(&through).faces, 6);
    let certified = certify_junction(through, &mut decider()).unwrap();
    assert_eq!(certified.junction().unwrap().saddles().len(), 2);

    let flipped = build(&parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
            json!({"id":"flipped","op":"mirror","body":"round",
                "plane":{"point":["20","0","0"],"normal":["3","4","0"]},"keep_original":false}),
        ],
        "flipped",
    ));
    assert_eq!(written(&flipped).edges, 8);
    let certified = certify_junction(flipped, &mut decider()).unwrap();
    assert_eq!(certified.junction().unwrap().saddles().len(), 1);

    let lean = cylinder("wall", "4", "30", LEAN, ["1", "0", "0"]);
    let leaning = build(&recipe(
        [wide(), lean],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(written(&leaning).edges, 8);
    certify_junction(leaning, &mut decider()).unwrap();
}

/// A narrow cylinder that runs past the side of the wide one, so the two cross in
/// one closed edge that winds round neither.
fn past(x: &str) -> Value {
    cylinder("wall", "4", "60", LEVEL, [x, "0", "-30"])
}
fn over(kind: &str, order: [&str; 2], radius: &str) -> Recipe {
    recipe([wide(), past("8")], kind, order, Some(radius))
}

/// A ball on those two cylinders, from its two distance conditions and sharing
/// no code with the kernel. The wide axis is the y axis and the narrow one runs
/// along z through `(8, 0)`. A side of `1` puts the ball outside that cylinder.
#[derive(Clone, Copy)]
struct Ball {
    radius: f64,
    sides: [f64; 2],
}
impl Ball {
    fn offsets(&self) -> [f64; 2] {
        [
            10.0 + self.sides[0] * self.radius,
            4.0 + self.sides[1] * self.radius,
        ]
    }
    fn distances(p: [f64; 3]) -> [f64; 2] {
        [p[0].hypot(p[2]), (p[0] - 8.0).hypot(p[1])]
    }
    // The centre at angle `t` round the narrow axis, above or below the plane
    // of the wide axis and the common perpendicular.
    fn centre(&self, t: f64, up: f64) -> Option<[f64; 3]> {
        let [big, small] = self.offsets();
        let (x, y) = (8.0 + small * t.cos(), small * t.sin());
        (x.abs() < big).then(|| [x, y, up * (big * big - x * x).sqrt()])
    }
    // A unit ray from a centre, halfway between its two contacts.
    fn ray(&self, c: [f64; 3]) -> [f64; 3] {
        let [big, small] = self.offsets();
        let on_wide = [c[0] * 10.0 / big, c[1], c[2] * 10.0 / big];
        let on_narrow = [8.0 + (c[0] - 8.0) * 4.0 / small, c[1] * 4.0 / small, c[2]];
        let ray: [f64; 3] = std::array::from_fn(|i| on_wide[i] + on_narrow[i] - 2.0 * c[i]);
        ray.map(|value| value / dot(ray, ray).sqrt())
    }
    fn allowed(&self, p: [f64; 3]) -> bool {
        let [big, small] = self.offsets();
        let [to_wide, to_narrow] = Self::distances(p);
        self.sides[0] * (to_wide - big) >= 0.0 && self.sides[1] * (to_narrow - small) >= 0.0
    }
    // Whether some ball that fits covers a point: one on either wall, or one on
    // the edge between them. `None` is too close to call.
    fn covers(&self, p: [f64; 3]) -> Option<bool> {
        let [big, small] = self.offsets();
        let [to_wide, to_narrow] = Self::distances(p);
        let on_wide = [p[0] * big / to_wide, p[1], p[2] * big / to_wide];
        let on_narrow = [
            8.0 + (p[0] - 8.0) * small / to_narrow,
            p[1] * small / to_narrow,
            p[2],
        ];
        let mut nearest = f64::INFINITY;
        if self.allowed(on_wide) {
            nearest = nearest.min((to_wide - big).abs());
        }
        if self.allowed(on_narrow) {
            nearest = nearest.min((to_narrow - small).abs());
        }
        for i in 0..2880 {
            let t = 2.0 * PI * f64::from(i) / 2880.0;
            for up in [1.0, -1.0] {
                if let Some(c) = self.centre(t, up) {
                    let d: [f64; 3] = std::array::from_fn(|k| p[k] - c[k]);
                    nearest = nearest.min(dot(d, d).sqrt());
                }
            }
        }
        ((nearest - self.radius).abs() > 1e-2).then_some(nearest < self.radius)
    }
    // Points round the loop on each side of the band, each with whether it stays
    // in the solid. Past the band a point is in the corner the ball cannot reach.
    fn witnesses(&self) -> Vec<([f64; 3], bool)> {
        let solid_at_ball = self.sides.iter().any(|side| *side < 0.0);
        let [outside_wide, outside_narrow] = self.sides.map(|side| side > 0.0);
        let mut found = Vec::new();
        for i in 0..8 {
            let t = 2.0 * PI * (f64::from(i) + 0.5) / 8.0;
            for up in [1.0, -1.0] {
                let Some(c) = self.centre(t, up) else {
                    continue;
                };
                let ray = self.ray(c);
                for share in [0.9, 1.04] {
                    let p: [f64; 3] = std::array::from_fn(|k| c[k] + share * self.radius * ray[k]);
                    let [to_wide, to_narrow] = Self::distances(p);
                    // Only the corner's own cell says anything about the ball.
                    if (to_wide > 10.0) != outside_wide || (to_narrow > 4.0) != outside_narrow {
                        continue;
                    }
                    if let Some(covered) = self.covers(p) {
                        found.push((p, covered == solid_at_ball));
                    }
                }
            }
        }
        found
    }
}

// Each witness on `built`, whose x axis `flip` may have turned over.
fn assert_witnesses(built: &mut Built, ball: Ball, flip: f64) {
    let witnesses = ball.witnesses();
    let solid = witnesses.iter().filter(|(_, inside)| *inside).count();
    assert!(
        solid >= 4 && witnesses.len() - solid >= 4,
        "{solid} of {}",
        witnesses.len()
    );
    for ([x, y, z], inside) in witnesses {
        let text = [flip * x, y, z].map(|value| format!("{}/4096", (value * 4096.0).round()));
        let want = if inside {
            Containment::Inside
        } else {
            Containment::Outside
        };
        classify(built, &[text.each_ref().map(String::as_str)], want);
    }
}

// Rounding a closed solid moves its surface by a little, so a trim that kept the
// wrong side of a loop would change the size by far more than this.
fn assert_close_to(built: &Built, plain: &Recipe) {
    let size = |mesh: &Mesh| -> f64 {
        mesh.triangles
            .iter()
            .map(|triangle| {
                let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
                (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]))
                    / 6.0
            })
            .sum()
    };
    let mesh = fine(built);
    assert_closed(&mesh);
    let [rounded, sharp] = [size(&mesh), size(&fine(&build(plain)))];
    assert!(
        (rounded - sharp).abs() < sharp * 2e-2,
        "rounded {rounded}, sharp {sharp}"
    );
}

#[test]
fn a_loop_rounds_each_way_two_cylinders_can_meet() {
    for (kind, order, radius, sides, faces) in [
        ("union", ["wide", "wall"], "1", [1.0, 1.0], 7),
        ("subtract", ["wide", "wall"], "1", [-1.0, 1.0], 5),
        ("intersect", ["wide", "wall"], "1", [-1.0, -1.0], 3),
        ("subtract", ["wall", "wide"], "1/2", [1.0, -1.0], 5),
    ] {
        let mut built = build(&over(kind, order, radius));
        assert_eq!(built.model.faces_of_body(built.result).len(), faces);
        assert_eq!(bands(&built).len(), 1, "{kind}");
        assert_close_to(&built, &recipe([wide(), past("8")], kind, order, None));
        let ball = Ball {
            radius: q(radius).to_f64_nearest(),
            sides,
        };
        assert_witnesses(&mut built, ball, 1.0);
    }
}

#[test]
fn a_larger_ball_and_either_face_order_round_the_same_loop() {
    let ball = Ball {
        radius: 3.0,
        sides: [-1.0, 1.0],
    };
    assert_witnesses(
        &mut build(&over("subtract", ["wide", "wall"], "3")),
        ball,
        1.0,
    );
    let swapped = recipe([past("8"), wide()], "union", ["wall", "wide"], Some("1"));
    let ball = Ball {
        radius: 1.0,
        sides: [1.0, 1.0],
    };
    assert_witnesses(&mut build(&swapped), ball, 1.0);
}

// A loop on the other side, a lean, and one almost past the wide cylinder.
#[test]
fn a_loop_meshes_closed_wherever_it_sits() {
    let lean = cylinder("wall", "4", "60", LEAN, ["8", "-18", "-24"]);
    let equal = cylinder("wall", "10", "60", LEVEL, ["1", "0", "-30"]);
    for (narrow, kind) in [
        (past("-8"), "union"),
        (lean.clone(), "union"),
        (lean, "subtract"),
        (equal, "union"),
        (past("27/2"), "union"),
    ] {
        let parts = [wide(), narrow];
        let built = build(&recipe(parts.clone(), kind, ["wide", "wall"], Some("1")));
        assert_eq!(bands(&built).len(), 1);
        assert_close_to(&built, &recipe(parts, kind, ["wide", "wall"], None));
    }
}

#[test]
fn a_mirrored_loop_keeps_its_band() {
    let mut steps = vec![
        wide(),
        past("8"),
        json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
        json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
    ];
    steps.push(json!({"id":"flipped","op":"mirror","body":"round",
        "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
    let mut flipped = build(&parse(steps, "flipped"));
    assert_eq!(bands(&flipped).len(), 1);
    assert_closed(&fine(&flipped));
    let ball = Ball {
        radius: 1.0,
        sides: [1.0, 1.0],
    };
    assert_witnesses(&mut flipped, ball, -1.0);
}

#[test]
fn a_loop_past_the_end_of_a_cylinder_refuses() {
    // The narrow rail climbs to about 10.58 over the wide cylinder.
    let short = cylinder("wall", "4", "21", LEVEL, ["8", "0", "-21/2"]);
    let (code, _) = refusal(&recipe(
        [wide(), short],
        "union",
        ["wide", "wall"],
        Some("1"),
    ));
    assert_eq!(code, "fillet.unsupported-overflow");
}

// At this radius the offsets' shadows touch, though the cylinders' own cross.
#[test]
fn a_ball_that_makes_a_loops_offsets_tangent_refuses() {
    let (code, _) = refusal(&over("subtract", ["wall", "wide"], "1"));
    assert_eq!(code, "fillet.unsupported-contact");
}

#[test]
fn a_loops_box_comes_from_its_walls_and_band() {
    let accuracy = Accuracy {
        relative: q("1/1000"),
        absolute: q("1/1000"),
    };
    let found = |recipe: &Recipe| {
        let built = build(recipe);
        let bounds = bounding_box(&built.model, built.result, &accuracy, &mut decider())
            .unwrap()
            .value
            .unwrap();
        let read = |side: &[Rational; 3]| -> [f64; 3] {
            std::array::from_fn(|i| side[i].to_f64_nearest())
        };
        [read(&bounds.min), read(&bounds.max)]
    };
    let near = |a: f64, c: f64| assert!((a - c).abs() < 1e-2, "{a} {c}");
    let [low, high] = found(&over("union", ["wide", "wall"], "1"));
    near(low[0], -10.0);
    near(high[0], 12.0);
    // A lens's walls are bounded by the rails alone. Generators of each cylinder
    // stay on it, so three sides do not move, and rounding lowers its top.
    let [low, high] = found(&over("intersect", ["wide", "wall"], "1"));
    near(low[0], 4.0);
    near(high[1], 4.0);
    near(low[1], -4.0);
    near(high[0], 10.0);
    assert!(high[2] < 84f64.sqrt() && low[2] > -(84f64.sqrt()));
}

// Each of these reads a ring that winds round the narrow cylinder.
#[test]
fn a_loops_measures_file_and_certificate_stay_refused() {
    let built = build(&over("union", ["wide", "wall"], "1"));
    let accuracy = Accuracy {
        relative: q("1/1000"),
        absolute: q("1"),
    };
    let found = volume(&built.model, built.result, &accuracy, &mut decider());
    assert_eq!(found.unwrap_err().code, "measure.unsupported-surface");
    let found = mass_properties(&built.model, built.result, &accuracy, &mut decider());
    assert_eq!(found.unwrap_err().code, "measure.unsupported-surface");
    let [band] = bands(&built)[..] else {
        panic!("one band");
    };
    let found = face_area(&built.model, band, &accuracy, &mut decider());
    assert_eq!(found.unwrap_err().code, "measure.unsupported-surface");
    let rail = catallax_topology::saddle_rows(&built.model, band).unwrap()[0].0;
    let found = edge_length(&built.model, rail, &accuracy, &mut decider());
    assert_eq!(found.unwrap_err().code, "measure.unsupported-surface");
    let mut bytes = Vec::new();
    let file = write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut decider(),
        &mut bytes,
    );
    assert_eq!(
        file.map(|_| ()).unwrap_err().code,
        "step.unsupported-geometry"
    );
    let certified = certify_junction(built, &mut decider());
    assert_eq!(
        certified.map(|_| ()).unwrap_err().code,
        "numbers.geometric-gap-unproven"
    );
    let field = Plan::compile_checked(&over("union", ["wide", "wall"], "1"), &mut decider());
    assert_eq!(field.map(|_| ()).unwrap_err().code, "field.fillet-contact");
}
