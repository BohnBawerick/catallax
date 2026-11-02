//! Rolling-ball fillets on the ring where two nonparallel cylinders cross.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::Plan;
use catallax_measure::{Accuracy, volume};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, with_operation_log};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{EnvelopeGrid, RPoint3, SplineSupport, Surface};
use catallax_step::write_step;
use catallax_topology::FaceId;
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
        const ANGLES: usize = 720;
        const RADII: usize = 480;
        let width = 2.0 * self.ball / RADII as f64;
        let turn = 2.0 * PI / ANGLES as f64;
        let mut sum = 0.0;
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
                sum += height * radial * width * turn;
            }
        }
        sum
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
fn a_leaning_boss_has_no_plane_of_symmetry_to_lean_on() {
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
fn a_cached_ring_records_its_phases() {
    let joined = recipe([wide(), wall()], "union", ["wide", "wall"], None);
    let prefix = evaluate_incremental(&joined, None, &mut decider()).unwrap();
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&boss("1"), Some(prefix), &mut dec));
    let built = result.unwrap().into_built();
    for name in ["fillet.contact", "fillet.fitting", "fillet.trimming"] {
        let measured: Vec<_> = phases.iter().filter(|p| p.operation == name).collect();
        assert!(!measured.is_empty(), "{name} was not recorded");
        let total: f64 = measured.iter().map(|p| p.wall_ms).sum();
        eprintln!("{name}={total:.3}ms over {} spans", measured.len());
    }
    let start = std::time::Instant::now();
    assert_eq!(built.model.check(built.result), vec![]);
    eprintln!("certification={:?}", start.elapsed());
    let start = std::time::Instant::now();
    fine(&built);
    eprintln!("display={:?}", start.elapsed());
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
    // A copy that keeps its original needs a body box, which waits for measurement.
    let copies = parse(
        vec![
            wide(),
            wall(),
            json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
            json!({"id":"round","op":"fillet","body":"joined","radius":"1","edges":ring()}),
            json!({"id":"copies","op":"pattern","body":"round",
                "layout":{"linear":{"vector":["100","0","0"],"count":2}}}),
        ],
        "copies",
    );
    assert_eq!(refusal(&copies).0, "measure.unsupported-surface");
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

// A partial overlap has one ring, and no section turns all the way round.
#[test]
fn a_partial_overlap_refuses() {
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

#[test]
fn volume_step_and_the_field_stay_refused() {
    let built = build(&boss("1"));
    let refused = volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: q("1/1000"),
            absolute: q("1"),
        },
        &mut decider(),
    )
    .unwrap_err();
    assert_eq!(refused.code, "measure.unsupported-surface");
    let step = write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut decider(),
        &mut Vec::new(),
    )
    .unwrap_err();
    assert_eq!(step.code, "step.unsupported-geometry");
    let field = Plan::compile_checked(&boss("1"), &mut decider()).unwrap_err();
    assert_eq!(field.code, "field.fillet-contact");
}
