//! 4.2a step 3: quadric pairs. The coaxial file's rounds, the saddle file's
//! rings and loops, the sphere-ring file and the equal-crossing file.
//!
//! The field compiles these fillets from its own path pieces. Its membership
//! must match the witnesses of each test file (checked against the B-rep, and
//! for the corpus rows against OpenCASCADE). Where the B-rep builds and
//! displays, every display vertex of a B-rep blend face must lie within the
//! field's guaranteed distance range of the field's own face, at the 1 nm
//! own-part tolerance, with the same records at 1, 2 and 4 workers.
use catallax_features::evaluate;
use catallax_field::rolling_region::piece::Piece;
use catallax_field::rolling_region::{Batch, Range, distance_batch};
use catallax_field::{Location, Plan, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, OperationAllowance, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RPoint3, SplineSupport, Surface};
use serde_json::{Value, json};
use std::collections::BTreeMap;

const UP: [&str; 4] = ["1", "0", "0", "0"];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn parse(value: Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
/// A binary64 coordinate as an exact rational, to nine decimals as the B-rep
/// test files write their witnesses.
fn exact(x: f64) -> Rational {
    q(&format!("{x:.9}"))
}

/// One recipe with points whose membership a test file fixed.
struct Case {
    name: String,
    recipe: Recipe,
    inside: Vec<[Rational; 3]>,
    outside: Vec<[Rational; 3]>,
    /// The B-rep builds and displays it and no earlier case has the same
    /// recipe, so its display vertices are checked.
    brep: bool,
    /// Rolling pieces the field compiles.
    pieces: usize,
}

// The coaxial file: rounds on circles where two surfaces of revolution about
// one axis meet. Each ball is worked by hand in the half-plane through the
// axis, as in that file.
type Pair = [f64; 2];
fn unit(v: Pair) -> Pair {
    let length = (v[0] * v[0] + v[1] * v[1]).sqrt();
    [v[0] / length, v[1] / length]
}
/// Where two lines `n . p = d` meet.
fn meet(a: Pair, d: f64, b: Pair, e: f64) -> Pair {
    let det = a[0] * b[1] - a[1] * b[0];
    [(d * b[1] - e * a[1]) / det, (a[0] * e - b[0] * d) / det]
}
/// A point inside the ball, and one in the corner between it and the edge.
fn ball_witnesses(edge: Pair, centre: Pair, radius: f64) -> [Pair; 2] {
    let (dx, dz) = (edge[0] - centre[0], edge[1] - centre[1]);
    let length = (dx * dx + dz * dz).sqrt();
    let near = [
        centre[0] + dx * radius / length,
        centre[1] + dz * radius / length,
    ];
    let half = radius / 2.0 / length;
    let inside = [centre[0] + dx * half, centre[1] + dz * half];
    [
        inside,
        [(near[0] + edge[0]) / 2.0, (near[1] + edge[1]) / 2.0],
    ]
}
fn pose(z: &str) -> Value {
    json!({"rotation": UP, "translation": ["0", "0", z]})
}
fn face(step: &str, role: &str) -> Value {
    json!({"step": step, "role": role})
}
fn side(step: &str, segment: &str) -> Value {
    json!({"step": step, "role": {"side": {"segment": segment}}})
}
fn wall(step: &str) -> Value {
    side(step, "wall")
}
fn block(id: &str, z: &str, height: &str) -> Vec<Value> {
    vec![
        json!({"id": format!("{id}-section"), "op": "sketch",
        "plane": {"rotation": UP, "translation": ["0", "0", z]},
        "loops": [[
            {"kind": "line", "id": "s", "to": ["20", "-20"]},
            {"kind": "line", "id": "e", "to": ["20", "20"]},
            {"kind": "line", "id": "n", "to": ["-20", "20"]},
            {"kind": "line", "id": "w", "to": ["-20", "-20"]}
        ]]}),
        json!({"id": id, "op": "extrude", "profile": format!("{id}-section"), "distance": height}),
    ]
}
fn cone(base: &str, top: &str, height: &str) -> Vec<Value> {
    vec![
        json!({"id": "c", "op": "cone", "base_radius": base, "top_radius": top,
        "height": height, "pose": pose("0")}),
    ]
}
fn chamfered() -> Vec<Value> {
    vec![
        json!({"id": "c", "op": "cylinder", "radius": "5", "height": "10", "pose": pose("0")}),
        json!({"id": "bevel", "op": "chamfer", "body": "c", "distance": "1",
            "edges": [{"between": [face("c", "end"), wall("c")]}]}),
    ]
}
fn dome() -> Vec<Value> {
    let mut steps = vec![json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("0")})];
    steps.extend(block("b", "2", "10"));
    steps
        .push(json!({"id": "j", "op": "boolean", "kind": "intersect", "target": "s", "tool": "b"}));
    steps
}
fn seat() -> Vec<Value> {
    let mut steps = block("b", "0", "20");
    steps.push(json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("23")}));
    steps.push(json!({"id": "j", "op": "boolean", "kind": "union", "target": "b", "tool": "s"}));
    steps
}
fn countersunk() -> Vec<Value> {
    let mut steps = block("b", "0", "20");
    steps.push(
        json!({"id": "h", "op": "hole", "target": "b", "at": [{"id": "p",
        "point": ["0", "0", "20"], "direction": ["0", "0", "-1"]}],
        "diameter": "6", "depth": "through", "tip": "flat",
        "head": {"countersink": {"diameter": "12", "angle": "90"}}}),
    );
    steps
}
fn hole(part: &str) -> Value {
    json!({"step": "h", "role": {"hole": {"at": "p", "part": part}}})
}
fn capped_rod() -> Vec<Value> {
    vec![
        json!({"id": "section", "op": "sketch",
            "plane": {"rotation": UP, "translation": ["0", "0", "0"]}, "loops": [[
            {"kind": "line", "id": "base", "to": ["6", "0"]},
            {"kind": "line", "id": "rod", "to": ["6", "8"]},
            {"kind": "arc", "id": "cap", "centre": ["0", "0"], "to": ["0", "10"], "turn": "ccw"},
            {"kind": "line", "id": "axis", "to": ["0", "0"]}
        ]]}),
        json!({"id": "r", "op": "revolve", "profile": "section",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "1", "0"]}, "angle": "360"}),
    ]
}
fn rounded(steps: &[Value], body: &str, radius: &str, edges: &[[Value; 2]]) -> Recipe {
    let mut steps = steps.to_vec();
    let edges: Vec<_> = edges.iter().map(|e| json!({"between": e})).collect();
    steps.push(
        json!({"id": "round", "op": "fillet", "body": body, "radius": radius, "edges": edges}),
    );
    parse(
        json!({"format": "catallax-recipe/3", "name": "coaxial round", "steps": steps,
        "result": "round"}),
    )
}

/// One coaxial file case: its balls as `(edge, centre)` in the half-plane, and
/// whether rounding removes material.
#[allow(clippy::too_many_arguments)]
fn coaxial(
    name: &str,
    steps: Vec<Value>,
    body: &str,
    radius: &str,
    edges: &[[Value; 2]],
    balls: &[(Pair, Pair)],
    convex: bool,
    about_y: bool,
) -> Case {
    let r: f64 = Rational::parse(radius).unwrap().to_f64_nearest();
    let at = |p: Pair| -> [Rational; 3] {
        if about_y {
            [exact(p[0]), exact(p[1]), Rational::zero()]
        } else {
            [exact(p[0]), Rational::zero(), exact(p[1])]
        }
    };
    let mut inside = Vec::new();
    let mut outside = Vec::new();
    for (edge, centre) in balls {
        let [ball, corner] = ball_witnesses(*edge, *centre, r);
        let (solid, void) = if convex {
            (ball, corner)
        } else {
            (corner, ball)
        };
        inside.push(at(solid));
        outside.push(at(void));
    }
    Case {
        name: name.into(),
        recipe: rounded(&steps, body, radius, edges),
        inside,
        outside,
        brep: true,
        pieces: balls.len(),
    }
}

fn coaxial_cases() -> Vec<Case> {
    let mut cases = Vec::new();
    for (name, base, top, n, d, edge_role, edge, level) in [
        (
            "frustum top",
            "5",
            "3",
            unit([3.0, 1.0]),
            15.0 / 10f64.sqrt(),
            "end",
            [3.0, 6.0],
            5.0,
        ),
        (
            "frustum base",
            "5",
            "3",
            unit([3.0, 1.0]),
            15.0 / 10f64.sqrt(),
            "start",
            [5.0, 0.0],
            1.0,
        ),
        (
            "flared top",
            "3",
            "5",
            unit([3.0, -1.0]),
            9.0 / 10f64.sqrt(),
            "end",
            [5.0, 6.0],
            5.0,
        ),
        (
            "narrow top",
            "4",
            "1",
            unit([2.0, 1.0]),
            8.0 / 5f64.sqrt(),
            "end",
            [1.0, 6.0],
            5.0,
        ),
    ] {
        cases.push(coaxial(
            name,
            cone(base, top, "6"),
            "c",
            "1",
            &[[face("c", edge_role), wall("c")]],
            &[(edge, meet(n, d - 1.0, [0.0, 1.0], level))],
            true,
            false,
        ));
    }
    let n = unit([1.0, 1.0]);
    let d = 14.0 / 2f64.sqrt();
    let chamfer = |r: f64| {
        [
            ([4.0, 10.0], meet(n, d - r, [0.0, 1.0], 10.0 - r)),
            ([5.0, 9.0], meet(n, d - r, [1.0, 0.0], 5.0 - r)),
        ]
    };
    let [top, flank] = [
        [face("c", "end"), json!({"step": "bevel"})],
        [wall("c"), json!({"step": "bevel"})],
    ];
    cases.push(coaxial(
        "chamfer top",
        chamfered(),
        "bevel",
        "1/2",
        std::slice::from_ref(&top),
        &chamfer(0.5)[..1],
        true,
        false,
    ));
    cases.push(coaxial(
        "chamfer side",
        chamfered(),
        "bevel",
        "1/2",
        std::slice::from_ref(&flank),
        &chamfer(0.5)[1..],
        true,
        false,
    ));
    cases.push(coaxial(
        "chamfer both",
        chamfered(),
        "bevel",
        "1/4",
        &[top, flank],
        &chamfer(0.25),
        true,
        false,
    ));
    cases.push(coaxial(
        "dome rim",
        dome(),
        "j",
        "1",
        &[[face("b", "start"), wall("s")]],
        &[([21f64.sqrt(), 2.0], [(16.0f64 - 9.0).sqrt(), 3.0])],
        true,
        false,
    ));
    cases.push(coaxial(
        "sphere seat",
        seat(),
        "j",
        "1",
        &[[face("b", "end"), wall("s")]],
        &[([4.0, 20.0], [32f64.sqrt(), 21.0])],
        false,
        false,
    ));
    let n = unit([1.0, -1.0]);
    let d = -14.0 / 2f64.sqrt();
    // The mouth's hole entry cylinder passes through the mouth circle into
    // the region, which only its face's extent keeps clear; the mouth is
    // among the refusals below.
    cases.push(coaxial(
        "countersink throat",
        countersunk(),
        "h",
        "1",
        &[[hole("countersink"), hole("wall")]],
        &[([3.0, 17.0], meet(n, d + 1.0, [1.0, 0.0], 4.0))],
        true,
        false,
    ));
    cases.push(coaxial(
        "rod cap",
        capped_rod(),
        "r",
        "1",
        &[[side("r", "rod"), side("r", "cap")]],
        &[([6.0, 8.0], [5.0, 56f64.sqrt()])],
        true,
        true,
    ));
    cases
}

// The saddle file: a wide cylinder along y and a narrow one that crosses it.
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
fn narrow() -> Value {
    cylinder("wall", "4", "30", LEVEL, ["1", "0", "0"])
}
/// A narrow cylinder through the whole height, past the wide one's side.
fn past(x: &str) -> Value {
    cylinder("wall", "4", "60", LEVEL, [x, "0", "-30"])
}
fn ring() -> Value {
    json!([{"between":[{"step":"wide","role":{"side":{"segment":"wall"}}},
        {"step":"wall","role":{"side":{"segment":"wall"}}}]}])
}
fn saddle_steps(parts: [Value; 2], kind: &str, order: [&str; 2], radius: &str) -> Vec<Value> {
    let [first, second] = parts;
    vec![
        first,
        second,
        json!({"id":"joined","op":"boolean","kind":kind,"target":order[0],"tool":order[1]}),
        json!({"id":"round","op":"fillet","body":"joined","radius":radius,"edges":ring()}),
    ]
}
fn saddle_recipe(steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    parse(json!({"format":"catallax-recipe/3","name":"ring","steps":steps,"result":result}))
}

/// A ball on the wide cylinder and a narrow one along z through `(8, 0)`, from
/// the two distance conditions alone, as the saddle file models it. A side of
/// `1` puts the ball outside that cylinder.
#[derive(Clone, Copy)]
struct LoopBall {
    radius: f64,
    sides: [f64; 2],
}
fn dot3(a: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * c[0] + a[1] * c[1] + a[2] * c[2]
}
impl LoopBall {
    fn offsets(&self) -> [f64; 2] {
        [
            10.0 + self.sides[0] * self.radius,
            4.0 + self.sides[1] * self.radius,
        ]
    }
    fn distances(p: [f64; 3]) -> [f64; 2] {
        [p[0].hypot(p[2]), (p[0] - 8.0).hypot(p[1])]
    }
    fn centre(&self, t: f64, up: f64) -> Option<[f64; 3]> {
        let [big, small] = self.offsets();
        let (x, y) = (8.0 + small * t.cos(), small * t.sin());
        (x.abs() < big).then(|| [x, y, up * (big * big - x * x).sqrt()])
    }
    fn ray(&self, c: [f64; 3]) -> [f64; 3] {
        let [big, small] = self.offsets();
        let on_wide = [c[0] * 10.0 / big, c[1], c[2] * 10.0 / big];
        let on_narrow = [8.0 + (c[0] - 8.0) * 4.0 / small, c[1] * 4.0 / small, c[2]];
        let ray: [f64; 3] = std::array::from_fn(|i| on_wide[i] + on_narrow[i] - 2.0 * c[i]);
        ray.map(|value| value / dot3(ray, ray).sqrt())
    }
    fn allowed(&self, p: [f64; 3]) -> bool {
        let [big, small] = self.offsets();
        let [to_wide, to_narrow] = Self::distances(p);
        self.sides[0] * (to_wide - big) >= -1e-9 && self.sides[1] * (to_narrow - small) >= -1e-9
    }
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
            let t = std::f64::consts::TAU * f64::from(i) / 2880.0;
            for up in [1.0, -1.0] {
                if let Some(c) = self.centre(t, up) {
                    let d: [f64; 3] = std::array::from_fn(|k| p[k] - c[k]);
                    nearest = nearest.min(dot3(d, d).sqrt());
                }
            }
        }
        ((nearest - self.radius).abs() > 1e-2).then_some(nearest < self.radius)
    }
    /// Points round the loop on each side of the band, each with whether it
    /// stays in the solid, in `flip`'s mirror of x.
    fn witnesses(&self, flip: f64) -> (Vec<[Rational; 3]>, Vec<[Rational; 3]>) {
        let solid_at_ball = self.sides.iter().any(|side| *side < 0.0);
        let [outside_wide, outside_narrow] = self.sides.map(|side| side > 0.0);
        let (mut inside, mut outside) = (Vec::new(), Vec::new());
        for i in 0..8 {
            let t = std::f64::consts::TAU * (f64::from(i) + 0.5) / 8.0;
            for up in [1.0, -1.0] {
                let Some(c) = self.centre(t, up) else {
                    continue;
                };
                let ray = self.ray(c);
                for share in [0.9, 1.04] {
                    let p: [f64; 3] = std::array::from_fn(|k| c[k] + share * self.radius * ray[k]);
                    let [to_wide, to_narrow] = Self::distances(p);
                    if (to_wide > 10.0) != outside_wide || (to_narrow > 4.0) != outside_narrow {
                        continue;
                    }
                    if let Some(covered) = self.covers(p) {
                        let point = [flip * p[0], p[1], p[2]]
                            .map(|v| q(&format!("{}/4096", (v * 4096.0).round())));
                        if covered == solid_at_ball {
                            inside.push(point);
                        } else {
                            outside.push(point);
                        }
                    }
                }
            }
        }
        (inside, outside)
    }
}

fn saddle_case(name: &str, steps: Vec<Value>, inside: &[[&str; 3]], outside: &[[&str; 3]]) -> Case {
    Case {
        name: name.into(),
        recipe: saddle_recipe(steps),
        inside: inside.iter().map(|p| p.map(q)).collect(),
        outside: outside.iter().map(|p| p.map(q)).collect(),
        brep: true,
        pieces: 1,
    }
}

fn saddle_cases() -> Vec<Case> {
    let mut cases = vec![
        saddle_case(
            "saddle boss",
            saddle_steps([wide(), narrow()], "union", ["wide", "wall"], "1"),
            &[
                ["1", "81/20", "10"],
                ["1", "-81/20", "10"],
                ["101/20", "0", "87/10"],
                ["-61/20", "0", "48/5"],
                ["1", "0", "20"],
                ["0", "20", "0"],
            ],
            &[
                ["1", "5", "219/20"],
                ["1", "-5", "219/20"],
                ["6", "0", "46/5"],
                ["-4", "0", "41/4"],
                ["20", "0", "0"],
            ],
        ),
        saddle_case(
            "saddle bore",
            saddle_steps([wide(), narrow()], "subtract", ["wide", "wall"], "1"),
            &[["1", "5", "89/10"], ["6", "0", "71/10"], ["-6", "0", "0"]],
            &[
                ["1", "81/20", "99/10"],
                ["1", "-81/20", "99/10"],
                ["101/20", "0", "43/5"],
            ],
        ),
        saddle_case(
            "saddle boss r=2, wall first",
            saddle_steps([wide(), narrow()], "union", ["wall", "wide"], "2"),
            &[],
            &[],
        ),
        saddle_case(
            "saddle boss r=3/2",
            saddle_steps([wide(), narrow()], "union", ["wide", "wall"], "3/2"),
            &[],
            &[],
        ),
        saddle_case(
            "leaning boss",
            saddle_steps(
                [wide(), cylinder("wall", "4", "30", LEAN, ["1", "0", "0"])],
                "union",
                ["wide", "wall"],
                "1",
            ),
            &[],
            &[],
        ),
    ];
    let mut through = saddle_case(
        "through bore",
        saddle_steps(
            [
                wide(),
                cylinder("wall", "4", "60", LEVEL, ["1", "0", "-30"]),
            ],
            "subtract",
            ["wide", "wall"],
            "1",
        ),
        &[["1", "5", "89/10"], ["1", "5", "-89/10"]],
        &[["1", "81/20", "99/10"], ["1", "81/20", "-99/10"]],
    );
    through.pieces = 2;
    cases.push(through);
    let mut steps = saddle_steps([wide(), narrow()], "union", ["wide", "wall"], "1");
    steps.push(json!({"id":"flipped","op":"mirror","body":"round",
        "plane":{"point":["20","0","0"],"normal":["3","4","0"]},"keep_original":false}));
    let mut mirrored = saddle_case(
        "mirrored boss",
        steps,
        &[["5396/500", "8553/500", "10"]],
        &[["247/25", "421/25", "219/20"]],
    );
    mirrored.pieces = 2;
    cases.push(mirrored);
    for (name, kind, order, radius, sides, parts) in [
        (
            "loop boss",
            "union",
            ["wide", "wall"],
            "1",
            [1.0, 1.0],
            [wide(), past("8")],
        ),
        (
            "loop notch",
            "subtract",
            ["wide", "wall"],
            "1",
            [-1.0, 1.0],
            [wide(), past("8")],
        ),
        (
            "loop lens",
            "intersect",
            ["wide", "wall"],
            "1",
            [-1.0, -1.0],
            [wide(), past("8")],
        ),
        (
            "loop bite",
            "subtract",
            ["wall", "wide"],
            "1/2",
            [1.0, -1.0],
            [wide(), past("8")],
        ),
        (
            "loop notch r=3",
            "subtract",
            ["wide", "wall"],
            "3",
            [-1.0, 1.0],
            [wide(), past("8")],
        ),
        (
            "loop boss, wall first",
            "union",
            ["wall", "wide"],
            "1",
            [1.0, 1.0],
            [past("8"), wide()],
        ),
    ] {
        let ball = LoopBall {
            radius: Rational::parse(radius).unwrap().to_f64_nearest(),
            sides,
        };
        let (inside, outside) = ball.witnesses(1.0);
        cases.push(Case {
            name: name.into(),
            recipe: saddle_recipe(saddle_steps(parts, kind, order, radius)),
            inside,
            outside,
            brep: true,
            pieces: 1,
        });
    }
    let lean = cylinder("wall", "4", "60", LEAN, ["8", "-18", "-24"]);
    for (name, narrow, kind) in [
        ("loop on the far side", past("-8"), "union"),
        ("leaning loop boss", lean.clone(), "union"),
        ("leaning loop notch", lean, "subtract"),
        (
            "equal loop",
            cylinder("wall", "10", "60", LEVEL, ["1", "0", "-30"]),
            "union",
        ),
        ("loop almost past", past("27/2"), "union"),
    ] {
        cases.push(saddle_case(
            name,
            saddle_steps([wide(), narrow], kind, ["wide", "wall"], "1"),
            &[],
            &[],
        ));
    }
    let mut steps = saddle_steps([wide(), past("8")], "union", ["wide", "wall"], "1");
    steps.push(json!({"id":"flipped","op":"mirror","body":"round",
        "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
    let (inside, outside) = LoopBall {
        radius: 1.0,
        sides: [1.0, 1.0],
    }
    .witnesses(-1.0);
    cases.push(Case {
        name: "mirrored loop".into(),
        recipe: saddle_recipe(steps),
        inside,
        outside,
        brep: true,
        pieces: 2,
    });
    cases
}

// The equal-crossing file: cylinders of radius 8 along x and along z.
fn crossing_steps(kind: &str, radius: &str, second: &str, gap: &str) -> Vec<Value> {
    vec![
        json!({"id":"a","op":"cylinder","radius":"8","height":"32",
            "pose":{"rotation":["1","0","1","0"],"translation":["-16","0","0"]}}),
        json!({"id":"b","op":"cylinder","radius":second,"height":"32",
            "pose":{"rotation":["1","0","0","0"],"translation":["0",gap,"-16"]}}),
        json!({"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"}),
        json!({"id":"round","op":"fillet","body":"joined","radius":radius,"edges":[
            {"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
                        {"step":"b","role":{"side":{"segment":"wall"}}}]}]}),
    ]
}
fn corpus(id: &str, pieces: usize) -> Case {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../parts/fillets/cases.jsonl"
    ))
    .unwrap();
    let row: Value = text
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap();
    let mut inside = Vec::new();
    let mut outside = Vec::new();
    for witness in row["witnesses"].as_array().into_iter().flatten() {
        let point = std::array::from_fn(|i| q(witness["point"][i].as_str().unwrap()));
        match witness["state"].as_str().unwrap() {
            "TopAbs_IN" => inside.push(point),
            "TopAbs_OUT" => outside.push(point),
            other => panic!("{other}"),
        }
    }
    Case {
        name: id.into(),
        recipe: parse(row["recipe"].clone()),
        inside,
        outside,
        brep: row["expected"]["status"] == "built",
        pieces,
    }
}

fn crossing_cases() -> Vec<Case> {
    let mut cases = vec![
        saddle_case(
            "crossing union",
            crossing_steps("union", "1", "8", "0"),
            &[["81/10", "0", "81/10"]],
            &[["17/2", "0", "17/2"]],
        ),
        saddle_case(
            "crossing intersect",
            crossing_steps("intersect", "1", "8", "0"),
            &[["15/2", "0", "15/2"]],
            &[["79/10", "0", "79/10"]],
        ),
    ];
    for (name, rotation, translation, inside) in [
        (
            "oblique crossing",
            ["3", "0", "1", "0"],
            ["-96/5", "0", "-128/5"],
            &[][..],
        ),
        (
            "spatial crossing",
            ["2", "1", "1", "0"],
            ["-64/3", "64/3", "-32/3"],
            &[["2401/300", "1/150", "481/60"]][..],
        ),
    ] {
        let mut steps = crossing_steps("union", "3/2", "8", "0");
        steps[0]["height"] = json!("64");
        steps[0]["pose"]["translation"] = json!(["-32", "0", "0"]);
        steps[1]["height"] = json!("64");
        steps[1]["pose"] = json!({"rotation": rotation, "translation": translation});
        cases.push(saddle_case(name, steps, inside, &[]));
    }
    let mut steps = crossing_steps("union", "1", "8", "0");
    steps.push(json!({"id":"reflected","op":"mirror","body":"round",
        "plane":{"point":["5","0","0"],"normal":["1","0","0"]},"keep_original":false}));
    let mut reflected = saddle_case(
        "reflected crossing",
        steps,
        &[["197/100", "1", "803/100"]],
        &[["3/2", "1", "17/2"]],
    );
    reflected.pieces = 2;
    cases.push(reflected);
    // Offset equal axes cross in one loop; a narrower cylinder passes through
    // the wider one in two rings.
    for (name, second, gap, pieces) in [
        ("offset equal axes", "8", "1", 1),
        ("unequal crossing", "4", "0", 2),
    ] {
        let mut case = saddle_case(
            name,
            crossing_steps("union", "1", second, gap),
            &[["12", "0", "0"], ["0", "1", "12"]],
            &[["12", "12", "12"]],
        );
        case.pieces = pieces;
        cases.push(case);
    }
    for (id, pieces) in [("cross-equal", 1), ("cross-unequal", 2)] {
        cases.push(corpus(id, pieces));
    }
    cases
}

// The sphere-ring file: a radius-5 sphere and a radius-2 cylinder through it,
// its axis 1 off the centre, rounded at both rings.
fn sphere_cases() -> Vec<Case> {
    let mut ring = corpus("sphere-cylinder", 2);
    let mut mirrored = corpus("sphere-cylinder", 4);
    mirrored.name = "mirrored sphere ring".into();
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../parts/fillets/cases.jsonl"
    ))
    .unwrap();
    let mut recipe: Value = text
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == "sphere-cylinder")
        .unwrap()["recipe"]
        .clone();
    recipe["steps"].as_array_mut().unwrap().push(json!({
        "id": "flipped", "op": "mirror", "body": "round",
        "plane": {"point": ["0", "0", "0"], "normal": ["0", "1", "0"]},
        "keep_original": false
    }));
    recipe["result"] = "flipped".into();
    mirrored.recipe = parse(recipe);
    for p in [
        ["49/10", "0", "0"],
        ["-49/10", "0", "0"],
        ["0", "49/10", "0"],
    ] {
        mirrored.inside.push(p.map(q));
    }
    for p in [["51/10", "0", "0"], ["0", "-51/10", "0"]] {
        mirrored.outside.push(p.map(q));
    }
    // The B-rep builds both rings but cannot display their bands, so only
    // membership is checked.
    ring.brep = false;
    mirrored.brep = false;
    vec![ring, mirrored]
}

fn cases() -> Vec<Case> {
    let mut cases = coaxial_cases();
    cases.extend(saddle_cases());
    cases.extend(crossing_cases());
    cases.extend(sphere_cases());
    cases
}

fn location(plan: &Plan, p: &[Rational; 3]) -> Location {
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    plan.proven_checked(&point, &mut Workspace::default(), &mut decider())
        .unwrap()
        .location
}

#[test]
fn the_field_builds_every_step_three_fillet_and_agrees_on_membership() {
    let mut failures = Vec::new();
    for case in cases() {
        let start = std::time::Instant::now();
        let plan = match Plan::compile_checked(&case.recipe, &mut decider()) {
            Ok(plan) => plan,
            Err(e) => {
                failures.push(format!("{}: {e:?}", case.name));
                continue;
            }
        };
        let pieces = plan.rolling_pieces().count();
        if pieces != case.pieces {
            failures.push(format!("{}: {pieces} rolling pieces", case.name));
        }
        for (points, want) in [
            (&case.inside, Location::Inside),
            (&case.outside, Location::Outside),
        ] {
            for p in points {
                let found = location(&plan, p);
                if found != want {
                    failures.push(format!("{} at {p:?}: {found:?}", case.name));
                }
            }
        }
        eprintln!(
            "{}: {} witnesses in {:?}",
            case.name,
            case.inside.len() + case.outside.len(),
            start.elapsed()
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The cases whose B-rep blend vertices are checked, one test per group so
/// they run side by side.
const VERTEX_GROUPS: &[&[&str]] = &[
    &["frustum top", "frustum base", "flared top", "narrow top"],
    &["chamfer top", "chamfer side", "chamfer both"],
    &["dome rim", "sphere seat", "rod cap"],
    &["countersink throat"],
    &["saddle boss", "saddle bore"],
    &["saddle boss r=2, wall first", "saddle boss r=3/2"],
    &["leaning boss", "through bore", "mirrored boss"],
    &["loop boss", "loop notch", "loop lens"],
    &["loop bite", "loop notch r=3", "loop boss, wall first"],
    &[
        "loop on the far side",
        "leaning loop boss",
        "leaning loop notch",
    ],
    &["equal loop", "loop almost past", "mirrored loop"],
    &["crossing union", "crossing intersect", "reflected crossing"],
    &["oblique crossing", "spatial crossing"],
    &["offset equal axes", "unequal crossing"],
];
macro_rules! vertex_tests {
    ($($test:ident = $group:literal;)+) => {
        $(
            #[test]
            fn $test() {
                blend_vertices_lie_on_the_field_face(VERTEX_GROUPS[$group]);
            }
        )+
    };
}
vertex_tests! {
    frustum_blend_vertices_lie_on_the_field_face = 0;
    chamfer_blend_vertices_lie_on_the_field_face = 1;
    sphere_blend_vertices_lie_on_the_field_face = 2;
    countersink_blend_vertices_lie_on_the_field_face = 3;
    saddle_ring_blend_vertices_lie_on_the_field_face = 4;
    saddle_wide_ball_blend_vertices_lie_on_the_field_face = 5;
    saddle_posed_blend_vertices_lie_on_the_field_face = 6;
    loop_blend_vertices_lie_on_the_field_face = 7;
    loop_sides_blend_vertices_lie_on_the_field_face = 8;
    loop_pose_blend_vertices_lie_on_the_field_face = 9;
    loop_shape_blend_vertices_lie_on_the_field_face = 10;
    crossing_blend_vertices_lie_on_the_field_face = 11;
    oblique_crossing_blend_vertices_lie_on_the_field_face = 12;
    crossing_saddle_blend_vertices_lie_on_the_field_face = 13;
}

#[test]
fn every_built_case_has_its_vertices_checked() {
    let listed: std::collections::BTreeSet<String> = VERTEX_GROUPS
        .iter()
        .flat_map(|group| group.iter().map(|name| name.to_string()))
        .collect();
    let wanted: std::collections::BTreeSet<String> = cases()
        .into_iter()
        .filter(|c| c.brep)
        .map(|c| c.name)
        .collect();
    assert_eq!(listed, wanted);
}

/// Every B-rep display vertex of a blend face lies within the field's
/// guaranteed range of a field face, at the 1 nm own-part tolerance.
fn blend_vertices_lie_on_the_field_face(names: &[&str]) {
    // The vertex lies within the field's bound once a face witness within the
    // tolerance is proven, so the range need be no narrower than that.
    let own = q("1/1000000");
    for case in cases()
        .into_iter()
        .filter(|c| names.contains(&c.name.as_str()))
    {
        assert!(case.brep, "{}", case.name);
        let start = std::time::Instant::now();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let compiled = start.elapsed();
        let built = evaluate(&case.recipe, &mut decider())
            .unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        let brep = start.elapsed();
        let mesh = mesh_body(
            &built.model,
            built.result,
            &MeshOptions {
                chord_tolerance: q("1/20"),
                ..MeshOptions::default()
            },
            &mut decider(),
        )
        .unwrap()
        .value;
        // A fitted band carries its certified bound; its vertices lie within
        // the part's tolerance of the exact face, the larger of 1 nm and that
        // bound (rules T0 to T8).
        let mut indices: BTreeMap<usize, Rational> = BTreeMap::new();
        for face in &mesh.faces {
            if !built.model.face_origins(face.face).iter().any(|o| {
                let mut role = &o.role;
                while let Role::Copy { of, .. } = role {
                    role = &of.role;
                }
                matches!(role, Role::Blend { .. })
            }) {
                continue;
            }
            let bound = match built
                .model
                .geom()
                .surface(built.model.face(face.face).surface)
            {
                Surface::BSpline(SplineSupport::RollingSaddle(band)) => band.bound().upper.clone(),
                Surface::BSpline(SplineSupport::RollingSphere(band)) => band.bound().upper.clone(),
                _ => Rational::zero(),
            };
            let tolerance = own.clone().max(bound);
            for triangle in &mesh.triangles[face.triangles.clone()] {
                for v in triangle.vertices {
                    let entry = indices.entry(v).or_insert_with(|| tolerance.clone());
                    *entry = entry.clone().max(tolerance.clone());
                }
            }
        }
        let tolerances: Vec<Rational> = indices.values().cloned().collect();
        let points: Vec<Point3> = indices
            .keys()
            .map(|&v| {
                let [x, y, z] = mesh.positions[v].map(|c| Rational::from_f64_exact(c).unwrap());
                Point3::new(x.into(), y.into(), z.into())
            })
            .collect();
        assert!(!points.is_empty(), "{}: no blend vertices", case.name);
        let pieces: Vec<_> = plan.rolling_pieces().collect();
        let mut nearest = vec![None::<Rational>; points.len()];
        // Every vertex once on all workers; every 128th also at 1, 2 and 4
        // workers, whose records and work must agree. A range at a far piece
        // needs deep refinement, so each query has a work cap; an exhausted
        // query only yields no witness, never a wrong one.
        const CAP: u64 = 20_000;
        let allowance = |n: usize| OperationAllowance::reserve(n as u64 * CAP).unwrap();
        let mut used = 0;
        // A piece whose box lies farther than the tolerance from a vertex
        // cannot hold its witness.
        let near = |piece: &Piece, p: &Point3, tolerance: &Rational| {
            let Some(bounds) = piece.bounds() else {
                return true;
            };
            let [x, y, z] = [&p.x, &p.y, &p.z].map(|c| c.rational_ref().unwrap().clone());
            let gap = [x, y, z]
                .iter()
                .zip(bounds)
                .map(|(c, [lo, hi])| {
                    let g = if c < lo {
                        lo - c
                    } else if c > hi {
                        c - hi
                    } else {
                        Rational::zero()
                    };
                    &g * &g
                })
                .fold(Rational::zero(), |sum, g| sum + g);
            gap <= tolerance * tolerance
        };
        // Newest piece first: a mirror or pattern copy is compiled after its
        // original, and a vertex leaves the search once it has a witness.
        let mut widths: Vec<Rational> = tolerances.clone();
        widths.sort();
        widths.dedup();
        for (piece, width) in pieces
            .iter()
            .rev()
            .flat_map(|piece| widths.iter().map(move |w| (piece, w)))
        {
            let index: Vec<usize> = (0..points.len())
                .filter(|&i| {
                    &tolerances[i] == width
                        && nearest[i].as_ref().is_none_or(|u| u > width)
                        && near(piece, &points[i], width)
                })
                .collect();
            if index.is_empty() {
                continue;
            }
            let mine: Vec<Point3> = index.iter().map(|&i| points[i].clone()).collect();
            let sample: Vec<Point3> = mine.iter().step_by(128).cloned().collect();
            let all = distance_batch(
                *piece,
                &mine,
                width,
                None,
                allowance(mine.len()),
                &mut decider(),
            )
            .unwrap();
            let mut batches: Vec<Batch<Range>> = Vec::new();
            for workers in [1, 2, 4] {
                batches.push(
                    distance_batch(
                        *piece,
                        &sample,
                        width,
                        Some(workers),
                        allowance(sample.len()),
                        &mut decider(),
                    )
                    .unwrap(),
                );
            }
            assert_eq!(batches[0].records, batches[1].records, "{}", case.name);
            assert_eq!(batches[0].records, batches[2].records, "{}", case.name);
            for (&i, record) in index.iter().zip(&all.records) {
                let Ok(range) = &record.result else {
                    assert!(record.exhausted, "{}: {record:?}", case.name);
                    continue;
                };
                if range.upper.as_ref().is_some_and(|u| u <= width) {
                    used = used.max(record.work);
                }
                if let Some(upper) = &range.upper
                    && nearest[i].as_ref().is_none_or(|best| upper < best)
                {
                    nearest[i] = Some(upper.clone());
                }
            }
        }
        let queried = start.elapsed();
        eprintln!(
            "{}: field compile {compiled:?}, B-rep and mesh {:?}, queries {:?}, \
             most work for a witness {used}",
            case.name,
            brep - compiled,
            queried - brep
        );
        let mut worst = Rational::zero();
        let mut loosest = own.clone();
        for (found, tolerance) in nearest.iter().zip(&tolerances) {
            let found = found
                .clone()
                .unwrap_or_else(|| panic!("{}: no witness", case.name));
            assert!(
                &found <= tolerance,
                "{}: {found:?} over {tolerance:?}",
                case.name
            );
            worst = worst.max(found);
            loosest = loosest.max(tolerance.clone());
        }
        eprintln!(
            "{}: {} blend vertices, worst field distance {}, largest tolerance {}",
            case.name,
            points.len(),
            worst.to_f64_nearest(),
            loosest.to_f64_nearest()
        );
    }
}

/// The step 3 files' recipes that the field refuses, or that the B-rep
/// refuses.
fn refusals() -> Vec<(String, Recipe)> {
    let mut found = Vec::new();
    let mut add = |name: &str, steps: Vec<Value>| {
        found.push((name.to_string(), saddle_recipe(steps)));
    };
    // The saddle file's refusals.
    let bore =
        |radius: &str| saddle_steps([wide(), narrow()], "subtract", ["wide", "wall"], radius);
    add("bore r=3", bore("3"));
    add("bore r=10", bore("10"));
    add(
        "short narrow",
        saddle_steps(
            [wide(), cylinder("wall", "4", "11", LEVEL, ["1", "0", "0"])],
            "union",
            ["wide", "wall"],
            "1",
        ),
    );
    add(
        "short wide",
        saddle_steps(
            [
                cylinder("wide", "10", "12", FLAT, ["0", "6", "0"]),
                narrow(),
            ],
            "union",
            ["wide", "wall"],
            "3",
        ),
    );
    add(
        "clipped equal crossing",
        saddle_steps(
            [wide(), cylinder("wall", "10", "30", LEVEL, ["0", "0", "0"])],
            "union",
            ["wide", "wall"],
            "1",
        ),
    );
    add(
        "open overlap",
        saddle_steps(
            [wide(), cylinder("wall", "4", "30", LEVEL, ["8", "0", "0"])],
            "union",
            ["wide", "wall"],
            "1",
        ),
    );
    for (name, x, radius) in [
        ("neighbour in the way", "7", "7/4"),
        ("neighbour clear", "15/2", "2"),
    ] {
        add(
            name,
            vec![
                wide(),
                narrow(),
                json!({"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}),
                cylinder("next", radius, "30", LEVEL, [x, "0", "0"]),
                json!({"id":"both","op":"boolean","kind":"union","target":"joined","tool":"next"}),
                json!({"id":"round","op":"fillet","body":"both","radius":"1","edges":ring()}),
            ],
        );
    }
    add(
        "loop past the end",
        saddle_steps(
            [
                wide(),
                cylinder("wall", "4", "21", LEVEL, ["8", "0", "-21/2"]),
            ],
            "union",
            ["wide", "wall"],
            "1",
        ),
    );
    add(
        "loop offsets tangent",
        saddle_steps([wide(), past("8")], "subtract", ["wall", "wide"], "1"),
    );
    // The equal-crossing file's refusals.
    let mut clipped = crossing_steps("union", "3/2", "8", "0");
    clipped[1]["pose"] = json!({"rotation":["3","0","1","0"],"translation":["-48/5","0","-64/5"]});
    add("crossing clipped at its ends", clipped);
    add(
        "crossing with different ball sides",
        crossing_steps("subtract", "1", "8", "0"),
    );
    let mut window = crossing_steps("union", "1", "8", "0");
    let fillet = window.pop().unwrap();
    window.extend([
        json!({"id":"window","op":"block","size":["40","40","2"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-20","-20","2"]}}),
        json!({"id":"trimmed","op":"boolean","kind":"intersect","target":"joined","tool":"window"}),
    ]);
    let mut fillet = fillet;
    fillet["body"] = json!("trimmed");
    window.push(fillet);
    add("crossing clipped by a window", window);
    // The coaxial file's refusals and its torus groove.
    let groove = vec![
        json!({"id": "section", "op": "sketch",
            "plane": {"rotation": UP, "translation": ["0", "0", "0"]}, "loops": [[
            {"kind": "line", "id": "floor", "to": ["6", "0"]},
            {"kind": "arc", "id": "groove", "centre": ["10", "0"], "to": ["10", "4"], "turn": "cw"},
            {"kind": "line", "id": "rim", "to": ["10", "6"]},
            {"kind": "line", "id": "top", "to": ["0", "6"]},
            {"kind": "line", "id": "axis", "to": ["0", "0"]}
        ]]}),
        json!({"id": "r", "op": "revolve", "profile": "section",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "1", "0"]}, "angle": "360"}),
    ];
    let mut spike = block("b", "-20", "20");
    spike.push(
        json!({"id": "c", "op": "cone", "base_radius": "3", "top_radius": "0",
        "height": "4", "pose": pose("0")}),
    );
    spike.push(json!({"id": "j", "op": "boolean", "kind": "union", "target": "b", "tool": "c"}));
    let mut torus =
        vec![json!({"id": "t", "op": "torus", "major": "6", "minor": "2", "pose": pose("0")})];
    torus.extend(block("b", "1", "10"));
    torus
        .push(json!({"id": "j", "op": "boolean", "kind": "intersect", "target": "t", "tool": "b"}));
    let mouth = [face("b", "end"), hole("countersink")];
    let throat = [hole("countersink"), hole("wall")];
    for (name, steps, body, radius, edges) in [
        (
            "groove",
            groove,
            "r",
            "1",
            vec![
                [side("r", "floor"), side("r", "groove")],
                [side("r", "groove"), side("r", "rim")],
            ],
        ),
        (
            "countersink mouth",
            countersunk(),
            "h",
            "1",
            vec![mouth.clone()],
        ),
        (
            "countersink both",
            countersunk(),
            "h",
            "1",
            vec![mouth, throat],
        ),
        (
            "dome r=2",
            dome(),
            "j",
            "2",
            vec![[face("b", "start"), wall("s")]],
        ),
        (
            "dome r=3/2",
            dome(),
            "j",
            "3/2",
            vec![[face("b", "start"), wall("s")]],
        ),
        (
            "narrow frustum r=2",
            cone("5", "1", "6"),
            "c",
            "2",
            vec![[face("c", "end"), wall("c")]],
        ),
        (
            "spike r=12",
            spike,
            "j",
            "12",
            vec![[face("b", "end"), wall("c")]],
        ),
        (
            "torus slice r=3/5",
            torus,
            "j",
            "3/5",
            vec![[face("b", "start"), wall("t")]],
        ),
    ] {
        found.push((name.to_string(), rounded(&steps, body, radius, &edges)));
    }
    found
}

/// The field's answer on each of these, by name. Where the B-rep builds
/// (`neighbour clear`, the countersink's mouth) the field refusal is one of
/// step 3's remaining refusals; elsewhere both kernels refuse, each with its
/// own reason.
const REFUSED: &[(&str, &str)] = &[
    ("bore r=3", "field.rolling-path-reach"),
    ("bore r=10", "field.rolling-path-reach"),
    ("short narrow", "fillet.unsupported-overflow"),
    ("short wide", "fillet.unsupported-overflow"),
    ("clipped equal crossing", "field.fillet-contact"),
    ("open overlap", "fillet.lost-edge"),
    ("neighbour in the way", "fillet.unsupported-overflow"),
    ("neighbour clear", "fillet.unsupported-overflow"),
    ("loop past the end", "fillet.unsupported-overflow"),
    ("loop offsets tangent", "field.rolling-path-reach"),
    ("crossing clipped at its ends", "field.fillet-contact"),
    (
        "crossing with different ball sides",
        "field.rolling-path-reach",
    ),
    ("crossing clipped by a window", "fillet.lost-edge"),
    ("groove", "field.fillet-contact"),
    ("countersink mouth", "fillet.unsupported-overflow"),
    ("countersink both", "fillet.unsupported-overflow"),
    ("dome r=2", "field.rolling-path-reach"),
    ("dome r=3/2", "field.rolling-path-reach"),
    ("narrow frustum r=2", "field.rolling-path-reach"),
    ("spike r=12", "field.rolling-path-reach"),
    // Step 6 takes a plane with a torus. No ball wider than 1/2 fits this
    // section (the B-rep refuses with that radius limit), and the field's
    // offset plane misses the offset torus.
    ("torus slice r=3/5", "field.rolling-path-reach"),
];

#[test]
fn the_field_names_what_it_refuses() {
    let cases = refusals();
    assert_eq!(cases.len(), REFUSED.len());
    let mut wrong = Vec::new();
    for ((name, recipe), (listed, want)) in cases.iter().zip(REFUSED) {
        assert_eq!(name, listed);
        match Plan::compile_checked(recipe, &mut decider()) {
            Ok(plan) => wrong.push(format!(
                "{name}: builds, {} pieces",
                plan.rolling_pieces().count()
            )),
            Err(e) if e.code != *want => wrong.push(format!("{name}: {} ({})", e.code, e.detail)),
            Err(e) => eprintln!("{name}: {} ({})", e.code, e.detail),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
