//! Fillets on circles where two surfaces of revolution about one axis meet,
//! beyond a plane and a cylinder: frustum rims, chamfer and countersink
//! edges, a dome's rim, a seated sphere and a revolved groove. Each strip is
//! an exact torus. Its volume is checked against an independent section
//! integral turned about the axis, its sides against exact witnesses, and its
//! display, names and STEP file as for the other elementary strips.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use catallax_boolean::{Containment, classify_point};
use catallax_features::{
    BuildError, Built, evaluate, evaluate_incremental,
    fillet::{RADIUS_TOO_LARGE, radius_limit},
};
use catallax_field::Plan;
use catallax_measure::{Accuracy, face_area, mass_properties, volume};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Refusal, with_operation_log};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use catallax_step::write_step;
use serde_json::{Value, json};

const UP: [&str; 4] = ["1", "0", "0", "0"];

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
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
/// A cylinder of radius 5 and height 10 with its top rim chamfered by 1.
fn chamfered() -> Vec<Value> {
    vec![
        json!({"id": "c", "op": "cylinder", "radius": "5", "height": "10", "pose": pose("0")}),
        json!({"id": "bevel", "op": "chamfer", "body": "c", "distance": "1",
            "edges": [{"between": [face("c", "end"), wall("c")]}]}),
    ]
}
/// The sphere of radius 5 above the plane z = 2.
fn dome() -> Vec<Value> {
    let mut steps = vec![json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("0")})];
    steps.extend(block("b", "2", "10"));
    steps
        .push(json!({"id": "j", "op": "boolean", "kind": "intersect", "target": "s", "tool": "b"}));
    steps
}
/// A sphere of radius 5 centred 3 above a block's top.
fn seat() -> Vec<Value> {
    let mut steps = block("b", "0", "20");
    steps.push(json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("23")}));
    steps.push(json!({"id": "j", "op": "boolean", "kind": "union", "target": "b", "tool": "s"}));
    steps
}
/// A through hole of diameter 6 with a 90 degree countersink 12 across.
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
fn revolved(loop_: Value) -> Vec<Value> {
    vec![
        json!({"id": "section", "op": "sketch",
            "plane": {"rotation": UP, "translation": ["0", "0", "0"]}, "loops": [loop_]}),
        json!({"id": "r", "op": "revolve", "profile": "section",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "1", "0"]}, "angle": "360"}),
    ]
}
/// A disc of radius 10 and height 6 about the y axis, with a quarter-round
/// groove of radius 4 cut round its lower rim.
fn groove() -> Vec<Value> {
    revolved(json!([
        {"kind": "line", "id": "floor", "to": ["6", "0"]},
        {"kind": "arc", "id": "groove", "centre": ["10", "0"], "to": ["10", "4"], "turn": "cw"},
        {"kind": "line", "id": "rim", "to": ["10", "6"]},
        {"kind": "line", "id": "top", "to": ["0", "6"]},
        {"kind": "line", "id": "axis", "to": ["0", "0"]}
    ]))
}
/// A rod of radius 6 and height 8 about the y axis, capped by the sphere of
/// radius 10 about the origin.
fn capped_rod() -> Vec<Value> {
    revolved(json!([
        {"kind": "line", "id": "base", "to": ["6", "0"]},
        {"kind": "line", "id": "rod", "to": ["6", "8"]},
        {"kind": "arc", "id": "cap", "centre": ["0", "0"], "to": ["0", "10"], "turn": "ccw"},
        {"kind": "line", "id": "axis", "to": ["0", "0"]}
    ]))
}

fn recipe(steps: &[Value], body: &str, round: Option<(&str, Vec<[Value; 2]>)>) -> Recipe {
    let mut steps = steps.to_vec();
    let mut result = body;
    if let Some((radius, edges)) = round {
        let edges: Vec<_> = edges.into_iter().map(|e| json!({"between": e})).collect();
        steps.push(
            json!({"id": "round", "op": "fillet", "body": body, "radius": radius,
            "edges": edges}),
        );
        result = "round";
    }
    Recipe::parse(
        &json!({"format": "catallax-recipe/3", "name": "coaxial round", "steps": steps,
            "result": result})
        .to_string(),
    )
    .unwrap()
}
fn build(recipe: &Recipe) -> Result<Built, Refusal> {
    evaluate(recipe, &mut decider()).map_err(|error| match error {
        BuildError::Geometry { refusal, .. } => refusal,
        BuildError::Recipe(error) => panic!("invalid test recipe: {error:?}"),
    })
}

// The half-plane through the axis: distance from the axis, then height.
type Pair = [f64; 2];

/// A support's profile there: a line `n . p = d` with unit `n`, or a circle.
#[derive(Clone, Copy)]
enum Profile {
    Line(Pair, f64),
    Circle(Pair, f64),
}

/// `∫ ρ²/2 dz` along a profile piece from `p` to `q`, the shorter way round a
/// circle. Summed round a closed boundary it is the region's first moment
/// about the axis, so turning the region sweeps `2 π` times it.
fn moment(profile: Profile, p: Pair, q: Pair) -> f64 {
    match profile {
        Profile::Line(..) => {
            let (dr, dz) = (q[0] - p[0], q[1] - p[1]);
            dz / 2.0 * (p[0] * p[0] + p[0] * dr + dr * dr / 3.0)
        }
        Profile::Circle(c, radius) => {
            let angle = |p: Pair| (p[1] - c[1]).atan2(p[0] - c[0]);
            let start = angle(p);
            let turn = (angle(q) - start + PI).rem_euclid(2.0 * PI) - PI;
            let f = |t: f64| {
                radius / 2.0
                    * (c[0] * c[0] * t.sin()
                        + 2.0 * c[0] * radius * (t / 2.0 + (2.0 * t).sin() / 4.0)
                        + radius * radius * (t.sin() - t.sin().powi(3) / 3.0))
            };
            f(start + turn) - f(start)
        }
    }
}

/// A rolling ball worked by hand in the half-plane: the edge point, the two
/// profiles, the ball centre and radius. Feet follow from the profiles.
struct Ball {
    edge: Pair,
    profiles: [Profile; 2],
    centre: Pair,
    radius: f64,
}
impl Ball {
    fn foot(&self, index: usize) -> Pair {
        let c = self.centre;
        match self.profiles[index] {
            Profile::Line(n, d) => {
                let off = n[0] * c[0] + n[1] * c[1] - d;
                [c[0] - off * n[0], c[1] - off * n[1]]
            }
            Profile::Circle(m, radius) => {
                let (dx, dz) = (c[0] - m[0], c[1] - m[1]);
                let length = (dx * dx + dz * dz).sqrt();
                [m[0] + dx * radius / length, m[1] + dz * radius / length]
            }
        }
    }
    /// The volume between the edge and the ball, turned about the axis.
    fn region(&self) -> f64 {
        let [a, b] = [self.foot(0), self.foot(1)];
        let arc = Profile::Circle(self.centre, self.radius);
        let total = moment(self.profiles[0], self.edge, a)
            + moment(arc, a, b)
            + moment(self.profiles[1], b, self.edge);
        2.0 * PI * total.abs()
    }
    /// The area the shorter arc between the feet sweeps about the axis.
    fn band_area(&self) -> f64 {
        let [a, b] = [self.foot(0), self.foot(1)];
        let c = self.centre;
        let start = (a[1] - c[1]).atan2(a[0] - c[0]);
        let turn = ((b[1] - c[1]).atan2(b[0] - c[0]) - start + PI).rem_euclid(2.0 * PI) - PI;
        let (low, high) = if turn > 0.0 {
            (start, start + turn)
        } else {
            (start + turn, start)
        };
        2.0 * PI * self.radius * (c[0] * (high - low) + self.radius * (high.sin() - low.sin()))
    }
    /// A point inside the ball, and one in the corner between it and the edge.
    fn witnesses(&self) -> [Pair; 2] {
        let (c, e) = (self.centre, self.edge);
        let (dx, dz) = (e[0] - c[0], e[1] - c[1]);
        let length = (dx * dx + dz * dz).sqrt();
        let near = [
            c[0] + dx * self.radius / length,
            c[1] + dz * self.radius / length,
        ];
        let half = self.radius / 2.0 / length;
        let inside = [c[0] + dx * half, c[1] + dz * half];
        [inside, [(near[0] + e[0]) / 2.0, (near[1] + e[1]) / 2.0]]
    }
}

/// Where two lines `n . p = d` meet.
fn meet(a: Pair, d: f64, b: Pair, e: f64) -> Pair {
    let det = a[0] * b[1] - a[1] * b[0];
    [(d * b[1] - e * a[1]) / det, (a[0] * e - b[0] * d) / det]
}
fn unit(v: Pair) -> Pair {
    let length = (v[0] * v[0] + v[1] * v[1]).sqrt();
    [v[0] / length, v[1] / length]
}

/// One selected edge and the ball that rounds it.
struct Strip {
    edge: [Value; 2],
    ball: Ball,
}
/// A body, its rounded edges, its exact volume and whether rounding removes
/// material, with the half-plane's axis: the z axis or the y axis.
struct Case {
    name: &'static str,
    steps: Vec<Value>,
    body: &'static str,
    radius: &'static str,
    strips: Vec<Strip>,
    volume: f64,
    convex: bool,
    about_y: bool,
}

fn frustum_top() -> Case {
    let n = unit([3.0, 1.0]);
    let d = 15.0 / 10f64.sqrt();
    Case {
        name: "frustum top",
        steps: cone("5", "3", "6"),
        body: "c",
        radius: "1",
        strips: vec![Strip {
            edge: [face("c", "end"), wall("c")],
            ball: Ball {
                edge: [3.0, 6.0],
                profiles: [Profile::Line([0.0, 1.0], 6.0), Profile::Line(n, d)],
                centre: meet(n, d - 1.0, [0.0, 1.0], 5.0),
                radius: 1.0,
            },
        }],
        volume: 98.0 * PI,
        convex: true,
        about_y: false,
    }
}
// The acute base rim: its band runs across the torus's outer equator.
fn frustum_base() -> Case {
    let n = unit([3.0, 1.0]);
    let d = 15.0 / 10f64.sqrt();
    Case {
        name: "frustum base",
        steps: cone("5", "3", "6"),
        body: "c",
        radius: "1",
        strips: vec![Strip {
            edge: [face("c", "start"), wall("c")],
            ball: Ball {
                edge: [5.0, 0.0],
                profiles: [Profile::Line([0.0, 1.0], 0.0), Profile::Line(n, d)],
                centre: meet(n, d - 1.0, [0.0, 1.0], 1.0),
                radius: 1.0,
            },
        }],
        volume: 98.0 * PI,
        convex: true,
        about_y: false,
    }
}
fn flared_top() -> Case {
    let n = unit([3.0, -1.0]);
    let d = 9.0 / 10f64.sqrt();
    Case {
        name: "flared top",
        steps: cone("3", "5", "6"),
        body: "c",
        radius: "1",
        strips: vec![Strip {
            edge: [face("c", "end"), wall("c")],
            ball: Ball {
                edge: [5.0, 6.0],
                profiles: [Profile::Line([0.0, 1.0], 6.0), Profile::Line(n, d)],
                centre: meet(n, d - 1.0, [0.0, 1.0], 5.0),
                radius: 1.0,
            },
        }],
        volume: 98.0 * PI,
        convex: true,
        about_y: false,
    }
}
// The circle of ball centres is narrower than the ball: an outer torus sheet.
fn narrow_top() -> Case {
    let n = unit([2.0, 1.0]);
    let d = 8.0 / 5f64.sqrt();
    Case {
        name: "narrow top",
        steps: cone("4", "1", "6"),
        body: "c",
        radius: "1",
        strips: vec![Strip {
            edge: [face("c", "end"), wall("c")],
            ball: Ball {
                edge: [1.0, 6.0],
                profiles: [Profile::Line([0.0, 1.0], 6.0), Profile::Line(n, d)],
                centre: meet(n, d - 1.0, [0.0, 1.0], 5.0),
                radius: 1.0,
            },
        }],
        volume: 42.0 * PI,
        convex: true,
        about_y: false,
    }
}
// Both rounds at once take the smaller ball: at radius 1/2 their complete
// tori touch, although their bands lie apart on the chamfer.
fn chamfer_strips(top: bool, side: bool) -> Case {
    let n = unit([1.0, 1.0]);
    let d = 14.0 / 2f64.sqrt();
    let r = if top && side { 0.25 } else { 0.5 };
    let mut strips = Vec::new();
    if top {
        strips.push(Strip {
            edge: [face("c", "end"), json!({"step": "bevel"})],
            ball: Ball {
                edge: [4.0, 10.0],
                profiles: [Profile::Line([0.0, 1.0], 10.0), Profile::Line(n, d)],
                centre: meet(n, d - r, [0.0, 1.0], 10.0 - r),
                radius: r,
            },
        });
    }
    if side {
        strips.push(Strip {
            edge: [wall("c"), json!({"step": "bevel"})],
            ball: Ball {
                edge: [5.0, 9.0],
                profiles: [Profile::Line([1.0, 0.0], 5.0), Profile::Line(n, d)],
                centre: meet(n, d - r, [1.0, 0.0], 5.0 - r),
                radius: r,
            },
        });
    }
    Case {
        name: match (top, side) {
            (true, false) => "chamfer top",
            (false, true) => "chamfer side",
            _ => "chamfer both",
        },
        steps: chamfered(),
        body: "bevel",
        radius: if top && side { "1/4" } else { "1/2" },
        strips,
        volume: 250.0 * PI - 14.0 * PI / 3.0,
        convex: true,
        about_y: false,
    }
}
// A dome's acute rim: its band also runs across the torus's outer equator.
fn dome_rim(radius: &'static str) -> Case {
    let r: f64 = radius.parse().unwrap();
    let height = 2.0 + r;
    Case {
        name: "dome rim",
        steps: dome(),
        body: "j",
        radius,
        strips: vec![Strip {
            edge: [face("b", "start"), wall("s")],
            ball: Ball {
                edge: [21f64.sqrt(), 2.0],
                profiles: [
                    Profile::Line([0.0, 1.0], 2.0),
                    Profile::Circle([0.0, 0.0], 5.0),
                ],
                centre: [((5.0 - r).powi(2) - height * height).sqrt(), height],
                radius: r,
            },
        }],
        volume: 36.0 * PI,
        convex: true,
        about_y: false,
    }
}
fn sphere_seat() -> Case {
    Case {
        name: "sphere seat",
        steps: seat(),
        body: "j",
        radius: "1",
        strips: vec![Strip {
            edge: [face("b", "end"), wall("s")],
            ball: Ball {
                edge: [4.0, 20.0],
                profiles: [
                    Profile::Line([0.0, 1.0], 20.0),
                    Profile::Circle([0.0, 23.0], 5.0),
                ],
                centre: [32f64.sqrt(), 21.0],
                radius: 1.0,
            },
        }],
        volume: 32000.0 + 448.0 * PI / 3.0,
        convex: false,
        about_y: false,
    }
}
fn countersink(mouth: bool, throat: bool) -> Case {
    // The countersink cone rho - z = -14, from the mouth at (6, 20) down to
    // the wall of radius 3 at height 17. Its material lies away from the axis.
    let n = unit([1.0, -1.0]);
    let d = -14.0 / 2f64.sqrt();
    let mut strips = Vec::new();
    if mouth {
        strips.push(Strip {
            edge: [face("b", "end"), hole("countersink")],
            ball: Ball {
                edge: [6.0, 20.0],
                profiles: [Profile::Line([0.0, 1.0], 20.0), Profile::Line(n, d)],
                centre: meet(n, d + 1.0, [0.0, 1.0], 19.0),
                radius: 1.0,
            },
        });
    }
    if throat {
        strips.push(Strip {
            edge: [hole("countersink"), hole("wall")],
            ball: Ball {
                edge: [3.0, 17.0],
                profiles: [Profile::Line(n, d), Profile::Line([1.0, 0.0], 3.0)],
                centre: meet(n, d + 1.0, [1.0, 0.0], 4.0),
                radius: 1.0,
            },
        });
    }
    Case {
        name: match (mouth, throat) {
            (true, false) => "countersink mouth",
            (false, true) => "countersink throat",
            _ => "countersink both",
        },
        steps: countersunk(),
        body: "h",
        radius: "1",
        strips,
        volume: 32000.0 - 216.0 * PI,
        convex: true,
        about_y: false,
    }
}
fn groove_strips() -> Case {
    let tube = Profile::Circle([10.0, 0.0], 4.0);
    Case {
        name: "groove",
        steps: groove(),
        body: "r",
        radius: "1",
        strips: vec![
            Strip {
                edge: [side("r", "floor"), side("r", "groove")],
                ball: Ball {
                    edge: [6.0, 0.0],
                    profiles: [Profile::Line([0.0, 1.0], 0.0), tube],
                    centre: [10.0 - 24f64.sqrt(), 1.0],
                    radius: 1.0,
                },
            },
            Strip {
                edge: [side("r", "groove"), side("r", "rim")],
                ball: Ball {
                    edge: [10.0, 4.0],
                    profiles: [tube, Profile::Line([1.0, 0.0], 10.0)],
                    centre: [9.0, 24f64.sqrt()],
                    radius: 1.0,
                },
            },
        ],
        volume: 600.0 * PI - 80.0 * PI * PI + 128.0 * PI / 3.0,
        convex: true,
        about_y: true,
    }
}
fn rod_cap() -> Case {
    Case {
        name: "rod cap",
        steps: capped_rod(),
        body: "r",
        radius: "1",
        strips: vec![Strip {
            edge: [side("r", "rod"), side("r", "cap")],
            ball: Ball {
                edge: [6.0, 8.0],
                profiles: [
                    Profile::Line([1.0, 0.0], 6.0),
                    Profile::Circle([0.0, 0.0], 10.0),
                ],
                centre: [5.0, 56f64.sqrt()],
                radius: 1.0,
            },
        }],
        volume: 288.0 * PI + 112.0 * PI / 3.0,
        convex: true,
        about_y: true,
    }
}
fn cases() -> Vec<Case> {
    vec![
        frustum_top(),
        frustum_base(),
        flared_top(),
        narrow_top(),
        chamfer_strips(true, false),
        chamfer_strips(false, true),
        chamfer_strips(true, true),
        dome_rim("1"),
        sphere_seat(),
        countersink(true, false),
        countersink(false, true),
        countersink(true, true),
        groove_strips(),
        rod_cap(),
    ]
}
fn rounded(case: &Case) -> Recipe {
    let edges = case.strips.iter().map(|s| s.edge.clone()).collect();
    recipe(&case.steps, case.body, Some((case.radius, edges)))
}
fn built(case: &Case) -> Built {
    let built = build(&rounded(case)).unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
    assert_eq!(built.model.check(built.result), vec![], "{}", case.name);
    assert_eq!(
        built.model.body(built.result).lumps.len(),
        1,
        "{}",
        case.name
    );
    assert_eq!(built.fillets.len(), case.strips.len(), "{}", case.name);
    built
}
fn fine() -> Accuracy {
    Accuracy {
        relative: q("1/1000000000"),
        absolute: q("1/1000000000"),
    }
}
fn at(case: &Case, p: Pair) -> RPoint3 {
    let rational = |x: f64| q(&format!("{x:.9}"));
    if case.about_y {
        RPoint3::new(rational(p[0]), rational(p[1]), Rational::zero())
    } else {
        RPoint3::new(rational(p[0]), Rational::zero(), rational(p[1]))
    }
}

#[test]
fn each_strip_is_the_exact_torus_of_its_ball() {
    for case in cases() {
        let built = built(&case);
        for (strip, definition) in case.strips.iter().zip(&built.fillets) {
            assert_eq!(definition.convex, case.convex, "{}", case.name);
            let (centre, major, minor, ring) = match &definition.surface {
                Surface::Torus(t) => (&t.center, &t.major, &t.minor, true),
                Surface::TorusSheet(t) => (&t.center, &t.major, &t.minor, false),
                other => panic!("{}: carrier {:?}", case.name, other.kind()),
            };
            let ball = &strip.ball;
            let height = if case.about_y {
                centre.y.to_f64()
            } else {
                centre.z.to_f64()
            };
            assert!((height - ball.centre[1]).abs() < 1e-12, "{}", case.name);
            assert!(
                (major.to_f64() - ball.centre[0]).abs() < 1e-12,
                "{}",
                case.name
            );
            assert!(
                (minor.to_f64() - ball.radius).abs() < 1e-12,
                "{}",
                case.name
            );
            assert_eq!(ring, ball.centre[0] > ball.radius, "{}", case.name);
            // The published face carries this definition's carrier.
            let face = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .find(|f| built.model.face_origins(*f).contains(&definition.origin))
                .unwrap();
            assert!(
                built
                    .model
                    .geom()
                    .surface(built.model.face(face).surface)
                    .same_definition(&definition.surface),
                "{}",
                case.name
            );
        }
    }
}

#[test]
fn each_round_changes_the_volume_by_its_turned_section() {
    for case in cases() {
        let built = built(&case);
        let change: f64 = case.strips.iter().map(|s| s.ball.region()).sum();
        let want = if case.convex {
            case.volume - change
        } else {
            case.volume + change
        };
        let found = volume(&built.model, built.result, &fine(), &mut decider()).unwrap();
        assert!(found.marks.is_empty(), "{}", case.name);
        let (lo, hi) = (
            found.value.lo.to_f64_nearest(),
            found.value.hi.to_f64_nearest(),
        );
        assert!(
            lo - 1e-7 <= want && want <= hi + 1e-7 && hi - lo < 1e-6,
            "{}: [{lo}, {hi}] does not hold {want}",
            case.name
        );
        // Each band's area, and the body's centroid on its axis.
        for (strip, definition) in case.strips.iter().zip(&built.fillets) {
            let face = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .find(|f| built.model.face_origins(*f).contains(&definition.origin))
                .unwrap();
            let area = face_area(&built.model, face, &fine(), &mut decider()).unwrap();
            assert!(area.marks.is_empty(), "{}", case.name);
            let mid = ((&area.value.lo + &area.value.hi) * q("1/2")).to_f64_nearest();
            let want = strip.ball.band_area();
            assert!(
                (mid - want).abs() < 1e-6,
                "{}: area {mid} against {want}",
                case.name
            );
        }
        let properties = mass_properties(&built.model, built.result, &fine(), &mut decider())
            .unwrap()
            .value;
        let centre = properties.centroid.unwrap();
        let off_axis = if case.about_y { [0, 2] } else { [0, 1] };
        for axis in off_axis {
            let mid = ((&centre[axis].lo + &centre[axis].hi) * q("1/2")).to_f64_nearest();
            assert!(mid.abs() < 1e-9, "{}: centroid {axis} at {mid}", case.name);
        }
    }
}

#[test]
fn witnesses_lie_on_the_material_side_of_each_round() {
    for case in cases() {
        let mut built = built(&case);
        for strip in &case.strips {
            let [ball, corner] = strip.ball.witnesses();
            let (solid, void) = if case.convex {
                (ball, corner)
            } else {
                (corner, ball)
            };
            for (p, want) in [(solid, Containment::Inside), (void, Containment::Outside)] {
                let point = at(&case, p);
                let found = classify_point(&mut built.model, built.result, &point, &mut decider());
                assert!(
                    matches!(found, Outcome::Proven(c) if c == want),
                    "{} at {p:?}: wanted {want:?}, got {found:?}",
                    case.name
                );
            }
        }
    }
}

fn assert_closed(name: &str, mesh: &Mesh, options: &MeshOptions) {
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(
        incidence.values().all(|v| *v == (2, 0)),
        "{name}: open or misoriented mesh"
    );
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
            assert_eq!(
                incidence.get(&key),
                Some(&(2, 0)),
                "{name}: rail not shared"
            );
        }
    }
    assert!(
        mesh.faces
            .iter()
            .all(|f| f.max_chord_error <= options.chord_tolerance),
        "{name}: chord bound"
    );
    if let Some(limit) = &options.angular_deviation {
        assert!(
            mesh.faces
                .iter()
                .all(|f| f.max_angular_deviation.as_ref().is_some_and(|a| a <= limit)),
            "{name}: angular bound"
        );
    }
}

// At kernel defaults, and at SolidDraft's automatic viewport levels: 15
// degrees and `0.002 * diagonal * factor`, halved for binary32 rounding.
#[test]
fn each_rounded_body_displays_closed_at_every_viewport_level() {
    for case in cases() {
        let built = built(&case);
        let defaults = MeshOptions::default();
        let mesh = mesh_body(&built.model, built.result, &defaults, &mut decider())
            .unwrap_or_else(|e| panic!("{}: {e:?}", case.name))
            .value;
        assert_closed(case.name, &mesh, &defaults);
        let extent = |axis: usize| {
            let values = mesh.positions.iter().map(|p| p[axis]);
            values.clone().fold(f64::MIN, f64::max) - values.fold(f64::MAX, f64::min)
        };
        let diagonal = (0..3).map(|a| extent(a).powi(2)).sum::<f64>().sqrt();
        for factor in [4.0, 1.0, 0.25] {
            let chord = 0.002 * diagonal * factor / 2.0;
            let options = MeshOptions {
                chord_tolerance: q(&format!("{chord:.9}")),
                angular_deviation: Some(q("0.2617993877991494")),
                ..MeshOptions::default()
            };
            let mesh = mesh_body(&built.model, built.result, &options, &mut decider())
                .unwrap_or_else(|e| panic!("{} at {factor}: {e:?}", case.name))
                .value;
            assert_closed(case.name, &mesh, &options);
        }
    }
}

fn blend_names(built: &Built) -> Vec<catallax_recipe::Origin> {
    let mut names: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f).to_vec())
        .filter(|o| matches!(o.role, Role::Blend { .. }))
        .collect();
    names.sort();
    names
}

#[test]
fn radius_and_upstream_edits_keep_the_blend_names() {
    // A radius edit on the chamfer's two rounds.
    let case = chamfer_strips(true, true);
    let first = blend_names(&built(&case));
    assert_eq!(first.len(), 2);
    let edges: Vec<_> = case.strips.iter().map(|s| s.edge.clone()).collect();
    let smaller = build(&recipe(
        &case.steps,
        case.body,
        Some(("1/5", edges.clone())),
    ))
    .unwrap();
    assert_eq!(first, blend_names(&smaller));
    // A wider chamfer moves both edges; the rounds keep their names.
    let mut steps = case.steps.clone();
    steps[1]["distance"] = json!("3/2");
    let wider = build(&recipe(&steps, case.body, Some(("1/4", edges)))).unwrap();
    assert_eq!(wider.model.check(wider.result), vec![]);
    assert_eq!(first, blend_names(&wider));
    // A narrower frustum top.
    let case = frustum_top();
    let first = blend_names(&built(&case));
    let edges: Vec<_> = case.strips.iter().map(|s| s.edge.clone()).collect();
    let narrower = build(&recipe(&cone("5", "5/2", "6"), "c", Some(("1", edges)))).unwrap();
    assert_eq!(first, blend_names(&narrower));
}

#[test]
fn a_cached_round_records_its_contact_phase() {
    let case = frustum_top();
    let prefix =
        evaluate_incremental(&recipe(&case.steps, case.body, None), None, &mut decider()).unwrap();
    let mut dec = decider();
    let (result, phases) =
        with_operation_log(|| evaluate_incremental(&rounded(&case), Some(prefix), &mut dec));
    let built = result.unwrap().into_built();
    assert_eq!(built.model.check(built.result), vec![]);
    let contact: Vec<_> = phases
        .iter()
        .filter(|p| p.operation == "fillet.contact")
        .collect();
    assert!(!contact.is_empty(), "fillet.contact was not recorded");
    let total: f64 = contact.iter().map(|p| p.wall_ms).sum();
    eprintln!("fillet.contact={total:.3}ms over {} spans", contact.len());
}

// Bands beside a cone are written on their own seams, and a band across its
// torus's outer equator is cut there. The files of spherical and toroidal
// supports fit spline pcurves for several seconds each, so the STEP
// regressions read those back instead, with every other case, raw.
#[test]
fn rounded_bodies_write_step_with_no_marks() {
    for case in cases() {
        if !case
            .steps
            .iter()
            .all(|step| step["op"] != "sphere" && step["op"] != "revolve")
        {
            continue;
        }
        let built = built(&case);
        let mut bytes = Vec::new();
        let answer = write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut decider(),
            &mut bytes,
        )
        .unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        if case.name == "narrow top" {
            let report = &answer.value;
            assert_eq!(
                (
                    report.solids,
                    report.faces,
                    report.edges,
                    report.vertices,
                    report.approximated_curves
                ),
                (1, 4, 5, 3, 0)
            );
            assert!(report.regularizations.is_empty());
        }
        assert!(answer.marks.is_empty(), "{}", case.name);
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("TOROIDAL_SURFACE"), "{}", case.name);
        assert!(text.contains("SEAM_CURVE"), "{}", case.name);
    }
}

#[test]
fn the_field_still_refuses_curved_rounds_beyond_the_cylinder() {
    for case in [frustum_top(), dome_rim("1"), groove_strips()] {
        let refusal = Plan::compile_checked(&rounded(&case), &mut decider()).unwrap_err();
        assert_eq!(refusal.code, "field.fillet-contact", "{}", case.name);
    }
}

fn refusal(steps: &[Value], body: &str, radius: &str, edge: [Value; 2]) -> Refusal {
    build(&recipe(steps, body, Some((radius, vec![edge])))).expect_err("expected a refusal")
}

#[test]
fn unsupported_curved_rounds_refuse_with_their_own_codes() {
    // The dome is 3 high, so no ball of radius 2 rests on both its faces.
    let found = refusal(&dome(), "j", "2", [face("b", "start"), wall("s")]);
    assert_eq!(found.code, RADIUS_TOO_LARGE);
    assert_eq!(radius_limit(&found), Some(q("3/2")));
    // At the limit itself the two offsets only touch.
    let found = refusal(&dome(), "j", "3/2", [face("b", "start"), wall("s")]);
    assert_eq!(found.code, RADIUS_TOO_LARGE);
    assert!(found.detail.contains("regular contact limit"), "{found:?}");
    assert_eq!(radius_limit(&found), Some(q("3/2")));
    // A frustum too narrow for the ball: its centre would cross the axis.
    let found = refusal(
        &cone("5", "1", "6"),
        "c",
        "2",
        [face("c", "end"), wall("c")],
    );
    assert_eq!(found.code, "fillet.unsupported-contact");
    assert!(
        found.detail.contains("centre reaches the axis"),
        "{found:?}"
    );
    // A pointed cone on a block: past its apex the cone contact would cross
    // the axis.
    let mut spike = block("b", "-20", "20");
    spike.push(
        json!({"id": "c", "op": "cone", "base_radius": "3", "top_radius": "0",
        "height": "4", "pose": pose("0")}),
    );
    spike.push(json!({"id": "j", "op": "boolean", "kind": "union", "target": "b", "tool": "c"}));
    let found = refusal(&spike, "j", "12", [face("b", "end"), wall("c")]);
    assert_eq!(found.code, "fillet.unsupported-contact");
    assert!(
        found.detail.contains("contact reaches the axis"),
        "{found:?}"
    );
    // A round wide enough that its complete torus reaches the cylinder below
    // the chamfer: the guard checks the carrier, not yet the band.
    let found = refusal(
        &chamfered(),
        "bevel",
        "2",
        [face("c", "end"), json!({"step": "bevel"})],
    );
    assert_eq!(found.code, "fillet.unsupported-overflow");
    // A plane tilted through a cone meets it off the axis.
    let mut tilted = cone("5", "3", "6");
    tilted.push(json!({"id": "cut-section", "op": "sketch",
        "plane": {"rotation": ["24/25", "0", "7/25", "0"], "translation": ["0", "0", "4"]},
        "loops": [[{"kind": "line", "id": "s", "to": ["30", "-30"]},
            {"kind": "line", "id": "e", "to": ["30", "30"]},
            {"kind": "line", "id": "n", "to": ["-30", "30"]},
            {"kind": "line", "id": "w", "to": ["-30", "-30"]}]]}));
    tilted.push(json!({"id": "cut", "op": "extrude", "profile": "cut-section", "distance": "30"}));
    tilted.push(
        json!({"id": "j", "op": "boolean", "kind": "subtract", "target": "c", "tool": "cut"}),
    );
    let found = refusal(&tilted, "j", "1/2", [face("cut", "start"), wall("c")]);
    assert_eq!(found.code, "fillet.unsupported-contact");
    // A plane through a torus meets it on circles kept as intersection curves.
    let mut ring =
        vec![json!({"id": "t", "op": "torus", "major": "6", "minor": "2", "pose": pose("0")})];
    ring.extend(block("b", "1", "10"));
    ring.push(json!({"id": "j", "op": "boolean", "kind": "intersect", "target": "t", "tool": "b"}));
    let found = refusal(&ring, "j", "1/2", [face("b", "start"), wall("t")]);
    assert_eq!(found.code, "fillet.unsupported-contact");
}
