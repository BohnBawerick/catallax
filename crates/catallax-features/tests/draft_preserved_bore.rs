//! Planar draft around separated axial bores, checked in both kernels.
//!
//! Every expected value comes from the part's definition. No expectation reads
//! a kernel result, and neither kernel's answer supplies the other's.
use catallax_boolean::{Containment, classify_point, face_contains};
use catallax_features::{BuildError, Built, ConsumedFace, draft::bored, evaluate};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Integer, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{ElementId, HolePart, Origin, Recipe, Role, StepId};
use catallax_shapes::{Curve, Cylinder, Plane, Point3, RPoint3, Surface, Vector3};
use catallax_topology::{BodyKind, Loop, Sense};
use serde_json::{Value as Json, json};
use std::collections::{BTreeMap, BTreeSet};

/// One ordinary application request's allowance.
const WORK: u64 = 20_000_000;

fn q(x: &str) -> Rational {
    Rational::parse(x).unwrap()
}
fn fixture(name: &str) -> Json {
    let path = format!(
        "{}/tests/fixtures/draft-preserved-bore/{name}.json",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}
fn recipe(source: &Json) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}

// Outward rational interval arithmetic for the hand formulas.
#[derive(Clone, Debug)]
struct Bound {
    lo: Rational,
    hi: Rational,
}
impl Bound {
    fn new(lo: Rational, hi: Rational) -> Self {
        let grid = q("1e70");
        let round = |value: Integer| Rational::new(value, Integer::one()).unwrap() / grid.clone();
        Self {
            lo: round((&lo * &grid).floor()),
            hi: round((&hi * &grid).ceil()),
        }
    }
    fn exact(value: &str) -> Self {
        Self::new(q(value), q(value))
    }
    fn add(&self, other: &Self) -> Self {
        Self::new(&self.lo + &other.lo, &self.hi + &other.hi)
    }
    fn sub(&self, other: &Self) -> Self {
        Self::new(&self.lo - &other.hi, &self.hi - &other.lo)
    }
    fn mul(&self, other: &Self) -> Self {
        let values = [
            &self.lo * &other.lo,
            &self.lo * &other.hi,
            &self.hi * &other.lo,
            &self.hi * &other.hi,
        ];
        let lo = values.iter().min().unwrap().clone();
        let hi = values.iter().max().unwrap().clone();
        Self::new(lo, hi)
    }
    fn times(&self, value: &str) -> Self {
        self.mul(&Self::exact(value))
    }
    fn inverse(&self) -> Self {
        assert!(self.lo > Rational::zero());
        Self::new(self.hi.recip().unwrap(), self.lo.recip().unwrap())
    }
    fn square(&self) -> Self {
        self.mul(self)
    }
}

fn arctangent(x: &Rational, terms: u32) -> Bound {
    let square = x * x;
    let (mut sum, mut power) = (Rational::zero(), x.clone());
    for k in 0..terms {
        let term = &power / &Rational::from_i64(i64::from(2 * k + 1));
        sum = if k % 2 == 0 { sum + term } else { sum - term };
        power = &power * &square;
    }
    let omitted = &power / &Rational::from_i64(i64::from(2 * terms + 1));
    Bound::new(&sum - &omitted, &sum + &omitted)
}
/// Machin's identity, bounded by the first omitted alternating terms.
fn pi() -> Bound {
    arctangent(&q("1/5"), 80)
        .times("16")
        .sub(&arctangent(&q("1/239"), 25).times("4"))
}
/// Alternating series and its first omitted term, with power `first` of x.
fn series(x: &Rational, first: u32) -> [Rational; 2] {
    let square = x * x;
    let (mut sum, mut term) = (
        Rational::zero(),
        if first == 0 {
            Rational::one()
        } else {
            x.clone()
        },
    );
    for k in 0..40u32 {
        sum = if k % 2 == 0 {
            sum + term.clone()
        } else {
            sum - term.clone()
        };
        let n = 2 * k + first;
        term = &(&term * &square) / &Rational::from_i64(i64::from((n + 1) * (n + 2)));
    }
    [sum, term]
}
/// Tangent and secant of a small positive angle in degrees.
fn slope(degrees: &str) -> [Bound; 2] {
    let x = pi().times(degrees).times("1/180");
    let [low, low_omitted] = series(&x.lo, 1);
    let [high, high_omitted] = series(&x.hi, 1);
    let sine = Bound::new(&low - &low_omitted, &high + &high_omitted);
    let [low, low_omitted] = series(&x.hi, 0);
    let [high, high_omitted] = series(&x.lo, 0);
    let cosine = Bound::new(&low - &low_omitted, &high + &high_omitted);
    let secant = cosine.inverse();
    [sine.mul(&secant), secant]
}

struct Reference {
    volume: Bound,
    area: Bound,
    centroid: [Bound; 3],
    inertia: [[Bound; 3]; 3],
}
impl Reference {
    /// Section area `30(40-2zt)` for z in [0,20] and its moments, less a full
    /// cylinder of each `(radius, x)` bore on the local x axis.
    fn local(degrees: &str, bores: &[(&str, &str)]) -> Self {
        let pi = pi();
        let [t, secant] = slope(degrees);
        let t2 = t.square();
        let t3 = t2.mul(&t);
        let sum = |term: &dyn Fn(&Rational, &Rational) -> Rational| {
            let total = bores
                .iter()
                .fold(Rational::zero(), |sum, (r, x)| sum + term(&q(r), &q(x)));
            pi.mul(&Bound::new(total.clone(), total))
        };
        let r1 = sum(&|r, _| r.clone());
        let r2 = sum(&|r, _| r * r);
        let r4 = sum(&|r, _| &(r * r) * &(r * r));
        let r2x = sum(&|r, x| &(r * r) * x);
        let r2x2 = sum(&|r, x| &(&(r * r) * x) * x);
        let volume = Bound::exact("24000")
            .sub(&t.times("12000"))
            .sub(&r2.times("20"));
        let area = Bound::exact("4000")
            .sub(&t.times("2000"))
            .add(&secant.times("1200"))
            .add(&r1.times("40"))
            .sub(&r2.times("2"));
        let x = r2x.times("-20");
        let z = Bound::exact("240000")
            .sub(&t.times("160000"))
            .sub(&r2.times("200"));
        let x2 = Bound::exact("3200000")
            .sub(&t.times("4800000"))
            .add(&t2.times("3200000"))
            .sub(&t3.times("800000"))
            .sub(&r2x2.times("20"))
            .sub(&r4.times("5"));
        let y2 = Bound::exact("1800000")
            .sub(&t.times("900000"))
            .sub(&r4.times("5"));
        let z2 = Bound::exact("3200000")
            .sub(&t.times("2400000"))
            .sub(&r2.times("8000/3"));
        let xz = r2x.times("-200");
        let per_volume = volume.inverse();
        let xx = x.square().mul(&per_volume);
        let zz = z.square().mul(&per_volume);
        let cross = xz.times("-1").add(&x.mul(&z).mul(&per_volume));
        let zero = Bound::exact("0");
        Self {
            centroid: [x.mul(&per_volume), zero.clone(), z.mul(&per_volume)],
            volume,
            area,
            inertia: [
                [y2.add(&z2).sub(&zz), zero.clone(), cross.clone()],
                [zero.clone(), x2.add(&z2).sub(&xx.add(&zz)), zero.clone()],
                [cross, zero.clone(), x2.add(&y2).sub(&xx)],
            ],
        }
    }
    fn scaled(&self, factor: &str) -> Self {
        let cube = |x: &Bound| x.times(factor).times(factor).times(factor);
        Self {
            volume: cube(&self.volume),
            area: self.area.times(factor).times(factor),
            centroid: self.centroid.clone().map(|x| x.times(factor)),
            inertia: self
                .inertia
                .clone()
                .map(|row| row.map(|x| cube(&x).times(factor).times(factor))),
        }
    }
    fn placed(&self, frame: &Frame) -> Self {
        let m = frame
            .rotation
            .clone()
            .map(|row| row.map(|x| Bound::new(x.clone(), x)));
        let centroid = std::array::from_fn(|i| {
            (0..3).fold(
                Bound::new(frame.translation[i].clone(), frame.translation[i].clone()),
                |sum, j| sum.add(&m[i][j].mul(&self.centroid[j])),
            )
        });
        let inertia = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Bound::exact("0"), |sum, a| {
                    (0..3).fold(sum, |sum, b| {
                        sum.add(&m[i][a].mul(&self.inertia[a][b]).mul(&m[j][b]))
                    })
                })
            })
        });
        Self {
            volume: self.volume.clone(),
            area: self.area.clone(),
            centroid,
            inertia,
        }
    }
}

#[derive(Clone)]
struct Frame {
    rotation: [[Rational; 3]; 3],
    translation: [Rational; 3],
}

/// The exact proper rotation from quaternion (1,2,3,4), and its translation.
fn oblique() -> Frame {
    Frame {
        rotation: [
            [q("-2/3"), q("2/15"), q("11/15")],
            [q("2/3"), q("-1/3"), q("2/3")],
            [q("1/3"), q("14/15"), q("2/15")],
        ],
        translation: [q("7"), q("11"), q("-13")],
    }
}

/// One expected bore in local coordinates on the x axis.
struct Bored {
    wall: Origin,
    radius: &'static str,
    x: &'static str,
    /// A point of this wall only.
    on_wall: [&'static str; 3],
}

struct Case {
    name: &'static str,
    degrees: &'static str,
    scale: &'static str,
    frame: Option<Frame>,
    bores: Vec<Bored>,
    /// A point of the retained end cap outside every hole.
    on_cap: [&'static str; 3],
    material: &'static [Row],
}

fn side(step: &str, segment: &str) -> Origin {
    Origin {
        step: StepId::new(step).unwrap(),
        role: Role::Side {
            segment: ElementId::new(segment).unwrap(),
            path: None,
        },
    }
}
fn hole_wall(step: &str, at: &str) -> Origin {
    Origin {
        step: StepId::new(step).unwrap(),
        role: Role::Hole {
            at: ElementId::new(at).unwrap(),
            part: HolePart::Wall,
        },
    }
}
fn one_bore(wall: Origin) -> Vec<Bored> {
    vec![Bored {
        wall,
        radius: "5",
        x: "0",
        on_wall: ["5", "0", "10"],
    }]
}
fn two_bores(left: Origin, right: Origin) -> Vec<Bored> {
    vec![
        Bored {
            wall: left,
            radius: "3",
            x: "-8",
            on_wall: ["-5", "0", "10"],
        },
        Bored {
            wall: right,
            radius: "4",
            x: "8",
            on_wall: ["4", "0", "10"],
        },
    ]
}
fn native_two_bores() -> Vec<Bored> {
    two_bores(
        hole_wall("hole-left", "left"),
        hole_wall("hole-right", "right"),
    )
}
fn tool_two_bores() -> Vec<Bored> {
    two_bores(side("tool-left", "wall"), side("tool-right", "wall"))
}
fn one_bore_cases() -> Vec<Case> {
    let tool = || one_bore(side("tool", "wall"));
    let case = |name, degrees, scale, frame, bores| Case {
        name,
        degrees,
        scale,
        frame,
        bores,
        on_cap: ["10", "0", "20"],
        material: &MATERIAL,
    };
    vec![
        case("K-H", "3", "1", None, one_bore(hole_wall("hole", "centre"))),
        case("K-E", "5", "1", None, tool()),
        case("K-R", "3", "1", Some(oblique()), tool()),
        case("K-S", "5", "1/2", None, tool()),
    ]
}
fn two_bore_cases() -> Vec<Case> {
    let case = |name, degrees, scale, frame, bores| Case {
        name,
        degrees,
        scale,
        frame,
        bores,
        on_cap: ["0", "0", "20"],
        material: &TWO_BORE_MATERIAL,
    };
    vec![
        case("L-N", "3", "1", None, native_two_bores()),
        case("L-B", "5", "1", None, tool_two_bores()),
        case("L-R", "3", "1", Some(oblique()), native_two_bores()),
        case("L-S", "5", "1/2", None, tool_two_bores()),
    ]
}
fn cases() -> Vec<Case> {
    let mut all = one_bore_cases();
    all.extend(two_bore_cases());
    all
}

impl Case {
    fn point(&self, p: [&str; 3]) -> RPoint3 {
        let local = p.map(|x| &q(x) * &q(self.scale));
        let [x, y, z] = match &self.frame {
            None => local,
            Some(frame) => std::array::from_fn(|i| {
                (0..3).fold(frame.translation[i].clone(), |sum, j| {
                    sum + &frame.rotation[i][j] * &local[j]
                })
            }),
        };
        RPoint3::new(x, y, z)
    }
    fn local(&self, p: [f64; 3], translate: bool) -> [f64; 3] {
        let scale = q(self.scale).to_f64_nearest();
        let Some(frame) = &self.frame else {
            return if translate { p.map(|x| x / scale) } else { p };
        };
        let t = frame.translation.each_ref().map(Rational::to_f64_nearest);
        let d: [f64; 3] = std::array::from_fn(|i| p[i] - if translate { t[i] } else { 0.0 });
        std::array::from_fn(|i| {
            let along: f64 = (0..3)
                .map(|j| frame.rotation[j][i].to_f64_nearest() * d[j])
                .sum();
            if translate { along / scale } else { along }
        })
    }
    fn reference(&self) -> Reference {
        let bores: Vec<_> = self.bores.iter().map(|b| (b.radius, b.x)).collect();
        let reference = Reference::local(self.degrees, &bores).scaled(self.scale);
        match &self.frame {
            Some(frame) => reference.placed(frame),
            None => reference,
        }
    }
}

fn certified(x: &Real) -> f64 {
    Decider::default()
        .enclose(x, &q("1e-13"))
        .unwrap()
        .mid()
        .to_f64_nearest()
}
fn xyz(p: &Point3) -> [f64; 3] {
    [&p.x, &p.y, &p.z].map(certified)
}
fn vector(v: &Vector3) -> [f64; 3] {
    [&v.x, &v.y, &v.z].map(certified)
}
fn close(a: [f64; 3], b: [f64; 3], what: &str) {
    for i in 0..3 {
        assert!((a[i] - b[i]).abs() < 1e-9, "{what}: {a:?} vs {b:?}");
    }
}

/// Exact supports, senses, loops, adjacency and names of the admitted topology.
fn check_topology(case: &Case, built: &Built) {
    let model = &built.model;
    assert!(model.check(built.result).is_empty(), "{}", case.name);
    let counts = model.counts(built.result);
    let holes = case.bores.len();
    assert_eq!(
        (counts.lumps, counts.shells, counts.faces, counts.loops),
        (1, 1, 6 + holes, 6 + 4 * holes)
    );
    assert_eq!(
        (counts.edges, counts.ring_edges, counts.vertices),
        (12 + 2 * holes, 2 * holes, 8)
    );
    assert!(built.consumed.is_empty());
    let t = (q(case.degrees).to_f64_nearest() * std::f64::consts::PI / 180.0).tan();
    let base = |role: Role| Origin {
        step: StepId::new("base").unwrap(),
        role,
    };
    let planes = [
        (base(Role::Start), [0.0, 0.0, -1.0], 0.0),
        (base(Role::End), [0.0, 0.0, 1.0], 20.0),
        (side("base", "bottom"), [0.0, -1.0, 0.0], 15.0),
        (side("base", "top"), [0.0, 1.0, 0.0], 15.0),
        (side("base", "right"), [1.0, 0.0, t], 20.0),
        (side("base", "left"), [-1.0, 0.0, t], 20.0),
    ];
    let mut wall_rings = Vec::new();
    let mut cap_rings = Vec::new();
    let mut disks = 0;
    let mut seen = 0;
    for face in model.faces_of_body(built.result) {
        let f = model.face(face);
        let origins = model.face_origins(face);
        assert_eq!(origins.len(), 1, "{origins:?}");
        let loops: Vec<_> = f
            .loops
            .iter()
            .map(|lp| match model.loop_(*lp) {
                Loop::Coedges(cs) => cs.clone(),
                Loop::Vertex(_) => panic!("no vertex loops"),
            })
            .collect();
        seen += 1;
        match model.geom().surface(f.surface) {
            Surface::Plane(plane) => {
                let (_, normal, offset) = planes
                    .iter()
                    .find(|(origin, ..)| *origin == origins[0])
                    .unwrap_or_else(|| panic!("{origins:?}"));
                let length = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
                let sign = if f.sense == Sense::Same { 1.0 } else { -1.0 };
                let n = vector(&plane.normal);
                let n = case.local(n, false).map(|x| x * sign);
                let n_length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                close(
                    n.map(|x| x / n_length),
                    normal.map(|x| x / length),
                    "support normal",
                );
                let p = case.local(xyz(&plane.origin), true);
                let value: f64 = (0..3).map(|i| p[i] * normal[i]).sum();
                assert!(
                    (value - offset).abs() < 1e-9,
                    "{origins:?} {value} {offset}"
                );
                if origins[0].step.as_str() == "base"
                    && matches!(origins[0].role, Role::Start | Role::End)
                {
                    let [outer, rings @ ..] = loops.as_slice() else {
                        panic!("a cap keeps its outer rectangle and holes");
                    };
                    assert_eq!(outer.len(), 4);
                    assert_eq!(rings.len(), holes, "one exclusion per bore");
                    for ring in rings {
                        let [ring] = ring.as_slice() else {
                            panic!("one circular hole edge");
                        };
                        cap_rings.push(model.coedge(*ring).edge);
                    }
                } else {
                    assert_eq!(loops.len(), 1);
                    assert_eq!(loops[0].len(), 4);
                    disks += 1;
                }
            }
            Surface::Cylinder(cylinder) => {
                let bore = case
                    .bores
                    .iter()
                    .find(|bore| bore.wall == origins[0])
                    .unwrap_or_else(|| panic!("{origins:?}"));
                let x = q(bore.x).to_f64_nearest();
                assert_eq!(f.sense, Sense::Reversed);
                let scaled = &q(bore.radius) * &q(case.scale);
                let radius = &scaled * &scaled;
                assert_eq!(cylinder.radius_sq.rational_ref(), Some(&radius));
                let axis = case.local(vector(&cylinder.axis), false);
                assert!(axis[0].abs() < 1e-9 && axis[1].abs() < 1e-9 && axis[2] != 0.0);
                let origin = case.local(xyz(&cylinder.origin), true);
                assert!((origin[0] - x).abs() < 1e-9 && origin[1].abs() < 1e-9);
                let mut heights = Vec::new();
                for cycle in &loops {
                    let [coedge] = cycle.as_slice() else {
                        panic!("one complete ring per loop");
                    };
                    let edge = model.coedge(*coedge).edge;
                    wall_rings.push(edge);
                    let e = model.edge(edge);
                    assert!(e.start.is_none() && e.end.is_none());
                    let Curve::Circle(circle) = model.geom().curve(e.curve) else {
                        panic!("rings stay circles");
                    };
                    assert_eq!(circle.radius_sq.rational_ref(), Some(&radius));
                    let centre = case.local(xyz(&circle.center), true);
                    assert!((centre[0] - x).abs() < 1e-9 && centre[1].abs() < 1e-9);
                    heights.push((centre[2] * 1e6).round() / 1e6);
                }
                heights.sort_by(f64::total_cmp);
                assert_eq!(heights, [0.0, 20.0]);
            }
            surface => panic!("unexpected {surface:?}"),
        }
    }
    assert_eq!(seen, 6 + holes);
    // Each bore wall is an annulus and each cap has Euler characteristic 1-n, with
    // ring edges contributing nothing. The four side disks close the line graph,
    // so the shell's characteristic is 2-2n: genus n, not a filled body.
    assert_eq!(disks, 4);
    let n = holes as i64;
    assert_eq!(
        counts.vertices as i64 - (counts.edges - counts.ring_edges) as i64 + disks + 2 * (1 - n),
        2 - 2 * n
    );
    wall_rings.sort();
    cap_rings.sort();
    assert_eq!(wall_rings, cap_rings);
}

type Row = ([&'static str; 3], Location, Location);
const MATERIAL: [Row; 7] = [
    (["10", "0", "10"], Location::Inside, Location::Inside),
    (["0", "0", "10"], Location::Outside, Location::Outside),
    (["19", "0", "19"], Location::Inside, Location::Outside),
    (["21", "0", "10"], Location::Outside, Location::Outside),
    (["5", "0", "10"], Location::Boundary, Location::Boundary),
    (["0", "0", "20"], Location::Outside, Location::Outside),
    (["10", "0", "20"], Location::Boundary, Location::Boundary),
];
const TWO_BORE_MATERIAL: [Row; 10] = [
    (["0", "0", "10"], Location::Inside, Location::Inside),
    (["-8", "0", "10"], Location::Outside, Location::Outside),
    (["8", "0", "10"], Location::Outside, Location::Outside),
    (["19", "0", "19"], Location::Inside, Location::Outside),
    (["21", "0", "10"], Location::Outside, Location::Outside),
    (["-8", "0", "20"], Location::Outside, Location::Outside),
    (["8", "0", "20"], Location::Outside, Location::Outside),
    (["-5", "0", "10"], Location::Boundary, Location::Boundary),
    (["4", "0", "10"], Location::Boundary, Location::Boundary),
    (["0", "0", "20"], Location::Boundary, Location::Boundary),
];

#[test]
fn both_kernels_draft_the_walls_around_an_unchanged_bore() {
    for case in cases() {
        let source = recipe(&fixture(case.name));
        let (built, _) = with_operation_budget(WORK, || evaluate(&source, &mut Decider::default()));
        let mut built = built.unwrap();
        check_topology(&case, &built);
        let (plan, _) = with_operation_budget(WORK, || {
            Plan::compile_checked(&source, &mut Decider::default())
        });
        let plan = plan.unwrap();
        assert!(plan.consumed_faces().is_empty());
        let mut storage = Workspace::default();
        for (point, three, five) in case.material.iter().cloned() {
            let expected = if case.degrees == "3" { three } else { five };
            let at = case.point(point);
            let (field, work) = with_operation_budget(WORK, || {
                plan.proven_checked(&at, &mut storage, &mut Decider::default())
            });
            let field = field.unwrap_or_else(|e| panic!("{} {point:?} {e:?} {work:?}", case.name));
            assert_eq!(field.location, expected, "{} field {point:?}", case.name);
            let (native, _) = with_operation_budget(WORK, || {
                classify_point(&mut built.model, built.result, &at, &mut Decider::default())
            });
            assert!(
                match expected {
                    Location::Inside => matches!(native, Outcome::Proven(Containment::Inside)),
                    Location::Outside => matches!(native, Outcome::Proven(Containment::Outside)),
                    _ => matches!(native, Outcome::ProvenCoincident(_)),
                },
                "{} native {point:?}: {native:?}",
                case.name
            );
        }
        // Names resolve to retained regions, and the cap's support plane is not
        // its face: the hole stays excluded from the named cap.
        let end = Origin {
            step: StepId::new("base").unwrap(),
            role: Role::End,
        };
        let state = |origin: &Origin, point: [&str; 3], storage: &mut Workspace| {
            plan.face_query(plan.select(origin).unwrap(), &case.point(point), storage)
                .unwrap()
                .state
        };
        assert_eq!(state(&end, case.on_cap, &mut storage), FaceState::OnFace);
        let cap = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .find(|face| built.model.face_origins(*face).contains(&end))
            .unwrap();
        let on_face = |point: &RPoint3, built: &mut Built| {
            face_contains(&mut built.model, cap, point, &mut Decider::default()).unwrap()
        };
        assert!(
            on_face(&case.point(case.on_cap), &mut built),
            "{}",
            case.name
        );
        for bore in &case.bores {
            assert_eq!(
                state(&bore.wall, bore.on_wall, &mut storage),
                FaceState::OnFace,
                "{}",
                case.name
            );
            let centre = [bore.x, "0", "20"];
            assert_eq!(state(&end, centre, &mut storage), FaceState::NotOnFace);
            assert!(!on_face(&case.point(centre), &mut built), "{}", case.name);
            // A surviving label on the other bore names the wrong wall.
            for other in &case.bores {
                if other.wall != bore.wall {
                    assert_eq!(
                        state(&bore.wall, other.on_wall, &mut storage),
                        FaceState::NotOnFace
                    );
                }
            }
        }
    }
}

fn overlap(e: &Enclosure, reference: &Bound, what: &str) {
    assert!(
        e.lo <= reference.hi && reference.lo <= e.hi,
        "{what}: {e:?} misses {reference:?}"
    );
}

#[test]
fn requested_properties_overlap_independent_enclosures() {
    for case in cases() {
        let built = evaluate(&recipe(&fixture(case.name)), &mut Decider::default()).unwrap();
        let (answer, _) = with_operation_budget(WORK, || {
            catallax_measure::mass_properties(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
        });
        let answer = answer.unwrap();
        assert!(answer.marks.is_empty(), "{}", case.name);
        let properties = answer.value;
        let reference = case.reference();
        overlap(&properties.volume, &reference.volume, "volume");
        overlap(&properties.area, &reference.area, "area");
        assert!(properties.volume.width() <= q("1e-4"), "{}", case.name);
        assert!(properties.area.width() <= q("1e-5"), "{}", case.name);
        let centroid = properties.centroid.unwrap();
        let inertia = properties.inertia.unwrap();
        for (i, row) in inertia.iter().enumerate() {
            overlap(&centroid[i], &reference.centroid[i], "centroid");
            for (j, entry) in row.iter().enumerate() {
                overlap(entry, &reference.inertia[i][j], "inertia");
            }
        }
    }
}

#[test]
fn bounds_and_published_table_rows_match_the_definition() {
    // The rounded report rows must agree with these formulas to their digits.
    const ONE: [(&str, &str); 1] = [("5", "0")];
    const TWO: [(&str, &str); 2] = [("3", "-8"), ("4", "8")];
    for (degrees, volume, area, height) in [
        (
            "3",
            "21800.310321808609",
            "5568.070154669892",
            "9.903840305923",
        ),
        (
            "5",
            "21379.339710894015",
            "5500.845376038638",
            "9.836311757596",
        ),
    ] {
        let reference = Reference::local(degrees, &ONE);
        for (bound, row) in [
            (&reference.volume, volume),
            (&reference.area, area),
            (&reference.centroid[2], height),
        ] {
            agree(bound, row);
        }
    }
    for (degrees, volume, area, x, z, ixx, iyy, izz, ixz) in [
        (
            "3",
            "21800.310321808609",
            "5819.397566957075",
            "-0.161400627793",
            "9.903840305923",
            "2474014.845109994",
            "3577199.526884613",
            "4598263.511822547",
            "338.345939101",
        ),
        (
            "5",
            "21379.339710894015",
            "5752.172788325821",
            "-0.164578692308",
            "9.836311757596",
            "2428038.441612656",
            "3409680.527797197",
            "4413575.324595280",
            "575.950793394",
        ),
    ] {
        let reference = Reference::local(degrees, &TWO);
        for (bound, row) in [
            (&reference.volume, volume),
            (&reference.area, area),
            (&reference.centroid[0], x),
            (&reference.centroid[2], z),
            (&reference.inertia[0][0], ixx),
            (&reference.inertia[1][1], iyy),
            (&reference.inertia[2][2], izz),
            (&reference.inertia[0][2], ixz),
            (&reference.inertia[2][0], ixz),
        ] {
            agree(bound, row);
        }
    }
    // Bottom edges retain every extremum, so drafting keeps the local box.
    for name in ["K-H", "L-N"] {
        let built = evaluate(&recipe(&fixture(name)), &mut Decider::default()).unwrap();
        let boxed = catallax_measure::bounding_box(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value
        .unwrap();
        for (axis, (lo, hi)) in [("-20", "20"), ("-15", "15"), ("0", "20")]
            .iter()
            .enumerate()
        {
            assert!((&boxed.min[axis] - &q(lo)).abs() <= q("0.001"));
            assert!((&boxed.max[axis] - &q(hi)).abs() <= q("0.001"));
        }
    }
}

fn agree(bound: &Bound, row: &str) {
    let middle = (&bound.lo + &bound.hi).mul_pow2(-1);
    assert!((&middle - &q(row)).abs() < q("1e-8"), "{row}");
}

/// Closed, selectable display whose triangles stay out of every `(radius, x)` hole.
fn check_display(name: &str, options: &catallax_mesh::MeshOptions, bores: &[(f64, f64)]) {
    let source = recipe(&fixture(name));
    let built = evaluate(&source, &mut Decider::default()).unwrap();
    let mesh =
        catallax_mesh::mesh_body(&built.model, built.result, options, &mut Decider::default())
            .unwrap()
            .value;
    let faces = 6 + bores.len();
    assert_eq!(mesh.faces.len(), faces);
    assert!(mesh.unavailable_controls.is_empty());
    assert!(
        mesh.faces
            .iter()
            .all(|face| face.max_chord_error <= options.chord_tolerance)
    );
    assert!(
        mesh.edges
            .iter()
            .all(|edge| edge.max_chord_error <= options.chord_tolerance)
    );
    let chord = options.chord_tolerance.to_f64_nearest();
    let mut regions = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for triangle in &mesh.triangles {
        *regions.entry(triangle.face).or_insert(0) += 1;
        let face = built.model.face(triangle.face);
        if let Surface::Plane(plane) = built.model.geom().surface(face.surface)
            && plane.normal.x.rational_ref().is_some_and(Rational::is_zero)
            && plane.normal.y.rational_ref().is_some_and(Rational::is_zero)
        {
            let p = triangle.vertices.map(|i| mesh.positions[i]);
            let x = (p[0][0] + p[1][0] + p[2][0]) / 3.0;
            let y = (p[0][1] + p[1][1] + p[2][1]) / 3.0;
            for (radius, centre) in bores {
                assert!(
                    (x - centre).hypot(y) >= radius - chord,
                    "filled cap hole at {x} {y}"
                );
            }
        }
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
            counts[usize::from(a > b)] += 1;
        }
    }
    assert_eq!(regions.len(), faces);
    assert!(edges.values().all(|counts| *counts == [1, 1]));
}

#[test]
fn certified_display_is_closed_selectable_and_keeps_the_cap_holes_open() {
    let finer = catallax_mesh::MeshOptions {
        chord_tolerance: q("0.005"),
        ..Default::default()
    };
    check_display("K-H", &finer, &[(5.0, 0.0)]);
    check_display("L-N", &Default::default(), &[(3.0, -8.0), (4.0, 8.0)]);
}

fn refusal(source: &Json) -> (String, String) {
    let parsed = recipe(source);
    let Err(BuildError::Geometry { refusal, .. }) = evaluate(&parsed, &mut Decider::default())
    else {
        panic!("B-rep must refuse {source}");
    };
    let field = Plan::compile_checked(&parsed, &mut Decider::default())
        .err()
        .unwrap_or_else(|| panic!("field must refuse {source}"));
    (refusal.code.to_string(), field.code.to_string())
}

#[test]
fn contact_crossing_and_invalid_angle_keep_their_typed_refusals() {
    for name in ["K-T", "K-X"] {
        let (brep, field) = refusal(&fixture(name));
        assert_eq!(brep, "draft.unsupported-arrangement", "{name}");
        assert_eq!(field, "draft.unsupported-source-boundary", "{name}");
    }
    assert_eq!(
        refusal(&fixture("K-BAD")),
        ("draft.angle".into(), "draft.angle".into())
    );
}

#[test]
fn two_bores_keep_outer_wall_contact_crossing_and_angle_refusals() {
    for name in ["L-T", "L-X"] {
        let (brep, field) = refusal(&fixture(name));
        assert_eq!(brep, "draft.unsupported-arrangement", "{name}");
        assert_eq!(field, "draft.unsupported-source-boundary", "{name}");
    }
    assert_eq!(
        refusal(&fixture("L-BAD")),
        ("draft.angle".into(), "draft.angle".into())
    );
}

fn plane(point: [&str; 3], normal: [&str; 3]) -> Plane {
    let vector = |[x, y, z]: [&str; 3]| Vector3::new(q(x).into(), q(y).into(), q(z).into());
    let [x, y, z] = point.map(|x| q(x).into());
    Plane {
        origin: Point3::new(x, y, z),
        normal: vector(normal),
        ref_dir: vector(if normal[0] == "0" && normal[1] == "0" {
            ["1", "0", "0"]
        } else {
            ["0", "0", "1"]
        }),
    }
}
fn bore(x: &str, radius: &str, axis: &str) -> Cylinder {
    Cylinder {
        origin: Point3::new(q(x).into(), Real::from_i64(0), Real::from_i64(0)),
        axis: Vector3::new(Real::from_i64(0), Real::from_i64(0), q(axis).into()),
        radius_sq: (&q(radius) * &q(radius)).into(),
        ref_dir: Vector3::new(Real::from_i64(1), Real::from_i64(0), Real::from_i64(0)),
    }
}

#[test]
fn bore_pairs_are_admitted_only_at_a_strict_gap() {
    let envelope = [
        plane(["0", "0", "0"], ["0", "0", "-1"]),
        plane(["0", "0", "20"], ["0", "0", "1"]),
        plane(["-40", "0", "0"], ["-1", "0", "0"]),
        plane(["40", "0", "0"], ["1", "0", "0"]),
        plane(["0", "-40", "0"], ["0", "-1", "0"]),
        plane(["0", "40", "0"], ["0", "1", "0"]),
    ];
    // Radii 3 and 4 touch at 7, overlap at 13/2, and leave a 1e-7 mm gap above 7.
    for (separation, admitted) in [("7", false), ("13/2", false), ("70000001/10000000", true)] {
        for axis in ["1", "-1"] {
            let pair = [&bore("0", "3", "1"), &bore(separation, "4", axis)];
            let caps = bored::certify_bores(envelope.iter(), &pair, &mut Decider::default());
            match caps {
                Ok(caps) => assert!(admitted && caps == [0, 1], "{separation}"),
                Err(refusal) => {
                    assert!(!admitted, "{separation} {refusal:?}");
                    assert_eq!(refusal.code, "draft.unsupported-arrangement");
                }
            }
        }
    }
}

#[test]
fn the_new_branch_admits_only_a_separated_unselected_axial_through_bore() {
    let edited = |edit: &dyn Fn(&mut Json)| {
        let mut source = fixture("K-E");
        edit(&mut source);
        source
    };
    let draft = 4;
    for source in [
        // The bore itself is selected.
        edited(&|s| {
            s["steps"][draft]["faces"] =
                json!([{"face":{"step":"tool","role":{"side":{"segment":"wall"}}}}]);
        }),
        // A cap carrying a ring is selected.
        edited(&|s| {
            s["steps"][draft]["faces"] = json!([{"face":{"step":"base","role":"end"}}]);
        }),
        // The tool stops inside the part.
        edited(&|s| {
            s["steps"][2]["height"] = json!("15");
            s["steps"][2]["pose"]["translation"][2] = json!("10");
        }),
        // Pull and neutral each leave the bore axis.
        edited(&|s| s["steps"][draft]["pull"] = json!(["0", "1", "20"])),
        edited(&|s| s["steps"][draft]["neutral"]["normal"] = json!(["0", "1", "20"])),
    ] {
        assert_eq!(
            refusal(&source),
            (
                "draft.unsupported-arrangement".into(),
                "draft.unsupported-source-boundary".into()
            ),
            "{source}"
        );
    }
    // Saved intent inside the bore names no material, on the wall neither.
    for point in [["0", "0", "10"], ["5", "0", "10"]] {
        let source = edited(&|s| s["steps"][draft]["keep_regions"] = json!([point]));
        assert_eq!(
            refusal(&source),
            (
                "draft.invalid-keep-region".into(),
                "draft.invalid-keep-region".into()
            )
        );
    }
    let kept = recipe(&edited(&|s| {
        s["steps"][draft]["keep_regions"] = json!([["10", "0", "10"]]);
    }));
    let built = evaluate(&kept, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 7);
    Plan::compile_checked(&kept, &mut Decider::default()).unwrap();
}

#[test]
fn every_bore_of_a_collection_stays_unselected_axial_and_through() {
    let edited = |edit: &dyn Fn(&mut Json)| {
        let mut source = fixture("L-B");
        edit(&mut source);
        source
    };
    let draft = 6;
    for source in [
        // The second bore is selected.
        edited(&|s| {
            s["steps"][draft]["faces"] =
                json!([{"face":{"step":"tool-right","role":{"side":{"segment":"wall"}}}}]);
        }),
        // A cap carrying both rings is selected.
        edited(&|s| {
            s["steps"][draft]["faces"] = json!([{"face":{"step":"base","role":"end"}}]);
        }),
        // The second tool stops inside the part.
        edited(&|s| {
            s["steps"][4]["height"] = json!("15");
            s["steps"][4]["pose"]["translation"][2] = json!("10");
        }),
    ] {
        assert_eq!(
            refusal(&source),
            (
                "draft.unsupported-arrangement".into(),
                "draft.unsupported-source-boundary".into()
            ),
            "{source}"
        );
    }
    // Saved intent inside either bore or on its wall names no material.
    for point in [["-8", "0", "10"], ["8", "0", "10"], ["4", "0", "10"]] {
        let source = edited(&|s| s["steps"][draft]["keep_regions"] = json!([point]));
        assert_eq!(
            refusal(&source),
            (
                "draft.invalid-keep-region".into(),
                "draft.invalid-keep-region".into()
            ),
            "{point:?}"
        );
    }
    let kept = recipe(&edited(&|s| {
        s["steps"][draft]["keep_regions"] = json!([["0", "0", "10"]]);
    }));
    let built = evaluate(&kept, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 8);
    Plan::compile_checked(&kept, &mut Decider::default()).unwrap();
}

#[test]
fn one_tool_step_may_cut_several_bores() {
    // Equal radius 3 bores at x=-8 and x=8, from one native step with two
    // placements, one Boolean step with a union tool and one two-circle cut.
    let mut placed = fixture("L-N");
    let steps = placed["steps"].as_array_mut().unwrap();
    steps[2]["at"] = json!([
        {"id":"left","point":["-8","0","20"],"direction":["0","0","-1"]},
        {"id":"right","point":["8","0","20"],"direction":["0","0","-1"]},
    ]);
    steps.remove(3);
    steps[3]["body"] = json!("hole-left");
    let mut union = fixture("L-B");
    let steps = union["steps"].as_array_mut().unwrap();
    steps[3] = steps[4].clone();
    steps[3]["radius"] = json!("3");
    steps[4] = json!({"id":"tools","op":"boolean","kind":"union","target":"tool-left",
                      "tool":"tool-right"});
    steps[5] = json!({"id":"cut","op":"boolean","kind":"subtract","target":"base",
                      "tool":"tools"});
    steps[6]["body"] = json!("cut");
    let mut cut = fixture("L-N");
    let steps = cut["steps"].as_array_mut().unwrap();
    steps[2] = json!({"id":"drill","op":"sketch","plane":steps[0]["plane"],"loops":[
        [{"kind":"circle","id":"left","centre":["-8","0"],"radius":"3"}],
        [{"kind":"circle","id":"right","centre":["8","0"],"radius":"3"}]]});
    steps[3] = json!({"id":"cut","op":"extrude_through","target":"base","profile":"drill"});
    steps[4]["body"] = json!("cut");
    for source in [placed, union, cut] {
        let source = recipe(&source);
        let mut built = evaluate(&source, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        assert_eq!(built.model.counts(built.result).faces, 8);
        let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
        let mut storage = Workspace::default();
        for (point, expected) in [
            (["0", "0", "10"], Location::Inside),
            (["-8", "0", "10"], Location::Outside),
            (["8", "0", "10"], Location::Outside),
            (["-5", "0", "10"], Location::Boundary),
            (["5", "0", "10"], Location::Boundary),
        ] {
            let at = RPoint3::new(q(point[0]), q(point[1]), q(point[2]));
            let field = plan
                .proven_checked(&at, &mut storage, &mut Decider::default())
                .unwrap();
            assert_eq!(field.location, expected, "field {point:?}");
            let native =
                classify_point(&mut built.model, built.result, &at, &mut Decider::default());
            assert!(
                match expected {
                    Location::Inside => matches!(native, Outcome::Proven(Containment::Inside)),
                    Location::Outside => matches!(native, Outcome::Proven(Containment::Outside)),
                    _ => matches!(native, Outcome::ProvenCoincident(_)),
                },
                "native {point:?}: {native:?}"
            );
        }
    }
}

/// Comparable `(source, by)` removal history, independent of list order.
fn history(consumed: &[ConsumedFace]) -> BTreeSet<(BTreeSet<Origin>, BTreeSet<Origin>)> {
    consumed
        .iter()
        .map(|face| {
            (
                face.source.iter().cloned().collect(),
                face.by.iter().cloned().collect(),
            )
        })
        .collect()
}

#[test]
fn an_empty_result_reports_the_same_removed_faces_in_both_kernels() {
    for (name, draft, walls) in [
        ("K-E", 4, vec![side("tool", "wall")]),
        ("K-H", 3, vec![hole_wall("hole", "centre")]),
        (
            "L-B",
            6,
            vec![side("tool-left", "wall"), side("tool-right", "wall")],
        ),
        (
            "L-N",
            4,
            vec![
                hole_wall("hole-left", "left"),
                hole_wall("hole-right", "right"),
            ],
        ),
    ] {
        let neutral: &dyn Fn(&mut Json) =
            &|s| s["steps"][draft]["neutral"]["point"][2] = json!("100");
        let keep: &dyn Fn(&mut Json) = &|s| s["steps"][draft]["keep_regions"] = json!([]);
        for edit in [neutral, keep] {
            let mut source = fixture(name);
            edit(&mut source);
            let source = recipe(&source);
            let built = evaluate(&source, &mut Decider::default()).unwrap();
            assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
            let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
            assert!(plan.is_proven_empty());
            let native = history(&built.consumed);
            assert_eq!(native.len(), 6 + walls.len(), "{name}");
            for wall in &walls {
                assert!(native.iter().any(|(removed, _)| removed.contains(wall)));
            }
            assert_eq!(history(plan.consumed_faces()), native, "{name}");
        }
    }
}
