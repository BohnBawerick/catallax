//! Ordinary-angle planar draft of a concave L bracket, checked in both kernels.
//!
//! Every expected value comes from the part's definition. No expectation reads
//! a kernel result, and neither kernel's answer supplies the other's.
use catallax_boolean::{Containment, classify_point, face_contains};
use catallax_features::{BuildError, Built, evaluate};
use catallax_field::{FaceState, Location, Plan, Value, Workspace};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Integer, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{EdgeName, ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::{Point3, RPoint3, Surface, Vector3};
use catallax_topology::{Loop, Sense};
use serde_json::{Value as Json, json};
use std::collections::BTreeMap;

/// One ordinary application request's allowance.
const WORK: u64 = 20_000_000;
/// The section's vertices, each with the segment that leaves it.
const SIDES: [(&str, [i32; 2]); 6] = [
    ("bottom", [-20, -15]),
    ("right", [20, -15]),
    ("notch", [20, -5]),
    ("inner", [0, -5]),
    ("top", [0, 15]),
    ("left", [-20, 15]),
];

fn q(x: &str) -> Rational {
    Rational::parse(x).unwrap()
}
fn fixture(name: &str) -> Json {
    let path = format!(
        "{}/tests/fixtures/concave-draft/{name}.json",
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
    fn nearest(&self) -> f64 {
        (&self.lo + &self.hi).mul_pow2(-1).to_f64_nearest()
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
    /// The two arms' rectangular moments, integrated by hand for z in [0,20].
    fn local(degrees: &str) -> Self {
        let [t, secant] = slope(degrees);
        let (t2, zero) = (t.square(), Bound::exact("0"));
        let t3 = t2.mul(&t);
        let volume = Bound::exact("16000").sub(&t.times("12000"));
        let area = Bound::exact("3200")
            .sub(&t.times("2000"))
            .add(&secant.times("1200"));
        let moment = [
            t.times("80000").sub(&Bound::exact("80000")),
            Bound::exact("-40000"),
            Bound::exact("160000").sub(&t.times("160000")),
        ];
        let square = [
            Bound::exact("6400000/3")
                .sub(&t.times("3200000"))
                .add(&t2.times("6400000/3"))
                .sub(&t3.times("800000")),
            Bound::exact("4000000/3").sub(&t.times("900000")),
            Bound::exact("6400000/3").sub(&t.times("2400000")),
        ];
        let product = |i: usize, j: usize| match i + j {
            1 => t.times("400000").sub(&Bound::exact("400000")),
            2 => t.times("3200000/3").sub(&Bound::exact("800000")),
            _ => Bound::exact("-400000"),
        };
        let per_volume = volume.inverse();
        let inertia = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    (0..3).filter(|k| *k != i).fold(zero.clone(), |sum, k| {
                        sum.add(&square[k])
                            .sub(&moment[k].square().mul(&per_volume))
                    })
                } else {
                    moment[i]
                        .mul(&moment[j])
                        .mul(&per_volume)
                        .sub(&product(i, j))
                }
            })
        });
        Self {
            centroid: moment.clone().map(|m| m.mul(&per_volume)),
            volume,
            area,
            inertia,
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

struct Case {
    name: &'static str,
    degrees: &'static str,
    scale: &'static str,
    frame: Option<Frame>,
}

fn cases() -> [Case; 3] {
    [
        Case {
            name: "R-P3",
            degrees: "3",
            scale: "1",
            frame: None,
        },
        Case {
            name: "R-P5",
            degrees: "5",
            scale: "1",
            frame: None,
        },
        Case {
            name: "R-PRS",
            degrees: "5",
            scale: "1/2",
            frame: Some(oblique()),
        },
    ]
}

fn origin(role: Role) -> Origin {
    Origin {
        step: StepId::new("base").unwrap(),
        role,
    }
}
fn side(segment: &str) -> Origin {
    origin(Role::Side {
        segment: ElementId::new(segment).unwrap(),
        path: None,
    })
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
        let reference = Reference::local(self.degrees).scaled(self.scale);
        match &self.frame {
            Some(frame) => reference.placed(frame),
            None => reference,
        }
    }
    /// The drafted corner of one section vertex: the three walls move in x only.
    fn corner(&self, index: usize, top: bool) -> [f64; 3] {
        let [x, y] = SIDES[index].1.map(f64::from);
        let moved = if !top {
            0.0
        } else {
            let run = 20.0 * slope(self.degrees)[0].nearest();
            match index {
                1..=4 => -run,
                _ => run,
            }
        };
        [x + moved, y, if top { 20.0 } else { 0.0 }]
    }
    fn built(&self) -> (Recipe, Built) {
        let source = recipe(&fixture(self.name));
        let (built, work) =
            with_operation_budget(WORK, || evaluate(&source, &mut Decider::default()));
        assert!(!work.exhausted, "{}: {work:?}", self.name);
        (
            source,
            built.unwrap_or_else(|e| panic!("{}: {e:?}", self.name)),
        )
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

/// Eight disk faces, eighteen lines and twelve vertices of a genus-zero shell.
fn check_topology(case: &Case, built: &Built) {
    let model = &built.model;
    assert!(model.check(built.result).is_empty(), "{}", case.name);
    let counts = model.counts(built.result);
    assert_eq!(
        (counts.lumps, counts.shells, counts.faces, counts.loops),
        (1, 1, 8, 8),
        "{}",
        case.name
    );
    assert_eq!(
        (counts.edges, counts.ring_edges, counts.vertices),
        (18, 0, 12),
        "{}",
        case.name
    );
    assert_eq!(counts.vertices as i64 - counts.edges as i64 + 8, 2);
    assert!(built.consumed.is_empty());
    let t = slope(case.degrees)[0].nearest();
    let planes = [
        (origin(Role::Start), [0.0, 0.0, -1.0], 0.0),
        (origin(Role::End), [0.0, 0.0, 1.0], 20.0),
        (side("bottom"), [0.0, -1.0, 0.0], 15.0),
        (side("right"), [1.0, 0.0, t], 20.0),
        (side("notch"), [0.0, 1.0, 0.0], -5.0),
        (side("inner"), [1.0, 0.0, t], 0.0),
        (side("top"), [0.0, 1.0, 0.0], 15.0),
        (side("left"), [-1.0, 0.0, t], 20.0),
    ];
    let mut seen = Vec::new();
    for face in model.faces_of_body(built.result) {
        let f = model.face(face);
        let origins = model.face_origins(face);
        assert_eq!(origins.len(), 1, "{origins:?}");
        seen.push(origins[0].clone());
        let [wire] = f.loops.as_slice() else {
            panic!("{}: {origins:?} keeps one wire", case.name);
        };
        let Loop::Coedges(wire) = model.loop_(*wire) else {
            panic!("no vertex loops");
        };
        let cap = matches!(origins[0].role, Role::Start | Role::End);
        assert_eq!(wire.len(), if cap { 6 } else { 4 }, "{origins:?}");
        let Surface::Plane(plane) = model.geom().surface(f.surface) else {
            panic!("every face stays planar");
        };
        let (_, normal, offset) = planes
            .iter()
            .find(|(origin, ..)| *origin == origins[0])
            .unwrap_or_else(|| panic!("{origins:?}"));
        let length = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
        let sign = if f.sense == Sense::Same { 1.0 } else { -1.0 };
        let n = case.local(vector(&plane.normal), false).map(|x| x * sign);
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
    }
    seen.sort_by_key(|origin| planes.iter().position(|(o, ..)| o == origin));
    assert_eq!(seen, planes.map(|(origin, ..)| origin));
}

/// Each semantic name binds its one current edge, with these drafted end points.
fn check_edge_names(case: &Case, built: &Built) {
    let model = &built.model;
    let ends = |name: EdgeName, expected: [[f64; 3]; 2]| {
        let [edge] = model.edges_with_name(built.result, &name)[..] else {
            panic!("{} resolves {name:?} once", case.name);
        };
        let edge = model.edge(edge);
        let found = [edge.start, edge.end].map(|vertex| {
            let point = model.vertex(vertex.expect("lines keep end points")).point;
            case.local(xyz(&model.geom().position(point)), true)
        });
        // An edge's stored direction is not part of its name.
        let swapped = (0..3).any(|i| (found[0][i] - expected[0][i]).abs() >= 1e-9);
        close(found[0], expected[usize::from(swapped)], "edge start");
        close(found[1], expected[usize::from(!swapped)], "edge end");
    };
    for (i, (segment, _)) in SIDES.iter().enumerate() {
        let next = (i + 1) % 6;
        for (cap, top) in [(Role::Start, false), (Role::End, true)] {
            ends(
                EdgeName::new(origin(cap), side(segment)).unwrap(),
                [case.corner(i, top), case.corner(next, top)],
            );
        }
        ends(
            EdgeName::new(side(segment), side(SIDES[next].0)).unwrap(),
            [case.corner(next, false), case.corner(next, true)],
        );
    }
}

type Row = ([&'static str; 3], Location, Location);
const MATERIAL: [Row; 8] = [
    (["-10", "-10", "10"], Location::Inside, Location::Inside),
    (["-10", "5", "10"], Location::Inside, Location::Inside),
    (["-10", "-5", "10"], Location::Inside, Location::Inside),
    (["5", "5", "10"], Location::Outside, Location::Outside),
    (["21", "-10", "10"], Location::Outside, Location::Outside),
    (["19", "-10", "19"], Location::Inside, Location::Outside),
    (["-1", "0", "19"], Location::Inside, Location::Outside),
    (["-10", "-10", "20"], Location::Boundary, Location::Boundary),
];

/// One strict or exact point, proved separately by each kernel.
fn check_location(
    case: &Case,
    point: [&str; 3],
    expected: Location,
    (plan, storage): (&Plan, &mut Workspace),
    built: &mut Built,
) {
    let at = case.point(point);
    let (field, work) = with_operation_budget(WORK, || {
        plan.proven_checked(&at, storage, &mut Decider::default())
    });
    let field = field.unwrap_or_else(|e| panic!("{} {point:?} {e:?} {work:?}", case.name));
    assert_eq!(field.location, expected, "{} field {point:?}", case.name);
    let (native, work) = with_operation_budget(WORK, || {
        classify_point(&mut built.model, built.result, &at, &mut Decider::default())
    });
    assert!(!work.exhausted, "{} {point:?}", case.name);
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

#[test]
fn both_kernels_draft_the_concave_bracket_within_the_application_allowance() {
    for case in cases() {
        let (source, mut built) = case.built();
        check_topology(&case, &built);
        check_edge_names(&case, &built);
        let (plan, work) = with_operation_budget(WORK, || {
            Plan::compile_checked(&source, &mut Decider::default())
        });
        assert!(!work.exhausted, "{}: {work:?}", case.name);
        let plan = plan.unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        assert!(plan.consumed_faces().is_empty());
        let mut storage = Workspace::default();
        for (point, three, five) in MATERIAL {
            let expected = if case.degrees == "3" { three } else { five };
            check_location(&case, point, expected, (&plan, &mut storage), &mut built);
        }
    }
}

/// These points tie between the two x walls and between the caps. The seam point's
/// foot also meets the inner wall's trim, and it lies on the notch's carrier, which
/// is not its face. Each stays (1 - t) cos(a) * 10 from those walls, at scale.
#[test]
fn a_tied_distance_is_the_wall_distance_and_not_a_carrier() {
    for case in cases() {
        let source = recipe(&fixture(case.name));
        let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
        let mut storage = Workspace::default();
        let [t, secant] = slope(case.degrees);
        let wall = Bound::exact("1")
            .sub(&t)
            .mul(&secant.inverse())
            .times("10")
            .times(case.scale);
        for point in [["-10", "5", "10"], ["-10", "-5", "10"]] {
            let at = case.point(point);
            let (answer, work) = with_operation_budget(WORK, || {
                plan.proven_checked(&at, &mut storage, &mut Decider::default())
            });
            let answer = answer.unwrap_or_else(|e| panic!("{} {e:?} {work:?}", case.name));
            assert_eq!(answer.location, Location::Inside, "{}", case.name);
            let Value::TrueDistance(value) = answer.value else {
                panic!("{} {point:?} keeps a true distance", case.name);
            };
            assert!(
                value.lo() <= -wall.nearest() + 1e-9 && -wall.nearest() - 1e-9 <= value.hi(),
                "{} {point:?}: {value:?} misses {}",
                case.name,
                wall.nearest()
            );
            assert!(value.hi() - value.lo() < 1e-9, "{} {point:?}", case.name);
        }
    }
}

#[test]
fn trimmed_faces_are_finite_and_a_support_plane_is_not_its_face() {
    let case = &cases()[0];
    let (source, mut built) = case.built();
    let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
    let mut storage = Workspace::default();
    for (face, point, on_face) in [
        (origin(Role::End), ["-10", "5", "20"], true),
        (origin(Role::End), ["5", "5", "20"], false),
        (side("notch"), ["10", "-5", "10"], true),
        (side("notch"), ["-10", "-5", "10"], false),
    ] {
        let at = case.point(point);
        let id = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .find(|id| built.model.face_origins(*id).contains(&face))
            .unwrap();
        let native = face_contains(&mut built.model, id, &at, &mut Decider::default()).unwrap();
        assert_eq!(native, on_face, "{face:?} {point:?}");
        let state = plan
            .face_query(plan.select(&face).unwrap(), &at, &mut storage)
            .unwrap()
            .state;
        let expected = if on_face {
            FaceState::OnFace
        } else {
            FaceState::NotOnFace
        };
        assert_eq!(state, expected, "{face:?} {point:?}");
    }
}

fn overlap(e: &Enclosure, reference: &Bound, what: &str) {
    assert!(
        e.lo <= reference.hi && reference.lo <= e.hi,
        "{what}: {e:?} misses {reference:?}"
    );
}

#[test]
fn full_mass_properties_and_boxes_match_the_hand_derivation() {
    for case in cases() {
        let (_, built) = case.built();
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
        // Every extremum of the frame is one of the twelve drafted corners.
        let corners: Vec<_> = (0..6)
            .flat_map(|i| [false, true].map(|top| (i, top)))
            .map(|(i, top)| {
                let [x, y, z] = case.corner(i, top).map(|v| v.to_string());
                let p = case.point([x.as_str(), y.as_str(), z.as_str()]);
                [&p.x, &p.y, &p.z].map(Rational::to_f64_nearest)
            })
            .collect();
        let boxed = catallax_measure::bounding_box(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value
        .unwrap();
        for axis in 0..3 {
            let values = corners.iter().map(|corner| corner[axis]);
            let low = values.clone().fold(f64::INFINITY, f64::min);
            let high = values.fold(f64::NEG_INFINITY, f64::max);
            assert!((boxed.min[axis].to_f64_nearest() - low).abs() <= 1e-3);
            assert!((boxed.max[axis].to_f64_nearest() - high).abs() <= 1e-3);
        }
    }
}

#[test]
fn published_table_rows_are_the_hand_formulas() {
    fn agree(bound: &Bound, row: &str) {
        let middle = (&bound.lo + &bound.hi).mul_pow2(-1);
        assert!((&middle - &q(row)).abs() < q("1e-8"), "{row}");
    }
    for (degrees, volume, area, centroid, inertia) in [
        (
            "3",
            "15371.106648603506",
            "4296.831256631423",
            ["-4.931810011496", "-2.602284982756", "9.863620022992"],
            [
                "1694159.259419132",
                "2109589.389474316",
                "2779579.995391907",
                "576309.288746627",
                "-3636.799386876",
                "5455.199080314",
            ],
        ),
        (
            "5",
            "14950.136037688912",
            "4229.606478000169",
            ["-4.882959374677", "-2.675560937985", "9.765918749354"],
            [
                "1645089.788332550",
                "2010221.209052094",
                "2660273.618000760",
                "560322.909576701",
                "-6242.166683906",
                "9363.250025858",
            ],
        ),
    ] {
        let reference = Reference::local(degrees);
        agree(&reference.volume, volume);
        agree(&reference.area, area);
        for (bound, row) in reference.centroid.iter().zip(centroid) {
            agree(bound, row);
        }
        for ((i, j), row) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .zip(inertia)
        {
            agree(&reference.inertia[i][j], row);
            agree(&reference.inertia[j][i], row);
        }
    }
}

#[test]
fn certified_display_is_closed_and_selectable() {
    let case = &cases()[0];
    let (_, built) = case.built();
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: q("0.005"),
        ..Default::default()
    };
    let (mesh, work) = with_operation_budget(WORK, || {
        catallax_mesh::mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted, "{work:?}");
    let mesh = mesh.unwrap().value;
    assert_eq!(mesh.faces.len(), 8);
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
    let mut regions = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for triangle in &mesh.triangles {
        *regions.entry(triangle.face).or_insert(0) += 1;
        // No triangle of the bracket's material fills the notch.
        let p = triangle.vertices.map(|i| mesh.positions[i]);
        let centroid: [f64; 3] = std::array::from_fn(|i| (p[0][i] + p[1][i] + p[2][i]) / 3.0);
        assert!(centroid[0] < 0.1 || centroid[1] < -4.9, "{centroid:?}");
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
            counts[usize::from(a > b)] += 1;
        }
    }
    assert_eq!(regions.len(), 8);
    assert!(edges.values().all(|counts| *counts == [1, 1]));
}

fn with_region(point: [&str; 3]) -> Json {
    let mut source = fixture("R-P3");
    source["steps"][2]["keep_regions"] = json!([point]);
    source
}

#[test]
fn a_saved_region_keeps_the_whole_connected_bracket_and_never_a_tile() {
    let case = &cases()[0];
    let picked = recipe(&with_region(["-10", "-10", "10"]));
    let (built, work) = with_operation_budget(WORK, || evaluate(&picked, &mut Decider::default()));
    assert!(!work.exhausted, "{work:?}");
    let mut built = built.unwrap();
    check_topology(case, &built);
    check_edge_names(case, &built);
    let plan = Plan::compile_checked(&picked, &mut Decider::default()).unwrap();
    let mut storage = Workspace::default();
    for (point, expected, _) in MATERIAL {
        check_location(case, point, expected, (&plan, &mut storage), &mut built);
    }
    let properties = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    overlap(&properties.volume, &case.reference().volume, "volume");
}

#[test]
fn a_notch_or_boundary_region_is_refused_by_both_kernels() {
    for point in [["5", "5", "10"], ["20", "-15", "0"]] {
        let source = recipe(&with_region(point));
        let (built, _) = with_operation_budget(WORK, || evaluate(&source, &mut Decider::default()));
        let Err(BuildError::Geometry { refusal, .. }) = built else {
            panic!("B-rep must refuse {point:?}");
        };
        assert_eq!(refusal.code, "draft.invalid-keep-region", "{point:?}");
        let (plan, _) = with_operation_budget(WORK, || {
            Plan::compile_checked(&source, &mut Decider::default())
        });
        let field = plan
            .err()
            .unwrap_or_else(|| panic!("field must refuse {point:?}"));
        assert_eq!(field.code, "draft.invalid-keep-region", "{point:?}");
    }
}
