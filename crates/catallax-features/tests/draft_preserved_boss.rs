//! Planar draft under a preserved raised boss, checked in both kernels.
//!
//! Every expected value comes from the part's definition. No expectation reads
//! a kernel result, and neither kernel's answer supplies the other's.
use catallax_boolean::{Containment, classify_point, face_contains};
use catallax_features::{BuildError, Built, ConsumedFace, evaluate};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Integer, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::{Curve, Point3, RPoint3, Surface, Vector3};
use catallax_topology::{BodyKind, FaceId, Loop, Sense};
use serde_json::{Value as Json, json};
use std::collections::{BTreeMap, BTreeSet};

/// One ordinary application request's allowance.
const WORK: u64 = 20_000_000;

fn q(x: &str) -> Rational {
    Rational::parse(x).unwrap()
}
fn fixture(name: &str) -> Json {
    let path = format!(
        "{}/tests/fixtures/draft-preserved-boss/{name}.json",
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
    fn exact(value: &Rational) -> Self {
        Self::new(value.clone(), value.clone())
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
        self.mul(&Self::exact(&q(value)))
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
/// Tangent and secant of a positive angle in degrees. At 45 degrees the tangent
/// is exactly one and the secant lies between a rational `s` below the square
/// root of two and `2/s` above it.
fn slope(degrees: &str) -> [Bound; 2] {
    if degrees == "45" {
        let s = q("1.4142135623730950488016887242096980785696718753769480731766797379907");
        assert!(&s * &s < q("2"));
        let two = q("2");
        return [Bound::exact(&q("1")), Bound::new(s.clone(), &two / &s)];
    }
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
    /// The base section `30(2W-2zt)` over 0 <= z <= 20 and the radius-4, height-5
    /// boss about (6,2,45/2), with radial variance 4 and axial variance 25/12. The
    /// boss's end replaces the shoulder disk it covers, so only its side adds area.
    fn local(half_width: &str, degrees: &str) -> Self {
        let [t, secant] = slope(degrees);
        let w = q(half_width);
        let (w2, w3) = (&w * &w, &(&w * &w) * &w);
        let scaled = |x: &Rational, by: &str| Bound::exact(&(x * &q(by)));
        let (t2, t3) = (t.square(), t.square().mul(&t));
        let boss = pi().times("80");
        let base = scaled(&w, "1200").sub(&t.times("12000"));
        let volume = base.add(&boss);
        let area = scaled(&w, "200")
            .sub(&t.times("2000"))
            .add(&secant.times("1200"))
            .add(&pi().times("40"));
        let first = [
            boss.times("6"),
            boss.times("2"),
            scaled(&w, "12000")
                .sub(&t.times("160000"))
                .add(&boss.times("45/2")),
        ];
        let second = [
            scaled(&w3, "400")
                .sub(&t.mul(&scaled(&w2, "12000")))
                .add(&t2.mul(&scaled(&w, "160000")))
                .sub(&t3.times("800000"))
                .add(&boss.times("40")),
            base.times("75").add(&boss.times("8")),
            scaled(&w, "160000")
                .sub(&t.times("2400000"))
                .add(&boss.times("1525/3")),
        ];
        let cross = [boss.times("12"), boss.times("135"), boss.times("45")];
        let per_volume = volume.inverse();
        let squared = |i: usize| first[i].square().mul(&per_volume);
        let diagonal =
            |a: usize, b: usize| second[a].add(&second[b]).sub(&squared(a).add(&squared(b)));
        let off = |a: usize, b: usize, q: &Bound| first[a].mul(&first[b]).mul(&per_volume).sub(q);
        let (xy, xz, yz) = (
            off(0, 1, &cross[0]),
            off(0, 2, &cross[1]),
            off(1, 2, &cross[2]),
        );
        Self {
            centroid: first.clone().map(|m| m.mul(&per_volume)),
            volume,
            area,
            inertia: [
                [diagonal(1, 2), xy.clone(), xz.clone()],
                [xy, diagonal(0, 2), yz.clone()],
                [xz, yz, diagonal(0, 1)],
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
            .map(|row| row.map(|x| Bound::exact(&x)));
        let centroid = std::array::from_fn(|i| {
            (0..3).fold(Bound::exact(&frame.translation[i]), |sum, j| {
                sum.add(&m[i][j].mul(&self.centroid[j]))
            })
        });
        let inertia = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Bound::exact(&Rational::zero()), |sum, a| {
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
    half_width: &'static str,
    degrees: &'static str,
    scale: &'static str,
    frame: Option<Frame>,
    /// Whether (19,0,19) keeps inside this slope: `19t < W-19`.
    edit_point: Location,
}

fn cases() -> Vec<Case> {
    let case = |name, half_width, degrees, scale, frame, edit_point| Case {
        name,
        half_width,
        degrees,
        scale,
        frame,
        edit_point,
    };
    vec![
        case("M-N", "20", "3", "1", None, Location::Inside),
        case("M-E", "20", "5", "1", None, Location::Outside),
        case("M-R", "20", "3", "1", Some(oblique()), Location::Inside),
        case("M-S", "31", "45", "1/2", None, Location::Outside),
    ]
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
        let reference = Reference::local(self.half_width, self.degrees).scaled(self.scale);
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

fn origin(step: &str, role: Role) -> Origin {
    Origin {
        step: StepId::new(step).unwrap(),
        role,
    }
}
fn side(step: &str, segment: &str) -> Origin {
    origin(
        step,
        Role::Side {
            segment: ElementId::new(segment).unwrap(),
            path: None,
        },
    )
}
fn shoulder() -> Origin {
    origin("base", Role::End)
}
fn top() -> Origin {
    origin("boss", Role::End)
}
fn named(built: &Built, name: &Origin) -> FaceId {
    let faces: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|face| built.model.face_origins(*face).contains(name))
        .collect();
    let [face] = faces.as_slice() else {
        panic!("{name:?} names {faces:?}");
    };
    *face
}

/// Exact supports, senses, rings, adjacency and names of the admitted topology.
fn check_topology(case: &Case, built: &Built) {
    let model = &built.model;
    assert!(model.check(built.result).is_empty(), "{}", case.name);
    let counts = model.counts(built.result);
    assert_eq!(
        (counts.lumps, counts.shells, counts.faces, counts.loops),
        (1, 1, 8, 10),
        "{}",
        case.name
    );
    assert_eq!(
        (counts.edges, counts.ring_edges, counts.vertices),
        (14, 2, 8)
    );
    assert!(built.consumed.is_empty());
    let degrees = q(case.degrees).to_f64_nearest();
    let t = (degrees * std::f64::consts::PI / 180.0).tan();
    let w = q(case.half_width).to_f64_nearest();
    let planes = [
        (origin("base", Role::Start), [0.0, 0.0, -1.0], 0.0, 1),
        (shoulder(), [0.0, 0.0, 1.0], 20.0, 2),
        (top(), [0.0, 0.0, 1.0], 25.0, 1),
        (side("base", "bottom"), [0.0, -1.0, 0.0], 15.0, 1),
        (side("base", "top"), [0.0, 1.0, 0.0], 15.0, 1),
        (side("base", "right"), [1.0, 0.0, t], w, 1),
        (side("base", "left"), [-1.0, 0.0, t], w, 1),
    ];
    let mut wall_rings = Vec::new();
    let mut plane_rings = Vec::new();
    let mut names = BTreeSet::new();
    let mut disks = 0;
    for face in model.faces_of_body(built.result) {
        let f = model.face(face);
        let origins = model.face_origins(face);
        assert_eq!(origins.len(), 1, "{origins:?}");
        assert!(names.insert(origins[0].clone()), "{origins:?}");
        let loops: Vec<_> = f
            .loops
            .iter()
            .map(|lp| match model.loop_(*lp) {
                Loop::Coedges(cs) => cs.clone(),
                Loop::Vertex(_) => panic!("no vertex loops"),
            })
            .collect();
        match model.geom().surface(f.surface) {
            Surface::Plane(plane) => {
                let (_, normal, offset, wires) = planes
                    .iter()
                    .find(|(name, ..)| *name == origins[0])
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
                assert_eq!(loops.len(), *wires, "{origins:?}");
                if *wires == 1 {
                    disks += 1;
                }
                for wire in &loops {
                    match wire.as_slice() {
                        [ring] => plane_rings.push(model.coedge(*ring).edge),
                        wire => assert_eq!(wire.len(), 4, "{origins:?}"),
                    }
                }
            }
            Surface::Cylinder(cylinder) => {
                assert_eq!(origins[0], side("boss", "wall"));
                assert_eq!(f.sense, Sense::Same, "the wall bounds material on its axis");
                let scaled = &q("4") * &q(case.scale);
                let radius = &scaled * &scaled;
                assert_eq!(cylinder.radius_sq.rational_ref(), Some(&radius));
                let axis = case.local(vector(&cylinder.axis), false);
                assert!(axis[0].abs() < 1e-9 && axis[1].abs() < 1e-9 && axis[2] != 0.0);
                let at = case.local(xyz(&cylinder.origin), true);
                assert!((at[0] - 6.0).abs() < 1e-9 && (at[1] - 2.0).abs() < 1e-9);
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
                    assert!((centre[0] - 6.0).abs() < 1e-9 && (centre[1] - 2.0).abs() < 1e-9);
                    heights.push((centre[2] * 1e6).round() / 1e6);
                }
                heights.sort_by(f64::total_cmp);
                assert_eq!(heights, [20.0, 25.0]);
            }
            surface => panic!("unexpected {surface:?}"),
        }
    }
    // Eight distinct names, none of them the buried `boss.Start`.
    assert_eq!(names.len(), 8);
    assert!(!names.contains(&origin("boss", Role::Start)));
    wall_rings.sort();
    plane_rings.sort();
    assert_eq!(wall_rings, plane_rings);
    // Six disks, one annulus and one band: rings add nothing, so V-E+F is 2 and
    // the shell has genus zero. A disk left in the joint would break the count.
    assert_eq!(disks, 6);
    assert_eq!(
        counts.vertices as i64 - (counts.edges - counts.ring_edges) as i64 + disks,
        2
    );
}

/// Every row but the two that depend on the slope and the part's width.
const MATERIAL: [([&str; 3], Location); 7] = [
    (["0", "0", "10"], Location::Inside),
    (["6", "2", "22"], Location::Inside),
    (["11", "2", "22"], Location::Outside),
    (["6", "2", "20"], Location::Inside),
    (["10", "2", "22"], Location::Boundary),
    (["6", "2", "25"], Location::Boundary),
    (["-10", "0", "20"], Location::Boundary),
];

#[test]
fn both_kernels_draft_the_base_under_an_unchanged_boss() {
    for case in cases() {
        let source = recipe(&fixture(case.name));
        let (built, _) = with_operation_budget(WORK, || evaluate(&source, &mut Decider::default()));
        let mut built = built.unwrap_or_else(|e| panic!("{} {e:?}", case.name));
        check_topology(&case, &built);
        let (plan, _) = with_operation_budget(WORK, || {
            Plan::compile_checked(&source, &mut Decider::default())
        });
        let plan = plan.unwrap_or_else(|e| panic!("{} {e:?}", case.name));
        assert!(plan.consumed_faces().is_empty());
        let beyond = (&q(case.half_width) + &q("1")).to_string();
        let mut rows = MATERIAL.to_vec();
        rows.push((["19", "0", "19"], case.edit_point));
        rows.push(([beyond.as_str(), "0", "10"], Location::Outside));
        let mut storage = Workspace::default();
        for (point, expected) in rows {
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
        // A carrier is not a face: the shoulder stops at the joint, the top at
        // the boss's rim.
        let (shoulder_face, top_face) = (named(&built, &shoulder()), named(&built, &top()));
        for (name, face, point, on) in [
            (shoulder(), shoulder_face, ["6", "2", "20"], false),
            (shoulder(), shoulder_face, ["-10", "0", "20"], true),
            (top(), top_face, ["6", "2", "25"], true),
            (top(), top_face, ["11", "2", "25"], false),
        ] {
            let at = case.point(point);
            let contains =
                face_contains(&mut built.model, face, &at, &mut Decider::default()).unwrap();
            assert_eq!(contains, on, "{} native {name:?} {point:?}", case.name);
            let state = plan
                .face_query(plan.select(&name).unwrap(), &at, &mut storage)
                .unwrap()
                .state;
            let expected = if on {
                FaceState::OnFace
            } else {
                FaceState::NotOnFace
            };
            assert_eq!(state, expected, "{} field {name:?} {point:?}", case.name);
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

fn agree(bound: &Bound, row: &str) {
    let middle = (&bound.lo + &bound.hi).mul_pow2(-1);
    assert!((&middle - &q(row)).abs() < q("1e-8"), "{row}");
}

#[test]
fn bounds_and_published_table_rows_match_the_definition() {
    // The rounded report rows must agree with these formulas to their digits.
    for (w, degrees, scale, row) in [
        (
            "20",
            "3",
            "1",
            [
                "23622.434060890689",
                "5222.494962775014",
                "0.063836117389",
                "0.021278705796",
                "10.044249524819",
                "2573617.064122793",
                "3785857.522882187",
                "4721906.302958257",
                "-2983.841415059",
                "-18782.829210132",
                "-6260.943070044",
            ],
        ),
        (
            "20",
            "5",
            "1",
            [
                "23201.463449976095",
                "5155.270184143761",
                "0.064994368867",
                "0.021664789622",
                "9.984571926326",
                "2528052.451210093",
                "3618759.944130454",
                "4537227.357348239",
                "-2983.259214365",
                "-18872.820908537",
                "-6290.940302846",
            ],
        ),
        (
            "31",
            "45",
            "1/2",
            [
                "3181.415926535898",
                "1505.679995247826",
                "0.029624475952",
                "0.009874825317",
                "4.275904027532",
                "84950.504698886",
                "168149.350296879",
                "201448.888850564",
                "-93.317099248",
                "-657.293060176",
                "-219.097686725",
            ],
        ),
    ] {
        let reference = Reference::local(w, degrees).scaled(scale);
        let [x, y, z] = &reference.centroid;
        let i = &reference.inertia;
        for (bound, row) in [
            &reference.volume,
            &reference.area,
            x,
            y,
            z,
            &i[0][0],
            &i[1][1],
            &i[2][2],
            &i[0][1],
            &i[0][2],
            &i[1][2],
        ]
        .into_iter()
        .zip(row)
        {
            agree(bound, row);
        }
    }
    // Bottom edges and the boss's end retain every extremum of the local box.
    for case in cases().iter().filter(|case| case.frame.is_none()) {
        let built = evaluate(&recipe(&fixture(case.name)), &mut Decider::default()).unwrap();
        let boxed = catallax_measure::bounding_box(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value
        .unwrap();
        let w = q(case.half_width);
        let scale = q(case.scale);
        for (axis, (lo, hi)) in [(-&w, w.clone()), (q("-15"), q("15")), (q("0"), q("25"))]
            .iter()
            .enumerate()
        {
            assert!((&boxed.min[axis] - &(lo * &scale)).abs() <= q("0.001"));
            assert!((&boxed.max[axis] - &(hi * &scale)).abs() <= q("0.001"));
        }
    }
}

#[test]
fn certified_display_is_closed_selectable_and_keeps_the_joint_open() {
    let source = recipe(&fixture("M-N"));
    let built = evaluate(&source, &mut Decider::default()).unwrap();
    let options = catallax_mesh::MeshOptions::default();
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
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
    let chord = options.chord_tolerance.to_f64_nearest();
    let shoulder_face = named(&built, &shoulder());
    let mut regions = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for triangle in &mesh.triangles {
        *regions.entry(triangle.face).or_insert(0) += 1;
        if triangle.face == shoulder_face {
            let p = triangle.vertices.map(|i| mesh.positions[i]);
            let x = (p[0][0] + p[1][0] + p[2][0]) / 3.0;
            let y = (p[0][1] + p[1][1] + p[2][1]) / 3.0;
            assert!((x - 6.0).hypot(y - 2.0) >= 4.0 - chord, "joint at {x} {y}");
        }
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
            counts[usize::from(a > b)] += 1;
        }
    }
    let names: BTreeSet<_> = regions
        .keys()
        .flat_map(|face| built.model.face_origins(*face))
        .collect();
    assert_eq!((regions.len(), names.len()), (8, 8));
    assert!(edges.values().all(|counts| *counts == [1, 1]));
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

fn coverage() -> (String, String) {
    (
        "draft.unsupported-arrangement".into(),
        "draft.unsupported-source-boundary".into(),
    )
}

#[test]
fn contact_crossing_and_invalid_angle_keep_their_typed_refusals() {
    for name in ["M-T", "M-X"] {
        assert_eq!(refusal(&fixture(name)), coverage(), "{name}");
    }
    assert_eq!(
        refusal(&fixture("M-BAD")),
        ("draft.angle".into(), "draft.angle".into())
    );
}

const DRAFT: usize = 5;

fn edited(edit: &dyn Fn(&mut Json)) -> Json {
    let mut source = fixture("M-N");
    edit(&mut source);
    source
}

#[test]
fn the_new_branch_admits_only_an_unselected_axial_boss() {
    let with = |face: Json| {
        edited(&|s| {
            s["steps"][DRAFT]["faces"]
                .as_array_mut()
                .unwrap()
                .push(json!({"face": face}));
        })
    };
    for source in [
        with(json!({"step":"boss","role":{"side":{"segment":"wall"}}})),
        with(json!({"step":"boss","role":"end"})),
        with(json!({"step":"base","role":"end"})),
        // Pull and neutral each leave the boss's axis.
        edited(&|s| s["steps"][DRAFT]["pull"] = json!(["0", "1", "20"])),
        edited(&|s| s["steps"][DRAFT]["neutral"]["normal"] = json!(["0", "1", "20"])),
    ] {
        assert_eq!(refusal(&source), coverage(), "{source}");
    }
}

#[test]
fn saved_regions_name_the_boss_or_nothing() {
    let kept = recipe(&edited(&|s| {
        s["steps"][DRAFT]["keep_regions"] = json!([["6", "2", "22"]]);
    }));
    let built = evaluate(&kept, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 8);
    let plan = Plan::compile_checked(&kept, &mut Decider::default()).unwrap();
    let at = RPoint3::new(q("0"), q("0"), q("10"));
    let field = plan
        .proven_checked(&at, &mut Workspace::default(), &mut Decider::default())
        .unwrap();
    assert_eq!(field.location, Location::Inside, "the base stays attached");
    let outside = edited(&|s| {
        s["steps"][DRAFT]["keep_regions"] = json!([["11", "2", "22"]]);
    });
    assert_eq!(
        refusal(&outside),
        (
            "draft.invalid-keep-region".into(),
            "draft.invalid-keep-region".into()
        )
    );
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
fn a_neutral_plane_through_the_boss_anchors_and_one_past_it_empties_both_kernels() {
    let across = recipe(&edited(&|s| {
        s["steps"][DRAFT]["neutral"]["point"][2] = json!("22");
    }));
    let built = evaluate(&across, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 8);
    assert!(
        !Plan::compile_checked(&across, &mut Decider::default())
            .unwrap()
            .is_proven_empty()
    );
    let past: &dyn Fn(&mut Json) = &|s| s["steps"][DRAFT]["neutral"]["point"][2] = json!("100");
    let none: &dyn Fn(&mut Json) = &|s| s["steps"][DRAFT]["keep_regions"] = json!([]);
    for edit in [past, none] {
        let source = recipe(&edited(edit));
        let built = evaluate(&source, &mut Decider::default()).unwrap();
        assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
        let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
        assert!(plan.is_proven_empty());
        let native = history(&built.consumed);
        assert_eq!(native.len(), 8);
        assert_eq!(history(plan.consumed_faces()), native);
    }
}

#[test]
fn union_operand_order_does_not_change_the_answer() {
    let source = recipe(&edited(&|s| {
        s["steps"][4]["target"] = json!("boss");
        s["steps"][4]["tool"] = json!("base");
    }));
    let built = evaluate(&source, &mut Decider::default()).unwrap();
    check_topology(&cases().remove(0), &built);
    Plan::compile_checked(&source, &mut Decider::default()).unwrap();
}

#[test]
fn a_boss_under_the_base_stands_on_the_neutral_cap() {
    // The sketch turns over, so the boss runs from z=0 down to z=-5.
    let source = recipe(&edited(&|s| {
        s["steps"][2]["plane"]["rotation"] = json!(["0", "1", "0", "0"]);
        s["steps"][2]["plane"]["translation"][2] = json!("0");
    }));
    let mut built = evaluate(&source, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 8);
    let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
    let mut storage = Workspace::default();
    for (point, expected) in [
        (["6", "2", "-3"], Location::Inside),
        (["6", "2", "0"], Location::Inside),
        (["6", "2", "22"], Location::Outside),
    ] {
        let at = RPoint3::new(q(point[0]), q(point[1]), q(point[2]));
        let field = plan
            .proven_checked(&at, &mut storage, &mut Decider::default())
            .unwrap();
        assert_eq!(field.location, expected, "field {point:?}");
        let native = classify_point(&mut built.model, built.result, &at, &mut Decider::default());
        assert!(
            match expected {
                Location::Inside => matches!(native, Outcome::Proven(Containment::Inside)),
                _ => matches!(native, Outcome::Proven(Containment::Outside)),
            },
            "native {point:?}: {native:?}"
        );
    }
}

#[test]
fn a_tool_sunk_into_its_target_keeps_the_field_refusal() {
    // The B-rep sees the same ring and answers. The field holds the whole tool, so
    // it would have to prove the buried length clear of the drafted walls.
    let source = recipe(&edited(&|s| {
        s["steps"][2]["plane"]["translation"][2] = json!("18");
        s["steps"][3]["distance"] = json!("7");
    }));
    let built = evaluate(&source, &mut Decider::default()).unwrap();
    assert_eq!(built.model.counts(built.result).faces, 8);
    let field = Plan::compile_checked(&source, &mut Decider::default())
        .err()
        .unwrap();
    assert_eq!(field.code, "draft.unsupported-source-boundary");
}
