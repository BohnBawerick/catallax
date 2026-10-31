//! Planar draft around one separated axial bore, checked in both kernels.
//!
//! Every expected value comes from the part's definition. No expectation reads
//! a kernel result, and neither kernel's answer supplies the other's.
use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, Built, ConsumedFace, evaluate};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Integer, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{ElementId, HolePart, Origin, Recipe, Role, StepId};
use catallax_shapes::{Curve, Point3, RPoint3, Surface, Vector3};
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
    /// Section area `30(40-2zt)-25pi` for z in [0,20] and its moments.
    fn local(degrees: &str) -> Self {
        let pi = pi();
        let [t, secant] = slope(degrees);
        let t2 = t.square();
        let t3 = t2.mul(&t);
        let volume = Bound::exact("24000")
            .sub(&t.times("12000"))
            .sub(&pi.times("500"));
        let area = Bound::exact("4000")
            .sub(&t.times("2000"))
            .add(&secant.times("1200"))
            .add(&pi.times("150"));
        let z = Bound::exact("240000")
            .sub(&t.times("160000"))
            .sub(&pi.times("5000"));
        let x2 = Bound::exact("3200000")
            .sub(&t.times("4800000"))
            .add(&t2.times("3200000"))
            .sub(&t3.times("800000"))
            .sub(&pi.times("3125"));
        let y2 = Bound::exact("1800000")
            .sub(&t.times("900000"))
            .sub(&pi.times("3125"));
        let z2 = Bound::exact("3200000")
            .sub(&t.times("2400000"))
            .sub(&pi.times("200000/3"));
        let shift = z.square().mul(&volume.inverse());
        let zero = Bound::exact("0");
        Self {
            centroid: [zero.clone(), zero.clone(), z.mul(&volume.inverse())],
            volume,
            area,
            inertia: [
                [y2.add(&z2).sub(&shift), zero.clone(), zero.clone()],
                [zero.clone(), x2.add(&z2).sub(&shift), zero.clone()],
                [zero.clone(), zero.clone(), x2.add(&y2)],
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

/// K-R's exact proper rotation from quaternion (1,2,3,4), and its translation.
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
    wall: Origin,
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
fn cases() -> Vec<Case> {
    let native = Origin {
        step: StepId::new("hole").unwrap(),
        role: Role::Hole {
            at: ElementId::new("centre").unwrap(),
            part: HolePart::Wall,
        },
    };
    let tool = || side("tool", "wall");
    vec![
        Case {
            name: "K-H",
            degrees: "3",
            scale: "1",
            frame: None,
            wall: native,
        },
        Case {
            name: "K-E",
            degrees: "5",
            scale: "1",
            frame: None,
            wall: tool(),
        },
        Case {
            name: "K-R",
            degrees: "3",
            scale: "1",
            frame: Some(oblique()),
            wall: tool(),
        },
        Case {
            name: "K-S",
            degrees: "5",
            scale: "1/2",
            frame: None,
            wall: tool(),
        },
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
        let reference = Reference::local(self.degrees).scaled(self.scale);
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
    assert_eq!(
        (counts.lumps, counts.shells, counts.faces, counts.loops),
        (1, 1, 7, 10)
    );
    assert_eq!(
        (counts.edges, counts.ring_edges, counts.vertices),
        (14, 2, 8)
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
                    let [outer, ring] = loops.as_slice() else {
                        panic!("a cap keeps its outer rectangle and hole");
                    };
                    assert_eq!(outer.len(), 4);
                    let [ring] = ring.as_slice() else {
                        panic!("one circular hole edge");
                    };
                    cap_rings.push(model.coedge(*ring).edge);
                } else {
                    assert_eq!(loops.len(), 1);
                    assert_eq!(loops[0].len(), 4);
                    disks += 1;
                }
            }
            Surface::Cylinder(cylinder) => {
                assert_eq!(origins[0], case.wall);
                assert_eq!(f.sense, Sense::Reversed);
                let radius = &q("25") * &(&q(case.scale) * &q(case.scale));
                assert_eq!(cylinder.radius_sq.rational_ref(), Some(&radius));
                let axis = case.local(vector(&cylinder.axis), false);
                assert!(axis[0].abs() < 1e-9 && axis[1].abs() < 1e-9 && axis[2] != 0.0);
                let origin = case.local(xyz(&cylinder.origin), true);
                assert!(origin[0].abs() < 1e-9 && origin[1].abs() < 1e-9);
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
                    assert!(centre[0].abs() < 1e-9 && centre[1].abs() < 1e-9);
                    heights.push((centre[2] * 1e6).round() / 1e6);
                }
                heights.sort_by(f64::total_cmp);
                assert_eq!(heights, [0.0, 20.0]);
            }
            surface => panic!("unexpected {surface:?}"),
        }
    }
    assert_eq!(seen, 7);
    // Caps and wall are annuli. The four side disks close `V-E+F=0` on the line
    // graph, with ring edges contributing nothing: genus one, not a filled body.
    assert_eq!(disks, 4);
    assert_eq!(
        counts.vertices as i64 - (counts.edges - counts.ring_edges) as i64 + disks,
        0
    );
    wall_rings.sort();
    cap_rings.sort();
    assert_eq!(wall_rings, cap_rings);
}

const MATERIAL: [([&str; 3], Location, Location); 7] = [
    (["10", "0", "10"], Location::Inside, Location::Inside),
    (["0", "0", "10"], Location::Outside, Location::Outside),
    (["19", "0", "19"], Location::Inside, Location::Outside),
    (["21", "0", "10"], Location::Outside, Location::Outside),
    (["5", "0", "10"], Location::Boundary, Location::Boundary),
    (["0", "0", "20"], Location::Outside, Location::Outside),
    (["10", "0", "20"], Location::Boundary, Location::Boundary),
];

#[test]
fn both_kernels_draft_the_walls_around_an_unchanged_bore() {
    for case in cases() {
        let source = recipe(&fixture(case.name));
        let mut built = evaluate(&source, &mut Decider::default()).unwrap();
        check_topology(&case, &built);
        let plan = Plan::compile_checked(&source, &mut Decider::default()).unwrap();
        assert!(plan.consumed_faces().is_empty());
        let mut storage = Workspace::default();
        for (point, three, five) in MATERIAL {
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
        assert_eq!(
            state(&case.wall, ["5", "0", "10"], &mut storage),
            FaceState::OnFace
        );
        assert_eq!(
            state(&end, ["10", "0", "20"], &mut storage),
            FaceState::OnFace
        );
        assert_eq!(
            state(&end, ["0", "0", "20"], &mut storage),
            FaceState::NotOnFace
        );
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
        let answer = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
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
        let reference = Reference::local(degrees);
        for (bound, row) in [
            (&reference.volume, volume),
            (&reference.area, area),
            (&reference.centroid[2], height),
        ] {
            let middle = (&bound.lo + &bound.hi).mul_pow2(-1);
            assert!((&middle - &q(row)).abs() < q("1e-9"), "{degrees} {row}");
        }
    }
    // Bottom edges retain every extremum, so drafting keeps the local box.
    let built = evaluate(&recipe(&fixture("K-H")), &mut Decider::default()).unwrap();
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

#[test]
fn certified_display_is_closed_selectable_and_keeps_the_cap_holes_open() {
    let source = recipe(&fixture("K-H"));
    let built = evaluate(&source, &mut Decider::default()).unwrap();
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: q("0.005"),
        ..Default::default()
    };
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(mesh.faces.len(), 7);
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
            assert!(x.hypot(y) >= 5.0 - 0.005, "filled cap hole at {x} {y}");
        }
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
            counts[usize::from(a > b)] += 1;
        }
    }
    assert_eq!(regions.len(), 7);
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
    for (name, draft, wall) in [
        ("K-E", 4, side("tool", "wall")),
        (
            "K-H",
            3,
            Origin {
                step: StepId::new("hole").unwrap(),
                role: Role::Hole {
                    at: ElementId::new("centre").unwrap(),
                    part: HolePart::Wall,
                },
            },
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
            assert_eq!(native.len(), 7, "{name}");
            assert!(native.iter().any(|(removed, _)| removed.contains(&wall)));
            assert_eq!(history(plan.consumed_faces()), native, "{name}");
        }
    }
}
