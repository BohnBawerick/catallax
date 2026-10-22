//! Timings and observable stored-number growth, never performance gates.

mod common;

use catallax_boolean::{
    BoolOp, Stage, boolean, classify, classify_point, imprint, intersect, keep, sew,
};
use catallax_numbers::{Decider, DecisionStats, Rational, Real};
use catallax_shapes::{
    Curve, CurveId, MeetOf, Point3, PointDef, PointId, Pose, Position, RPoint3, RVector3, Rotation,
    Surface, SurfaceId, Vector3,
};
use catallax_topology::{BodyId, Model};
use common::{Case, Primitive, ball_at, block_at, cylinder_at, r};
use std::collections::BTreeSet;
use std::time::Instant;

#[derive(Default)]
struct NumberSize {
    rational_bits: u64,
    numerator_bits: u64,
    denominator_bits: u64,
    opaque_real_coordinates: usize,
}
impl NumberSize {
    fn rational(&mut self, r: &Rational) {
        let n = r.numerator().bit_length();
        let d = r.denominator().bit_length();
        self.rational_bits = self.rational_bits.max(n + d);
        self.numerator_bits = self.numerator_bits.max(n);
        self.denominator_bits = self.denominator_bits.max(d);
    }
    fn real(&mut self, value: &Real) {
        if let Some(value) = value.rational_ref() {
            self.rational(value);
        } else {
            self.opaque_real_coordinates += 1;
        }
    }
    fn point(&mut self, p: &Point3) {
        for value in [&p.x, &p.y, &p.z] {
            self.real(value);
        }
    }
    fn vector(&mut self, v: &Vector3) {
        for value in [&v.x, &v.y, &v.z] {
            self.real(value);
        }
    }
    fn polynomial(&mut self, p: &catallax_numbers::RealPoly) {
        for (_, coefficient) in p.terms() {
            self.real(&coefficient);
        }
    }
    fn rational_polynomial(&mut self, p: &catallax_numbers::MPoly) {
        for (_, coefficient) in p.terms() {
            self.rational(&coefficient);
        }
    }
}
fn inventory(model: &Model, body: BodyId) -> NumberSize {
    let mut result = NumberSize::default();
    let mut surfaces = Vec::<SurfaceId>::new();
    let mut curves = Vec::<CurveId>::new();
    let mut points = Vec::<PointId>::new();
    for f in model.faces_of_body(body) {
        surfaces.push(model.face(f).surface);
    }
    for e in model.edges_of_body(body) {
        curves.push(model.edge(e).curve);
    }
    for v in model.vertices_of_body(body) {
        points.push(model.vertex(v).point);
    }
    for mark in model.body(body).marks.iter() {
        result.rational(&mark.bound);
    }
    let mut seen_surfaces = BTreeSet::new();
    let mut seen_curves = BTreeSet::new();
    let mut seen_points = BTreeSet::new();
    while !surfaces.is_empty() || !curves.is_empty() || !points.is_empty() {
        while let Some(surface) = surfaces.pop() {
            if !seen_surfaces.insert(surface) {
                continue;
            }
            result.polynomial(&model.geom().implicit(surface));
            match model.geom().surface(surface) {
                Surface::Plane(s) => {
                    result.point(&s.origin);
                    result.vector(&s.normal);
                    result.vector(&s.ref_dir);
                }
                Surface::Cylinder(s) => {
                    result.point(&s.origin);
                    result.vector(&s.axis);
                    result.vector(&s.ref_dir);
                    result.real(&s.radius);
                }
                Surface::Sphere(s) => {
                    result.point(&s.center);
                    result.vector(&s.axis);
                    result.vector(&s.ref_dir);
                    result.real(&s.radius);
                }
                Surface::Cone(s) => {
                    result.point(&s.apex);
                    result.vector(&s.axis);
                    result.vector(&s.ref_dir);
                    result.real(&s.slope);
                }
                Surface::Torus(s) => {
                    result.point(&s.center);
                    result.vector(&s.axis);
                    result.vector(&s.ref_dir);
                    result.real(&s.major);
                    result.real(&s.minor);
                }
                Surface::EllipticCylinder(s) => {
                    result.point(&s.origin);
                    result.vector(&s.axis);
                    result.vector(&s.major_dir);
                    result.real(&s.major_sq);
                    result.real(&s.minor_sq);
                }
                _ => {}
            }
        }
        while let Some(curve) = curves.pop() {
            if !seen_curves.insert(curve) {
                continue;
            }
            for equation in model.geom().curve(curve).equations(model.geom()) {
                result.polynomial(&equation);
            }
            match model.geom().curve(curve) {
                Curve::Line(c) => {
                    result.point(&c.origin);
                    result.vector(&c.dir);
                }
                Curve::Circle(c) => {
                    result.point(&c.center);
                    result.vector(&c.normal);
                    result.vector(&c.ref_dir);
                    result.real(&c.radius_sq);
                }
                Curve::Ellipse(c) => {
                    result.point(&c.center);
                    result.vector(&c.normal);
                    result.vector(&c.major_dir);
                    result.real(&c.major_sq);
                    result.real(&c.minor_sq);
                }
                Curve::Intersection(c) => {
                    surfaces.extend(c.surfaces);
                    points.push(c.witness);
                }
                _ => {}
            }
        }
        while let Some(point) = points.pop() {
            if !seen_points.insert(point) {
                continue;
            }
            match model.geom().point(point) {
                PointDef::Explicit(p) => result.point(p),
                PointDef::Meet(m) => {
                    match &m.of {
                        MeetOf::Surfaces(s) => surfaces.extend(s),
                        MeetOf::CurveSurface(c, s) => {
                            curves.push(*c);
                            surfaces.push(*s);
                        }
                        MeetOf::Curves(a, b) => curves.extend([a, b]),
                    }
                    match &m.position {
                        Position::Closed(p) => result.point(p),
                        Position::Isolated(root) => {
                            for r in root.root_box().lo().iter().chain(root.root_box().hi()) {
                                result.rational(r);
                            }
                            for p in root.system().equations() {
                                result.rational_polynomial(p);
                            }
                        }
                    }
                }
            }
        }
    }
    result
}
fn growth(label: &str, model: &Model, body: BodyId, dec: &Decider) {
    let size = inventory(model, body);
    println!(
        "growth/{label}: max_precision_bits={} max_observable_rational_bits={} max_numerator_bits={} max_denominator_bits={} opaque_real_coordinates={}",
        dec.stats().max_bits,
        size.rational_bits,
        size.numerator_bits,
        size.denominator_bits,
        size.opaque_real_coordinates
    );
}
fn stage(
    label: &str,
    stage: Stage,
    start: Instant,
    before: &DecisionStats,
    dec: &Decider,
    marks: usize,
) {
    let after = dec.stats();
    let decisions = after.decisions - before.decisions;
    let filter = after.filter - before.filter;
    let facts = after.coincident_fact - before.coincident_fact;
    let denominator = decisions.max(1) as f64;
    println!(
        "stage/{label}/{stage:?}: ns={} decisions={} facts={} facts_share={:.6} filter={} filter_share={:.6} marks={} max_bits={}",
        start.elapsed().as_nanos(),
        decisions,
        facts,
        facts as f64 / denominator,
        filter,
        filter as f64 / denominator,
        marks,
        after.max_bits
    );
}

fn operation(case: &Case, op: BoolOp, placement: &Pose, label: &str) {
    let a = case.a.transformed(&Rational::one(), placement);
    let b = case.b.transformed(&Rational::one(), placement);
    let mut model = Model::new();
    let a_id = a.build(&mut model).unwrap();
    let b_id = b.build(&mut model).unwrap();
    let mut dec = Decider::default();
    let total = Instant::now();
    let label = format!("case{}-{}/{op:?}/{label}", case.number, case.variant);
    macro_rules! timed {
        ($which:expr,$call:expr) => {{
            let before = dec.stats().clone();
            let start = Instant::now();
            match $call {
                Ok(answer) => {
                    stage(&label, $which, start, &before, &dec, answer.marks.len());
                    answer.value
                }
                Err(e) => {
                    println!(
                        "refusal/{label}/{:?}: {} elapsed_ns={}",
                        $which,
                        e.code,
                        total.elapsed().as_nanos()
                    );
                    return;
                }
            }
        }};
    }
    let graph = timed!(
        Stage::Intersect,
        intersect(&mut model, a_id, b_id, &mut dec)
    );
    let imprinted = timed!(Stage::Imprint, imprint(&mut model, &graph, &mut dec));
    let classes = timed!(Stage::Classify, classify(&mut model, &imprinted, &mut dec));
    let before = dec.stats().clone();
    let start = Instant::now();
    let selection = keep(&classes, &imprinted, op);
    stage(&label, Stage::Keep, start, &before, &dec, 0);
    let body = timed!(
        Stage::Sew,
        sew(&mut model, &imprinted, &selection, &mut dec)
    );
    println!(
        "operation/{label}: ns={} faces={} expected_counts={:?}",
        total.elapsed().as_nanos(),
        model.counts(body).faces,
        case.expected[match op {
            BoolOp::Unite => 0,
            BoolOp::Intersect => 1,
            _ => 2,
        }]
    );
    growth(&label, &model, body, &dec);
    let start = Instant::now();
    let mut answered = 0;
    let mut refused = 0;
    let mut mismatches = 0;
    let (al, ah) = a.bounds();
    let (bl, bh) = b.bounds();
    for k in 0..64 {
        let fraction = Rational::from_ratio(k, 63).unwrap();
        let xyz: [Rational; 3] = std::array::from_fn(|i| {
            let lo = al[i].clone().min(bl[i].clone());
            let hi = ah[i].clone().max(bh[i].clone());
            &lo + &((hi - lo.clone()) * fraction.clone())
        });
        let [x, y, z] = xyz;
        let p = RPoint3::new(x, y, z);
        let expected = (a.classify(&p), b.classify(&p));
        let outcome = classify_point(&mut model, body, &p, &mut dec);
        if let catallax_numbers::Outcome::Proven(actual) = outcome {
            answered += 1;
            let (ai, bi) = (
                expected.0 == common::Position::Inside,
                expected.1 == common::Position::Inside,
            );
            if expected.0 != common::Position::Boundary && expected.1 != common::Position::Boundary
            {
                let inside = match op {
                    BoolOp::Unite => ai || bi,
                    BoolOp::Intersect => ai && bi,
                    _ => ai && !bi,
                };
                if (actual == catallax_boolean::Containment::Inside) != inside {
                    mismatches += 1;
                }
            }
        } else if outcome.is_refused() {
            refused += 1;
        }
        let _ = std::hint::black_box(outcome);
    }
    println!(
        "containment/{label}: faces={} ns_per_query={} answered={} refused={} oracle_mismatches={}",
        model.counts(body).faces,
        start.elapsed().as_nanos() / 64,
        answered,
        refused,
        mismatches
    );
}
fn chain(label: &str, base: Primitive, steps: Vec<(Primitive, BoolOp)>) {
    let mut model = Model::new();
    let mut body = base.build(&mut model).unwrap();
    let mut dec = Decider::default();
    growth(&format!("{label}/0"), &model, body, &dec);
    for (i, (operand, op)) in steps.into_iter().enumerate() {
        let other = operand.build(&mut model).unwrap();
        let start = Instant::now();
        match boolean(&mut model, body, other, op, &mut dec) {
            Ok(answer) => {
                body = answer.body;
                println!(
                    "chain/{label}/{}: operation={op:?} ns={} faces={}",
                    i + 1,
                    start.elapsed().as_nanos(),
                    model.counts(body).faces
                );
                growth(&format!("{label}/{}", i + 1), &model, body, &dec);
            }
            Err(e) => {
                println!(
                    "chain-refusal/{label}/{}: stage={:?} code={}",
                    i + 1,
                    e.stage,
                    e.refusal.code
                );
                break;
            }
        }
    }
}
fn chains(selected: Option<u8>) {
    if selected.is_none() || selected == Some(41) {
        chain(
            "case41",
            block_at([10, 10, 10], [0, 0, 0]),
            vec![
                (cylinder_at(2, 12, [5, 5, -1], false), BoolOp::Subtract),
                (cylinder_at(2, 12, [-1, 5, 5], true), BoolOp::Subtract),
            ],
        );
    }
    if selected.is_none() || selected == Some(42) {
        chain(
            "case42",
            block_at([3, 1, 1], [0, 0, 0]),
            vec![
                (block_at([1, 3, 1], [2, 0, 0]), BoolOp::Unite),
                (block_at([3, 1, 1], [0, 2, 0]), BoolOp::Unite),
                (block_at([1, 3, 1], [0, 0, 0]), BoolOp::Unite),
            ],
        );
    }
    if selected.is_none() || selected == Some(43) {
        chain(
            "case43",
            ball_at(5, [0, 0, 0]),
            vec![
                (ball_at(5, [6, 0, 0]), BoolOp::Subtract),
                (ball_at(5, [-6, 0, 0]), BoolOp::Subtract),
            ],
        );
    }
    if selected.is_none() || selected == Some(44) {
        chain(
            "case44",
            block_at([2, 1, 1], [0, 0, 0]),
            (1..=5)
                .map(|k| (block_at([2, 1, 1], [k, 0, 0]), BoolOp::Unite))
                .collect(),
        );
    }
    if selected.is_none() || selected == Some(0) {
        let mut seed = 0x6581_4027u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed % 201) as i64 - 100
        };
        let operands = [
            cylinder_at(3, 12, [5, 5, -1], false),
            ball_at(3, [8, 5, 5]),
            block_at([4, 4, 4], [8, 4, 3]),
            cylinder_at(2, 16, [7, 7, -3], false),
            ball_at(4, [8, 8, 8]),
        ];
        let steps = operands
            .into_iter()
            .enumerate()
            .map(|(i, p)| {
                let rotation =
                    Rotation::from_quaternion(r(101), r(next()), r(next()), r(next())).unwrap();
                let pose = Pose::new(
                    rotation,
                    (RVector3::new(r(next() % 5), r(next() % 5), r(next() % 5))).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap();
                (
                    p.transformed(&Rational::one(), &pose),
                    if i % 2 == 0 {
                        BoolOp::Unite
                    } else {
                        BoolOp::Subtract
                    },
                )
            })
            .collect();
        chain(
            "alternating-rational-poses",
            block_at([10, 10, 10], [0, 0, 0]),
            steps,
        );
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let selected = args
        .windows(2)
        .find(|w| w[0] == "--case")
        .map(|w| w[1].parse::<u8>().expect("--case needs a case number"));
    let placed = Pose::new(
        Rotation::from_quaternion(r(17), r(2), r(-3), r(5)).unwrap(),
        (RVector3::new(r(1000000), r(-1000000), r(300000))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    for case in common::cases() {
        if case.number > 40 || selected.is_some_and(|n| n != case.number) {
            continue;
        }
        for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
            operation(&case, op, &Pose::identity(), "identity");
            operation(&case, op, &placed, "rational-pose");
        }
    }
    chains(selected);
}
