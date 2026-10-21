//! Opt-in measurements. No elapsed-time assertion is a correctness condition.
#[path = "../tests/common/mod.rs"]
mod common;
use common::*;
use std::time::Instant;

fn report(label: &str, stage: &str, start: Instant, stats: &DecisionStats, marks: usize) {
    let denominator = stats.decisions.max(1) as f64;
    println!(
        "{label},stage={stage},ns={},decisions={},facts={},filter={},fact_share={},filter_share={},work={},max_bits={},marks={marks}",
        start.elapsed().as_nanos(),
        stats.decisions,
        stats.coincident_fact,
        stats.filter,
        stats.coincident_fact as f64 / denominator,
        stats.filter as f64 / denominator,
        stats.work,
        stats.max_bits
    );
}
fn stage_benchmark(c: &Case, pose: &Pose, pose_name: &str, op: BoolOp) {
    let label = format!("case={},pose={pose_name},op={op:?}", c.id);
    let mut m = Model::new();
    let a = c.a.transformed(&r(1), pose).build(&mut m);
    let b = c.b.transformed(&r(1), pose).build(&mut m);
    let mut dec = Decider::default();
    macro_rules! stage {
        ($name:literal,$expr:expr) => {{
            let before = dec.stats().clone();
            let start = Instant::now();
            let result = $expr;
            let stats = stats_delta(&before, dec.stats());
            match result {
                Ok(answer) => {
                    report(&label, $name, start, &stats, answer.marks.len());
                    answer.value
                }
                Err(e) => {
                    report(&label, $name, start, &stats, 0);
                    println!("{label},refused_stage={},code={}", $name, e.code);
                    return;
                }
            }
        }};
    }
    let graph = stage!("intersect", intersect(&mut m, a, b, &mut dec));
    let imprinted = stage!("imprint", imprint(&mut m, &graph, &mut dec));
    let classification = stage!("classify", classify(&mut m, &imprinted, &mut dec));
    let before = dec.stats().clone();
    let start = Instant::now();
    let selection = keep(&classification, &imprinted, op);
    report(&label, "keep", start, &stats_delta(&before, dec.stats()), 0);
    let body = stage!("sew", sew(&mut m, &imprinted, &selection, &mut dec));
    let faces = m.counts(body).faces;
    let point = pose.apply_point(&p("1/7", "2/7", "3/7"));
    let mut outcomes = [0usize; 4];
    let before = dec.stats().clone();
    let start = Instant::now();
    for _ in 0..32 {
        let index = match std::hint::black_box(classify_point(&mut m, body, &point, &mut dec)) {
            Outcome::Proven(_) => 0,
            Outcome::ProvenCoincident(_) => 1,
            Outcome::NotExact(_) => 2,
            Outcome::Refused(_) => 3,
        };
        outcomes[index] += 1;
    }
    println!(
        "{label},query=classify_point,faces={faces},samples=32,outcomes_proven_coincident_notexact_refused={outcomes:?},total_ns={},stats={:?}",
        start.elapsed().as_nanos(),
        stats_delta(&before, dec.stats())
    );
}

// Public stored-definition metric. No private representation or cached f64 is
// read. Opaque non-rational Real expression trees are counted separately.
#[derive(Default)]
struct Size {
    largest_numerator_bits: u64,
    largest_denominator_bits: u64,
    opaque_reals: usize,
    seen: BTreeSet<GeomKey>,
}
impl Size {
    fn rational(&mut self, q: &Rational) {
        self.largest_numerator_bits = self.largest_numerator_bits.max(q.numerator().bit_length());
        self.largest_denominator_bits = self
            .largest_denominator_bits
            .max(q.denominator().bit_length());
    }
    fn point(&mut self, p: &RPoint3) {
        for q in [&p.x, &p.y, &p.z] {
            self.rational(q);
        }
    }
    fn vector(&mut self, p: &RVector3) {
        for q in [&p.x, &p.y, &p.z] {
            self.rational(q);
        }
    }
    fn real_point(&mut self, p: &Point3) {
        for x in [&p.x, &p.y, &p.z] {
            if let Some(q) = x.as_rational() {
                self.rational(&q);
            } else {
                self.opaque_reals += 1;
            }
        }
    }
    fn surface(&mut self, m: &Model, id: SurfaceId) {
        if !self.seen.insert(GeomKey::Surface(id)) {
            return;
        }
        match m.geom().surface(id) {
            Surface::Plane(s) => {
                self.point(&s.origin);
                self.vector(&s.normal);
                self.vector(&s.ref_dir);
            }
            Surface::Cylinder(s) => {
                self.point(&s.origin);
                self.vector(&s.axis);
                self.vector(&s.ref_dir);
                self.rational(&s.radius);
            }
            Surface::Sphere(s) => {
                self.point(&s.center);
                self.vector(&s.axis);
                self.vector(&s.ref_dir);
                self.rational(&s.radius);
            }
            _ => panic!("benchmark uses only supported surfaces"),
        }
        for (_, q) in m.geom().implicit(id).terms() {
            self.rational(&q);
        }
    }
    fn curve(&mut self, m: &Model, id: CurveId) {
        if !self.seen.insert(GeomKey::Curve(id)) {
            return;
        }
        match m.geom().curve(id) {
            Curve::Line(c) => {
                self.point(&c.origin);
                self.vector(&c.dir);
            }
            Curve::Circle(c) => {
                self.point(&c.center);
                self.vector(&c.normal);
                self.vector(&c.ref_dir);
                self.rational(&c.radius_sq);
            }
            Curve::Ellipse(c) => {
                self.point(&c.center);
                self.vector(&c.normal);
                self.vector(&c.major_dir);
                self.rational(&c.major_sq);
                self.rational(&c.minor_sq);
            }
            Curve::Intersection(c) => {
                for &s in &c.surfaces {
                    self.surface(m, s);
                }
                self.stored_point(m, c.witness);
            }
            _ => panic!("new curve definition"),
        }
        for equation in m.geom().curve(id).equations(m.geom()) {
            for (_, q) in equation.terms() {
                self.rational(&q);
            }
        }
    }
    fn stored_point(&mut self, m: &Model, id: PointId) {
        if !self.seen.insert(GeomKey::Point(id)) {
            return;
        }
        match m.geom().point(id) {
            PointDef::Explicit(p) => self.point(p),
            PointDef::Meet(meet) => {
                match &meet.of {
                    MeetOf::Surfaces(ss) => {
                        for &s in ss {
                            self.surface(m, s);
                        }
                    }
                    MeetOf::CurveSurface(c, s) => {
                        self.curve(m, *c);
                        self.surface(m, *s);
                    }
                    MeetOf::Curves(a, b) => {
                        self.curve(m, *a);
                        self.curve(m, *b);
                    }
                }
                match &meet.position {
                    Position::Closed(p) => self.real_point(p),
                    Position::Isolated(root) => {
                        for q in root.root_box().lo().iter().chain(root.root_box().hi()) {
                            self.rational(q);
                        }
                        for polynomial in root.system().equations() {
                            for (_, q) in polynomial.terms() {
                                self.rational(&q);
                            }
                        }
                    }
                }
            }
        }
    }
    fn body(m: &Model, b: BodyId) -> Self {
        let mut size = Self::default();
        for f in m.faces_of_body(b) {
            size.surface(m, m.face(f).surface);
        }
        for e in m.edges_of_body(b) {
            size.curve(m, m.edge(e).curve);
        }
        for v in m.vertices_of_body(b) {
            size.stored_point(m, m.vertex(v).point);
        }
        for mark in m.body(b).marks.iter() {
            size.rational(&mark.bound);
        }
        size
    }
}
fn growth(label: &str, m: &Model, b: BodyId, step: usize, stats: &DecisionStats) {
    let size = Size::body(m, b);
    println!(
        "chain={label},step={step},max_precision_bits={},largest_observable_stored_rational_numerator_bits={},largest_observable_stored_rational_denominator_bits={},opaque_reals={}",
        stats.max_bits,
        size.largest_numerator_bits,
        size.largest_denominator_bits,
        size.opaque_reals
    );
}
fn chains() {
    for n in 41..=44 {
        let spec = chain(n);
        let mut m = Model::new();
        let mut bodies: Vec<_> = spec.operands.iter().map(|p| p.build(&mut m)).collect();
        let mut dec = Decider::default();
        for (step, (a, b, op, _)) in spec.steps.iter().enumerate() {
            match boolean(&mut m, bodies[*a], bodies[*b], *op, &mut dec) {
                Ok(answer) => {
                    growth(&n.to_string(), &m, answer.body, step + 1, dec.stats());
                    bodies.push(answer.body);
                }
                Err(e) => {
                    println!("chain={n},step={},refusal={e:?}", step + 1);
                    break;
                }
            }
        }
    }
    let mut rng = Rng(0x5a17e);
    let mut m = Model::new();
    let mut current = block_at([10; 3], [0; 3]).build(&mut m);
    let mut dec = Decider::default();
    for step in 0..5 {
        // Small translations keep the mixed operands in the same neighbourhood;
        // quaternion components are the same fixed-seed rational pose family.
        let mut pose = rng.pose();
        pose.translation = pose.translation.scaled(&q("1e-6"));
        let primitive = match step % 3 {
            0 => cyl_at(3, 12, [5, 5, -1], false),
            1 => ball_at(4, [5; 3]),
            _ => block_at([6; 3], [2; 3]),
        }
        .transformed(&r(1), &pose);
        let next = primitive.build(&mut m);
        let op = if step % 2 == 0 {
            BoolOp::Unite
        } else {
            BoolOp::Subtract
        };
        match boolean(&mut m, current, next, op, &mut dec) {
            Ok(answer) => {
                current = answer.body;
                growth("alternating-mixed", &m, current, step + 1, dec.stats());
            }
            Err(e) => {
                println!("chain=alternating-mixed,step={},refusal={e:?}", step + 1);
                break;
            }
        }
    }
}
fn main() {
    if std::env::var_os("CATALLAX_BOOLEAN_RUN_BENCH").as_deref() != Some(std::ffi::OsStr::new("1"))
    {
        println!(
            "Enable after pieces 4 and 5 are implemented: CATALLAX_BOOLEAN_RUN_BENCH=1 cargo bench -p catallax-boolean --bench speed"
        );
        return;
    }
    println!("boolean benchmark: time is measured, never gated");
    let moved = Rng(0x5beac4).pose();
    for &id in IDS {
        let c = case(id);
        for (pose_name, pose) in [("identity", Pose::identity()), ("rational", moved.clone())] {
            for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
                let mut m = Model::new();
                let a = c.a.transformed(&r(1), &pose).build(&mut m);
                let b = c.b.transformed(&r(1), &pose).build(&mut m);
                let mut dec = Decider::default();
                let start = Instant::now();
                let result = boolean(&mut m, a, b, op, &mut dec);
                report(
                    &format!("case={id},pose={pose_name},op={op:?}"),
                    "whole",
                    start,
                    dec.stats(),
                    result.as_ref().map_or(0, |a| a.marks.len()),
                );
                if let Err(e) = result {
                    println!("case={id},whole_refusal={e:?}");
                }
                stage_benchmark(&c, &pose, pose_name, op);
            }
        }
    }
    chains();
}
