#![allow(dead_code, unused_imports, unused_macros)]
pub use catallax_numbers::{
    Answer, Ball, BoundKind, BudgetPolicy, Decider, Fact, FactKind, FactOrigin, Mark, MarkLog,
    Rational,
};
pub use catallax_shapes::*;
pub use catallax_step::{StepOptions, StepReport, write_step};
pub use catallax_topology::*;
pub use std::collections::{BTreeMap, BTreeSet};
pub mod budget;
pub mod fixtures;
pub mod reader;
pub mod spline;
pub mod topology;
pub use budget::*;
pub use fixtures::*;
pub use reader::*;
pub use spline::*;
pub use topology::*;

macro_rules! blind_test {
    ($(#[$attr:meta])* $name:ident, $body:block) => {
        #[test]
        $(#[$attr])*
        fn $name() $body
    };
}
pub(crate) use blind_test;
pub fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
pub fn n(i: i64) -> Rational {
    Rational::from_i64(i)
}
pub fn xyz(p: &RPoint3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
pub fn vec_xyz(p: &RVector3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
pub fn real_xyz(p: &Point3) -> [Rational; 3] {
    p.rational_components()
        .expect("rational STEP fixture")
        .map(Clone::clone)
}
pub fn real_vec_xyz(p: &Vector3) -> [Rational; 3] {
    p.rational_components()
        .expect("rational STEP fixture")
        .map(Clone::clone)
}
pub fn point(p: [Rational; 3]) -> RPoint3 {
    let [x, y, z] = p;
    RPoint3::new(x, y, z)
}
pub fn dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    (0..3).fold(n(0), |s, i| s + &a[i] * &b[i])
}
pub fn sub(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
pub fn sq(q: &Rational) -> Rational {
    q * q
}
pub fn distance_sq(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    let v = sub(a, b);
    dot(&v, &v)
}
pub fn options(digits: u32) -> StepOptions {
    StepOptions {
        product_name: "catallax".into(),
        curve_tolerance: r("1e-7"),
        significant_digits: digits,
    }
}
pub fn export(m: &Model, bs: &[BodyId], o: &StepOptions) -> (Vec<u8>, Answer<StepReport>, Step) {
    // Every default-budget fixture, including all generated cases 13..15,
    // also exercises the architect's zero-work variant. Default refusal is
    // still a failure; the zero-work run alone may refuse atomically.
    let result = export_with(m, bs, o, &mut Decider::default());
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    });
    let before = b"zero-work output prefix".to_vec();
    let mut out = before.clone();
    match write_step(m, bs, o, &mut dec, &mut out) {
        Ok(answer) => {
            assert_eq!(&out[..before.len()], before);
            check_export(m, bs, o, out.split_off(before.len()), answer);
        }
        Err(refusal) => assert_budget_refusal(refusal.code, &before, &out),
    }
    result
}
pub fn export_with(
    m: &Model,
    bs: &[BodyId],
    o: &StepOptions,
    dec: &mut Decider,
) -> (Vec<u8>, Answer<StepReport>, Step) {
    for &b in bs {
        assert_eq!(m.check(b), vec![], "invalid test fixture");
    }
    let mut bytes = Vec::new();
    let answer = write_step(m, bs, o, dec, &mut bytes).expect("export must answer");
    check_export(m, bs, o, bytes, answer)
}
pub fn check_export(
    m: &Model,
    bs: &[BodyId],
    o: &StepOptions,
    bytes: Vec<u8>,
    answer: Answer<StepReport>,
) -> (Vec<u8>, Answer<StepReport>, Step) {
    let step = Step::parse(&bytes, &o.product_name);
    let report = &answer.value;
    assert!(report.max_deviation >= n(0));
    assert_eq!(report.solids, step.items().len() - 1);
    assert_eq!(report.faces, step.count("ADVANCED_FACE"));
    assert_eq!(report.edges, step.count("EDGE_CURVE"));
    assert_eq!(report.vertices, step.count("VERTEX_POINT"));
    // Bounded chart pcurves are 2D splines inside a definitional
    // representation. Only 3D splines are approximated space curves.
    let pcurves: BTreeSet<usize> = step
        .ids("DEFINITIONAL_REPRESENTATION")
        .into_iter()
        .flat_map(|id| list(&step.entities[&id].args[1]).into_iter().map(reference))
        .collect();
    assert_eq!(
        report.approximated_curves,
        step.ids("B_SPLINE_CURVE_WITH_KNOTS")
            .into_iter()
            .filter(|id| !pcurves.contains(id))
            .count()
    );
    assert_eq!(
        answer.marks.len(),
        usize::from(report.max_deviation > o.curve_tolerance)
    );
    for mark in answer.marks.iter() {
        assert_eq!(mark.code, "step.tolerance-not-reached");
        assert_eq!(mark.kind, BoundKind::Distance);
        assert_eq!(mark.bound, report.max_deviation);
    }
    let uncertainty = report.max_deviation.clone().max(o.curve_tolerance.clone());
    step.fixed(
        &o.product_name,
        &decimal(&uncertainty, o.significant_digits, true),
    );
    assert_topology(m, bs, &step, o, &report.max_deviation);
    (bytes, answer, step)
}
pub fn counts(s: &Step, expected: &[(&str, usize)]) {
    for &(kind, count) in expected {
        assert_eq!(s.count(kind), count, "{kind}");
    }
}
pub fn report(a: &Answer<StepReport>, values: [usize; 7], exact: bool) {
    let p = &a.value;
    assert_eq!(
        [
            p.solids,
            p.faces,
            p.edges,
            p.vertices,
            p.added_vertices,
            p.added_edges,
            p.approximated_curves
        ],
        values
    );
    if exact {
        assert_eq!(p.max_deviation, n(0));
        assert!(a.marks.is_empty());
    }
}
pub fn vertex_set(s: &Step) -> BTreeSet<[Rational; 3]> {
    s.ids("VERTEX_POINT")
        .iter()
        .map(|&id| s.xyz(s.vertex(id)))
        .collect()
}
// |sqrt(q)-radius| <= d, squared without an irrational oracle or f64.
pub fn radial_distance(q: Rational, radius: &Rational, d: &Rational) {
    assert!(q >= n(0));
    assert!(q <= sq(&(radius + d)), "radial upper distance");
    if radius > d {
        assert!(q >= sq(&(radius - d)), "radial lower distance");
    }
}
pub fn on_written_surface(s: &Step, surface: usize, p: &[Rational; 3], d: &Rational) {
    let entity = &s.entities[&surface];
    let (origin, axis, _) = s.placement(reference(&entity.args[1]));
    let v = sub(p, &origin);
    let q = dot(&v, &v);
    let a2 = dot(&axis, &axis);
    assert!(a2 > n(0));
    match entity.kind.as_str() {
        "PLANE" => assert!(
            sq(&dot(&v, &axis)) <= sq(d) * a2,
            "vertex outside plane deviation"
        ),
        "SPHERICAL_SURFACE" => radial_distance(q, &number(&entity.args[2]), d),
        "CYLINDRICAL_SURFACE" => {
            radial_distance(q - sq(&dot(&v, &axis)) / a2, &number(&entity.args[2]), d)
        }
        "CONICAL_SURFACE" | "TOROIDAL_SURFACE" => {
            let axial2 = sq(&dot(&v, &axis)) / a2.clone();
            let radial2 = &q - &axial2;
            let ball = |v: &Rational| Ball::from_rational(v, 512);
            let rho = ball(&radial2).sqrt().unwrap();
            let distance = if entity.kind == "TOROIDAL_SURFACE" {
                let major = number(&entity.args[2]);
                let minor = number(&entity.args[3]);
                let residual = &q + &sq(&major) - sq(&minor);
                // On a ring torus, the squared implicit equation is exact and
                // residual=2 R rho is nonnegative. This also proves d=0 cases.
                if residual >= n(0) && sq(&residual) == n(4) * sq(&major) * radial2 {
                    return;
                }
                let dr = &rho - &ball(&major);
                (&dr * &dr + ball(&axial2)).sqrt().unwrap() - ball(&minor)
            } else {
                if q.is_zero() {
                    return;
                }
                let h = ball(&dot(&v, &axis))
                    .checked_div(&ball(&a2).sqrt().unwrap())
                    .unwrap();
                let angle = ball(&number(&entity.args[3]));
                let sin = angle.sin();
                let cos = angle.cos();
                let along = &rho * &sin + &h * &cos;
                if along.hi() < n(0) {
                    ball(&q).sqrt().unwrap()
                } else {
                    &rho * &cos - &h * &sin
                }
            };
            assert!(
                distance.lo().abs().max(distance.hi().abs()) <= *d,
                "vertex outside curved surface deviation"
            );
        }
        k => panic!("unexpected surface {k}"),
    }
}
// Quaternion rotation is computed here from its formula, independently of Pose.
pub fn rotation_oracle(q: [i64; 4]) -> [[Rational; 3]; 3] {
    let [w, x, y, z] = q.map(n);
    let den = &w * &w + &x * &x + &y * &y + &z * &z;
    [
        [
            &w * &w + &x * &x - &y * &y - &z * &z,
            n(2) * (&x * &y - &w * &z),
            n(2) * (&x * &z + &w * &y),
        ],
        [
            n(2) * (&x * &y + &w * &z),
            &w * &w - &x * &x + &y * &y - &z * &z,
            n(2) * (&y * &z - &w * &x),
        ],
        [
            n(2) * (&x * &z - &w * &y),
            n(2) * (&y * &z + &w * &x),
            &w * &w - &x * &x - &y * &y + &z * &z,
        ],
    ]
    .map(|row| row.map(|v| &v / &den))
}
pub fn apply(matrix: &[[Rational; 3]; 3], t: &[Rational; 3], p: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| &dot(&matrix[i], p) + &t[i])
}
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn integer(&mut self, span: u64) -> i64 {
        (self.next() % (2 * span + 1)) as i64 - span as i64
    }
    pub fn positive(&mut self) -> Rational {
        Rational::from_ratio(
            (self.next() % 100 + 1) as i64,
            (self.next() % 29 + 1) as i64,
        )
        .unwrap()
    }
    pub fn pose(&mut self) -> (Pose, [[Rational; 3]; 3], [Rational; 3]) {
        let q = [
            (self.next() % 100 + 1) as i64,
            self.integer(100),
            self.integer(100),
            self.integer(100),
        ];
        let t = std::array::from_fn(|_| {
            Rational::from_ratio(self.integer(1_000_000), (self.next() % 29 + 1) as i64).unwrap()
        });
        let [w, x, y, z] = q.map(n);
        (
            Pose::new(
                Rotation::from_quaternion(w, x, y, z).unwrap(),
                (RVector3::new(t[0].clone(), t[1].clone(), t[2].clone())).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
            rotation_oracle(q),
            t,
        )
    }
}
