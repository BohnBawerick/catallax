#![allow(dead_code, unused_macros, unused_imports)]

use catallax_numbers::*;
use catallax_shapes::*;
pub use rug::Rational as Q;
use std::sync::Arc;

pub mod cases;

macro_rules! blind_test {
    ($name:ident, $body:block) => {
        #[test]
        fn $name() $body
    };
}
pub(crate) use blind_test;

pub fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
pub fn real(s: &str) -> Real {
    Real::from_rational(r(s))
}
pub fn p(x: &str, y: &str, z: &str) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
pub fn v(x: &str, y: &str, z: &str) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
pub fn rp(a: [Rational; 3]) -> RPoint3 {
    let [x, y, z] = a;
    RPoint3::new(x, y, z)
}
pub fn xyz(p: &RPoint3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
pub fn reals(p: &Point3) -> [Real; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
pub fn vx(v: &RVector3) -> [Rational; 3] {
    [v.x.clone(), v.y.clone(), v.z.clone()]
}
pub fn q(n: i64, d: i64) -> Q {
    Q::from((n, d))
}
pub fn oq(x: &Rational) -> Q {
    Q::from((
        x.numerator().to_string().parse::<rug::Integer>().unwrap(),
        x.denominator().to_string().parse::<rug::Integer>().unwrap(),
    ))
}
pub fn eq(x: &Q) -> Rational {
    r(&x.to_string())
}
pub fn qp(p: &RPoint3) -> [Q; 3] {
    xyz(p).map(|x| oq(&x))
}
pub fn qv(v: &RVector3) -> [Q; 3] {
    vx(v).map(|x| oq(&x))
}
pub fn dot(a: &[Q; 3], b: &[Q; 3]) -> Q {
    (0..3).map(|i| a[i].clone() * &b[i]).sum()
}
pub fn cross(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| {
        a[(i + 1) % 3].clone() * &b[(i + 2) % 3] - a[(i + 2) % 3].clone() * &b[(i + 1) % 3]
    })
}
pub fn sub(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| a[i].clone() - &b[i])
}
pub fn sq(x: &Q) -> Q {
    x.clone() * x
}
pub fn sign(x: &Q) -> i8 {
    if x > &0 {
        1
    } else if x < &0 {
        -1
    } else {
        0
    }
}

// Table formulas, evaluated in GMP rationals. No catallax-shapes arithmetic is used.
pub fn implicit_oracle(s: &Surface, point: &RPoint3) -> Q {
    let point = qp(point);
    match s {
        Surface::Plane(s) => dot(
            &sub(
                &point,
                &qp(&s.origin.as_rational().expect("rational oracle fixture")),
            ),
            &qv(&s.normal.as_rational().expect("rational oracle fixture")),
        ),
        Surface::Cylinder(s) => {
            let d = sub(
                &point,
                &qp(&s.origin.as_rational().expect("rational oracle fixture")),
            );
            let a = qv(&s.axis.as_rational().expect("rational oracle fixture"));
            let c = cross(&d, &a);
            dot(&c, &c)
                - sq(&oq(&s
                    .radius
                    .as_rational()
                    .expect("rational oracle fixture")))
                    * dot(&a, &a)
        }
        Surface::Cone(s) => {
            let d = sub(
                &point,
                &qp(&s.apex.as_rational().expect("rational oracle fixture")),
            );
            let a = qv(&s.axis.as_rational().expect("rational oracle fixture"));
            let c = cross(&d, &a);
            dot(&c, &c)
                - sq(&oq(&s
                    .slope
                    .as_rational()
                    .expect("rational oracle fixture")))
                    * sq(&dot(&d, &a))
        }
        Surface::Sphere(s) => {
            let d = sub(
                &point,
                &qp(&s.center.as_rational().expect("rational oracle fixture")),
            );
            dot(&d, &d)
                - sq(&oq(&s
                    .radius
                    .as_rational()
                    .expect("rational oracle fixture")))
        }
        Surface::Torus(s) => {
            let d = sub(
                &point,
                &qp(&s.center.as_rational().expect("rational oracle fixture")),
            );
            let a = qv(&s.axis.as_rational().expect("rational oracle fixture"));
            let norm = dot(&d, &d);
            let major = sq(&oq(&s
                .major
                .as_rational()
                .expect("rational oracle fixture")));
            let minor = sq(&oq(&s
                .minor
                .as_rational()
                .expect("rational oracle fixture")));
            sq(&(norm.clone() + &major - minor))
                - q(4, 1) * major * (norm - sq(&dot(&d, &a)) / dot(&a, &a))
        }
        _ => panic!("add an oracle for the new surface kind"),
    }
}

pub fn plane() -> Surface {
    Surface::Plane(
        Plane::new(
            (p("0", "0", "1")).into(),
            (v("0", "0", "2")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn cylinder(radius: &str) -> Surface {
    Surface::Cylinder(
        Cylinder::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r(radius)).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn cone() -> Surface {
    Surface::Cone(
        Cone::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r("1/2")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn sphere(radius: &str) -> Surface {
    Surface::Sphere(
        Sphere::new(
            (p("0", "0", "0")).into(),
            (r(radius)).into(),
            (v("0", "0", "1")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn torus() -> Surface {
    Surface::Torus(
        Torus::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r("3")).into(),
            (r("1")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn surfaces() -> Vec<Surface> {
    vec![plane(), cylinder("5"), cone(), sphere("3"), torus()]
}
pub fn line() -> Curve {
    Curve::Line(
        Line::new(
            (p("1", "2", "3")).into(),
            (v("2", "-1", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn circle(radius_sq: &str) -> Curve {
    Curve::Circle(
        Circle::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r(radius_sq)).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn ellipse() -> Curve {
    Curve::Ellipse(
        Ellipse::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (v("2", "0", "0")).into(),
            (r("25")).into(),
            (r("9")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
pub fn curves() -> Vec<Curve> {
    vec![line(), circle("25"), ellipse()]
}
pub fn rotation(a: [&str; 4]) -> Rotation {
    let [w, x, y, z] = a.map(r);
    Rotation::from_quaternion(w, x, y, z).unwrap()
}
pub fn poly(terms: &[([u32; 3], &str)]) -> MPoly {
    MPoly::from_terms(3, terms.iter().map(|(e, c)| (e.to_vec(), r(c))).collect())
}
pub fn ball(x: &str, bits: u32) -> Ball {
    Ball::from_rational(&r(x), bits)
}
pub fn balls(p: &RPoint3, bits: u32) -> [Ball; 3] {
    xyz(p).map(|x| Ball::from_rational(&x, bits))
}
pub fn contains(a: &[Ball; 3], expected: &RPoint3) {
    for (b, x) in a.iter().zip(xyz(expected)) {
        assert!(b.contains(&x), "{b:?} excludes {x}");
    }
}
pub fn overlap(a: &Ball, b: &Ball) {
    assert!(
        a.lo() <= b.hi() && b.lo() <= a.hi(),
        "disjoint {a:?}, {b:?}"
    );
}
pub fn real_eq(a: &Real, b: &Real) {
    assert!(
        Decider::default().compare(a, b).is_coincident(),
        "unequal exact reals {a:?}, {b:?}"
    );
}
pub fn point_eq(a: &Point3, b: &Point3) {
    for (a, b) in reals(a).iter().zip(reals(b)) {
        real_eq(a, &b);
    }
}
pub fn vector_eq(a: &Vector3, b: &RVector3) {
    for (a, b) in [&a.x, &a.y, &a.z].into_iter().zip(vx(b)) {
        real_eq(a, &Real::from_rational(b));
    }
}
pub fn coincident<T: std::fmt::Debug>(out: Outcome<T>, proof: Option<ZeroProof>) {
    match out {
        Outcome::ProvenCoincident(c) => {
            if let Some(proof) = proof {
                assert_eq!(c.proof, proof);
            }
        }
        other => panic!("expected coincidence, got {other:?}"),
    }
}
pub fn side_answer(out: Outcome<Side>, expected: i8) {
    match expected {
        0 => coincident(out, None),
        1 => assert_eq!(out, Outcome::Proven(Side::Front)),
        -1 => assert_eq!(out, Outcome::Proven(Side::Back)),
        _ => unreachable!(),
    }
}
pub fn incidence(out: Outcome<Apart>, on: bool) {
    if on {
        coincident(out, None)
    } else {
        assert_eq!(out, Outcome::Proven(Apart))
    }
}
pub fn refusal<T: std::fmt::Debug>(out: Outcome<T>, code: &str) {
    match out {
        Outcome::Refused(e) => assert_eq!(e.code, code),
        other => panic!("expected {code}, got {other:?}"),
    }
}
pub fn rejected<T: std::fmt::Debug>(out: Result<T, Refusal>, code: &str) {
    assert_eq!(out.unwrap_err().code, code);
}
pub fn tiny() -> Decider {
    Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    })
}
pub fn expired() -> Decider {
    let mut d = Decider::default();
    d.set_deadline(Some(
        std::time::Instant::now() - std::time::Duration::from_secs(1),
    ));
    d
}
pub fn explicit(store: &mut GeomStore, point: RPoint3) -> PointId {
    store.add_point(PointDef::Explicit((point).into()))
}
pub fn closed_on(store: &mut GeomStore, point: Point3, supports: Vec<SurfaceId>) -> PointId {
    store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(supports),
        position: Position::Closed(point),
    }))
}
pub fn fact(kind: FactKind, a: GeomKey, b: GeomKey) -> Fact<GeomKey> {
    Fact {
        kind,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("test.exact-construction"),
    }
}

pub fn sphere_root(positive: bool) -> IsolatedRoot {
    // Subtraction gives 12x=36 and 12y=36; z^2=25-18=7.
    let base = poly(&[
        ([2, 0, 0], "1"),
        ([0, 2, 0], "1"),
        ([0, 0, 2], "1"),
        ([0, 0, 0], "-25"),
    ]);
    let sx = &base + &poly(&[([1, 0, 0], "-12"), ([0, 0, 0], "36")]);
    let sy = &base + &poly(&[([0, 1, 0], "-12"), ([0, 0, 0], "36")]);
    let system = Arc::new(PolySystem::new(vec![base, sx, sy]).unwrap());
    let box_ = RootBox::new(
        vec![r("2"), r("2"), r(if positive { "2" } else { "-3" })],
        vec![r("4"), r("4"), r(if positive { "3" } else { "-2" })],
    )
    .unwrap();
    match Decider::default().certify_root(&system, &box_) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        other => panic!("piece 1 root prerequisite: {other:?}"),
    }
}
pub fn circle_root() -> IsolatedRoot {
    // x=y>0, z=0, 2x^2=1 isolates the positive sqrt(1/2).
    let system = Arc::new(
        PolySystem::new(vec![
            poly(&[
                ([2, 0, 0], "1"),
                ([0, 2, 0], "1"),
                ([0, 0, 2], "1"),
                ([0, 0, 0], "-1"),
            ]),
            poly(&[([1, 0, 0], "1"), ([0, 1, 0], "-1")]),
            poly(&[([0, 0, 1], "1")]),
        ])
        .unwrap(),
    );
    let box_ = RootBox::new(
        vec![r("1/2"), r("1/2"), r("-1/2")],
        vec![r("1"), r("1"), r("1/2")],
    )
    .unwrap();
    match Decider::default().certify_root(&system, &box_) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        other => panic!("piece 1 root prerequisite: {other:?}"),
    }
}

pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    pub fn small(&mut self) -> i64 {
        (self.next() % 201) as i64 - 100
    }
    pub fn rational(&mut self) -> Rational {
        Rational::from_ratio(self.small(), (self.next() % 31 + 1) as i64).unwrap()
    }
    pub fn pose(&mut self) -> Pose {
        let w = Rational::from_i64((self.next() % 100 + 1) as i64);
        let [x, y, z] = std::array::from_fn(|_| Rational::from_i64(self.small()));
        let rotation = Rotation::from_quaternion(w, x, y, z).unwrap();
        let [x, y, z] = std::array::from_fn(|_| {
            Rational::from_i64((self.next() % 2_000_001) as i64 - 1_000_000)
        });
        Pose::new(
            rotation,
            (RVector3::new(x, y, z)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap()
    }
}

pub fn quaternion_matrix(a: &[Q; 4]) -> [[Q; 3]; 3] {
    let [w, x, y, z] = a;
    let n = sq(w) + sq(x) + sq(y) + sq(z);
    let twice = |x: Q| q(2, 1) * x;
    [
        [
            (sq(w) + sq(x) - sq(y) - sq(z)) / &n,
            twice(x.clone() * y - w.clone() * z) / &n,
            twice(x.clone() * z + w.clone() * y) / &n,
        ],
        [
            twice(x.clone() * y + w.clone() * z) / &n,
            (sq(w) - sq(x) + sq(y) - sq(z)) / &n,
            twice(y.clone() * z - w.clone() * x) / &n,
        ],
        [
            twice(x.clone() * z - w.clone() * y) / &n,
            twice(y.clone() * z + w.clone() * x) / &n,
            (sq(w) - sq(x) - sq(y) + sq(z)) / &n,
        ],
    ]
}
pub fn hamilton(a: &[Q; 4], b: &[Q; 4]) -> [Q; 4] {
    let av = [a[1].clone(), a[2].clone(), a[3].clone()];
    let bv = [b[1].clone(), b[2].clone(), b[3].clone()];
    let c = cross(&av, &bv);
    [
        a[0].clone() * &b[0] - dot(&av, &bv),
        a[0].clone() * &b[1] + b[0].clone() * &a[1] + &c[0],
        a[0].clone() * &b[2] + b[0].clone() * &a[2] + &c[1],
        a[0].clone() * &b[3] + b[0].clone() * &a[3] + &c[2],
    ]
}

pub fn scale_point(p: &RPoint3, s: &Rational) -> RPoint3 {
    rp(xyz(p).map(|x| &x * s))
}
fn scale_definition_point(p: &Point3, s: &Real) -> Point3 {
    Point3::new(&p.x * s, &p.y * s, &p.z * s)
}
pub fn scale_surface(surface: &Surface, s: &Rational) -> Surface {
    let scale = Real::from(s.clone());
    match surface {
        Surface::Plane(a) => Surface::Plane(Plane {
            origin: scale_definition_point(&a.origin, &scale),
            ..a.clone()
        }),
        Surface::Cylinder(a) => Surface::Cylinder(Cylinder {
            origin: scale_definition_point(&a.origin, &scale),
            radius: &a.radius * &scale,
            ..a.clone()
        }),
        Surface::Cone(a) => Surface::Cone(Cone {
            apex: scale_definition_point(&a.apex, &scale),
            ..a.clone()
        }),
        Surface::Sphere(a) => Surface::Sphere(Sphere {
            center: scale_definition_point(&a.center, &scale),
            radius: &a.radius * &scale,
            ..a.clone()
        }),
        Surface::Torus(a) => Surface::Torus(Torus {
            center: scale_definition_point(&a.center, &scale),
            major: &a.major * &scale,
            minor: &a.minor * &scale,
            ..a.clone()
        }),
        _ => unreachable!(),
    }
}
pub fn scale_curve(curve: &Curve, s: &Rational) -> Curve {
    let scale = Real::from(s.clone());
    let squared = Real::from(s.pow(2).unwrap());
    match curve {
        Curve::Line(a) => Curve::Line(Line {
            origin: scale_definition_point(&a.origin, &scale),
            ..a.clone()
        }),
        Curve::Circle(a) => Curve::Circle(Circle {
            center: scale_definition_point(&a.center, &scale),
            radius_sq: &a.radius_sq * &squared,
            ..a.clone()
        }),
        Curve::Ellipse(a) => Curve::Ellipse(Ellipse {
            center: scale_definition_point(&a.center, &scale),
            major_sq: &a.major_sq * &squared,
            minor_sq: &a.minor_sq * &squared,
            ..a.clone()
        }),
        _ => unreachable!(),
    }
}
