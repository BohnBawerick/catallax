#![allow(dead_code, unused_imports, unused_macros)]

pub use catallax_intersect::*;
pub use catallax_numbers::*;
pub use catallax_shapes::*;
use std::collections::BTreeSet;

pub mod cases;
pub mod curves;
pub mod oracle;

// Passing blind contracts run in CI; comparison exclusions live at their tests.
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
pub fn n(x: i64) -> Rational {
    Rational::from_i64(x)
}
pub fn rr(x: i64) -> Real {
    Real::from_i64(x)
}
pub fn p(a: [&str; 3]) -> RPoint3 {
    let [x, y, z] = a.map(r);
    RPoint3::new(x, y, z)
}
pub fn v(a: [&str; 3]) -> RVector3 {
    let [x, y, z] = a.map(r);
    RVector3::new(x, y, z)
}
pub fn pi(a: [i64; 3]) -> RPoint3 {
    let [x, y, z] = a.map(n);
    RPoint3::new(x, y, z)
}
pub fn vi(a: [i64; 3]) -> RVector3 {
    let [x, y, z] = a.map(n);
    RVector3::new(x, y, z)
}
pub fn xyz(p: &Point3) -> [Real; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
pub fn rational_point(p: &Point3) -> Option<RPoint3> {
    Some(RPoint3::new(
        p.x.as_rational()?,
        p.y.as_rational()?,
        p.z.as_rational()?,
    ))
}
pub fn square_root(q: &Rational) -> Option<Rational> {
    let a: rug::Integer = q.numerator().to_string().parse().unwrap();
    let b: rug::Integer = q.denominator().to_string().parse().unwrap();
    if a < 0 {
        return None;
    }
    let x = a.clone().sqrt();
    let y = b.clone().sqrt();
    if x.clone() * &x == a && y.clone() * &y == b {
        Some(r(&format!("{x}/{y}")))
    } else {
        None
    }
}
pub fn root(q: &Rational) -> Real {
    square_root(q)
        .map(Real::from_rational)
        .unwrap_or_else(|| Real::from_rational(q.clone()).sqrt())
}
pub fn plus(o: &RPoint3, d: &RVector3, t: &Real) -> Point3 {
    o.to_real().plus(&d.to_real().scaled(t))
}
pub fn perp(w: &RVector3) -> RVector3 {
    let a = [w.x.abs(), w.y.abs(), w.z.abs()];
    let k = (0..3).min_by_key(|&i| a[i].clone()).unwrap();
    w.cross(&vi(std::array::from_fn(|i| i64::from(i == k))))
}
pub fn plane(o: RPoint3, normal: RVector3) -> Surface {
    let reference = if normal.x.is_zero() && normal.y.is_zero() {
        vi([1, 0, 0])
    } else if normal.x.is_zero() && normal.z.is_zero() {
        vi([0, 0, 1])
    } else if normal.y.is_zero() && normal.z.is_zero() {
        vi([0, 1, 0])
    } else {
        perp(&normal)
    };
    Surface::Plane(Plane::new(o, normal, reference).unwrap())
}
pub fn axis_plane(axis: usize, k: Rational) -> Surface {
    let mut a = [n(0), n(0), n(0)];
    a[axis] = k;
    let [x, y, z] = a;
    plane(
        RPoint3::new(x, y, z),
        vi(std::array::from_fn(|i| i64::from(i == axis))),
    )
}
pub fn cyl(o: RPoint3, axis: RVector3, radius: Rational) -> Surface {
    let reference = if axis.x.is_zero() {
        vi([1, 0, 0])
    } else if axis.y.is_zero() && axis.z.is_zero() {
        vi([0, 0, 1])
    } else if axis.y.is_zero() {
        vi([0, 1, 0])
    } else {
        vi([0, 0, 1])
    };
    Surface::Cylinder(Cylinder::new(o, axis, radius, reference).unwrap())
}
pub fn zc(radius: i64) -> Surface {
    cyl(pi([0, 0, 0]), vi([0, 0, 1]), n(radius))
}
pub fn sphere(o: RPoint3, radius: Rational) -> Surface {
    Surface::Sphere(Sphere::new(o, radius, vi([0, 0, 1]), vi([1, 0, 0])).unwrap())
}
pub fn k5() -> Curve {
    Curve::Circle(Circle::new(pi([0, 0, 0]), vi([0, 0, 1]), n(25), vi([1, 0, 0])).unwrap())
}
pub fn line28() -> Curve {
    Curve::Line(Line::new(pi([-10, 0, 0]), vi([1, 0, 0])).unwrap())
}
pub fn ellipse4() -> Curve {
    Curve::Ellipse(
        Ellipse::new(
            pi([0, 0, 0]),
            vi([0, 3, 4]),
            v(["0", "-12/25", "9/25"]),
            r("625/16"),
            n(25),
        )
        .unwrap(),
    )
}
pub fn closed(g: &mut GeomStore, p: Point3, supports: &[SurfaceId]) -> PointId {
    g.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(supports.to_vec()),
        position: Position::Closed(p),
    }))
}
pub fn explicit(g: &mut GeomStore, p: RPoint3) -> PointId {
    g.add_point(PointDef::Explicit(p))
}
pub fn exact<T: std::fmt::Debug>(a: Answered<T>) -> T {
    let a = a.unwrap();
    assert!(a.marks.is_empty(), "unexpected marks: {:?}", a.marks);
    a.value
}
pub fn coincident<T: std::fmt::Debug>(a: Outcome<T>) {
    assert!(a.is_coincident(), "expected coincidence: {a:?}");
}
pub fn real_eq(a: &Real, b: &Real) {
    coincident(Decider::default().compare(a, b));
}
pub fn point_eq(a: &Point3, b: &Point3) {
    for (a, b) in xyz(a).iter().zip(xyz(b)) {
        real_eq(a, &b);
    }
}
pub fn vector_eq(a: &Vector3, b: &Vector3) {
    real_eq(&a.x, &b.x);
    real_eq(&a.y, &b.y);
    real_eq(&a.z, &b.z);
}
pub fn positive(x: &Real) {
    assert_eq!(
        Decider::default().sign(x),
        Outcome::Proven(StrictSign::Positive)
    );
}
pub fn positive_multiple(a: &Vector3, b: &Vector3) {
    vector_eq(&a.cross(b), &vi([0, 0, 0]).to_real());
    positive(&a.dot(b));
}
pub fn expired() -> Decider {
    let mut d = Decider::default();
    d.set_deadline(Some(
        std::time::Instant::now() - std::time::Duration::from_secs(1),
    ));
    d
}
pub fn tiny() -> Decider {
    Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    })
}
pub fn starved() -> Decider {
    Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    })
}
pub fn fingerprint(g: &GeomStore) -> String {
    format!("{g:?}")
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
    pub fn pose(&mut self) -> Pose {
        let w = n((self.next() % 100 + 1) as i64);
        let [x, y, z] = std::array::from_fn(|_| n((self.next() % 201) as i64 - 100));
        let rotation = Rotation::from_quaternion(w, x, y, z).unwrap();
        let [x, y, z] = std::array::from_fn(|_| n((self.next() % 2_000_001) as i64 - 1_000_000));
        Pose::new(rotation, RVector3::new(x, y, z))
    }
}
pub fn scaled_surface(s: &Surface, k: &Rational) -> Surface {
    let point = |p: &RPoint3| RPoint3::new(&p.x * k, &p.y * k, &p.z * k);
    match s {
        Surface::Plane(s) => Surface::Plane(Plane {
            origin: point(&s.origin),
            ..s.clone()
        }),
        Surface::Cylinder(s) => Surface::Cylinder(Cylinder {
            origin: point(&s.origin),
            radius: &s.radius * k,
            ..s.clone()
        }),
        Surface::Sphere(s) => Surface::Sphere(Sphere {
            center: point(&s.center),
            radius: &s.radius * k,
            ..s.clone()
        }),
        _ => panic!("unsupported fixture"),
    }
}
pub fn scaled_curve(c: &Curve, k: &Rational) -> Curve {
    let point = |p: &RPoint3| RPoint3::new(&p.x * k, &p.y * k, &p.z * k);
    match c {
        Curve::Line(c) => Curve::Line(Line {
            origin: point(&c.origin),
            ..c.clone()
        }),
        Curve::Circle(c) => Curve::Circle(Circle {
            center: point(&c.center),
            radius_sq: &(&c.radius_sq * k) * k,
            ..c.clone()
        }),
        Curve::Ellipse(c) => Curve::Ellipse(Ellipse {
            center: point(&c.center),
            major_sq: &(&c.major_sq * k) * k,
            minor_sq: &(&c.minor_sq * k) * k,
            ..c.clone()
        }),
        _ => panic!("analytic fixture required"),
    }
}
