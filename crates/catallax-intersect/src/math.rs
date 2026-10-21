use catallax_numbers::{
    Coincidence, Decider, Integer, Outcome, Rational, Real, Refusal, Sign, StrictSign, ZeroProof,
};
use catallax_shapes::{Point3, RPoint3, RVector3, Vector3};

#[derive(Debug)]
pub(super) enum Failure {
    Uncertain(catallax_numbers::NotExact),
    Refused(Refusal),
}
impl From<Refusal> for Failure {
    fn from(e: Refusal) -> Self {
        Self::Refused(e)
    }
}
impl From<Failure> for Refusal {
    fn from(e: Failure) -> Self {
        match e {
            Failure::Uncertain(n) => Refusal::new(n.code, n.detail),
            Failure::Refused(r) => r,
        }
    }
}
impl Failure {
    pub fn outcome<T>(self) -> Outcome<T> {
        match self {
            Self::Uncertain(n) => Outcome::NotExact(n),
            Self::Refused(r) => Outcome::Refused(r),
        }
    }
}
pub(super) fn sq(x: &Rational) -> Rational {
    x * x
}
pub(super) fn r(x: i64) -> Rational {
    Rational::from_i64(x)
}
pub(super) fn real(x: &Rational) -> Real {
    x.clone().into()
}
pub(super) fn zero<T>() -> Outcome<T> {
    Outcome::ProvenCoincident(Coincidence {
        proof: ZeroProof::ExactArithmetic,
    })
}
pub(super) fn fact<T>(dec: &mut Decider) -> Outcome<T> {
    dec.note_fact_decision();
    Outcome::ProvenCoincident(Coincidence {
        proof: ZeroProof::ConstructionFact,
    })
}
pub(super) fn refused_curve() -> Refusal {
    Refusal::new(
        "intersect.unsupported-curve",
        "curve is outside piece 4 coverage",
    )
}
pub(super) fn refused_pair() -> Refusal {
    Refusal::new(
        "intersect.unsupported-pair",
        "surface pair is outside piece 4 coverage",
    )
}
pub(super) fn sign(x: &Real, dec: &mut Decider) -> Result<Sign, Failure> {
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Failure::Uncertain(n)),
        Outcome::Refused(e) => Err(Failure::Refused(e)),
    }
}
fn isqrt(n: &Integer) -> Integer {
    assert!(n.sign() != Sign::Negative);
    if n.is_zero() {
        return Integer::zero();
    }
    let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
    loop {
        let y = (&x + &n.div_rem_floor(&x).unwrap().0)
            .div_rem_floor(&Integer::from_i64(2))
            .unwrap()
            .0;
        if y >= x {
            return x;
        }
        x = y;
    }
}
pub(super) fn rational_sqrt(x: &Rational) -> Option<Rational> {
    if x.sign() == Sign::Negative {
        return None;
    }
    let n = isqrt(x.numerator());
    let d = isqrt(x.denominator());
    if &n * &n == *x.numerator() && &d * &d == *x.denominator() {
        Rational::new(n, d)
    } else {
        None
    }
}
pub(super) fn root(x: &Rational) -> Real {
    rational_sqrt(x).map_or_else(|| real(x).sqrt(), Real::from_rational)
}
pub(super) fn point_rational(p: &Point3) -> Option<RPoint3> {
    Some(RPoint3::new(
        p.x.as_rational()?,
        p.y.as_rational()?,
        p.z.as_rational()?,
    ))
}
pub(super) fn perp(v: &RVector3) -> RVector3 {
    let a = [v.x.abs(), v.y.abs(), v.z.abs()];
    let k = (0..3).min_by_key(|&k| &a[k]).unwrap();
    let mut e = [r(0), r(0), r(0)];
    e[k] = r(1);
    v.cross(&RVector3::new(e[0].clone(), e[1].clone(), e[2].clone()))
}
pub(super) fn unit(v: &RVector3) -> Vector3 {
    v.to_real().scaled(&(Real::one() / root(&v.norm_squared())))
}
pub(super) fn midpoint(a: &Point3, b: &Point3) -> Point3 {
    a.plus(
        &b.minus(a)
            .scaled(&Real::from_rational(Rational::from_ratio(1, 2).unwrap())),
    )
}
