use catallax_numbers::{
    Coincidence, Decider, Integer, Outcome, Rational, Real, Refusal, Sign, StrictSign, ZeroProof,
};
use catallax_shapes::{Point3, RVector3, Vector3};

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
pub(super) fn r(x: i64) -> Rational {
    Rational::from_i64(x)
}
pub(super) fn real<T: Clone + Into<Real>>(x: &T) -> Real {
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
    if let Some(value) = x.rational_ref() {
        return Ok(value.sign());
    }
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
pub(super) fn root_real(x: &Real) -> Real {
    x.rational_ref().map_or_else(|| x.sqrt(), root)
}
pub(super) fn perp(v: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    if let Some([x, y, z]) = v.rational_components() {
        let a = [x.abs(), y.abs(), z.abs()];
        let k = (0..3).min_by_key(|&k| &a[k]).unwrap();
        let direction = match k {
            0 => RVector3::new(r(0), z.clone(), -y),
            1 => RVector3::new(-z, r(0), x.clone()),
            _ => RVector3::new(y.clone(), -x, r(0)),
        };
        return Ok(direction.into());
    }
    for axis in 0..3 {
        let components = std::array::from_fn::<_, 3, _>(|i| Real::from_i64(i64::from(i == axis)));
        let [x, y, z] = components;
        let direction = v.cross(&Vector3::new(x, y, z));
        if sign(&direction.norm_squared(), dec)? != Sign::Zero {
            return Ok(direction);
        }
    }
    Err(Refusal::new(
        "shapes.zero-vector",
        "perpendicular direction requires a nonzero vector",
    ))
}
pub(super) fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(Real::one() / root_real(&v.norm_squared())))
}
pub(super) fn midpoint(a: &Point3, b: &Point3) -> Point3 {
    a.plus(
        &b.minus(a)
            .scaled(&Real::from_rational(Rational::from_ratio(1, 2).unwrap())),
    )
}
