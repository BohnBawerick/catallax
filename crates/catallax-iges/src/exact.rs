//! Exact rational recovery and decimal encodings. Candidates are never evidence:
//! the shared decider must prove equality against the retained Real definition.
use crate::IgesError;
use catallax_numbers::{Decider, Integer, Outcome, Rational, Real};
use catallax_shapes::{Point3, Vector3};

pub(crate) fn rational(x: &Real, dec: &mut Decider) -> Result<Rational, IgesError> {
    // Also enclose rational leaves so caller work/deadline policies are observed.
    let ball = dec.enclose_at(x, 128)?;
    if let Some(q) = x.as_rational() {
        return Ok(q);
    }
    // Recover rational expressions such as sqrt(4/9) without replacing an
    // algebraic value by a nearby decimal. Continued fractions propose small
    // rationals; only an exact equality proof can accept one.
    let mid = ball.mid();
    let negative = mid < Rational::zero();
    let mut numerator = mid.numerator().abs();
    let mut denominator = mid.denominator().clone();
    let (mut p0, mut p1) = (Integer::zero(), Integer::one());
    let (mut q0, mut q1) = (Integer::one(), Integer::zero());
    for _ in 0..128 {
        if denominator.is_zero() {
            break;
        }
        let (a, remainder) = numerator.div_rem_floor(&denominator).unwrap();
        let p = &(&a * &p1) + &p0;
        let q = &(&a * &q1) + &q0;
        if q.bit_length() > 64 {
            break;
        }
        let candidate =
            Rational::new(if negative { -p.clone() } else { p.clone() }, q.clone()).unwrap();
        if ball.contains(&candidate) {
            match dec.compare(x, &Real::from(candidate.clone())) {
                Outcome::ProvenCoincident(_) => return Ok(candidate),
                Outcome::Refused(r) => return Err(r.into()),
                _ => {}
            }
        }
        (p0, p1) = (p1, p);
        (q0, q1) = (q1, q);
        numerator = denominator;
        denominator = remainder;
    }
    for digits in [17, 34] {
        let (lo, hi) = ball.mid().to_decimal_bounds(digits);
        for text in [lo, hi] {
            let q = Rational::parse(&text)
                .map_err(|_| IgesError::UnrepresentableDefinition("rational candidate"))?;
            match dec.compare(x, &Real::from(q.clone())) {
                Outcome::ProvenCoincident(_) => return Ok(q),
                Outcome::Refused(r) => return Err(r.into()),
                _ => {}
            }
        }
    }
    Err(IgesError::UnrepresentableDefinition(
        "non-rational numeric definition",
    ))
}
pub(crate) fn point(p: &Point3, dec: &mut Decider) -> Result<[Rational; 3], IgesError> {
    Ok([
        rational(&p.x, dec)?,
        rational(&p.y, dec)?,
        rational(&p.z, dec)?,
    ])
}
pub(crate) fn vector(p: &Vector3, dec: &mut Decider) -> Result<[Rational; 3], IgesError> {
    point(&Point3::new(p.x.clone(), p.y.clone(), p.z.clone()), dec)
}
pub(crate) fn unit(p: &Vector3, dec: &mut Decider) -> Result<[Rational; 3], IgesError> {
    vector(&p.scaled(&(Real::one() / p.norm_squared().sqrt())), dec)
}
pub(crate) fn lcm(a: &Integer, b: &Integer) -> Integer {
    &a.div_rem_floor(&a.gcd(b)).unwrap().0 * b
}
pub(crate) fn denominator(values: impl IntoIterator<Item = Rational>) -> Integer {
    values
        .into_iter()
        .fold(Integer::one(), |d, q| lcm(&d, q.denominator()))
}
pub(crate) fn integer_vector(v: [Rational; 3]) -> [Rational; 3] {
    let d = Rational::from_integer(denominator(v.clone()));
    v.map(|x| &x * &d)
}
pub(crate) fn dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    (0..3)
        .map(|i| &a[i] * &b[i])
        .fold(Rational::zero(), |a, b| a + b)
}
pub(crate) fn cross(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
pub(crate) fn sub(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
pub(crate) fn bound(p: &[Rational; 3]) -> Rational {
    p.iter().fold(Rational::zero(), |a, b| a + b.abs())
}
// Only terminating decimals reach this helper. The explicit limit fails closed.
pub(crate) fn decimal(q: &Rational) -> Result<String, IgesError> {
    for digits in [17, 34, 100, 512, 2048, 8192] {
        let (lo, hi) = q.to_decimal_bounds(digits);
        if lo == hi {
            return Ok(lo.replace('e', "E"));
        }
    }
    Err(IgesError::UnrepresentableDefinition(
        "decimal field exceeds 8192 significant digits",
    ))
}
pub(crate) fn scale(d: Integer, digits: u32) -> Result<Rational, IgesError> {
    let reciprocal = Rational::new(Integer::one(), d.clone()).unwrap();
    let (a, b) = reciprocal.to_decimal_bounds(digits);
    let a = Rational::parse(&a).map_err(|_| {
        IgesError::UnrepresentableDefinition("scale exponent exceeds numeric parser range")
    })?;
    let b = Rational::parse(&b).map_err(|_| {
        IgesError::UnrepresentableDefinition("scale exponent exceeds numeric parser range")
    })?;
    let q = if (&reciprocal - &a).abs() <= (&b - &reciprocal).abs() {
        a
    } else {
        b
    };
    Ok(Rational::from_integer(d) * q)
}
