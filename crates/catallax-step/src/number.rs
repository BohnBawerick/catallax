//! Decimal rounding and rational error accounting. No platform floating point.

use catallax_numbers::{
    Ball, BudgetPolicy, Decider, Integer, Outcome, Rational, Real, Refusal, StrictOrder,
};
use catallax_shapes::{Point3, Vector3};

#[derive(Clone)]
pub(crate) struct Number {
    pub(crate) value: Rational,
    pub(crate) text: String,
    pub(crate) error: Rational,
}

#[derive(Clone)]
pub(crate) struct Point {
    pub(crate) xyz: [Number; 3],
}

impl Point {
    pub(crate) fn error(&self) -> Rational {
        self.xyz.iter().fold(Rational::zero(), |s, n| &s + &n.error)
    }

    pub(crate) fn extent(&self, origin: &[Rational; 3]) -> Rational {
        self.xyz
            .iter()
            .zip(origin)
            .fold(Rational::zero(), |s, (n, o)| {
                &(s + (&n.value - o).abs()) + &n.error
            })
    }
}

pub(crate) fn decimal(s: &str) -> String {
    if let Some((mantissa, exponent)) = s.split_once('e') {
        let mantissa = if mantissa.contains('.') {
            mantissa.to_owned()
        } else {
            format!("{mantissa}.")
        };
        format!("{mantissa}E{exponent}")
    } else if s.contains('.') {
        s.to_owned()
    } else {
        format!("{s}.")
    }
}

pub(crate) fn rounded(lo: &Rational, hi: &Rational, digits: u32) -> Number {
    let mid = (lo + hi).mul_pow2(-1);
    let (lower, upper) = mid.to_decimal_bounds(digits.max(1));
    let a = formatted_value(&lower);
    let b = formatted_value(&upper);
    let (value, text) = if (&mid - &a).abs() <= (&b - &mid).abs() {
        (a, lower)
    } else {
        (b, upper)
    };
    let error = (&value - lo).abs().max((&value - hi).abs());
    Number {
        value,
        text: decimal(&text),
        error,
    }
}

pub(crate) fn rational(value: &Rational, digits: u32) -> Number {
    rounded(value, value, digits)
}

fn request(value: &Real, digits: u32) -> Rational {
    let scale = value.as_rational().unwrap_or_else(Rational::one).abs();
    let scale = if scale.is_zero() {
        Rational::one()
    } else {
        scale
    };
    scale
        * Rational::new(
            Integer::one(),
            Integer::from_i64(10).pow(digits.max(1).saturating_add(8)),
        )
        .unwrap()
}

pub(crate) fn real(value: &Real, digits: u32, dec: &mut Decider) -> Result<Number, Refusal> {
    // Even rational coordinates use enclose, so the deadline is observed here.
    let b = dec.enclose(value, &request(value, digits))?;
    enclosed(value, b, digits, dec)
}

pub(crate) fn enclosed(
    value: &Real,
    b: Ball,
    digits: u32,
    dec: &mut Decider,
) -> Result<Number, Refusal> {
    enclosed_by(value, b, digits, dec, |number, _, dec| {
        dec.compare(value, &Real::from_rational(number.clone()))
    })
}

fn enclosed_by(
    value: &Real,
    mut b: Ball,
    digits: u32,
    dec: &mut Decider,
    mut equal: impl FnMut(&Rational, &Ball, &mut Decider) -> Outcome<StrictOrder>,
) -> Result<Number, Refusal> {
    if let Some(q) = value.as_rational() {
        return Ok(rational(&q, digits));
    }
    loop {
        let mut number = rounded(&b.lo(), &b.hi(), digits);
        let scale = b.mid().abs().max(b.rad());
        let requested = scale
            * Rational::new(
                Integer::one(),
                Integer::from_i64(10).pow(digits.max(1).saturating_add(8)),
            )
            .unwrap();
        // Perfect-square expressions may still have a nonzero binary radius.
        // Certify equality only when the written decimal is inside the ball. A
        // tighter enclosure usually excludes it; only a surviving candidate pays
        // for a proof, whose zero test can be algebraically huge.
        if !b.rad().is_zero()
            && b.contains(&number.value)
            && dec
                .enclose(value, &requested)
                .map_or(true, |tight| tight.contains(&number.value))
        {
            match equal(&number.value, &b, dec) {
                Outcome::ProvenCoincident(_) => number.error = Rational::zero(),
                Outcome::Refused(refusal) => return Err(refusal),
                _ => {}
            }
        }
        if let Ok(number) = certified(number, digits) {
            return Ok(number);
        }
        // Repeat a relative request when cancellation hid the magnitude.
        // A stalled budget-limited enclosure cannot certify the text.
        if requested.is_zero() {
            return Err(budget("no decimal accuracy certificate"));
        }
        let next = dec.enclose(value, &requested)?;
        if next.rad() >= b.rad() {
            // A nonnegative enclosure of an exact algebraic zero can keep a
            // positive midpoint at every precision. A stalled enclosure is not
            // a reason to choose zero; an exact coincidence proof is.
            if b.contains(&Rational::zero()) {
                match dec.sign(value) {
                    Outcome::ProvenCoincident(_) => {
                        return Ok(rational(&Rational::zero(), digits));
                    }
                    Outcome::Refused(refusal) => return Err(refusal),
                    _ => {}
                }
            }
            return Err(budget("enclosure too wide to certify the decimal"));
        }
        b = next;
    }
}

pub(crate) fn point(value: &Point3, digits: u32, dec: &mut Decider) -> Result<Point, Refusal> {
    Ok(Point {
        xyz: [
            real(&value.x, digits, dec)?,
            real(&value.y, digits, dec)?,
            real(&value.z, digits, dec)?,
        ],
    })
}

pub(crate) fn vertex(
    model: &catallax_topology::Model,
    id: catallax_topology::VertexId,
    options: &crate::StepOptions,
    dec: &mut Decider,
) -> Result<Point, Refusal> {
    let point_id = model.vertex(id).point;
    // Observe the public enclosure/deadline contract even for rational points.
    let [x, y, z] = model.geom().enclose_point(
        point_id,
        &options
            .curve_tolerance
            .abs()
            .max(Rational::from_ratio(1, 1_000_000_000).unwrap()),
        dec,
    )?;
    let position = model.geom().position(point_id);
    Ok(Point {
        xyz: [
            enclosed(&position.x, x, options.significant_digits.max(1), dec)?,
            enclosed(&position.y, y, options.significant_digits.max(1), dec)?,
            enclosed(&position.z, z, options.significant_digits.max(1), dec)?,
        ],
    })
}

pub(crate) fn sqrt(value: &Rational) -> Real {
    let n = integer_sqrt(value.numerator());
    let d = integer_sqrt(value.denominator());
    if &n * &n == *value.numerator() && &d * &d == *value.denominator() {
        Real::from_rational(Rational::new(n, d).unwrap())
    } else {
        Real::from_rational(value.clone()).sqrt()
    }
}
pub(crate) fn sqrt_real(value: &Real) -> Real {
    value
        .as_rational()
        .map_or_else(|| value.sqrt(), |q| sqrt(&q))
}

fn integer_sqrt(value: &Integer) -> Integer {
    if value.is_zero() {
        return Integer::zero();
    }
    let mut x = Integer::one().mul_pow2(value.bit_length().div_ceil(2));
    loop {
        let q = value.div_rem_floor(&x).unwrap().0;
        let y = (&x + &q).div_rem_floor(&Integer::from_i64(2)).unwrap().0;
        if y >= x {
            return x;
        }
        x = y;
    }
}

pub(crate) fn unit(value: &Vector3) -> Vector3 {
    value.scaled(&(Real::one() / sqrt_real(&value.norm_squared())))
}

pub(crate) fn direction(value: &Vector3, digits: u32, dec: &mut Decider) -> Result<Point, Refusal> {
    let v = unit(value);
    let norm_squared = value.norm_squared();
    Ok(Point {
        xyz: [
            component(&v.x, &value.x, &norm_squared, digits, dec)?,
            component(&v.y, &value.y, &norm_squared, digits, dec)?,
            component(&v.z, &value.z, &norm_squared, digits, dec)?,
        ],
    })
}

/// One component of `source / sqrt(norm_squared)`. A ball holding the component
/// and a nonzero decimal gives them one sign, so `source^2 = decimal^2 *
/// norm_squared` proves equality without the radical.
fn component(
    value: &Real,
    source: &Real,
    norm_squared: &Real,
    digits: u32,
    dec: &mut Decider,
) -> Result<Number, Refusal> {
    let b = dec.enclose(value, &request(value, digits))?;
    enclosed_by(value, b, digits, dec, |decimal, b, dec| {
        let decimal = Real::from_rational(decimal.clone());
        if !b.contains_zero() {
            match dec.compare(&source.square(), &(&decimal.square() * norm_squared)) {
                Outcome::NotExact(_) => {}
                settled => return settled,
            }
        }
        dec.compare(value, &decimal)
    })
}

pub(crate) fn ball_number(value: &Ball, digits: u32) -> Number {
    rounded(&value.lo(), &value.hi(), digits)
}

pub(crate) fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

pub(crate) fn point_plus(origin: &Point3, vector: &Vector3, amount: &Real) -> Point3 {
    origin.plus(&vector.scaled(amount))
}

pub(crate) fn upper(value: &Real, dec: &Decider) -> Result<Rational, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.clone());
    }
    Ok(dec.enclose_at(value, 64)?.hi())
}

pub(crate) fn budget(detail: impl Into<String>) -> Refusal {
    Refusal::new("step.budget-exhausted", detail)
}

pub(crate) fn certified(number: Number, digits: u32) -> Result<Number, Refusal> {
    if number.error.is_zero() {
        return Ok(number);
    }
    if !number.value.is_zero() {
        let lower = number.value.abs().to_decimal_bounds(1).0;
        let leading = lower.bytes().find(|c| matches!(c, b'1'..=b'9')).unwrap();
        let decade = formatted_value(&lower) / Rational::from_i64(i64::from(leading - b'0'));
        let unit = decade / Rational::from_integer(Integer::from_i64(10).pow(digits.max(1) - 1));
        if number.error <= unit {
            return Ok(number);
        }
    }
    Err(budget("enclosure too wide to certify the decimal"))
}

pub(crate) fn work_precision(policy: &BudgetPolicy, desired: u32, operations: u64) -> Option<u32> {
    assert!(operations > 0);
    let limbs = (policy.max_work / operations).isqrt();
    let precision = limbs
        .saturating_mul(64)
        .min(u64::from(policy.max_precision_bits))
        .min(u64::from(desired.max(1))) as u32;
    (precision != 0).then_some(precision)
}

pub(crate) fn floor_ball(value: &Ball, digits: u32) -> Number {
    let lower = value.lo();
    let upper = value.hi();
    let text = lower.to_decimal_bounds(digits.max(1)).0;
    let number = formatted_value(&text);
    let error = (&number - &lower).abs().max((&number - &upper).abs());
    Number {
        value: number,
        text: decimal(&text),
        error,
    }
}

pub(crate) fn formatted_value(text: &str) -> Rational {
    // The input parser limits exponents, but arithmetic can cross that limit.
    // Only text produced by to_decimal_bounds reaches this conversion.
    let Some((mantissa, exponent)) = text.split_once(['e', 'E']) else {
        return Rational::parse(text).expect("decimal formatter mantissa");
    };
    let value = Rational::parse(mantissa).expect("decimal formatter mantissa");
    let exponent: i64 = exponent.parse().expect("decimal formatter exponent");
    let scale = Rational::from_integer(
        Integer::from_i64(10).pow(
            exponent
                .unsigned_abs()
                .try_into()
                .expect("stored decimal exponent"),
        ),
    );
    if exponent < 0 {
        value / scale
    } else {
        value * scale
    }
}
