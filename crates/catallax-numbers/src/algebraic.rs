//! Real algebraic numbers and filtered Descartes root isolation.

use crate::decider::{BudgetPolicy, Work};
use crate::{Ball, Integer, Rational, Sign, UPoly};
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};

const NEWTON_GUARD_BITS: i64 = 64;

/// Exact rational, or a single real root in an open rational interval.
#[derive(Clone, Debug)]
pub struct Algebraic {
    pub(crate) repr: AlgebraicRepr,
}
#[derive(Clone, Debug)]
pub(crate) enum AlgebraicRepr {
    Rational(Rational),
    Root(Arc<AlgebraicRoot>),
}

#[derive(Debug)]
pub(crate) struct AlgebraicRoot {
    pub polynomial: UPoly,
    pub lo: Rational,
    pub hi: Rational,
    pub enclosures: Mutex<BTreeMap<u32, RootEnclosure>>,
}
#[derive(Debug)]
pub(crate) struct RootEnclosure {
    pub ball: Ball,
    pub cost: u64,
    pub lo: Rational,
    pub hi: Rational,
    pub lower_sign: Option<Sign>,
}
impl AlgebraicRoot {
    fn new(polynomial: UPoly, lo: Rational, hi: Rational) -> Self {
        Self {
            polynomial,
            lo,
            hi,
            enclosures: Mutex::new(BTreeMap::new()),
        }
    }
}
impl Clone for AlgebraicRoot {
    fn clone(&self) -> Self {
        Self::new(self.polynomial.clone(), self.lo.clone(), self.hi.clone())
    }
}

impl Algebraic {
    pub(crate) fn isolated(polynomial: UPoly, lo: Rational, hi: Rational) -> Self {
        Self {
            repr: AlgebraicRepr::Root(Arc::new(AlgebraicRoot::new(polynomial, lo, hi))),
        }
    }

    /// Constructs an exact rational algebraic number.
    pub fn from_rational(r: Rational) -> Self {
        Self {
            repr: AlgebraicRepr::Rational(r),
        }
    }
    /// Exact rational value, when applicable.
    pub fn as_rational(&self) -> Option<&Rational> {
        if let AlgebraicRepr::Rational(r) = &self.repr {
            Some(r)
        } else {
            None
        }
    }
    /// Square-free primitive integer defining polynomial with positive leading coefficient.
    pub fn defining_polynomial(&self) -> UPoly {
        match &self.repr {
            AlgebraicRepr::Rational(r) => UPoly::new(vec![
                Rational::from_integer(-r.numerator()),
                Rational::from_integer(r.denominator().clone()),
            ]),
            AlgebraicRepr::Root(root) => root.polynomial.clone(),
        }
    }
    /// Open isolating interval, or the exact rational twice.
    pub fn isolating_interval(&self) -> (Rational, Rational) {
        let (lo, hi) = self.interval_refs();
        (lo.clone(), hi.clone())
    }
    fn interval_refs(&self) -> (&Rational, &Rational) {
        match &self.repr {
            AlgebraicRepr::Rational(r) => (r, r),
            AlgebraicRepr::Root(root) => (&root.lo, &root.hi),
        }
    }
    /// Refines this same number to a positive requested width.
    pub fn refined(&self, max_width: &Rational) -> Self {
        let mut work = Work::new(
            BudgetPolicy {
                max_precision_bits: u32::MAX,
                max_work: u64::MAX,
            },
            None,
        );
        self.refined_checked(max_width, &mut work)
            .expect("unbounded algebraic refinement")
    }
    fn refined_checked(&self, max_width: &Rational, work: &mut Work) -> Result<Self, &'static str> {
        work.check_deadline()?;
        assert!(max_width.sign() == Sign::Positive);
        let AlgebraicRepr::Root(root) = &self.repr else {
            return Ok(self.clone());
        };
        let AlgebraicRoot {
            polynomial, lo, hi, ..
        } = root.as_ref();
        if hi - lo <= *max_width {
            return Ok(self.clone());
        }
        let mut lo = lo.clone();
        let mut hi = hi.clone();
        let derivative = polynomial.derivative();
        let s = filtered_sign(polynomial, &lo, work)?;
        while &hi - &lo > *max_width {
            work.exact(&lo, &hi, false)?;
            let mid = (&lo + &hi).mul_pow2(-1);
            let t = filtered_sign(polynomial, &mid, work)?;
            if t == Sign::Zero {
                return Ok(Self::from_rational(mid));
            }
            if let Some(step) =
                newton_bracket(polynomial, &derivative, &lo, &hi, &mid, max_width, work)?
            {
                match step {
                    Narrowing::Exact(root) => return Ok(Self::from_rational(root)),
                    Narrowing::Bracket(a, b) => {
                        lo = a;
                        hi = b;
                        continue;
                    }
                }
            }
            if s == t {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(Self::isolated(polynomial.clone(), lo, hi))
    }
    /// Proven ball enclosure.
    pub fn to_ball(&self, precision: u32) -> Ball {
        if let Some(r) = self.as_rational() {
            return Ball::from_rational(r, precision);
        }
        let refined = self.refined(&Rational::one().mul_pow2(-i64::from(precision.max(1))));
        let (lo, hi) = refined.isolating_interval();
        Ball::bounds(lo, hi, precision)
    }
}
/// One real root with its multiplicity in the input polynomial.
#[derive(Clone, Debug)]
pub struct RealRoot {
    /// The exact root.
    pub value: Algebraic,
    /// Multiplicity.
    pub multiplicity: u32,
}

pub(crate) fn isolate(p: &UPoly, work: &mut Work) -> Result<Vec<RealRoot>, &'static str> {
    let factors = if p.degree() == Some(2) {
        // A nonzero quadratic discriminant is a complete square-free
        // certificate. Avoid a generic polynomial GCD reservation when the
        // exact quadratic path can isolate these roots directly.
        let [c, b, a] = p.coefficients() else {
            unreachable!("quadratic degree has three coefficients")
        };
        work.exact(b, b, true)?;
        let squared = b * b;
        work.exact(a, c, true)?;
        let product = (a * c).mul_pow2(2);
        work.exact(&squared, &product, false)?;
        let discriminant = squared - product;
        if discriminant.is_zero() {
            let numerator = -b;
            let denominator = a.mul_pow2(1);
            work.exact(&numerator, &denominator, true)?;
            let root = numerator / denominator;
            return Ok(vec![RealRoot {
                value: Algebraic::from_rational(root),
                multiplicity: 2,
            }]);
        }
        vec![(p.clone(), 1)]
    } else {
        work.poly(p, p)?;
        p.square_free_factors_checked(&|| work.check_deadline())?
    };
    let mut roots = Vec::new();
    for (f, multiplicity) in factors {
        for value in isolate_square_free(&f, work)? {
            roots.push(RealRoot {
                value,
                multiplicity,
            });
        }
    }
    // Disjoint intervals are required across factors, not just within each factor.
    let mut width = Rational::one();
    loop {
        work.check_deadline()?;
        roots.sort_by(|a, b| a.value.interval_refs().0.cmp(b.value.interval_refs().0));
        if roots
            .windows(2)
            .all(|w| w[0].value.interval_refs().1 < w[1].value.interval_refs().0)
        {
            break;
        }
        width = width.mul_pow2(-1);
        let widest = roots
            .iter()
            .map(|root| {
                let (lo, hi) = root.value.interval_refs();
                hi - lo
            })
            .max()
            .expect("overlapping roots are nonempty");
        if widest <= width {
            // Skip only targets that would leave every interval unchanged.
            // Keep the original power-of-two refinement grid and its brackets.
            width = Rational::one().mul_pow2(widest.binary_exponent());
            if width == widest {
                width = width.mul_pow2(-1);
            }
        }
        for root in &mut roots {
            if let AlgebraicRepr::Root(data) = &root.value.repr {
                if &data.hi - &data.lo <= width {
                    continue;
                }
                work.poly(&data.polynomial, &data.polynomial)?;
                root.value = root.value.refined_checked(&width, work)?;
            }
        }
    }
    Ok(roots)
}

// For (x²+u x+v)(x²+w x+t), y=v+t is a root of
// y³-b y²+(a c-4d)y+4bd-a²d-c². A rational factorization
// has rational y and u-w; if u=w, recover v-t instead.
fn split_quartic(f: &UPoly, work: &mut Work) -> Result<Option<[UPoly; 2]>, &'static str> {
    let c = f.coefficients();
    let (a, b, c, d) = (&c[3], &c[2], &c[1], &c[0]);
    debug_assert_eq!(f.leading(), Some(&Rational::one()));
    let half_a = a.mul_pow2(-1);
    work.exact(&half_a, &half_a, true)?;
    let square = &half_a * &half_a;
    work.exact(b, &square, false)?;
    let y = b - &square;
    work.exact(&half_a, &y, true)?;
    if &half_a * &y == *c
        && let Some(factors) = quadratics_for_sum(a, b, c, d, &y, work)?
    {
        return Ok(Some(factors));
    }
    work.exact(a, c, true)?;
    let ac = a * c;
    work.exact(&ac, &d.mul_pow2(2), false)?;
    let linear = &ac - &d.mul_pow2(2);
    work.exact(b, d, true)?;
    let four_bd = (b * d).mul_pow2(2);
    work.exact(a, a, true)?;
    let a2 = a * a;
    work.exact(&a2, d, true)?;
    let a2d = &a2 * d;
    work.exact(c, c, true)?;
    let c2 = c * c;
    work.exact(&four_bd, &a2d, false)?;
    let constant = &four_bd - &a2d;
    work.exact(&constant, &c2, false)?;
    let resolvent = UPoly::new(vec![constant - c2, linear, -b, Rational::one()]);
    if !rational_root_possible(&resolvent, work)? {
        return Ok(None);
    }
    let coefficients = resolvent.primitive_integer_checked(&|| work.check_deadline())?;
    if let Some(candidates) = rational_candidates(&coefficients, work)? {
        for y in candidates {
            if filtered_sign(&resolvent, &y, work)? == Sign::Zero
                && let Some(factors) = quadratics_for_sum(a, b, c, d, &y, work)?
            {
                return Ok(Some(factors));
            }
        }
        return Ok(None);
    }
    // The optional divisor search declined before completeness. Preserve the
    // exact continued-fraction fallback used by the full resolver isolator.
    for root in isolate_square_free(&resolvent, work)? {
        if let Some(y) = root.as_rational()
            && let Some(factors) = quadratics_for_sum(a, b, c, d, y, work)?
        {
            return Ok(Some(factors));
        }
    }
    Ok(None)
}

// A rational root of a monic cubic reduces modulo any prime avoiding its
// coefficient denominators. A rootless reduction excludes every rational root.
pub(crate) fn rational_root_possible(
    polynomial: &UPoly,
    work: &mut Work,
) -> Result<bool, &'static str> {
    debug_assert_eq!(polynomial.degree(), Some(3));
    debug_assert_eq!(polynomial.coefficients()[3], Rational::one());
    'prime: for prime in [
        3u32, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    ] {
        work.charge(u64::from(prime) * 16)?;
        let modulo = |value: &Integer| {
            Integer(&value.0 % i64::from(prime))
                .to_i64()
                .unwrap()
                .rem_euclid(i64::from(prime)) as u32
        };
        let mut coefficients = [0u32; 4];
        for (out, coefficient) in coefficients.iter_mut().zip(polynomial.coefficients()) {
            work.exact(coefficient, coefficient, false)?;
            let denominator = modulo(coefficient.denominator());
            if denominator == 0 {
                continue 'prime;
            }
            let inverse = (1..prime).find(|x| (x * denominator) % prime == 1).unwrap();
            *out = (modulo(coefficient.numerator()) * inverse) % prime;
        }
        if !(0..prime).any(|x| {
            coefficients
                .iter()
                .rev()
                .fold(0, |value, coefficient| (value * x + coefficient) % prime)
                == 0
        }) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn quadratics_for_sum(
    a: &Rational,
    b: &Rational,
    c: &Rational,
    d: &Rational,
    y: &Rational,
    work: &mut Work,
) -> Result<Option<[UPoly; 2]>, &'static str> {
    work.exact(a, a, true)?;
    let a2 = a * a;
    work.exact(b, y, false)?;
    let four_uw = (b - y).mul_pow2(2);
    work.exact(&a2, &four_uw, false)?;
    let Some(s) = crate::quadratic::rational_sqrt(&(&a2 - &four_uw), work)? else {
        return Ok(None);
    };
    work.exact(a, &s, false)?;
    let u = (a + &s).mul_pow2(-1);
    work.exact(a, &s, false)?;
    let w = (a - &s).mul_pow2(-1);
    let (v, t) = if s.is_zero() {
        work.exact(y, y, true)?;
        let y2 = y * y;
        work.exact(&y2, &d.mul_pow2(2), false)?;
        let Some(delta) = crate::quadratic::rational_sqrt(&(&y2 - &d.mul_pow2(2)), work)? else {
            return Ok(None);
        };
        work.exact(y, &delta, false)?;
        let v = (y + &delta).mul_pow2(-1);
        work.exact(y, &delta, false)?;
        (v, (y - &delta).mul_pow2(-1))
    } else {
        work.exact(&u, y, true)?;
        let uy = &u * y;
        work.exact(&uy, c, false)?;
        let numerator = &uy - c;
        work.exact(&numerator, &s, true)?;
        let v = &numerator / &s;
        work.exact(y, &v, false)?;
        let t = y - &v;
        (v, t)
    };
    Ok(Some([
        UPoly::new(vec![v, u, Rational::one()]),
        UPoly::new(vec![t, w, Rational::one()]),
    ]))
}
enum QuadraticIsolation {
    Complete,
    NoRationalRoots,
}

fn isolate_quadratic(
    polynomial: &UPoly,
    target: &Rational,
    output: &mut Vec<Algebraic>,
    work: &mut Work,
) -> Result<QuadraticIsolation, &'static str> {
    let [c, b, a] = polynomial.coefficients() else {
        unreachable!("quadratic isolator requires degree two");
    };
    work.exact(b, b, true)?;
    let squared = b * b;
    work.exact(a, c, true)?;
    let product = (a * c).mul_pow2(2);
    work.exact(&squared, &product, false)?;
    let discriminant = squared - product;
    if discriminant.sign() == Sign::Negative {
        return Ok(QuadraticIsolation::Complete);
    }
    let center = -b;
    let denominator = a.mul_pow2(1);
    if let Some(root) = crate::quadratic::rational_sqrt(&discriminant, work)? {
        work.exact(&center, &root, false)?;
        let lower = &center - &root;
        work.exact(&lower, &denominator, true)?;
        let lower = Algebraic::from_rational(&lower / &denominator);
        if root.is_zero() {
            output.push(lower);
            return Ok(QuadraticIsolation::Complete);
        }
        work.exact(&center, &root, false)?;
        let upper = &center + &root;
        work.exact(&upper, &denominator, true)?;
        output.extend([lower, Algebraic::from_rational(&upper / &denominator)]);
        return Ok(QuadraticIsolation::Complete);
    }
    if work.policy.max_precision_bits == 0 {
        return Ok(QuadraticIsolation::NoRationalRoots);
    }
    // Try a directly certified closed-form bracket at the ordinary
    // rational-reconstruction width. Exact endpoint signs and the width check
    // below gate every enclosure; an inconclusive ball falls back to Descartes.
    let magnitude = (polynomial
        .coefficients()
        .iter()
        .filter(|x| !x.is_zero())
        .map(Rational::binary_exponent)
        .max()
        .unwrap()
        - a.binary_exponent()
        + 2)
    .max(1);
    let precision = (2 * a.numerator().bit_length() as i64 + magnitude + 16)
        .clamp(1, i64::from(work.policy.max_precision_bits)) as u32;
    work.ball(precision, true)?;
    let discriminant = Ball::from_rational(&discriminant, precision);
    work.ball(precision, true)?;
    let Some(root) = discriminant.sqrt() else {
        return Ok(QuadraticIsolation::NoRationalRoots);
    };
    work.ball(precision, true)?;
    let center = Ball::from_rational(&center, precision);
    work.ball(precision, true)?;
    let denominator = Ball::from_rational(&denominator, precision);
    let mut roots = [None, None];
    for (index, subtract) in [true, false].into_iter().enumerate() {
        work.ball(precision, false)?;
        let numerator = if subtract {
            &center - &root
        } else {
            &center + &root
        };
        work.ball(precision, true)?;
        work.ball(precision, true)?;
        let Some(enclosure) = numerator.checked_div(&denominator) else {
            return Ok(QuadraticIsolation::NoRationalRoots);
        };
        let lo = enclosure.lo();
        let hi = enclosure.hi();
        work.exact(&lo, &hi, false)?;
        if lo >= hi || &hi - &lo >= *target {
            return Ok(QuadraticIsolation::NoRationalRoots);
        }
        let lower_sign = filtered_sign(polynomial, &lo, work)?;
        let upper_sign = filtered_sign(polynomial, &hi, work)?;
        if lower_sign == Sign::Zero || upper_sign == Sign::Zero || lower_sign == upper_sign {
            return Ok(QuadraticIsolation::NoRationalRoots);
        }
        roots[index] = Some(Algebraic::isolated(polynomial.clone(), lo, hi));
    }
    output.extend(roots.map(Option::unwrap));
    Ok(QuadraticIsolation::Complete)
}

fn isolate_square_free(f: &UPoly, work: &mut Work) -> Result<Vec<Algebraic>, &'static str> {
    if f.degree() == Some(1) {
        let a = &f.coefficients()[0];
        let b = &f.coefficients()[1];
        work.exact(a, b, true)?;
        return Ok(vec![Algebraic::from_rational((-a) / b.clone())]);
    }
    let mut remaining = f.clone();
    let mut result = Vec::new();
    if remaining.coefficients()[0].is_zero() {
        result.push(Algebraic::from_rational(Rational::zero()));
        remaining = UPoly::new(remaining.coefficients()[1..].to_vec());
    }
    let integers = remaining.primitive_integer_checked(&|| work.check_deadline())?;
    let candidate_primitive = UPoly::new(
        integers
            .iter()
            .cloned()
            .map(Rational::from_integer)
            .collect(),
    );
    let candidate_quadratic = candidate_primitive.degree() == Some(2);
    let mut quadratic_proves_no_rational_roots = false;
    if candidate_quadratic {
        let q = integers.last().unwrap().abs();
        let target = Rational::new(Integer::one(), (&q * &q).mul_pow2(2)).unwrap();
        match isolate_quadratic(&candidate_primitive, &target, &mut result, work)? {
            QuadraticIsolation::Complete => return Ok(result),
            QuadraticIsolation::NoRationalRoots => {
                // rational_sqrt proved the discriminant nonsquare, so the
                // optional rational-root theorem search cannot find anything.
                quadratic_proves_no_rational_roots = true;
            }
        }
    }
    let candidates = if quadratic_proves_no_rational_roots {
        None
    } else {
        rational_candidates(&integers, work)?
    };
    let mut rational_complete = quadratic_proves_no_rational_roots || candidates.is_some();
    if let Some(candidates) = candidates {
        for r in candidates {
            if filtered_sign(&remaining, &r, work)? == Sign::Zero {
                work.poly(&remaining, &remaining)?;
                remaining = remaining
                    .div_rem_checked(&UPoly::new(vec![-&r, Rational::one()]), &|| {
                        work.check_deadline()
                    })?
                    .unwrap()
                    .0;
                result.push(Algebraic::from_rational(r));
            }
        }
    }
    if remaining.degree().is_none_or(|d| d == 0) {
        return Ok(result);
    }
    if remaining.degree() == Some(4)
        && let Some(factors) = split_quartic(&remaining, work)?
    {
        for factor in factors {
            result.extend(isolate_square_free(&factor, work)?);
        }
        return Ok(result);
    }
    let integers = remaining.primitive_integer_checked(&|| work.check_deadline())?;
    let primitive = UPoly::new(
        integers
            .iter()
            .cloned()
            .map(Rational::from_integer)
            .collect(),
    );
    let q = integers.last().unwrap().abs();
    let target = Rational::new(Integer::one(), (&q * &q).mul_pow2(2)).unwrap();
    if !candidate_quadratic && remaining.degree() == Some(2) && !rational_complete {
        match isolate_quadratic(&primitive, &target, &mut result, work)? {
            QuadraticIsolation::Complete => return Ok(result),
            QuadraticIsolation::NoRationalRoots => rational_complete = true,
        }
    }
    let irrational_start = result.len();
    let mut derivative = None;
    let bound = root_bound(&primitive, work)?;
    // Zero was removed above. All other roots lie strictly inside (-B, B).
    let mut pending = vec![(-&bound, Rational::zero()), (Rational::zero(), bound)];
    while let Some((lo, hi)) = pending.pop() {
        let variations = descartes_variations(&primitive, &lo, &hi, work)?;
        if variations == 0 {
            continue;
        }
        if variations == 1 {
            let mut lo = lo;
            let mut hi = hi;
            let s = filtered_sign(&primitive, &lo, work)?;
            let value = loop {
                if rational_complete {
                    break Algebraic::isolated(primitive.clone(), lo, hi);
                }
                work.exact(&lo, &hi, false)?;
                let mid = (&lo + &hi).mul_pow2(-1);
                let t = filtered_sign(&primitive, &mid, work)?;
                if t == Sign::Zero {
                    break Algebraic::from_rational(mid);
                }
                if &hi - &lo < target {
                    if let Some(r) = rational_root(&primitive, &mid, &q, &lo, &hi, work)? {
                        break Algebraic::from_rational(r);
                    }
                    break Algebraic::isolated(primitive.clone(), lo, hi);
                }
                if let Some(step) = newton_bracket(
                    &primitive,
                    derivative.get_or_insert_with(|| primitive.derivative()),
                    &lo,
                    &hi,
                    &mid,
                    &target,
                    work,
                )? {
                    match step {
                        Narrowing::Exact(root) => break Algebraic::from_rational(root),
                        Narrowing::Bracket(a, b) => {
                            lo = a;
                            hi = b;
                            continue;
                        }
                    }
                }
                if s == t {
                    lo = mid;
                } else {
                    hi = mid;
                }
            };
            result.push(value);
            continue;
        }
        work.exact(&lo, &hi, false)?;
        let mut mid = (&lo + &hi).mul_pow2(-1);
        // Never put a root at a child endpoint. Successive quarter splits stay
        // dyadic, and a nonzero polynomial has only finitely many such roots.
        while filtered_sign(&primitive, &mid, work)? == Sign::Zero {
            work.exact(&lo, &mid, false)?;
            mid = (&lo + &mid).mul_pow2(-1);
        }
        pending.push((mid.clone(), hi));
        pending.push((lo, mid));
    }
    // Continued-fraction recovery can discover rational roots that bounded
    // divisor enumeration missed. Remove them from the other roots' defining
    // polynomial too, so degree <= 4 still certifies an irreducible factor.
    if !rational_complete && primitive.degree().is_some_and(|d| d <= 4) {
        let mut factor = primitive.clone();
        for root in &result[irrational_start..] {
            if let Some(r) = root.as_rational() {
                work.poly(&factor, &factor)?;
                factor = factor
                    .div_rem_checked(&UPoly::new(vec![-r, Rational::one()]), &|| {
                        work.check_deadline()
                    })?
                    .unwrap()
                    .0;
            }
        }
        if factor.degree() != primitive.degree() {
            let integers = factor.primitive_integer_checked(&|| work.check_deadline())?;
            let factor = UPoly::new(integers.into_iter().map(Rational::from_integer).collect());
            for root in &mut result[irrational_start..] {
                if let AlgebraicRepr::Root(root) = &mut root.repr {
                    let root = Arc::make_mut(root);
                    root.polynomial = factor.clone();
                    root.enclosures.get_mut().unwrap().clear();
                }
            }
        }
    }
    result.sort_by(|a, b| a.isolating_interval().0.cmp(&b.isolating_interval().0));
    Ok(result)
}

enum Narrowing {
    Exact(Rational),
    Bracket(Rational, Rational),
}

// Newton only proposes an interval. Opposite certified endpoint signs, inside
// an existing one-root interval, are the proof; a failed proposal changes nothing.
fn newton_bracket(
    p: &UPoly,
    derivative: &UPoly,
    lo: &Rational,
    hi: &Rational,
    mid: &Rational,
    target: &Rational,
    work: &mut Work,
) -> Result<Option<Narrowing>, &'static str> {
    work.check_deadline()?;
    if work.policy.max_precision_bits == 0 {
        return Ok(None);
    }
    let width = hi - lo;
    let width_bits = -width.binary_exponent();
    let exponent = if mid.is_zero() {
        0
    } else {
        mid.binary_exponent()
    };
    // Near a simple root of magnitude 2^e the Newton error is O(width^2/2^e).
    // Leave guard bits for conditioning; exact signs reject overbold proposals.
    // Absolute dyadic bits may be negative for large-scale roots.
    let bits = (2 * width_bits + exponent - 12)
        .max(width_bits + 4)
        .min(-target.binary_exponent() + 3);
    let precision = (bits + exponent + NEWTON_GUARD_BITS)
        .clamp(1, i64::from(work.policy.max_precision_bits)) as u32;
    work.ball(precision, true)?;
    let x = Ball::from_rational(mid, precision);
    let value = eval_ball_checked(p, &x, work)?;
    let slope = eval_ball_checked(derivative, &x, work)?;
    // Ball division is a reciprocal and a multiplication.
    work.ball(precision, true)?;
    work.ball(precision, true)?;
    let Some(correction) = value.checked_div(&slope) else {
        return Ok(None);
    };
    let correction = correction.mid();
    work.exact(mid, &correction, false)?;
    let proposal = mid - &correction;
    let scaled = proposal.mul_pow2(bits);
    work.exact(&scaled, &scaled, true)?;
    let center = Rational::from_integer(scaled.floor()).mul_pow2(-bits);
    let radius = Rational::one().mul_pow2(1 - bits);
    work.exact(&center, &radius, false)?;
    let a = (&center - &radius).max(lo.clone());
    work.exact(&center, &radius, false)?;
    let b = (&center + &radius).min(hi.clone());
    if a >= b || (&b - &a).mul_pow2(1) > width {
        return interval_newton_bracket(p, derivative, (lo, hi), &x, &value, &width, work);
    }
    let sa = filtered_sign(p, &a, work)?;
    if sa == Sign::Zero {
        return Ok(Some(Narrowing::Exact(a)));
    }
    let sb = filtered_sign(p, &b, work)?;
    if sb == Sign::Zero {
        return Ok(Some(Narrowing::Exact(b)));
    }
    if sa != sb {
        return Ok(Some(Narrowing::Bracket(a, b)));
    }
    interval_newton_bracket(p, derivative, (lo, hi), &x, &value, &width, work)
}

// Clustered roots are conditioned by their separation, not their distance from
// zero. The derivative over the whole bracket bounds the Newton correction
// without that assumption. Exact endpoint signs still certify every new bracket.
fn interval_newton_bracket(
    p: &UPoly,
    derivative: &UPoly,
    (lo, hi): (&Rational, &Rational),
    point: &Ball,
    value: &Ball,
    width: &Rational,
    work: &mut Work,
) -> Result<Option<Narrowing>, &'static str> {
    let precision = point.precision();
    work.ball(precision, true)?;
    let range = Ball::bounds(lo.clone(), hi.clone(), precision);
    let slope = eval_ball_checked(derivative, &range, work)?;
    work.ball(precision, true)?;
    work.ball(precision, true)?;
    let Some(correction) = value.checked_div(&slope) else {
        return Ok(None);
    };
    work.ball(precision, false)?;
    let enclosure = point - &correction;
    let a = enclosure.lo().max(lo.clone());
    let b = enclosure.hi().min(hi.clone());
    work.exact(&a, &b, false)?;
    if a >= b || (&b - &a).mul_pow2(1) > *width {
        return Ok(None);
    }
    let sa = filtered_sign(p, &a, work)?;
    if sa == Sign::Zero {
        return Ok(Some(Narrowing::Exact(a)));
    }
    let sb = filtered_sign(p, &b, work)?;
    if sb == Sign::Zero {
        return Ok(Some(Narrowing::Exact(b)));
    }
    Ok((sa != sb).then_some(Narrowing::Bracket(a, b)))
}

pub(crate) fn eval_ball_checked(
    p: &UPoly,
    x: &Ball,
    work: &mut Work,
) -> Result<Ball, &'static str> {
    let precision = x.precision();
    let mut value = Ball::from_rational(&Rational::zero(), precision);
    for c in p.coefficients().iter().rev() {
        work.ball(precision, true)?;
        let c = Ball::from_rational(c, precision);
        work.ball(precision, true)?;
        value = &value * x;
        work.ball(precision, false)?;
        value = &value + &c;
    }
    Ok(value)
}

// If |a_i/a_n| <= 2^(h(n-i)), then for |x| >= 2^(h+1) the
// nonleading terms have relative sum at most sum_{j=1}^n 2^-j < 1.
// This dyadic Fujiwara-style bound avoids a large Cauchy search interval.
fn root_bound(p: &UPoly, work: &mut Work) -> Result<Rational, &'static str> {
    let n = p.degree().unwrap();
    let leading = p.leading().unwrap();
    let mut h = i64::MIN;
    for (i, c) in p.coefficients().iter().enumerate().take(n) {
        if c.is_zero() {
            continue;
        }
        work.exact(c, leading, true)?;
        let ratio = (c / leading).abs();
        let e = ratio.binary_exponent() + 1;
        let d = (n - i) as i64;
        let k = e.div_euclid(d) + i64::from(e.rem_euclid(d) != 0);
        h = h.max(k);
    }
    Ok(Rational::one().mul_pow2(if h == i64::MIN { 0 } else { h + 1 }))
}

// Filtered subdivision and Taylor shifts follow the bitstream Descartes
// principle of Eigenwillig, Kettner, Krandick, Mehlhorn, Schmitt and Wolpert,
// "A Descartes Algorithm for Polynomials with Bit-Stream Coefficients", 2005.
// An uncertain coefficient is never assigned a sign: double precision, then
// use exact coefficients at the cap. Descartes gives 0 roots for V=0 and
// exactly one for V=1; larger variations only request another subdivision.
fn descartes_variations(
    p: &UPoly,
    lo: &Rational,
    hi: &Rational,
    work: &mut Work,
) -> Result<usize, &'static str> {
    let mut precision = 53u32.min(work.policy.max_precision_bits);
    while precision != 0 {
        let mut coefficients = Vec::with_capacity(p.coefficients().len());
        for c in p.coefficients() {
            work.ball(precision, true)?;
            coefficients.push(Ball::from_rational(c, precision));
        }
        work.ball(precision, true)?;
        let shift = Ball::from_rational(lo, precision);
        work.exact(hi, lo, false)?;
        work.ball(precision, true)?;
        let width = Ball::from_rational(&(hi - lo), precision);
        taylor_shift(&mut coefficients, &shift, precision, work)?;
        let mut power = Ball::from_rational(&Rational::one(), precision);
        for c in &mut coefficients {
            work.ball(precision, true)?;
            *c = &*c * &power;
            work.ball(precision, true)?;
            power = &power * &width;
        }
        coefficients.reverse();
        taylor_shift(
            &mut coefficients,
            &Ball::from_rational(&Rational::one(), precision),
            precision,
            work,
        )?;
        let signs: Option<Vec<_>> = coefficients.iter().map(ball_sign).collect();
        if let Some(signs) = signs {
            return Ok(sign_variations(signs));
        }
        if precision == 53u32.min(work.policy.max_precision_bits)
            && let Some(variations) =
                bounded_exact(work, |work| exact_descartes_variations(p, lo, hi, work))?
        {
            return Ok(variations);
        }
        if precision == work.policy.max_precision_bits {
            break;
        }
        precision = precision
            .saturating_mul(2)
            .min(work.policy.max_precision_bits);
    }
    exact_descartes_variations(p, lo, hi, work)
}
fn exact_descartes_variations(
    p: &UPoly,
    lo: &Rational,
    hi: &Rational,
    work: &mut Work,
) -> Result<usize, &'static str> {
    let mut coefficients = p.coefficients().to_vec();
    exact_taylor_shift(&mut coefficients, lo, work)?;
    work.exact(hi, lo, false)?;
    let width = hi - lo;
    let mut power = Rational::one();
    for c in &mut coefficients {
        work.exact(c, &power, true)?;
        *c = &*c * &power;
        work.exact(&power, &width, true)?;
        power = &power * &width;
    }
    coefficients.reverse();
    exact_taylor_shift(&mut coefficients, &Rational::one(), work)?;
    Ok(sign_variations(coefficients.iter().map(Rational::sign)))
}
fn taylor_shift(
    coefficients: &mut [Ball],
    shift: &Ball,
    precision: u32,
    work: &mut Work,
) -> Result<(), &'static str> {
    for i in (0..coefficients.len() - 1).rev() {
        for j in i..coefficients.len() - 1 {
            work.ball(precision, true)?;
            let term = shift * &coefficients[j + 1];
            work.ball(precision, false)?;
            coefficients[j] = &coefficients[j] + &term;
        }
    }
    Ok(())
}
fn exact_taylor_shift(
    coefficients: &mut [Rational],
    shift: &Rational,
    work: &mut Work,
) -> Result<(), &'static str> {
    for i in (0..coefficients.len() - 1).rev() {
        for j in i..coefficients.len() - 1 {
            work.exact(shift, &coefficients[j + 1], true)?;
            let term = shift * &coefficients[j + 1];
            work.exact(&coefficients[j], &term, false)?;
            coefficients[j] = &coefficients[j] + &term;
        }
    }
    Ok(())
}
fn filtered_sign(p: &UPoly, x: &Rational, work: &mut Work) -> Result<Sign, &'static str> {
    let mut precision = 53u32.min(work.policy.max_precision_bits);
    while precision != 0 {
        work.ball(precision, true)?;
        let argument = Ball::from_rational(x, precision);
        let value = eval_ball_checked(p, &argument, work)?;
        if let Some(sign) = ball_sign(&value) {
            return Ok(sign);
        }
        if precision == 53u32.min(work.policy.max_precision_bits)
            && let Some(sign) = bounded_exact(work, |work| Ok(work.eval_poly(p, x)?.sign()))?
        {
            return Ok(sign);
        }
        if precision == work.policy.max_precision_bits {
            break;
        }
        precision = precision
            .saturating_mul(2)
            .min(work.policy.max_precision_bits);
    }
    Ok(work.eval_poly(p, x)?.sign())
}

// Exact rational signs can be much cheaper than refining a ball around zero.
// Reserve only a fraction of the remaining budget for the early attempt, so
// large coefficients still leave room for the interval filter to refine.
fn bounded_exact<T>(
    work: &mut Work,
    attempt: impl FnOnce(&mut Work) -> Result<T, &'static str>,
) -> Result<Option<T>, &'static str> {
    let limit = work.policy.max_work;
    let allowance = 1_000_000.min(limit.saturating_sub(work.used) / 8);
    work.policy.max_work = work.used + allowance;
    let result = attempt(work);
    work.policy.max_work = limit;
    match result {
        Ok(value) => Ok(Some(value)),
        Err("numbers.budget-exhausted") => Ok(None),
        Err(code) => Err(code),
    }
}
fn ball_sign(value: &Ball) -> Option<Sign> {
    match value.sign() {
        Some(crate::StrictSign::Negative) => Some(Sign::Negative),
        Some(crate::StrictSign::Positive) => Some(Sign::Positive),
        None if value.is_exact_zero() => Some(Sign::Zero),
        None => None,
    }
}
fn sign_variations(signs: impl IntoIterator<Item = Sign>) -> usize {
    let mut previous = Sign::Zero;
    let mut variations = 0;
    for sign in signs {
        if sign != Sign::Zero {
            if previous != Sign::Zero && previous != sign {
                variations += 1;
            }
            previous = sign;
        }
    }
    variations
}
fn rational_root(
    p: &UPoly,
    x: &Rational,
    max_den: &Integer,
    lo: &Rational,
    hi: &Rational,
    work: &mut Work,
) -> Result<Option<Rational>, &'static str> {
    // Legendre's continued-fraction criterion: an approximation closer than
    // 1/(2 q^2) has p/q among its convergents. The rational-root theorem bounds
    // q by the primitive leading coefficient; target width is 1/(4 Q^2).
    let mut n = x.numerator().clone();
    let mut d = x.denominator().clone();
    let mut p0 = Integer::zero();
    let mut p1 = Integer::one();
    let mut q0 = Integer::one();
    let mut q1 = Integer::zero();
    while !d.is_zero() {
        work.exact(
            &Rational::from_integer(n.clone()),
            &Rational::from_integer(d.clone()),
            true,
        )?;
        let (a, r) = n.div_rem_floor(&d).unwrap();
        let pn = &(&a * &p1) + &p0;
        let qn = &(&a * &q1) + &q0;
        if qn > *max_den {
            break;
        }
        let value = Rational::new(pn.clone(), qn.clone()).unwrap();
        if value > *lo && value < *hi && work.eval_poly(p, &value)?.is_zero() {
            return Ok(Some(value));
        }
        p0 = p1;
        p1 = pn;
        q0 = q1;
        q1 = qn;
        n = d;
        d = r;
    }
    Ok(None)
}

// Complete bounded divisor enumeration is an optional rational-root-theorem
// certificate. Failure to factor or a large candidate set uses continued
// fractions instead; it never asserts that unexamined roots are irrational.
fn rational_candidates(
    coefficients: &[Integer],
    work: &mut Work,
) -> Result<Option<Vec<Rational>>, &'static str> {
    // The rational-root theorem is only a convenience before exact
    // continued-fraction isolation. Keep this optional search bounded so a
    // large coefficient cannot consume the request allowance before the
    // complete fallback gets a chance to run.
    const OPTIONAL_CANDIDATE_WORK: u64 = 100_000;
    let used_at_start = work.used;
    let ceiling = (work.policy.max_work / 16).min(OPTIONAL_CANDIDATE_WORK);
    let Some(numerators) = divisors(&coefficients[0].abs(), work, used_at_start, ceiling)? else {
        return Ok(None);
    };
    let Some(denominators) = divisors(
        &coefficients.last().unwrap().abs(),
        work,
        used_at_start,
        ceiling,
    )?
    else {
        return Ok(None);
    };
    if numerators.len().saturating_mul(denominators.len()) > 4096 {
        return Ok(None);
    }
    let mut candidates = std::collections::BTreeSet::new();
    for n in numerators {
        for d in &denominators {
            if !work.exact_within(
                &Rational::from_integer(n.clone()),
                &Rational::from_integer(d.clone()),
                true,
                used_at_start,
                ceiling,
            )? {
                return Ok(None);
            }
            let r = Rational::new(n.clone(), d.clone()).unwrap();
            candidates.insert(-&r);
            candidates.insert(r);
        }
    }
    Ok(Some(candidates.into_iter().collect()))
}
fn divisors(
    n: &Integer,
    work: &mut Work,
    used_at_start: u64,
    ceiling: u64,
) -> Result<Option<Vec<Integer>>, &'static str> {
    let mut n = n.clone();
    let mut factors = Vec::new();
    let mut trial = 2i64;
    while n > Integer::one() && trial <= 1000 {
        let d = Integer::from_i64(trial);
        if &d * &d > n {
            factors.push((n.clone(), 1u32));
            n = Integer::one();
            break;
        }
        let mut exponent = 0;
        loop {
            if !work.exact_within(
                &Rational::from_integer(n.clone()),
                &Rational::from_integer(d.clone()),
                true,
                used_at_start,
                ceiling,
            )? {
                return Ok(None);
            }
            let (q, r) = n.div_rem_floor(&d).unwrap();
            if !r.is_zero() {
                break;
            }
            n = q;
            exponent += 1;
        }
        if exponent != 0 {
            factors.push((d, exponent));
        }
        trial += 1;
    }
    if n != Integer::one() {
        return Ok(None);
    }
    let mut result = vec![Integer::one()];
    for (prime, exponent) in factors {
        let old = result.clone();
        let mut power = Integer::one();
        for _ in 0..exponent {
            power = &power * &prime;
            for d in &old {
                if result.len() >= 4096 {
                    return Ok(None);
                }
                result.push(d * &power);
            }
        }
    }
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_rational_root_search_leaves_room_for_exact_fallback() {
        let policy = crate::BudgetPolicy {
            max_precision_bits: 128,
            max_work: 10_000,
        };
        let polynomial = UPoly::new(vec![
            Rational::from_i64(2),
            Rational::from_i64(-2018),
            Rational::from_i64(-1),
            Rational::from_i64(1009),
        ]);
        let mut work = Work::new(policy, None);
        let roots = isolate(&polynomial, &mut work).unwrap();
        assert_eq!(roots.len(), 3);
        assert!(roots.iter().any(|root| root.value.as_rational()
            == Rational::new(Integer::one(), Integer::from_i64(1009)).as_ref()));
        assert!(work.used < policy.max_work);
    }

    #[test]
    fn rational_quadratic_roots_use_the_closed_form_before_divisor_search() {
        let polynomial = UPoly::new(vec![
            Rational::from_i64(2),
            Rational::from_i64(-3027),
            Rational::from_i64(1_018_081),
        ]);
        let policy = crate::BudgetPolicy {
            max_precision_bits: 128,
            max_work: 10_000,
        };
        let mut work = Work::new(policy, None);
        let roots = isolate(&polynomial, &mut work).unwrap();
        assert_eq!(
            roots
                .iter()
                .filter(|root| root.value.as_rational().is_some())
                .count(),
            2
        );
        assert!(work.used < 1_000, "closed form used {} work", work.used);
    }

    #[test]
    fn nonsquare_quadratic_roots_use_checked_closed_form_brackets() {
        let polynomial = UPoly::new(vec![
            Rational::from_i64(-2),
            Rational::zero(),
            Rational::one(),
        ]);
        let mut work = Work::new(crate::BudgetPolicy::default(), None);
        let roots = isolate(&polynomial, &mut work).unwrap();
        assert_eq!(roots.len(), 2);
        assert!(roots.iter().all(|root| root.value.as_rational().is_none()));
        let intervals: Vec<_> = roots
            .iter()
            .map(|root| root.value.isolating_interval())
            .collect();
        assert!(intervals[0].0 < intervals[0].1);
        assert!(intervals[0].1 < intervals[1].0);
        assert!(intervals[1].0 < intervals[1].1);
        assert!(work.used < 1_000, "closed form used {} work", work.used);
    }

    #[test]
    fn optional_exact_step_does_not_exhaust_tight_request_allowance() {
        let policy = crate::BudgetPolicy {
            max_precision_bits: 128,
            max_work: 100,
        };
        let mut work = Work::new(policy, None);
        let one = Rational::one();
        let (accepted, observed) = crate::with_operation_budget(4, || {
            work.exact_within(&one, &one, true, 0, 100).unwrap()
        });
        assert!(!accepted);
        assert_eq!(observed.used, 0);
        assert!(!observed.exhausted);

        let mut work = Work::new(policy, None);
        let (accepted, observed) = crate::with_operation_budget(5, || {
            work.exact_within(&one, &one, true, 0, 100).unwrap()
        });
        assert!(accepted);
        assert_eq!(observed.used, 5);
        assert!(!observed.exhausted);
    }
}
