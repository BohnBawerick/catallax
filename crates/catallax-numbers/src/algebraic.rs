//! Real algebraic numbers and filtered Descartes root isolation.

use crate::decider::{BudgetPolicy, Work};
use crate::{Ball, Integer, Rational, Sign, UPoly};
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};

const NEWTON_GUARD_BITS: i64 = 64;

// After removing smaller prime powers, composite trial divisors cannot divide.
const TRIAL_PRIMES: [i64; 168] = {
    let mut primes = [0; 168];
    let mut count = 0;
    let mut candidate = 2;
    while candidate <= 1000 {
        let mut divisor = 2;
        let mut prime = true;
        while divisor * divisor <= candidate {
            if candidate % divisor == 0 {
                prime = false;
                break;
            }
            divisor += 1;
        }
        if prime {
            primes[count] = candidate;
            count += 1;
        }
        candidate += 1;
    }
    assert!(count == primes.len());
    primes
};

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
    if p.degree() == Some(0) {
        work.check_deadline()?;
        return Ok(Vec::new());
    }
    if let [constant, linear] = p.coefficients() {
        work.exact(constant, linear, true)?;
        return Ok(vec![RealRoot {
            value: Algebraic::from_rational((-constant) / linear.clone()),
            multiplicity: 1,
        }]);
    }
    // A quadratic's discriminant certifies multiplicity and rational roots
    // directly, without a dense square-free GCD or divisor enumeration.
    if p.degree() == Some(2) {
        match quadratic_roots(p, work)? {
            QuadraticRoots::Nonreal => return Ok(Vec::new()),
            QuadraticRoots::Repeated(root) => {
                return Ok(vec![RealRoot {
                    value: Algebraic::from_rational(root),
                    multiplicity: 2,
                }]);
            }
            QuadraticRoots::Rational(lower, upper) => {
                return Ok(vec![
                    RealRoot {
                        value: Algebraic::from_rational(lower),
                        multiplicity: 1,
                    },
                    RealRoot {
                        value: Algebraic::from_rational(upper),
                        multiplicity: 1,
                    },
                ]);
            }
            QuadraticRoots::Irrational(_) => {}
        }
    }
    let mut roots = Vec::new();
    if p.degree() == Some(2) {
        // A nonsquare discriminant proves both roots simple without a GCD.
        roots.extend(
            isolate_square_free(p, work)?
                .into_iter()
                .map(|value| RealRoot {
                    value,
                    multiplicity: 1,
                }),
        );
    } else {
        // Rootless intersection probes need no square-free or cubic factorization.
        // An inconclusive exclusion leaves the ordinary isolator responsible.
        if p.degree() == Some(4)
            && !p.coefficients()[0].is_zero()
            && bounded_exact(work, |work| quartic_has_no_real_roots(p, work))? == Some(true)
        {
            return Ok(Vec::new());
        }
        work.poly(p, p)?;
        let factors = p.square_free_factors_checked(&|| work.check_deadline())?;
        for (f, multiplicity) in factors {
            for value in isolate_square_free(&f, work)? {
                roots.push(RealRoot {
                    value,
                    multiplicity,
                });
            }
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
enum QuadraticRoots {
    Nonreal,
    Repeated(Rational),
    Rational(Rational, Rational),
    Irrational(Rational),
}

fn quadratic_roots(polynomial: &UPoly, work: &mut Work) -> Result<QuadraticRoots, &'static str> {
    let [c, b, a] = polynomial.coefficients() else {
        unreachable!("quadratic solver requires degree two");
    };
    work.exact(b, b, true)?;
    let squared = b * b;
    work.exact(a, c, true)?;
    let product = (a * c).mul_pow2(2);
    work.exact(&squared, &product, false)?;
    let discriminant = squared - product;
    if discriminant.sign() == Sign::Negative {
        return Ok(QuadraticRoots::Nonreal);
    }
    let center = -b;
    let denominator = a.mul_pow2(1);
    if discriminant.is_zero() {
        work.exact(&center, &denominator, true)?;
        return Ok(QuadraticRoots::Repeated(&center / &denominator));
    }
    let Some(root) = crate::quadratic::rational_sqrt(&discriminant, work)? else {
        return Ok(QuadraticRoots::Irrational(discriminant));
    };
    work.exact(&center, &root, false)?;
    let lower = &center - &root;
    work.exact(&lower, &denominator, true)?;
    let mut lower = &lower / &denominator;
    work.exact(&center, &root, false)?;
    let upper = &center + &root;
    work.exact(&upper, &denominator, true)?;
    let mut upper = &upper / &denominator;
    if lower > upper {
        std::mem::swap(&mut lower, &mut upper);
    }
    Ok(QuadraticRoots::Rational(lower, upper))
}

fn guarded_quadratic(
    polynomial: &UPoly,
    target: &Rational,
    output: &mut Vec<Algebraic>,
    work: &mut Work,
) -> Result<bool, &'static str> {
    let [c, b, a] = polynomial.coefficients() else {
        unreachable!("quadratic isolator requires degree two");
    };
    let QuadraticRoots::Irrational(discriminant) = quadratic_roots(polynomial, work)? else {
        return Ok(false);
    };
    let cap = work.policy.max_precision_bits;
    if cap == 0 {
        return Ok(false);
    }
    // Midpoint and residual radius each need guard bits. Leave additional room
    // before testing that both endpoints round to the same guarded ball.
    let guard = 256.min(cap.saturating_sub(64) / 2).max(1);
    let mut precision = (2 * guard + 64).min(cap);
    loop {
        let intervals = 'attempt: {
            work.ball(precision, true)?;
            let d = Ball::from_rational(&discriminant, precision);
            work.ball(precision, true)?;
            let Some(root) = d.sqrt() else {
                break 'attempt None;
            };
            work.ball(precision, true)?;
            let center = Ball::from_rational(&(-b), precision);
            work.ball(precision, true)?;
            let denominator = Ball::from_rational(&a.mul_pow2(1), precision);
            work.ball(precision, true)?;
            let product = Ball::from_rational(&c.mul_pow2(1), precision);
            // Select the noncancelling branch, then use the exact root product
            // for the other branch. Tiny positive roots must not become zero.
            work.ball(precision, false)?;
            let numerator = if b.sign() == Sign::Negative {
                &center + &root
            } else {
                &center - &root
            };
            work.ball(precision, true)?;
            work.ball(precision, true)?;
            let Some(first) = numerator.checked_div(&denominator) else {
                break 'attempt None;
            };
            work.ball(precision, true)?;
            work.ball(precision, true)?;
            let Some(second) = product.checked_div(&numerator) else {
                break 'attempt None;
            };
            let mut intervals = [(first.lo(), first.hi()), (second.lo(), second.hi())];
            if intervals[0].0 > intervals[1].0 {
                intervals.swap(0, 1);
            }
            if intervals[0].1 >= intervals[1].0 {
                break 'attempt None;
            }
            for (lo, hi) in &intervals {
                work.exact(lo, hi, false)?;
                if lo >= hi {
                    break 'attempt None;
                }
                // A coarser established target is cheap and may carry a
                // deliberately wider export bound. Keep its original bracket.
                if hi - lo <= *target {
                    return Ok(false);
                }
                let lower = filtered_sign(polynomial, lo, work)?;
                let upper = filtered_sign(polynomial, hi, work)?;
                if lower == Sign::Zero || upper == Sign::Zero || lower == upper {
                    break 'attempt None;
                }
                // Preserve the unique outward-rounded enclosure at consumer
                // precisions, including roots arbitrarily close to dyadic ties.
                work.ball(guard, true)?;
                let lower = Ball::from_rational(lo, guard);
                work.ball(guard, true)?;
                let upper = Ball::from_rational(hi, guard);
                if lower.lo() != upper.lo() || lower.hi() != upper.hi() {
                    break 'attempt None;
                }
            }
            Some(intervals)
        };
        if let Some(intervals) = intervals {
            output
                .extend(intervals.map(|(lo, hi)| Algebraic::isolated(polynomial.clone(), lo, hi)));
            return Ok(true);
        }
        if precision == cap {
            return Ok(false);
        }
        precision = precision.saturating_mul(2).min(cap);
    }
}

// Retain the established target-width brackets for clustered roots when the
// target is already coarser than guarded construction would produce.
fn target_quadratic(
    polynomial: &UPoly,
    target: &Rational,
    output: &mut Vec<Algebraic>,
    work: &mut Work,
) -> Result<bool, &'static str> {
    let [_, b, a] = polynomial.coefficients() else {
        unreachable!("quadratic isolator requires degree two");
    };
    let QuadraticRoots::Irrational(discriminant) = quadratic_roots(polynomial, work)? else {
        return Ok(false);
    };
    let center = -b;
    let denominator = a.mul_pow2(1);
    if work.policy.max_precision_bits == 0 {
        return Ok(false);
    }
    let scale = &center.abs() + &denominator.abs();
    work.exact(&scale, &scale, true)?;
    let guard = NEWTON_GUARD_BITS.min(i64::from(work.policy.max_precision_bits));
    if discriminant >= (&scale * &scale).mul_pow2(-2 * guard) {
        return Ok(false);
    }
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
        return Ok(false);
    };
    work.ball(precision, true)?;
    let center = Ball::from_rational(&center, precision);
    work.ball(precision, true)?;
    let denominator = Ball::from_rational(&denominator, precision);
    let mut intervals = [None, None];
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
            return Ok(false);
        };
        let lo = enclosure.lo();
        let hi = enclosure.hi();
        work.exact(&lo, &hi, false)?;
        if lo >= hi || &hi - &lo >= *target {
            return Ok(false);
        }
        let lower = filtered_sign(polynomial, &lo, work)?;
        let upper = filtered_sign(polynomial, &hi, work)?;
        if lower == Sign::Zero || upper == Sign::Zero || lower == upper {
            return Ok(false);
        }
        intervals[index] = Some((lo, hi));
    }
    output.extend(intervals.map(|interval| {
        let (lo, hi) = interval.unwrap();
        Algebraic::isolated(polynomial.clone(), lo, hi)
    }));
    Ok(true)
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
    let irrational_quadratic = if remaining.degree() == Some(2) {
        match quadratic_roots(&remaining, work)? {
            QuadraticRoots::Nonreal => return Ok(result),
            QuadraticRoots::Repeated(root) => {
                result.push(Algebraic::from_rational(root));
                return Ok(result);
            }
            QuadraticRoots::Rational(lower, upper) => {
                result.extend([
                    Algebraic::from_rational(lower),
                    Algebraic::from_rational(upper),
                ]);
                return Ok(result);
            }
            QuadraticRoots::Irrational(_) => true,
        }
    } else {
        false
    };
    let integers = remaining.primitive_integer_checked(&|| work.check_deadline())?;
    // Retain the established coarse brackets when bounded divisor enumeration
    // was complete. Only formerly reconstructed quadratics need guarded bounds.
    let (candidates, enumeration_complete) = if irrational_quadratic {
        (None, rational_search_complete(&integers, work)?)
    } else {
        let candidates = rational_candidates(&integers, work)?;
        let complete = candidates.is_some();
        (candidates, complete)
    };
    let rational_complete = enumeration_complete;
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
    let integers = if irrational_quadratic {
        integers
    } else {
        remaining.primitive_integer_checked(&|| work.check_deadline())?
    };
    let primitive = UPoly::new(
        integers
            .iter()
            .cloned()
            .map(Rational::from_integer)
            .collect(),
    );
    let q = integers.last().unwrap().abs();
    let target = Rational::new(Integer::one(), (&q * &q).mul_pow2(2)).unwrap();
    if irrational_quadratic
        && !enumeration_complete
        && guarded_quadratic(&primitive, &target, &mut result, work)?
    {
        return Ok(result);
    }
    if remaining.degree() == Some(2)
        && !rational_complete
        && target_quadratic(&primitive, &target, &mut result, work)?
    {
        return Ok(result);
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
                    if !irrational_quadratic
                        && let Some(r) = rational_root(&primitive, &mid, &q, &lo, &hi, work)?
                    {
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

// Descartes' zero-variation certificate also excludes repeated real roots.
// The caller's numeric suballowance bounds this optional filter.
fn quartic_has_no_real_roots(p: &UPoly, work: &mut Work) -> Result<bool, &'static str> {
    let bound = root_bound(p, work)?;
    let mut pending = vec![(-&bound, Rational::zero()), (Rational::zero(), bound)];
    let mut splits = 0;
    while let Some((lo, hi)) = pending.pop() {
        let variations = descartes_variations(p, &lo, &hi, work)?;
        if variations == 0 {
            continue;
        }
        if variations == 1 || splits == 16 {
            return Ok(false);
        }
        work.exact(&lo, &hi, false)?;
        let mid = (&lo + &hi).mul_pow2(-1);
        if filtered_sign(p, &mid, work)? == Sign::Zero {
            return Ok(false);
        }
        splits += 1;
        pending.push((mid.clone(), hi));
        pending.push((lo, mid));
    }
    Ok(true)
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
    let Some(numerators) = divisors(&coefficients[0].abs(), work)? else {
        return Ok(None);
    };
    let Some(denominators) = divisors(&coefficients.last().unwrap().abs(), work)? else {
        return Ok(None);
    };
    if numerators.len().saturating_mul(denominators.len()) > 4096 {
        return Ok(None);
    }
    let mut candidates = std::collections::BTreeSet::new();
    for n in numerators {
        for d in &denominators {
            work.exact(
                &Rational::from_integer(n.clone()),
                &Rational::from_integer(d.clone()),
                true,
            )?;
            let r = Rational::new(n.clone(), d.clone()).unwrap();
            candidates.insert(-&r);
            candidates.insert(r);
        }
    }
    Ok(Some(candidates.into_iter().collect()))
}
fn small_prime_factors(
    n: &Integer,
    work: &mut Work,
    mut factor: impl FnMut(Integer, u32),
) -> Result<bool, &'static str> {
    if n.bit_length() <= 1 {
        return Ok(!n.is_zero());
    }
    let mut n = n.clone();
    for &trial in &TRIAL_PRIMES {
        if n.bit_length() == 1 {
            return Ok(true);
        }
        if n.to_i64().is_some_and(|value| value < trial * trial) {
            factor(n, 1);
            return Ok(true);
        }
        let d = Integer::from_i64(trial);
        let mut exponent = 0;
        loop {
            // The same integer/word division charge as exact rational arithmetic,
            // without cloning a large integer only to compute its cost.
            let limbs = n.bit_length().div_ceil(64).max(1).saturating_add(1);
            work.charge(limbs.saturating_mul(2))?;
            let (q, r) = n.div_rem_floor(&d).unwrap();
            if !r.is_zero() {
                break;
            }
            n = q;
            exponent += 1;
        }
        if exponent != 0 {
            factor(d, exponent);
        }
    }
    if n.bit_length() > 1 {
        // All primes <= 1000 were removed. A remaining value <= 1000² is prime.
        // This retains the old trial loop's final composite-index square check.
        if n.to_i64().is_none_or(|value| value > 1_000_000) {
            return Ok(false);
        }
        factor(n, 1);
    }
    Ok(true)
}

fn rational_search_complete(
    coefficients: &[Integer],
    work: &mut Work,
) -> Result<bool, &'static str> {
    let mut candidates = 1usize;
    for coefficient in [&coefficients[0], coefficients.last().unwrap()] {
        let complete = small_prime_factors(&coefficient.abs(), work, |_, exponent| {
            candidates = candidates.saturating_mul((exponent as usize).saturating_add(1));
        })?;
        if !complete || candidates > 4096 {
            return Ok(false);
        }
    }
    Ok(true)
}

fn divisors(n: &Integer, work: &mut Work) -> Result<Option<Vec<Integer>>, &'static str> {
    let mut factors = Vec::new();
    if !small_prime_factors(n, work, |prime, exponent| factors.push((prime, exponent)))? {
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
