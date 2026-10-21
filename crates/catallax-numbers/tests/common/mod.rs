#![allow(dead_code)]

use catallax_numbers::*;
use rug::{Float, Integer as Z, Rational as Q, float::Round};
use std::sync::Arc;

#[allow(unused_macros)]
macro_rules! blind_test {
    ($name:ident, $body:block) => {
        #[test]
        fn $name() $body
    };
}
#[allow(unused_imports)]
pub(crate) use blind_test;

pub fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
pub fn real(s: &str) -> Real {
    Real::from_rational(r(s))
}
pub fn q(n: i64, d: i64) -> Q {
    Q::from((n, d))
}
pub fn z(s: &str) -> Z {
    s.parse().unwrap()
}
pub fn eq(q: &Q) -> Rational {
    Rational::new(
        Integer::parse(&q.numer().to_string()).unwrap(),
        Integer::parse(&q.denom().to_string()).unwrap(),
    )
    .unwrap()
}
pub fn oq(r: &Rational) -> Q {
    Q::from((
        z(&r.numerator().to_string()),
        z(&r.denominator().to_string()),
    ))
}
pub fn assert_q(actual: &Rational, expected: &Q) {
    assert_eq!(actual.numerator().to_string(), expected.numer().to_string());
    assert_eq!(
        actual.denominator().to_string(),
        expected.denom().to_string()
    );
}
pub fn interval_contains(i: Interval, x: &Q) {
    assert!(!i.lo().is_nan() && !i.hi().is_nan());
    assert!(i.lo() <= i.hi());
    assert!(
        i.lo() == f64::NEG_INFINITY || Q::from_f64(i.lo()).unwrap() <= *x,
        "lower endpoint {i:?}, value {x}"
    );
    assert!(
        i.hi() == f64::INFINITY || Q::from_f64(i.hi()).unwrap() >= *x,
        "upper endpoint {i:?}, value {x}"
    );
}
pub fn ball_contains(b: &Ball, x: &Q) {
    assert!(oq(&b.rad()) >= 0);
    assert!(
        oq(&b.lo()) <= *x && *x <= oq(&b.hi()),
        "ball {b:?} excludes {x}"
    );
}
pub fn overlap(a: &Ball, b: &Ball) {
    assert!(oq(&a.lo()) <= oq(&b.hi()) && oq(&b.lo()) <= oq(&a.hi()));
}
pub fn proven_sign(d: &mut Decider, x: &Real, expected: i8) {
    let out = d.sign(x);
    match (expected, &out) {
        (-1, Outcome::Proven(StrictSign::Negative))
        | (1, Outcome::Proven(StrictSign::Positive)) => (),
        (0, Outcome::ProvenCoincident(c)) => assert_ne!(c.proof, ZeroProof::ConstructionFact),
        _ => panic!("expected sign {expected}, got {out:?}"),
    }
}
pub fn refused<T: std::fmt::Debug>(out: Outcome<T>, code: &str) {
    match out {
        Outcome::Refused(e) => assert_eq!(e.code, code),
        _ => panic!("expected {code}, got {out:?}"),
    }
}
pub fn stats_sum(s: &DecisionStats) {
    assert_eq!(
        s.decisions,
        s.filter
            + s.refined
            + s.exact
            + s.coincident_exact
            + s.coincident_separation
            + s.coincident_fact
            + s.not_exact
            + s.refused
    );
}
pub fn roots(d: &mut Decider, p: &UPoly) -> Vec<RealRoot> {
    match d.real_roots(p) {
        Outcome::Proven(v) => v,
        out => panic!("expected complete roots, got {out:?}"),
    }
}
pub fn poly(coeffs: &[i64]) -> UPoly {
    UPoly::new(coeffs.iter().map(|&x| Rational::from_i64(x)).collect())
}
pub fn chebyshev(n: usize) -> UPoly {
    let mut a = poly(&[1]);
    if n == 0 {
        return a;
    }
    let mut b = UPoly::x();
    for _ in 1..n {
        let c = &(&poly(&[0, 2]) * &b) - &a;
        a = b;
        b = c;
    }
    b
}
pub fn mignotte(n: usize) -> UPoly {
    let mut c = vec![Rational::zero(); n + 1];
    c[0] = r("-2");
    c[1] = r("40");
    c[2] = r("-200");
    c[n] = r("1");
    UPoly::new(c)
}
pub fn sphere_system() -> Arc<PolySystem> {
    // Subtracting sphere equations gives x=y=3, then z^2=25-9-9=7.
    let sphere = MPoly::from_terms(
        3,
        vec![
            (vec![2, 0, 0], r("1")),
            (vec![0, 2, 0], r("1")),
            (vec![0, 0, 2], r("1")),
            (vec![0, 0, 0], r("-25")),
        ],
    );
    let sx =
        &sphere + &MPoly::from_terms(3, vec![(vec![1, 0, 0], r("-12")), (vec![0, 0, 0], r("36"))]);
    let sy =
        &sphere + &MPoly::from_terms(3, vec![(vec![0, 1, 0], r("-12")), (vec![0, 0, 0], r("36"))]);
    Arc::new(PolySystem::new(vec![sphere, sx, sy]).unwrap())
}
pub fn sphere_box(positive: bool) -> RootBox {
    RootBox::new(
        vec![r("5/2"), r("5/2"), r(if positive { "2" } else { "-3" })],
        vec![r("7/2"), r("7/2"), r(if positive { "3" } else { "-2" })],
    )
    .unwrap()
}
pub fn hidden_zeros() -> Vec<Real> {
    // sqrt(2)^2=2. Both principal roots are nonnegative.
    // (sqrt(2)+sqrt(3))^2=5+2sqrt(6), so the positive square roots agree.
    // (1+sqrt(2))^2=3+2sqrt(2), again with positive sides.
    // For g=(1+sqrt(5))/2, g^2=(3+sqrt(5))/2=g+1.
    let s2 = real("2").sqrt();
    let s3 = real("3").sqrt();
    let g = (real("1") + real("5").sqrt()) / real("2");
    vec![
        &s2 * &s2 - real("2"),
        &s2 + &s3 - (real("5") + real("2") * real("6").sqrt()).sqrt(),
        (real("3") + real("2") * s2.clone()).sqrt() - real("1") - s2,
        g.square() - g - real("1"),
    ]
}
pub fn near_tie() -> Real {
    real("10000000000000000000000000000000000000001").sqrt() - real("100000000000000000000")
}

// SplitMix64, fixed seed per family. No dependency on platform RNG or scheduling.
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
    pub fn rational(&mut self) -> Q {
        q(self.small(), (self.next() % 100 + 1) as i64)
    }
    pub fn wide(&mut self) -> Q {
        let mut n = Z::from(0);
        let mut d = Z::from(0);
        for _ in 0..4 {
            n <<= 64;
            n += self.next();
            d <<= 64;
            d += self.next();
        }
        if d == 0 {
            d = Z::from(1);
        }
        if self.next() & 1 != 0 {
            n = -n;
        }
        Q::from((n, d))
    }
}

// MPFR is a test-only oracle. Directed rounding certifies both ends; no
// catallax-numbers arithmetic contributes to the expected value.
#[derive(Clone, Debug)]
pub struct Oracle {
    pub lo: Float,
    pub hi: Float,
}
impl Oracle {
    pub fn rational(x: &Q, p: u32) -> Self {
        Self {
            lo: Float::with_val_round(p, x, Round::Down).0,
            hi: Float::with_val_round(p, x, Round::Up).0,
        }
    }
    pub fn pi(p: u32) -> Self {
        Self {
            lo: Float::with_val_round(p, rug::float::Constant::Pi, Round::Down).0,
            hi: Float::with_val_round(p, rug::float::Constant::Pi, Round::Up).0,
        }
    }
    pub fn unary(&self, f: fn(&mut Float, Round) -> std::cmp::Ordering, increasing: bool) -> Self {
        let mut lo = if increasing {
            self.lo.clone()
        } else {
            self.hi.clone()
        };
        let mut hi = if increasing {
            self.hi.clone()
        } else {
            self.lo.clone()
        };
        f(&mut lo, Round::Down);
        f(&mut hi, Round::Up);
        Self { lo, hi }
    }
    pub fn add(&self, b: &Self) -> Self {
        let p = self.lo.prec();
        Self {
            lo: Float::with_val_round(p, &self.lo + &b.lo, Round::Down).0,
            hi: Float::with_val_round(p, &self.hi + &b.hi, Round::Up).0,
        }
    }
    pub fn neg(&self) -> Self {
        Self {
            lo: -self.hi.clone(),
            hi: -self.lo.clone(),
        }
    }
    pub fn sub(&self, b: &Self) -> Self {
        self.add(&b.neg())
    }
    pub fn mul(&self, b: &Self) -> Self {
        let p = self.lo.prec();
        let mut lows = Vec::new();
        let mut highs = Vec::new();
        for a in [&self.lo, &self.hi] {
            for b in [&b.lo, &b.hi] {
                lows.push(Float::with_val_round(p, a * b, Round::Down).0);
                highs.push(Float::with_val_round(p, a * b, Round::Up).0);
            }
        }
        Self {
            lo: lows
                .into_iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap(),
            hi: highs
                .into_iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap(),
        }
    }
    pub fn div(&self, b: &Self) -> Option<Self> {
        if b.lo <= 0 && b.hi >= 0 {
            return None;
        }
        let p = self.lo.prec();
        let inv = Self {
            lo: Float::with_val_round(p, 1 / &b.hi, Round::Down).0,
            hi: Float::with_val_round(p, 1 / &b.lo, Round::Up).0,
        };
        Some(self.mul(&inv))
    }
    pub fn sqrt(&self) -> Option<Self> {
        if self.lo < 0 {
            None
        } else {
            Some(self.unary(Float::sqrt_round, true))
        }
    }
    pub fn assert_ball(&self, b: &Ball) {
        // Sufficient containment proof, using an oracle at much higher precision.
        assert!(oq(&b.lo()) <= self.lo.to_rational().unwrap());
        assert!(oq(&b.hi()) >= self.hi.to_rational().unwrap());
    }
    pub fn assert_bound(&self, bound: &Rational) {
        let b = oq(bound);
        assert!(b >= 0);
        assert!(self.lo.to_rational().unwrap() >= -b.clone());
        assert!(self.hi.to_rational().unwrap() <= b);
    }
}

// Test-only Sturm sequence over GMP rationals. This checks the certificate's
// defining polynomial independently of the kernel's polynomial algorithms.
pub fn oracle_eval(coefficients: &[Q], x: &Q) -> Q {
    coefficients.iter().rev().fold(q(0, 1), |v, c| v * x + c)
}
fn trim(p: &mut Vec<Q>) {
    while p.last().is_some_and(|x| *x == 0) {
        p.pop();
    }
}
fn remainder(mut a: Vec<Q>, b: &[Q]) -> Vec<Q> {
    while a.len() >= b.len() && !a.is_empty() {
        let shift = a.len() - b.len();
        let factor = a.last().unwrap().clone() / b.last().unwrap();
        for (i, c) in b.iter().enumerate() {
            a[i + shift] -= factor.clone() * c;
        }
        trim(&mut a);
    }
    a
}
pub fn sturm(p: &[Q]) -> Vec<Vec<Q>> {
    let derivative: Vec<_> = p
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.clone() * i)
        .collect();
    let mut sequence = vec![p.to_vec(), derivative];
    while !sequence.last().unwrap().is_empty() {
        let n = sequence.len();
        let rem = remainder(sequence[n - 2].clone(), &sequence[n - 1]);
        sequence.push(rem.into_iter().map(|x| -x).collect());
    }
    sequence.pop();
    sequence
}
pub fn variations(sequence: &[Vec<Q>], x: &Q) -> usize {
    let signs: Vec<_> = sequence
        .iter()
        .map(|p| oracle_eval(p, x).cmp(&q(0, 1)))
        .filter(|s| *s != std::cmp::Ordering::Equal)
        .collect();
    signs.windows(2).filter(|s| s[0] != s[1]).count()
}
pub fn check_irrational_certificate(value: &Algebraic) {
    let coefficients: Vec<_> = value
        .defining_polynomial()
        .coefficients()
        .iter()
        .map(oq)
        .collect();
    assert!(coefficients.last().unwrap() > &q(0, 1));
    let mut gcd = Z::from(0);
    for c in &coefficients {
        assert_eq!(c.denom(), &Z::from(1));
        gcd = gcd.gcd(c.numer());
    }
    assert_eq!(gcd, 1);
    let sequence = sturm(&coefficients);
    assert_eq!(
        sequence.last().unwrap().len(),
        1,
        "defining polynomial must be square-free"
    );
    let (lo, hi) = value.isolating_interval();
    let (lo, hi) = (oq(&lo), oq(&hi));
    assert!(lo < hi);
    assert_ne!(oracle_eval(&coefficients, &lo), 0);
    assert_ne!(oracle_eval(&coefficients, &hi), 0);
    assert_eq!(variations(&sequence, &lo) - variations(&sequence, &hi), 1);
}
