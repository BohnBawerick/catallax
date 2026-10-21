//! Correlated polynomial arithmetic in one algebraic generator.

use crate::{
    Algebraic, MPoly, Rational, Real, UPoly, algebraic::AlgebraicRepr, decider::Work, real::Expr,
};
use std::{collections::HashMap, sync::Arc};

#[derive(Clone)]
struct Value<'a> {
    generator: Option<&'a Algebraic>,
    polynomial: UPoly,
}
fn defining(a: &Algebraic) -> &UPoly {
    match &a.repr {
        AlgebraicRepr::Root(root) => &root.polynomial,
        AlgebraicRepr::Rational(_) => unreachable!("rational leaves have no generator"),
    }
}
fn same(a: &Algebraic, b: &Algebraic) -> bool {
    match (&a.repr, &b.repr) {
        (AlgebraicRepr::Root(a), AlgebraicRepr::Root(b)) => {
            Arc::ptr_eq(a, b) || (a.polynomial == b.polynomial && a.lo == b.lo && a.hi == b.hi)
        }
        _ => false,
    }
}
fn remainder(p: UPoly, modulus: &UPoly, work: &mut Work) -> Result<UPoly, &'static str> {
    if p.degree() < modulus.degree() {
        return Ok(p);
    }
    work.poly(&p, modulus)?;
    Ok(p.div_rem_checked(modulus, &|| work.check_deadline())?
        .unwrap()
        .1)
}
fn inverse(p: &UPoly, modulus: &UPoly, work: &mut Work) -> Result<Option<UPoly>, &'static str> {
    let (mut a, mut b) = (modulus.clone(), p.clone());
    let (mut s, mut t) = (UPoly::zero(), UPoly::constant(Rational::one()));
    while !b.is_zero() {
        work.poly(&a, &b)?;
        let (q, r) = a.div_rem_checked(&b, &|| work.check_deadline())?.unwrap();
        work.poly(&q, &t)?;
        let qt = &q * &t;
        work.poly(&s, &qt)?;
        let next = &s - &qt;
        (a, b) = (b, r);
        (s, t) = (t, next);
    }
    // The defining polynomial need not be irreducible. Invert only units of
    // Q[t]/(P): gcd(p,P) must be a nonzero constant. This gives an identity
    // valid at EVERY root of P, without assuming which factor is minimal.
    if a.degree() != Some(0) {
        return Ok(None);
    }
    work.poly(&s, &a)?;
    Ok(Some(remainder(
        s.scale(&a.coefficients()[0].recip().unwrap()),
        modulus,
        work,
    )?))
}
impl<'a> Value<'a> {
    fn constant(r: Rational) -> Self {
        Self {
            generator: None,
            polynomial: UPoly::constant(r),
        }
    }
    fn combine(&self, rhs: &Self, op: u8, work: &mut Work) -> Result<Option<Self>, &'static str> {
        let generator = match (self.generator, rhs.generator) {
            (Some(a), Some(b)) if !same(a, b) => return Ok(None),
            (Some(a), _) | (_, Some(a)) => Some(a),
            _ => None,
        };
        let a = &self.polynomial;
        let b = &rhs.polynomial;
        work.poly(a, b)?;
        let p = match op {
            0 => a + b,
            1 => a - b,
            2 => a * b,
            _ => {
                let reciprocal = if let Some(g) = generator {
                    let Some(v) = inverse(b, defining(g), work)? else {
                        return Ok(None);
                    };
                    v
                } else {
                    let Some(r) = b.coefficients().first().and_then(Rational::recip) else {
                        return Ok(None);
                    };
                    UPoly::constant(r)
                };
                work.poly(a, &reciprocal)?;
                a * &reciprocal
            }
        };
        let polynomial = if let Some(g) = generator {
            remainder(p, defining(g), work)?
        } else {
            p
        };
        Ok(Some(Self {
            generator,
            polynomial,
        }))
    }
}
fn visit<'a>(
    x: &'a Real,
    work: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Option<Value<'a>>>,
) -> Result<Option<Value<'a>>, &'static str> {
    work.charge(1)?;
    let key = Arc::as_ptr(&x.node);
    if let Some(v) = cache.get(&key) {
        return Ok(v.clone());
    }
    let result = if let Some(r) = &x.node.rational {
        work.exact(r, r, false)?;
        Some(Value::constant(r.clone()))
    } else {
        match &x.node.expr {
            Expr::Rational(_) => unreachable!(),
            Expr::Algebraic(a) => Some(Value {
                generator: Some(a),
                polynomial: UPoly::x(),
            }),
            Expr::Coordinate(..) => None,
            Expr::Neg(a) => visit(a, work, cache)?.map(|v| Value {
                generator: v.generator,
                polynomial: -&v.polynomial,
            }),
            Expr::Sqrt(a) => match visit(a, work, cache)? {
                Some(a) if a.polynomial.degree().is_none_or(|d| d == 0) => {
                    let r = a
                        .polynomial
                        .coefficients()
                        .first()
                        .cloned()
                        .unwrap_or_default();
                    crate::quadratic::rational_sqrt(&r, work)?.map(Value::constant)
                }
                _ => None,
            },
            Expr::Square(a) | Expr::Pow(a, _) => {
                if let Some(mut a) = visit(a, work, cache)? {
                    let mut n = match &x.node.expr {
                        Expr::Pow(_, n) => *n,
                        _ => 2,
                    };
                    let mut result = Value::constant(Rational::one());
                    while n != 0 {
                        if n & 1 != 0 {
                            result = result.combine(&a, 2, work)?.unwrap();
                        }
                        n >>= 1;
                        if n != 0 {
                            a = a.combine(&a, 2, work)?.unwrap();
                        }
                    }
                    Some(result)
                } else {
                    None
                }
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                let a = visit(a, work, cache)?;
                let b = visit(b, work, cache)?;
                match a.zip(b) {
                    Some((a, b)) => a.combine(
                        &b,
                        match &x.node.expr {
                            Expr::Add(..) => 0,
                            Expr::Sub(..) => 1,
                            Expr::Mul(..) => 2,
                            _ => 3,
                        },
                        work,
                    )?,
                    None => None,
                }
            }
        }
    };
    cache.insert(key, result.clone());
    Ok(result)
}

// Reduction modulo any P(alpha)=0 proves a ZERO remainder is zero. A nonzero
// remainder does not prove nonzero unless P is minimal, so we never make that
// inference. A constant remainder is an exact rational identity at every root.
// All callers establish definedness before using these identities. Caches are
// confined to this call. Polynomial arithmetic reserves Work::poly before every
// operation, using the existing deterministic coefficient/limb cost model.
// Every visit, including a memoized visit, additionally costs one work unit.
pub(crate) fn rational(x: &Real, work: &mut Work) -> Result<Option<Rational>, &'static str> {
    Ok(
        visit(x, work, &mut HashMap::new())?.and_then(|v| match v.polynomial.degree() {
            None => Some(Rational::zero()),
            Some(0) => Some(v.polynomial.coefficients()[0].clone()),
            _ => None,
        }),
    )
}
pub(crate) fn annihilator(x: &Real, work: &mut Work) -> Result<Option<UPoly>, &'static str> {
    let Some(value) = visit(x, work, &mut HashMap::new())? else {
        return Ok(None);
    };
    if value.polynomial.degree().is_none_or(|d| d == 0) {
        let r = value
            .polynomial
            .coefficients()
            .first()
            .cloned()
            .unwrap_or_default();
        return Ok(Some(UPoly::new(vec![-r, Rational::one()])));
    }
    let Some(generator) = value.generator else {
        return Ok(None);
    };
    let lift = |p: &UPoly| {
        MPoly::from_terms(
            2,
            p.coefficients()
                .iter()
                .enumerate()
                .map(|(i, c)| (vec![0, i as u32], c.clone()))
                .collect(),
        )
    };
    let f = lift(defining(generator));
    let g = &MPoly::var(2, 0) - &lift(&value.polynomial);
    work.mpoly(&f, &g)?;
    // Res_t(P(t), z-V(t)) vanishes at z=V(alpha). Its degree is at most
    // deg(P), retaining all uses of alpha as one correlated generator rather
    // than multiplying independent conjugates at each arithmetic node.
    Ok(f.resultant_checked(&g, 1, &|| work.check_deadline())?
        .to_upoly(0))
}
