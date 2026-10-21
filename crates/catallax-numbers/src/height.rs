//! Resultant-free nonzero bounds from absolute logarithmic heights.

use crate::{Rational, Real, UPoly, decider::Work, real::Expr, system};
use std::{collections::HashMap, sync::Arc};

// Write h₂ for the absolute logarithmic Weil height divided by log(2).
// The product formula gives |a| >= 2^(-[Q(a):Q] h₂(a)) for a != 0.
// Also h₂(a±b) <= h₂(a)+h₂(b)+1, h₂(ab), h₂(a/b) <=
// h₂(a)+h₂(b), h₂(a^n)=n h₂(a), and h₂(sqrt(a))=h₂(a)/2.
// These follow place by place from the triangle inequality and the product
// formula; see Bombieri and Gubler, Heights in Diophantine Geometry, §1.5.
// We round UP at every step using integer bits, never floating point.
// Each summary represents x = scale*y with h₂(y) <= bits. Extracting rational
// scales before propagating heights avoids charging a homothety at every
// arithmetic node. The final nonzero bound is |scale|*2^(-degree*bits).
// A field containing y is obtained by adjoining each distinct algebraic leaf
// and sqrt node once. Tower degrees multiply, including redundant generators.
// System coordinates use actual elimination polynomials, not an assumed
// zero-dimensional Bézout bound. If elimination fails, this method declines.
// Each visit costs one work unit, including memoized visits. Rational scale
// arithmetic uses Work::exact; polynomial leaves reserve Work::poly.
#[derive(Clone)]
struct Height {
    bits: u64,
    scale: Rational,
}
struct Bound<'a> {
    work: &'a mut Work,
    cache: HashMap<*const crate::real::Node, Option<Height>>,
    degree: u64,
}
fn rational_bits(r: &Rational) -> u64 {
    if r.is_zero() || (r.numerator().bit_length() == 1 && r.denominator().bit_length() == 1) {
        return 0;
    }
    r.numerator().bit_length().max(r.denominator().bit_length())
}
pub(crate) fn separation(x: &Real, work: &mut Work) -> Result<Option<Rational>, &'static str> {
    let mut bound = Bound {
        work,
        cache: HashMap::new(),
        degree: 1,
    };
    let Some(height) = bound.visit(x)? else {
        return Ok(None);
    };
    let bits = height.bits.saturating_mul(bound.degree);
    // An unusable bound is not a proof and must not allocate an enormous integer.
    if bits > u64::from(bound.work.policy.max_precision_bits) {
        return Ok(None);
    }
    bound.work.charge(bits.div_ceil(64).max(1))?;
    let dyadic = Rational::one().mul_pow2(-(bits as i64));
    bound.work.exact(&height.scale, &dyadic, true)?;
    Ok(Some(if height.scale.is_zero() {
        Rational::one()
    } else {
        &height.scale.abs() * &dyadic
    }))
}
impl Bound<'_> {
    fn polynomial(&mut self, p: &UPoly) -> Option<Height> {
        let degree = p.degree()? as u64;
        if degree == 0 {
            return None;
        }
        // M(P) <= ||P||₁ for an integer polynomial. The minimal polynomial
        // divides primitive(P) over Z and integer factors have measure >=1.
        // Clearing denominators by their product gives a valid upper bound
        // without a content GCD or an assumption that P is irreducible.
        let denominator_bits = p
            .coefficients()
            .iter()
            .map(|c| {
                if c.denominator().bit_length() == 1 {
                    0
                } else {
                    c.denominator().bit_length()
                }
            })
            .fold(0u64, u64::saturating_add);
        let numerator_bits = p
            .coefficients()
            .iter()
            .map(|c| c.numerator().bit_length())
            .max()
            .unwrap_or(0);
        self.degree = self.degree.saturating_mul(degree);
        Some(Height {
            bits: denominator_bits
                .saturating_add(numerator_bits)
                .saturating_add((p.coefficients().len() as u64).ilog2() as u64 + 1),
            scale: Rational::one(),
        })
    }
    fn visit(&mut self, x: &Real) -> Result<Option<Height>, &'static str> {
        self.work.charge(1)?;
        let key = Arc::as_ptr(&x.node);
        if let Some(h) = self.cache.get(&key) {
            return Ok(h.clone());
        }
        let value = if let Some(r) = &x.node.rational {
            self.work.exact(r, r, false)?;
            Some(Height {
                bits: 0,
                scale: r.clone(),
            })
        } else {
            match &x.node.expr {
                Expr::Rational(_) => unreachable!(),
                Expr::Algebraic(a) => {
                    let p = a.defining_polynomial();
                    self.work.poly(&p, &p)?;
                    self.polynomial(&p)
                }
                Expr::Coordinate(root, axis) => {
                    match system::coordinate_polynomial(root.system(), *axis, self.work)? {
                        Some(p) => self.polynomial(&p),
                        None => None,
                    }
                }
                Expr::Neg(a) => self.visit(a)?.map(|h| Height {
                    bits: h.bits,
                    scale: -h.scale,
                }),
                Expr::Square(a) | Expr::Pow(a, _) => {
                    if let Some(h) = self.visit(a)? {
                        let n = match &x.node.expr {
                            Expr::Pow(_, n) => *n,
                            _ => 2,
                        };
                        let mut scale = Rational::one();
                        let mut factor = h.scale;
                        let mut exponent = n;
                        while exponent != 0 {
                            if exponent & 1 != 0 {
                                self.work.exact(&scale, &factor, true)?;
                                scale = &scale * &factor;
                            }
                            exponent >>= 1;
                            if exponent != 0 {
                                self.work.exact(&factor, &factor, true)?;
                                factor = &factor * &factor;
                            }
                        }
                        Some(Height {
                            bits: h.bits.saturating_mul(u64::from(n)),
                            scale,
                        })
                    } else {
                        None
                    }
                }
                Expr::Sqrt(a) => {
                    let h = self.visit(a)?;
                    self.degree = self.degree.saturating_mul(2);
                    if let Some(h) = h {
                        match crate::quadratic::rational_sqrt(&h.scale.abs(), self.work)? {
                            Some(scale) => Some(Height {
                                bits: h.bits.div_ceil(2),
                                scale,
                            }),
                            None => Some(Height {
                                bits: h.bits.saturating_add(rational_bits(&h.scale)).div_ceil(2),
                                scale: Rational::one(),
                            }),
                        }
                    } else {
                        None
                    }
                }
                Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                    let a = self.visit(a)?;
                    let b = self.visit(b)?;
                    if let Some((a, b)) = a.zip(b) {
                        self.work.exact(&a.scale, &b.scale, true)?;
                        match &x.node.expr {
                            Expr::Add(..) | Expr::Sub(..) => {
                                if a.scale.is_zero() {
                                    Some(b)
                                } else if b.scale.is_zero() {
                                    Some(a)
                                } else {
                                    let ratio = &b.scale / &a.scale;
                                    Some(Height {
                                        bits: a
                                            .bits
                                            .saturating_add(b.bits)
                                            .saturating_add(rational_bits(&ratio))
                                            .saturating_add(1),
                                        scale: a.scale,
                                    })
                                }
                            }
                            Expr::Mul(..) => Some(Height {
                                bits: a.bits.saturating_add(b.bits),
                                scale: &a.scale * &b.scale,
                            }),
                            _ if !b.scale.is_zero() => Some(Height {
                                bits: a.bits.saturating_add(b.bits),
                                scale: &a.scale / &b.scale,
                            }),
                            _ => None,
                        }
                    } else {
                        None
                    }
                }
            }
        };
        self.cache.insert(key, value.clone());
        Ok(value)
    }
}
