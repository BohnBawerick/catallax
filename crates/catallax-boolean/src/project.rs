//! Rational charts. The cylindrical coordinate is a lift of a diamond,
//! not an angle or a transcendental parameter.

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{RPoint3, RVector3, Rotation, Surface};

use crate::geometry::{r, refusal, sign, unit};

#[derive(Clone)]
pub(crate) enum Projection {
    Plane {
        origin: RPoint3,
        u: RVector3,
        v: RVector3,
    },
    Cylinder {
        origin: RPoint3,
        u: RVector3,
        v: RVector3,
        axis: RVector3,
    },
    Sphere {
        origin: RPoint3,
        radius: Rational,
        u: RVector3,
        v: RVector3,
        pole: RVector3,
    },
}

pub(crate) fn ball(lo: &Rational, hi: &Rational, bits: u32) -> Ball {
    Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits))
}

fn dot(a: &[Ball; 3], b: &RVector3) -> Ball {
    let p = a[0].precision();
    &(&a[0] * &Ball::from_rational(&b.x, p))
        + &(&(&a[1] * &Ball::from_rational(&b.y, p)) + &(&a[2] * &Ball::from_rational(&b.z, p)))
}

fn displacement(p: &[Ball; 3], origin: &RPoint3) -> [Ball; 3] {
    let bits = p[0].precision();
    let c = [&origin.x, &origin.y, &origin.z];
    std::array::from_fn(|i| &p[i] - &Ball::from_rational(c[i], bits))
}

fn phase(x: &Rational, y: &Rational) -> Rational {
    if y.sign() != Sign::Negative {
        if x.sign() != Sign::Negative {
            y / &(x + y)
        } else {
            r(1) + (-x) / (-x + y.clone())
        }
    } else if x.sign() == Sign::Negative {
        r(2) + (-y) / (-x - y.clone())
    } else {
        r(3) + x / &(x - y)
    }
}

fn phase_box(x: &Ball, y: &Ball) -> Option<Ball> {
    if x.contains_zero() && y.contains_zero() {
        return None;
    }
    let mut corners = Vec::with_capacity(4);
    for a in [x.lo(), x.hi()] {
        for b in [y.lo(), y.hi()] {
            if a.is_zero() && b.is_zero() {
                return None;
            }
            corners.push(phase(&a, &b));
        }
    }
    let anchor = corners[0].clone();
    for c in &mut corners {
        if &*c - &anchor > r(2) {
            *c = &*c - &r(4);
        }
        if &anchor - &*c > r(2) {
            *c = &*c + &r(4);
        }
    }
    let lo = corners.iter().min()?.clone();
    let hi = corners.iter().max()?.clone();
    if &hi - &lo >= r(2) {
        return None;
    }
    let result = ball(&lo, &hi, x.precision());
    (result.hi() - result.lo() < r(2)).then_some(result)
}

impl Projection {
    pub(crate) fn new(surface: &Surface, choice: i64) -> Result<Self, Refusal> {
        Ok(match surface {
            Surface::Plane(s) => Self::Plane {
                origin: s.origin.clone(),
                u: s.ref_dir.clone(),
                v: s.normal.cross(&s.ref_dir),
            },
            Surface::Cylinder(s) => Self::Cylinder {
                origin: s.origin.clone(),
                u: s.ref_dir.clone(),
                v: s.axis.cross(&s.ref_dir),
                axis: s.axis.clone(),
            },
            Surface::Sphere(s) => {
                let k = r(choice);
                let rotation = Rotation::from_quaternion(
                    r(1),
                    k.clone(),
                    &k * &k,
                    k.pow(4).expect("nonnegative power"),
                )?;
                Self::Sphere {
                    origin: s.center.clone(),
                    radius: s.radius.clone(),
                    u: rotation.apply(&unit(0)),
                    v: rotation.apply(&unit(1)),
                    pole: rotation.apply(&unit(2)),
                }
            }
            _ => {
                return Err(refusal(
                    "intersect.unsupported-pair",
                    "surface chart is outside intersection coverage",
                ));
            }
        })
    }

    pub(crate) fn periodic(&self) -> bool {
        matches!(self, Self::Cylinder { .. })
    }
    pub(crate) fn pole_point(&self) -> Option<RPoint3> {
        match self {
            Self::Sphere {
                origin,
                radius,
                pole,
                ..
            } => Some(origin.plus(&pole.scaled(radius))),
            _ => None,
        }
    }

    pub(crate) fn project_certified(
        &self,
        p: &[Ball; 3],
        dec: &mut Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        let excludes_zero = |value: &Ball, dec: &mut Decider| -> Result<bool, Refusal> {
            if sign(&Real::from(value.lo()), dec, "boolean.classify-ambiguous")? == Sign::Positive {
                return Ok(true);
            }
            Ok(sign(&Real::from(value.hi()), dec, "boolean.classify-ambiguous")? == Sign::Negative)
        };
        Ok(match self {
            Self::Plane { origin, u, v } => {
                let d = displacement(p, origin);
                Some([dot(&d, u), dot(&d, v)])
            }
            Self::Cylinder { origin, u, v, axis } => {
                let d = displacement(p, origin);
                let x = dot(&d, u);
                let y = dot(&d, v);
                if !excludes_zero(&x, dec)? && !excludes_zero(&y, dec)? {
                    return Ok(None);
                }
                let Some(phase) = phase_box(&x, &y) else {
                    return Ok(None);
                };
                Some([phase, dot(&d, axis)])
            }
            Self::Sphere {
                origin,
                radius,
                u,
                v,
                pole,
            } => {
                let d = displacement(p, origin);
                let denominator = Ball::from_rational(radius, p[0].precision()) - dot(&d, pole);
                if !excludes_zero(&denominator, dec)? {
                    return Ok(None);
                }
                let (Some(x), Some(y)) = (
                    dot(&d, u).checked_div(&denominator),
                    (-dot(&d, v)).checked_div(&denominator),
                ) else {
                    return Ok(None);
                };
                Some([x, y])
            }
        })
    }
}

pub(crate) type Box2 = [Ball; 2];
pub(crate) type Point2 = [Rational; 2];

pub(crate) fn shift(b: &mut Box2, amount: &Rational) {
    b[0] = &b[0] + &Ball::from_rational(amount, b[0].precision());
}

pub(crate) fn align(b: &mut Box2, reference: &Box2) -> Rational {
    let turns = ((reference[0].mid() - b[0].mid() + r(2)) / r(4)).floor();
    let amount = Rational::from_integer(turns) * r(4);
    shift(b, &amount);
    amount
}

pub(crate) fn separated(a: &Box2, b: &Box2, dec: &mut Decider) -> Result<bool, Refusal> {
    for i in 0..2 {
        for gap in [b[i].lo() - a[i].hi(), a[i].lo() - b[i].hi()] {
            if sign(&Real::from(gap), dec, "boolean.classify-ambiguous")? == Sign::Positive {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

pub(crate) fn common_point(boxes: &[&Box2]) -> Option<Point2> {
    let mut p = [Rational::zero(), Rational::zero()];
    for (i, x) in p.iter_mut().enumerate() {
        let lo = boxes.iter().map(|b| b[i].lo()).max()?;
        let hi = boxes.iter().map(|b| b[i].hi()).min()?;
        if lo > hi {
            return None;
        }
        *x = (lo + hi).mul_pow2(-1);
    }
    Some(p)
}

pub(crate) fn winding(
    segments: &[(Point2, Point2)],
    q: &Point2,
    periodic: bool,
    dec: &mut Decider,
) -> Result<i64, Refusal> {
    let mut result = 0i64;
    for (a, b) in segments {
        let (lo, hi) = if periodic {
            let min = a[0].clone().min(b[0].clone());
            let max = a[0].clone().max(b[0].clone());
            let lo = ((&min - &q[0]) / r(4)).floor().to_i64();
            let hi = ((&max - &q[0]) / r(4)).ceil().to_i64();
            (lo, hi)
        } else {
            (Some(0), Some(0))
        };
        let (Some(lo), Some(hi)) = (lo, hi) else {
            return Err(refusal(
                "boolean.classify-ambiguous",
                "cylindrical winding exceeds representable traversal",
            ));
        };
        for k in lo..=hi {
            let x = &q[0] + &(r(4) * r(k));
            let left = sign(&Real::from(&x - &a[0]), dec, "boolean.classify-ambiguous")?;
            let right = sign(&Real::from(&x - &b[0]), dec, "boolean.classify-ambiguous")?;
            let forward = left != Sign::Negative && right == Sign::Negative;
            let backward = right != Sign::Negative && left == Sign::Negative;
            if !forward && !backward {
                continue;
            }
            let cross = (&b[0] - &a[0]) * (&q[1] - &a[1]) - (&b[1] - &a[1]) * (&x - &a[0]);
            let crossing = sign(&Real::from(cross), dec, "boolean.classify-ambiguous")?;
            if crossing == Sign::Zero {
                return Err(refusal(
                    "boolean.classify-ambiguous",
                    "query meets a certified boundary polygon",
                ));
            }
            if forward && crossing == Sign::Negative {
                result -= 1;
            }
            if backward && crossing == Sign::Positive {
                result += 1;
            }
        }
    }
    Ok(result)
}
