//! Rational charts. The cylindrical coordinate is a lift of a diamond,
//! not an angle or a transcendental parameter.

use catallax_numbers::{Ball, Decider, Integer, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, Rotation, Surface, Vector3};

use crate::geometry::{r, refusal, sign, unit};

#[derive(Clone)]
pub(crate) enum Projection {
    Plane {
        origin: Point3,
        u: Vector3,
        v: Vector3,
    },
    Cylinder {
        origin: Point3,
        u: Vector3,
        v: Vector3,
        axis: Vector3,
    },
    Sphere {
        origin: Point3,
        radius: Real,
        u: Vector3,
        v: Vector3,
        pole: Vector3,
    },
}

pub(crate) fn ball(lo: &Rational, hi: &Rational, bits: u32) -> Ball {
    Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits))
}

fn dot(a: &[Ball; 3], b: &Vector3, dec: &Decider) -> Result<Ball, Refusal> {
    let p = a[0].precision();
    Ok(&(&a[0] * &dec.enclose_at(&b.x, p)?)
        + &(&(&a[1] * &dec.enclose_at(&b.y, p)?) + &(&a[2] * &dec.enclose_at(&b.z, p)?)))
}

fn displacement(p: &[Ball; 3], origin: &Point3, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    let bits = p[0].precision();
    Ok([
        &p[0] - &dec.enclose_at(&origin.x, bits)?,
        &p[1] - &dec.enclose_at(&origin.y, bits)?,
        &p[2] - &dec.enclose_at(&origin.z, bits)?,
    ])
}

// Corner coordinates are dyadic. Keep their phase homogeneous until the
// final binary enclosure; normalizing each corner only adds general GCDs.
struct Phase {
    numerator: Integer,
    denominator: Integer,
}
impl Phase {
    fn new(x: &Rational, y: &Rational) -> Self {
        let bits = x
            .denominator()
            .bit_length()
            .max(y.denominator().bit_length());
        let x = x.numerator().mul_pow2(bits - x.denominator().bit_length());
        let y = y.numerator().mul_pow2(bits - y.denominator().bit_length());
        let denominator = x.abs() + y.abs();
        let numerator = if y.sign() != Sign::Negative {
            if x.sign() != Sign::Negative {
                y
            } else {
                &denominator - &x
            }
        } else if x.sign() == Sign::Negative {
            denominator.mul_pow2(1) - y
        } else {
            &denominator * &Integer::from_i64(3) + x
        };
        Self {
            numerator,
            denominator,
        }
    }

    fn shifted(&self, amount: i64) -> Self {
        Self {
            numerator: &self.numerator + &(&self.denominator * &Integer::from_i64(amount)),
            denominator: self.denominator.clone(),
        }
    }

    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (&self.numerator * &other.denominator).cmp(&(&other.numerator * &self.denominator))
    }

    fn enclose(&self, bits: u32) -> Ball {
        Ball::from_ratio(&self.numerator, &self.denominator, bits)
            .expect("phase corners exclude the pole")
    }
}

fn phase_box(x: &Ball, y: &Ball) -> Option<Ball> {
    if x.contains_zero() && y.contains_zero() {
        return None;
    }
    // The pole guard excludes (0, 0) at every corner.
    let (x0, x1, y0, y1) = (x.lo(), x.hi(), y.lo(), y.hi());
    let mut corners = [
        Phase::new(&x0, &y0),
        Phase::new(&x0, &y1),
        Phase::new(&x1, &y0),
        Phase::new(&x1, &y1),
    ];
    let lower = corners[0].shifted(-2);
    let upper = corners[0].shifted(2);
    for c in &mut corners {
        if c.cmp(&upper).is_gt() {
            *c = c.shifted(-4);
        }
        if c.cmp(&lower).is_lt() {
            *c = c.shifted(4);
        }
    }
    let lo = corners.iter().min_by(|a, b| a.cmp(b))?;
    let hi = corners.iter().max_by(|a, b| a.cmp(b))?;
    if hi.cmp(&lo.shifted(2)).is_ge() {
        return None;
    }
    let result = lo.enclose(x.precision()).hull(&hi.enclose(x.precision()));
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
            // Radial projection is a global cone chart, including the apex.
            // Its disk orientation is opposite the nappe's axis.
            Surface::Cone(s) => Self::Plane {
                origin: s.apex.clone(),
                u: s.ref_dir.clone(),
                v: s.ref_dir.cross(&s.axis),
            },
            Surface::EllipticCylinder(s) => Self::Cylinder {
                origin: s.origin.clone(),
                u: s.major_dir.clone(),
                v: s.axis.cross(&s.major_dir),
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
            Surface::Torus(_) => {
                return Err(refusal(
                    "boolean.torus-chart-required",
                    "non-circular torus boundaries require a certified periodic chart",
                ));
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
    pub(crate) fn pole_point(&self) -> Option<Point3> {
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
                let d = displacement(p, origin, dec)?;
                Some([dot(&d, u, dec)?, dot(&d, v, dec)?])
            }
            Self::Cylinder { origin, u, v, axis } => {
                let d = displacement(p, origin, dec)?;
                let x = dot(&d, u, dec)?;
                let y = dot(&d, v, dec)?;
                if !excludes_zero(&x, dec)? && !excludes_zero(&y, dec)? {
                    return Ok(None);
                }
                let Some(phase) = phase_box(&x, &y) else {
                    return Ok(None);
                };
                Some([phase, dot(&d, axis, dec)?])
            }
            Self::Sphere {
                origin,
                radius,
                u,
                v,
                pole,
            } => {
                let d = displacement(p, origin, dec)?;
                let denominator = dec.enclose_at(radius, p[0].precision())? - dot(&d, pole, dec)?;
                if !excludes_zero(&denominator, dec)? {
                    return Ok(None);
                }
                let (Some(x), Some(y)) = (
                    dot(&d, u, dec)?.checked_div(&denominator),
                    (-dot(&d, v, dec)?).checked_div(&denominator),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diamond_projection_keeps_quadrants_seams_and_poles_distinct() {
        let epsilon = Rational::parse("1e-30").unwrap();
        for (x0, x1, y0, y1, lo, hi) in [
            ("1", "2", "1", "2", "1/3", "2/3"),
            ("-2", "-1", "1", "2", "4/3", "5/3"),
            ("-2", "-1", "-2", "-1", "7/3", "8/3"),
            ("1", "2", "-2", "-1", "10/3", "11/3"),
            ("1", "2", "-1/8", "1/8", "35/9", "37/9"),
            ("-2", "-1", "0", "0", "2", "2"),
        ] {
            let q = |s: &str| Rational::parse(s).unwrap();
            let result = phase_box(&ball(&q(x0), &q(x1), 128), &ball(&q(y0), &q(y1), 128))
                .expect("a rectangle avoiding the pole has a diamond lift");
            let (lo, hi) = (q(lo), q(hi));
            assert!(result.contains(&lo) && result.contains(&hi));
            assert!(result.lo() > &lo - &epsilon);
            assert!(result.hi() < &hi + &epsilon);
        }
        assert!(phase_box(&ball(&r(-1), &r(1), 128), &ball(&r(-1), &r(1), 128)).is_none());
    }

    #[test]
    fn diamond_point_bounds_keep_unequal_dyadic_scales() {
        for (x, y, phase) in [
            ("1", "1/8", "1/9"),
            ("-1", "1/8", "17/9"),
            ("-1", "-1/8", "19/9"),
            ("1", "-1/8", "35/9"),
            ("1", "0", "0"),
            ("0", "1", "1"),
            ("-1", "0", "2"),
            ("0", "-1", "3"),
        ] {
            for power in [-4096, 0, 4096] {
                let x = Rational::parse(x).unwrap().mul_pow2(power);
                let y = Rational::parse(y).unwrap().mul_pow2(power);
                for bits in [8, 64, 128, 512] {
                    let result = phase_box(
                        &Ball::from_rational(&x, bits),
                        &Ball::from_rational(&y, bits),
                    )
                    .unwrap();
                    let expected = Ball::from_rational(&Rational::parse(phase).unwrap(), bits);
                    assert_eq!(result.mid(), expected.mid());
                    assert_eq!(result.rad(), expected.rad());
                }
            }
        }
    }
}
