//! Polynomial carrier integrals with exact, possibly algebraic, trim endpoints.

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{Curve, Point3};

use super::{Polynomial, Program, Work, scale, sqrt, trigonometric_inputs};
use crate::integrate::Checked;
use crate::planar::{components, direction};

struct Angle {
    cosine: Real,
    sine: Real,
    wrap: bool,
}

pub(super) struct Integral {
    angle: Option<Angle>,
    values: Vec<(Rational, Real)>,
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Option<Sign>, Refusal> {
    Ok(match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Some(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Some(Sign::Negative),
        Outcome::ProvenCoincident(_) => Some(Sign::Zero),
        Outcome::NotExact(_) => None,
        Outcome::Refused(error) => return Err(error),
    })
}

fn monomial(from: &[Real; 2], to: &[Real; 2], cosine: u32, sine: u32) -> (Rational, Real) {
    let [c0, s0] = from;
    let [c1, s1] = to;
    if cosine == 1 {
        let power = sine + 1;
        return (
            Rational::zero(),
            (s1.pow_u32(power) - s0.pow_u32(power)) / Real::from_i64(i64::from(power)),
        );
    }
    match sine {
        0 => (Rational::one(), Real::zero()),
        1 => (Rational::zero(), c0 - c1),
        k => {
            let divisor = Rational::from_i64(i64::from(k));
            let factor = Rational::from_i64(i64::from(k - 1)) / divisor.clone();
            let (angle, endpoint) = monomial(from, to, 0, k - 2);
            (
                angle * factor.clone(),
                endpoint * Real::from(factor)
                    - (&s1.pow_u32(k - 1) * c1 - &s0.pow_u32(k - 1) * c0) / Real::from(divisor),
            )
        }
    }
}

impl Integral {
    pub(super) fn evaluate(&self, work: &mut Work) -> Checked<Vec<Ball>> {
        let angle = if self
            .values
            .iter()
            .any(|(coefficient, _)| !coefficient.is_zero())
        {
            let angle = self.angle.as_ref().unwrap();
            let cosine = work.real(&angle.cosine)?;
            let sine = work.real(&angle.sine)?;
            let swept = crate::integrate::angle(&cosine, &sine, work)?;
            Some(if angle.wrap {
                let pi = work.pi()?;
                let two = work.rational(&Rational::from_i64(2))?;
                let period = work.mul(&pi, &two)?;
                work.add(&swept, &period)?
            } else {
                swept
            })
        } else {
            None
        };
        self.values
            .iter()
            .map(|(coefficient, endpoint)| {
                let endpoint = work.real(endpoint)?;
                if coefficient.is_zero() {
                    Ok(endpoint)
                } else {
                    let coefficient = work.rational(coefficient)?;
                    let swept = work.mul(&coefficient, angle.as_ref().unwrap())?;
                    work.add(&swept, &endpoint)
                }
            })
            .collect()
    }
}

impl Program {
    pub(super) fn algebraic_endpoints(
        &self,
        curve: &Curve,
        from: &Point3,
        to: &Point3,
        constants: &[Option<Rational>],
        dec: &mut Decider,
    ) -> Result<Option<Integral>, Refusal> {
        let mut inputs = Vec::with_capacity(6);
        let mut angle = None;
        let (start, end) = match curve {
            Curve::Line(line) => {
                let origin = line.origin.as_rational().unwrap();
                let delta = line.dir.as_rational().unwrap();
                let offset = [from, to].map(|p| p.minus(&line.origin));
                for offset in &offset {
                    if sign(&offset.cross(&line.dir).norm_squared(), dec)? != Some(Sign::Zero) {
                        return Ok(None);
                    }
                }
                for (origin, delta) in components(&origin).into_iter().zip(direction(&delta)) {
                    let mut polynomial = Polynomial::constant(origin.clone());
                    polynomial.add_term((0, 1), delta.clone());
                    inputs.push(polynomial);
                }
                inputs.extend(direction(&delta).map(|d| Polynomial::constant(d.clone())));
                let norm = Real::from(delta.norm_squared());
                let [start, end] = offset.map(|v| &v.dot(&line.dir) / &norm);
                ([Real::zero(), start], [Real::zero(), end])
            }
            Curve::Circle(_) | Curve::Ellipse(_) => {
                let (center, normal, reference, a2, b2) = match curve {
                    Curve::Circle(c) => {
                        (&c.center, &c.normal, &c.ref_dir, &c.radius_sq, &c.radius_sq)
                    }
                    Curve::Ellipse(c) => {
                        (&c.center, &c.normal, &c.major_dir, &c.major_sq, &c.minor_sq)
                    }
                    _ => unreachable!(),
                };
                let center_r = center.as_rational().unwrap();
                let normal_r = normal.as_rational().unwrap();
                let reference_r = reference.as_rational().unwrap();
                let Some(a) = sqrt(&(a2.as_rational().unwrap() / reference_r.norm_squared()))
                else {
                    return Ok(None);
                };
                let second = normal_r.cross(&reference_r);
                let Some(bb) = sqrt(&(b2.as_rational().unwrap() / second.norm_squared())) else {
                    return Ok(None);
                };
                let u = scale(&reference_r, &a);
                let v = scale(&second, &bb);
                trigonometric_inputs(
                    &center_r,
                    &u,
                    &v,
                    &center_r.plus(&u),
                    &center_r.plus(&v),
                    &mut inputs,
                )
                .unwrap();
                let offset = [from, to].map(|p| p.minus(center));
                let uv = [u.to_real(), v.to_real()];
                let coordinates = offset.each_ref().map(|q| {
                    std::array::from_fn::<_, 2, _>(|i| q.dot(&uv[i]) / uv[i].norm_squared())
                });
                for (q, coordinates) in offset.iter().zip(&coordinates) {
                    if sign(&q.dot(normal), dec)? != Some(Sign::Zero)
                        || sign(
                            &(coordinates[0].square() + coordinates[1].square() - Real::one()),
                            dec,
                        )? != Some(Sign::Zero)
                    {
                        return Ok(None);
                    }
                }
                let [start, end] = coordinates;
                let cosine = &start[0] * &end[0] + &start[1] * &end[1];
                let sine = &start[0] * &end[1] - &start[1] * &end[0];
                let Some(cosine_sign) = sign(&cosine, dec)? else {
                    return Ok(None);
                };
                let Some(sine_sign) = sign(&sine, dec)? else {
                    return Ok(None);
                };
                if cosine_sign == Sign::Positive && sine_sign == Sign::Zero {
                    return Ok(None);
                }
                angle = Some(Angle {
                    cosine,
                    sine,
                    wrap: cosine_sign != Sign::Negative && sine_sign == Sign::Negative,
                });
                (start, end)
            }
            _ => return Ok(None),
        };
        // Variable divisors, including a sphere's gauge, do not become
        // polynomial merely because its endpoints are algebraic.
        let Some(polynomials) = self.substitute(&inputs, constants, None) else {
            return Ok(None);
        };
        let values = self
            .outputs
            .iter()
            .map(|index| {
                let mut swept = Rational::zero();
                let mut endpoint = Real::zero();
                for ((c, s), coefficient) in &polynomials[*index].0 {
                    let (a, e) = if angle.is_some() {
                        monomial(&start, &end, *c, *s)
                    } else {
                        let power = s + 1;
                        (
                            Rational::zero(),
                            (end[1].pow_u32(power) - start[1].pow_u32(power))
                                / Real::from_i64(i64::from(power)),
                        )
                    };
                    swept = swept + coefficient * &a;
                    endpoint = endpoint + Real::from(coefficient.clone()) * e;
                }
                (swept, endpoint)
            })
            .collect();
        Ok(Some(Integral { angle, values }))
    }
}
