//! A probe is adjacent to a boundary arc only after its whole chart path
//! is certified free of other boundary arcs; displacement alone proves nothing.

use crate::geometry::{natural_ends, on_arc, r, sign, unit};
use catallax_numbers::{Decider, MPoly, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{CurveId, Point3, PointId, RPoint3, RVector3, Rotation, Surface, Vector3};
use catallax_topology::{EdgeId, Model};

fn coefficient_sign(x: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if x.as_rational().is_none() && x.approx().sign().is_none() {
        // Try a modest enclosure before an algebraic zero proof. Near-contact
        // coefficients often need more than the sign filter's first 64 bits,
        // but do not need elimination. Only a zero-excluding ball decides.
        let bits = dec.policy().max_precision_bits.min(256);
        let width = Rational::one().mul_pow2(-i64::from(bits));
        if let Some(value) = dec.enclose(x, &width)?.sign() {
            return Ok(match value {
                StrictSign::Negative => Sign::Negative,
                StrictSign::Positive => Sign::Positive,
            });
        }
    }
    sign(x, dec, "boolean.imprint-ambiguous")
}

#[derive(Clone)]
struct P(Vec<Real>);
impl P {
    fn constant(x: Real) -> Self {
        Self(vec![x])
    }
    fn line(a: Real, b: Real) -> Self {
        Self(vec![a.clone(), b - a])
    }
    fn coefficient(&self, i: usize) -> Real {
        self.0.get(i).cloned().unwrap_or_else(Real::zero)
    }
    fn add(&self, b: &Self) -> Self {
        Self(
            (0..self.0.len().max(b.0.len()))
                .map(|i| self.coefficient(i) + b.coefficient(i))
                .collect(),
        )
    }
    fn neg(&self) -> Self {
        Self(self.0.iter().map(|v| -v).collect())
    }
    fn sub(&self, b: &Self) -> Self {
        self.add(&b.neg())
    }
    fn scale(&self, k: Real) -> Self {
        Self(self.0.iter().map(|x| x * &k).collect())
    }
    fn mul(&self, b: &Self) -> Self {
        let mut result = vec![Real::zero(); self.0.len() + b.0.len() - 1];
        for (i, x) in self.0.iter().enumerate() {
            for (j, y) in b.0.iter().enumerate() {
                result[i + j] = &result[i + j] + &(x * y);
            }
        }
        Self(result)
    }
    fn pow(&self, n: u32) -> Self {
        (0..n).fold(Self::constant(Real::one()), |p, _| p.mul(self))
    }
    fn is_zero(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        for x in &self.0 {
            if coefficient_sign(x, dec)? != Sign::Zero {
                return Ok(false);
            }
        }
        Ok(true)
    }
    fn excludes_zero(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut first = 0;
        let expected = loop {
            let Some(coefficient) = self.0.get(first) else {
                return Ok(false);
            };
            let value = coefficient_sign(coefficient, dec)?;
            if value != Sign::Zero {
                break value;
            }
            first += 1;
        };
        // Factor the permitted contact at t=0. On [0,1], Bernstein basis
        // functions are nonnegative. Coefficients of one weak sign, with
        // strict endpoint coefficients, therefore exclude every later hit.
        // This avoids repeatedly enclosing a deep expression tree merely
        // to prove a low-degree polynomial has a sign.
        let coefficients = &self.0[first..];
        let degree = coefficients.len() - 1;
        for k in 1..=degree {
            let mut value = coefficients[0].clone();
            let mut ratio = Rational::one();
            for (i, coefficient) in coefficients.iter().enumerate().take(k + 1).skip(1) {
                // C(k,i)/C(degree,i); supported path degrees are at most four.
                ratio = ratio * (r((k + 1 - i) as i64) / r((degree + 1 - i) as i64));
                value = &value + &(coefficient * &Real::from(ratio.clone()));
            }
            let value = coefficient_sign(&value, dec)?;
            if value != expected && (value != Sign::Zero || k == degree) {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

struct Path {
    numerators: [P; 3],
    denominator: P,
}
fn vector(v: &RVector3) -> [Rational; 3] {
    [v.x.clone(), v.y.clone(), v.z.clone()]
}
fn point(p: &RPoint3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
fn linear_map(origin: &RPoint3, terms: &[(&RVector3, P)], denominator: P) -> Path {
    let o = point(origin);
    let numerators = std::array::from_fn(|i| {
        terms
            .iter()
            .fold(denominator.scale(Real::from(o[i].clone())), |p, (v, t)| {
                p.add(&t.scale(Real::from(vector(v)[i].clone())))
            })
    });
    Path {
        numerators,
        denominator,
    }
}
fn parameter(point: &Point3, origin: &RPoint3, direction: &RVector3) -> Real {
    point.minus(&origin.to_real()).dot(&direction.to_real()) / Real::from(direction.norm_squared())
}

fn paths(
    surface: &Surface,
    p: &Point3,
    q: &RPoint3,
    dec: &mut Decider,
) -> Result<Vec<Path>, Refusal> {
    let one = P::constant(Real::one());
    let result = match surface {
        Surface::Plane(s) => {
            let u = s.ref_dir.clone();
            let v = s.normal.cross(&u);
            let a = P::line(
                parameter(p, &s.origin, &u),
                parameter(&q.to_real(), &s.origin, &u),
            );
            let b = P::line(
                parameter(p, &s.origin, &v),
                parameter(&q.to_real(), &s.origin, &v),
            );
            vec![linear_map(&s.origin, &[(&u, a), (&v, b)], one)]
        }
        Surface::Sphere(s) => {
            let mut result = Vec::new();
            for k in 0..3 {
                let rotation = Rotation::from_quaternion(r(1), r(k), r(k * k), r(k * k * k))?;
                let u = rotation.apply(&unit(0));
                let v = rotation.apply(&unit(1));
                let pole = rotation.apply(&unit(2));
                let coordinates =
                    |point: &Point3, dec: &mut Decider| -> Result<Option<(Real, Real)>, Refusal> {
                        let d = point.minus(&s.center.to_real());
                        let denominator = Real::from(s.radius.clone()) - d.dot(&pole.to_real());
                        if sign(&denominator, dec, "boolean.imprint-ambiguous")? == Sign::Zero {
                            return Ok(None);
                        }
                        Ok(Some((
                            &d.dot(&u.to_real()) / &denominator,
                            &d.dot(&v.to_real()) / &denominator,
                        )))
                    };
                let (Some((a, b)), Some((c, d))) =
                    (coordinates(p, dec)?, coordinates(&q.to_real(), dec)?)
                else {
                    continue;
                };
                let a = P::line(a, c);
                let b = P::line(b, d);
                let norm = a.mul(&a).add(&b.mul(&b));
                let denominator = one.add(&norm);
                let rr = Real::from(s.radius.clone());
                result.push(linear_map(
                    &s.center,
                    &[
                        (&u, a.scale(Real::from_i64(2) * rr.clone())),
                        (&v, b.scale(Real::from_i64(2) * rr.clone())),
                        (&pole, norm.sub(&one).scale(rr)),
                    ],
                    denominator,
                ));
            }
            result
        }
        Surface::Cylinder(s) => {
            let Some(length) = crate::points::sqrt_rational(&s.ref_dir.norm_squared()) else {
                return Ok(Vec::new());
            };
            let base = s.ref_dir.scaled(&(&s.radius / &length));
            let mut result = Vec::new();
            for b in [base.clone(), base.scaled(&r(-1))] {
                let v = s.axis.cross(&b);
                let rr = b.norm_squared();
                let vv = v.norm_squared();
                let coordinates =
                    |point: &Point3, dec: &mut Decider| -> Result<Option<(Real, Real)>, Refusal> {
                        let d = point.minus(&s.origin.to_real());
                        let denominator = Real::from(rr.clone()) + d.dot(&b.to_real());
                        if sign(&denominator, dec, "boolean.imprint-ambiguous")? == Sign::Zero {
                            return Ok(None);
                        }
                        Ok(Some((
                            &d.dot(&v.to_real()) / &denominator * Real::from(&rr / &vv),
                            parameter(point, &s.origin, &s.axis),
                        )))
                    };
                let (Some((a, h)), Some((c, j))) =
                    (coordinates(p, dec)?, coordinates(&q.to_real(), dec)?)
                else {
                    continue;
                };
                let a = P::line(a, c);
                let h = P::line(h, j);
                let square = a.mul(&a).scale(Real::from(vv));
                let denominator = P::constant(Real::from(rr.clone())).add(&square);
                let radial = P::constant(Real::from(rr.clone())).sub(&square);
                let cross = a.scale(Real::from(rr) * Real::from_i64(2));
                result.push(linear_map(
                    &s.origin,
                    &[(&b, radial), (&v, cross), (&s.axis, h.mul(&denominator))],
                    denominator,
                ));
            }
            result
        }
        _ => Vec::new(),
    };
    Ok(result)
}

fn numerator(equation: &MPoly, path: &Path) -> P {
    let degree = equation.total_degree().unwrap_or(0);
    equation
        .terms()
        .into_iter()
        .fold(P::constant(Real::zero()), |sum, (monomial, k)| {
            let mut term = P::constant(Real::from(k));
            for (i, &power) in monomial.exponents().iter().enumerate() {
                term = term.mul(&path.numerators[i].pow(power));
            }
            term = term.mul(&path.denominator.pow(degree - monomial.degree()));
            sum.add(&term)
        })
}

/// The source position is an exact arc probe on the supplied source curve.
pub(crate) fn adjacent(
    model: &Model,
    surface: &Surface,
    source: (PointId, &Point3, CurveId),
    q: &RPoint3,
    left: &Vector3,
    barriers: &[EdgeId],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (source, p, source_curve) = source;
    let carrier = surface.implicit();
    for path in paths(surface, p, q, dec)? {
        let d0 = path.denominator.coefficient(0);
        let d1 = path.denominator.coefficient(1);
        let derivatives: [Real; 3] = std::array::from_fn(|i| {
            path.numerators[i].coefficient(1) * d0.clone()
                - path.numerators[i].coefficient(0) * d1.clone()
        });
        let [x, y, z] = derivatives;
        if sign(
            &left.dot(&Vector3::new(x, y, z)),
            dec,
            "boolean.imprint-ambiguous",
        )? != Sign::Positive
        {
            continue;
        }
        let mut clear = true;
        for &edge_id in barriers {
            let edge = model.edge(edge_id);
            let equations: Vec<_> = model
                .geom()
                .curve(edge.curve)
                .equations(model.geom())
                .iter()
                .map(|e| {
                    // Every path is a parametrisation of this surface.
                    // Re-expanding its defining zero as algebraic Real
                    // coefficients would discard that construction proof.
                    let mut polynomial = if *e == carrier {
                        P::constant(Real::zero())
                    } else {
                        numerator(e, &path)
                    };
                    if edge.curve == source_curve {
                        // The exact arc probe is on this carrier, and path(0)=p.
                        // Re-expanding that known zero loses shared-root
                        // information and can trigger needless elimination.
                        polynomial.0[0] = Real::zero();
                    }
                    polynomial
                })
                .collect();
            let mut separated = false;
            for equation in &equations {
                if equation.excludes_zero(dec)? {
                    separated = true;
                    break;
                }
            }
            if !separated {
                let mut contained = true;
                for equation in &equations {
                    if !equation.is_zero(dec)? {
                        contained = false;
                        break;
                    }
                }
                let (from, to) = natural_ends(model, edge_id);
                if !contained
                    || from.is_none()
                    || on_arc(
                        model,
                        edge.curve,
                        from,
                        to,
                        source,
                        dec,
                        "boolean.imprint-ambiguous",
                    )?
                {
                    clear = false;
                    break;
                }
                let mut avoids_ends = true;
                for endpoint in [from, to].into_iter().flatten() {
                    let position = model.geom().position(endpoint);
                    let coordinates = [position.x, position.y, position.z];
                    let mut avoids = false;
                    for (axis, coordinate) in coordinates.into_iter().enumerate() {
                        if path.numerators[axis]
                            .sub(&path.denominator.scale(coordinate))
                            .excludes_zero(dec)?
                        {
                            avoids = true;
                            break;
                        }
                    }
                    if !avoids {
                        avoids_ends = false;
                        break;
                    }
                }
                if !avoids_ends {
                    clear = false;
                    break;
                }
            }
        }
        if clear {
            return Ok(true);
        }
    }
    Ok(false)
}
