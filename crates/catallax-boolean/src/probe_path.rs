//! A probe is adjacent to a boundary arc only after its whole chart path
//! is certified free of other boundary arcs; displacement alone proves nothing.

use crate::geometry::{natural_ends, on_arc, r, sign, unit};
use catallax_numbers::{Decider, Rational, Real, RealPoly, Refusal, Sign, StrictSign};
use catallax_shapes::{CurveId, Point3, PointId, Rotation, Surface, Vector3};
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
                // C(k,i)/C(degree,i); the identity applies at every path degree.
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
fn vector(v: &Vector3) -> [&Real; 3] {
    [&v.x, &v.y, &v.z]
}
fn point(p: &Point3) -> [&Real; 3] {
    [&p.x, &p.y, &p.z]
}
fn linear_map(origin: &Point3, terms: &[(&Vector3, P)], denominator: P) -> Path {
    let o = point(origin);
    let numerators = std::array::from_fn(|i| {
        terms
            .iter()
            .fold(denominator.scale(o[i].clone()), |p, (v, t)| {
                p.add(&t.scale(vector(v)[i].clone()))
            })
    });
    Path {
        numerators,
        denominator,
    }
}
fn parameter(point: &Point3, origin: &Point3, direction: &Vector3) -> Real {
    point.minus(origin).dot(direction) / direction.norm_squared()
}

fn paths(
    surface: &Surface,
    p: &Point3,
    q: &Point3,
    dec: &mut Decider,
) -> Result<Vec<Path>, Refusal> {
    let one = P::constant(Real::one());
    let result = match surface {
        Surface::Plane(s) => {
            let u = s.ref_dir.clone();
            let v = s.normal.cross(&u);
            let a = P::line(parameter(p, &s.origin, &u), parameter(q, &s.origin, &u));
            let b = P::line(parameter(p, &s.origin, &v), parameter(q, &s.origin, &v));
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
                        let d = point.minus(&s.center);
                        let denominator = &s.radius - &d.dot(&pole);
                        if sign(&denominator, dec, "boolean.imprint-ambiguous")? == Sign::Zero {
                            return Ok(None);
                        }
                        Ok(Some((&d.dot(&u) / &denominator, &d.dot(&v) / &denominator)))
                    };
                let (Some((a, b)), Some((c, d))) = (coordinates(p, dec)?, coordinates(q, dec)?)
                else {
                    continue;
                };
                let a = P::line(a, c);
                let b = P::line(b, d);
                let norm = a.mul(&a).add(&b.mul(&b));
                let denominator = one.add(&norm);
                let rr = s.radius.clone();
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
            let length = crate::points::root(&s.ref_dir.norm_squared());
            let base = s.ref_dir.scaled(&(&s.radius / &length));
            let mut result = Vec::new();
            for b in [base.clone(), base.scaled(&Real::from_i64(-1))] {
                let v = s.axis.cross(&b);
                let rr = b.norm_squared();
                let vv = v.norm_squared();
                let coordinates =
                    |point: &Point3, dec: &mut Decider| -> Result<Option<(Real, Real)>, Refusal> {
                        let d = point.minus(&s.origin);
                        let denominator = &rr + &d.dot(&b);
                        if sign(&denominator, dec, "boolean.imprint-ambiguous")? == Sign::Zero {
                            return Ok(None);
                        }
                        Ok(Some((
                            &d.dot(&v) / &denominator * (&rr / &vv),
                            parameter(point, &s.origin, &s.axis),
                        )))
                    };
                let (Some((a, h)), Some((c, j))) = (coordinates(p, dec)?, coordinates(q, dec)?)
                else {
                    continue;
                };
                let a = P::line(a, c);
                let h = P::line(h, j);
                let square = a.mul(&a).scale(vv);
                let denominator = P::constant(rr.clone()).add(&square);
                let radial = P::constant(rr.clone()).sub(&square);
                let cross = a.scale(rr * Real::from_i64(2));
                result.push(linear_map(
                    &s.origin,
                    &[(&b, radial), (&v, cross), (&s.axis, h.mul(&denominator))],
                    denominator,
                ));
            }
            result
        }
        Surface::Cone(s) => {
            let Some([u, v, axis]) = crate::points::analytic_frame(&s.axis, &s.ref_dir) else {
                return Ok(Vec::new());
            };
            let coordinates = |point: &Point3| {
                let d = point.minus(&s.apex);
                [d.dot(&u), d.dot(&v), d.dot(&axis)]
            };
            let a = coordinates(p);
            let b = coordinates(q);
            if sign(&a[2], dec, "boolean.imprint-ambiguous")? != Sign::Positive
                || sign(&b[2], dec, "boolean.imprint-ambiguous")? != Sign::Positive
            {
                return Ok(Vec::new());
            }
            let radius_a = &s.slope * &a[2];
            let radius_b = &s.slope * &b[2];
            let h = P::line(a[2].clone(), b[2].clone());
            let mut result = Vec::new();
            for side in [1, -1] {
                let (Some(from), Some(to)) = (
                    crate::points::stereographic(&a[0], &a[1], &radius_a, side, dec)?,
                    crate::points::stereographic(&b[0], &b[1], &radius_b, side, dec)?,
                ) else {
                    continue;
                };
                let angular = P::line(from, to);
                let square = angular.mul(&angular);
                let denominator = one.add(&square);
                let radial = h.mul(&one.sub(&square)).scale(s.slope.clone());
                let cross = h.mul(&angular).scale(&Real::from_i64(2) * &s.slope);
                result.push(linear_map(
                    &s.apex,
                    &[
                        (&u.scaled(&Real::from_i64(side)), radial),
                        (&v.scaled(&Real::from_i64(side)), cross),
                        (&axis, h.mul(&denominator)),
                    ],
                    denominator,
                ));
            }
            result
        }
        Surface::EllipticCylinder(s) => {
            let first = s
                .major_dir
                .scaled(&(&Real::one() / &crate::points::root(&s.major_dir.norm_squared())));
            let axis = s
                .axis
                .scaled(&(&Real::one() / &crate::points::root(&s.axis.norm_squared())));
            let second = axis.cross(&first);
            let major = crate::points::root(&s.major_sq);
            let minor = crate::points::root(&s.minor_sq);
            let mut result = Vec::new();
            for sense in [1, -1] {
                let first = first.scaled(&Real::from_i64(sense));
                let second = second.scaled(&Real::from_i64(sense));
                let coordinates =
                    |point: &Point3, dec: &mut Decider| -> Result<Option<(Real, Real)>, Refusal> {
                        let d = point.minus(&s.origin);
                        let denominator = &Real::one() + &(&d.dot(&first) / &major);
                        if sign(&denominator, dec, "boolean.imprint-ambiguous")? == Sign::Zero {
                            return Ok(None);
                        }
                        Ok(Some((
                            (&d.dot(&second) / &minor) / denominator,
                            parameter(point, &s.origin, &s.axis),
                        )))
                    };
                let (Some((a, h)), Some((c, j))) = (coordinates(p, dec)?, coordinates(q, dec)?)
                else {
                    continue;
                };
                let t = P::line(a, c);
                let h = P::line(h, j);
                let square = t.mul(&t);
                let denominator = one.add(&square);
                result.push(linear_map(
                    &s.origin,
                    &[
                        (&first, one.sub(&square).scale(major.clone())),
                        (&second, t.scale(Real::from_i64(2) * minor.clone())),
                        (&s.axis, h.mul(&denominator)),
                    ],
                    denominator,
                ));
            }
            result
        }
        Surface::Torus(s) => {
            let Some([u, v, axis]) = crate::points::analytic_frame(&s.axis, &s.ref_dir) else {
                return Ok(Vec::new());
            };
            let coordinates = |point: &Point3| {
                let d = point.minus(&s.center);
                let x = d.dot(&u);
                let y = d.dot(&v);
                let radius = (x.square() + y.square()).sqrt();
                [x, y, d.dot(&axis), radius]
            };
            let a = coordinates(p);
            let b = coordinates(q);
            let major = &s.major;
            let minor = &s.minor;
            let tube_a = &a[3] - major;
            let tube_b = &b[3] - major;
            let mut result = Vec::new();
            for radial_side in [1, -1] {
                let (Some(from), Some(to)) = (
                    crate::points::stereographic(&a[0], &a[1], &a[3], radial_side, dec)?,
                    crate::points::stereographic(&b[0], &b[1], &b[3], radial_side, dec)?,
                ) else {
                    continue;
                };
                let angular = P::line(from, to);
                let square = angular.mul(&angular);
                let du = one.add(&square);
                let cos_u = one.sub(&square);
                let sin_u = angular.scale(Real::from_i64(2));
                for tube_side in [1, -1] {
                    let (Some(from), Some(to)) = (
                        crate::points::stereographic(&tube_a, &a[2], minor, tube_side, dec)?,
                        crate::points::stereographic(&tube_b, &b[2], minor, tube_side, dec)?,
                    ) else {
                        continue;
                    };
                    let tube = P::line(from, to);
                    let square = tube.mul(&tube);
                    let dv = one.add(&square);
                    let signed_minor = &Real::from_i64(tube_side) * minor;
                    let radius = dv
                        .scale(major.clone())
                        .add(&one.sub(&square).scale(signed_minor.clone()));
                    let height = tube.scale(Real::from_i64(2) * signed_minor).mul(&du);
                    result.push(linear_map(
                        &s.center,
                        &[
                            (&u.scaled(&Real::from_i64(radial_side)), radius.mul(&cos_u)),
                            (&v.scaled(&Real::from_i64(radial_side)), radius.mul(&sin_u)),
                            (&axis, height),
                        ],
                        du.mul(&dv),
                    ));
                }
            }
            result
        }
        _ => Vec::new(),
    };
    Ok(result)
}

fn numerator(equation: &RealPoly, path: &Path, dec: &mut Decider) -> Result<P, Refusal> {
    let degree = equation.total_degree(dec)?.unwrap_or(0);
    let mut sum = P::constant(Real::zero());
    for (monomial, coefficient) in equation.terms() {
        let zero = if let Some(value) = coefficient.rational_ref() {
            value.is_zero()
        } else {
            coefficient_sign(&coefficient, dec)? == Sign::Zero
        };
        if zero {
            continue;
        }
        let mut term = P::constant(coefficient);
        for (i, &power) in monomial.exponents().iter().enumerate() {
            term = term.mul(&path.numerators[i].pow(power));
        }
        term = term.mul(&path.denominator.pow(degree - monomial.degree()));
        sum = sum.add(&term);
    }
    Ok(sum)
}

/// The source position is an exact arc probe on the supplied source curve.
pub(crate) fn adjacent(
    model: &Model,
    surface: &Surface,
    source: (PointId, &Point3, CurveId),
    q: &Point3,
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
                .map(|e| -> Result<P, Refusal> {
                    // Every path is a parametrisation of this surface.
                    // Re-expanding its defining zero as algebraic Real
                    // coefficients would discard that construction proof.
                    let mut polynomial = if e.same_definition(&carrier) {
                        P::constant(Real::zero())
                    } else {
                        numerator(e, &path, dec)?
                    };
                    if edge.curve == source_curve {
                        // The exact arc probe is on this carrier, and path(0)=p.
                        // Re-expanding that known zero loses shared-root
                        // information and can trigger needless elimination.
                        polynomial.0[0] = Real::zero();
                    }
                    Ok(polynomial)
                })
                .collect::<Result<_, _>>()?;
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
