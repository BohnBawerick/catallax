//! Integrate polynomial boundary forms on rationally framed analytic arcs.
//!
//! Substitution and antiderivatives are exact. Only the swept angle needs Ball
//! arithmetic. Unsupported radicals or variable divisors keep Taylor quadrature.

use std::collections::BTreeMap;

use catallax_shapes::{Curve, RPoint3, RVector3};

use super::{Checked, Node, Program, Work};
use crate::planar::{components, direction, sqrt};
use catallax_numbers::{Ball, Rational, Sign};

mod endpoints;
mod reduced;

// cos(theta)^2 = 1 - sin(theta)^2, so the cosine degree is at most one.
#[derive(Clone, Default)]
struct Polynomial(BTreeMap<(u32, u32), Rational>);

impl Polynomial {
    fn constant(value: Rational) -> Self {
        let mut result = Self::default();
        result.add_term((0, 0), value);
        result
    }

    fn add_term(&mut self, power: (u32, u32), value: Rational) {
        let sum = self.0.entry(power).or_insert_with(Rational::zero);
        *sum = &*sum + &value;
        if sum.is_zero() {
            self.0.remove(&power);
        }
    }

    fn add(&self, other: &Self) -> Self {
        let mut result = self.clone();
        for (power, value) in &other.0 {
            result.add_term(*power, value.clone());
        }
        result
    }

    fn scale(&self, factor: &Rational) -> Self {
        let mut result = Self::default();
        for (power, value) in &self.0 {
            result.add_term(*power, value * factor);
        }
        result
    }

    fn mul(&self, other: &Self) -> Self {
        let mut result = Self::default();
        for ((c, s), a) in &self.0 {
            for ((d, t), b) in &other.0 {
                let value = a * b;
                if c + d == 2 {
                    result.add_term((0, s + t), value.clone());
                    result.add_term((0, s + t + 2), -value);
                } else {
                    result.add_term((c + d, s + t), value);
                }
            }
        }
        result
    }

    fn as_constant(&self) -> Option<Rational> {
        if self.0.is_empty() {
            Some(Rational::zero())
        } else if self.0.len() == 1 {
            self.0.get(&(0, 0)).cloned()
        } else {
            None
        }
    }
}

struct Sweep {
    from: [Rational; 2],
    to: [Rational; 2],
}

impl Sweep {
    fn contains(&self, point: &[Rational; 2]) -> bool {
        fn ordered(a: &[Rational; 2], b: &[Rational; 2]) -> bool {
            let half = |p: &[Rational; 2]| {
                p[1] > Rational::zero() || (p[1].is_zero() && p[0] >= Rational::zero())
            };
            if half(a) != half(b) {
                half(a)
            } else {
                &a[0] * &b[1] - &a[1] * &b[0] >= Rational::zero()
            }
        }
        if ordered(&self.from, &self.to) {
            ordered(&self.from, point) && ordered(point, &self.to)
        } else {
            ordered(&self.from, point) || ordered(point, &self.to)
        }
    }

    fn positive(&self, polynomial: &Polynomial) -> bool {
        // The sphere gauge denominator is A + B cos(theta) + C sin(theta).
        // A >= sqrt(B^2+C^2) proves nonnegativity on the whole carrier. At
        // equality its sole zero must lie outside the closed integration arc.
        if polynomial
            .0
            .keys()
            .any(|p| !matches!(p, (0, 0) | (1, 0) | (0, 1)))
        {
            return false;
        }
        let get = |power| {
            polynomial
                .0
                .get(&power)
                .cloned()
                .unwrap_or_else(Rational::zero)
        };
        let a = get((0, 0));
        let b = get((1, 0));
        let c = get((0, 1));
        if a <= Rational::zero() {
            return false;
        }
        let square = &a * &a;
        let radial = &b * &b + &c * &c;
        square > radial || (square == radial && !self.contains(&[-b / a.clone(), -c / a]))
    }

    fn angle(&self, work: &mut Work) -> Checked<Ball> {
        let dot = &self.from[0] * &self.to[0] + &self.from[1] * &self.to[1];
        let cross = &self.from[0] * &self.to[1] - &self.from[1] * &self.to[0];
        let cosine = work.rational(&dot)?;
        let sine = work.rational(&cross)?;
        let angle = super::angle(&cosine, &sine, work)?;
        if dot >= Rational::zero() && cross.sign() == Sign::Negative {
            let pi = work.pi()?;
            let two = work.rational(&Rational::from_i64(2))?;
            let period = work.mul(&pi, &two)?;
            work.add(&angle, &period)
        } else {
            Ok(angle)
        }
    }

    // Integral = coefficient * swept_angle + rational endpoint term.
    fn monomial(&self, cosine: u32, sine: u32) -> [Rational; 2] {
        let [c0, s0] = &self.from;
        let [c1, s1] = &self.to;
        if cosine == 1 {
            let power = i32::try_from(sine + 1).unwrap();
            return [
                Rational::zero(),
                (s1.pow(power).unwrap() - s0.pow(power).unwrap())
                    / Rational::from_i64(i64::from(power)),
            ];
        }
        match sine {
            0 => [Rational::one(), Rational::zero()],
            1 => [Rational::zero(), c0 - c1],
            k => {
                let divisor = Rational::from_i64(i64::from(k));
                let factor = Rational::from_i64(i64::from(k - 1)) / divisor.clone();
                let [angle, endpoint] = self.monomial(0, k - 2);
                let power = i32::try_from(k - 1).unwrap();
                [
                    angle * factor.clone(),
                    endpoint * factor
                        - (&s1.pow(power).unwrap() * c1 - &s0.pow(power).unwrap() * c0) / divisor,
                ]
            }
        }
    }
}

pub(crate) struct Analytic(Form);

enum Form {
    Rational {
        sweep: Option<Sweep>,
        integrals: Vec<[Rational; 2]>,
    },
    Algebraic(endpoints::Integral),
}

impl Analytic {
    pub(crate) fn evaluate(&self, work: &mut Work) -> Checked<Vec<Ball>> {
        let (sweep, integrals) = match &self.0 {
            Form::Rational { sweep, integrals } => (sweep, integrals),
            Form::Algebraic(integral) => return integral.evaluate(work),
        };
        let angle = if integrals.iter().any(|i| !i[0].is_zero()) {
            Some(sweep.as_ref().unwrap().angle(work)?)
        } else {
            None
        };
        integrals
            .iter()
            .map(|[coefficient, endpoint]| {
                let endpoint = work.rational(endpoint)?;
                if coefficient.is_zero() {
                    Ok(endpoint)
                } else {
                    let coefficient = work.rational(coefficient)?;
                    let term = work.mul(&coefficient, angle.as_ref().unwrap())?;
                    work.add(&term, &endpoint)
                }
            })
            .collect()
    }
}

impl Program {
    pub(crate) fn analytic(&self, curve: &Curve, from: &RPoint3, to: &RPoint3) -> Option<Analytic> {
        self.analytic_rational(curve, from, to, &[])
    }

    fn analytic_rational(
        &self,
        curve: &Curve,
        from: &RPoint3,
        to: &RPoint3,
        constants: &[Option<Rational>],
    ) -> Option<Analytic> {
        if from == to {
            return None;
        }
        let mut inputs = Vec::with_capacity(6);
        let sweep = match curve {
            Curve::Line(line) => {
                let origin = line.origin.as_rational()?;
                let line_direction = line.dir.as_rational()?;
                if [from, to].iter().any(|point| {
                    !point
                        .minus(&origin)
                        .cross(&line_direction)
                        .norm_squared()
                        .is_zero()
                }) {
                    return None;
                }
                let delta = to.minus(from);
                for (origin, delta) in components(from).into_iter().zip(direction(&delta)) {
                    let mut polynomial = Polynomial::constant(origin.clone());
                    polynomial.add_term((0, 1), delta.clone());
                    inputs.push(polynomial);
                }
                inputs.extend(direction(&delta).map(|d| Polynomial::constant(d.clone())));
                None
            }
            Curve::Circle(circle) => {
                let center = circle.center.as_rational()?;
                let normal = circle.normal.as_rational()?;
                let reference = circle.ref_dir.as_rational()?;
                let radius_sq = circle.radius_sq.as_rational()?;
                let u = scale(
                    &reference,
                    &sqrt(&(&radius_sq / &reference.norm_squared()))?,
                );
                let second = normal.cross(&reference);
                let v = scale(&second, &sqrt(&(&radius_sq / &second.norm_squared()))?);
                Some(trigonometric_inputs(
                    &center,
                    &u,
                    &v,
                    from,
                    to,
                    &mut inputs,
                )?)
            }
            Curve::Ellipse(ellipse) => {
                let center = ellipse.center.as_rational()?;
                let normal = ellipse.normal.as_rational()?;
                let major_dir = ellipse.major_dir.as_rational()?;
                let major_sq = ellipse.major_sq.as_rational()?;
                let minor_sq = ellipse.minor_sq.as_rational()?;
                let u = scale(&major_dir, &sqrt(&(&major_sq / &major_dir.norm_squared()))?);
                let second = normal.cross(&major_dir);
                let v = scale(&second, &sqrt(&(&minor_sq / &second.norm_squared()))?);
                Some(trigonometric_inputs(
                    &center,
                    &u,
                    &v,
                    from,
                    to,
                    &mut inputs,
                )?)
            }
            _ => return None,
        };
        let values = self.substitute(&inputs, constants, sweep.as_ref())?;
        let integrals = self
            .outputs
            .iter()
            .map(|i| {
                let mut result = [Rational::zero(), Rational::zero()];
                for ((cosine, sine), coefficient) in &values[*i].0 {
                    let integral = if let Some(sweep) = &sweep {
                        sweep.monomial(*cosine, *sine)
                    } else {
                        [
                            Rational::zero(),
                            Rational::one() / Rational::from_i64(i64::from(*sine) + 1),
                        ]
                    };
                    for (sum, value) in result.iter_mut().zip(integral) {
                        *sum = &*sum + &(coefficient * &value);
                    }
                }
                result
            })
            .collect();
        Some(Analytic(Form::Rational { sweep, integrals }))
    }

    fn substitute(
        &self,
        inputs: &[Polynomial],
        constants: &[Option<Rational>],
        sweep: Option<&Sweep>,
    ) -> Option<Vec<Polynomial>> {
        let mut values: Vec<Polynomial> = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            let value = match node {
                Node::Constant(r) => Polynomial::constant(r.clone()),
                Node::Exact(i) => Polynomial::constant(constants.get(*i)?.as_ref()?.clone()),
                Node::Input(i) => inputs.get(*i)?.clone(),
                Node::Add(a, b) => values[*a].add(&values[*b]),
                Node::Mul(a, b) => values[*a].mul(&values[*b]),
                Node::Neg(a) => values[*a].scale(&Rational::from_i64(-1)),
                Node::Div(a, b) | Node::GaugeDiv(a, b) => {
                    if let Some(divisor) = values[*b].as_constant() {
                        values[*a].scale(&divisor.recip()?)
                    } else if values[*a].0.is_empty()
                        && sweep.is_some_and(|s| s.positive(&values[*b]))
                    {
                        Polynomial::default()
                    } else {
                        return None;
                    }
                }
                Node::Sqrt(a) => Polynomial::constant(sqrt(&values[*a].as_constant()?)?),
            };
            values.push(value);
        }
        Some(values)
    }
}

fn scale(vector: &RVector3, factor: &Rational) -> RVector3 {
    RVector3::new(&vector.x * factor, &vector.y * factor, &vector.z * factor)
}

fn trigonometric_inputs(
    center: &RPoint3,
    u: &RVector3,
    v: &RVector3,
    from: &RPoint3,
    to: &RPoint3,
    inputs: &mut Vec<Polynomial>,
) -> Option<Sweep> {
    let coordinates = |point: &RPoint3| {
        let delta = point.minus(center);
        [
            delta.dot(u) / u.norm_squared(),
            delta.dot(v) / v.norm_squared(),
        ]
    };
    let sweep = Sweep {
        from: coordinates(from),
        to: coordinates(to),
    };
    if [from, to]
        .iter()
        .any(|point| !point.minus(center).dot(&u.cross(v)).is_zero())
    {
        return None;
    }
    if [&sweep.from, &sweep.to]
        .iter()
        .any(|p| &p[0] * &p[0] + &p[1] * &p[1] != Rational::one())
    {
        return None;
    }
    for ((center, u), v) in components(center)
        .into_iter()
        .zip(direction(u))
        .zip(direction(v))
    {
        let mut position = Polynomial::constant(center.clone());
        position.add_term((1, 0), u.clone());
        position.add_term((0, 1), v.clone());
        inputs.push(position);
    }
    for (u, v) in direction(u).into_iter().zip(direction(v)) {
        let mut derivative = Polynomial::default();
        derivative.add_term((0, 1), -u);
        derivative.add_term((1, 0), v.clone());
        inputs.push(derivative);
    }
    Some(sweep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface::{self, Property};
    use catallax_numbers::{BudgetPolicy, Decider};
    use catallax_shapes::{Circle, Plane, Sphere, Surface};

    fn r(value: i64) -> Rational {
        Rational::from_i64(value)
    }

    fn point(x: i64, y: i64, z: i64) -> RPoint3 {
        RPoint3::new(r(x), r(y), r(z))
    }

    fn vector(x: i64, y: i64, z: i64) -> RVector3 {
        RVector3::new(r(x), r(y), r(z))
    }

    #[test]
    fn circular_sweeps_preserve_long_arcs_across_reference_seams() {
        let surface = Surface::Plane(
            Plane::new(
                (point(0, 0, 0)).into(),
                (vector(0, 0, 1)).into(),
                (vector(1, 0, 0)).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        let program = surface::compile(
            &surface,
            &point(0, 0, 0),
            Property::Area,
            false,
            &vector(0, 0, 1).to_real(),
        )
        .unwrap()
        .boundary;
        for reference in [vector(1, 0, 0), vector(3, 4, 0), vector(-3, -4, 0)] {
            let curve = Curve::Circle(
                Circle::new(
                    (point(0, 0, 0)).into(),
                    (vector(0, 0, 1)).into(),
                    (r(4)).into(),
                    (reference).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            );
            for (from, to, multiple) in [
                (point(2, 0, 0), point(0, 2, 0), 1),
                (point(0, 2, 0), point(2, 0, 0), 3),
            ] {
                let integral = program.analytic(&curve, &from, &to).unwrap();
                let policy = BudgetPolicy::default();
                let dec = Decider::default();
                let mut work = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
                let value = integral.evaluate(&mut work).unwrap().remove(0);
                let expected = &Ball::pi(512) * &Ball::from_rational(&r(multiple), 512);
                assert!(value.lo() <= expected.lo() && value.hi() >= expected.hi());
                assert!(value.hi() - value.lo() <= Rational::parse("1e-30").unwrap());
            }
        }
    }

    #[test]
    fn zero_meridian_form_never_erases_a_gauge_pole_on_the_arc() {
        let surface = Surface::Sphere(
            Sphere::new(
                (point(0, 0, 0)).into(),
                (r(2)).into(),
                (vector(0, 0, 1)).into(),
                (vector(1, 0, 0)).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        let program = surface::compile(
            &surface,
            &point(0, 0, 0),
            Property::Area,
            false,
            &vector(0, 0, 1).to_real(),
        )
        .unwrap()
        .boundary;
        let curve = Curve::Circle(
            Circle::new(
                (point(0, 0, 0)).into(),
                (vector(1, 0, 0)).into(),
                (r(4)).into(),
                (vector(0, 1, 0)).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        // +y to +z is a quarter turn away from the south gauge pole.
        let integral = program
            .analytic(&curve, &point(0, 2, 0), &point(0, 0, 2))
            .unwrap();
        let policy = BudgetPolicy::default();
        let dec = Decider::default();
        let mut work = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
        let value = integral.evaluate(&mut work).unwrap().remove(0);
        assert!(value.lo().is_zero() && value.hi().is_zero());
        // The complementary three-quarter turn passes through the pole.
        assert!(
            program
                .analytic(&curve, &point(0, 0, 2), &point(0, 2, 0))
                .is_none()
        );
        // A pole at either closed endpoint also requires the existing fallback.
        assert!(
            program
                .analytic(&curve, &point(0, 2, 0), &point(0, 0, -2))
                .is_none()
        );
    }
}
