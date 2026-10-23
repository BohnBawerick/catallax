//! Dyadic turns make translated periodic chart copies bit-identical.

use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal};
use catallax_shapes::{Surface, Vector3};

use crate::failure;

pub(crate) const BITS: u32 = 128;
pub(crate) type Uv = [f64; 2];

pub(crate) fn q(x: f64) -> Rational {
    Rational::from_f64_exact(x).expect("finite mesh coordinate")
}

pub(crate) fn midpoint(b: &Ball) -> f64 {
    b.mid().to_f64_nearest()
}

pub(crate) fn position(b: &[Ball; 3]) -> Result<[f64; 3], Refusal> {
    let p = b.each_ref().map(midpoint);
    if p.iter().all(|x| x.is_finite()) {
        Ok(p.map(|x| if x == 0.0 { 0.0 } else { x }))
    } else {
        Err(failure("geometry exceeds finite binary64 coordinates"))
    }
}

pub(crate) fn residual(p: &[f64; 3], b: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        sum + (q(p[i]) - b[i].lo()).abs().max((q(p[i]) - b[i].hi()).abs())
    })
}

fn unit(v: &Vector3, dec: &Decider) -> Result<[f64; 3], Refusal> {
    if let Some(v) = v.rational_components() {
        let norm = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let length = Ball::from_rational(&norm, BITS)
            .sqrt()
            .expect("nonzero definition axis");
        return Ok(v.map(|x| midpoint(&Ball::from_rational(x, BITS).checked_div(&length).unwrap())));
    }
    let length = v.norm_squared().sqrt();
    Ok([
        midpoint(&dec.enclose_at(&(&v.x / &length), BITS)?),
        midpoint(&dec.enclose_at(&(&v.y / &length), BITS)?),
        midpoint(&dec.enclose_at(&(&v.z / &length), BITS)?),
    ])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A chart uses turns for angles, and millimetres for linear parameters.
pub(crate) struct Chart {
    pub(crate) surface: Surface,
    pub(crate) periodic: [bool; 2],
    origin: [f64; 3],
    frame: [[f64; 3]; 3],
    turn: Option<Ball>,
    evaluator: Option<crate::evaluate::Evaluator>,
    ellipse_radii: Option<[f64; 2]>,
}

impl Chart {
    pub(crate) fn new(surface: &Surface, dec: &Decider) -> Result<Self, Refusal> {
        let (origin, axis, reference, periodic) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal, &s.ref_dir, [false, false]),
            Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir, [true, false]),
            Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir, [true, false]),
            Surface::Sphere(s) => (&s.center, &s.axis, &s.ref_dir, [true, false]),
            Surface::Torus(s) => (&s.center, &s.axis, &s.ref_dir, [true, true]),
            Surface::EllipticCylinder(s) => (&s.origin, &s.axis, &s.major_dir, [true, false]),
            _ => return Err(failure("unknown analytic surface definition")),
        };
        let scalar = |value: &Real| -> Result<f64, Refusal> {
            if let Some(value) = value.rational_ref() {
                Ok(value.to_f64_nearest())
            } else {
                Ok(midpoint(&dec.enclose_at(value, BITS)?))
            }
        };
        let origin = [scalar(&origin.x)?, scalar(&origin.y)?, scalar(&origin.z)?];
        if !origin.iter().all(|v| v.is_finite()) {
            return Err(failure("surface origin exceeds binary64"));
        }
        let turn = if matches!(surface, Surface::Plane(_)) {
            None
        } else {
            Some(Ball::pi_checked(BITS, dec)? * Ball::from_rational(&Rational::from_i64(2), BITS))
        };
        Ok(Self {
            surface: surface.clone(),
            periodic,
            origin,
            turn,
            evaluator: if matches!(surface, Surface::Plane(_)) {
                None
            } else {
                Some(crate::evaluate::Evaluator::new(surface, dec)?)
            },
            ellipse_radii: if let Surface::EllipticCylinder(s) = surface {
                Some([
                    midpoint(&dec.enclose_at(&s.major_sq.sqrt(), BITS)?),
                    midpoint(&dec.enclose_at(&s.minor_sq.sqrt(), BITS)?),
                ])
            } else {
                None
            },
            frame: [
                unit(reference, dec)?,
                unit(&axis.cross(reference), dec)?,
                unit(axis, dec)?,
            ],
        })
    }

    pub(crate) fn project(&self, p: [f64; 3]) -> Uv {
        let d = std::array::from_fn(|i| p[i] - self.origin[i]);
        let [x, y, z] = self.frame.map(|axis| dot(d, axis));
        if matches!(self.surface, Surface::Plane(_)) {
            return [x, y];
        }
        let angle = |y, x| {
            let t = libm::atan2(y, x) / std::f64::consts::TAU;
            // Leave enough mantissa bits for exact translations by -2..=2.
            libm::round(t * 281_474_976_710_656.0) / 281_474_976_710_656.0
        };
        let u = if let Some([major, minor]) = self.ellipse_radii {
            angle(y / minor, x / major)
        } else {
            angle(y, x)
        };
        let v = match &self.surface {
            Surface::Sphere(_) => angle(z, libm::sqrt(x * x + y * y)),
            Surface::Torus(s) => angle(z, libm::sqrt(x * x + y * y) - s.major.to_f64()),
            _ => z,
        };
        [u, v]
    }

    pub(crate) fn canonical(&self, mut uv: Uv) -> Uv {
        for (x, periodic) in uv.iter_mut().zip(self.periodic) {
            if periodic {
                *x -= libm::floor(*x);
            }
            if *x == 0.0 {
                *x = 0.0;
            }
        }
        // Every longitude denotes the same pole/apex vertex.
        if (matches!(self.surface, Surface::Sphere(_))
            && uv[1].abs().to_bits() == 0.25f64.to_bits())
            || (matches!(self.surface, Surface::Cone(_)) && uv[1] == 0.0)
        {
            uv[0] = 0.0;
        }
        uv
    }

    pub(crate) fn evaluate(&self, uv: Uv, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        let mut parameters = uv.map(|x| Ball::from_rational(&q(x), BITS));
        if let Some(turn) = &self.turn {
            parameters[0] = &parameters[0] * turn;
            if matches!(self.surface, Surface::Sphere(_) | Surface::Torus(_)) {
                parameters[1] = &parameters[1] * turn;
            }
        }
        self.surface.evaluate(&parameters[0], &parameters[1], dec)
    }

    pub(crate) fn sample(
        &self,
        uv: Uv,
        boundary: Option<[f64; 3]>,
        tolerance: &Rational,
        dec: &Decider,
    ) -> Result<([f64; 3], Rational), Refusal> {
        dec.check_deadline()?;
        if let Some(values) = self
            .evaluator
            .as_ref()
            .and_then(|e| e.evaluate(&self.surface, uv))
        {
            let position = boundary.unwrap_or_else(|| values.map(|v| v.mid()));
            let deltas = std::array::from_fn::<_, 3, _>(|i| {
                catallax_numbers::Interval::point(position[i]).unwrap() - values[i]
            });
            let error = deltas
                .iter()
                .all(catallax_numbers::Interval::is_finite)
                .then(|| {
                    deltas.iter().fold(Rational::zero(), |sum, delta| {
                        sum + q(delta.lo().abs().max(delta.hi().abs()))
                    })
                });
            if let Some(error) = error
                && error <= tolerance.mul_pow2(-4)
            {
                return Ok((position, error));
            }
        }
        let balls = self.evaluate(uv, dec)?;
        let position = match boundary {
            Some(position) => position,
            None => position(&balls)?,
        };
        let error = residual(&position, &balls);
        Ok((position, error))
    }

    /// Euclidean upper bounds on S_uu, S_uv and S_vv in this chart.
    pub(crate) fn hessian(&self, vmax: f64, dec: &Decider) -> Result<[Rational; 3], Refusal> {
        if matches!(self.surface, Surface::Plane(_)) {
            return Ok(std::array::from_fn(|_| Rational::zero()));
        }
        let tau = self.turn.as_ref().unwrap().hi();
        let tau2 = &tau * &tau;
        let zero = Rational::zero();
        let upper = |value: &Real| -> Result<Rational, Refusal> {
            if let Some(value) = value.rational_ref() {
                return Ok(value.clone());
            }
            Ok(dec.enclose_at(value, BITS)?.hi())
        };
        Ok(match &self.surface {
            Surface::Plane(_) => [zero.clone(), zero.clone(), zero],
            Surface::Cylinder(s) => [&upper(&s.radius())? * &tau2, zero.clone(), zero],
            Surface::Cone(s) => {
                let slope = upper(&s.slope)?;
                [&(&slope * &q(vmax.abs())) * &tau2, &slope * &tau, zero]
            }
            Surface::Sphere(s) => {
                let radius = upper(&s.radius)?;
                std::array::from_fn(|_| &radius * &tau2)
            }
            Surface::Torus(s) => {
                let minor = &upper(&s.minor)? * &tau2;
                [
                    &upper(&(&s.major + &s.minor))? * &tau2,
                    minor.clone(),
                    minor,
                ]
            }
            Surface::EllipticCylinder(s) => {
                [&upper(&s.major_sq.sqrt())? * &tau2, zero.clone(), zero]
            }
            _ => unreachable!(),
        })
    }
}

pub(crate) fn align(mut b: Uv, a: Uv, periodic: [bool; 2]) -> Uv {
    for i in 0..2 {
        if periodic[i] {
            b[i] += libm::round(a[i] - b[i]);
        }
    }
    b
}

pub(crate) fn filtered_chord_bound(uvs: &[Uv], hessian: &[Interval; 3]) -> Interval {
    let [u, v] = std::array::from_fn::<_, 2, _>(|i| {
        let lo = uvs.iter().map(|uv| uv[i]).fold(f64::INFINITY, f64::min);
        let hi = uvs.iter().map(|uv| uv[i]).fold(f64::NEG_INFINITY, f64::max);
        Interval::point(hi).unwrap() - Interval::point(lo).unwrap()
    });
    (hessian[0] * u.square()
        + hessian[1] * u * v * Interval::point(2.0).unwrap()
        + hessian[2] * v.square())
        * Interval::point(0.125).unwrap()
}

/// Taylor's theorem bounds interpolation error by half a Hessian-weighted
/// variance. Each parameter variance is at most its squared range / 4;
/// Cauchy-Schwarz bounds the mixed term. Endpoint residuals interpolate convexly.
pub(crate) fn chord_bound(uvs: &[Uv], hessian: &[Rational; 3]) -> Rational {
    let span: [Rational; 2] = std::array::from_fn(|i| {
        let lo = uvs.iter().map(|uv| uv[i]).fold(f64::INFINITY, f64::min);
        let hi = uvs.iter().map(|uv| uv[i]).fold(f64::NEG_INFINITY, f64::max);
        q(hi) - q(lo)
    });
    let [u, v] = &span;
    (&(&hessian[0] * u) * u
        + &(&hessian[1] * u) * v * Rational::from_i64(2)
        + &(&hessian[2] * v) * v)
        .mul_pow2(-3)
}

pub(crate) fn key(uv: Uv) -> [u64; 2] {
    uv.map(|value| if value == 0.0 { 0 } else { value.to_bits() })
}
