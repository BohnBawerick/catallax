//! Dyadic turns make translated periodic chart copies bit-identical.

use catallax_numbers::{Ball, Rational, Refusal};
use catallax_shapes::{RVector3, Surface};

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

fn unit(v: &RVector3) -> [f64; 3] {
    let length = Ball::from_rational(&v.norm_squared(), BITS)
        .sqrt()
        .expect("nonzero definition axis");
    [&v.x, &v.y, &v.z]
        .map(|x| midpoint(&Ball::from_rational(x, BITS).checked_div(&length).unwrap()))
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
}

impl Chart {
    pub(crate) fn new(surface: &Surface) -> Result<Self, Refusal> {
        let (origin, axis, reference, periodic) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal, &s.ref_dir, [false, false]),
            Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir, [true, false]),
            Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir, [true, false]),
            Surface::Sphere(s) => (&s.center, &s.axis, &s.ref_dir, [true, false]),
            Surface::Torus(s) => (&s.center, &s.axis, &s.ref_dir, [true, true]),
            _ => return Err(failure("unknown analytic surface definition")),
        };
        let origin = [&origin.x, &origin.y, &origin.z].map(Rational::to_f64_nearest);
        if !origin.iter().all(|v| v.is_finite()) {
            return Err(failure("surface origin exceeds binary64"));
        }
        Ok(Self {
            surface: surface.clone(),
            periodic,
            origin,
            frame: [unit(reference), unit(&axis.cross(reference)), unit(axis)],
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
        let u = angle(y, x);
        let v = match &self.surface {
            Surface::Sphere(_) => angle(z, libm::sqrt(x * x + y * y)),
            Surface::Torus(s) => angle(z, libm::sqrt(x * x + y * y) - s.major.to_f64_nearest()),
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

    pub(crate) fn evaluate(&self, uv: Uv) -> [Ball; 3] {
        let turn = Ball::pi(BITS) * Ball::from_rational(&Rational::from_i64(2), BITS);
        let mut parameters = uv.map(|x| Ball::from_rational(&q(x), BITS));
        if !matches!(self.surface, Surface::Plane(_)) {
            parameters[0] = &parameters[0] * &turn;
        }
        if matches!(self.surface, Surface::Sphere(_) | Surface::Torus(_)) {
            parameters[1] = &parameters[1] * &turn;
        }
        self.surface.evaluate(&parameters[0], &parameters[1])
    }

    /// Euclidean upper bounds on S_uu, S_uv and S_vv in this chart.
    pub(crate) fn hessian(&self, vmax: f64) -> [Rational; 3] {
        let tau = Ball::pi(BITS).hi() * Rational::from_i64(2);
        let tau2 = &tau * &tau;
        let zero = Rational::zero();
        match &self.surface {
            Surface::Plane(_) => [zero.clone(), zero.clone(), zero],
            Surface::Cylinder(s) => [&s.radius * &tau2, zero.clone(), zero],
            Surface::Cone(s) => [&(&s.slope * &q(vmax.abs())) * &tau2, &s.slope * &tau, zero],
            Surface::Sphere(s) => std::array::from_fn(|_| &s.radius * &tau2),
            Surface::Torus(s) => [
                &(&s.major + &s.minor) * &tau2,
                &s.minor * &tau2,
                &s.minor * &tau2,
            ],
            _ => unreachable!(),
        }
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

/// Taylor's theorem bounds interpolation error by half a Hessian-weighted
/// variance. Each parameter variance is at most its squared range / 4;
/// Cauchy-Schwarz bounds the mixed term. Endpoint residuals interpolate convexly.
pub(crate) fn chord_bound(uvs: &[Uv], hessian: &[Rational; 3]) -> Rational {
    let span: [Rational; 2] = std::array::from_fn(|i| {
        let lo = uvs.iter().map(|uv| q(uv[i])).min().unwrap();
        let hi = uvs.iter().map(|uv| q(uv[i])).max().unwrap();
        hi - lo
    });
    let [u, v] = &span;
    (&(&hessian[0] * u) * u
        + &(&hessian[1] * u) * v * Rational::from_i64(2)
        + &(&hessian[2] * v) * v)
        .mul_pow2(-3)
}

pub(crate) fn key(uv: Uv) -> [u64; 2] {
    uv.map(f64::to_bits)
}
