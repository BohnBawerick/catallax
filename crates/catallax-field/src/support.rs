//! Signed 1-Lipschitz residuals of untrimmed analytic supports.

use crate::metric::{PRECISION, Scalar};
use catallax_numbers::{Ball, Rational};
use catallax_shapes::{RPoint3, Surface};

#[derive(Clone, Copy, Debug)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
}
#[derive(Clone, Debug)]
pub(crate) struct Support {
    kind: Kind,
    balls: [Ball; 10],
    fast: [f64; 10],
}
impl Support {
    pub fn new(surface: &Surface) -> Self {
        let (kind, origin, axis, a, b) = match surface {
            Surface::Plane(s) => (
                Kind::Plane,
                &s.origin,
                &s.normal,
                Rational::zero(),
                Rational::zero(),
            ),
            Surface::Cylinder(s) => (
                Kind::Cylinder,
                &s.origin,
                &s.axis,
                s.radius.clone(),
                Rational::zero(),
            ),
            Surface::Cone(s) => (
                Kind::Cone,
                &s.apex,
                &s.axis,
                s.slope.clone(),
                Rational::zero(),
            ),
            Surface::Sphere(s) => (
                Kind::Sphere,
                &s.center,
                &s.axis,
                s.radius.clone(),
                Rational::zero(),
            ),
            Surface::Torus(s) => (
                Kind::Torus,
                &s.center,
                &s.axis,
                s.major.clone(),
                s.minor.clone(),
            ),
            _ => unreachable!("feature analytic support"),
        };
        let q = [
            origin.x.clone(),
            origin.y.clone(),
            origin.z.clone(),
            axis.x.clone(),
            axis.y.clone(),
            axis.z.clone(),
            axis.norm_squared(),
            a.clone(),
            b,
            &a * &a + Rational::one(),
        ];
        Self {
            kind,
            balls: q.each_ref().map(|q| Ball::from_rational(q, PRECISION)),
            fast: q.each_ref().map(Rational::to_f64_nearest),
        }
    }
    pub fn proven(&self, p: &RPoint3) -> Ball {
        self.evaluate(
            &self.balls,
            &[&p.x, &p.y, &p.z].map(|q| Ball::from_rational(q, PRECISION)),
        )
    }
    pub fn fast(&self, p: &[f64; 3]) -> f64 {
        self.evaluate(&self.fast, p)
    }
    fn evaluate<T: Scalar>(&self, q: &[T; 10], p: &[T; 3]) -> T {
        let d: [T; 3] = std::array::from_fn(|i| p[i].sub(&q[i]));
        let norm = |p: &[T; 3]| p[0].square().add(&p[1].square()).add(&p[2].square()).sqrt();
        let abs = |x: T| x.max(&x.neg());
        if matches!(self.kind, Kind::Sphere) {
            return abs(norm(&d).sub(&q[7]));
        }
        let length = q[6].sqrt();
        let z = d[0]
            .mul(&q[3])
            .add(&d[1].mul(&q[4]))
            .add(&d[2].mul(&q[5]))
            .div(&length);
        if matches!(self.kind, Kind::Plane) {
            return abs(z);
        }
        let cross = [
            d[1].mul(&q[5]).sub(&d[2].mul(&q[4])),
            d[2].mul(&q[3]).sub(&d[0].mul(&q[5])),
            d[0].mul(&q[4]).sub(&d[1].mul(&q[3])),
        ];
        let rho = norm(&cross).div(&length);
        match self.kind {
            Kind::Cylinder => abs(rho.sub(&q[7])),
            Kind::Cone => abs(rho.sub(&q[7].mul(&abs(z)))).div(&q[9].sqrt()),
            Kind::Torus => abs(rho.sub(&q[7]).square().add(&z.square()).sqrt().sub(&q[8])),
            _ => unreachable!(),
        }
    }
}
