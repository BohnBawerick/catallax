//! Signed 1-Lipschitz residuals of untrimmed analytic supports.

use crate::metric::{PRECISION, Scalar};
use catallax_numbers::{Ball, Decider, Real, Refusal};
use catallax_shapes::{Point3, Pose, Surface, Vector3};

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
    exact: [Real; 10],
    balls: Option<[Ball; 10]>,
    fast: [f64; 10],
}
impl Support {
    pub fn new(surface: &Surface) -> Self {
        let (kind, origin, axis, a, b) = match surface {
            Surface::Plane(s) => (
                Kind::Plane,
                &s.origin,
                &s.normal,
                Real::zero(),
                Real::zero(),
            ),
            Surface::Cylinder(s) => (
                Kind::Cylinder,
                &s.origin,
                &s.axis,
                s.radius.clone(),
                Real::zero(),
            ),
            Surface::Cone(s) => (Kind::Cone, &s.apex, &s.axis, s.slope.clone(), Real::zero()),
            Surface::Sphere(s) => (
                Kind::Sphere,
                &s.center,
                &s.axis,
                s.radius.clone(),
                Real::zero(),
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
        let factor = &a * &a + Real::one();
        let exact = [
            origin.x.clone(),
            origin.y.clone(),
            origin.z.clone(),
            axis.x.clone(),
            axis.y.clone(),
            axis.z.clone(),
            axis.norm_squared(),
            a,
            b,
            factor,
        ];
        Self {
            kind,
            balls: rational_balls(&exact),
            fast: exact.each_ref().map(Real::to_f64),
            exact,
        }
    }
    pub fn transform(&mut self, pose: &Pose) {
        let origin = pose.apply_point(&Point3::new(
            self.exact[0].clone(),
            self.exact[1].clone(),
            self.exact[2].clone(),
        ));
        let axis = pose.apply_vector(&Vector3::new(
            self.exact[3].clone(),
            self.exact[4].clone(),
            self.exact[5].clone(),
        ));
        self.exact[0] = origin.x;
        self.exact[1] = origin.y;
        self.exact[2] = origin.z;
        self.exact[3] = axis.x;
        self.exact[4] = axis.y;
        self.exact[5] = axis.z;
        // A certified rigid motion preserves the defining axis norm and dimensions.
        self.balls = rational_balls(&self.exact);
        self.fast = self.exact.each_ref().map(Real::to_f64);
    }
    pub fn proven(&self, p: &[Ball; 3], dec: &mut Decider) -> Result<Ball, Refusal> {
        if let Some(q) = &self.balls {
            return Ok(self.evaluate(q, p));
        }
        let q = [
            crate::geometry::enclose(&self.exact[0], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[1], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[2], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[3], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[4], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[5], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[6], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[7], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[8], PRECISION, dec)?,
            crate::geometry::enclose(&self.exact[9], PRECISION, dec)?,
        ];
        Ok(self.evaluate(&q, p))
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

fn rational_balls(exact: &[Real; 10]) -> Option<[Ball; 10]> {
    exact.iter().all(|q| q.rational_ref().is_some()).then(|| {
        exact
            .each_ref()
            .map(|q| Ball::from_rational(q.rational_ref().expect("rational support"), PRECISION))
    })
}
