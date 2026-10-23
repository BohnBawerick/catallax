//! Signed 1-Lipschitz residuals of untrimmed analytic supports.

use crate::metric::{PRECISION, Scalar};
use catallax_numbers::{Ball, Decider, Real, Refusal};
use catallax_shapes::{Isometry, Point3, Surface, Vector3};

#[derive(Clone, Copy, Debug)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    EllipticCylinder,
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
                s.radius_sq.clone(),
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
            Surface::EllipticCylinder(s) => return Self::elliptic(s),
            _ => unreachable!("feature analytic support"),
        };
        let factor = if matches!(kind, Kind::Cone) {
            &a * &a + Real::one()
        } else {
            Real::one()
        };
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
        Self::from_exact(kind, exact)
    }
    fn from_exact(kind: Kind, exact: [Real; 10]) -> Self {
        let mut fast = exact.each_ref().map(Real::to_f64);
        if matches!(kind, Kind::Cylinder) {
            fast[7] = exact[7].sqrt().to_f64();
        }
        Self {
            kind,
            balls: rational_balls(kind, &exact),
            fast,
            exact,
        }
    }
    fn elliptic(surface: &catallax_shapes::EllipticCylinder) -> Self {
        // Pack the two scaled transverse vectors into the existing ten-value
        // cache. b*(||[x/a,y/b]||-1) is 1-Lipschitz for a >= b, so its
        // absolute value is a lower bound, never an affine-stretched distance.
        let major = surface.major_dir.scaled(
            &(Real::one() / (&surface.major_dir.norm_squared() * &surface.major_sq).sqrt()),
        );
        let minor = surface.axis.cross(&surface.major_dir);
        let minor =
            minor.scaled(&(Real::one() / (&minor.norm_squared() * &surface.minor_sq).sqrt()));
        let exact = [
            surface.origin.x.clone(),
            surface.origin.y.clone(),
            surface.origin.z.clone(),
            major.x,
            major.y,
            major.z,
            minor.x,
            minor.y,
            minor.z,
            surface.minor_sq.sqrt(),
        ];
        Self {
            kind: Kind::EllipticCylinder,
            balls: rational_balls(Kind::EllipticCylinder, &exact),
            fast: exact.each_ref().map(Real::to_f64),
            exact,
        }
    }
    pub fn transform(&mut self, pose: &Isometry) {
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
        if matches!(self.kind, Kind::EllipticCylinder) {
            let minor = pose.apply_vector(&Vector3::new(
                self.exact[6].clone(),
                self.exact[7].clone(),
                self.exact[8].clone(),
            ));
            self.exact[6] = minor.x;
            self.exact[7] = minor.y;
            self.exact[8] = minor.z;
        }
        // A certified rigid motion preserves vector norms and dimension scalars.
        self.balls = rational_balls(self.kind, &self.exact);
        let radius = self.fast[7];
        self.fast = self.exact.each_ref().map(Real::to_f64);
        if matches!(self.kind, Kind::Cylinder) {
            self.fast[7] = radius;
        }
    }
    pub fn proven(&self, p: &[Ball; 3], dec: &mut Decider) -> Result<Ball, Refusal> {
        if let Some(q) = &self.balls {
            return Ok(self.evaluate(q, p));
        }
        let mut q = [
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
        if matches!(self.kind, Kind::Cylinder) {
            q[7] = q[7].sqrt().expect("validated positive squared radius");
        }
        Ok(self.evaluate(&q, p))
    }
    pub fn fast(&self, p: &[f64; 3]) -> f64 {
        self.evaluate(&self.fast, p)
    }
    fn evaluate<T: Scalar>(&self, q: &[T; 10], p: &[T; 3]) -> T {
        let d: [T; 3] = std::array::from_fn(|i| p[i].sub(&q[i]));
        let norm = |p: &[T; 3]| p[0].square().add(&p[1].square()).add(&p[2].square()).sqrt();
        let abs = |x: T| x.max(&x.neg());
        if matches!(self.kind, Kind::EllipticCylinder) {
            let dot = |offset| {
                d[0].mul(&q[offset])
                    .add(&d[1].mul(&q[offset + 1]))
                    .add(&d[2].mul(&q[offset + 2]))
            };
            let rho = dot(3).square().add(&dot(6).square()).sqrt();
            return abs(rho.sub(&T::integer(1))).mul(&q[9]);
        }
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

fn rational_balls(kind: Kind, exact: &[Real; 10]) -> Option<[Ball; 10]> {
    exact.iter().all(|q| q.rational_ref().is_some()).then(|| {
        let mut balls = exact
            .each_ref()
            .map(|q| Ball::from_rational(q.rational_ref().expect("rational support"), PRECISION));
        if matches!(kind, Kind::Cylinder) {
            balls[7] = if let Some(radius) = crate::profile::sqrt(exact[7].rational_ref().unwrap())
            {
                Ball::from_rational(&radius, PRECISION)
            } else {
                balls[7].sqrt().expect("validated positive squared radius")
            };
        }
        balls
    })
}
