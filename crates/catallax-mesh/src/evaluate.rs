//! Allocation-free display enclosures, with the exact chart as the fallback.

use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal};
use catallax_shapes::{Surface, Vector3};

use crate::chart::{BITS, Uv};

fn point(x: f64) -> Interval {
    Interval::point(x).expect("finite chart coefficient")
}

fn unit(v: &Vector3, dec: &Decider) -> Result<[Interval; 3], Refusal> {
    let length = v.norm_squared().sqrt();
    let values = [&v.x, &v.y, &v.z].map(|x| {
        dec.enclose_at(&(x / &length), BITS)
            .map(|b| b.to_interval())
    });
    let [x, y, z] = values;
    Ok([x?, y?, z?])
}

/// Fixed definition coefficients only. No rounded geometry becomes a definition.
pub(crate) struct Evaluator {
    origin: [Interval; 3],
    frame: [[Interval; 3]; 3],
    dimensions: [Interval; 2],
    tau: Interval,
}

impl Evaluator {
    pub(crate) fn new(surface: &Surface, dec: &Decider) -> Result<Self, Refusal> {
        let (origin, axis, reference) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal, &s.ref_dir),
            Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
            Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
            Surface::Sphere(s) => (&s.center, &s.axis, &s.ref_dir),
            Surface::Torus(s) => (&s.center, &s.axis, &s.ref_dir),
            Surface::EllipticCylinder(s) => (&s.origin, &s.axis, &s.major_dir),
            _ => unreachable!("unknown analytic surface"),
        };
        let enclose = |x: &Real| dec.enclose_at(x, BITS).map(|b| b.to_interval());
        let [x, y, z] = [&origin.x, &origin.y, &origin.z].map(enclose);
        let dimensions = match surface {
            Surface::Plane(_) => [point(0.0); 2],
            Surface::Cylinder(s) => [enclose(&s.radius)?, point(0.0)],
            Surface::Cone(s) => [enclose(&s.slope)?, point(0.0)],
            Surface::Sphere(s) => [enclose(&s.radius)?, point(0.0)],
            Surface::Torus(s) => [enclose(&s.major)?, enclose(&s.minor)?],
            Surface::EllipticCylinder(s) => {
                [enclose(&s.major_sq.sqrt())?, enclose(&s.minor_sq.sqrt())?]
            }
            _ => unreachable!("unknown analytic surface"),
        };
        Ok(Self {
            origin: [x?, y?, z?],
            frame: [
                unit(reference, dec)?,
                unit(&axis.cross(reference), dec)?,
                unit(axis, dec)?,
            ],
            dimensions,
            tau: (Ball::pi_checked(BITS, dec)? * Ball::from_rational(&Rational::from_i64(2), BITS))
                .to_interval(),
        })
    }

    pub(crate) fn evaluate(&self, surface: &Surface, uv: Uv) -> Option<[Interval; 3]> {
        let v = point(uv[1]);
        let [a, b] = self.dimensions;
        let coefficients = if matches!(surface, Surface::Plane(_)) {
            [point(uv[0]), v, point(0.0)]
        } else {
            let [sin, cos] = self.trig(uv[0])?;
            match surface {
                Surface::Cylinder(_) => [a * cos, a * sin, v],
                Surface::Cone(_) => [v * a * cos, v * a * sin, v],
                Surface::Sphere(_) => {
                    let [sv, cv] = self.trig(uv[1])?;
                    [a * cv * cos, a * cv * sin, a * sv]
                }
                Surface::Torus(_) => {
                    let [sv, cv] = self.trig(uv[1])?;
                    let radial = a + b * cv;
                    [radial * cos, radial * sin, b * sv]
                }
                Surface::EllipticCylinder(_) => [a * cos, b * sin, v],
                _ => unreachable!("unknown analytic surface"),
            }
        };
        let result = std::array::from_fn(|i| {
            self.origin[i]
                + (self.frame[0][i] * coefficients[0] + self.frame[1][i] * coefficients[1])
                + self.frame[2][i] * coefficients[2]
        });
        result.iter().all(Interval::is_finite).then_some(result)
    }

    fn trig(&self, turns: f64) -> Option<[Interval; 2]> {
        // Reduction uses an exact quarter turn, not an approximate multiple of pi.
        // Every arithmetic operation, including the reduction, rounds outwards.
        let quadrant = libm::round(turns * 4.0);
        let x = (point(turns) - point(quadrant * 0.25)) * self.tau;
        if x.lo() < -1.0 || x.hi() > 1.0 {
            return None;
        }
        let square = -x.square();
        let mut sine_term = x;
        let mut cosine_term = point(1.0);
        let mut sine = sine_term;
        let mut cosine = cosine_term;
        for k in 1..=10 {
            let n = f64::from(2 * k);
            sine_term = sine_term * square / point(n * (n + 1.0));
            cosine_term = cosine_term * square / point(n * (n - 1.0));
            sine = sine + sine_term;
            cosine = cosine + cosine_term;
        }
        // On |x| <= 1 the alternating Taylor tails are bounded by 1/23!
        // and 1/22!, both strictly less than 2^-69.
        let tail = Interval::new(-1.6940658945086007e-21, 1.6940658945086007e-21).unwrap();
        sine = sine + tail;
        cosine = cosine + tail;
        Some(match (quadrant as i64).rem_euclid(4) {
            0 => [sine, cosine],
            1 => [cosine, -sine],
            2 => [-sine, -cosine],
            _ => [-cosine, sine],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{EllipticCylinder, Pose, RPoint3, RVector3, Rotation};
    use catallax_topology::{self as topology, Model};

    #[test]
    fn interval_charts_enclose_high_precision_definitions_across_quadrants() {
        let mut dec = Decider::default();
        let rotation = Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap();
        let pose = Pose::new(
            rotation,
            RVector3::new(7.into(), (-3).into(), 2.into()).into(),
            &mut dec,
        )
        .unwrap();
        let mut model = Model::new();
        let bodies = [
            topology::cylinder(&mut model, 3.into(), 5.into(), &pose).unwrap(),
            topology::cone(&mut model, 3.into(), 1.into(), 5.into(), &pose).unwrap(),
            topology::sphere(&mut model, 3.into(), &pose).unwrap(),
            topology::torus(&mut model, 5.into(), 2.into(), &pose).unwrap(),
        ];
        let mut surfaces: Vec<_> = bodies
            .into_iter()
            .flat_map(|b| model.faces_of_body(b))
            .map(|f| model.geom().surface(model.face(f).surface).clone())
            .collect();
        surfaces.push(
            Surface::EllipticCylinder(
                EllipticCylinder::new(
                    RPoint3::origin().into(),
                    RVector3::new(0.into(), 0.into(), 1.into()).into(),
                    RVector3::new(1.into(), 0.into(), 0.into()).into(),
                    9.into(),
                    2.into(),
                    &mut dec,
                )
                .unwrap(),
            )
            .transformed(&pose),
        );
        let tau = Ball::pi(256) * Ball::from_rational(&Rational::from_i64(2), 256);
        for surface in surfaces {
            let evaluator = Evaluator::new(&surface, &dec).unwrap();
            for u in [
                -0.5, -0.375, -0.25, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.875, 1.0,
            ] {
                for v in [-0.25, 0.0, 0.125, 0.25] {
                    let filtered = evaluator.evaluate(&surface, [u, v]).unwrap();
                    let mut parameters =
                        [u, v].map(|x| Ball::from_rational(&crate::chart::q(x), 256));
                    if !matches!(surface, Surface::Plane(_)) {
                        parameters[0] = &parameters[0] * &tau;
                    }
                    if matches!(surface, Surface::Sphere(_) | Surface::Torus(_)) {
                        parameters[1] = &parameters[1] * &tau;
                    }
                    let exact = surface
                        .evaluate(&parameters[0], &parameters[1], &dec)
                        .unwrap();
                    for i in 0..3 {
                        assert!(crate::chart::q(filtered[i].lo()) <= exact[i].lo());
                        assert!(crate::chart::q(filtered[i].hi()) >= exact[i].hi());
                    }
                }
            }
        }
    }
}
