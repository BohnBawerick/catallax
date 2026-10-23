//! Rational charts, exact witnesses and uniform interpolation bounds.

use catallax_numbers::{Ball, Rational, Real, Refusal};
use catallax_shapes::{RPoint3, RVector3, Surface, Vector3};

use super::{Cell, failure};

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn rational(x: &Real) -> Result<Rational, Refusal> {
    x.as_rational().ok_or_else(|| {
        failure(
            "field.display-unsupported",
            "display chart needs rational frame and dimensions",
        )
    })
}
fn unit(v: &Vector3) -> Result<[Rational; 3], Refusal> {
    let v = [rational(&v.x)?, rational(&v.y)?, rational(&v.z)?];
    let squared = v.iter().fold(q(0), |s, x| s + x * x);
    let length = crate::profile::sqrt(&squared)
        .filter(|n| !n.is_zero())
        .ok_or_else(|| {
            failure(
                "field.display-unsupported",
                "display chart needs a rational unit frame",
            )
        })?;
    Ok(v.map(|x| &x / &length))
}
#[derive(Clone, Copy)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
}

pub(super) struct Chart {
    kind: Kind,
    origin: [Rational; 3],
    frame: [[Rational; 3]; 3],
    radius: Rational,
    minor: Rational,
    quadrant: [usize; 2],
    pub(super) domain: Cell,
    pub(super) lipschitz: [Rational; 2],
    pub(super) hessian: [Rational; 3],
}
impl Chart {
    pub(super) fn all(surface: &Surface, bound: &Rational) -> Result<Vec<Self>, Refusal> {
        let cylinder_radius = match surface {
            Surface::Cylinder(s) => Some(s.radius()),
            _ => None,
        };
        let radial_frame = if let Surface::Cylinder(s) = surface
            && cylinder_radius.as_ref().unwrap().rational_ref().is_none()
        {
            let radial = s.radial_reference().as_rational().ok_or_else(|| {
                failure(
                    "field.display-unsupported",
                    "cylinder has no rational point-on-circle frame",
                )
            })?;
            let axis = unit(&s.axis)?;
            let tangent =
                RVector3::new(axis[0].clone(), axis[1].clone(), axis[2].clone()).cross(&radial);
            Some([
                [radial.x, radial.y, radial.z],
                [tangent.x, tangent.y, tangent.z],
                axis,
            ])
        } else {
            None
        };
        let radial_bound = radial_frame
            .as_ref()
            .map(|f| f[0].iter().fold(q(0), |sum, x| sum + x.abs()));
        let (kind, origin, axis, reference, radius, minor) = match surface {
            Surface::Plane(s) => (Kind::Plane, &s.origin, &s.normal, &s.ref_dir, q(0), q(0)),
            Surface::Cylinder(s) => (
                Kind::Cylinder,
                &s.origin,
                &s.axis,
                &s.ref_dir,
                if radial_frame.is_some() {
                    q(1)
                } else {
                    rational(cylinder_radius.as_ref().unwrap())?
                },
                q(0),
            ),
            Surface::Cone(s) => (
                Kind::Cone,
                &s.apex,
                &s.axis,
                &s.ref_dir,
                rational(&s.slope)?,
                q(0),
            ),
            Surface::Sphere(s) => (
                Kind::Sphere,
                &s.center,
                &s.axis,
                &s.ref_dir,
                rational(&s.radius)?,
                q(0),
            ),
            Surface::Torus(s) => (
                Kind::Torus,
                &s.center,
                &s.axis,
                &s.ref_dir,
                rational(&s.major)?,
                rational(&s.minor)?,
            ),
            _ => {
                return Err(failure(
                    "field.display-unsupported",
                    "no certified native display chart for this supporting surface",
                ));
            }
        };
        let origin = [
            rational(&origin.x)?,
            rational(&origin.y)?,
            rational(&origin.z)?,
        ];
        let frame = if let Some(frame) = radial_frame {
            frame
        } else {
            [unit(reference)?, unit(&axis.cross(reference))?, unit(axis)?]
        };
        let extent = bound + &origin.iter().fold(q(0), |sum, x| sum + x.abs());
        let mut b = q(1);
        while b < extent {
            b = b.mul_pow2(1);
        }
        let [lu, lv, huu, huv, hvv] = match kind {
            Kind::Plane => [q(1), q(1), q(0), q(0), q(0)],
            Kind::Cylinder => {
                let bound = radial_bound.as_ref().unwrap_or(&radius);
                [bound.mul_pow2(1), q(1), bound.mul_pow2(3), q(0), q(0)]
            }
            Kind::Cone => [
                (&radius * &b).mul_pow2(1),
                &q(1) + &radius,
                (&radius * &b).mul_pow2(3),
                radius.mul_pow2(1),
                q(0),
            ],
            Kind::Sphere => [
                radius.mul_pow2(1),
                radius.mul_pow2(1),
                radius.mul_pow2(3),
                radius.mul_pow2(2),
                radius.mul_pow2(3),
            ],
            Kind::Torus => [
                (&radius + &minor).mul_pow2(1),
                minor.mul_pow2(1),
                (&radius + &minor).mul_pow2(3),
                minor.mul_pow2(2),
                minor.mul_pow2(3),
            ],
        };
        let u = if matches!(kind, Kind::Plane) {
            [-&b, b.clone()]
        } else {
            [q(0), q(1)]
        };
        let v = match kind {
            Kind::Sphere => [q(-1), q(1)],
            Kind::Torus => [q(0), q(1)],
            _ => [-&b, b],
        };
        let mut charts = Vec::new();
        for i in 0..if matches!(kind, Kind::Plane) { 1 } else { 4 } {
            for j in 0..if matches!(kind, Kind::Torus) { 4 } else { 1 } {
                charts.push(Self {
                    kind,
                    origin: origin.clone(),
                    frame: frame.clone(),
                    radius: radius.clone(),
                    minor: minor.clone(),
                    quadrant: [i, j],
                    domain: Cell {
                        range: [u.clone(), v.clone()],
                        depth: 0,
                    },
                    lipschitz: [lu.clone(), lv.clone()],
                    hessian: [huu.clone(), huv.clone(), hvv.clone()],
                });
            }
        }
        Ok(charts)
    }
    pub(super) fn point(&self, uv: &[Rational; 2]) -> RPoint3 {
        let circle = |t: &Rational, quadrant: usize| {
            let t2 = t * t;
            let denominator = &q(1) + &t2;
            let c = (&q(1) - &t2) / denominator.clone();
            let s = t.mul_pow2(1) / denominator;
            match quadrant {
                0 => [c, s],
                1 => [-s, c],
                2 => [-c, -s],
                _ => [s, -c],
            }
        };
        let [u, v] = uv;
        let local = if matches!(self.kind, Kind::Plane) {
            [u.clone(), v.clone(), q(0)]
        } else {
            let [c, s] = circle(u, self.quadrant[0]);
            match self.kind {
                Kind::Cylinder => [&c * &self.radius, &s * &self.radius, v.clone()],
                Kind::Cone => [
                    &(&c * &self.radius) * v,
                    &(&s * &self.radius) * v,
                    v.clone(),
                ],
                Kind::Sphere => {
                    let [cv, sv] = circle(v, 0);
                    [
                        &(&c * &cv) * &self.radius,
                        &(&s * &cv) * &self.radius,
                        &sv * &self.radius,
                    ]
                }
                Kind::Torus => {
                    let [cv, sv] = circle(v, self.quadrant[1]);
                    let r = &self.radius + &(&cv * &self.minor);
                    [&c * &r, &s * &r, &sv * &self.minor]
                }
                Kind::Plane => unreachable!(),
            }
        };
        let p: [Rational; 3] = std::array::from_fn(|i| {
            (0..3).fold(self.origin[i].clone(), |sum, j| {
                sum + &self.frame[j][i] * &local[j]
            })
        });
        let [x, y, z] = p;
        RPoint3::new(x, y, z)
    }
    /// A linear constraint is affine in height. Exact Bernstein signs of its
    /// angular derivative prove when corner values bound the entire cell.
    pub(super) fn corner_signs_suffice(&self, normal: &[Rational; 3], cell: &Cell) -> bool {
        match self.kind {
            Kind::Plane => true,
            Kind::Cylinder => {
                let dot =
                    |axis: usize| (0..3).fold(q(0), |s, i| s + &self.frame[axis][i] * &normal[i]);
                let (a, b) = (dot(0), dot(1));
                if a.is_zero() || b.is_zero() {
                    return true;
                }
                let (a, b) = match self.quadrant[0] {
                    0 => (a, b),
                    1 => (b, -a),
                    2 => (-a, -b),
                    _ => (-b, a),
                };
                let [lo, hi] = &cell.range[0];
                // d/dt of a*cos(t)+b*sin(t) has the sign of
                // b*(1-t²)-2*a*t in the rational quarter-circle chart.
                let coefficients = [
                    &b * &(&q(1) - &(lo * lo)) - (&a * lo).mul_pow2(1),
                    &b * &(&q(1) - &(lo * hi)) - &a * &(lo + hi),
                    &b * &(&q(1) - &(hi * hi)) - (&a * hi).mul_pow2(1),
                ];
                coefficients.iter().all(|x| x >= &q(0)) || coefficients.iter().all(|x| x <= &q(0))
            }
            _ => false,
        }
    }
    pub(super) fn cover(&self, cell: &Cell, p: &RPoint3) -> [Ball; 3] {
        // Coordinate-wise derivative bounds avoid widening directions that the
        // chart does not move, especially plane/plane and cylinder/cap trims.
        let span = cell.span();
        let mut radii = [q(0), q(0), q(0)];
        for (i, radius) in radii.iter_mut().enumerate() {
            let (lu, lv) = match self.kind {
                Kind::Plane => (self.frame[0][i].abs(), self.frame[1][i].abs()),
                Kind::Cylinder => (
                    (&self.radius * &(self.frame[0][i].abs() + self.frame[1][i].abs())).mul_pow2(1),
                    self.frame[2][i].abs(),
                ),
                _ => (self.lipschitz[0].clone(), self.lipschitz[1].clone()),
            };
            *radius = (&lu * &span[0] + &lv * &span[1]).mul_pow2(-1);
        }
        let coordinates = [&p.x, &p.y, &p.z];
        std::array::from_fn(|i| {
            Ball::from_rational(&(coordinates[i] - &radii[i]), 128)
                .hull(&Ball::from_rational(&(coordinates[i] + &radii[i]), 128))
        })
    }
    pub(super) fn chord(&self, cell: &Cell) -> Rational {
        let [u, v] = cell.span();
        (&(&self.hessian[0] * &u) * &u
            + (&(&self.hessian[1] * &u) * &v).mul_pow2(1)
            + &(&self.hessian[2] * &v) * &v)
            .mul_pow2(-3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Cylinder, Point3};

    #[test]
    fn irrational_cylinder_chart_does_not_assume_monotonicity_across_an_extremum() {
        let surface = Surface::Cylinder(Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.into(), 0.into(), 1.into()),
            radius_sq: 2.into(),
            ref_dir: Vector3::new(1.into(), 1.into(), 0.into()),
        });
        let charts = Chart::all(&surface, &q(4)).unwrap();
        let normal = [q(1), q(0), q(0)];
        assert!(charts[0].corner_signs_suffice(&normal, &charts[0].domain));
        assert!(!charts[1].corner_signs_suffice(&normal, &charts[1].domain));
        for range in [
            [q(0), Rational::from_ratio(1, 4).unwrap()],
            [Rational::from_ratio(3, 4).unwrap(), q(1)],
        ] {
            let cell = Cell {
                range: [range, [q(-1), q(1)]],
                depth: 2,
            };
            assert!(charts[1].corner_signs_suffice(&normal, &cell));
        }
        let mut dec = catallax_numbers::Decider::default();
        for chart in charts {
            for uv in [
                [q(0), q(0)],
                [Rational::from_ratio(1, 3).unwrap(), q(2)],
                [q(1), q(-2)],
            ] {
                let p = chart.point(&uv);
                assert!(
                    dec.sign(&surface.implicit().eval_real(&[
                        Real::from(p.x),
                        Real::from(p.y),
                        Real::from(p.z)
                    ],))
                        .is_coincident()
                );
            }
        }
    }
}
