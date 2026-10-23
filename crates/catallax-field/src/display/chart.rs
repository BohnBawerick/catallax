//! Rational charts, exact witnesses and uniform interpolation bounds.

use catallax_numbers::{Interval, Rational, Real, Refusal};
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
    Disk,
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
    cached_origin: [Interval; 3],
    cached_frame: [[Interval; 3]; 3],
    cached_radius: Interval,
    cached_minor: Interval,
    pub(super) domain: Cell,
    pub(super) lipschitz: [Rational; 2],
    pub(super) hessian: [Rational; 3],
}
impl Chart {
    pub(super) fn is_disk(&self) -> bool {
        matches!(self.kind, Kind::Disk)
    }

    /// A nonzero affine pullback can vanish only on a set with no open planar
    /// patch. No point sample or coefficient enclosure is an identity proof.
    pub(super) fn has_nonzero_planar_restriction(
        &self,
        polynomial: &catallax_numbers::MPoly,
    ) -> bool {
        use catallax_numbers::MPoly;
        if !matches!(self.kind, Kind::Plane | Kind::Disk) {
            return false;
        }
        let coordinates: [MPoly; 3] = std::array::from_fn(|i| {
            MPoly::from_terms(
                2,
                vec![
                    (vec![0, 0], self.origin[i].clone()),
                    (vec![1, 0], self.frame[0][i].clone()),
                    (vec![0, 1], self.frame[1][i].clone()),
                ],
            )
        });
        let mut restricted = MPoly::zero(2);
        for (monomial, coefficient) in polynomial.terms() {
            let mut term = MPoly::constant(2, coefficient);
            for (axis, &power) in monomial.exponents().iter().enumerate() {
                for _ in 0..power {
                    term = &term * &coordinates[axis];
                }
            }
            restricted = &restricted + &term;
        }
        !restricted.is_zero()
    }

    pub(super) fn disk(
        plane: &catallax_shapes::Plane,
        cylinder: &catallax_shapes::Cylinder,
        boxes: &[&[Interval; 3]],
        bound: &Rational,
    ) -> Result<Vec<Self>, Refusal> {
        let axis = unit(&cylinder.axis)?;
        let normal = unit(&plane.normal)?;
        let origin = [
            rational(&cylinder.origin.x)?,
            rational(&cylinder.origin.y)?,
            rational(&cylinder.origin.z)?,
        ];
        let plane_origin = [
            rational(&plane.origin.x)?,
            rational(&plane.origin.y)?,
            rational(&plane.origin.z)?,
        ];
        let offset = (0..3).fold(q(0), |s, i| {
            s + &normal[i] * &(&plane_origin[i] - &origin[i])
        });
        let divisor = (0..3).fold(q(0), |s, i| s + &normal[i] * &axis[i]);
        let centre: [Rational; 3] =
            std::array::from_fn(|i| &origin[i] + &(&axis[i] * &(&offset / &divisor)));
        let radius = boxes
            .iter()
            .map(|b| {
                (0..3).fold(q(0), |sum, i| {
                    let lo = (&Rational::from_f64_exact(b[i].lo()).unwrap() - &centre[i]).abs();
                    let hi = (&Rational::from_f64_exact(b[i].hi()).unwrap() - &centre[i]).abs();
                    sum + lo.max(hi)
                })
            })
            .max()
            .expect("source boxes");
        let support = Surface::Cylinder(catallax_shapes::Cylinder {
            origin: RPoint3::new(centre[0].clone(), centre[1].clone(), centre[2].clone()).into(),
            ..cylinder.clone()
        });
        let mut charts = Self::all(&support, bound)?;
        // Quarter charts share the radial frame. A squared-radius cylinder
        // may use a rational, nonunit frame: bound its parameter extent and
        // physical distances once for the entire cap.
        let radial = &charts[0].frame[0];
        let squared = radial.iter().fold(q(0), |s, x| s + x * x);
        let (lower, upper) = if let Some(length) = crate::profile::sqrt(&squared) {
            (length.clone(), length)
        } else {
            let upper = radial.iter().fold(q(0), |s, x| s + x.abs());
            (&squared / &upper, upper)
        };
        let extent = &radius / &lower;
        let physical = &extent * &upper;
        let domain = [q(0), extent];
        let lipschitz = [physical.mul_pow2(1), upper.clone()];
        let hessian = [physical.mul_pow2(3), upper.mul_pow2(1), q(0)];
        for chart in &mut charts {
            chart.kind = Kind::Disk;
            chart.domain.range[1] = domain.clone();
            chart.lipschitz = lipschitz.clone();
            chart.hessian = hessian.clone();
        }
        Ok(charts)
    }

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
            Kind::Plane | Kind::Disk => [q(1), q(1), q(0), q(0), q(0)],
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
                    cached_origin: origin.each_ref().map(Interval::from_rational),
                    cached_frame: frame
                        .each_ref()
                        .map(|v| v.each_ref().map(Interval::from_rational)),
                    cached_radius: Interval::from_rational(&radius),
                    cached_minor: Interval::from_rational(&minor),
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
    /// Project every source box into linear chart coordinates. The union of
    /// those projections covers every source face on this support.
    pub(super) fn restrict(&mut self, boxes: &[&[Interval; 3]]) {
        if boxes.is_empty() || boxes.iter().any(|b| b.iter().any(|v| !v.is_finite())) {
            return;
        }
        for axis in 0..2 {
            let frame = match self.kind {
                Kind::Plane => &self.frame[axis],
                Kind::Cylinder if axis == 1 => &self.frame[2],
                _ => continue,
            };
            let mut bounds: Option<[Rational; 2]> = None;
            for b in boxes {
                let mut range = [q(0), q(0)];
                for i in 0..3 {
                    let a = &frame[i]
                        * &(&Rational::from_f64_exact(b[i].lo()).unwrap() - &self.origin[i]);
                    let z = &frame[i]
                        * &(&Rational::from_f64_exact(b[i].hi()).unwrap() - &self.origin[i]);
                    range[0] = &range[0] + (&a).min(&z);
                    range[1] = &range[1] + (&a).max(&z);
                }
                bounds = Some(match bounds {
                    None => range,
                    Some([lo, hi]) => [lo.min(range[0].clone()), hi.max(range[1].clone())],
                });
            }
            let [lo, hi] = bounds.expect("nonempty source boxes");
            self.domain.range[axis][0] = self.domain.range[axis][0].clone().max(lo);
            self.domain.range[axis][1] = self.domain.range[axis][1].clone().min(hi);
        }
    }

    /// Exact parameter cuts where an axis-aligned planar trim changes sign.
    pub(super) fn cuts(&self, constraints: &[Vec<super::Constraint>]) -> [Vec<Rational>; 2] {
        let mut cuts = [Vec::new(), Vec::new()];
        for c in constraints.iter().flatten() {
            let Some((normal, constant)) = &c.linear else {
                continue;
            };
            let dot = |v: &[Rational; 3]| (0..3).fold(q(0), |s, i| s + &v[i] * &normal[i]);
            let offset = &dot(&self.origin) + constant;
            let coefficients = match self.kind {
                Kind::Plane => [dot(&self.frame[0]), dot(&self.frame[1])],
                Kind::Cylinder
                    if dot(&self.frame[0]).is_zero() && dot(&self.frame[1]).is_zero() =>
                {
                    [q(0), dot(&self.frame[2])]
                }
                _ => continue,
            };
            for axis in 0..2 {
                if !coefficients[axis].is_zero() && coefficients[1 - axis].is_zero() {
                    let root = -&offset / coefficients[axis].clone();
                    if root > self.domain.range[axis][0] && root < self.domain.range[axis][1] {
                        cuts[axis].push(root);
                    }
                }
            }
        }
        for axis in &mut cuts {
            axis.sort();
            axis.dedup();
        }
        cuts
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
                Kind::Disk => [&c * v, &s * v, q(0)],
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
            Kind::Cylinder | Kind::Disk => {
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
    pub(super) fn cover(&self, cell: &Cell) -> [Interval; 3] {
        let range = |axis: usize| {
            Interval::new(
                Interval::from_rational(&cell.range[axis][0]).lo(),
                Interval::from_rational(&cell.range[axis][1]).hi(),
            )
            .unwrap()
        };
        let zero = Interval::point(0.0).unwrap();
        let one = Interval::point(1.0).unwrap();
        let circle = |axis: usize, quadrant: usize| {
            let at = |t: &Rational| {
                let t = Interval::from_rational(t);
                let square = t.square();
                [(one - square) / (one + square), (t + t) / (one + square)]
            };
            let a = at(&cell.range[axis][0]);
            let b = at(&cell.range[axis][1]);
            // Cosine is even and decreases with |t|. Sine increases on [-1, 1].
            let c = Interval::new(
                a[0].lo().min(b[0].lo()),
                if cell.range[axis][0] <= q(0) && cell.range[axis][1] >= q(0) {
                    1.0
                } else {
                    a[0].hi().max(b[0].hi())
                },
            )
            .unwrap();
            let s = Interval::new(a[1].lo(), b[1].hi()).unwrap();
            match quadrant {
                0 => [c, s],
                1 => [-s, c],
                2 => [-c, -s],
                _ => [s, -c],
            }
        };
        let local = if matches!(self.kind, Kind::Plane) {
            [range(0), range(1), zero]
        } else {
            let [c, s] = circle(0, self.quadrant[0]);
            match self.kind {
                Kind::Disk => [c * range(1), s * range(1), zero],
                Kind::Cylinder => [c * self.cached_radius, s * self.cached_radius, range(1)],
                Kind::Cone => [
                    c * self.cached_radius * range(1),
                    s * self.cached_radius * range(1),
                    range(1),
                ],
                Kind::Sphere => {
                    let [cv, sv] = circle(1, 0);
                    [
                        c * cv * self.cached_radius,
                        s * cv * self.cached_radius,
                        sv * self.cached_radius,
                    ]
                }
                Kind::Torus => {
                    let [cv, sv] = circle(1, self.quadrant[1]);
                    let radius = self.cached_radius + cv * self.cached_minor;
                    [c * radius, s * radius, sv * self.cached_minor]
                }
                Kind::Plane => unreachable!(),
            }
        };
        std::array::from_fn(|i| {
            (0..3).fold(self.cached_origin[i], |sum, j| {
                sum + self.cached_frame[j][i] * local[j]
            })
        })
    }
    pub(super) fn chord(&self, cell: &Cell) -> Rational {
        let [u, v] = cell.span();
        if self.is_disk() {
            // Each radial line is exact. The two triangles form the annular
            // sector's chord trapezoid, whose inner and outer circular arcs
            // deviate by at most max_physical_radius * du^2. This is a set-distance
            // bound, not interpolation at matching (angle, radius) parameters.
            let radius = &cell.range[1][1] * &self.lipschitz[1];
            return &(&radius * &u) * &u;
        }
        (&(&self.hessian[0] * &u) * &u
            + (&(&self.hessian[1] * &u) * &v).mul_pow2(1)
            + &(&self.hessian[2] * &v) * &v)
            .mul_pow2(-3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Cone, Cylinder, Plane, Point3, Sphere, Torus};

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
        let implicit = surface.implicit().expect("analytic fixture support");
        for chart in charts {
            for uv in [
                [q(0), q(0)],
                [Rational::from_ratio(1, 3).unwrap(), q(2)],
                [q(1), q(-2)],
            ] {
                let p = chart.point(&uv);
                assert!(
                    dec.sign(&implicit.eval_real(&[
                        Real::from(p.x),
                        Real::from(p.y),
                        Real::from(p.z)
                    ],))
                        .is_coincident()
                );
            }
        }
    }

    #[test]
    fn scaled_radial_cap_chord_bounds_physical_distance() {
        for scale in [q(10), Rational::from_ratio(1, 10).unwrap()] {
            let cylinder = Cylinder {
                origin: Point3::origin(),
                axis: Vector3::new(0.into(), 0.into(), 1.into()),
                radius_sq: (&scale * &scale).mul_pow2(1).into(),
                ref_dir: Vector3::new(scale.clone().into(), scale.clone().into(), 0.into()),
            };
            let plane = Plane {
                origin: Point3::origin(),
                normal: cylinder.axis.clone(),
                ref_dir: cylinder.ref_dir.clone(),
            };
            let extent = (&scale * &Rational::from_ratio(3, 2).unwrap()).to_f64_nearest();
            let bounds = [Interval::new(-extent, extent).unwrap(); 3];
            let charts = Chart::disk(&plane, &cylinder, &[&bounds], &q(32)).unwrap();
            for chart in charts {
                // Parameter radius one is on the source circle at both scales.
                assert!(chart.domain.range[1][1] >= q(1), "scale={scale}");
                let cell = Cell {
                    range: [[q(0), Rational::from_ratio(1, 4).unwrap()], [q(0), q(1)]],
                    depth: 2,
                };
                let a = chart.point(&[q(0), q(1)]);
                let b = chart.point(&[Rational::from_ratio(1, 4).unwrap(), q(1)]);
                let p = chart.point(&[Rational::from_ratio(1, 8).unwrap(), q(1)]);
                let chord = b.minus(&a);
                let offset = p.minus(&a);
                let projection = offset.dot(&chord);
                assert!(projection >= q(0) && projection <= chord.norm_squared());
                let distance_squared = offset.cross(&chord).norm_squared() / chord.norm_squared();
                let bound = chart.chord(&cell);
                assert!(distance_squared <= &bound * &bound, "scale={scale}");
            }
        }
    }

    #[test]
    fn interval_covers_contain_exact_chart_points_in_every_quadrant() {
        let origin = Point3::new(
            Rational::from_ratio(7_000_001, 7).unwrap().into(),
            Rational::from_ratio(-11, 13).unwrap().into(),
            q(17).into(),
        );
        let axis = Vector3::new(0.into(), 0.into(), 1.into());
        let reference = Vector3::new(
            Rational::from_ratio(3, 5).unwrap().into(),
            Rational::from_ratio(4, 5).unwrap().into(),
            0.into(),
        );
        let surfaces = [
            Surface::Plane(Plane {
                origin: origin.clone(),
                normal: axis.clone(),
                ref_dir: reference.clone(),
            }),
            Surface::Cylinder(Cylinder {
                origin: origin.clone(),
                axis: axis.clone(),
                ref_dir: reference.clone(),
                radius_sq: q(9).into(),
            }),
            Surface::Cylinder(Cylinder {
                origin: origin.clone(),
                axis: axis.clone(),
                ref_dir: Vector3::new(10.into(), 10.into(), 0.into()),
                radius_sq: q(200).into(),
            }),
            Surface::Cone(Cone {
                apex: origin.clone(),
                axis: axis.clone(),
                ref_dir: reference.clone(),
                slope: Rational::from_ratio(2, 3).unwrap().into(),
            }),
            Surface::Sphere(Sphere {
                center: origin.clone(),
                axis: axis.clone(),
                ref_dir: reference.clone(),
                radius: q(3).into(),
            }),
            Surface::Torus(Torus {
                center: origin,
                axis,
                ref_dir: reference,
                major: q(5).into(),
                minor: q(2).into(),
            }),
        ];
        for surface in surfaces {
            let mut charts = Chart::all(&surface, &q(32)).unwrap();
            if let Surface::Cylinder(cylinder) = &surface {
                let plane = Plane {
                    origin: cylinder.origin.clone(),
                    normal: cylinder.axis.clone(),
                    ref_dir: cylinder.ref_dir.clone(),
                };
                let bounds = [&plane.origin.x, &plane.origin.y, &plane.origin.z].map(|v| {
                    let v = rational(v).unwrap();
                    Interval::from_rational(&(&v - &q(32)))
                        .hull(&Interval::from_rational(&(&v + &q(32))))
                });
                charts.extend(Chart::disk(&plane, cylinder, &[&bounds], &q(32)).unwrap());
            }
            for chart in charts {
                let mut cells = vec![chart.domain.clone()];
                for axis in 0..2 {
                    cells = cells
                        .into_iter()
                        .flat_map(|cell| cell.split(axis, &[]))
                        .collect();
                }
                for cell in cells {
                    let cover = chart.cover(&cell);
                    for u in 0..=4 {
                        for v in 0..=4 {
                            let uv = [u, v].map(|i| q(i).mul_pow2(-2));
                            let p = chart.point(&std::array::from_fn(|i| {
                                &cell.range[i][0] + &(&cell.span()[i] * &uv[i])
                            }));
                            for (coordinate, interval) in [p.x, p.y, p.z].iter().zip(cover) {
                                assert!(
                                    Rational::from_f64_exact(interval.lo()).unwrap() <= *coordinate
                                        && Rational::from_f64_exact(interval.hi()).unwrap()
                                            >= *coordinate,
                                    "{surface:?}: {coordinate} outside {interval:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn planar_side_sharing_checks_a_complete_polynomial_identity() {
        use catallax_numbers::MPoly;
        let plane = Surface::Plane(Plane {
            origin: Point3::origin(),
            normal: Vector3::new(0.into(), 0.into(), 1.into()),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        });
        let chart = Chart::all(&plane, &q(4)).unwrap().remove(0);
        let x = MPoly::var(3, 0);
        let y = MPoly::var(3, 1);
        let z = MPoly::var(3, 2);
        assert!(!chart.has_nonzero_planar_restriction(&(&x * &z)));
        assert!(chart.has_nonzero_planar_restriction(&(&x + &y)));
    }
}
