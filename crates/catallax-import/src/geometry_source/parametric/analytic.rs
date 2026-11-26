use super::{Expr, Frame, MappingError};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use std::f64::consts::PI;

#[derive(Debug)]
pub(super) enum Kind {
    Plane,
    Cylinder(Rational),
    Cone(Rational, Rational),
    Sphere(Rational),
    Torus(Rational, Rational),
}

/// The unwrapped source chart. Radii and frame origin are in millimetres.
#[derive(Debug)]
pub(super) struct AnalyticChart {
    pub(super) frame: Frame,
    pub(super) length: Rational,
    pub(super) angle: Expr,
    pub(super) kind: Kind,
}

impl AnalyticChart {
    pub(super) fn sample(
        &self,
        [u, v]: &[Ball; 2],
        bits: u32,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        dec.charge_work(1)?;
        let number = |value| Ball::from_rational(value, bits);
        let height = || v * &number(&self.length);
        let local = if matches!(self.kind, Kind::Plane) {
            [
                u * &number(&self.length),
                height(),
                number(&Rational::zero()),
            ]
        } else {
            let angle = self.angle.evaluate(&[], bits, dec)?;
            let azimuth = u * &angle;
            let (sin_u, cos_u) = (azimuth.sin(), azimuth.cos());
            let (radial, z) = match &self.kind {
                Kind::Plane => unreachable!(),
                Kind::Cylinder(radius) => (number(radius), height()),
                Kind::Cone(radius, semi_angle) => {
                    let slope = cone_slope(&number(semi_angle), &angle).ok_or_else(|| {
                        Refusal::new(
                            "import.geometry-parameter-domain",
                            "parameter interval contains a singularity; subdivide the source chart",
                        )
                    })?;
                    let z = height();
                    (number(radius) + &z * &slope, z)
                }
                Kind::Sphere(radius) => {
                    let latitude = v * &angle;
                    let radius = number(radius);
                    (&radius * &latitude.cos(), &radius * &latitude.sin())
                }
                Kind::Torus(major, minor) => {
                    let meridian = v * &angle;
                    let minor = number(minor);
                    (
                        number(major) + &minor * &meridian.cos(),
                        &minor * &meridian.sin(),
                    )
                }
            };
            [&radial * &cos_u, &radial * &sin_u, z]
        };
        let point = [
            &self.frame.origin.x,
            &self.frame.origin.y,
            &self.frame.origin.z,
        ];
        let [x, y, z] = self
            .frame
            .axes
            .each_ref()
            .map(|axis| [&axis.x, &axis.y, &axis.z]);
        let coordinate = |i: usize| -> Result<Ball, MappingError> {
            // Preserve the expression chart's accumulator order and bounds.
            let mut result = dec.enclose_at(point[i], bits)?;
            for (axis, value) in [x[i], y[i], z[i]].into_iter().zip(&local) {
                result = result + &dec.enclose_at(axis, bits)? * value;
            }
            Ok(result)
        };
        Ok([coordinate(0)?, coordinate(1)?, coordinate(2)?])
    }

    pub(super) fn projection(
        &self,
        bits: u32,
        dec: &mut Decider,
    ) -> Result<Option<AnalyticProjection>, MappingError> {
        let length = self.length.to_f64_nearest();
        let angle = self.angle.evaluate(&[], bits, dec)?;
        let radians = angle.mid().to_f64_nearest();
        if !length.is_finite() || length <= 0.0 || !radians.is_finite() || radians == 0.0 {
            return Ok(None);
        }
        let radius = Rational::to_f64_nearest;
        let kind = match &self.kind {
            Kind::Plane => ProjectionKind::Plane,
            Kind::Cylinder(r) => ProjectionKind::Cylinder(radius(r)),
            Kind::Cone(r, a) => {
                let Some(slope) = cone_slope(&Ball::from_rational(a, bits), &angle) else {
                    return Ok(None);
                };
                ProjectionKind::Cone(radius(r), slope.mid().to_f64_nearest())
            }
            Kind::Sphere(r) => ProjectionKind::Sphere(radius(r)),
            Kind::Torus(r, t) => ProjectionKind::Torus(radius(r), radius(t)),
        };
        let nearest = |x: &Real| -> Result<f64, MappingError> {
            Ok(dec.enclose_at(x, bits)?.mid().to_f64_nearest())
        };
        let centre = [
            nearest(&self.frame.origin.x)?,
            nearest(&self.frame.origin.y)?,
            nearest(&self.frame.origin.z)?,
        ];
        let axis = |v: &catallax_shapes::Vector3| -> Result<_, MappingError> {
            Ok([nearest(&v.x)?, nearest(&v.y)?, nearest(&v.z)?])
        };
        let [x, y, z] = &self.frame.axes;
        Ok(Some(AnalyticProjection {
            centre,
            axes: [axis(x)?, axis(y)?, axis(z)?],
            length,
            radians,
            kind,
        }))
    }
}

fn cone_slope(angle: &Ball, scale: &Ball) -> Option<Ball> {
    let angle = angle * scale;
    angle.sin().checked_div(&angle.cos())
}

#[derive(Debug)]
enum ProjectionKind {
    Plane,
    Cylinder(f64),
    Cone(f64, f64),
    Sphere(f64),
    Torus(f64, f64),
}

/// Rounded inputs propose a location only. The caller verifies its residual
/// with the enclosing sample before accepting it.
#[derive(Debug)]
pub(crate) struct AnalyticProjection {
    centre: [f64; 3],
    axes: [[f64; 3]; 3],
    length: f64,
    radians: f64,
    kind: ProjectionKind,
}

impl AnalyticProjection {
    pub(crate) fn parameters(&self, target: [f64; 3], period: f64) -> Option<[f64; 2]> {
        let delta = std::array::from_fn::<_, 3, _>(|i| target[i] - self.centre[i]);
        let [x, y, z] = std::array::from_fn(|i| {
            let [a, b, c] = self.axes[i];
            a * delta[0] + b * delta[1] + c * delta[2]
        });
        if matches!(self.kind, ProjectionKind::Plane) {
            return Some([x / self.length, y / self.length]);
        }
        let rho = x.hypot(y);
        if rho == 0.0 || !rho.is_finite() || !z.is_finite() {
            return None;
        }
        let phi = y.atan2(x);
        let angular = |a: f64| (a / self.radians).rem_euclid(period);
        let shift = |sign: f64| if sign < 0.0 { PI } else { 0.0 };
        Some(match self.kind {
            ProjectionKind::Plane => unreachable!(),
            ProjectionKind::Cylinder(r) if r.is_finite() && r != 0.0 => {
                [angular(phi + shift(r)), z / self.length]
            }
            ProjectionKind::Cone(r, slope) if r.is_finite() && slope.is_finite() => {
                let inv_norm = 1.0 / slope.hypot(1.0);
                let sheet = |sign: f64| {
                    let h = (z * inv_norm + (sign * rho - r) * (slope * inv_norm)) * inv_norm;
                    let radial = r + slope * h;
                    (sign * radial > 0.0).then(|| {
                        (
                            [angular(phi + shift(sign)), h / self.length],
                            (rho - sign * radial).hypot(z - h),
                        )
                    })
                };
                closer(sheet(1.0), sheet(-1.0))?
            }
            ProjectionKind::Sphere(r) if r.is_finite() && r != 0.0 => [
                angular(phi + shift(r)),
                z.atan2(rho) / self.radians * r.signum(),
            ],
            ProjectionKind::Torus(major, minor)
                if major > 0.0 && minor > 0.0 && minor.is_finite() && major.is_finite() =>
            {
                let meridian = |sign: f64| {
                    let radial = sign * rho - major;
                    let norm = radial.hypot(z);
                    if norm == 0.0
                        || !norm.is_finite()
                        || sign * (major + minor * (radial / norm)) <= 0.0
                    {
                        return None;
                    }
                    Some((
                        [angular(phi + shift(sign)), angular(z.atan2(radial))],
                        (norm - minor).abs(),
                    ))
                };
                closer(
                    meridian(1.0),
                    if minor > major { meridian(-1.0) } else { None },
                )?
            }
            _ => return None,
        })
    }
}

fn closer(a: Option<([f64; 2], f64)>, b: Option<([f64; 2], f64)>) -> Option<[f64; 2]> {
    match (a, b) {
        (Some((a, d)), Some((b, e))) if d != e => Some(if d < e { a } else { b }),
        (Some((a, _)), None) | (None, Some((a, _))) => Some(a),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::super::{AngleUnit, ParsedStep};
    use super::*;
    use crate::tolerant::native::NativeSurface;

    #[test]
    fn analytic_samples_keep_enclosures_and_invert_in_source_units() {
        let bytes = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('analytic projection'),'2;1');
FILE_NAME('analytic.step','',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(4.1,-3.2,7.3));
#2=DIRECTION('',(1.,1.,0.));
#3=DIRECTION('',(0.,0.,1.));
#4=AXIS2_PLACEMENT_3D('',#1,#3,#2);
#5=PLANE('',#4);
#6=CYLINDRICAL_SURFACE('',#4,2.5);
#7=CONICAL_SURFACE('',#4,2.5,30.);
#8=SPHERICAL_SURFACE('',#4,2.5);
#9=TOROIDAL_SURFACE('',#4,4.,1.);
#10=TOROIDAL_SURFACE('',#4,1.,2.);
ENDSEC;
END-ISO-10303-21;
";
        let mut dec = Decider::default();
        let source = ParsedStep::read(bytes, &mut dec).unwrap();
        let length = Rational::parse("1.25").unwrap();
        for (id, uv) in [
            (5, [0.125, -1.5]),
            (6, [312.0, -1.5]),
            (7, [39.0, -8.0]),
            (8, [70.0, -22.0]),
            (9, [223.0, 262.0]),
            (10, [44.0, 175.0]),
        ] {
            let charts = source
                .surface_charts(id, &length, AngleUnit::Degrees, &mut dec)
                .unwrap();
            let box_uv = uv.map(|u| {
                [u - 1e-5, u + 1e-5].map(|value| Rational::from_f64_exact(value).unwrap())
            });
            let expected = charts[0].evaluate_box(&box_uv, 96, &mut dec).unwrap();
            let actual = charts[0].sample_box(&box_uv, 96, &mut dec).unwrap();
            for (a, b) in actual.iter().zip(expected) {
                assert_eq!([a.lo(), a.hi()], [b.lo(), b.hi()], "surface #{id}");
            }

            let surface =
                NativeSurface::from_source(&source, id, &length, AngleUnit::Degrees, &mut dec)
                    .unwrap();
            let sample = surface
                .sample(surface.location(uv).unwrap(), &mut dec)
                .unwrap();
            let fit = surface.project(sample.point_mm, None, &mut dec).unwrap();
            assert!(
                fit.converged && fit.iterations == 0,
                "surface #{id}: {fit:?}"
            );
            for (value, expected) in fit.location.parameters.into_iter().zip(uv) {
                assert!((value - expected).abs() < 1e-9, "surface #{id}: {fit:?}");
            }
        }
    }
}
