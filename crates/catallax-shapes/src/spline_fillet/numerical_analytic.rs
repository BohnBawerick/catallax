//! Binary64 proposals in the same native frames as Surface::evaluate.
//! These values carry no geometric authority; accepted stations use Surface's
//! interval evaluator and interval native derivatives independently.
//!
//! Sines, cosines, arctangents and hypotenuses in binary64 go through the
//! `libm` crate (the musl library's routines in Rust), not the methods of
//! `f64`, which call the platform's library and differ in the last bit from
//! one machine to the next. The work that follows from these values is
//! counted, and the predictor's done-when lines ask for the same work on two
//! machines (Checkpoint 2 decisions, section 4.5).

use super::*;
use crate::Surface;

#[derive(Debug)]
enum Profile {
    Plane,
    Cylinder(f64),
    Cone {
        slope: f64,
        inverse_normal_length: f64,
    },
    Sphere(f64),
    Torus {
        major: f64,
        minor: f64,
    },
    Ellipse {
        major: f64,
        minor: f64,
    },
}

#[derive(Debug)]
pub(super) struct Analytic64 {
    origin: [f64; 3],
    frame: [[f64; 3]; 3],
    profile: Profile,
}

fn scalar(value: &Real, dec: &mut Decider) -> Result<f64, Refusal> {
    let value = if let Some(value) = value.as_rational() {
        value.to_f64_nearest()
    } else {
        dec.enclose_budgeted_at(value, 128)?.mid().to_f64_nearest()
    };
    if !value.is_finite() {
        return Err(bad("native rolling-ball scalar is outside binary64 range"));
    }
    Ok(value)
}

fn unit(value: [f64; 3]) -> Result<[f64; 3], Refusal> {
    let length = norm(value);
    if !length.is_finite() || length <= 0.0 {
        return Err(bad("native rolling-ball frame is outside binary64 range"));
    }
    Ok(value.map(|coordinate| coordinate / length))
}

impl Analytic64 {
    pub(super) fn new(surface: &Surface, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let profile = match surface {
            Surface::Plane(_) => Profile::Plane,
            Surface::Cylinder(value) => Profile::Cylinder(scalar(&value.radius(), dec)?),
            Surface::Cone(value) => {
                let slope = scalar(&value.slope, dec)?;
                Profile::Cone {
                    slope,
                    inverse_normal_length: 1.0 / (1.0 + slope * slope).sqrt(),
                }
            }
            Surface::Sphere(value) => Profile::Sphere(scalar(&value.radius, dec)?),
            Surface::Torus(value) => Profile::Torus {
                major: scalar(&value.major, dec)?,
                minor: scalar(&value.minor, dec)?,
            },
            Surface::TorusSheet(value) => Profile::Torus {
                major: scalar(&value.major, dec)?,
                minor: scalar(&value.minor, dec)?,
            },
            Surface::EllipticCylinder(value) => Profile::Ellipse {
                major: scalar(&value.major_sq.sqrt(), dec)?,
                minor: scalar(&value.minor_sq.sqrt(), dec)?,
            },
            Surface::BSpline(_) => return Ok(None),
        };
        dec.charge_work(128)?;
        let (origin, axis, reference) = surface
            .definition_frame()
            .ok_or_else(|| bad("native support has no analytic frame"))?;
        let e3 = unit([
            scalar(&axis.x, dec)?,
            scalar(&axis.y, dec)?,
            scalar(&axis.z, dec)?,
        ])?;
        let e1 = unit([
            scalar(&reference.x, dec)?,
            scalar(&reference.y, dec)?,
            scalar(&reference.z, dec)?,
        ])?;
        let e2 = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            e3[j] * e1[k] - e3[k] * e1[j]
        });
        let origin = [
            scalar(&origin.x, dec)?,
            scalar(&origin.y, dec)?,
            scalar(&origin.z, dec)?,
        ];
        Ok(Some(Self {
            origin,
            frame: [e1, e2, e3],
            profile,
        }))
    }

    pub(super) fn contact(&self, [u, v]: [f64; 2]) -> Result<([f64; 3], [f64; 3]), Refusal> {
        let (point, normal) = match self.profile {
            Profile::Plane => ([u, v, 0.0], [0.0, 0.0, 1.0]),
            Profile::Cylinder(radius) => {
                let (sin, cos) = libm::sincos(u);
                ([radius * cos, radius * sin, v], [cos, sin, 0.0])
            }
            Profile::Cone {
                slope,
                inverse_normal_length: inverse,
            } => {
                if v <= 0.0 {
                    return Err(bad(
                        "rolling-ball contact leaves the cone's positive-height chart",
                    ));
                }
                let (sin, cos) = libm::sincos(u);
                (
                    [v * slope * cos, v * slope * sin, v],
                    [cos * inverse, sin * inverse, -slope * inverse],
                )
            }
            Profile::Sphere(radius) => {
                if v.abs() >= std::f64::consts::FRAC_PI_2 {
                    return Err(bad("rolling-ball contact leaves the regular sphere chart"));
                }
                let (su, cu) = libm::sincos(u);
                let (sv, cv) = libm::sincos(v);
                let normal = [cv * cu, cv * su, sv];
                (normal.map(|value| radius * value), normal)
            }
            Profile::Torus { major, minor } => {
                let (su, cu) = libm::sincos(u);
                let (sv, cv) = libm::sincos(v);
                let radial = major + minor * cv;
                if radial <= 0.0 {
                    return Err(bad(
                        "rolling-ball contact leaves the torus's regular outer sheet",
                    ));
                }
                (
                    [radial * cu, radial * su, minor * sv],
                    [cv * cu, cv * su, sv],
                )
            }
            Profile::Ellipse { major, minor } => {
                let (sin, cos) = libm::sincos(u);
                (
                    [major * cos, minor * sin, v],
                    unit([cos / major, sin / minor, 0.0])?,
                )
            }
        };
        let point = std::array::from_fn(|i| {
            self.origin[i] + (0..3).map(|j| point[j] * self.frame[j][i]).sum::<f64>()
        });
        let normal =
            std::array::from_fn(|i| (0..3).map(|j| normal[j] * self.frame[j][i]).sum::<f64>());
        if !point.iter().chain(&normal).all(|value| value.is_finite()) {
            return Err(bad("non-finite native rolling-ball proposal"));
        }
        Ok((point, normal))
    }
}
