//! Binary64 proposals in the same native frames as Surface::evaluate.
//! These values carry no geometric authority; accepted stations use Surface's
//! interval evaluator and interval native derivatives independently.

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
            Surface::Cylinder(value) => Profile::Cylinder(value.radius().to_f64()),
            Surface::Cone(value) => {
                let slope = value.slope.to_f64();
                Profile::Cone {
                    slope,
                    inverse_normal_length: 1.0 / (1.0 + slope * slope).sqrt(),
                }
            }
            Surface::Sphere(value) => Profile::Sphere(value.radius.to_f64()),
            Surface::Torus(value) => Profile::Torus {
                major: value.major.to_f64(),
                minor: value.minor.to_f64(),
            },
            Surface::TorusSheet(value) => Profile::Torus {
                major: value.major.to_f64(),
                minor: value.minor.to_f64(),
            },
            Surface::EllipticCylinder(value) => Profile::Ellipse {
                major: value.major_sq.sqrt().to_f64(),
                minor: value.minor_sq.sqrt().to_f64(),
            },
            Surface::BSpline(_) => return Ok(None),
        };
        dec.charge_work(128)?;
        let (origin, axis, reference) = surface
            .definition_frame()
            .ok_or_else(|| bad("native support has no analytic frame"))?;
        let e3 = unit([axis.x.to_f64(), axis.y.to_f64(), axis.z.to_f64()])?;
        let e1 = unit([
            reference.x.to_f64(),
            reference.y.to_f64(),
            reference.z.to_f64(),
        ])?;
        let e2 = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            e3[j] * e1[k] - e3[k] * e1[j]
        });
        let origin = [origin.x.to_f64(), origin.y.to_f64(), origin.z.to_f64()];
        if !origin.iter().all(|value| value.is_finite()) {
            return Err(bad("native rolling-ball origin is outside binary64 range"));
        }
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
                let (sin, cos) = u.sin_cos();
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
                let (sin, cos) = u.sin_cos();
                (
                    [v * slope * cos, v * slope * sin, v],
                    [cos * inverse, sin * inverse, -slope * inverse],
                )
            }
            Profile::Sphere(radius) => {
                if v.abs() >= std::f64::consts::FRAC_PI_2 {
                    return Err(bad("rolling-ball contact leaves the regular sphere chart"));
                }
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let normal = [cv * cu, cv * su, sv];
                (normal.map(|value| radius * value), normal)
            }
            Profile::Torus { major, minor } => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
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
                let (sin, cos) = u.sin_cos();
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
