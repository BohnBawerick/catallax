//! Exact parallel supports. Trimming, envelopes and material selection belong to features.

use catallax_numbers::{Decider, Real, Refusal, Sign};

use crate::{Surface, math};

impl Surface {
    /// Moves a regular analytic support by `distance` along its outward unit normal.
    ///
    /// This constructs a support, not a trimmed offset solid. A missing support
    /// means the cylinder radius, sphere radius or torus tube was consumed.
    /// A cone's translated support still needs its original normal-image domain.
    /// Curved joins and singular envelopes must be constructed separately.
    pub fn parallel_support(
        &self,
        distance: &Real,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        math::validate(distance, dec)?;
        if math::sign(distance, dec)? == Sign::Zero {
            return Ok(Some(self.clone()));
        }
        let positive =
            |value: &Real, dec: &mut Decider| math::sign(value, dec).map(|s| s == Sign::Positive);
        Ok(Some(match self {
            Self::Plane(plane) => {
                let mut plane = plane.clone();
                let length = math::exact_sqrt(&plane.normal.norm_squared());
                plane.origin = plane
                    .origin
                    .plus(&plane.normal.scaled(&(distance / &length)));
                Self::Plane(plane)
            }
            Self::Cylinder(cylinder) => {
                let radius = &cylinder.radius() + distance;
                if !positive(&radius, dec)? {
                    return Ok(None);
                }
                let mut cylinder = cylinder.clone();
                cylinder.radius_sq = radius.square();
                Self::Cylinder(cylinder)
            }
            Self::Sphere(sphere) => {
                let radius = &sphere.radius + distance;
                if !positive(&radius, dec)? {
                    return Ok(None);
                }
                let mut sphere = sphere.clone();
                sphere.radius = radius;
                Self::Sphere(sphere)
            }
            Self::Cone(cone) => {
                let mut cone = cone.clone();
                let radial = distance * &math::exact_sqrt(&(&Real::one() + &cone.slope.square()));
                let axial = -(&radial / &cone.slope);
                let length = math::exact_sqrt(&cone.axis.norm_squared());
                cone.apex = cone.apex.plus(&cone.axis.scaled(&(&axial / &length)));
                Self::Cone(cone)
            }
            Self::Torus(torus) => {
                let minor = &torus.minor + distance;
                if !positive(&minor, dec)? {
                    return Ok(None);
                }
                if !positive(&(&torus.major - &minor), dec)? {
                    return Err(Refusal::new(
                        "offset.spindle-carrier",
                        "the distance envelope needs a horn or spindle torus domain, not a ring support",
                    ));
                }
                let mut torus = torus.clone();
                torus.minor = minor;
                Self::Torus(torus)
            }
            Self::TorusSheet(_) => {
                return Err(Refusal::new(
                    "offset.spindle-carrier",
                    "a parallel of a torus outer sheet needs its own bounded sheet domain",
                ));
            }
            Self::EllipticCylinder(_) | Self::BSpline(_) => {
                return Err(Refusal::new(
                    "offset.derived-carrier",
                    "this parallel needs a certified spline carrier and its retained definition",
                ));
            }
        }))
    }
}
