//! Direct point distances to complete circular cylinder bands and perpendicular
//! circular caps. The finite closed face is measured, never its infinite support.
//!
//! Recognition is exact and conservative: anything that is not a full band between
//! two perpendicular coaxial rim circles, or a full disc bounded by one circle
//! perpendicular to its plane, is left to the general complete search.

use super::{Bounds, Options, validate, witness_distance};
use crate::{Circle, Point3, Surface, Vector3, math};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};

/// Which nearest-point case produced a witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// Interior of the band wall or of the cap disc.
    Interior,
    /// A rim circle, reached from outside the finite extent or from beyond the edge.
    Rim,
    /// The point lies on the cylinder axis, so every point of the nearest circle
    /// ties and the seam point is reported.
    Axis,
}

/// A recognized finite analytic face.
#[derive(Clone, Debug)]
pub struct Finite(Shape);
#[derive(Clone, Debug)]
enum Shape {
    /// A full cylinder wall between two axial heights.
    Band {
        origin: Point3,
        axis: Vector3,
        reference: Vector3,
        radius: Real,
        heights: [Real; 2],
    },
    /// A full disc in a plane.
    Cap {
        center: Point3,
        normal: Vector3,
        radius_sq: Real,
        radius: Real,
    },
}

fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(Real::one() / v.norm_squared().sqrt()))
}
fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(math::sign(value, dec)? == Sign::Zero)
}
/// The rim circle lies in a plane perpendicular to `axis`, centred on the line
/// through `origin` along `axis`.
fn coaxial(
    rim: &Circle,
    origin: &Point3,
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(zero(&rim.normal.cross(axis).norm_squared(), dec)?
        && zero(&rim.center.minus(origin).cross(axis).norm_squared(), dec)?)
}

impl Finite {
    /// Recognizes `surface` bounded exactly by `rims`, the complete set of one
    /// closed circle per boundary loop. Returns `None` for every other trim.
    pub fn recognize(
        surface: &Surface,
        rims: &[&Circle],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        match (surface, rims) {
            (Surface::Cylinder(s), [a, b]) => {
                for rim in [a, b] {
                    if !coaxial(rim, &s.origin, &s.axis, dec)?
                        || !zero(&(&rim.radius_sq - &s.radius_sq), dec)?
                    {
                        return Ok(None);
                    }
                }
                let axis = unit(&s.axis);
                let h = |rim: &Circle| rim.center.minus(&s.origin).dot(&axis);
                let heights = [h(a), h(b)];
                let heights = match math::sign(&(&heights[1] - &heights[0]), dec)? {
                    Sign::Zero => return Ok(None),
                    Sign::Positive => heights,
                    Sign::Negative => [heights[1].clone(), heights[0].clone()],
                };
                Ok(Some(Self(Shape::Band {
                    origin: s.origin.clone(),
                    reference: unit(&s.ref_dir),
                    radius: s.radius(),
                    axis,
                    heights,
                })))
            }
            (Surface::Plane(s), [rim]) => {
                if !zero(&rim.normal.cross(&s.normal).norm_squared(), dec)?
                    || !zero(&rim.center.minus(&s.origin).dot(&s.normal), dec)?
                {
                    return Ok(None);
                }
                Ok(Some(Self(Shape::Cap {
                    center: rim.center.clone(),
                    normal: unit(&s.normal),
                    radius: math::exact_sqrt(&rim.radius_sq),
                    radius_sq: rim.radius_sq.clone(),
                })))
            }
            _ => Ok(None),
        }
    }

    /// Nearest point of the closed finite face to `p` as exact coordinates, the
    /// squared distance, and the branch taken.
    fn nearest(&self, p: &Point3, dec: &mut Decider) -> Result<(Point3, Real, Branch), Refusal> {
        match &self.0 {
            Shape::Band {
                origin,
                axis,
                reference,
                radius,
                heights,
            } => {
                let d = p.minus(origin);
                let h = d.dot(axis);
                let planar = d.minus(&axis.scaled(&h));
                let rho_sq = planar.norm_squared();
                let on_axis = zero(&rho_sq, dec)?;
                let (clamped, wall) = if math::sign(&(&h - &heights[0]), dec)? == Sign::Negative {
                    (heights[0].clone(), false)
                } else if math::sign(&(&h - &heights[1]), dec)? == Sign::Positive {
                    (heights[1].clone(), false)
                } else {
                    (h.clone(), true)
                };
                let (radial, rho) = if on_axis {
                    (reference.clone(), Real::zero())
                } else {
                    let rho = rho_sq.sqrt();
                    (planar.scaled(&(Real::one() / rho.clone())), rho)
                };
                let witness = origin
                    .plus(&axis.scaled(&clamped))
                    .plus(&radial.scaled(radius));
                let squared = (&rho - radius).square() + (&h - &clamped).square();
                let branch = if on_axis {
                    Branch::Axis
                } else if wall {
                    Branch::Interior
                } else {
                    Branch::Rim
                };
                Ok((witness, squared, branch))
            }
            Shape::Cap {
                center,
                normal,
                radius_sq,
                radius,
            } => {
                let d = p.minus(center);
                let h = d.dot(normal);
                let planar = d.minus(&normal.scaled(&h));
                let rho_sq = planar.norm_squared();
                if math::sign(&(&rho_sq - radius_sq), dec)? != Sign::Positive {
                    return Ok((p.plus(&normal.scaled(&(-&h))), h.square(), Branch::Interior));
                }
                let rho = rho_sq.sqrt();
                let witness = center.plus(&planar.scaled(&(radius / &rho)));
                Ok((witness, (&rho - radius).square() + h.square(), Branch::Rim))
            }
        }
    }

    /// Bounds the distance from `p` to the finite face directly. Returns `None`
    /// when the requested width is not reached, so the caller can fall back to the
    /// general search. The lower bound encloses the exact distance and the upper
    /// bound is measured to an exact witness on the closed face.
    pub fn point_distance(
        &self,
        p: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Option<(Bounds, Branch)>, Refusal> {
        validate(options)?;
        let (witness, squared, branch) = self.nearest(p, dec)?;
        let lower = dec
            .enclose_budgeted_at(&squared, super::BITS)?
            .sqrt()
            .unwrap()
            .lo()
            .max(Rational::zero());
        let upper = witness_distance(p, &witness, dec)?;
        if upper < lower || &upper - &lower > options.width {
            return Ok(None);
        }
        Ok(Some((
            Bounds {
                lower,
                upper: Some(upper),
                witnesses: Some([p.clone(), witness]),
                converged: true,
                cells: 0,
            },
            branch,
        )))
    }
}
