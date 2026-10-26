//! Direct point distances to complete circular cylinder bands and perpendicular
//! circular discs, spherical caps and coaxial spherical bands. The finite closed
//! face is measured, never its infinite support.
//!
//! Recognition proves complete circular boundaries and, for spheres, their
//! oriented retained sides. Every unrecognized trim keeps the complete search.

use super::{Bounds, Options, validate, witness_distance};
use crate::{Circle, Point3, Sphere, Surface, Vector3, math};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};

/// Which nearest-point case produced a witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// The radial sphere point, band wall or cap disc, including boundary ties.
    Interior,
    /// A rim circle, reached from outside the finite extent or from beyond the edge.
    Rim,
    /// The query lies on the symmetry axis (including the sphere centre).
    /// A nearest circle tie uses its stored reference direction.
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
    /// A closed spherical height interval proven by oriented full circles.
    Sphere {
        center: Point3,
        axis: Vector3,
        radius: Real,
        heights: [Real; 2],
        rims: Vec<SphericalRim>,
    },
}

#[derive(Clone, Debug)]
struct SphericalRim {
    height: Real,
    radius: Real,
    reference: Vector3,
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
    /// Recognizes complete circular boundaries with retained-side metadata.
    /// Each boolean is true when travel follows the circle's natural direction
    /// relative to the carrier's outward normal. A topology adapter must combine
    /// edge, coedge and face senses. On a sphere the retained side is then toward
    /// the circle normal; reversing that boolean selects the complementary cap.
    /// Only one cap or a compatible two-rim coaxial band is recognized. Other
    /// spherical trims return `None`; other carriers use [`Self::recognize`].
    pub fn recognize_oriented(
        surface: &Surface,
        rims: &[(&Circle, bool)],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if let Surface::Sphere(s) = surface {
            return Self::spherical(s, rims, dec);
        }
        Self::recognize(
            surface,
            &rims.iter().map(|(rim, _)| *rim).collect::<Vec<_>>(),
            dec,
        )
    }

    fn spherical(
        sphere: &Sphere,
        rims: &[(&Circle, bool)],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if !(1..=2).contains(&rims.len()) || zero(&rims[0].0.normal.norm_squared(), dec)? {
            return Ok(None);
        }
        let axis = unit(&rims[0].0.normal);
        let mut heights = [-&sphere.radius, sphere.radius.clone()];
        let mut lower = false;
        let mut upper = false;
        let mut circles = Vec::new();
        for &(rim, forward) in rims {
            if zero(&rim.normal.norm_squared(), dec)?
                || zero(&rim.ref_dir.norm_squared(), dec)?
                || math::sign(&rim.radius_sq, dec)? != Sign::Positive
                || !coaxial(rim, &sphere.center, &axis, dec)?
                || !zero(&rim.ref_dir.dot(&axis), dec)?
            {
                return Ok(None);
            }
            let height = rim.center.minus(&sphere.center).dot(&axis);
            if !zero(
                &(&rim.radius_sq + &height.square() - sphere.radius.square()),
                dec,
            )? {
                return Ok(None);
            }
            // With outward sphere normal n and natural tangent t, n cross t
            // points toward the circle normal. The effective travel fixes which
            // of the two spherical components is retained.
            let increasing =
                (math::sign(&rim.normal.dot(&axis), dec)? == Sign::Positive) == forward;
            if increasing {
                if lower {
                    return Ok(None);
                }
                lower = true;
                heights[0] = height.clone();
            } else {
                if upper {
                    return Ok(None);
                }
                upper = true;
                heights[1] = height.clone();
            }
            circles.push(SphericalRim {
                height,
                radius: math::exact_sqrt(&rim.radius_sq),
                reference: unit(&rim.ref_dir),
            });
        }
        if math::sign(&(&heights[1] - &heights[0]), dec)? != Sign::Positive {
            return Ok(None);
        }
        Ok(Some(Self(Shape::Sphere {
            center: sphere.center.clone(),
            axis,
            radius: sphere.radius.clone(),
            heights,
            rims: circles,
        })))
    }

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
            Shape::Sphere {
                center,
                axis,
                radius,
                heights,
                rims,
            } => {
                let d = p.minus(center);
                let norm_sq = d.norm_squared();
                if zero(&norm_sq, dec)? {
                    let rim = &rims[0];
                    let witness = center
                        .plus(&axis.scaled(&rim.height))
                        .plus(&rim.reference.scaled(&rim.radius));
                    return Ok((witness, radius.square(), Branch::Axis));
                }
                let norm = math::exact_sqrt(&norm_sq);
                let radial = d.scaled(&(radius / &norm));
                let height = radial.dot(axis);
                let z = d.dot(axis);
                let planar = d.minus(&axis.scaled(&z));
                let rho_sq = planar.norm_squared();
                let on_axis = zero(&rho_sq, dec)?;
                let branch = if on_axis {
                    Branch::Axis
                } else {
                    Branch::Interior
                };
                if math::sign(&(&height - &heights[0]), dec)? != Sign::Negative
                    && math::sign(&(&height - &heights[1]), dec)? != Sign::Positive
                {
                    return Ok((center.plus(&radial), (&norm - radius).square(), branch));
                }
                // Squared distance is |d|²+R²-2*d·x. At each retained height,
                // azimuth aligns with planar(d). Its maximum over height is the
                // admitted radial point or a finite rim. On-axis it is linear;
                // at the centre every admitted point ties (handled above).
                let rho = math::exact_sqrt(&rho_sq);
                let mut best: Option<(Point3, Real, Branch)> = None;
                for rim in rims {
                    let direction = if on_axis {
                        rim.reference.clone()
                    } else {
                        planar.scaled(&(&Real::one() / &rho))
                    };
                    let witness = center
                        .plus(&axis.scaled(&rim.height))
                        .plus(&direction.scaled(&rim.radius));
                    let squared = (&rho - &rim.radius).square() + (&z - &rim.height).square();
                    if let Some((_, old, _)) = &best
                        && math::sign(&(&squared - old), dec)? != Sign::Negative
                    {
                        continue;
                    }
                    best = Some((
                        witness,
                        squared,
                        if on_axis { Branch::Axis } else { Branch::Rim },
                    ));
                }
                Ok(best.expect("a recognized spherical region has at least one rim"))
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
