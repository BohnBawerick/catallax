//! One unchanged pull-parallel bore strictly inside a convex planar draft cell.
//!
//! Both kernels share this definition math. Neither passes topology through it.

use crate::draft::{DraftDefinition, cell::PlanarCell, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Plane, Point3, Vector3};

/// A convex planar envelope with one finite bore removed from cap to cap.
/// The bore's support and both ring trims never change.
#[derive(Clone, Debug)]
pub struct BoredCell {
    /// Filled outer material, before the bore is removed.
    pub outer: PlanarCell,
    /// Unchanged bore wall carrier.
    pub cylinder: Cylinder,
    /// Complete source names of the bore wall.
    pub origins: Vec<Origin>,
    /// Lower and upper cap support indices along the cylinder axis.
    pub caps: [usize; 2],
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-arrangement", detail)
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}

/// The cylinder axis point on a plane perpendicular to that axis.
pub fn axis_point(cylinder: &Cylinder, plane: &Plane) -> Point3 {
    let height = -&(value(plane, &cylinder.origin) / plane.normal.dot(&cylinder.axis));
    cylinder.origin.plus(&cylinder.axis.scaled(&height))
}

/// An unselected wall has no neutral curve, but this bounded family keeps its
/// axial mouths parallel to the neutral section: pull, axis and normal agree.
pub fn check_axial(
    cylinder: &Cylinder,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let parallel = |a: &Vector3, dec: &mut Decider| {
        Ok::<_, Refusal>(sign(&a.cross(&definition.pull).norm_squared(), dec)? == Sign::Zero)
    };
    if !parallel(&cylinder.axis, dec)? || !parallel(&definition.neutral.normal, dec)? {
        return Err(unsupported(
            "a preserved bore needs its axis and the neutral normal along pull",
        ));
    }
    Ok(())
}

/// Finds the two opposed caps and proves the finite cylinder between them is
/// strictly inside every other halfspace. A wall value is affine along the axis
/// and its circular maximum is the axis value plus radius times the normal's
/// radial length, so strict inequalities at both rings cover the whole height.
pub fn certify_clear<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinder: &Cylinder,
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
    let axis_sq = cylinder.axis.norm_squared();
    let mut lower = None;
    let mut upper = None;
    for (index, plane) in planes.clone().enumerate() {
        if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
            continue;
        }
        let side = if sign(&plane.normal.dot(&cylinder.axis), dec)? == Sign::Positive {
            &mut upper
        } else {
            &mut lower
        };
        if side.replace(index).is_some() {
            return Err(unsupported(
                "a preserved bore needs exactly two opposed caps",
            ));
        }
    }
    let (Some(lower), Some(upper)) = (lower, upper) else {
        return Err(unsupported(
            "a preserved bore needs exactly two opposed caps",
        ));
    };
    let caps = [lower, upper];
    let plane = |index: usize| planes.clone().nth(index).expect("indexed above");
    let centres = caps.map(|i| axis_point(cylinder, plane(i)));
    if sign(&centres[1].minus(&centres[0]).dot(&cylinder.axis), dec)? != Sign::Positive {
        return Err(unsupported("preserved bore caps have no axial interval"));
    }
    for (index, wall) in planes.enumerate() {
        if caps.contains(&index) {
            continue;
        }
        let along = wall.normal.dot(&cylinder.axis);
        let radial = wall.normal.norm_squared() - along.square() / axis_sq.clone();
        for centre in &centres {
            let offset = value(wall, centre);
            if sign(&offset, dec)? != Sign::Negative
                || sign(&(offset.square() - &cylinder.radius_sq * &radial), dec)? != Sign::Positive
            {
                return Err(unsupported(
                    "the bore touches or crosses a wall of its planar envelope",
                ));
            }
        }
    }
    Ok(caps)
}

impl BoredCell {
    /// Certifies the bore against the replaced envelope before any publication.
    pub fn new(
        outer: PlanarCell,
        cylinder: Cylinder,
        origins: Vec<Origin>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        check_axial(&cylinder, definition, dec)?;
        let caps = certify_clear(outer.supports.iter().map(|s| &s.plane), &cylinder, dec)?;
        Ok(Self {
            outer,
            cylinder,
            origins,
            caps,
        })
    }

    /// Exact strict material membership, excluding the bore and every boundary.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.outer.contains_interior(point, dec)?
            && sign(
                &(point
                    .minus(&self.cylinder.origin)
                    .cross(&self.cylinder.axis)
                    .norm_squared()
                    - &self.cylinder.radius_sq * &self.cylinder.axis.norm_squared()),
                dec,
            )? == Sign::Positive)
    }
}
