//! Unchanged pull-parallel bores strictly inside a convex planar draft cell.
//!
//! Both kernels share this definition math. Neither passes topology through it.

use crate::draft::{DraftDefinition, cell::PlanarCell, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Plane, Point3, Vector3};

/// One finite bore's unchanged wall.
#[derive(Clone, Debug)]
pub struct Bore {
    /// Unchanged bore wall carrier.
    pub cylinder: Cylinder,
    /// Complete source names of the bore wall.
    pub origins: Vec<Origin>,
}

/// A convex planar envelope with finite bores removed from cap to cap.
/// Every bore's support and both ring trims never change.
#[derive(Clone, Debug)]
pub struct BoredCell {
    /// Filled outer material, before any bore is removed.
    pub outer: PlanarCell,
    /// Mutually separated bores, all between the same two caps.
    pub bores: Vec<Bore>,
    /// The two cap support indices, in ascending order.
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

/// Whether a circle about the axis at `centre` lies strictly inside a halfspace.
/// The circular maximum is the axis value plus radius times the normal's radial
/// length, compared on squares so the proof needs no root.
pub fn ring_clear(
    wall: &Plane,
    cylinder: &Cylinder,
    centre: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let along = wall.normal.dot(&cylinder.axis);
    let radial = wall.normal.norm_squared() - along.square() / cylinder.axis.norm_squared();
    let offset = value(wall, centre);
    Ok(sign(&offset, dec)? == Sign::Negative
        && sign(&(offset.square() - &cylinder.radius_sq * &radial), dec)? == Sign::Positive)
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
            "a preserved wall needs its axis and the neutral normal along pull",
        ));
    }
    Ok(())
}

/// Finds the two opposed caps and proves the finite cylinder between them is
/// strictly inside every other halfspace. A wall value is affine along the axis,
/// so strict inequalities at both rings cover the whole height.
pub fn certify_clear<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinder: &Cylinder,
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
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
        for centre in &centres {
            if !ring_clear(wall, cylinder, centre, dec)? {
                return Err(unsupported(
                    "the bore touches or crosses a wall of its planar envelope",
                ));
            }
        }
    }
    Ok(caps)
}

/// Proves two walls strictly separated. Parallel walls are nearest along their
/// common normal, so the axis distance `d` must exceed the radii's sum. Squaring
/// `d^2-a^2-b^2 > 2ab` twice keeps the proof exact without a root. An exact
/// contact refuses, and a gap of any positive size stands unchanged.
pub fn certify_disjoint(a: &Cylinder, b: &Cylinder, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(&a.axis.cross(&b.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(unsupported("preserved walls need parallel axes"));
    }
    let offset = b.origin.minus(&a.origin).cross(&a.axis).norm_squared() / a.axis.norm_squared();
    let spare = offset - a.radius_sq.clone() - b.radius_sq.clone();
    let reach = Real::from_i64(4) * (&a.radius_sq * &b.radius_sq);
    if sign(&spare, dec)? != Sign::Positive
        || sign(&(spare.square() - reach), dec)? != Sign::Positive
    {
        return Err(unsupported("preserved walls touch or cross each other"));
    }
    Ok(())
}

/// Proves every bore clear of every other halfspace, between the same two caps
/// and strictly apart from each other. Opposite axis charts of one direction
/// agree, because their caps are compared as a set. Returns those caps ascending.
pub fn certify_bores<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinders: &[&Cylinder],
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
    let mut shared = None;
    for (index, cylinder) in cylinders.iter().enumerate() {
        let mut caps = certify_clear(planes.clone(), cylinder, dec)?;
        caps.sort_unstable();
        if shared.replace(caps).is_some_and(|other| other != caps) {
            return Err(unsupported("preserved bores need the same two caps"));
        }
        for other in &cylinders[..index] {
            certify_disjoint(other, cylinder, dec)?;
        }
    }
    shared.ok_or_else(|| unsupported("a preserved bore needs a wall"))
}

impl BoredCell {
    /// Certifies the bores against the replaced envelope before any publication.
    pub fn new(
        outer: PlanarCell,
        bores: Vec<Bore>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        for bore in &bores {
            check_axial(&bore.cylinder, definition, dec)?;
        }
        let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
        let caps = certify_bores(outer.supports.iter().map(|s| &s.plane), &cylinders, dec)?;
        Ok(Self { outer, bores, caps })
    }

    /// Exact strict material membership, excluding every bore and boundary.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if !self.outer.contains_interior(point, dec)? {
            return Ok(false);
        }
        for Bore { cylinder, .. } in &self.bores {
            let radial = point.minus(&cylinder.origin).cross(&cylinder.axis);
            let excess =
                radial.norm_squared() - &cylinder.radius_sq * &cylinder.axis.norm_squared();
            if sign(&excess, dec)? != Sign::Positive {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
