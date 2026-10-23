//! Exact single-nappe cells between two planes perpendicular to pull.

use crate::draft::{DraftDefinition, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Plane, Point3, Surface, Vector3};

/// An axial cylindrical draft, with neighbours trimmed at the apex when needed.
#[derive(Clone, Debug)]
pub struct RoundCell {
    /// Axis reference point of the input cylinder.
    pub origin: Point3,
    /// Unit axis in the pull direction.
    pub axis: Vector3,
    /// Exact carrier, including its chosen nappe.
    pub surface: Surface,
    /// Side-face origin set, unchanged by draft.
    pub origins: Vec<Origin>,
    /// Axial coordinates of the surviving lower and upper boundaries.
    pub heights: [Real; 2],
    /// Radii at those boundaries. Zero denotes an apex, not a cap.
    pub radii: [Real; 2],
    /// Surviving cap definitions and origins. An apex has no cap.
    pub caps: [Option<(Plane, Vec<Origin>)>; 2],
    /// Perpendicular seam direction, with its original exact scale.
    pub reference: Vector3,
}

impl RoundCell {
    /// Constructs the exact bounded material on the neutral-circle nappe.
    /// Caps may be swallowed; no replacement cap is invented at an apex.
    pub fn new(
        cylinder: &Cylinder,
        origins: Vec<Origin>,
        caps: [(Plane, Vec<Origin>); 2],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let support = definition.support(&Surface::Cylinder(cylinder.clone()), true, dec)?;
        let axis = definition.pull.clone();
        let mut heights = [Real::zero(), Real::zero()];
        let mut ordered: [Option<(Plane, Vec<Origin>)>; 2] = [None, None];
        for (plane, origins) in caps {
            if sign(&plane.normal.cross(&axis).norm_squared(), dec)? != Sign::Zero {
                return Err(Refusal::new(
                    "draft.unsupported-round-cap",
                    "non-axial neighbouring planes need conic trimmed-cell construction",
                ));
            }
            let dot = plane.normal.dot(&axis);
            let index = if sign(&dot, dec)? == Sign::Negative {
                0
            } else {
                1
            };
            if ordered[index].is_some() {
                return Err(Refusal::new(
                    "draft.nonclosing-boundary",
                    "round cell has no pair of opposing caps",
                ));
            }
            heights[index] = plane.origin.minus(&cylinder.origin).dot(&plane.normal) / dot;
            ordered[index] = Some((plane, origins));
        }
        if ordered.iter().any(Option::is_none)
            || sign(&(&heights[1] - &heights[0]), dec)? != Sign::Positive
        {
            return Err(Refusal::new(
                "draft.reversed-solid",
                "round cell has no positive axial interior",
            ));
        }
        let neutral_height = definition
            .neutral
            .origin
            .minus(&cylinder.origin)
            .dot(&definition.neutral.normal)
            / definition.pull.dot(&definition.neutral.normal);
        let radius = cylinder.radius_sq.sqrt();
        let tangent = &definition.sine / &definition.cosine;
        let mut radii = heights
            .each_ref()
            .map(|h| &radius - &(&tangent * &(h - &neutral_height)));
        let signs = [sign(&radii[0], dec)?, sign(&radii[1], dec)?];
        if signs.iter().all(|s| *s != Sign::Positive) {
            return Err(Refusal::new(
                "draft.reversed-solid",
                "both round-cell ends lie beyond the neutral-circle nappe",
            ));
        }
        for i in 0..2 {
            if signs[i] != Sign::Positive {
                heights[i] = &neutral_height + &(&radius / &tangent);
                radii[i] = Real::zero();
                ordered[i] = None;
            }
        }
        Ok(Self {
            origin: cylinder.origin.clone(),
            axis,
            surface: support.surface,
            origins,
            heights,
            radii,
            caps: ordered,
            reference: cylinder.ref_dir.clone(),
        })
    }

    /// Centre of a boundary circle or apex.
    pub fn centre(&self, index: usize) -> Point3 {
        self.origin.plus(&self.axis.scaled(&self.heights[index]))
    }

    /// Exact strict membership. A cone's other nappe is never included.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = point.minus(&self.origin);
        let h = d.dot(&self.axis);
        if sign(&(&h - &self.heights[0]), dec)? != Sign::Positive
            || sign(&(&self.heights[1] - &h), dec)? != Sign::Positive
        {
            return Ok(false);
        }
        let t = (&h - &self.heights[0]) / (&self.heights[1] - &self.heights[0]);
        let radius = &self.radii[0] + &(t * (&self.radii[1] - &self.radii[0]));
        Ok(sign(&(d.cross(&self.axis).norm_squared() - radius.square()), dec)? == Sign::Negative)
    }

    /// The neutral plane must cut a positive-radius section, not only an apex.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let h = neutral.origin.minus(&self.origin).dot(&neutral.normal)
            / self.axis.dot(&neutral.normal);
        if sign(&(&h - &self.heights[0]), dec)? == Sign::Negative
            || sign(&(&self.heights[1] - &h), dec)? == Sign::Negative
        {
            return Ok(false);
        }
        let t = (&h - &self.heights[0]) / (&self.heights[1] - &self.heights[0]);
        let radius = &self.radii[0] + &(t * (&self.radii[1] - &self.radii[0]));
        Ok(sign(&radius, dec)? == Sign::Positive)
    }
}
