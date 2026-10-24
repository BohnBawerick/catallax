//! Exact single-nappe cells between two planes perpendicular to pull.

use crate::draft::{DraftDefinition, reduced, sign};
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
            heights[index] = reduced(
                plane.origin.minus(&cylinder.origin).dot(&plane.normal) / dot,
                dec,
            )?;
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
        let neutral_height = reduced(
            definition
                .neutral
                .origin
                .minus(&cylinder.origin)
                .dot(&definition.neutral.normal)
                / definition.pull.dot(&definition.neutral.normal),
            dec,
        )?;
        let radius = reduced(cylinder.radius_sq.sqrt(), dec)?;
        let tangent = reduced(&definition.sine / &definition.cosine, dec)?;
        let mut radius_at = |h: &Real| reduced(&radius - &(&tangent * &(h - &neutral_height)), dec);
        let mut radii = [radius_at(&heights[0])?, radius_at(&heights[1])?];
        let signs = [sign(&radii[0], dec)?, sign(&radii[1], dec)?];
        if signs.iter().all(|s| *s != Sign::Positive) {
            return Err(Refusal::new(
                "draft.reversed-solid",
                "both round-cell ends lie beyond the neutral-circle nappe",
            ));
        }
        for i in 0..2 {
            if signs[i] != Sign::Positive {
                heights[i] = reduced(&neutral_height + &(&radius / &tangent), dec)?;
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

    /// Proves closed separation for cells expressed along the same pull axis.
    /// A false result requires a curved arrangement, not an untrimmed assembly.
    pub fn certified_disjoint(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        if sign(&self.axis.cross(&other.axis).norm_squared(), dec)? != Sign::Zero
            || sign(&self.axis.dot(&other.axis), dec)? != Sign::Positive
        {
            return Ok(false);
        }
        let delta = other.origin.minus(&self.origin);
        let shift = delta.dot(&self.axis);
        let other_lower = &other.heights[0] + &shift;
        let lower = if sign(&(&other_lower - &self.heights[0]), dec)? == Sign::Positive {
            other_lower
        } else {
            self.heights[0].clone()
        };
        let other_upper = &other.heights[1] + &shift;
        let upper = if sign(&(&other_upper - &self.heights[1]), dec)? == Sign::Negative {
            other_upper
        } else {
            self.heights[1].clone()
        };
        let span = &upper - &lower;
        if sign(&span, dec)? == Sign::Negative {
            return Ok(true);
        }
        let a_rate = (&self.radii[1] - &self.radii[0]) / (&self.heights[1] - &self.heights[0]);
        let b_rate = (&other.radii[1] - &other.radii[0]) / (&other.heights[1] - &other.heights[0]);
        let a_radius = &self.radii[0] + &(&a_rate * &(&lower - &self.heights[0]));
        let b_height = &(&lower - &shift) - &other.heights[0];
        let b_radius = &other.radii[0] + &(&b_rate * &b_height);
        let lower_sum = &a_radius + &b_radius;
        let upper_sum = &lower_sum + &((&a_rate + &b_rate) * span);
        let maximum = if sign(&(&upper_sum - &lower_sum), dec)? == Sign::Positive {
            upper_sum
        } else {
            lower_sum
        };
        let lateral = delta.minus(&self.axis.scaled(&shift));
        Ok(sign(&(lateral.norm_squared() - maximum.square()), dec)? == Sign::Positive)
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

    /// A neutral section has area if it cuts the interior, or is a nonzero cap.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let normal_rate = self.axis.dot(&neutral.normal);
        let lateral = neutral.normal.minus(&self.axis.scaled(&normal_rate));
        let lateral_squared = lateral.norm_squared();
        if sign(&lateral_squared, dec)? != Sign::Zero {
            // The extrema of n·(p-N) on each end disk are centre ± |n_perp|r.
            // Both are affine in height, so their global extrema occur at ends.
            // A non-axial tangency contributes only a point or generator, no area.
            let spread_rate = lateral_squared.sqrt();
            let offset = self.origin.minus(&neutral.origin).dot(&neutral.normal);
            let mut below = false;
            let mut above = false;
            for i in 0..2 {
                let centre = &offset + &(&normal_rate * &self.heights[i]);
                let spread = &spread_rate * &self.radii[i];
                below |= sign(&(&centre - &spread), dec)? == Sign::Negative;
                above |= sign(&(&centre + &spread), dec)? == Sign::Positive;
            }
            return Ok(below && above);
        }
        let h = neutral.origin.minus(&self.origin).dot(&neutral.normal) / normal_rate;
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
