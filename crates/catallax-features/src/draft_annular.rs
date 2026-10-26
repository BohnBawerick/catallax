//! A finite hollow axial cell; both radial graphs stay strictly separated.

use crate::draft::{DraftDefinition, reduced, round::RoundCell, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Plane, Point3, Surface};

/// One complete source cylinder, with its material sense and saved selection.
#[derive(Clone, Debug)]
pub struct AnnularWall {
    /// Exact source carrier.
    pub cylinder: Cylinder,
    /// Complete source face names.
    pub origins: Vec<Origin>,
    /// True for the exterior wall, false for the wall of the hole.
    pub same_sense: bool,
    /// Whether this carrier is replaced by the draft law.
    pub selected: bool,
}

/// One connected annulus between two perpendicular planes.
/// The inner round cell describes the excluded hole, never material.
#[derive(Clone, Debug)]
pub struct AnnularCell {
    /// Complete finite exterior radial graph.
    pub outer: RoundCell,
    /// Complete finite hole radial graph, expressed in the same axial frame.
    pub inner: RoundCell,
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-arrangement", detail)
}

impl AnnularCell {
    /// Certifies the bounded annular slice before publishing any topology.
    /// Folds, contacts, crossings and swallowed ends need a curved arrangement.
    pub fn new(
        mut walls: [AnnularWall; 2],
        caps: [(Plane, Vec<Origin>); 2],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if !walls[0].same_sense {
            walls.swap(0, 1);
        }
        if !walls[0].same_sense || walls[1].same_sense {
            return Err(unsupported("annular walls need opposite material senses"));
        }
        let origin = walls[0].cylinder.origin.clone();
        for wall in &mut walls {
            if sign(
                &wall.cylinder.axis.cross(&definition.pull).norm_squared(),
                dec,
            )? != Sign::Zero
                || sign(
                    &wall
                        .cylinder
                        .origin
                        .minus(&origin)
                        .cross(&definition.pull)
                        .norm_squared(),
                    dec,
                )? != Sign::Zero
            {
                return Err(unsupported(
                    "annular walls must be coaxial and parallel to pull",
                ));
            }
            // An axial origin shift does not change a cylinder. Use one frame
            // so all endpoint and separation proofs have identical heights.
            wall.cylinder.origin = origin.clone();
        }
        if sign(&walls[1].cylinder.radius_sq, dec)? != Sign::Positive
            || sign(
                &(&walls[0].cylinder.radius_sq - &walls[1].cylinder.radius_sq),
                dec,
            )? != Sign::Positive
        {
            return Err(unsupported(
                "source annular radii must be positive and separated",
            ));
        }
        for (plane, _) in &caps {
            if sign(&plane.normal.cross(&definition.pull).norm_squared(), dec)? != Sign::Zero {
                return Err(unsupported("annular caps must be perpendicular to pull"));
            }
        }
        let mut build = |wall: AnnularWall| -> Result<RoundCell, Refusal> {
            // An unchanged wall must not invent a neutral-curve requirement.
            let fixed = DraftDefinition {
                pull: definition.pull.clone(),
                neutral: Plane {
                    origin: origin.clone(),
                    normal: definition.pull.clone(),
                    ref_dir: wall.cylinder.ref_dir.clone(),
                },
                sine: Real::zero(),
                cosine: Real::one(),
            };
            let mut law = if wall.selected {
                definition.clone()
            } else {
                fixed.clone()
            };
            if !wall.same_sense {
                // The outward material normal at the hole points toward the
                // axis: a positive authored angle increases its radius.
                law.sine = -law.sine;
            }
            law.support(&Surface::Cylinder(wall.cylinder.clone()), true, dec)?;
            let original = RoundCell::new(
                &wall.cylinder,
                wall.origins.clone(),
                caps.clone(),
                &fixed,
                dec,
            )?;
            let neutral = reduced(
                law.neutral.origin.minus(&origin).dot(&law.neutral.normal)
                    / law.pull.dot(&law.neutral.normal),
                dec,
            )?;
            let tangent = reduced(&law.sine / &law.cosine, dec)?;
            let radius = reduced(wall.cylinder.radius_sq.sqrt(), dec)?;
            for height in &original.heights {
                let r = reduced(&radius - &(&tangent * &(height - &neutral)), dec)?;
                if sign(&r, dec)? != Sign::Positive {
                    return Err(unsupported(
                        "annular wall reaches or crosses the axis within its finite height",
                    ));
                }
            }
            RoundCell::new(&wall.cylinder, wall.origins, caps.clone(), &law, dec)
        };
        let [outer, inner] = walls;
        let outer = build(outer)?;
        let inner = build(inner)?;
        for i in 0..2 {
            if sign(&(&outer.radii[i] - &inner.radii[i]), dec)? != Sign::Positive {
                return Err(unsupported(
                    "annular walls touch or cross within the finite height",
                ));
            }
        }
        // Both radii and their difference are affine in height. Strict endpoint
        // inequalities certify the complete interval, not sampled sections.
        Ok(Self { outer, inner })
    }

    /// Exact strict material membership, excluding the hole and every boundary.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if !self.outer.contains_interior(point, dec)? {
            return Ok(false);
        }
        let delta = point.minus(&self.inner.origin);
        let height = delta.dot(&self.inner.axis);
        let t =
            (&height - &self.inner.heights[0]) / (&self.inner.heights[1] - &self.inner.heights[0]);
        let radius = &self.inner.radii[0] + &(t * (&self.inner.radii[1] - &self.inner.radii[0]));
        Ok(sign(
            &(delta.cross(&self.inner.axis).norm_squared() - radius.square()),
            dec,
        )? == Sign::Positive)
    }

    /// Exact positive-area anchoring, including an annular cap section.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        // Removing a strictly smaller coaxial radial graph cannot remove all
        // relatively open area of any plane section of this outer radial graph.
        self.outer.anchored(neutral, dec)
    }

    /// A strict material point halfway along the height and radial interval.
    pub fn interior_witness(&self) -> Point3 {
        let height = (&self.outer.heights[0] + &self.outer.heights[1]) / Real::from_i64(2);
        let radius = self
            .outer
            .radii
            .iter()
            .chain(&self.inner.radii)
            .fold(Real::zero(), |sum, r| &sum + r)
            / Real::from_i64(4);
        self.outer
            .origin
            .plus(&self.outer.axis.scaled(&height))
            .plus(
                &self
                    .outer
                    .reference
                    .scaled(&(&radius / &self.outer.reference.norm_squared().sqrt())),
            )
    }
}
