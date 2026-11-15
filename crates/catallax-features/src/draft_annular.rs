//! Finite hollow axial cells and checked oblique radial-separation bounds.

use super::oblique::ObliqueCell;
use crate::draft::{DraftDefinition, reduced, round::RoundCell, sign};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{BSplineCurve, Cylinder, Plane, Point3, SplineBoundary, Surface};

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

/// Exact relation of the two radial laws over the complete finite height.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnularSeparation {
    /// Every outer radius strictly exceeds the inner radius.
    Separated,
    /// The minimum radial gap is zero, with no negative gap.
    Touching,
    /// The gap takes both positive and negative values.
    Crossing,
    /// The outer radius never exceeds the inner radius.
    Reversed,
}

/// Checked radial separation and explicit plane/wall intersections.
/// This is a bound on two finite supports, never a publishable body or a
/// shortest Euclidean wall-distance certificate. Cap topology remains unbuilt.
#[derive(Clone, Debug)]
pub struct ObliqueAnnularBound {
    walls: [ObliqueCell; 2],
    caps: [(Plane, [BSplineCurve; 2]); 2],
    radial_gap: [Rational; 2],
    separation: AnnularSeparation,
    source_detail: String,
}

impl ObliqueAnnularBound {
    /// Checks coaxial source walls, finite ruled laws and all four cap curves.
    /// Zero neutral-axis rate and unsupported folds retain typed refusals.
    pub fn new(
        walls: [AnnularWall; 2],
        caps: [(Plane, Vec<Origin>); 2],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut walls = checked_source(walls, &caps, definition, dec)?;
        let axis_rate = definition.neutral.normal.dot(&definition.pull);
        let zero_rate = sign(&axis_rate, dec)? == Sign::Zero;
        // A cylinder has no distinguished azimuth. Use one orthonormal frame
        // for the two polar laws even if the source reference directions differ.
        walls[1].cylinder.ref_dir = walls[0].cylinder.ref_dir.clone();
        let tangent = reduced(&definition.sine / &definition.cosine, dec)?;
        let span = caps[1]
            .0
            .origin
            .minus(&caps[0].0.origin)
            .dot(&definition.pull)
            .square()
            .sqrt();
        let mut enclose = |value: &Real| -> Result<String, Refusal> {
            let b = dec.enclose_budgeted_at(value, 128)?;
            Ok(format!("[{},{}]", b.lo(), b.hi()))
        };
        let tilt = if zero_rate {
            "undefined (neutral_axis_rate=0)".to_owned()
        } else {
            enclose(
                &(definition
                    .neutral
                    .normal
                    .cross(&definition.pull)
                    .norm_squared()
                    / axis_rate.square()),
            )?
        };
        let source_detail = format!(
            "outer_radius_squared_mm2={}; inner_radius_squared_mm2={}; selected_outer={}; selected_inner={}; tan_angle={}; cap_span_mm={}; neutral_transverse_ratio_squared={}",
            enclose(&walls[0].cylinder.radius_sq)?,
            enclose(&walls[1].cylinder.radius_sq)?,
            walls[0].selected,
            walls[1].selected,
            enclose(&tangent)?,
            enclose(&span)?,
            tilt,
        );
        if zero_rate {
            return Err(Refusal::new(
                "draft.oblique-annulus-zero-neutral-rate",
                format!(
                    "{source_detail}; no closed transverse neutral curve; certified_radial_gap=undefined; body_published=false"
                ),
            ));
        }
        let mut build = |wall: AnnularWall| {
            let mut law = definition.clone();
            if !wall.selected {
                law.sine = Real::zero();
                law.cosine = Real::one();
            } else if !wall.same_sense {
                law.sine = -law.sine;
            }
            ObliqueCell::new(&wall.cylinder, wall.origins, caps.clone(), &law, dec)
        };
        let [outer, inner] = walls;
        let walls = [build(outer)?, build(inner)?];
        let [outer, inner] = walls.each_ref().map(|wall| &wall.support);
        let differences = outer
            .heights()
            .each_ref()
            .map(|z| outer.amplitude(z) - inner.amplitude(z));
        let lean = reduced(
            ((&outer.lean()[0] - &inner.lean()[0]).square()
                + (&outer.lean()[1] - &inner.lean()[1]).square())
            .sqrt(),
            dec,
        )?;
        let order = sign(&(&differences[1] - &differences[0]), dec)?;
        let (lo, hi) = if order == Sign::Negative {
            (1, 0)
        } else {
            (0, 1)
        };
        // The gap is affine in z and equals dA(z)+db*cos(u)+dc*sin(u).
        // Its extrema on the entire closed domain are the two endpoint
        // amplitudes plus/minus sqrt(db^2+dc^2), attained by unit rays.
        let minimum = reduced(&differences[lo] - &lean, dec)?;
        let maximum = reduced(&differences[hi] + &lean, dec)?;
        let separation = match sign(&minimum, dec)? {
            Sign::Positive => AnnularSeparation::Separated,
            Sign::Zero => AnnularSeparation::Touching,
            Sign::Negative if sign(&maximum, dec)? == Sign::Positive => AnnularSeparation::Crossing,
            Sign::Negative => AnnularSeparation::Reversed,
        };
        let radial_gap = [
            dec.enclose_budgeted_at(&minimum, 128)?.lo(),
            dec.enclose_budgeted_at(&maximum, 128)?.hi(),
        ];
        let caps = std::array::from_fn(|i| {
            let boundary = [SplineBoundary::VMin, SplineBoundary::VMax][i];
            (
                walls[0].caps[i].0.clone(),
                std::array::from_fn(|j| {
                    BSplineCurve::boundary(walls[j].support.carrier(), boundary)
                }),
            )
        });
        for (plane, rings) in &caps {
            for ring in rings {
                if ring.is_closed() != Some(true) {
                    return Err(unsupported("oblique cap intersection is not closed"));
                }
                for point in ring.controls() {
                    if sign(&point.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero {
                        return Err(unsupported("oblique cap intersection leaves its plane"));
                    }
                }
            }
        }
        Ok(Self {
            walls,
            caps,
            radial_gap,
            separation,
            source_detail,
        })
    }

    /// Complete finite outer and inner laws in one axial and angular frame.
    pub fn walls(&self) -> &[ObliqueCell; 2] {
        &self.walls
    }

    /// Lower and upper cap planes, each with outer and inner exact intersections.
    /// These complete carrier rows are curves, not assembled trimmed cap faces.
    pub fn cap_intersections(&self) -> &[(Plane, [BSplineCurve; 2]); 2] {
        &self.caps
    }

    /// Outward enclosure of every signed radial gap, in millimetres.
    pub fn radial_gap(&self) -> &[Rational; 2] {
        &self.radial_gap
    }

    /// Relation decided from exact extrema, independently of enclosure width.
    pub fn separation(&self) -> AnnularSeparation {
        self.separation
    }

    fn publication_refusal(&self) -> Refusal {
        let code = match self.separation {
            AnnularSeparation::Separated => "draft.unsupported-oblique-annulus",
            AnnularSeparation::Touching => "draft.oblique-annulus-touching",
            AnnularSeparation::Crossing => "draft.oblique-annulus-crossing",
            AnnularSeparation::Reversed => "draft.oblique-annulus-reversed",
        };
        Refusal::new(
            code,
            format!(
                "{}; separation={:?}; certified_radial_gap_mm=[{},{}]; cap_intersections=4 exact closed planar carrier rows; shared_trimmed_cap_rings=unbuilt; body_published=false",
                self.source_detail, self.separation, self.radial_gap[0], self.radial_gap[1],
            ),
        )
    }
}

fn checked_source(
    mut walls: [AnnularWall; 2],
    caps: &[(Plane, Vec<Origin>); 2],
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<[AnnularWall; 2], Refusal> {
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
    for (plane, _) in caps {
        if sign(&plane.normal.cross(&definition.pull).norm_squared(), dec)? != Sign::Zero {
            return Err(unsupported("annular caps must be perpendicular to pull"));
        }
    }
    Ok(walls)
}

impl AnnularCell {
    /// Certifies the bounded annular slice before publishing any topology.
    /// Folds, contacts, crossings and swallowed ends need a curved arrangement.
    pub fn new(
        walls: [AnnularWall; 2],
        caps: [(Plane, Vec<Origin>); 2],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let walls = checked_source(walls, &caps, definition, dec)?;
        let origin = walls[0].cylinder.origin.clone();
        if walls.iter().any(|wall| wall.selected) && super::oblique::needed(definition, dec)? {
            let bound = ObliqueAnnularBound::new(walls, caps, definition, dec)?;
            return Err(bound.publication_refusal());
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
