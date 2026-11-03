//! Finite oblique cylindrical draft cells with two unchanged transverse caps.

use crate::draft::{DraftDefinition, reduced, sign};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, ObliqueDraft, Plane, Point3, SplineSupport, Surface};

/// One embedded ruled side and its complete planar cap disks.
#[derive(Clone, Debug)]
pub struct ObliqueCell {
    /// Checked exact ruled carrier, retaining its finite branch.
    pub support: ObliqueDraft,
    /// Source origins, unchanged by draft.
    pub origins: Vec<Origin>,
    /// Complete cap definitions and origins, ordered along pull.
    pub caps: [(Plane, Vec<Origin>); 2],
}
impl ObliqueCell {
    /// Rebuilds a complete circular source cell as an exact oblique draft.
    /// Partial caps, nontransverse neighbours and folds retain typed refusals.
    pub fn new(
        cylinder: &Cylinder,
        origins: Vec<Origin>,
        caps: [(Plane, Vec<Origin>); 2],
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut ordered = [None, None];
        let mut heights = [Real::zero(), Real::zero()];
        for (plane, names) in caps {
            if sign(&plane.normal.cross(&definition.pull).norm_squared(), dec)? != Sign::Zero {
                return Err(Refusal::new(
                    "draft.unsupported-oblique-cap",
                    "oblique cylindrical draft needs complete caps perpendicular to pull",
                ));
            }
            let rate = plane.normal.dot(&definition.pull);
            let i = match sign(&rate, dec)? {
                Sign::Negative => 0,
                Sign::Positive => 1,
                Sign::Zero => {
                    return Err(Refusal::new("draft.nonclosing-boundary", "zero cap normal"));
                }
            };
            if ordered[i].is_some() {
                return Err(Refusal::new(
                    "draft.nonclosing-boundary",
                    "oblique draft has no opposing cap pair",
                ));
            }
            heights[i] = reduced(
                &plane.origin.minus(&cylinder.origin).dot(&plane.normal) / &rate,
                dec,
            )?;
            ordered[i] = Some((plane, names));
        }
        let [Some(lo), Some(hi)] = ordered else {
            return Err(Refusal::new(
                "draft.nonclosing-boundary",
                "oblique draft is missing a complete cap",
            ));
        };
        let support = ObliqueDraft::new(
            cylinder,
            definition.pull.clone(),
            &definition.neutral,
            reduced(&definition.sine / &definition.cosine, dec)?,
            heights,
            dec,
        )?;
        Ok(Self {
            support,
            origins,
            caps: [lo, hi],
        })
    }
    /// Exact closed membership, with strict interior distinguished from boundary.
    pub fn locate(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let z = self.support.coordinates(p)[2].clone();
        let h = self.support.heights();
        let a = sign(&(&z - &h[0]), dec)?;
        let b = sign(&(&h[1] - &z), dec)?;
        if a == Sign::Negative || b == Sign::Negative {
            return Ok(Sign::Positive);
        }
        let side = self.support.side_at(p, dec)?;
        Ok(if side == Sign::Positive {
            side
        } else if a == Sign::Zero || b == Sign::Zero || side == Sign::Zero {
            Sign::Zero
        } else {
            Sign::Negative
        })
    }
    /// Exact strict membership for retained-region picks.
    pub fn contains_interior(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.locate(p, dec)? == Sign::Negative)
    }
    /// Proves a positive-area neutral section, or proves complete separation.
    /// An unresolved remote cut refuses rather than guessing retention.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let rate = self.support.axis().dot(&neutral.normal);
        let h = &neutral
            .origin
            .minus(self.support.origin())
            .dot(&neutral.normal)
            / &rate;
        let ends = self.support.heights();
        if sign(&(&h - &ends[0]), dec)? != Sign::Negative
            && sign(&(&ends[1] - &h), dec)? != Sign::Negative
        {
            return Ok(true);
        }
        let mut below = false;
        let mut above = false;
        // Positive rational weights contain the entire side and its cap disks
        // in the carrier-control hull. A separating plane proves no anchor.
        for p in self.support.carrier().controls() {
            match sign(&p.minus(&neutral.origin).dot(&neutral.normal), dec)? {
                Sign::Negative => below = true,
                Sign::Positive => above = true,
                Sign::Zero => (),
            }
        }
        if !below || !above {
            return Ok(false);
        }
        Err(Refusal::new(
            "draft.unsupported-oblique-retention",
            "remote neutral section needs an exact plane/ruled-cell extremum proof; use explicit retained-region points",
        ))
    }
    /// A model-origin radius enclosing the complete positive-weight carrier hull.
    pub fn bound_radius(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        self.support
            .carrier()
            .controls()
            .iter()
            .try_fold(Rational::zero(), |best, p| {
                let mut sum = Rational::zero();
                for x in [&p.x, &p.y, &p.z] {
                    let b = dec.enclose_at(x, 128)?;
                    sum = sum + b.lo().abs().max(b.hi().abs());
                }
                Ok(best.max(sum))
            })
    }
    /// Exact rigid transport of the finite law and both neighbouring caps.
    pub fn transform(&mut self, map: &catallax_shapes::Isometry) {
        let Surface::BSpline(SplineSupport::ObliqueDraft(support)) = map.surface(
            &Surface::BSpline(SplineSupport::ObliqueDraft(self.support.clone())),
        ) else {
            unreachable!("rigid transport preserves the checked chart family");
        };
        self.support = support;
        for (p, _) in &mut self.caps {
            p.origin = map.apply_point(&p.origin);
            p.normal = map.apply_vector(&p.normal);
            p.ref_dir = map.apply_vector(&p.ref_dir);
        }
    }
    /// An exact point strictly inside the finite cell, on its axis.
    pub fn interior_witness(&self) -> Point3 {
        let h = (&self.support.heights()[0] + &self.support.heights()[1]) / Real::from_i64(2);
        self.support.origin().plus(&self.support.axis().scaled(&h))
    }
    /// Complete outward support and cap surfaces, paired with their source names.
    pub fn surfaces(&self) -> [(catallax_shapes::Surface, Vec<Origin>); 3] {
        [
            (
                catallax_shapes::Surface::BSpline(catallax_shapes::SplineSupport::ObliqueDraft(
                    self.support.clone(),
                )),
                self.origins.clone(),
            ),
            (
                catallax_shapes::Surface::Plane(self.caps[0].0.clone()),
                self.caps[0].1.clone(),
            ),
            (
                catallax_shapes::Surface::Plane(self.caps[1].0.clone()),
                self.caps[1].1.clone(),
            ),
        ]
    }
}

/// Whether this selected cylindrical law requires the oblique ruled carrier.
pub fn needed(definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&definition.sine, dec)? != Sign::Zero
        && sign(
            &definition
                .neutral
                .normal
                .cross(&definition.pull)
                .norm_squared(),
            dec,
        )? != Sign::Zero)
}
