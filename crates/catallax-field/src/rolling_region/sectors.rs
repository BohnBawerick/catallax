//! Exact angular routing for an explicitly supplied set of finite mouths.
//! This proves a partition of directions, not incidence with a recipe solid.

use super::{Region, Station, sign};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::Vector3;

fn refusal(detail: &str) -> Refusal {
    Refusal::new("field.rolling-sectors-unresolved", detail)
}

/// End of a finite mouth chart in its own oriented frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Transition {
    /// Normalized station -1.
    Start,
    /// Normalized station +1.
    End,
}

/// Disjoint angular partition, including both closed-sector endpoints.
/// These labels are routing information, NOT solid membership or named faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sector {
    /// Strictly outside every supplied closed mouth sector.
    Torus,
    /// Strictly inside this supplied mouth sector (input index).
    Mouth(usize),
    /// Exactly on this supplied mouth's transition ray (input index).
    Transition(usize, Transition),
}

/// Immutable, pairwise strictly separated mouth sectors on one oriented ring.
///
/// All definitions must share the exact origin, axis, body and ball radius.
/// Mouth radii, offsets and normal cones may differ. Closed-sector touching,
/// overlap, unresolved comparisons and incompatible rings refuse. No recipe
/// discovery, material-equivalence or boundary-distance proof is implied.
#[derive(Clone, Debug)]
pub struct Sectors {
    regions: Vec<Region>,
}

impl Sectors {
    /// Certifies the entire angular domains, without sampled stations.
    pub fn new(regions: Vec<Region>, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let Some(first) = regions.first() else {
            return Err(refusal("at least one finite mouth is required"));
        };
        let ring = first.definition();
        for (i, region) in regions.iter().enumerate() {
            dec.charge_work(1)?;
            let d = region.definition();
            for difference in [
                d.origin().minus(ring.origin()).norm_squared(),
                d.axis().minus(ring.axis()).norm_squared(),
                d.body() - ring.body(),
                d.radius() - ring.radius(),
            ] {
                if sign(&difference, dec)? != Sign::Zero {
                    return Err(refusal("finite mouths must share one oriented ring"));
                }
            }
            for previous in &regions[..i] {
                dec.charge_work(1)?;
                // Each half-angle is strictly between 0 and pi/2 (Region::new).
                // Their sum is in (0,pi), where cosine strictly decreases.
                // The minor angle between unit mouth directions exceeds that
                // sum iff dot < cos(a+b). This excludes CLOSED-sector contact,
                // for every direction in either arc, including across the seam.
                let limit = &(&region.cosine * &previous.cosine) - &(&region.sine * &previous.sine);
                let dot = d.toward().dot(previous.definition().toward());
                if sign(&(&limit - &dot), dec)? != Sign::Positive {
                    return Err(refusal("closed finite mouth sectors touch or overlap"));
                }
            }
        }
        dec.check_deadline()?;
        Ok(Self { regions })
    }

    /// The validated definitions in original order; no mutable access.
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Routes a nonzero radial vector in the common plane. Scale is irrelevant.
    /// A refusal is retained if any exact comparison cannot be settled.
    pub fn locate(&self, direction: &Vector3, dec: &mut Decider) -> Result<Sector, Refusal> {
        dec.check_deadline()?;
        if sign(&direction.norm_squared(), dec)? != Sign::Positive
            || sign(&direction.dot(self.regions[0].definition().axis()), dec)? != Sign::Zero
        {
            return Err(refusal("direction must be nonzero and in the ring plane"));
        }
        for (i, region) in self.regions.iter().enumerate() {
            dec.charge_work(1)?;
            let d = region.definition();
            let x = direction.dot(d.toward());
            if sign(&x, dec)? != Sign::Positive {
                continue;
            }
            let z = direction.dot(d.across());
            let start = sign(&(&(&x * &region.sine) + &(&z * &region.cosine)), dec)?;
            let end = sign(&(&(&x * &region.sine) - &(&z * &region.cosine)), dec)?;
            if start == Sign::Negative || end == Sign::Negative {
                continue;
            }
            return Ok(if start == Sign::Zero {
                Sector::Transition(i, Transition::Start)
            } else if end == Sign::Zero {
                Sector::Transition(i, Transition::End)
            } else {
                Sector::Mouth(i)
            });
        }
        Ok(Sector::Torus)
    }

    /// Exact common rolling-ball contacts at a mouth/ordinary-torus transition.
    /// At either endpoint D=m, hence H=r: both centre and contacts agree with
    /// the ordinary torus. This is endpoint agreement, not a whole-seam proof
    /// for the lower envelope or an exposed boundary certificate.
    pub fn transition(
        &self,
        index: usize,
        end: Transition,
        dec: &mut Decider,
    ) -> Result<Station, Refusal> {
        dec.check_deadline()?;
        let region = self
            .regions
            .get(index)
            .ok_or_else(|| refusal("unknown mouth index"))?;
        let u = Rational::from_i64(match end {
            Transition::Start => -1,
            Transition::End => 1,
        });
        region.station(&u, dec)
    }
}
