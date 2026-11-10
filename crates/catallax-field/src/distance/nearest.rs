//! Proved nearest named face ("ownership") from any point.
//!
//! Every answer comes from the certified trimmed-face distance bounds of
//! [`Prepared::point_to_face`]. A named face owns a point only when its upper
//! bound is strictly below every other name's lower bound. Otherwise the query
//! returns [`Ownership::Refused`] with the unseparated names and their bounds.
//! No floating-point preview value, clock or worker count selects an answer.

use super::Prepared;
use crate::{Direction, FaceId, geometry};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    Point3, Surface,
    distance::{self, Membership, Options},
};

/// Controls for [`Prepared::nearest_face`].
///
/// Refinement runs in fixed stages of 2, 16, 128, ... cells per face, ending at
/// `max_cells`. The first stage that separates one name ends the query.
#[derive(Clone, Debug)]
pub struct NearestOptions {
    /// Requested enclosure width for each name's distance, in millimetres.
    /// It ends refinement early; it never widens a proved separation.
    pub width: Rational,
    /// Largest cell count per registry face at the last stage.
    pub max_cells: usize,
}

impl Default for NearestOptions {
    fn default() -> Self {
        Self {
            width: Rational::from_ratio(1, 1_000_000).unwrap(),
            max_cells: 1024,
        }
    }
}

/// Certified distance enclosure from the query point to one named face.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    /// Stable recipe name.
    pub origin: Origin,
    /// Surviving registry entries carrying this name, in registry order.
    pub faces: Vec<FaceId>,
    /// Proven lower bound on the distance to the union of `faces`, in millimetres.
    pub lower: Rational,
    /// Proven upper bound from an admitted witness, absent when none was found.
    pub upper: Option<Rational>,
}

/// A proved owner: the one named face strictly nearer than every other name.
#[derive(Clone, Debug)]
pub struct Owner {
    /// Stable recipe name of the nearest face.
    pub origin: Origin,
    /// Registry entry whose witness attains `upper`.
    pub face: FaceId,
    /// Proven lower bound on the distance to the owner, in millimetres.
    pub lower: Rational,
    /// Proven upper bound on the distance to the owner, in millimetres.
    pub upper: Rational,
    /// Exact admitted point on the owner's trimmed face at distance at most `upper`.
    pub witness: Point3,
    /// Unit vector from the query point toward `witness`, enclosed component by
    /// component. `Undefined` when the witness is the query point itself.
    /// It encloses the vector to the admitted witness only. It carries no
    /// uniqueness or angular-error guarantee: at the centre of a spherical
    /// face, opposite witnesses are equally near.
    pub direction: Direction,
    /// The excluded name with the smallest proven lower bound, with the bounds
    /// that excluded it. It limits `margin`. It is not a certified
    /// second-nearest name: another excluded name with a larger lower bound
    /// can be nearer, because exclusion stops refining each name as soon as
    /// its lower bound exceeds the owner's upper bound. Absent when the
    /// result has a single name.
    pub limiting: Option<Candidate>,
    /// Proven lower bound on the distance from the query point to the border
    /// between the owner's region and any other name's region, in millimetres:
    /// half of `limiting.lower - upper`. Every point closer than this has the
    /// same owner. Absent when the result has a single name.
    pub margin: Option<Rational>,
    /// Cells visited over all stages and faces.
    pub cells: usize,
}

/// Names whose distance bounds still overlap after the last stage.
#[derive(Clone, Debug, PartialEq)]
pub struct Ambiguity {
    /// Unseparated names in registry order. At least two, or one without any
    /// admitted witness.
    pub contenders: Vec<Candidate>,
    /// Proven bound on the difference between any two contenders' distances,
    /// in millimetres: the largest upper bound minus the smallest lower bound.
    /// Absent when some contender has no admitted witness. A point on an edge
    /// shared by two names, or in the region whose nearest point is that edge,
    /// has a spread that shrinks toward zero as both enclosures converge.
    pub spread: Option<Rational>,
    /// Cells visited over all stages and faces.
    pub cells: usize,
}

impl Ambiguity {
    /// The stable refusal for this answer.
    pub fn refusal(&self) -> Refusal {
        let names: Vec<String> = self
            .contenders
            .iter()
            .map(|c| format!("{:?}", c.origin))
            .collect();
        Refusal::new(
            "field.ownership-ambiguous",
            format!(
                "no single nearest named face: distance bounds of {} overlap after {} cells",
                names.join(", "),
                self.cells
            ),
        )
    }
}

/// Answer of [`Prepared::nearest_face`]. Both variants are deterministic.
#[derive(Clone, Debug)]
pub enum Ownership {
    /// Proved owner from certified distance bounds.
    Proved(Box<Owner>),
    /// Typed refusal: equidistant names, or names not separated within the
    /// requested cells. Never a guess.
    Refused(Ambiguity),
}

struct Entry {
    candidate: Candidate,
    witness: Option<(FaceId, Point3)>,
    live: bool,
}

impl Prepared<'_> {
    /// Prepares the exact planar trims of every surviving plane face. Ownership
    /// queries call it first, so their bounds never depend on earlier queries.
    /// Idempotent. Call it once before sharing this object between workers, so
    /// that each query's logical work is the same under any partition.
    pub fn prepare_nearest(&self, dec: &mut Decider) -> Result<(), Refusal> {
        for (index, domain) in self.domains.iter().enumerate() {
            if domain
                .as_ref()
                .is_some_and(|d| matches!(d.surface(), Surface::Plane(_)))
            {
                self.planar_boundary(index, dec)?;
            }
        }
        Ok(())
    }

    /// Proves which named face is nearest to `point`, with its distance bounds,
    /// a witness direction and the margin to the nearest ownership border.
    /// Every surviving named face takes part. Equidistant names, such as at a
    /// point on a shared edge, return [`Ownership::Refused`].
    pub fn nearest_face(
        &self,
        point: &Point3,
        options: &NearestOptions,
        dec: &mut Decider,
    ) -> Result<Ownership, Refusal> {
        if options.width <= Rational::zero() || options.max_cells < 2 {
            return Err(Refusal::new(
                "field.ownership-options",
                "positive width and at least two cells required",
            ));
        }
        self.prepare_nearest(dec)?;
        let mut entries: Vec<Entry> = Vec::new();
        for (index, domain) in self.domains.iter().enumerate() {
            if domain.is_none() {
                continue;
            }
            let origin = &self.plan.faces[index].origin;
            if let Some(entry) = entries.iter_mut().find(|e| &e.candidate.origin == origin) {
                entry.candidate.faces.push(FaceId(index));
            } else {
                entries.push(Entry {
                    candidate: Candidate {
                        origin: origin.clone(),
                        faces: vec![FaceId(index)],
                        lower: Rational::zero(),
                        upper: None,
                    },
                    witness: None,
                    live: true,
                });
            }
        }
        if entries.is_empty() {
            return Err(Refusal::new(
                "field.ownership-empty",
                "the result has no surviving named face",
            ));
        }
        let mut cells = 0usize;
        let mut stage = 2usize;
        loop {
            let controls = Options {
                width: options.width.clone(),
                max_cells: stage,
            };
            for entry in entries.iter_mut().filter(|e| e.live) {
                let mut lower: Option<Rational> = None;
                let mut upper: Option<Rational> = None;
                let mut witness = None;
                for &face in &entry.candidate.faces {
                    let mut bounds = self.point_to_face(face, point, &controls, dec)?;
                    cells += bounds.cells;
                    if bounds.upper.is_none() {
                        let (source, used) = self.source_lower(face, point, &controls, dec)?;
                        cells += used;
                        bounds.lower = bounds.lower.max(source);
                    }
                    lower = Some(lower.map_or(bounds.lower.clone(), |l| l.min(bounds.lower)));
                    if let Some(u) = bounds.upper
                        && upper.as_ref().is_none_or(|best| u < *best)
                    {
                        witness = bounds.witnesses.map(|[_, w]| (face, w));
                        upper = Some(u);
                    }
                }
                // Bounds from a later stage are recomputed from scratch; keep the
                // tighter proven side of each, which is still a proof.
                let lower = lower.expect("every name has a face");
                if lower > entry.candidate.lower {
                    entry.candidate.lower = lower;
                }
                if let Some(u) = upper
                    && entry.candidate.upper.as_ref().is_none_or(|best| u < *best)
                {
                    entry.candidate.upper = Some(u);
                    entry.witness = witness;
                }
            }
            let best = entries
                .iter()
                .filter(|e| e.live)
                .filter_map(|e| e.candidate.upper.clone())
                .min();
            if let Some(best) = &best {
                for entry in entries.iter_mut().filter(|e| e.live) {
                    if entry.candidate.lower > *best {
                        entry.live = false;
                    }
                }
            }
            let live = entries.iter().filter(|e| e.live).count();
            if live == 1 && best.is_some() {
                return Ok(Ownership::Proved(Box::new(owner(
                    point, &entries, cells, dec,
                )?)));
            }
            // Converged enclosures do not refine further, so a later stage
            // cannot separate them.
            let settled = entries.iter().filter(|e| e.live).all(|e| {
                e.candidate
                    .upper
                    .as_ref()
                    .is_some_and(|u| u - &e.candidate.lower <= options.width)
            });
            if stage >= options.max_cells || settled {
                let contenders: Vec<Candidate> = entries
                    .into_iter()
                    .filter(|e| e.live)
                    .map(|e| e.candidate)
                    .collect();
                let lowest = contenders.iter().map(|c| &c.lower).min().cloned();
                let highest = contenders
                    .iter()
                    .map(|c| c.upper.clone())
                    .collect::<Option<Vec<_>>>()
                    .and_then(|uppers| uppers.into_iter().max());
                return Ok(Ownership::Refused(Ambiguity {
                    spread: highest.zip(lowest).map(|(u, l)| u - l),
                    contenders,
                    cells,
                }));
            }
            stage = stage.saturating_mul(8).min(options.max_cells);
        }
    }
}

impl Prepared<'_> {
    // A face's surviving region lies inside its own primitive's face, so the
    // distance to that source face is a lower bound. Result trims can keep a
    // buried or tangent-touched face's lower bound at its carrier distance;
    // this bound uses only the face's own primitive constraints, so cells far
    // from its source face are excluded. It never supplies an upper bound.
    fn source_lower(
        &self,
        id: FaceId,
        point: &Point3,
        controls: &Options,
        dec: &mut Decider,
    ) -> Result<(Rational, usize), Refusal> {
        let face = &self.plan.faces[id.0];
        let leaf = &self.plan.leaves[face.leaf];
        let Some(domain) = self.domains[id.0].as_ref() else {
            return Ok((Rational::zero(), 0));
        };
        if !matches!(leaf, geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_)) {
            return Ok((Rational::zero(), 0));
        }
        let coincident = &self.coincident[id.0][face.leaf];
        let mut trim = |_: &Point3, cover: &[Ball; 3], dec: &mut Decider| {
            for (c, &same) in leaf.constraints().iter().zip(coincident) {
                if !same && c.polynomial.eval_ball(cover, dec)?.lo() > Rational::zero() {
                    return Ok(Membership {
                        witness: false,
                        excluded: true,
                    });
                }
            }
            Ok(Membership {
                witness: false,
                excluded: false,
            })
        };
        let bounds = distance::point_distance(domain, point, controls, &mut trim, dec)?;
        Ok((bounds.lower, bounds.cells))
    }
}

fn owner(
    point: &Point3,
    entries: &[Entry],
    cells: usize,
    dec: &mut Decider,
) -> Result<Owner, Refusal> {
    let winner = entries.iter().find(|e| e.live).expect("one live name");
    let (face, witness) = winner.witness.clone().expect("a live owner has a witness");
    let upper = winner.candidate.upper.clone().expect("owner upper");
    // The smallest excluding lower bound limits the margin. Fixed tie order:
    // the earliest registry name among equal lower bounds.
    let limiting = entries
        .iter()
        .filter(|e| !e.live)
        .min_by(|a, b| a.candidate.lower.cmp(&b.candidate.lower))
        .map(|e| e.candidate.clone());
    let margin = limiting.as_ref().map(|r| (&r.lower - &upper).mul_pow2(-1));
    let offset = witness.minus(point);
    let direction = if geometry::sign(&offset.norm_squared(), dec)? == Sign::Zero {
        Direction::Undefined
    } else {
        crate::unit_direction(offset, dec)?
    };
    Ok(Owner {
        origin: winner.candidate.origin.clone(),
        face,
        lower: winner.candidate.lower.clone(),
        upper,
        witness,
        direction,
        limiting,
        margin,
        cells,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_recipe::{Role, StepId};
    use catallax_shapes::RPoint3;

    fn entry(step: &str, lower: i64, upper: i64, live: bool) -> Entry {
        let origin = Origin {
            step: StepId::new(step).unwrap(),
            role: Role::Start,
        };
        let witness = RPoint3::new(
            Rational::from_i64(upper),
            Rational::zero(),
            Rational::zero(),
        );
        Entry {
            candidate: Candidate {
                origin,
                faces: vec![FaceId(0)],
                lower: Rational::from_i64(lower),
                upper: Some(Rational::from_i64(upper)),
            },
            witness: Some((FaceId(0), witness.to_real())),
            live,
        }
    }

    // Owner [1, 1]; A excluded with [2, 100] though its distance may be 50;
    // B excluded with [3, 3]. B is nearer than A can be shown to be, yet A's
    // lower bound is the one that limits the margin. The field names A as the
    // limiting name and never claims it is second nearest.
    #[test]
    fn limiting_name_is_the_smallest_excluding_lower_bound() {
        let entries = [
            entry("owner", 1, 1, true),
            entry("a", 2, 100, false),
            entry("b", 3, 3, false),
        ];
        let point = RPoint3::new(Rational::zero(), Rational::zero(), Rational::zero()).to_real();
        let owner = owner(&point, &entries, 0, &mut Decider::default()).unwrap();
        assert_eq!(owner.origin.step.as_str(), "owner");
        assert_eq!(owner.limiting.unwrap().origin.step.as_str(), "a");
        assert_eq!(owner.margin, Some(Rational::from_ratio(1, 2).unwrap()));
    }
}
