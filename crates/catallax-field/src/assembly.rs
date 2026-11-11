//! Bounded assembly queries on placed, complete field definitions.
//!
//! Distance lower bounds cover the Cartesian product of both solids. Upper
//! bounds use certified points in the closed solids. No residual, preview or
//! nearest-owner vector is used as a distance certificate.
use crate::{Location, Plan, Workspace, lattice::BoxState, mass, printing};
use catallax_measure::Enclosure;
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{Pose, RPoint3};

type Box3 = [[Rational; 2]; 3];
type Point = [Rational; 3];

/// A field with an additional rigid placement (recipe placements apply first).
#[derive(Clone, Debug)]
pub struct PlacedPart {
    plan: Plan,
}
impl PlacedPart {
    /// Applies a checked rigid placement; all subsequent coordinates are world mm.
    pub fn new(plan: &Plan, placement: &Pose, dec: &mut Decider) -> Result<Self, Refusal> {
        Ok(Self {
            plan: plan.transformed(placement, dec)?,
        })
    }
}

/// Explicit resolution and deterministic resource cap for a new assembly query.
#[derive(Clone, Debug)]
pub struct Options {
    /// Maximum width of the returned smallest-gap enclosure, in mm. Must be positive.
    pub gap_width: Rational,
    /// Maximum processed box pairs; exhaustion never publishes a partial answer.
    pub max_pairs: usize,
}

/// Proven relation, or a typed reason the requested relation remains unresolved.
#[derive(Clone, Debug)]
pub enum Relation {
    /// A world point certified strictly inside both parts (positive-volume clash).
    Clash {
        /// Common interior point in world mm.
        witness: Point,
    },
    /// Disjoint closed solids and an enclosure of their global smallest distance.
    Clear {
        /// Strictly positive lower bound; width no greater than the request.
        gap: Enclosure,
    },
    /// Contact, unsupported definitions, empty operands or resource exhaustion.
    Unresolved(Refusal),
}

/// Compare two placed parts. Touching is unresolved, never labelled clear or clash.
pub fn relation(a: &PlacedPart, b: &PlacedPart, options: &Options, dec: &mut Decider) -> Relation {
    match solve(&a.plan, &b.plan, options, dec) {
        Ok(result) => result,
        Err(reason) => Relation::Unresolved(reason),
    }
}

/// Uniform material mass and centre of gravity of a placed field.
#[derive(Clone, Debug)]
pub struct Mass {
    /// Mass in grams, with density supplied in grams per cubic millimetre.
    pub grams: Enclosure,
    /// World mm; absent for exactly zero volume.
    pub centre_of_gravity: Option<[Enclosure; 3]>,
}

/// Certified uniform-density mass. Density must be strictly positive and exact.
/// `accuracy` retains the existing unit-density volume/centroid/inertia contract;
/// the mass width is density times the volume width. No material is inferred.
/// Existing integration and definition refusals are preserved without fallback.
pub fn uniform_mass(
    part: &PlacedPart,
    density: &Rational,
    accuracy: &mass::Accuracy,
    options: &mass::Options,
    dec: &mut Decider,
) -> Result<Mass, Refusal> {
    if density <= &Rational::zero() {
        return Err(Refusal::new(
            "assembly.density",
            "uniform density must be positive in g/mm^3",
        ));
    }
    let properties = mass::properties(&part.plan, accuracy, options, dec)?;
    Ok(Mass {
        grams: Enclosure {
            lo: &properties.volume.lo * density,
            hi: &properties.volume.hi * density,
        },
        centre_of_gravity: properties.centroid,
    })
}

fn middle(b: &Box3) -> Point {
    std::array::from_fn(|i| (&b[i][0] + &b[i][1]).mul_pow2(-1))
}
fn point(p: &Point) -> RPoint3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone())
}
fn location(plan: &Plan, p: &Point, dec: &mut Decider) -> Result<Location, Refusal> {
    Ok(plan
        .proven_checked(&point(p), &mut Workspace::default(), dec)?
        .location)
}
fn material(l: Location) -> bool {
    matches!(l, Location::Inside | Location::Boundary)
}
fn square_distance(a: &Point, b: &Point) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        let d = &a[i] - &b[i];
        &sum + &(&d * &d)
    })
}
// Closest points of two closed boxes; their distance is a GLOBAL lower bound
// for the corresponding subset of solid x solid, never an upper bound by itself.
fn closest(a: &Box3, b: &Box3) -> (Point, Point) {
    let pairs: [(Rational, Rational); 3] = std::array::from_fn(|i| {
        if a[i][1] < b[i][0] {
            (a[i][1].clone(), b[i][0].clone())
        } else if b[i][1] < a[i][0] {
            (a[i][0].clone(), b[i][1].clone())
        } else {
            let v = (a[i][0].clone().max(b[i][0].clone()) + a[i][1].clone().min(b[i][1].clone()))
                .mul_pow2(-1);
            (v.clone(), v)
        }
    });
    (pairs.clone().map(|p| p.0), pairs.map(|p| p.1))
}
struct Pair {
    a: Box3,
    b: Box3,
    lower: Rational,
}
impl Pair {
    fn new(a: Box3, b: Box3) -> Self {
        let (p, q) = closest(&a, &b);
        Self {
            a,
            b,
            lower: square_distance(&p, &q),
        }
    }
}
fn root(value: &Rational) -> Enclosure {
    if let Some(r) = crate::profile::sqrt(value) {
        return Enclosure {
            lo: r.clone(),
            hi: r,
        };
    }
    // The lower endpoint is clamped at zero before sqrt: a tiny positive input's
    // outward ball may straddle zero. sqrt of the endpoints is monotone.
    let b = Ball::from_rational(value, 192);
    let lo = Ball::from_rational(&b.lo().max(Rational::zero()), 192).sqrt();
    let hi = Ball::from_rational(&b.hi(), 192).sqrt();
    Enclosure {
        lo: lo.map_or_else(Rational::zero, |v| v.lo().max(Rational::zero())),
        hi: hi.expect("positive upper endpoint").hi(),
    }
}
fn outside(plan: &Plan, b: &Box3, dec: &mut Decider) -> Result<bool, Refusal> {
    let lo = b.clone().map(|v| v[0].clone());
    let hi = b.clone().map(|v| v[1].clone());
    Ok(printing::classify_box(plan, &lo, &hi, dec)? == BoxState::Outside)
}
// Preserve exact axis-box endpoints instead of introducing a ball-width gap
// between a rational boundary and the proposed boundary witness.
fn bounds(plan: &Plan, dec: &mut Decider) -> Result<Box3, Refusal> {
    if let crate::Instruction::Leaf(id) = plan.instructions[plan.root]
        && let crate::geometry::Leaf::Primitive(primitive) = &plan.leaves[id]
        && let Some(bounds) = mass::axis_block(primitive)
    {
        return Ok(bounds);
    }
    printing::bounds(plan, dec)
}
// Keep top-level union pieces separate: their bounding-box hull would fill
// arbitrary voids and make the global product bound needlessly weak.
fn cover(plan: &Plan, limit: usize, dec: &mut Decider) -> Result<Vec<Box3>, Refusal> {
    let mut nodes = vec![plan.root];
    let mut result = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    while let Some(node) = nodes.pop() {
        dec.charge_work(1)?;
        if !seen.insert(node) {
            continue;
        }
        match plan.instructions[node] {
            crate::Instruction::Empty => (),
            crate::Instruction::Boolean {
                kind: catallax_recipe::BooleanKind::Union,
                a,
                b,
            } => {
                nodes.push(b);
                nodes.push(a);
            }
            _ => {
                let mut subtree = plan.clone();
                subtree.root = node;
                let bounds = bounds(&subtree, dec)?;
                if bounds.iter().all(|v| v[0] < v[1]) {
                    if result.len() == limit {
                        return Err(Refusal::new(
                            "assembly.pair-budget",
                            "initial union cover exceeds pair limit",
                        ));
                    }
                    result.push(bounds);
                }
            }
        }
    }
    Ok(result)
}
fn solve(a: &Plan, b: &Plan, options: &Options, dec: &mut Decider) -> Result<Relation, Refusal> {
    if options.gap_width <= Rational::zero() || options.max_pairs == 0 {
        return Err(Refusal::new(
            "assembly.options",
            "positive gap width and pair limit are required",
        ));
    }
    let ac = cover(a, options.max_pairs, dec)?;
    let bc = cover(b, options.max_pairs, dec)?;
    if ac.is_empty() || bc.is_empty() {
        return Err(Refusal::new(
            "assembly.empty",
            "empty operands have no finite smallest gap",
        ));
    }
    if ac
        .len()
        .checked_mul(bc.len())
        .is_none_or(|n| n > options.max_pairs)
    {
        return Err(Refusal::new(
            "assembly.pair-budget",
            "initial Cartesian cover exceeds pair limit",
        ));
    }
    let mut pending = Vec::new();
    for a in &ac {
        for b in &bc {
            pending.push(Pair::new(a.clone(), b.clone()));
        }
    }
    let mut upper: Option<Rational> = None;
    for _ in 0..options.max_pairs {
        dec.charge_work(1)?;
        let index = pending
            .iter()
            .enumerate()
            .min_by(|(_, x), (_, y)| x.lower.cmp(&y.lower))
            .map(|(i, _)| i);
        let Some(index) = index else {
            return Err(Refusal::new(
                "assembly.empty",
                "no nonempty pair remains certified",
            ));
        };
        let pair = pending.remove(index);
        if let Some(upper) = &upper {
            let lower = root(&pair.lower);
            let upper = root(upper);
            if lower.lo > Rational::zero() && &upper.hi - &lower.lo <= options.gap_width {
                return Ok(Relation::Clear {
                    gap: Enclosure {
                        lo: lower.lo,
                        hi: upper.hi,
                    },
                });
            }
        }
        let (p, q) = closest(&pair.a, &pair.b);
        for (p, q) in [(p, q), (middle(&pair.a), middle(&pair.b))] {
            let ap = location(a, &p, dec)?;
            let bp = location(b, &p, dec)?;
            let aq = location(a, &q, dec)?;
            let bq = location(b, &q, dec)?;
            for (witness, x, y) in [(&p, ap, bp), (&q, aq, bq)] {
                if x == Location::Inside && y == Location::Inside {
                    return Ok(Relation::Clash {
                        witness: witness.clone(),
                    });
                }
            }
            if material(ap) && material(bq) {
                let distance = square_distance(&p, &q);
                upper = Some(upper.map_or(distance.clone(), |old| old.min(distance)));
            }
        }
        // Re-check immediately, so an analytic box pair needs only one visit.
        if let Some(upper) = &upper {
            let lower_squared = pending
                .iter()
                .fold(pair.lower.clone(), |v, p| v.min(p.lower.clone()));
            let lower = root(&lower_squared);
            let upper = root(upper);
            if lower.lo > Rational::zero() && &upper.hi - &lower.lo <= options.gap_width {
                return Ok(Relation::Clear {
                    gap: Enclosure {
                        lo: lower.lo,
                        hi: upper.hi,
                    },
                });
            }
            if upper.hi.is_zero() {
                return Err(Refusal::new(
                    "assembly.contact",
                    "closed solids meet; positive-volume clash is not established",
                ));
            }
        }
        // Both tests must succeed before the pair can be subdivided. An unknown
        // classifier result keeps the entire closed box in the cover.
        if outside(a, &pair.a, dec)? || outside(b, &pair.b, dec)? {
            continue;
        }
        let mut side = 0;
        let mut axis = 0;
        let mut width = Rational::zero();
        for (s, bounds) in [&pair.a, &pair.b].iter().enumerate() {
            for (i, v) in bounds.iter().enumerate() {
                let w = &v[1] - &v[0];
                if w > width {
                    side = s;
                    axis = i;
                    width = w;
                }
            }
        }
        let mut first = [pair.a, pair.b];
        let mut second = first.clone();
        let mid = (&first[side][axis][0] + &first[side][axis][1]).mul_pow2(-1);
        first[side][axis][1] = mid.clone();
        second[side][axis][0] = mid;
        for [a, b] in [first, second] {
            pending.push(Pair::new(a, b));
        }
    }
    Err(Refusal::new(
        "assembly.pair-budget",
        "pair limit reached before a clash or requested global gap certificate",
    ))
}
