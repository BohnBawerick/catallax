//! Bounded assembly queries on placed, complete field definitions.
//!
//! Distance lower bounds cover the Cartesian product of both solids. Upper
//! bounds use certified points in the closed solids. No residual, preview or
//! nearest-owner vector is used as a distance certificate.
//!
//! Box pairs are evaluated in fixed rounds on every core; the request's thread
//! adopts each result in order, so answers, refusals and operation work are
//! the same at every thread count.
use crate::{Location, Plan, Workspace, lattice::BoxState, mass, printing};
use catallax_measure::Enclosure;
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{Pose, RPoint3};
use std::ops::ControlFlow;

mod ahead;

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
    /// The placed field, in world millimetres.
    pub fn plan(&self) -> &Plan {
        &self.plan
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

/// Compare two placed parts. Touching is unresolved, never labelled clear or clash:
/// a proven zero gap refuses `assembly.contact` once the pair limit is spent
/// without finding a point inside both.
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

/// Opt-in exact meridian integration for coaxial polygonal revolutions.
/// Density and property accuracies have the same units as [`uniform_mass`].
/// Unsupported meridians refuse without substituting their bounding envelopes.
pub fn revolved_uniform_mass(
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
    let properties = mass::revolved_properties(&part.plan, accuracy, options, dec)?;
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
// One pending subset of solid x solid. `order` breaks ties in creation order,
// so the search visits pairs in the same sequence at every thread count.
struct Pair {
    a: Box3,
    b: Box3,
    lower: Rational,
    order: u64,
}
impl Pair {
    fn new(a: Box3, b: Box3, order: &mut u64) -> Self {
        let (p, q) = closest(&a, &b);
        *order += 1;
        Self {
            lower: square_distance(&p, &q),
            a,
            b,
            order: *order,
        }
    }
}
impl PartialEq for Pair {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}
impl Eq for Pair {}
impl PartialOrd for Pair {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
// Reversed, so `BinaryHeap` pops the smallest lower bound, earliest first.
impl Ord for Pair {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (&other.lower, other.order).cmp(&(&self.lower, self.order))
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
/// Box pairs taken from the cover per round. A fixed count, independent of the
/// thread count, so every count visits the same pairs in the same order.
const ROUND: usize = 32;

/// Proven facts about one box pair, all decided from that pair alone.
struct Evaluation {
    /// Proposal pairs `(p, q)` with the locations of `p` and `q` in each solid:
    /// `[p in a, p in b, q in a, q in b]`.
    proposals: [(Point, Point, [Location; 4]); 2],
    /// Material points of each solid moved from the box centre toward the
    /// other box, when both centres are material.
    approach: Option<(Point, Point)>,
    /// One closed box is proven strictly outside its solid.
    outside: bool,
}
// Bisects from a certified material point toward `target`, keeping the material
// end, until the bracket is no longer than `width`. Every returned point is
// certified material; the search only proposes, it never bounds from below.
fn toward(
    plan: &Plan,
    material_point: &Point,
    target: &Point,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Point, Refusal> {
    let mut lo = material_point.clone();
    let mut hi = target.clone();
    if material(location(plan, &hi, dec)?) {
        return Ok(hi);
    }
    let width = width * width;
    while square_distance(&lo, &hi) > width {
        let mid = std::array::from_fn(|i| (&lo[i] + &hi[i]).mul_pow2(-1));
        if material(location(plan, &mid, dec)?) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(lo)
}
fn evaluate(
    a: &Plan,
    b: &Plan,
    pair: &Pair,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Evaluation, Refusal> {
    let (p, q) = closest(&pair.a, &pair.b);
    let mut proposals = Vec::with_capacity(2);
    for (p, q) in [(p.clone(), q.clone()), (middle(&pair.a), middle(&pair.b))] {
        let at = [
            location(a, &p, dec)?,
            location(b, &p, dec)?,
            location(a, &q, dec)?,
            location(b, &q, dec)?,
        ];
        proposals.push((p, q, at));
    }
    let (ma, mb, [ama, _, _, bmb]) = &proposals[1];
    let approach = if material(*ama) && material(*bmb) {
        // A quarter of the requested width keeps the proposal's own share of
        // the enclosure small; the lower bound comes only from the cover.
        let step = width.mul_pow2(-2);
        Some((
            toward(a, ma, &p, &step, dec)?,
            toward(b, mb, &q, &step, dec)?,
        ))
    } else {
        None
    };
    // An unknown classifier result keeps the entire closed box in the cover.
    let outside = outside(a, &pair.a, dec)? || outside(b, &pair.b, dec)?;
    Ok(Evaluation {
        proposals: proposals.try_into().unwrap_or_else(|_| unreachable!()),
        approach,
        outside,
    })
}
fn clear(lower_squared: &Rational, upper_squared: &Rational, width: &Rational) -> Option<Relation> {
    let lower = root(lower_squared);
    let upper = root(upper_squared);
    (lower.lo > Rational::zero() && &upper.hi - &lower.lo <= *width).then_some(Relation::Clear {
        gap: Enclosure {
            lo: lower.lo,
            hi: upper.hi,
        },
    })
}
/// The cover of solid x solid still pending, and the best certified upper bound.
struct Search<'w> {
    pending: std::collections::BinaryHeap<Pair>,
    order: u64,
    /// Squared distance between two certified material points, one per solid.
    upper: Option<Rational>,
    width: &'w Rational,
}
impl Search<'_> {
    /// Smallest lower bound of the cover while `pair` is still in it. Pairs
    /// leave a round in order and their bounds never decrease within it; a
    /// child's bound is at least its parent's. So the minimum is the pair's
    /// own, the round's next or the heap's.
    fn lower(&self, pair: &Pair, next: Option<&Pair>) -> Rational {
        [next, self.pending.peek()]
            .into_iter()
            .flatten()
            .fold(pair.lower.clone(), |v, p| v.min(p.lower.clone()))
    }
    fn clear(&self, pair: &Pair, next: Option<&Pair>) -> Option<Relation> {
        let upper = self.upper.as_ref()?;
        clear(&self.lower(pair, next), upper, self.width)
    }
    /// Applies one pair's evaluation: a clash or clear answer, or a split.
    fn visit(
        &mut self,
        pair: &Pair,
        next: Option<&Pair>,
        evaluation: Result<Evaluation, Refusal>,
    ) -> Result<Option<Relation>, Refusal> {
        if let Some(clear) = self.clear(pair, next) {
            return Ok(Some(clear));
        }
        let evaluation = evaluation?;
        let mut certified = Vec::with_capacity(3);
        for (p, q, [ap, bp, aq, bq]) in &evaluation.proposals {
            for (witness, x, y) in [(p, ap, bp), (q, aq, bq)] {
                if *x == Location::Inside && *y == Location::Inside {
                    return Ok(Some(Relation::Clash {
                        witness: witness.clone(),
                    }));
                }
            }
            if material(*ap) && material(*bq) {
                certified.push((p, q));
            }
        }
        certified.extend(evaluation.approach.as_ref().map(|(p, q)| (p, q)));
        for (p, q) in certified {
            let distance = square_distance(p, q);
            self.upper = Some(match self.upper.take() {
                Some(old) => old.min(distance),
                None => distance,
            });
        }
        // Re-check immediately, so an analytic box pair needs only one visit.
        if let Some(clear) = self.clear(pair, next) {
            return Ok(Some(clear));
        }
        if evaluation.outside {
            return Ok(None);
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
        let mut first = [pair.a.clone(), pair.b.clone()];
        let mut second = first.clone();
        let mid = (&first[side][axis][0] + &first[side][axis][1]).mul_pow2(-1);
        first[side][axis][1] = mid.clone();
        second[side][axis][0] = mid;
        for [a, b] in [first, second] {
            let pair = Pair::new(a, b, &mut self.order);
            self.pending.push(pair);
        }
        Ok(None)
    }
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
    let mut search = Search {
        pending: std::collections::BinaryHeap::new(),
        order: 0,
        upper: None,
        width: &options.gap_width,
    };
    for a in &ac {
        for b in &bc {
            let pair = Pair::new(a.clone(), b.clone(), &mut search.order);
            search.pending.push(pair);
        }
    }
    let mut processed = 0;
    while processed < options.max_pairs {
        let mut round = Vec::with_capacity(ROUND);
        while round.len() < ROUND.min(options.max_pairs - processed)
            && let Some(pair) = search.pending.pop()
        {
            round.push(pair);
        }
        if round.is_empty() {
            return Err(Refusal::new(
                "assembly.empty",
                "no nonempty pair remains certified",
            ));
        }
        let result = ahead::in_order(
            &round,
            |pair, dec| evaluate(a, b, pair, &options.gap_width, dec),
            dec,
            |index, evaluation, dec| {
                processed += 1;
                let visit = dec
                    .charge_work(1)
                    .and_then(|()| search.visit(&round[index], round.get(index + 1), evaluation));
                match visit {
                    Ok(None) => ControlFlow::Continue(()),
                    Ok(Some(relation)) => ControlFlow::Break(Ok(relation)),
                    Err(refusal) => ControlFlow::Break(Err(refusal)),
                }
            },
        );
        if let ControlFlow::Break(result) = result {
            return result;
        }
    }
    let Search { pending, upper, .. } = search;
    if upper.as_ref().is_some_and(Rational::is_zero) {
        return Err(Refusal::new(
            "assembly.contact",
            "closed solids meet at a certified common point; no point inside both was found within the pair limit",
        ));
    }
    // The reached enclosure is reported, never published as an answer: it is
    // wider than requested, and its lower end may be zero.
    // Six significant digits, rounded outward.
    let lower = pending.peek().map_or_else(
        || "0".to_owned(),
        |p| root(&p.lower).lo.to_decimal_bounds(6).0,
    );
    let reached = upper.map_or_else(
        || format!("smallest gap at least {lower} mm; no certified upper bound"),
        |upper| {
            let upper = root(&upper).hi.to_decimal_bounds(6).1;
            format!("smallest gap in [{lower}, {upper}] mm")
        },
    );
    Err(Refusal::new(
        "assembly.pair-budget",
        format!("pair limit reached before a clash or requested global gap certificate; {reached}"),
    ))
}
