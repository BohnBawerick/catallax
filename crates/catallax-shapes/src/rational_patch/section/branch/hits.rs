//! Where a section branch meets a third analytic surface.

use super::*;
use catallax_numbers::{BoxRoots, BudgetPolicy, Outcome, PolySystem, RootBox};

/// Deepest box split while isolating the meetings within one arc cell.
const MAX_HIT_DEPTH: u32 = 24;
/// Depth before the first root certificate is tried in a box.
const MIN_CERTIFY_DEPTH: u32 = 10;
/// Working precision of one root certificate. An inconclusive box splits
/// instead of escalating precision, which keeps failed boxes cheap.
const CERTIFY_BITS: u32 = 128;
/// Work cap of one root certificate.
const CERTIFY_WORK: u64 = 100_000;
/// Width of the rational box enclosing an exact meeting at a cell end.
const END_BOX_BITS: u32 = 64;

/// One meeting: its cell, its fraction along that cell, its parameters, and
/// a rational box holding it.
struct Hit {
    cell: usize,
    fraction: Real,
    at: [Real; 2],
    near: Domain,
}

/// Where a section branch meets a third analytic surface.
#[derive(Clone, Debug)]
pub enum SectionHits {
    /// The branch lies on the surface.
    Contained,
    /// Transverse crossings, in traversal order.
    Points(Vec<Point3>),
}

impl PatchSectionBranch {
    /// Every point where this branch meets `other`, in traversal order.
    ///
    /// Meetings at cell ends, including chart seams, are decided exactly by
    /// the sign of `other`'s cut equation there. Inside each arc cell the two
    /// cut equations form a square rational system. A box where either has a
    /// single-signed Bernstein net holds no meeting. Any other box either
    /// certifies one nonsingular root or splits. A box touching a known
    /// meeting that solves this cell's system is certified together with
    /// that meeting's box: one certified root there is the known one, so no
    /// root is counted twice and no equality of isolated roots is ever needed. A nonsingular root is a
    /// transverse crossing. Where the two cut equations are proportional the
    /// branch lies on `other`: over every cell that is `Contained`, over only
    /// some it refuses. Tangent meetings never certify and refuse with
    /// `intersect.patch-section-hit`.
    /// Algebraic cut coefficients use complete exact eliminants instead of
    /// rational root boxes, and check a nonzero Jacobian at every meeting.
    pub fn surface_hits(
        &self,
        other: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<SectionHits, Refusal> {
        let before_separation = dec.stats().coincident_separation;
        let mut algebraic = false;
        other.validate(dec)?;
        // Admit rational plane coefficients before combining their algebraic
        // definitions with an independent section-end root. The reductions
        // are exact numeric proofs; an inconclusive coefficient is unchanged.
        let reduced = match other {
            SectionSurface::Plane(plane)
                if [
                    &plane.origin.x,
                    &plane.origin.y,
                    &plane.origin.z,
                    &plane.normal.x,
                    &plane.normal.y,
                    &plane.normal.z,
                ]
                .iter()
                .any(|value| value.rational_ref().is_none()) =>
            {
                let mut plane = plane.clone();
                for value in [
                    &mut plane.origin.x,
                    &mut plane.origin.y,
                    &mut plane.origin.z,
                    &mut plane.normal.x,
                    &mut plane.normal.y,
                    &mut plane.normal.z,
                ] {
                    if value.rational_ref().is_none()
                        && let Some(rational) = dec.rational_value(value)?
                    {
                        *value = Real::from(rational);
                    }
                }
                Some(SectionSurface::Plane(plane))
            }
            _ => None,
        };
        let other = reduced.as_ref().unwrap_or(other);
        if matches!(other, SectionSurface::Cone(_)) {
            let mut seen = std::collections::BTreeSet::new();
            let mut retained = None;
            for cell in self.cells.iter() {
                if !seen.insert((cell.chart, cell.span)) {
                    continue;
                }
                let controls = self.charts[cell.chart]
                    .carrier()
                    .span_bezier_controls(self.span(cell).indices, dec)?;
                let inside = other.retains_span(&controls, dec)?;
                if retained.is_some_and(|old| old != inside) {
                    return Err(Refusal::new(
                        "intersect.patch-cone-nappe",
                        "third-surface meetings need one strict cone nappe over the branch",
                    ));
                }
                retained = Some(inside);
            }
            if retained == Some(false) {
                return Ok(SectionHits::Points(Vec::new()));
            }
        }
        let mut contained = 0;
        for cell in self.cells.iter() {
            let span = self.span(cell);
            if proportional(&self.surface.equation(span), &other.equation(span), dec)? {
                contained += 1;
            }
        }
        if contained == self.cells.len() {
            return Ok(SectionHits::Contained);
        }
        if contained > 0 {
            return Err(unresolved(
                "the branch lies on the third surface over part of its length",
            ));
        }
        let mut hits: Vec<Hit> = Vec::new();
        // Exact meetings at cell ends. Consecutive cells share their ends, so
        // each cell contributes its first end and the last cell its second.
        for (index, cell) in self.cells.iter().enumerate() {
            let last = index + 1 == self.cells.len();
            for end in 0..if last && !self.closed { 2 } else { 1 } {
                let at = &cell.ends[end];
                let cut = other.equation(self.span(cell));
                if exact_sign(&cut.eval_real(at), dec)? == Sign::Zero {
                    let mut near = cell.domain.clone();
                    for (axis, value) in at.iter().enumerate() {
                        let ball = dec.enclose_budgeted_at(value, END_BOX_BITS)?;
                        near[axis] = [ball.lo(), ball.hi()];
                    }
                    hits.push(Hit {
                        cell: index,
                        fraction: Real::from_i64(end as i64),
                        at: at.clone(),
                        near,
                    });
                }
            }
        }
        for (index, cell) in self.cells.iter().enumerate() {
            let patch = &self.charts[cell.chart];
            let span = self.span(cell);
            let equations = [self.surface.equation(span), other.equation(span)];
            let (Some(own), Some(cut)) = (equations[0].as_rational(), equations[1].as_rational())
            else {
                algebraic = true;
                let jacobian = &equations[0].partial(0) * &equations[1].partial(1)
                    - &equations[0].partial(1) * &equations[1].partial(0);
                let roots = match crate::rational_patch::query::linear_roots(
                    &equations,
                    &cell.domain,
                    dec,
                )? {
                    Some(roots) => roots,
                    None => crate::rational_patch::query::roots(&equations, &cell.domain, dec)?,
                };
                for at in roots {
                    if exact_sign(&jacobian.eval_real(&at), dec)? == Sign::Zero {
                        return Err(unresolved(
                            "an algebraic meeting has no strict transverse Jacobian",
                        ));
                    }
                    if same_uv(&at, &cell.ends[0], dec)? || same_uv(&at, &cell.ends[1], dec)? {
                        continue;
                    }
                    let free = 1 - cell.solved;
                    let fraction = (&at[free] - &cell.ends[0][free])
                        / (&cell.ends[1][free] - &cell.ends[0][free]);
                    let mut near = cell.domain.clone();
                    for (axis, value) in at.iter().enumerate() {
                        let bound = dec.enclose_budgeted_at(value, END_BOX_BITS)?;
                        near[axis] = [bound.lo(), bound.hi()];
                    }
                    hits.push(Hit {
                        cell: index,
                        fraction,
                        at,
                        near,
                    });
                }
                continue;
            };
            let system = Arc::new(PolySystem::new(vec![own, cut]).expect("two equations"));
            // Exact meeting coordinates, from each axis's eliminant, on demand.
            let mut eliminants: Option<[Vec<Real>; 2]> = None;
            let mut pending = vec![(cell.domain.clone(), 0u32)];
            while let Some((domain, depth)) = pending.pop() {
                dec.check_deadline()?;
                let controls = patch.carrier().span_bezier_controls_on(
                    span.indices,
                    domain.each_ref().map(|b| [&b[0], &b[1]]),
                    dec,
                )?;
                let degrees = patch.carrier().degrees();
                let (own_net, _) = self.surface.net(&controls, degrees, dec)?;
                let (cut_net, _) = other.net(&controls, degrees, dec)?;
                if uniform(&own_net, dec)?.is_some() || uniform(&cut_net, dec)?.is_some() {
                    continue;
                }
                // Krawczyk contracts slowly on wide boxes; split those first.
                if depth < MIN_CERTIFY_DEPTH {
                    let [low, high] = halves(domain);
                    pending.push((high, depth + 1));
                    pending.push((low, depth + 1));
                    continue;
                }
                // Inflated, so a root on this box's edge is interior to the
                // candidate; Krawczyk inclusion fails for a boundary root.
                let grown = inflated(&domain);
                let known = self.reusable(&hits, &grown, cell, &equations, dec)?;
                let candidate = match known {
                    Some(k) => inflated(&hull(&domain, &hits[k].near)),
                    None => grown,
                };
                let candidate = RootBox::new(
                    candidate.iter().map(|b| b[0].clone()).collect(),
                    candidate.iter().map(|b| b[1].clone()).collect(),
                )
                .expect("ordered box");
                // Krawczyk can creep, shrinking a box in tiny steps without
                // ever certifying. A per-certificate work cap makes such a box
                // inconclusive, so it splits. The request's operation
                // allowance still pays for every unit.
                let mut certifier = Decider::new(BudgetPolicy {
                    max_work: CERTIFY_WORK,
                    ..*dec.policy()
                });
                match certifier.certify_root_at(&system, &candidate, CERTIFY_BITS) {
                    // Unique in a box holding the known meeting: it is that one.
                    Outcome::Proven(BoxRoots::One(_)) if known.is_some() => continue,
                    Outcome::Proven(BoxRoots::One(root)) => {
                        let (lo, hi) = (root.root_box().lo(), root.root_box().hi());
                        let near: Domain = [
                            [lo[0].clone(), hi[0].clone()],
                            [lo[1].clone(), hi[1].clone()],
                        ];
                        let outside = (0..2).any(|a| {
                            near[a][1] < cell.domain[a][0] || near[a][0] > cell.domain[a][1]
                        });
                        if outside {
                            continue;
                        }
                        let interior = (0..2).all(|a| {
                            near[a][0] > cell.domain[a][0] && near[a][1] < cell.domain[a][1]
                        });
                        if interior {
                            if eliminants.is_none() {
                                eliminants = Some(axis_roots(&system, dec)?);
                            }
                            let roots = eliminants.as_ref().expect("computed above");
                            let Some(at) = exact_meeting(roots, &near, dec)? else {
                                // Another eliminant root shares this box; split.
                                if depth >= MAX_HIT_DEPTH {
                                    return Err(unresolved(
                                        "a certified meeting has no separated exact coordinates",
                                    ));
                                }
                                let [low, high] = halves(domain);
                                pending.push((high, depth + 1));
                                pending.push((low, depth + 1));
                                continue;
                            };
                            let free = 1 - cell.solved;
                            let fraction = (&at[free] - &cell.ends[0][free])
                                / (&cell.ends[1][free] - &cell.ends[0][free]);
                            hits.push(Hit {
                                cell: index,
                                fraction,
                                at,
                                near,
                            });
                            continue;
                        }
                        // The root box straddles the border. A root of both
                        // equations on an arc cell's border is a cell end,
                        // already known; one root in a box with that end is it.
                        // Otherwise the box splits until the root's place is
                        // decided by its certified box alone.
                        let around = inflated(&near);
                        if let Some(k) = self.reusable(&hits, &around, cell, &equations, dec)? {
                            let both = inflated(&hull(&near, &hits[k].near));
                            let both = RootBox::new(
                                both.iter().map(|b| b[0].clone()).collect(),
                                both.iter().map(|b| b[1].clone()).collect(),
                            )
                            .expect("ordered box");
                            let mut certifier = Decider::new(BudgetPolicy {
                                max_work: CERTIFY_WORK,
                                ..*dec.policy()
                            });
                            if let Outcome::Proven(BoxRoots::One(_)) =
                                certifier.certify_root_at(&system, &both, CERTIFY_BITS)
                            {
                                continue;
                            }
                        }
                    }
                    Outcome::Proven(BoxRoots::None) if known.is_none() => continue,
                    Outcome::Refused(error) if error.code == "numbers.operation-budget" => {
                        return Err(error);
                    }
                    _ => (),
                }
                if depth >= MAX_HIT_DEPTH {
                    let point = span.point(&midpoint(&domain).map(Real::from_rational));
                    return Err(unresolved(&format!(
                        "an uncertified tangent or contained meeting near ({:.6}, {:.6}, {:.6})",
                        point.x.to_f64(),
                        point.y.to_f64(),
                        point.z.to_f64()
                    )));
                }
                let [low, high] = halves(domain);
                pending.push((high, depth + 1));
                pending.push((low, depth + 1));
            }
        }
        // Traversal order: cell index, then the fraction along the cell.
        let mut ordered: Vec<Hit> = Vec::with_capacity(hits.len());
        for hit in hits {
            let mut position = ordered.len();
            for (i, old) in ordered.iter().enumerate() {
                if hit.cell < old.cell
                    || (hit.cell == old.cell
                        && exact_sign(&(&hit.fraction - &old.fraction), dec)? == Sign::Negative)
                {
                    position = i;
                    break;
                }
            }
            ordered.insert(position, hit);
        }
        if algebraic && dec.stats().coincident_separation != before_separation {
            return Err(unresolved(
                "algebraic meetings need exact root identities and cell ownership",
            ));
        }
        Ok(SectionHits::Points(
            ordered
                .into_iter()
                .map(|hit| self.span(&self.cells[hit.cell]).point(&hit.at))
                .collect(),
        ))
    }
}

/// The box grown by a quarter of its width on every side.
fn inflated(domain: &Domain) -> Domain {
    domain.clone().map(|[lo, hi]| {
        let margin = (&hi - &lo).mul_pow2(-2);
        [&lo - &margin, &hi + &margin]
    })
}

/// Real roots of each axis's eliminant: the resultant of the two equations
/// with the other parameter eliminated, as exact univariate algebraic numbers.
fn axis_roots(system: &PolySystem, dec: &mut Decider) -> Result<[Vec<Real>; 2], Refusal> {
    let [own, cut] = [&system.equations()[0], &system.equations()[1]];
    let mut roots: [Vec<Real>; 2] = [Vec::new(), Vec::new()];
    for (axis, slot) in roots.iter_mut().enumerate() {
        let eliminant = own.resultant(cut, 1 - axis);
        let terms = eliminant
            .terms()
            .into_iter()
            .map(|(m, c)| (vec![m.exponents()[axis]], Real::from_rational(c)))
            .collect();
        let eliminant = RealPoly::from_terms(1, terms);
        if eliminant.is_zero(dec)? {
            return Err(unresolved(
                "the cut equations share a common factor over this cell",
            ));
        }
        *slot = eliminant.real_roots(dec)?;
    }
    Ok(roots)
}

/// The exact meeting in a certified box: the one eliminant root on each axis
/// inside it. `None` when an axis has another root in the box too.
fn exact_meeting(
    roots: &[Vec<Real>; 2],
    near: &Domain,
    dec: &mut Decider,
) -> Result<Option<[Real; 2]>, Refusal> {
    let mut at: [Option<Real>; 2] = [None, None];
    for axis in 0..2 {
        for root in &roots[axis] {
            let [lo, hi] = &near[axis];
            if exact_sign(&(root - &Real::from_rational(lo.clone())), dec)? == Sign::Negative
                || exact_sign(&(&Real::from_rational(hi.clone()) - root), dec)? == Sign::Negative
            {
                continue;
            }
            if at[axis].is_some() {
                return Ok(None);
            }
            at[axis] = Some(root.clone());
        }
    }
    match at {
        [Some(u), Some(v)] => Ok(Some([u, v])),
        _ => Err(unresolved(
            "a certified meeting has no eliminant root in its box",
        )),
    }
}

/// Halves a box across its longer side.
fn halves(domain: Domain) -> [Domain; 2] {
    let axis = usize::from((&domain[0][1] - &domain[0][0]) < (&domain[1][1] - &domain[1][0]));
    let middle = (&domain[axis][0] + &domain[axis][1]).mul_pow2(-1);
    let mut low = domain.clone();
    let mut high = domain;
    low[axis][1] = middle.clone();
    high[axis][0] = middle;
    [low, high]
}

/// Whether `b` is a constant multiple of the nonzero polynomial `a`.
fn proportional(a: &RealPoly, b: &RealPoly, dec: &mut Decider) -> Result<bool, Refusal> {
    let terms = b.terms();
    for (monomial, ca) in a.terms() {
        if exact_sign(&ca, dec)? == Sign::Zero {
            continue;
        }
        let cb = terms
            .iter()
            .find(|(m, _)| m.exponents() == monomial.exponents())
            .map_or_else(Real::zero, |(_, c)| c.clone());
        return (&b.scale(&ca) - &a.scale(&cb)).is_zero(dec);
    }
    Ok(false)
}

fn overlaps(a: &Domain, b: &Domain) -> bool {
    (0..2).all(|axis| a[axis][0] <= b[axis][1] && b[axis][0] <= a[axis][1])
}

fn hull(a: &Domain, b: &Domain) -> Domain {
    std::array::from_fn(|axis| {
        [
            a[axis][0].clone().min(b[axis][0].clone()),
            a[axis][1].clone().max(b[axis][1].clone()),
        ]
    })
}

impl PatchSectionBranch {
    /// A known meeting near `domain` that is a root of the current cell's
    /// two equations. Only such a meeting may stand for a box's unique root.
    /// A meeting from the same chart span solves the same system. One from
    /// another span must be proved a root here: it need not solve these
    /// polynomials, and a distinct meeting here must not be discarded for it.
    fn reusable(
        &self,
        hits: &[Hit],
        domain: &Domain,
        cell: &ArcCell,
        equations: &[RealPoly; 2],
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        for (k, hit) in hits.iter().enumerate() {
            if !overlaps(&hit.near, domain) {
                continue;
            }
            let owner = &self.cells[hit.cell];
            if owner.chart == cell.chart && owner.span == cell.span {
                return Ok(Some(k));
            }
            let mut root = true;
            for equation in equations {
                if exact_sign(&equation.eval_real(&hit.at), dec)? != Sign::Zero {
                    root = false;
                    break;
                }
            }
            if root {
                return Ok(Some(k));
            }
        }
        Ok(None)
    }
}

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("intersect.patch-section-hit", detail)
}
