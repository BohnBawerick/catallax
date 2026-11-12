//! Face systems of a cell and its cover: exclusion, exact border points and
//! how the branch crosses each of them.

use super::*;
use catallax_numbers::{Outcome, StrictSign};
use net::{disjoint, less};
use std::collections::BTreeMap;

/// One end of an arc on the border of its unit.
#[derive(Clone, Debug)]
pub(super) struct Endpoint {
    pub(super) hit: Hit,
    /// The unit's border hyperplanes `(axis, value)` through the point.
    pub(super) active: Vec<(usize, Rational)>,
    /// Whether the arc starts here in the direction of first normal cross
    /// second normal; otherwise it ends here.
    pub(super) start: bool,
}

/// An exactly isolated meeting point on a face of a unit.
#[derive(Clone, Debug)]
pub(super) struct Hit {
    /// Exact coordinates `(axis, value)`, in axis order.
    pub(super) fixed: Vec<(usize, Rational)>,
    /// Axes of the root's variables, in order.
    pub(super) free: Vec<usize>,
    /// The root, in variables centred on `centre`.
    pub(super) root: IsolatedRoot,
    pub(super) centre: Vec<Rational>,
    /// The face equations in the free variables, before centring.
    pub(super) equations: Arc<Vec<MPoly>>,
    /// A box round the point in `(u, v, s, t)`.
    pub(super) bounds: [[Rational; 2]; 4],
    /// Boxes round the offset lengths `L`, one per offset side.
    pub(super) lengths: Vec<[Rational; 2]>,
    /// The face box in which the point is the only root of its equations.
    pub(super) face: [[Rational; 2]; 4],
    pub(super) spans: [usize; 2],
}

impl Hit {
    pub(super) fn parameters(&self) -> [Real; 4] {
        std::array::from_fn(|axis| {
            if let Some((_, value)) = self.fixed.iter().find(|(e, _)| *e == axis) {
                Real::from_rational(value.clone())
            } else {
                let index = self
                    .free
                    .iter()
                    .position(|e| *e == axis)
                    .expect("every axis is fixed or free");
                &Real::from_rational(self.centre[index].clone()) + &self.root.coordinate(index)
            }
        })
    }
}

/// Work cap of one root certificate, and its precision.
const CERTIFY_WORK: u64 = 50_000;
const CERTIFY_BITS: u32 = 64;
/// Krawczyk steps one root certificate of an offset face may take. A face
/// system with offset lengths has four or five unknowns, and one step on it
/// can cost more than the whole plain cap; its cap is this many of its own
/// steps instead. The request's allowance still pays every unit.
const OFFSET_STEPS: u64 = 16;

/// The work cap of one root certificate on `system`.
fn certify_cap(system: &PolySystem) -> u64 {
    let n = system.dimension() as u64;
    if n <= 3 {
        return CERTIFY_WORK;
    }
    let terms = system
        .equations()
        .iter()
        .map(|p| p.terms().len() as u64)
        .sum::<u64>();
    // One step's charge at 64 bits, as the certificate's decider counts it.
    let step = n
        .saturating_pow(3)
        .saturating_add(terms.saturating_mul(n).saturating_mul(8))
        .saturating_mul(16);
    step.saturating_mul(OFFSET_STEPS).max(CERTIFY_WORK)
}

pub(super) fn strict(interval: &[Rational; 2]) -> bool {
    interval[0].sign() == Sign::Positive || interval[1].sign() == Sign::Negative
}

pub(super) fn apart(
    a: &[Rational; 2],
    b: &[Rational; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(less(&a[1], &b[0], dec)? || less(&b[1], &a[0], dec)?)
}

/// Whether a cover cell can leave a first-patch cell's cover. Disjoint
/// boxes, or disjoint slabs along either cell's normal direction, prove that
/// the closed cells meet nowhere.
///
/// Boxes that touch exactly prove less. Every meeting point then lies on the
/// touching side, and where one patch reaches that side on one border of its
/// cell only, on that border. That is no meeting inside the pair, not none in
/// the closed cells. Each patch's cells own their points half-open: a cell
/// owns its low borders and leaves its high borders to the neighbour, except
/// on the physical border, which it owns. So the pair is dropped only when
/// that border is an internal high border. Its meeting points then belong to
/// the neighbouring cell across it, and a pair that owns a meeting point is
/// never dropped this way.
pub(super) fn excluded(
    ctx: &Context,
    a: &Cell,
    b: &Cell,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if disjoint(&a.bbox, &b.bbox, dec)? {
        return Ok(true);
    }
    // A binary64 estimate skips a slab test that clearly fails; only exact
    // arithmetic excludes.
    if (!a.clearly_meets(b, None, dec)? && apart(&a.slab(&b.frame[0], None, dec)?, &b.own[0], dec)?)
        || (!b.clearly_meets(a, None, dec)?
            && apart(&a.own[0], &b.slab(&a.frame[0], None, dec)?, dec)?)
    {
        return Ok(true);
    }
    // Boxes that touch exactly: only cells on knot lines, where the two
    // patches can share a value. The outward boxes then overlap by at most
    // their rounding, so the exact boxes decide. The touch rule reads exact
    // control nets, so it applies to rational patches only.
    if a.offset.is_some() || b.offset.is_some() {
        return Ok(false);
    }
    for axis in 0..3 {
        for high in [true, false] {
            let gap = net::sub(
                &a.bbox[axis][usize::from(high)],
                &b.bbox[axis][usize::from(!high)],
                dec,
            )?;
            let rounding = Rational::one().mul_pow2(-60);
            if less(&rounding, &gap.abs(), dec)? {
                continue;
            }
            let touch = a.exact_box(dec)?[axis][usize::from(high)].clone();
            let other = b.exact_box(dec)?[axis][usize::from(!high)].clone();
            dec.reserve_exact_arithmetic(&touch, &other, true)?;
            if touch == other
                && (a.extreme_on_border(axis, high, &passed(ctx, a, 0), dec)?
                    || b.extreme_on_border(axis, !high, &passed(ctx, b, 2), dec)?)
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// The borders whose points a cell leaves to its neighbour: its high `u` and
/// high `v` borders, edges 1 and 3, where they are not physical. `first` is 0
/// for a first-patch cell and 2 for a second-patch one.
fn passed(ctx: &Context, cell: &Cell, first: usize) -> Vec<usize> {
    let mut edges = Vec::with_capacity(2);
    for (axis, edge) in [(0, 1), (1, 3)] {
        if !ctx.physical(first + axis, &cell.rect[axis][1]) {
            edges.push(edge);
        }
    }
    edges
}

/// Exclusion of one face: the border curve against the other cell, by box
/// and by slab along the other cell's normal direction.
pub(super) fn face_excluded(
    border_cell: &Cell,
    edge: usize,
    other: &Cell,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if disjoint(&border_cell.edges[edge], &other.bbox, dec)? {
        return Ok(true);
    }
    if !border_cell.clearly_meets(other, Some(edge), dec)?
        && apart(
            &border_cell.slab(&other.frame[0], Some(edge), dec)?,
            &other.own[0],
            dec,
        )?
    {
        return Ok(true);
    }
    Ok(false)
}

/// The simplest rational in a closed interval of positive or negative
/// numbers: the continued-fraction rule. Integers come first.
pub(super) fn simplest(
    lo: &Rational,
    hi: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    dec.charge_work(4)?;
    if lo.sign() != Sign::Positive && hi.sign() != Sign::Negative {
        return Ok(Rational::zero());
    }
    if hi.sign() == Sign::Negative {
        return Ok(-simplest(&-hi.clone(), &-lo.clone(), dec)?);
    }
    let ceiling = Rational::from_integer(lo.ceil());
    if &ceiling <= hi {
        return Ok(ceiling);
    }
    let whole = Rational::from_integer(lo.floor());
    let low = (hi - &whole).recip().expect("positive remainder");
    let high = (lo - &whole).recip().expect("positive remainder");
    let inner = simplest(&low, &high, dec)?;
    Ok(&whole + &inner.recip().expect("positive inner"))
}

/// The tangent-contact refusal when the simplest rational parameters in the
/// pair are an exact meeting point with parallel normals. No certificate can
/// ever accept a pair holding such a point: both fields vanish there, the
/// normal cross product vanishes and no exclusion holds. So the pair holding
/// it refuses at once rather than at the split limit.
pub(super) fn contact(
    ctx: &Context,
    spans: [usize; 2],
    rect: &[[Rational; 2]; 4],
    describe: &str,
    dec: &mut Decider,
) -> Result<Option<Refusal>, Refusal> {
    // An offset point at rational parameters is not rational; offsets have no
    // exact contact test, and such a pair refuses at the split limit.
    if ctx.lengths() > 0 {
        return Ok(None);
    }
    let mut middle = Vec::with_capacity(4);
    for interval in rect {
        middle.push(simplest(&interval[0], &interval[1], dec)?);
    }
    let [sa, sb] = spans;
    let a = &ctx.sides[0].spans[sa];
    let b = &ctx.sides[1].spans[sb];
    let pa = [middle[0].clone(), middle[1].clone()];
    let pb = [middle[2].clone(), middle[3].clone()];
    let cost = a
        .homogeneous
        .iter()
        .chain(&b.homogeneous)
        .chain(&a.normal)
        .chain(&b.normal)
        .map(|p| p.terms().len() as u64)
        .sum::<u64>();
    dec.charge_work(cost.saturating_mul(16))?;
    let ha = a.homogeneous.each_ref().map(|p| p.eval(&pa));
    let hb = b.homogeneous.each_ref().map(|p| p.eval(&pb));
    if !(0..3).all(|k| &ha[k] * &hb[3] == &hb[k] * &ha[3]) || ha[3].is_zero() {
        return Ok(None);
    }
    let na = a.normal.each_ref().map(|p| p.eval(&pa));
    let nb = b.normal.each_ref().map(|p| p.eval(&pb));
    let parallel = (0..3).all(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &na[j] * &nb[k] == &na[k] * &nb[j]
    });
    if !parallel {
        return Ok(None);
    }
    // The point in the caller's frame: the rotation is orthogonal, so its
    // transpose undoes it.
    let turned: Vec<Rational> = (0..3).map(|k| &ha[k] / &ha[3]).collect();
    let mut point = Vec::with_capacity(3);
    for axis in 0..3 {
        let mut value = Rational::zero();
        for (row, coordinate) in ctx.rotation.iter().zip(&turned) {
            let term = net::mul(&row[axis], coordinate, dec)?;
            value = net::add(&value, &term, dec)?;
        }
        point.push(value.to_string());
    }
    Ok(Some(Refusal::new(
        "intersect.spline-pair-tangent",
        format!(
            "nonregular contact at ({}, {}, {}), parameters ({}, {}) and ({}, {}); {}",
            point[0], point[1], point[2], pa[0], pa[1], pb[0], pb[1], describe
        ),
    )))
}

pub(super) fn overlap(
    a: &[[Rational; 2]; 4],
    b: &[[Rational; 2]; 4],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for axis in 0..4 {
        if less(&a[axis][1], &b[axis][0], dec)? || less(&b[axis][1], &a[axis][0], dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Polynomials in `(u, v, s, t)` and any offset lengths, with some
/// coordinates fixed, in the remaining variables in axis order.
pub(super) fn restrict(polys: &[MPoly], fixed: &[(usize, Rational)]) -> Vec<MPoly> {
    polys
        .iter()
        .map(|poly| {
            MPoly::from_terms(
                poly.nvars() - fixed.len(),
                poly.terms()
                    .into_iter()
                    .map(|(monomial, coefficient)| {
                        let mut value = coefficient;
                        let mut exponents = Vec::with_capacity(poly.nvars());
                        for (axis, &e) in monomial.exponents().iter().enumerate() {
                            if let Some((_, at)) = fixed.iter().find(|(f, _)| *f == axis) {
                                for _ in 0..e {
                                    value = &value * at;
                                }
                            } else {
                                exponents.push(e);
                            }
                        }
                        (exponents, value)
                    })
                    .collect(),
            )
        })
        .collect()
}

fn charge_terms(polys: &[MPoly], dec: &mut Decider) -> Result<(), Refusal> {
    let terms = polys.iter().map(|p| p.terms().len() as u64).sum::<u64>();
    dec.charge_work(terms.saturating_mul(8))
}

/// The face's meeting point: `Some(None)` when it has none, `None` when the
/// face does not settle at this size. `lengths` encloses each offset side's
/// `L` over the face.
#[allow(clippy::option_option)]
pub(super) fn solve(
    ctx: &Context,
    spans: [usize; 2],
    rect: &[[Rational; 2]; 4],
    lengths: &[[Rational; 2]],
    axis: usize,
    value: &Rational,
    dec: &mut Decider,
) -> Result<Option<Option<Hit>>, Refusal> {
    let equations = ctx.equations(spans);
    let fixed = vec![(axis, value.clone())];
    if !lengths.is_empty() {
        // An offset face: the lifted system decides it, with each length in
        // a box of positive values. Corners and exactly dependent faces are
        // rational-patch cases; here a face that does not isolate splits.
        if lengths.iter().any(|l| l[0].sign() != Sign::Positive) {
            return Ok(None);
        }
        let reduced = restrict(equations, &fixed);
        return certify(reduced, fixed, rect, lengths, spans, dec);
    }
    let others: [usize; 2] = if axis < 2 { [2, 3] } else { [0, 1] };
    // Where a face of one patch's cell lies on a knot line of the other's,
    // as where both patches share a chart, the face system can vanish on a
    // whole corner of the box: some exact combination of its equations is
    // zero there. When that combination factors through the corner with a
    // quotient of one strict sign, every root of the face lies on the
    // corner, and the two other equations decide it.
    for other in others {
        for side in 0..2 {
            let corner = rect[other][side].clone();
            let mut pinned = vec![(axis, value.clone()), (other, corner.clone())];
            pinned.sort_by_key(|(e, _)| *e);
            let reduced = restrict(equations, &pinned);
            charge_terms(&reduced, dec)?;
            let Some(lambda) = dependence(&reduced, dec)? else {
                continue;
            };
            let face = restrict(equations, &fixed);
            let mut combined = MPoly::zero(3);
            for (k, weight) in lambda.iter().enumerate() {
                if !weight.is_zero() {
                    combined = &combined + &face[k].scale(weight);
                }
            }
            if !quotient_signed(&combined, other, &corner, rect, &fixed, dec)? {
                continue;
            }
            let k = (0..3)
                .rev()
                .find(|&k| !lambda[k].is_zero())
                .expect("nonzero");
            let rest: Vec<MPoly> = (0..3)
                .filter(|&j| j != k)
                .map(|j| reduced[j].clone())
                .collect();
            return certify(rest, pinned, rect, &[], spans, dec);
        }
    }
    let reduced = restrict(equations, &fixed);
    if reduce(&reduced, dec)?.1.len() < reduced.len() {
        return dependent_face(&reduced, rect, &fixed, spans, dec);
    }
    certify(reduced, fixed, rect, &[], spans, dec)
}

/// Exact row reduction of the monomial-by-equation matrix of `polys`: the
/// reduced rows and the pivots `(row, column)`. The number of pivots is the
/// rank of the polynomials as vectors of coefficients.
#[allow(clippy::type_complexity)]
fn reduce(
    polys: &[MPoly],
    dec: &mut Decider,
) -> Result<(Vec<Vec<Rational>>, Vec<(usize, usize)>), Refusal> {
    let n = polys.len();
    let mut monomials: BTreeMap<Vec<u32>, usize> = BTreeMap::new();
    let mut rows: Vec<Vec<Rational>> = Vec::new();
    for (k, poly) in polys.iter().enumerate() {
        for (monomial, coefficient) in poly.terms() {
            let row = *monomials
                .entry(monomial.exponents().to_vec())
                .or_insert_with(|| {
                    rows.push(vec![Rational::zero(); n]);
                    rows.len() - 1
                });
            rows[row][k] = coefficient;
        }
    }
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut next = 0;
    for column in 0..n {
        let Some(found) = (next..rows.len()).find(|&r| !rows[r][column].is_zero()) else {
            continue;
        };
        rows.swap(next, found);
        let pivot = rows[next][column].clone();
        for value in rows[next].iter_mut() {
            dec.reserve_exact_arithmetic(value, &pivot, true)?;
            *value = &*value / &pivot;
        }
        let top = rows[next].clone();
        for (r, row) in rows.iter_mut().enumerate() {
            if r == next || row[column].is_zero() {
                continue;
            }
            let factor = row[column].clone();
            for (value, above) in row.iter_mut().zip(&top) {
                let product = net::mul(&factor, above, dec)?;
                *value = net::sub(value, &product, dec)?;
            }
        }
        pivots.push((next, column));
        next += 1;
    }
    Ok((rows, pivots))
}

/// A nonzero exact combination of three polynomials that vanishes, when
/// their coefficient vectors span exactly two dimensions.
fn dependence(polys: &[MPoly], dec: &mut Decider) -> Result<Option<[Rational; 3]>, Refusal> {
    let (rows, pivots) = reduce(polys, dec)?;
    if pivots.len() != 2 {
        return Ok(None);
    }
    let free = (0..3)
        .find(|c| !pivots.iter().any(|(_, p)| p == c))
        .expect("one free column");
    let mut lambda: [Rational; 3] = std::array::from_fn(|_| Rational::zero());
    lambda[free] = Rational::one();
    for (row, column) in pivots {
        lambda[column] = -rows[row][free].clone();
    }
    Ok(Some(lambda))
}

/// Whether one of `polys`, in the free variables of `fixed`, keeps one strict
/// sign over the face box, so that the face holds no root.
fn some_strict(
    polys: &[MPoly],
    rect: &[[Rational; 2]; 4],
    fixed: &[(usize, Rational)],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut centre = Vec::with_capacity(3);
    let mut balls = Vec::with_capacity(3);
    for axis in (0..4).filter(|e| !fixed.iter().any(|(f, _)| f == e)) {
        let [lo, hi] = &rect[axis];
        let c = net::add(lo, hi, dec)?.mul_pow2(-1);
        let r = net::sub(hi, &c, dec)?;
        balls.push(Ball::from_rational(&-r.clone(), 128).hull(&Ball::from_rational(&r, 128)));
        centre.push(c);
    }
    for poly in polys {
        let moved = shifted(poly, &centre);
        charge_terms(std::slice::from_ref(&moved), dec)?;
        if moved.eval_ball(&balls).sign().is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A face whose equations are exactly dependent: its roots form a curve or
/// more, or there are none. With none proved, a branch may run along this
/// cell or patch border, which this stage does not certify; it refuses.
#[allow(clippy::option_option)]
fn dependent_face(
    polys: &[MPoly],
    rect: &[[Rational; 2]; 4],
    fixed: &[(usize, Rational)],
    spans: [usize; 2],
    dec: &mut Decider,
) -> Result<Option<Option<Hit>>, Refusal> {
    if some_strict(polys, rect, fixed, dec)? {
        return Ok(Some(None));
    }
    const NAMES: [&str; 4] = ["u", "v", "s", "t"];
    let face: Vec<String> = (0..4)
        .map(|axis| match fixed.iter().find(|(f, _)| *f == axis) {
            Some((_, value)) => format!("{} = {value}", NAMES[axis]),
            None => format!("{} in [{}, {}]", NAMES[axis], rect[axis][0], rect[axis][1]),
        })
        .collect();
    Err(Refusal::new(
        "intersect.spline-pair-boundary-overlap",
        format!(
            "the meeting may run along a cell or patch border: the equations of the face {} of first span {} against second span {} are exactly dependent, and no branch along a border is certified yet",
            face.join(", "),
            spans[0],
            spans[1]
        ),
    ))
}

/// Whether `face / (x_other - corner)` keeps one strict sign over the face
/// box; `face` is in the face's three free variables.
fn quotient_signed(
    face: &MPoly,
    other: usize,
    corner: &Rational,
    rect: &[[Rational; 2]; 4],
    fixed: &[(usize, Rational)],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let free: Vec<usize> = (0..4)
        .filter(|e| !fixed.iter().any(|(f, _)| f == e))
        .collect();
    let index = free.iter().position(|e| *e == other).expect("free axis");
    let n = face.nvars();
    let shift = &MPoly::var(n, index) + &MPoly::constant(n, corner.clone());
    let moved = face.substitute(index, &shift);
    let mut terms = Vec::new();
    for (monomial, coefficient) in moved.terms() {
        if coefficient.is_zero() {
            continue;
        }
        let mut exponents = monomial.exponents().to_vec();
        if exponents[index] == 0 {
            // The corner value is not a root after all.
            return Ok(false);
        }
        exponents[index] -= 1;
        terms.push((exponents, coefficient));
    }
    let quotient = MPoly::from_terms(n, terms);
    let mut centre = Vec::with_capacity(n);
    let mut radius = Vec::with_capacity(n);
    for (k, &e) in free.iter().enumerate() {
        let (lo, hi) = if k == index {
            (
                net::sub(&rect[e][0], corner, dec)?,
                net::sub(&rect[e][1], corner, dec)?,
            )
        } else {
            (rect[e][0].clone(), rect[e][1].clone())
        };
        let c = net::add(&lo, &hi, dec)?.mul_pow2(-1);
        radius.push(net::sub(&hi, &c, dec)?);
        centre.push(c);
    }
    let quotient = shifted(&quotient, &centre);
    charge_terms(std::slice::from_ref(&quotient), dec)?;
    let balls: Vec<Ball> = radius
        .iter()
        .map(|r| Ball::from_rational(&-r.clone(), 128).hull(&Ball::from_rational(r, 128)))
        .collect();
    Ok(quotient.eval_ball(&balls).sign().is_some())
}

/// Exact root certificate on the face box, then on a tight box round the
/// same root so later sign questions settle cheaply.
#[allow(clippy::option_option)]
fn certify(
    equations: Vec<MPoly>,
    fixed: Vec<(usize, Rational)>,
    rect: &[[Rational; 2]; 4],
    lengths: &[[Rational; 2]],
    spans: [usize; 2],
    dec: &mut Decider,
) -> Result<Option<Option<Hit>>, Refusal> {
    let free: Vec<usize> = (0..4 + lengths.len())
        .filter(|e| !fixed.iter().any(|(f, _)| f == e))
        .collect();
    let range = |e: usize| if e < 4 { &rect[e] } else { &lengths[e - 4] };
    let lo: Vec<Rational> = free.iter().map(|&e| range(e)[0].clone()).collect();
    let hi: Vec<Rational> = free.iter().map(|&e| range(e)[1].clone()).collect();
    let equations = Arc::new(equations);
    let found = isolate(&equations, &lo, &hi, dec)?;
    let Some(found) = found else {
        return Ok(None);
    };
    let Some((centre, mut root)) = found else {
        return Ok(Some(None));
    };
    let system = root.system().clone();
    let better = match tight_box(&system, &root, dec) {
        Some(tight) => root_in(&system, &tight, dec)?,
        None => None,
    };
    if let Some(Some(better)) = better {
        // A root in a sub-box of a box with exactly one root is that root.
        root = better;
    }
    let mut bounds: [[Rational; 2]; 4] = rect.clone();
    for (axis, value) in &fixed {
        bounds[*axis] = [value.clone(), value.clone()];
    }
    let face = bounds.clone();
    // `face` is the box, fixed coordinates pinned, in which the point is the
    // only root: a corner point is the only one on its corner, since every
    // root of that face lies on the corner.
    let mut boxes = Vec::with_capacity(lengths.len());
    for (k, &axis) in free.iter().enumerate() {
        let ends = [
            net::add(&centre[k], &root.root_box().lo()[k], dec)?,
            net::add(&centre[k], &root.root_box().hi()[k], dec)?,
        ];
        if axis < 4 {
            bounds[axis] = ends;
        } else {
            boxes.push(ends);
        }
    }
    Ok(Some(Some(Hit {
        fixed,
        free,
        root,
        centre,
        equations,
        bounds,
        lengths: boxes,
        face,
        spans,
    })))
}

/// The polynomial in variables shifted so that `centre` becomes the origin.
pub(super) fn shifted(poly: &MPoly, centre: &[Rational]) -> MPoly {
    let n = poly.nvars();
    let mut out = poly.clone();
    for (k, c) in centre.iter().enumerate() {
        if !c.is_zero() {
            out = out.substitute(k, &(&MPoly::var(n, k) + &MPoly::constant(n, c.clone())));
        }
    }
    out
}

/// The exact root certificate on the box `[lo, hi]`, for the equations
/// re-centred on the box so that ball evaluation stays tight. Returns the
/// centre with the root in centred variables.
#[allow(clippy::option_option, clippy::type_complexity)]
pub(super) fn isolate(
    equations: &[MPoly],
    lo: &[Rational],
    hi: &[Rational],
    dec: &mut Decider,
) -> Result<Option<Option<(Vec<Rational>, IsolatedRoot)>>, Refusal> {
    let mut centre = Vec::with_capacity(lo.len());
    let mut low = Vec::with_capacity(lo.len());
    let mut high = Vec::with_capacity(lo.len());
    for k in 0..lo.len() {
        let c = net::add(&lo[k], &hi[k], dec)?.mul_pow2(-1);
        low.push(net::sub(&lo[k], &c, dec)?);
        high.push(net::sub(&hi[k], &c, dec)?);
        centre.push(c);
    }
    let moved: Vec<MPoly> = equations.iter().map(|p| shifted(p, &centre)).collect();
    charge_terms(&moved, dec)?;
    let system = Arc::new(PolySystem::new(moved).expect("square face system"));
    let candidate = RootBox::new(low, high).expect("ordered face box");
    Ok(match root_in(&system, &candidate, dec)? {
        None => None,
        Some(None) => Some(None),
        Some(Some(root)) => Some(Some((centre, root))),
    })
}

/// `Some(Some(root))`, `Some(None)` for a proved empty box, `None` when the
/// certificate does not settle.
#[allow(clippy::option_option)]
pub(super) fn root_in(
    system: &Arc<PolySystem>,
    candidate: &RootBox,
    dec: &mut Decider,
) -> Result<Option<Option<IsolatedRoot>>, Refusal> {
    dec.check_deadline()?;
    // Krawczyk can creep, shrinking a box in tiny steps without ever
    // certifying. A per-certificate cap makes such a face inconclusive, so
    // its pair splits. The request's operation allowance pays every unit.
    let mut certifier = Decider::new(catallax_numbers::BudgetPolicy {
        max_work: certify_cap(system),
        ..*dec.policy()
    });
    certifier.set_deadline(dec.deadline());
    match certifier.certify_root_at(system, candidate, CERTIFY_BITS) {
        Outcome::Proven(BoxRoots::None) => Ok(Some(None)),
        Outcome::Proven(BoxRoots::One(root)) => Ok(Some(Some(root))),
        Outcome::Refused(refusal) if refusal.code == "numbers.root-not-certified" => {
            dec.check_deadline()?;
            Ok(None)
        }
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::ProvenCoincident(_) => Ok(None),
    }
}

/// A small dyadic box round a binary64 Newton estimate of the root. It only
/// proposes a candidate; the exact certificate decides.
fn tight_box(system: &PolySystem, root: &IsolatedRoot, dec: &mut Decider) -> Option<RootBox> {
    let n = system.dimension();
    let lo = root.root_box().lo();
    let hi = root.root_box().hi();
    let mut x: Vec<f64> = (0..n)
        .map(|k| ((&lo[k] + &hi[k]).mul_pow2(-1)).to_f64_nearest())
        .collect();
    let terms: Vec<Vec<(Vec<u32>, f64)>> = system
        .equations()
        .iter()
        .map(|p| {
            p.terms()
                .into_iter()
                .map(|(m, c)| (m.exponents().to_vec(), c.to_f64_nearest()))
                .collect()
        })
        .collect();
    let count = terms.iter().map(Vec::len).sum::<usize>() as u64;
    if dec
        .charge_work(count.saturating_mul(n as u64 + 1).saturating_mul(12))
        .is_err()
    {
        return None;
    }
    let mut step_size = f64::INFINITY;
    for _ in 0..12 {
        let mut value = vec![0.0; n];
        let mut jacobian = vec![vec![0.0; n]; n];
        for (row, poly) in terms.iter().enumerate() {
            for (exponents, c) in poly {
                let mut term = *c;
                for (k, &e) in exponents.iter().enumerate() {
                    for _ in 0..e {
                        term *= x[k];
                    }
                }
                value[row] += term;
                for (j, &e) in exponents.iter().enumerate() {
                    if e == 0 {
                        continue;
                    }
                    let mut partial = *c * f64::from(e);
                    for (k, &ek) in exponents.iter().enumerate() {
                        let power = if k == j { ek - 1 } else { ek };
                        for _ in 0..power {
                            partial *= x[k];
                        }
                    }
                    jacobian[row][j] += partial;
                }
            }
        }
        let step = solve_dense(jacobian, value)?;
        step_size = step.iter().fold(0.0f64, |m, s| m.max(s.abs()));
        for (xk, sk) in x.iter_mut().zip(&step) {
            *xk -= sk;
        }
        if !x.iter().all(|v| v.is_finite()) {
            return None;
        }
        if step_size == 0.0 {
            break;
        }
    }
    let scale = x.iter().fold(1.0f64, |m, v| m.max(v.abs()));
    let radius = (step_size * 1024.0).max(scale * f64::EPSILON * 1024.0);
    let mut new_lo = Vec::with_capacity(n);
    let mut new_hi = Vec::with_capacity(n);
    for k in 0..n {
        let a = Rational::from_f64_exact(x[k] - radius)?;
        let b = Rational::from_f64_exact(x[k] + radius)?;
        // Stay inside the certified box: a sub-box keeps the uniqueness.
        let a = a.max(lo[k].clone());
        let b = b.min(hi[k].clone());
        if a > b {
            return None;
        }
        new_lo.push(a);
        new_hi.push(b);
    }
    RootBox::new(new_lo, new_hi)
}

fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))?;
        if a[pivot][k] == 0.0 {
            return None;
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        let (top, rest) = a.split_at_mut(k + 1);
        let pivot_row = &top[k];
        for (offset, row) in rest.iter_mut().enumerate() {
            let f = row[k] / pivot_row[k];
            for (value, above) in row[k..n].iter_mut().zip(&pivot_row[k..n]) {
                *value -= f * above;
            }
            b[k + 1 + offset] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let mut s = b[k];
        for j in k + 1..n {
            s -= a[k][j] * x[j];
        }
        x[k] = s / a[k][k];
    }
    Some(x)
}

/// Whether the branch through a border point enters the pair there, leaves
/// it there, or only touches its corner. `None` when signs do not settle.
#[allow(clippy::option_option)]
pub(super) fn classify(
    ctx: &Context,
    rect: &[[Rational; 2]; 4],
    hit: Hit,
    dec: &mut Decider,
) -> Result<Option<Option<Endpoint>>, Refusal> {
    let mut active = hit.fixed.clone();
    for (k, &axis) in hit.free.iter().enumerate() {
        // Offset lengths are not borders of the unit.
        if axis >= 4 {
            continue;
        }
        for (side, bound) in rect[axis].iter().enumerate() {
            let inside = if side == 0 {
                less(bound, &hit.bounds[axis][0], dec)?
            } else {
                less(&hit.bounds[axis][1], bound, dec)?
            };
            if inside {
                continue;
            }
            let offset = &hit.root.coordinate(k)
                - &Real::from_rational(net::sub(bound, &hit.centre[k], dec)?);
            match dec.sign(&offset) {
                Outcome::Proven(StrictSign::Positive) if side == 0 => (),
                Outcome::Proven(StrictSign::Negative) if side == 1 => (),
                Outcome::ProvenCoincident(_) => active.push((axis, bound.clone())),
                Outcome::Refused(refusal)
                    if refusal.code == "numbers.operation-budget"
                        || refusal.code == "numbers.deadline" =>
                {
                    return Err(refusal);
                }
                _ => return Ok(None),
            }
        }
    }
    let Some(signs) = tangent_signs(ctx, &hit, dec)? else {
        return Ok(None);
    };
    let inward = |(axis, value): &(usize, Rational)| {
        if value == &rect[*axis][0] {
            Sign::Positive
        } else {
            Sign::Negative
        }
    };
    let mut forward = true;
    let mut backward = true;
    for plane in &active {
        let sign = signs[plane.0];
        if sign == Sign::Zero {
            return Ok(None);
        }
        let flipped = if sign == Sign::Positive {
            Sign::Negative
        } else {
            Sign::Positive
        };
        forward &= sign == inward(plane);
        backward &= flipped == inward(plane);
    }
    if !(forward || backward) {
        return Ok(Some(None));
    }
    Ok(Some(Some(Endpoint {
        hit,
        active,
        start: forward,
    })))
}

/// Signs of the branch tangent `(du, dv, ds, dt)` at the point, in the
/// direction of first normal cross second normal: `(-y, x, y', -x')`, with
/// `(x, y)` the field on the first patch and `(x', y')` the field on the
/// second. `Zero` where a sign does not settle; `None` when none settles.
fn tangent_signs(
    ctx: &Context,
    hit: &Hit,
    dec: &mut Decider,
) -> Result<Option<[Sign; 4]>, Refusal> {
    let a = &ctx.sides[0].spans[hit.spans[0]];
    let b = &ctx.sides[1].spans[hit.spans[1]];
    if ctx.lengths() > 0 {
        return offset_tangent_signs(a, b, hit, dec);
    }
    let all: Vec<&MPoly> = a
        .tangents
        .iter()
        .flatten()
        .chain(&a.normal)
        .chain(b.tangents.iter().flatten())
        .chain(&b.normal)
        .collect();
    let terms = all.iter().map(|p| p.terms().len() as u64).sum::<u64>();
    dec.charge_work(terms.saturating_mul(16))?;
    let balls: Vec<Ball> = hit
        .bounds
        .iter()
        .map(|[lo, hi]| Ball::from_rational(lo, 192).hull(&Ball::from_rational(hi, 192)))
        .collect();
    let pa = [balls[0].clone(), balls[1].clone()];
    let pb = [balls[2].clone(), balls[3].clone()];
    let eval = |polys: &[MPoly; 3], at: &[Ball; 2]| -> [Ball; 3] {
        std::array::from_fn(|i| polys[i].eval_ball(at))
    };
    let dot = |x: &[Ball; 3], y: &[Ball; 3]| -> Ball {
        &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
    };
    let na = eval(&a.normal, &pa);
    let nb = eval(&b.normal, &pb);
    let x = dot(&eval(&a.tangents[0], &pa), &nb);
    let y = dot(&eval(&a.tangents[1], &pa), &nb);
    let x2 = dot(&eval(&b.tangents[0], &pb), &na);
    let y2 = dot(&eval(&b.tangents[1], &pb), &na);
    let sign = |ball: &Ball, flip: bool| match ball.sign() {
        Some(StrictSign::Positive) if !flip => Sign::Positive,
        Some(StrictSign::Negative) if flip => Sign::Positive,
        Some(_) => Sign::Negative,
        None => Sign::Zero,
    };
    let signs = [
        sign(&y, true),
        sign(&x, false),
        sign(&y2, false),
        sign(&x2, true),
    ];
    if signs.iter().all(|s| *s == Sign::Zero) {
        return Ok(None);
    }
    Ok(Some(signs))
}

/// The same signs where a side is an offset: its tangents times `W^2 L^3`
/// are `L p1 + p0`, with `L` in the root's box, and its normal is a positive
/// multiple of `N` on the regular cells that hold the point.
fn offset_tangent_signs(
    a: &SideSpan,
    b: &SideSpan,
    hit: &Hit,
    dec: &mut Decider,
) -> Result<Option<[Sign; 4]>, Refusal> {
    let balls: Vec<Ball> = hit
        .bounds
        .iter()
        .map(|[lo, hi]| Ball::from_rational(lo, 192).hull(&Ball::from_rational(hi, 192)))
        .collect();
    let mut lengths = hit
        .lengths
        .iter()
        .map(|[lo, hi]| Ball::from_rational(lo, 192).hull(&Ball::from_rational(hi, 192)));
    let mut tangents = Vec::with_capacity(2);
    let mut normals = Vec::with_capacity(2);
    for (span, at) in [
        (a, [balls[0].clone(), balls[1].clone()]),
        (b, [balls[2].clone(), balls[3].clone()]),
    ] {
        let polys: Vec<&MPoly> = match &span.offset {
            Some(offset) => offset
                .p1
                .iter()
                .chain(&offset.p0)
                .flatten()
                .chain(&span.normal)
                .collect(),
            None => span.tangents.iter().flatten().chain(&span.normal).collect(),
        };
        let terms = polys.iter().map(|p| p.terms().len() as u64).sum::<u64>();
        dec.charge_work(terms.saturating_mul(16))?;
        let eval =
            |polys: &[MPoly; 3]| -> [Ball; 3] { std::array::from_fn(|i| polys[i].eval_ball(&at)) };
        let pair: [[Ball; 3]; 2] = match &span.offset {
            Some(offset) => {
                let length = lengths.next().expect("a length per offset side");
                std::array::from_fn(|axis| {
                    let [p1, p0] = [eval(&offset.p1[axis]), eval(&offset.p0[axis])];
                    std::array::from_fn(|i| &(&length * &p1[i]) + &p0[i])
                })
            }
            None => [eval(&span.tangents[0]), eval(&span.tangents[1])],
        };
        tangents.push(pair);
        normals.push(eval(&span.normal));
    }
    let dot = |x: &[Ball; 3], y: &[Ball; 3]| -> Ball {
        &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
    };
    let x = dot(&tangents[0][0], &normals[1]);
    let y = dot(&tangents[0][1], &normals[1]);
    let x2 = dot(&tangents[1][0], &normals[0]);
    let y2 = dot(&tangents[1][1], &normals[0]);
    let sign = |ball: &Ball, flip: bool| match ball.sign() {
        Some(StrictSign::Positive) if !flip => Sign::Positive,
        Some(StrictSign::Negative) if flip => Sign::Positive,
        Some(_) => Sign::Negative,
        None => Sign::Zero,
    };
    let signs = [
        sign(&y, true),
        sign(&x, false),
        sign(&y2, false),
        sign(&x2, true),
    ];
    if signs.iter().all(|s| *s == Sign::Zero) {
        return Ok(None);
    }
    Ok(Some(signs))
}
