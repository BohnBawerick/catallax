//! The complete meeting of two retained rational patches: M10's
//! spline-with-spline stage, on the 2026 field method of Yang and Jia with an
//! exact proof layer.
//!
//! The search splits the first patch's parameter domain into cells. Each cell
//! keeps a cover: the second patch's cells, split to its size, that exclusion
//! has not removed. Exclusion is exact: the boxes round the two control nets
//! are disjoint, or their slabs along either cell's normal direction are, or
//! the boxes touch only on a side that one patch reaches on its cell border
//! alone. A first-patch cell against the cover cells of one second-patch knot
//! span is a unit, accepted only with an exact certificate:
//!
//! - the field of Yang and Jia, `(Su . Nb, Sv . Nb)` on the first patch's
//!   parameters, or the same field on the second patch's parameters, keeps
//!   one component of one strict sign over the whole unit. Then one parameter
//!   of that patch is strictly monotone along every branch in the unit, so the
//!   unit holds no closed loop, no isolated contact and no critical point, and
//!   the winding number of the field round every closed curve in it is zero;
//! - the normal cone: one coordinate of `Na x Nb` keeps one strict sign, so
//!   that spatial coordinate is strictly monotone along every branch.
//!
//! The unit must also show zero border points, or exactly two that the branch
//! crosses transversally, one entering and one leaving, which bound its one
//! arc there. Every border point is an exact root of the polynomial system on
//! its face, isolated by an exact Krawczyk certificate. A cell with any other
//! unit splits into thirds; one still unresolved at the split limit refuses
//! and names its cells. Arcs join at their shared border points into open and
//! closed branches.
//!
//! Binary64 arithmetic only chooses: a rotation, a candidate box round a root
//! that the exact certificate then checks, and which exact slab tests to skip
//! because they would clearly fail. Every accepted fact is exact.

use super::*;
use catallax_numbers::{BoxRoots, IsolatedRoot, MPoly, PolySystem, RootBox};

mod face;
mod join;
mod net;
mod run;
mod unit;

use net::Cell;

/// Most first-patch cells one search examines before it refuses.
const MAX_CELLS: usize = 1 << 16;

/// The certificate that accepted a cell pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PairCertificate {
    /// The 2026 field on the first patch's parameters stays in the open
    /// half-plane `d . (x, y) > 0` for this direction `d`.
    FieldFirst([i8; 2]),
    /// The same field on the second patch's parameters.
    FieldSecond([i8; 2]),
    /// This coordinate of the cross product of the two normals keeps one
    /// strict sign.
    NormalCone(u8),
}

/// An exact meeting point on the border of a cell pair: parameters `u, v`
/// on the first patch and `s, t` on the second. Fixed coordinates are exact
/// rationals; the rest are coordinates of an exactly isolated system root.
#[derive(Clone, Debug)]
pub struct PairPoint {
    parameters: [Real; 4],
    enclosure: [[Rational; 2]; 4],
    first_span: usize,
}

impl PairPoint {
    /// Exact parameters `[u, v, s, t]`.
    pub fn parameters(&self) -> &[Real; 4] {
        &self.parameters
    }
    /// A rational box round the parameters, as certified by the search.
    pub fn enclosure(&self) -> &[[Rational; 2]; 4] {
        &self.enclosure
    }
}

/// One certified arc of a branch inside one cell pair.
#[derive(Clone, Debug)]
pub struct PairArc {
    /// Knot span of the first patch holding the first cell.
    pub first_span: usize,
    /// First patch's parameter cell.
    pub first_cell: [[Rational; 2]; 2],
    /// Knot span of the second patch holding the second cell.
    pub second_span: usize,
    /// Second patch's parameter cell.
    pub second_cell: [[Rational; 2]; 2],
    /// The certificate that proved the pair holds exactly this arc.
    pub certificate: PairCertificate,
    /// The arc's ends, in the direction of first normal cross second normal.
    pub ends: [PairPoint; 2],
}

/// One connected component of the meeting, as the ordered arcs that cover it.
#[derive(Clone, Debug)]
pub struct PatchPairBranch {
    arcs: Arc<[PairArc]>,
    closed: bool,
}

impl PatchPairBranch {
    /// The arcs in traversal order, in the direction of first normal cross
    /// second normal.
    pub fn arcs(&self) -> &[PairArc] {
        &self.arcs
    }
    /// Whether the branch is a closed loop.
    pub fn is_closed(&self) -> bool {
        self.closed
    }
    /// The two ends of an open branch, on a patch border; for a loop, its
    /// first join point twice.
    pub fn ends(&self) -> [&PairPoint; 2] {
        [
            &self.arcs[0].ends[0],
            &self.arcs[self.arcs.len() - 1].ends[1],
        ]
    }
}

/// Counts and work of one complete search, for the record.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PatchPairSearch {
    /// First-patch cells examined.
    pub cells: usize,
    /// Second-patch cover cells kept against them, summed.
    pub cover: usize,
    /// Cells whose whole cover was excluded.
    pub excluded: usize,
    /// Accepted units (a cell against one second-patch knot span), by the
    /// certificate that accepted them.
    pub field_first: usize,
    /// See `field_first`.
    pub field_second: usize,
    /// See `field_first`.
    pub normal_cone: usize,
    /// Accepted units holding one arc.
    pub arcs: usize,
    /// Cells split because no certificate held or the border did not settle.
    pub splits: usize,
    /// Face systems sent to the exact root certificate.
    pub face_systems: usize,
    /// Most splits of one first-patch cell.
    pub depth: u32,
    /// Rounds of the level-by-level search.
    pub levels: usize,
    /// Deterministic decision work of the whole search.
    pub work: u64,
}

/// The complete checked meeting of two rational patches.
#[derive(Clone, Debug)]
pub struct PatchPairIntersection {
    first: RationalPatch,
    branches: Vec<PatchPairBranch>,
    search: PatchPairSearch,
}

impl PatchPairIntersection {
    /// Every branch of the meeting.
    pub fn branches(&self) -> &[PatchPairBranch] {
        &self.branches
    }
    /// Counts and work.
    pub fn search(&self) -> &PatchPairSearch {
        &self.search
    }
    /// The exact meeting point, from the first patch.
    pub fn point(&self, at: &PairPoint) -> Point3 {
        let [u, v, _, _] = &at.parameters;
        self.first.data.spans[at.first_span].point(&[u.clone(), v.clone()])
    }
}

/// One patch's spans, in the forms the search reads.
struct Side {
    domain: Domain,
    spans: Vec<SideSpan>,
}

struct SideSpan {
    /// Global homogeneous coordinates.
    homogeneous: [MPoly; 4],
    /// Global homogeneous tangents `W P_u - W_u P` and `W P_v - W_v P`.
    tangents: [[MPoly; 3]; 2],
    /// Their cross product, a positive multiple of the normal.
    normal: [MPoly; 3],
    root: Arc<Cell>,
}

impl Side {
    fn new(
        patch: &RationalPatch,
        homogeneous: Vec<[MPoly; 4]>,
        rotation: &[[Rational; 3]; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut spans = Vec::with_capacity(homogeneous.len());
        for (index, (span, homogeneous)) in patch.data.spans.iter().zip(homogeneous).enumerate() {
            dec.check_deadline()?;
            let homogeneous = rotated(&homogeneous, rotation, dec)?;
            let mut tangents = Vec::with_capacity(2);
            for axis in 0..2 {
                let dw = homogeneous[3].partial(axis);
                let mut parts = Vec::with_capacity(3);
                for coordinate in &homogeneous[..3] {
                    let left = product(&coordinate.partial(axis), &homogeneous[3], dec)?;
                    let right = product(coordinate, &dw, dec)?;
                    parts.push(difference(&left, &right, dec)?);
                }
                tangents.push(parts.try_into().expect("three components"));
            }
            let tangents: [[MPoly; 3]; 2] = tangents.try_into().expect("two tangents");
            let [tu, tv] = &tangents;
            let mut normal = Vec::with_capacity(3);
            for i in 0..3 {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                let left = product(&tu[j], &tv[k], dec)?;
                let right = product(&tu[k], &tv[j], dec)?;
                normal.push(difference(&left, &right, dec)?);
            }
            let root = Arc::new(Cell::span(index, &span.domain, &homogeneous, dec)?);
            spans.push(SideSpan {
                homogeneous,
                tangents,
                normal: normal.try_into().expect("three components"),
                root,
            });
        }
        Ok(Self {
            domain: patch.data.domain.clone(),
            spans,
        })
    }
}

/// Each span's homogeneous coordinates as exact rational polynomials.
fn rational_spans(patch: &RationalPatch) -> Result<Vec<[MPoly; 4]>, Refusal> {
    patch
        .data
        .spans
        .iter()
        .map(|span| {
            let coordinates: Option<Vec<MPoly>> =
                span.coordinates.iter().map(RealPoly::as_rational).collect();
            coordinates
                .map(|c| c.try_into().expect("four coordinates"))
                .ok_or_else(|| {
                    Refusal::new(
                        "intersect.spline-pair-irrational",
                        "the spline pair search needs rational patch coefficients",
                    )
                })
        })
        .collect()
}

/// Refuses a fill of a charged shared cell on a run-ahead worker. Only the
/// request's own thread fills these, in a fixed order, so a worker never
/// should; if one does, `shared_caches` drops its record and the request's
/// thread proves that cell itself.
fn writer_only() -> Result<(), Refusal> {
    if catallax_numbers::shared_caches() {
        Ok(())
    } else {
        Err(Refusal::new(
            "intersect.spline-pair-worker-fill",
            "a run-ahead worker met a shared cell only the request's own thread fills",
        ))
    }
}

/// Exact sum of limbs over a polynomial's coefficients.
fn size(poly: &MPoly) -> u64 {
    poly.terms().iter().map(|(_, c)| net::limbs(c)).sum()
}

/// `a * b`, charged as one exact product per pair of coefficients.
fn product(a: &MPoly, b: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    dec.charge_work(size(a).saturating_mul(size(b)))?;
    Ok(a * b)
}

/// `a - b`, charged as one exact sum per coefficient.
fn difference(a: &MPoly, b: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    dec.charge_work(size(a).saturating_add(size(b)))?;
    Ok(a - b)
}

/// Everything a cell pair's proof reads, shared by every thread.
struct Context {
    sides: [Side; 2],
    /// The exact rotation applied to both patches; reports undo it.
    rotation: [[Rational; 3]; 3],
    /// `A_k B_w - B_k A_w` in `(u, v, s, t)`, per span pair, built once.
    equations: Vec<std::sync::OnceLock<[MPoly; 3]>>,
}

impl Context {
    fn index(&self, spans: [usize; 2]) -> usize {
        spans[0] * self.sides[1].spans.len() + spans[1]
    }

    /// Builds a span pair's equations. Only the request's own thread calls
    /// this, in a fixed order, so the charge does not depend on threads.
    fn ensure_equations(&self, spans: [usize; 2], dec: &mut Decider) -> Result<(), Refusal> {
        let slot = &self.equations[self.index(spans)];
        if slot.get().is_some() {
            return Ok(());
        }
        writer_only()?;
        let a = &self.sides[0].spans[spans[0]].homogeneous;
        let b = &self.sides[1].spans[spans[1]].homogeneous;
        let a = a.each_ref().map(|p| embed(p, [0, 1], 4));
        let b = b.each_ref().map(|p| embed(p, [2, 3], 4));
        let mut equations = Vec::with_capacity(3);
        for k in 0..3 {
            let left = product(&a[k], &b[3], dec)?;
            let right = product(&b[k], &a[3], dec)?;
            equations.push(difference(&left, &right, dec)?);
        }
        let _ = slot.set(equations.try_into().expect("three equations"));
        Ok(())
    }

    /// A span pair's equations; the writer built them before the level.
    fn equations(&self, spans: [usize; 2]) -> &[MPoly; 3] {
        self.equations[self.index(spans)]
            .get()
            .expect("the writer builds a span pair's equations before proving it")
    }

    /// Whether `value` on `axis` is a border of the physical domain.
    fn physical(&self, axis: usize, value: &Rational) -> bool {
        let domain = &self.sides[axis / 2].domain[axis % 2];
        &domain[0] == value || &domain[1] == value
    }
}

/// The exact rotation of a small integer quaternion that best brings the
/// first patch's mean normal over its span corners to the `z` axis and its
/// mean `u` tangent to the `x` axis. A binary64 search, charged a fixed
/// amount, only chooses the quaternion; the rotation is exact, and any
/// rotation keeps the meeting. The identity stays unless another quaternion
/// is clearly better.
fn upright(
    patch: &RationalPatch,
    spans: &[[MPoly; 4]],
    dec: &mut Decider,
) -> Result<[[Rational; 3]; 3], Refusal> {
    let terms = spans
        .iter()
        .flatten()
        .map(|p| p.terms().len() as u64)
        .sum::<u64>();
    dec.charge_work(terms.saturating_mul(16).saturating_add(4096))?;
    // Mean directions, each span's sum over its corners scaled to unit size.
    let mut means = [[0.0f64; 3]; 2];
    for (span, homogeneous) in patch.data.spans.iter().zip(spans) {
        let partials = [0, 1].map(|axis| homogeneous.each_ref().map(|p| p.partial(axis)));
        let mut totals = [[0.0f64; 3]; 2];
        for u in &span.domain[0] {
            for v in &span.domain[1] {
                let at = [u.to_f64_nearest(), v.to_f64_nearest()];
                let h = homogeneous.each_ref().map(|p| approximate(p, at));
                let [tu, tv] = partials.each_ref().map(|d| {
                    let d = d.each_ref().map(|p| approximate(p, at));
                    [0, 1, 2].map(|i| d[i] * h[3] - h[i] * d[3])
                });
                let normal = [0, 1, 2].map(|i| {
                    let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                    tu[j] * tv[k] - tu[k] * tv[j]
                });
                for (total, value) in totals.iter_mut().zip([normal, tu]) {
                    for axis in 0..3 {
                        total[axis] += value[axis];
                    }
                }
            }
        }
        for (mean, total) in means.iter_mut().zip(totals) {
            let size = total.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            if size > 0.0 && size.is_finite() {
                for (sum, value) in mean.iter_mut().zip(total) {
                    *sum += value / size;
                }
            }
        }
    }
    let unit = |v: [f64; 3]| {
        let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        (length > 0.0 && length.is_finite()).then(|| v.map(|x| x / length))
    };
    let mut best = [1i64, 0, 0, 0];
    if let (Some(n), Some(t)) = (unit(means[0]), unit(means[1])) {
        let score = |r: &[[f64; 3]; 3]| {
            let up = r[2][0] * n[0] + r[2][1] * n[1] + r[2][2] * n[2];
            let along = r[0][0] * t[0] + r[0][1] * t[1] + r[0][2] * t[2];
            up + along
        };
        let mut top = score(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        for w in 0..=8i64 {
            for x in -8..=8i64 {
                for y in -8..=8i64 {
                    for z in -8..=8i64 {
                        let q = [w, x, y, z];
                        let size = q.iter().map(|c| c * c).sum::<i64>();
                        if size == 0 {
                            continue;
                        }
                        let r = quaternion(q).map(|row| row.map(|c| c as f64 / size as f64));
                        let here = score(&r);
                        if here > top + 1e-3 {
                            top = here;
                            best = q;
                        }
                    }
                }
            }
        }
    }
    let size = best.iter().map(|c| c * c).sum::<i64>();
    Ok(quaternion(best)
        .map(|row| row.map(|c| Rational::from_ratio(c, size).expect("nonzero quaternion"))))
}

/// A binary64 value of a polynomial in two variables.
fn approximate(poly: &MPoly, at: [f64; 2]) -> f64 {
    poly.terms()
        .iter()
        .map(|(monomial, coefficient)| {
            let e = monomial.exponents();
            coefficient.to_f64_nearest() * at[0].powi(e[0] as i32) * at[1].powi(e[1] as i32)
        })
        .sum()
}

/// The integer rotation matrix of a quaternion, scaled by its squared norm.
fn quaternion([w, x, y, z]: [i64; 4]) -> [[i64; 3]; 3] {
    [
        [
            w * w + x * x - y * y - z * z,
            2 * (x * y - w * z),
            2 * (x * z + w * y),
        ],
        [
            2 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2 * (y * z - w * x),
        ],
        [
            2 * (x * z - w * y),
            2 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ]
}

/// Homogeneous coordinates turned by an exact rotation; the weight stays.
fn rotated(
    homogeneous: &[MPoly; 4],
    rotation: &[[Rational; 3]; 3],
    dec: &mut Decider,
) -> Result<[MPoly; 4], Refusal> {
    let identity =
        (0..3).all(|i| (0..3).all(|j| rotation[i][j] == Rational::from_i64(i64::from(i == j))));
    if identity {
        return Ok(homogeneous.clone());
    }
    let n = homogeneous[0].nvars();
    let mut out = Vec::with_capacity(4);
    for row in rotation {
        let mut sum = MPoly::zero(n);
        for (entry, coordinate) in row.iter().zip(homogeneous) {
            if entry.is_zero() {
                continue;
            }
            let scaled = product(coordinate, &MPoly::constant(n, entry.clone()), dec)?;
            dec.charge_work(size(&sum).saturating_add(size(&scaled)))?;
            sum = &sum + &scaled;
        }
        out.push(sum);
    }
    out.push(homogeneous[3].clone());
    Ok(out.try_into().expect("four coordinates"))
}

/// Re-indexes a polynomial's variables into a larger space.
fn embed(poly: &MPoly, targets: [usize; 2], nvars: usize) -> MPoly {
    MPoly::from_terms(
        nvars,
        poly.terms()
            .into_iter()
            .map(|(monomial, coefficient)| {
                let mut exponents = vec![0; nvars];
                for (k, &target) in targets.iter().enumerate() {
                    exponents[target] = monomial.exponents()[k];
                }
                (exponents, coefficient)
            })
            .collect(),
    )
}

impl RationalPatch {
    /// The complete exact meeting with another rational patch.
    ///
    /// Every branch is returned as an ordered cover of cell pairs, each
    /// accepted by an exact certificate and bounded by exactly isolated
    /// border points. A cell pair that no certificate settles within the
    /// split and work allowances refuses with
    /// `intersect.spline-pair-unresolved`, or `intersect.spline-pair-tangent`
    /// when an exact tangent contact sits at its centre; never a near-enough
    /// answer.
    pub fn patch_meeting(
        &self,
        other: &RationalPatch,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        // With an operation allowance active, its counter includes the capped
        // root certificates, which run on deciders of their own.
        let budget = catallax_numbers::operation_budget_active();
        let used = || catallax_numbers::current_operation_work().used;
        let start = if budget { used() } else { dec.stats().work };
        // A rotation changes no parameter of the meeting, and boxes along
        // the axes follow tilted surfaces once they are turned upright.
        let spans = [rational_spans(self)?, rational_spans(other)?];
        let rotation = upright(self, &spans[0], dec)?;
        let [first, second] = spans;
        let sides = [
            Side::new(self, first, &rotation, dec)?,
            Side::new(other, second, &rotation, dec)?,
        ];
        let pairs = sides[0].spans.len() * sides[1].spans.len();
        let context = Context {
            sides,
            rotation,
            equations: (0..pairs).map(|_| std::sync::OnceLock::new()).collect(),
        };
        let (arcs, mut search) = run::search(&context, dec)?;
        let branches = join::branches(&context, arcs, dec)?;
        search.work = if budget { used() } else { dec.stats().work } - start;
        Ok(PatchPairIntersection {
            first: self.clone(),
            branches,
            search,
        })
    }
}
