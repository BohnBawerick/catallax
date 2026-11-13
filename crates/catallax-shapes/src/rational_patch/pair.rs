//! The complete meeting of two retained rational patches: M10's
//! spline-with-spline stage, on the 2026 field method of Yang and Jia with an
//! exact proof layer.
//!
//! The meeting is that of the closed patches, borders included. The search
//! splits the first patch's parameter domain into cells. Each cell keeps a
//! cover: the second patch's cells, split to its size, that exclusion has not
//! removed. Exclusion is exact: the boxes round the two control nets are
//! disjoint, or their slabs along either cell's normal direction are. Boxes
//! that touch exactly drop a pair only when every meeting point in it lies on
//! an internal high border of one of its cells, which the neighbouring cell
//! owns; a point on a physical border is never dropped this way. A
//! first-patch cell against the cover cells of one second-patch knot span is
//! a unit, accepted only with an exact certificate:
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
//! its face, isolated by an exact Krawczyk certificate. A face whose equations
//! are exactly dependent, so that a branch may run along that border, refuses
//! unless one equation proves it empty. A cell with any other unit splits into
//! thirds; one still unresolved at the split limit refuses and names its
//! cells. Arcs join at their shared border points into open and closed
//! branches.
//!
//! Binary64 arithmetic only chooses: a rotation, a candidate box round a root
//! that the exact certificate then checks, and which exact slab tests to skip
//! because they would clearly fail. Every accepted fact is exact.
//!
//! Either side may instead be an exact offset `P + d N/|N|` of a rational
//! patch, with `N` its normal and `d` a signed rational distance. The offset
//! is not a rational surface, so the search carries the length `L = |N|` as
//! one more unknown with `L^2 = N . N` and `L > 0`, and keeps the support's
//! own parameter rectangle, borders included. Its position, its tangents
//! times `W^2 L^3` and its homogeneous coordinates are all linear in `L` with
//! polynomial coefficients, so cells bound them from exact Bernstein nets and
//! a rational enclosure of `L`. Exclusion widens each control box by the
//! box round `d N/|N|` over the cell. A unit is certified only where the
//! offset is proved regular with the support's orientation: the triple
//! product of its two tangents with `N` is positive over every cell. The
//! offset's normal is then a positive multiple of `N`, so the field and
//! normal-cone certificates read `N`. Face systems gain one length equation
//! per offset side. The exact tangent-contact test and the exact touch rule
//! apply to rational patches only; an offset pair with a tangent contact
//! refuses at the split limit instead.

use super::*;
use catallax_numbers::{BoxRoots, IsolatedRoot, MPoly, PolySystem, RootBox};

mod cover;
use cover::{RetainedSource, Source};
mod face;
mod join;
mod net;
mod original;
mod prepared;
mod reuse;
mod run;
mod trace;
mod unit;
pub use trace::{PairTracePiece, PatchPairCurve};

use net::Cell;
pub(crate) use net::PolynomialBounds;

// Positive spline weights keep every curve point in the convex hull
// of its Cartesian controls. Native curves use the same strict closed-box
// exclusion certificate as the patch-pair engine, without claiming a meeting
// or identity when the coarse hulls overlap.
pub(crate) fn curve_hulls_disjoint(
    first: &BSplineCurve,
    second: &BSplineCurve,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    fn control_box(
        curve: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<Option<net::Bounds3>, Refusal> {
        let mut controls = curve.controls();
        let Some(first) = controls
            .next()
            .expect("nonempty finite spline controls")
            .as_rational()
        else {
            return Ok(None);
        };
        let mut bounds = [
            [first.x.clone(), first.x],
            [first.y.clone(), first.y],
            [first.z.clone(), first.z],
        ];
        for point in controls {
            let Some(point) = point.as_rational() else {
                return Ok(None);
            };
            for (axis, value) in [point.x, point.y, point.z].into_iter().enumerate() {
                if net::less(&value, &bounds[axis][0], dec)? {
                    bounds[axis][0] = value.clone();
                }
                if net::less(&bounds[axis][1], &value, dec)? {
                    bounds[axis][1] = value;
                }
            }
        }
        Ok(Some(bounds))
    }
    let Some(first) = control_box(first, dec)? else {
        return Ok(false);
    };
    let Some(second) = control_box(second, dec)? else {
        return Ok(false);
    };
    net::disjoint(&first, &second, dec)
}

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
    first_chart: usize,
    first_span: usize,
    reconstruction: Arc<face::CoordinateWitness>,
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
    /// First retained source unit: a chart, or an original tensor span.
    pub first_chart: usize,
    /// Knot span within that unit; zero for original-span meetings.
    pub first_span: usize,
    /// First patch's parameter cell.
    pub first_cell: [[Rational; 2]; 2],
    /// Second retained source unit: a chart, or an original tensor span.
    pub second_chart: usize,
    /// Knot span within that unit; zero for original-span meetings.
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

/// The complete checked meeting of two rational patches, or of their exact
/// offsets.
#[derive(Clone, Debug)]
pub struct PatchPairIntersection {
    first: RetainedSource,
    /// The first side's offset distance, when it is an offset.
    distance: Option<Rational>,
    second: RetainedSource,
    second_distance: Option<Rational>,
    rotation: [[Rational; 3]; 3],
    branches: Vec<PatchPairBranch>,
    search: PatchPairSearch,
}

/// An exact offset of a rational patch: the surface `P + d N / |N|`, where
/// `N = P_u x P_v` and `d` is a signed rational distance. It is not a
/// rational surface in general; [`OffsetPatch::offset_meeting`] meets two
/// of them exactly over their supports' closed parameter domains.
#[derive(Clone, Debug)]
pub struct OffsetPatch {
    support: RationalPatch,
    distance: Rational,
}

impl OffsetPatch {
    /// The rational patch this offsets.
    pub fn support(&self) -> &RationalPatch {
        &self.support
    }
    /// The signed distance along `N / |N|`.
    pub fn distance(&self) -> &Rational {
        &self.distance
    }

    /// The complete exact meeting of two offsets, over both supports' closed
    /// parameter domains, borders included, by the same search and
    /// certificates as [`RationalPatch::patch_meeting`]. Parameters in the
    /// answer are the supports' parameters. A zero distance is the support
    /// itself.
    ///
    /// Besides that method's refusals, a cell whose offset is not proved
    /// regular with the support's orientation, such as one reaching a focal
    /// point, is never accepted; it splits and, at the split limit, refuses
    /// with `intersect.spline-pair-unresolved`.
    pub fn offset_meeting(
        &self,
        other: &OffsetPatch,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        let distance = |d: &Rational| (!d.is_zero()).then(|| d.clone());
        meeting(
            [Source::Patch(&self.support), Source::Patch(&other.support)],
            [distance(&self.distance), distance(&other.distance)],
            dec,
        )
    }
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
    /// The exact meeting point, from the first patch or its offset.
    pub fn point(&self, at: &PairPoint) -> Point3 {
        let [u, v, _, _] = &at.parameters;
        let span = self.first.span(at.first_chart, at.first_span);
        let at = [u.clone(), v.clone()];
        let point = span.point(&at);
        let Some(distance) = &self.distance else {
            return point;
        };
        let [du, dv] = span.derivatives(&at);
        let normal = du.cross(&dv);
        let length = normal.norm_squared().sqrt();
        point.plus(&normal.scaled(&(&Real::from_rational(distance.clone()) / &length)))
    }
}

/// One patch's spans, in the forms the search reads.
struct Side {
    periodic: [Option<[Rational; 2]>; 2],
    /// The offset distance, when this side is an offset of the patch.
    distance: Option<Rational>,
    spans: Vec<SideSpan>,
}

/// Original global-parameter polynomials, before chart-unit derivative scaling.
struct SpanPolynomials {
    homogeneous: [MPoly; 4],
    tangents: [[MPoly; 3]; 2],
    normal: [MPoly; 3],
}

// Keep legacy preparation inline; only explicit reuse allocates shared data.
#[allow(clippy::large_enum_variant)]
enum PolynomialStorage {
    Owned(SpanPolynomials),
    Shared(Arc<SpanPolynomials>),
}

struct SideSpan {
    chart: usize,
    local_span: usize,
    physical: [bool; 4],
    polynomials: PolynomialStorage,
    /// The offset's polynomials, when this side is an offset.
    offset: Option<OffsetSpan>,
    root: Arc<Cell>,
}

impl std::ops::Deref for SideSpan {
    type Target = SpanPolynomials;

    fn deref(&self) -> &Self::Target {
        match &self.polynomials {
            PolynomialStorage::Owned(value) => value,
            PolynomialStorage::Shared(value) => value,
        }
    }
}

/// One span of an offset `P + d N/|N|`, with `N` the normal above and
/// `L = |N|`. Each tangent of the offset times `W^2 L^3` is `L p1 + p0`.
pub(super) struct OffsetSpan {
    pub(super) distance: Rational,
    pub(super) normal: [MPoly; 3],
    /// `N . N`.
    pub(super) length2: MPoly,
    /// `(N . N)` times the homogeneous tangents.
    pub(super) p1: [[MPoly; 3]; 2],
    /// `d W^2 ((N . N) N_x - (N . N_x) N)` for `x` each parameter.
    pub(super) p0: [[MPoly; 3]; 2],
}

impl OffsetSpan {
    fn new(
        distance: &Rational,
        weight: &MPoly,
        tangents: &[[MPoly; 3]; 2],
        normal: &[MPoly; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let dot = |a: &[MPoly; 3], b: &[MPoly; 3], dec: &mut Decider| -> Result<MPoly, Refusal> {
            let mut total = product(&a[0], &b[0], dec)?;
            for i in 1..3 {
                total = sum(&total, &product(&a[i], &b[i], dec)?, dec)?;
            }
            Ok(total)
        };
        let length2 = dot(normal, normal, dec)?;
        let n = weight.nvars();
        let weight2 = product(weight, weight, dec)?;
        let scale = product(&weight2, &MPoly::constant(n, distance.clone()), dec)?;
        let mut p1 = Vec::with_capacity(2);
        let mut p0 = Vec::with_capacity(2);
        for (axis, tangent) in tangents.iter().enumerate() {
            let moved = normal.each_ref().map(|p| p.partial(axis));
            let along = dot(normal, &moved, dec)?;
            let mut first = Vec::with_capacity(3);
            let mut rest = Vec::with_capacity(3);
            for i in 0..3 {
                first.push(product(&length2, &tangent[i], dec)?);
                let turn = difference(
                    &product(&length2, &moved[i], dec)?,
                    &product(&along, &normal[i], dec)?,
                    dec,
                )?;
                rest.push(product(&scale, &turn, dec)?);
            }
            p1.push(first.try_into().expect("three components"));
            p0.push(rest.try_into().expect("three components"));
        }
        Ok(Self {
            distance: distance.clone(),
            normal: normal.clone(),
            length2,
            p1: p1.try_into().expect("two tangents"),
            p0: p0.try_into().expect("two tangents"),
        })
    }
}

impl Side {
    fn new(
        source: Source<'_>,
        homogeneous: Vec<[MPoly; 4]>,
        rotation: &[[Rational; 3]; 3],
        distance: Option<&Rational>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut spans = Vec::with_capacity(homogeneous.len());
        for (index, ((chart, local_span, span), homogeneous)) in
            source.spans().zip(homogeneous).enumerate()
        {
            dec.check_deadline()?;
            let polynomials = SpanPolynomials::new(&homogeneous, rotation, dec)?;
            let offset = match distance {
                Some(d) => Some(OffsetSpan::new(
                    d,
                    &polynomials.homogeneous[3],
                    &polynomials.tangents,
                    &polynomials.normal,
                    dec,
                )?),
                None => None,
            };
            let root = Arc::new(Cell::span(
                index,
                &span.domain,
                &polynomials.homogeneous,
                offset.as_ref(),
                dec,
            )?);
            spans.push(SideSpan {
                chart,
                local_span,
                physical: source.physical(chart, span),
                polynomials: PolynomialStorage::Owned(polynomials),
                offset,
                root,
            });
        }
        Ok(Self {
            periodic: source.periodic(),
            distance: distance.cloned(),
            spans,
        })
    }
}

/// Each span's homogeneous coordinates as exact rational polynomials.
fn rational_spans(source: Source<'_>) -> Result<Vec<[MPoly; 4]>, Refusal> {
    source
        .spans()
        .map(|(_, _, span)| {
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

/// `a + b`, charged as one exact sum per coefficient.
fn sum(a: &MPoly, b: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    dec.charge_work(size(a).saturating_add(size(b)))?;
    Ok(a + b)
}

/// Everything a cell pair's proof reads, shared by every thread.
struct Context {
    sides: [Side; 2],
    /// The exact rotation applied to both patches; reports undo it.
    rotation: [[Rational; 3]; 3],
    /// `A_k B_w - B_k A_w` in `(u, v, s, t)`, per span pair, built once.
    /// An offset side's homogeneous coordinates are `(L P + d W N, W L)`,
    /// with its length `L` one more variable after `(u, v, s, t)` and one
    /// more equation `L^2 - N . N`.
    equations: Vec<std::sync::OnceLock<Vec<MPoly>>>,
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
        let nvars = 4 + self.lengths();
        let mut next = 4;
        let mut sides = Vec::with_capacity(2);
        let mut extra = Vec::new();
        for (side, (span, targets)) in self
            .sides
            .iter()
            .zip(spans.into_iter().zip([[0, 1], [2, 3]]))
        {
            let span = &side.spans[span];
            let h = span
                .homogeneous
                .each_ref()
                .map(|p| embed(p, targets, nvars));
            let Some(offset) = &span.offset else {
                sides.push(h);
                continue;
            };
            let length = MPoly::var(nvars, next);
            next += 1;
            let scale = product(&h[3], &MPoly::constant(nvars, offset.distance.clone()), dec)?;
            let mut lifted = Vec::with_capacity(4);
            for (k, normal) in offset.normal.iter().enumerate() {
                let normal = embed(normal, targets, nvars);
                lifted.push(sum(
                    &product(&length, &h[k], dec)?,
                    &product(&scale, &normal, dec)?,
                    dec,
                )?);
            }
            lifted.push(product(&length, &h[3], dec)?);
            let square = product(&length, &length, dec)?;
            extra.push(difference(
                &square,
                &embed(&offset.length2, targets, nvars),
                dec,
            )?);
            sides.push(lifted.try_into().expect("four coordinates"));
        }
        let [a, b] = [&sides[0], &sides[1]];
        let mut equations = Vec::with_capacity(3 + extra.len());
        for k in 0..3 {
            let left = product(&a[k], &b[3], dec)?;
            let right = product(&b[k], &a[3], dec)?;
            equations.push(difference(&left, &right, dec)?);
        }
        equations.extend(extra);
        let _ = slot.set(equations);
        Ok(())
    }

    /// How many sides are offsets, each adding its length as an unknown.
    fn lengths(&self) -> usize {
        self.sides
            .iter()
            .filter(|side| side.distance.is_some())
            .count()
    }

    /// A span pair's equations; the writer built them before the level.
    fn equations(&self, spans: [usize; 2]) -> &[MPoly] {
        self.equations[self.index(spans)]
            .get()
            .expect("the writer builds a span pair's equations before proving it")
    }

    /// Whether `value` on `axis` is a border of the physical domain.
    fn physical(&self, axis: usize, value: &Rational, span: usize) -> bool {
        let span = &self.sides[axis / 2].spans[span];
        (0..2).any(|end| {
            span.physical[2 * (axis % 2) + end] && &span.root.rect[axis % 2][end] == value
        })
    }
}

/// The exact rotation of a small integer quaternion that best brings the
/// first patch's mean normal over its span corners to the `z` axis and its
/// mean `u` tangent to the `x` axis. A binary64 search, charged a fixed
/// amount, only chooses the quaternion; the rotation is exact, and any
/// rotation keeps the meeting. The identity stays unless another quaternion
/// is clearly better.
fn upright(
    source: Source<'_>,
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
    for ((_, _, span), homogeneous) in source.spans().zip(spans) {
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
    /// The complete exact meeting with another rational patch, over both
    /// closed parameter domains, borders included.
    ///
    /// Every branch is returned as an ordered cover of cell pairs, each
    /// accepted by an exact certificate and bounded by exactly isolated
    /// border points. No meeting point of the closed patches is dropped: a
    /// pair leaves the search only when it provably holds none, or when all
    /// it holds lie on an internal cell border that the neighbouring cell
    /// owns. An empty answer is therefore a proof that the closed patches do
    /// not meet.
    ///
    /// Where this cannot be certified, the search refuses; it never returns
    /// a near-enough answer. A face whose equations are exactly dependent,
    /// where a branch may run along a cell or patch border, refuses with
    /// `intersect.spline-pair-boundary-overlap`. An exact tangent contact
    /// refuses with `intersect.spline-pair-tangent`. Any other cell pair
    /// that no certificate settles within the split and work allowances,
    /// such as a contact on the border, refuses with
    /// `intersect.spline-pair-unresolved`.
    pub fn patch_meeting(
        &self,
        other: &RationalPatch,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        meeting(
            [Source::Patch(self), Source::Patch(other)],
            [None, None],
            dec,
        )
    }

    /// This patch's exact offset `P + d N / |N|` at signed distance `d`.
    pub fn offset(&self, distance: Rational) -> OffsetPatch {
        OffsetPatch {
            support: self.clone(),
            distance,
        }
    }
}

/// The meeting of two patches, each itself or its offset at a nonzero
/// distance.
fn meeting(
    patches: [Source<'_>; 2],
    distances: [Option<Rational>; 2],
    dec: &mut Decider,
) -> Result<PatchPairIntersection, Refusal> {
    let [this, other] = patches;
    // With an operation allowance active, its counter includes the capped
    // root certificates, which run on deciders of their own.
    let budget = catallax_numbers::operation_budget_active();
    let used = || catallax_numbers::current_operation_work().used;
    let start = if budget { used() } else { dec.stats().work };
    // A rotation changes no parameter of the meeting, and boxes along
    // the axes follow tilted surfaces once they are turned upright. An
    // offset of the turned patch is the turned offset.
    let spans = [rational_spans(this)?, rational_spans(other)?];
    let rotation = upright(this, &spans[0], dec)?;
    let [first, second] = spans;
    let sides = [
        Side::new(this, first, &rotation, distances[0].as_ref(), dec)?,
        Side::new(other, second, &rotation, distances[1].as_ref(), dec)?,
    ];
    finish(patches, distances, sides, rotation, (budget, start), dec)
}

fn finish(
    patches: [Source<'_>; 2],
    distances: [Option<Rational>; 2],
    sides: [Side; 2],
    rotation: [[Rational; 3]; 3],
    (budget, start): (bool, u64),
    dec: &mut Decider,
) -> Result<PatchPairIntersection, Refusal> {
    let pairs = sides[0].spans.len() * sides[1].spans.len();
    let context = Context {
        sides,
        rotation,
        equations: (0..pairs).map(|_| std::sync::OnceLock::new()).collect(),
    };
    let reuse = patches
        .iter()
        .any(|source| matches!(source, Source::Original(_)));
    let (arcs, mut search) = run::search(&context, reuse, dec)?;
    let branches = join::branches(&context, arcs, dec)?;
    search.work = if budget {
        catallax_numbers::current_operation_work().used
    } else {
        dec.stats().work
    } - start;
    Ok(PatchPairIntersection {
        first: patches[0].retained(),
        distance: distances[0].clone(),
        second: patches[1].retained(),
        second_distance: distances[1].clone(),
        rotation: context.rotation,
        branches,
        search,
    })
}
