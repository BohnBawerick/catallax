//! Finite rational charts with a checked global embedding projection.

use std::sync::{Arc, OnceLock};

use catallax_numbers::{
    Ball, Decider, MPoly, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof,
};

use crate::{BSplineCurve, BSplineSurface, Isometry, Point3, Pose, SplineBoundary, Vector3, math};

mod grid;
pub use grid::{RationalGridEmbedding, RationalPatchGrid};
mod pair;
pub use pair::{
    OffsetPatch, PairArc, PairCertificate, PairPoint, PatchPairBranch, PatchPairIntersection,
    PatchPairSearch,
};
mod query;
mod range;
mod section;
pub(crate) use section::{Cover, search_cover, vanishing_side_row};
pub use section::{
    PatchSectionBranch, RationalPatchSection, SectionBranchApproximation, SectionCandidate,
    SectionHits, SectionSearch, SectionSurface,
};
mod boundary;
mod join;

pub(crate) const MAX_PROOF_CELLS: usize = 16_384;
const MAX_PROOF_DEPTH: u32 = 24;

type Domain = [[Rational; 2]; 2];
type BezierSpanView<'a> = (&'a Domain, &'a [[Real; 4]]);

#[derive(Debug)]
struct Span {
    indices: [usize; 2],
    domain: Domain,
    coordinates: [RealPoly; 4],
    first: [[RealPoly; 4]; 2],
    normal: OnceLock<[RealPoly; 3]>,
}

impl Span {
    fn new(indices: [usize; 2], domain: Domain, coordinates: [RealPoly; 4]) -> Self {
        let first = std::array::from_fn(|axis| coordinates.each_ref().map(|p| p.partial(axis)));
        Self {
            indices,
            domain,
            coordinates,
            first,
            normal: OnceLock::new(),
        }
    }

    fn point(&self, at: &[Real; 2]) -> Point3 {
        let weight = self.coordinates[3].eval_real(at);
        Point3::new(
            &self.coordinates[0].eval_real(at) / &weight,
            &self.coordinates[1].eval_real(at) / &weight,
            &self.coordinates[2].eval_real(at) / &weight,
        )
    }

    fn derivatives(&self, at: &[Real; 2]) -> [Vector3; 2] {
        let weight = self.coordinates[3].eval_real(at);
        let point = self.point(at);
        let xyz = [&point.x, &point.y, &point.z];
        std::array::from_fn(|axis| {
            let dw = self.first[axis][3].eval_real(at);
            let value: [Real; 3] = std::array::from_fn(|i| {
                &(&self.first[axis][i].eval_real(at) - &(xyz[i] * &dw)) / &weight
            });
            let [x, y, z] = value;
            Vector3::new(x, y, z)
        })
    }

    fn normal_polynomials(&self) -> &[RealPoly; 3] {
        self.normal.get_or_init(|| {
            let tangents: [[RealPoly; 3]; 2] = std::array::from_fn(|axis| {
                std::array::from_fn(|i| {
                    &self.first[axis][i] * &self.coordinates[3]
                        - &self.coordinates[i] * &self.first[axis][3]
                })
            });
            let [u, v] = &tangents;
            [
                &u[1] * &v[2] - &u[2] * &v[1],
                &u[2] * &v[0] - &u[0] * &v[2],
                &u[0] * &v[1] - &u[1] * &v[0],
            ]
        })
    }
}

#[derive(Clone, Debug)]
struct ProjectionCell {
    domain: Domain,
    lower: Rational,
}

#[derive(Debug)]
struct BezierSpan {
    domain: Domain,
    controls: Vec<[Real; 4]>,
}

#[derive(Debug)]
struct Data {
    carrier: BSplineSurface,
    domain: Domain,
    widths: [Rational; 2],
    spans: Vec<Span>,
    projection: [Vector3; 2],
    cover: Arc<[ProjectionCell]>,
    lower: Rational,
    beziers: OnceLock<Vec<BezierSpan>>,
    // The same spans for run-ahead workers. Filling `beziers` is charged, so a
    // worker must not fill it for the writer; it fills this cell instead.
    ahead_beziers: OnceLock<Vec<BezierSpan>>,
}

/// An exact finite rational patch whose projected chart is strongly monotone.
/// The positive symmetric-Jacobian bound proves both regularity and injectivity
/// over the complete knot rectangle, including piecewise-smooth knot joins.
/// A failed proof is an implementation gap, never an impossible-geometry claim.
#[derive(Clone, Debug)]
pub struct RationalPatch {
    data: Arc<Data>,
}

fn unproven(detail: &str) -> Refusal {
    Refusal::new("shapes.rational-patch-unproven", detail)
}

fn domain_error(detail: &str) -> Refusal {
    Refusal::new("shapes.rational-patch-domain", detail)
}

// Work profiles attribute each proof to the line that asked for it.
#[track_caller]
pub(crate) fn exact_sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(unproven(
            "a patch proof requires exact equality, not bounded coincidence",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

fn midpoint(domain: &Domain) -> [Rational; 2] {
    std::array::from_fn(|axis| (&domain[axis][0] + &domain[axis][1]).mul_pow2(-1))
}

fn box_balls(domain: &Domain, bits: u32) -> [Ball; 2] {
    domain
        .each_ref()
        .map(|[lo, hi]| Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits)))
}

fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&a[0] * &b[0]) + &(&a[1] * &b[1]) + &a[2] * &b[2]
}

fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn split(domain: &Domain, widths: &[Rational; 2]) -> [Domain; 2] {
    let extent = domain.each_ref().map(|[lo, hi]| hi - lo);
    let axis = usize::from(&extent[1] * &widths[0] > &extent[0] * &widths[1]);
    let middle = (&domain[axis][0] + &domain[axis][1]).mul_pow2(-1);
    let mut low = domain.clone();
    let mut high = domain.clone();
    low[axis][1] = middle.clone();
    high[axis][0] = middle;
    [low, high]
}

fn insert(values: &mut Vec<Real>, value: Real, dec: &mut Decider) -> Result<(), Refusal> {
    for i in 0..values.len() {
        match exact_sign(&(&value - &values[i]), dec)? {
            Sign::Negative => {
                values.insert(i, value);
                return Ok(());
            }
            Sign::Zero => return Ok(()),
            Sign::Positive => (),
        }
    }
    values.push(value);
    Ok(())
}

fn projection_margin(
    derivatives: &[[Ball; 3]; 2],
    frame: &[[Ball; 3]; 2],
    scales: &[Ball; 2],
    bits: u32,
) -> Rational {
    let diagonal = [
        &dot(&frame[0], &derivatives[0]) * &scales[0],
        &dot(&frame[1], &derivatives[1]) * &scales[1],
    ];
    let off = (&(&dot(&frame[0], &derivatives[1]) * &scales[1])
        + &(&dot(&frame[1], &derivatives[0]) * &scales[0]))
        * Ball::from_rational(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
    diagonal[0].lo().min(diagonal[1].lo()) - magnitude(&off)
}

impl RationalPatch {
    /// Proves the complete chart using the dual of its centre tangents as a
    /// projection. No control is changed and no sample supplies the proof.
    /// Patches without a proved global projection retain a typed refusal.
    pub fn new(carrier: BSplineSurface, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::certified(carrier, MAX_PROOF_DEPTH, MAX_PROOF_CELLS, dec)
    }

    pub(crate) fn new_cell(carrier: BSplineSurface, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::certified(carrier, 0, 1, dec)
    }

    fn certified(
        carrier: BSplineSurface,
        max_depth: u32,
        max_cells: usize,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let domain = carrier.domain().map(|bounds| bounds.map(Clone::clone));
        let widths = domain.each_ref().map(|[lo, hi]| hi - lo);
        let knots = carrier.knots();
        let mut spans = Vec::new();
        for indices in carrier.span_indices() {
            let bounds = std::array::from_fn(|axis| {
                [
                    knots[axis][indices[axis]].clone(),
                    knots[axis][indices[axis] + 1].clone(),
                ]
            });
            spans.push(Span::new(
                indices,
                bounds,
                carrier.span_polynomials(indices, dec)?,
            ));
        }
        let center_rational = midpoint(&domain);
        let span = spans
            .iter()
            .find(|span| {
                center_rational
                    .iter()
                    .zip(&span.domain)
                    .all(|(value, [lo, hi])| value >= lo && value <= hi)
            })
            .expect("active spline rectangle has a centre span");
        let center = center_rational.map(Real::from_rational);
        let [a, b] = span.derivatives(&center);
        // Use the unit-domain metric, not the arbitrary units of authored knots.
        let a = a.scaled(&Real::from_rational(widths[0].clone()));
        let b = b.scaled(&Real::from_rational(widths[1].clone()));
        let aa = a.norm_squared();
        let bb = b.norm_squared();
        let ab = a.dot(&b);
        let determinant = &aa * &bb - ab.square();
        if exact_sign(&determinant, dec)? != Sign::Positive {
            return Err(unproven(
                "centre tangents do not provide an independent projection",
            ));
        }
        let projection = [
            a.scaled(&bb)
                .minus(&b.scaled(&ab))
                .scaled(&(&Real::one() / &determinant)),
            b.scaled(&aa)
                .minus(&a.scaled(&ab))
                .scaled(&(&Real::one() / &determinant)),
        ];
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(unproven(
                "no arithmetic precision is available for the embedding proof",
            ));
        }
        let frame: [[Ball; 3]; 2] = [
            math::budgeted_vector(&projection[0], bits, dec)?,
            math::budgeted_vector(&projection[1], bits, dec)?,
        ];
        let scales = widths
            .each_ref()
            .map(|width| Ball::from_rational(width, bits));
        let mut cover = Vec::new();
        let mut pending: Vec<_> = spans
            .iter()
            .rev()
            .map(|span| (span.indices, span.domain.clone(), 0u32))
            .collect();
        let mut workspace = carrier.span_workspace();
        let mut visited = 0usize;
        let mut lower: Option<Rational> = None;
        while let Some((indices, cell, depth)) = pending.pop() {
            dec.check_deadline()?;
            visited += 1;
            if visited > max_cells {
                return Err(unproven(
                    "the embedding proof exceeded its finite cell allowance",
                ));
            }
            let derivatives =
                carrier.derivatives_on_span(indices, &cell, bits, dec, &mut workspace)?;
            let margin = projection_margin(&derivatives, &frame, &scales, bits);
            if margin > Rational::zero() {
                lower = Some(lower.map_or_else(|| margin.clone(), |old| old.min(margin.clone())));
                cover.push(ProjectionCell {
                    domain: cell,
                    lower: margin,
                });
            } else {
                if depth == max_depth {
                    return Err(unproven(
                        "no positive whole-domain projected Jacobian bound was proved",
                    ));
                }
                let [low, high] = split(&cell, &widths);
                pending.push((indices, high, depth + 1));
                pending.push((indices, low, depth + 1));
            }
        }
        join::check(&carrier, &spans, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                carrier,
                domain,
                widths,
                spans,
                projection,
                cover: cover.into(),
                lower: lower.expect("nonempty certified parameter cover"),
                beziers: OnceLock::new(),
                ahead_beziers: OnceLock::new(),
            }),
        })
    }

    pub(crate) fn restricted(&self, bounds: Domain, dec: &Decider) -> Result<Self, Refusal> {
        assert_eq!(self.data.spans.len(), 1, "atlas cells have one Bezier span");
        assert_eq!(
            self.data.cover.len(),
            1,
            "atlas cells have one proof domain"
        );
        let parent = &self.data.spans[0];
        let carrier = self.carrier().bezier_patch(
            parent.indices,
            bounds.each_ref().map(|part| part.each_ref()),
            dec,
        )?;
        let widths = bounds.each_ref().map(|[lo, hi]| hi - lo);
        let ratios: [Rational; 2] =
            std::array::from_fn(|axis| &widths[axis] / &self.data.widths[axis]);
        let minimum = ratios[0].clone().min(ratios[1].clone());
        let minimum_squared = &minimum * &minimum;
        // Diagonal congruence preserves the parent's positive symmetric
        // Jacobian. Rescale so its certified unit-domain lower bound remains
        // valid without reevaluating an already proved cell.
        let projection = std::array::from_fn(|axis| {
            self.data.projection[axis]
                .scaled(&Real::from_rational(&ratios[axis] / &minimum_squared))
        });
        let span = Span {
            indices: parent.indices,
            domain: bounds.clone(),
            coordinates: parent.coordinates.clone(),
            first: parent.first.clone(),
            normal: parent.normal.clone(),
        };
        Ok(Self {
            data: Arc::new(Data {
                carrier,
                domain: bounds.clone(),
                widths,
                spans: vec![span],
                projection,
                cover: Arc::from([ProjectionCell {
                    domain: bounds,
                    lower: self.data.lower.clone(),
                }]),
                lower: self.data.lower.clone(),
                beziers: OnceLock::new(),
                ahead_beziers: OnceLock::new(),
            }),
        })
    }

    pub(crate) fn seam(
        &self,
        other: &Self,
        axis: usize,
        check_position: bool,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        assert_eq!(self.data.spans.len(), 1);
        assert_eq!(other.data.spans.len(), 1);
        let a = &self.data.spans[0];
        let b = &other.data.spans[0];
        let left_value = RealPoly::constant(2, Real::from_rational(a.domain[axis][1].clone()));
        let right_value = RealPoly::constant(2, Real::from_rational(b.domain[axis][0].clone()));
        let zero = |poly: &RealPoly, dec: &mut Decider| -> Result<bool, Refusal> {
            for (_, coefficient) in poly.terms() {
                if exact_sign(&coefficient, dec)? != Sign::Zero {
                    return Ok(false);
                }
            }
            Ok(true)
        };
        if check_position {
            let left = a
                .coordinates
                .each_ref()
                .map(|poly| poly.substitute(axis, &left_value));
            let right = b
                .coordinates
                .each_ref()
                .map(|poly| poly.substitute(axis, &right_value));
            for (a, b) in left[..3].iter().zip(&right[..3]) {
                if !zero(&(a * &right[3] - b * &left[3]), dec)? {
                    return Err(Refusal::new(
                        "shapes.rational-atlas-open",
                        "identified atlas rows have different exact positions",
                    ));
                }
            }
        }
        let left = a
            .normal_polynomials()
            .each_ref()
            .map(|poly| poly.substitute(axis, &left_value));
        let right = b
            .normal_polynomials()
            .each_ref()
            .map(|poly| poly.substitute(axis, &right_value));
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            if !zero(&(&left[i] * &right[j] - &left[j] * &right[i]), dec)? {
                return Ok(false);
            }
        }
        let mut at = [Real::zero(), Real::zero()];
        at[1 - axis] =
            Real::from_rational((&a.domain[1 - axis][0] + &a.domain[1 - axis][1]).mul_pow2(-1));
        let left = left.map(|poly| poly.eval_real(&at));
        let right = right.map(|poly| poly.eval_real(&at));
        let dot = left
            .iter()
            .zip(&right)
            .fold(Real::zero(), |sum, (a, b)| sum + a * b);
        if exact_sign(&dot, dec)? != Sign::Positive {
            return Err(Refusal::new(
                "shapes.rational-atlas-folded",
                "identified regular atlas rows have opposing normal orientations",
            ));
        }
        Ok(true)
    }

    /// The unchanged exact rational definition, including knots and weights.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }

    /// Each knot span's rectangle and exact homogeneous polynomials in the
    /// global parameters, or `None` when a coefficient is not rational.
    pub(crate) fn rational_span_polynomials(&self) -> Option<Vec<(Domain, [MPoly; 4])>> {
        self.data
            .spans
            .iter()
            .map(|span| {
                let coordinates: Option<Vec<MPoly>> =
                    span.coordinates.iter().map(RealPoly::as_rational).collect();
                Some((
                    span.domain.clone(),
                    coordinates?.try_into().expect("four coordinates"),
                ))
            })
            .collect()
    }

    /// The same rational map over a larger parameter rectangle, for a patch
    /// of one knot span. Its exact polynomials are written as Bernstein
    /// controls on `domain`, and the new patch is proved like any other. A
    /// spline fillet's supports are extended past their faces, so that the
    /// fillet can reach and be cut by the faces' neighbours.
    pub fn extended(&self, domain: [[Rational; 2]; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        let spans = self.rational_span_polynomials().ok_or_else(|| {
            domain_error("only a patch with rational coefficients can be extended")
        })?;
        let [(old, h)] = <[_; 1]>::try_from(spans)
            .map_err(|_| domain_error("only a patch of one knot span can be extended"))?;
        if (0..2).any(|i| domain[i][0] > old[i][0] || domain[i][1] < old[i][1]) {
            return Err(domain_error(
                "an extended domain must contain the patch's own domain",
            ));
        }
        let degree = self.carrier().degrees();
        // Each size-dependent exact operation is reserved before it runs.
        let mul = |a: &Rational, b: &Rational, dec: &mut Decider| -> Result<Rational, Refusal> {
            dec.reserve_exact_arithmetic(a, b, true)?;
            Ok(a * b)
        };
        let add = |a: &Rational, b: &Rational, dec: &mut Decider| -> Result<Rational, Refusal> {
            dec.reserve_exact_arithmetic(a, b, false)?;
            Ok(a + b)
        };
        // Each coordinate's coefficients of sigma^i omega^j, where
        // u = a + (b - a) sigma and v = c + (d - c) omega.
        let mut width = Vec::with_capacity(2);
        for side in &domain {
            width.push(add(&side[1], &-&side[0], dec)?);
        }
        let binomial = |n: usize, k: usize| -> Rational {
            (0..k).fold(Rational::one(), |out, i| {
                out * Rational::from_ratio((n - i) as i64, (i + 1) as i64).expect("nonzero")
            })
        };
        // Powers of each axis's start and width, up to its degree.
        let mut starts: Vec<Vec<Rational>> = Vec::with_capacity(2);
        let mut widths: Vec<Vec<Rational>> = Vec::with_capacity(2);
        for axis in 0..2 {
            let mut start = vec![Rational::one()];
            let mut span = vec![Rational::one()];
            for _ in 0..degree[axis] {
                dec.check_deadline()?;
                start.push(mul(start.last().expect("one"), &domain[axis][0], dec)?);
                span.push(mul(span.last().expect("one"), &width[axis], dec)?);
            }
            starts.push(start);
            widths.push(span);
        }
        // The coefficient of sigma^p in (a + w sigma)^i, per axis.
        let mut shifts: Vec<Vec<Vec<Rational>>> = Vec::with_capacity(2);
        for axis in 0..2 {
            let mut table = Vec::with_capacity(degree[axis] + 1);
            for i in 0..=degree[axis] {
                let mut row = Vec::with_capacity(i + 1);
                for p in 0..=i {
                    dec.check_deadline()?;
                    let power = mul(&starts[axis][i - p], &widths[axis][p], dec)?;
                    row.push(mul(&binomial(i, p), &power, dec)?);
                }
                table.push(row);
            }
            shifts.push(table);
        }
        let mut controls = Vec::with_capacity((degree[0] + 1) * (degree[1] + 1));
        let mut weights = Vec::with_capacity(controls.capacity());
        let mut grids: Vec<Vec<Vec<Rational>>> = Vec::with_capacity(4);
        for poly in &h {
            let mut grid = vec![vec![Rational::zero(); degree[1] + 1]; degree[0] + 1];
            for (monomial, coefficient) in poly.terms() {
                let e = monomial.exponents();
                let (i, j) = (e[0] as usize, e[1] as usize);
                if i > degree[0] || j > degree[1] {
                    return Err(domain_error("a span polynomial exceeds its degree"));
                }
                for (p, row) in grid.iter_mut().enumerate().take(i + 1) {
                    let along = mul(&coefficient, &shifts[0][i][p], dec)?;
                    for (q, slot) in row.iter_mut().enumerate().take(j + 1) {
                        dec.check_deadline()?;
                        let term = mul(&along, &shifts[1][j][q], dec)?;
                        *slot = add(slot, &term, dec)?;
                    }
                }
            }
            grids.push(grid);
        }
        for a in 0..=degree[0] {
            for b in 0..=degree[1] {
                let mut point: [Rational; 4] = std::array::from_fn(|_| Rational::zero());
                for (value, grid) in point.iter_mut().zip(&grids) {
                    for (i, row) in grid.iter().enumerate().take(a + 1) {
                        for (j, coefficient) in row.iter().enumerate().take(b + 1) {
                            dec.check_deadline()?;
                            let factor = binomial(a, i) / binomial(degree[0], i)
                                * (binomial(b, j) / binomial(degree[1], j));
                            let term = mul(&factor, coefficient, dec)?;
                            *value = add(value, &term, dec)?;
                        }
                    }
                }
                let [x, y, z, w] = point;
                if w.sign() != Sign::Positive {
                    return Err(domain_error(
                        "an extended patch needs positive weights over its whole domain",
                    ));
                }
                let mut cartesian = Vec::with_capacity(3);
                for value in [&x, &y, &z] {
                    // A division costs as a multiplication of its operands.
                    dec.reserve_exact_arithmetic(value, &w, true)?;
                    cartesian.push(Real::from_rational(value / &w));
                }
                let [x, y, z]: [Real; 3] = cartesian.try_into().expect("three coordinates");
                controls.push(Point3::new(x, y, z));
                weights.push(Real::from_rational(w));
            }
        }
        let knots = |axis: usize| -> Vec<Rational> {
            let mut out = vec![domain[axis][0].clone(); degree[axis] + 1];
            out.extend(std::iter::repeat_n(
                domain[axis][1].clone(),
                degree[axis] + 1,
            ));
            out
        };
        let carrier = BSplineSurface::new(
            degree,
            [knots(0), knots(1)],
            [degree[0] + 1, degree[1] + 1],
            controls,
            weights,
            dec,
        )?;
        Self::new(carrier, dec)
    }

    /// Closed complete parameter rectangle in original knot coordinates.
    pub fn domain(&self) -> &[[Rational; 2]; 2] {
        &self.data.domain
    }

    /// Positive lower bound on the symmetric projected Jacobian over the whole
    /// domain. It is a regularity/embedding proof, not a positional tolerance.
    pub fn embedding_lower_bound(&self) -> &Rational {
        &self.data.lower
    }

    /// The checked finite parameter cover and each cell's positive projected
    /// Jacobian margin. Together the cells cover the complete active domain.
    pub fn embedding_cover(&self) -> impl Iterator<Item = (&[[Rational; 2]; 2], &Rational)> {
        self.data
            .cover
            .iter()
            .map(|cell| (&cell.domain, &cell.lower))
    }

    /// Exact positive homogeneous Bézier controls on each complete original
    /// knot rectangle. The controls use unit-interval local u and v parameters.
    pub fn bezier_spans<'a>(
        &'a self,
        dec: &Decider,
    ) -> Result<impl Iterator<Item = BezierSpanView<'a>> + 'a + use<'a>, Refusal> {
        let beziers = if catallax_numbers::shared_caches() {
            &self.data.beziers
        } else {
            &self.data.ahead_beziers
        };
        // A filled cache skips the controls' charge, so a reader waits for
        // the fill.
        #[cfg(feature = "work-span")]
        let span_owner = std::ptr::from_ref(beziers) as usize as u64;
        #[cfg(feature = "work-span")]
        if beziers.get().is_some() {
            catallax_numbers::work_span::read_at("shared cache: patch spans", span_owner, 0);
        }
        if beziers.get().is_none() {
            let knots = self.carrier().knots();
            let mut pieces = Vec::with_capacity(self.data.spans.len());
            for span in self.carrier().span_indices() {
                dec.check_deadline()?;
                let domain = std::array::from_fn(|axis| {
                    [
                        knots[axis][span[axis]].clone(),
                        knots[axis][span[axis] + 1].clone(),
                    ]
                });
                pieces.push(BezierSpan {
                    domain,
                    controls: self.carrier().span_bezier_controls(span, dec)?,
                });
            }
            // Only successful, policy-independent exact definitions are cached.
            // A racing reader may finish first; neither reader caches a refusal.
            let _ = beziers.set(pieces);
            #[cfg(feature = "work-span")]
            catallax_numbers::work_span::write_at("shared cache: patch spans", span_owner, 0);
        }
        Ok(beziers
            .get()
            .expect("completed exact span cache")
            .iter()
            .map(|span| (&span.domain, span.controls.as_slice())))
    }

    /// A boundary view sharing the authoritative carrier's exact control table.
    pub fn boundary(&self, side: SplineBoundary) -> BSplineCurve {
        BSplineCurve::boundary(self.carrier(), side)
    }

    fn span_at(&self, at: &[Real; 2], dec: &mut Decider) -> Result<&Span, Refusal> {
        for span in &self.data.spans {
            if contains(&span.domain, at, dec)? {
                return Ok(span);
            }
        }
        Err(domain_error("the exact parameter leaves the finite patch"))
    }

    /// Exact algebraic point, without rounded intermediate coordinates.
    pub fn point_at(&self, at: &[Real; 2], dec: &mut Decider) -> Result<Point3, Refusal> {
        Ok(self.span_at(at, dec)?.point(at))
    }

    /// Exact first partials on the selected finite knot span.
    pub fn derivatives_at(
        &self,
        at: &[Real; 2],
        dec: &mut Decider,
    ) -> Result<[Vector3; 2], Refusal> {
        Ok(self.span_at(at, dec)?.derivatives(at))
    }

    /// Outward enclosures of the authoritative Euclidean controls and weights.
    pub fn enclosed_controls(&self, dec: &mut Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let bits = dec.policy().max_precision_bits.min(192);
        if bits == 0 {
            return Err(unproven("no control enclosure precision is available"));
        }
        self.carrier()
            .controls()
            .iter()
            .zip(self.carrier().weights())
            .map(|(p, w)| {
                Ok((
                    math::budgeted_point(p, bits, dec)?,
                    dec.enclose_budgeted_at(w, bits)?,
                ))
            })
            .collect()
    }

    /// Exact rigid transport. Dot products retain the same certified projection
    /// inequalities; reflecting the tensor chart reverses its natural normal.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }

    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        let carrier = self.carrier().mapped(map);
        let matrix = map.matrix();
        let translation = [&map.translation.x, &map.translation.y, &map.translation.z];
        let spans = self
            .data
            .spans
            .iter()
            .map(|span| {
                let coordinates = std::array::from_fn(|axis| {
                    if axis == 3 {
                        return span.coordinates[3].clone();
                    }
                    let mut value = span.coordinates[3].scale(translation[axis]);
                    for (coefficient, coordinate) in matrix[axis].iter().zip(&span.coordinates[..3])
                    {
                        value = value + coordinate.scale(coefficient);
                    }
                    value
                });
                Span::new(span.indices, span.domain.clone(), coordinates)
            })
            .collect();
        Self {
            data: Arc::new(Data {
                carrier,
                domain: self.data.domain.clone(),
                widths: self.data.widths.clone(),
                spans,
                projection: self.data.projection.each_ref().map(|p| map.apply_vector(p)),
                cover: self.data.cover.clone(),
                lower: self.data.lower.clone(),
                beziers: OnceLock::new(),
                ahead_beziers: OnceLock::new(),
            }),
        }
    }

    /// Positive identity of the complete retained definition and checked cover.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
}

fn contains(domain: &Domain, at: &[Real; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    for (value, [lo, hi]) in at.iter().zip(domain) {
        if exact_sign(&(value - &Real::from_rational(lo.clone())), dec)? == Sign::Negative
            || exact_sign(&(&Real::from_rational(hi.clone()) - value), dec)? == Sign::Negative
        {
            return Ok(false);
        }
    }
    Ok(true)
}
