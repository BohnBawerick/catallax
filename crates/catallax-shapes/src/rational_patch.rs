//! Finite rational charts with a checked global embedding projection.

use std::sync::{Arc, OnceLock};

use catallax_numbers::{
    Ball, Decider, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof,
};

use crate::{BSplineCurve, BSplineSurface, Isometry, Point3, Pose, SplineBoundary, Vector3, math};

mod grid;
pub use grid::{RationalGridEmbedding, RationalPatchGrid};
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
