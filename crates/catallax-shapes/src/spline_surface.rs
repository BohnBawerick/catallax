//! Finite regular spline supports with an independently checked exact target.

use crate::{
    BSplineCurve, BSplineSurface, EllipticOffsetDefinition, OffsetSweepProof, Point3, Pose,
    SplineBoundary, SplineCurveParameter, Vector3, math,
};
use catallax_numbers::{Ball, Decider, EvidenceBound, Rational, Real, Refusal, Sign};
use std::{
    collections::{BTreeMap, btree_map::Entry},
    sync::Arc,
};

#[derive(Debug)]
struct Data {
    carrier: BSplineSurface,
    definition: EllipticOffsetDefinition,
    axial: [Real; 2],
    height: Real,
    axis: Vector3,
    edge: BSplineCurve,
    bound: EvidenceBound,
    transport: Option<Transport>,
}

#[derive(Debug)]
struct Transport {
    source: SplineSurface,
    pose: Pose,
    inverse: Pose,
}

/// A checked finite, outward-oriented spline sweep. Unlike raw control data,
/// this support carries regularity, convex embedding and full target coverage.
/// It does not prove the validity or recipe identity of a body using the patch.
#[derive(Clone, Debug)]
pub struct SplineSurface {
    data: Arc<Data>,
}

/// A spline carrier with the retained definition and local proofs needed by
/// geometric queries. Each producer supplies its own checked chart family;
/// raw control data alone do not authorize surface-side decisions.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum SplineSupport {
    /// Regular convex normal offset of an elliptic cylinder, with finite caps.
    EllipticOffset(SplineSurface),
}

impl From<SplineSurface> for SplineSupport {
    fn from(surface: SplineSurface) -> Self {
        Self::EllipticOffset(surface)
    }
}
impl SplineSupport {
    /// The genuine finite B-spline carrier, independent of its exact target.
    pub fn carrier(&self) -> &BSplineSurface {
        match self {
            Self::EllipticOffset(surface) => surface.carrier(),
        }
    }
    /// Exact rigid transport of the carrier, retained law and local proofs.
    pub fn transformed(&self, pose: &Pose) -> Self {
        match self {
            Self::EllipticOffset(surface) => Self::EllipticOffset(surface.transformed(pose)),
        }
    }
    /// Positive expression-node identity for deterministic cache reuse.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::EllipticOffset(a), Self::EllipticOffset(b)) => a.shares_definition_nodes(b),
        }
    }
    /// Complete finite parameter cover for a spatial box. A returned cover
    /// does not prove incidence; separated branches and the periodic seam are
    /// reported as ambiguous rather than choosing one.
    pub fn parameters_of(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        match self {
            Self::EllipticOffset(surface) => surface.parameters_of_box(point, dec),
        }
    }

    /// Oriented side only where the finite axial chart supplies a cover.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        match self {
            Self::EllipticOffset(surface) => surface.side_at(point, dec),
        }
    }
    /// Outward normal at an actual finite carrier point; no global implicit
    /// gradient is fabricated for a parametric surface.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        match self {
            Self::EllipticOffset(surface) => surface
                .parameters_at(point, dec)?
                .ok_or_else(|| {
                    Refusal::new(
                        "shapes.bspline-domain",
                        "normal point is not on the finite spline support",
                    )
                })?
                .normal(dec),
        }
    }
}

/// Exact finite chart preimage, bound to its checked support and curve span.
#[derive(Debug)]
pub struct SplinePointParameters<'a> {
    surface: &'a SplineSurface,
    u: SplineCurveParameter<'a>,
    v: Real,
}
impl SplinePointParameters<'_> {
    /// Cross-section parameter in [0,4], with the seam represented first at zero.
    pub fn u(&self) -> &Real {
        self.u.value()
    }
    /// Axial parameter in [0,1].
    pub fn v(&self) -> &Real {
        &self.v
    }
    /// Exact represented point, not the ideal target point.
    pub fn point(&self, dec: &Decider) -> Result<Point3, Refusal> {
        let (source, pose) = self.surface.source_frame();
        let point = self
            .u
            .point(dec)?
            .plus(&source.axis().scaled(&(&self.v * source.height())));
        Ok(match pose {
            None => point,
            Some(pose) => pose.apply_point(&point),
        })
    }
    /// Exact chart tangents; neither is normalised.
    pub fn first(&self, dec: &Decider) -> Result<[Vector3; 2], Refusal> {
        let mut tangent = self.u.first(dec)?;
        if let (_, Some(pose)) = self.surface.source_frame() {
            tangent = pose.apply_vector(&tangent);
        }
        Ok([tangent, self.surface.axis().scaled(self.surface.height())])
    }
    /// Outward normal covector at this actual spline point, not an implicit gradient.
    pub fn normal(&self, dec: &Decider) -> Result<Vector3, Refusal> {
        let (source, pose) = self.surface.source_frame();
        let normal = self.u.first(dec)?.cross(source.axis());
        Ok(match pose {
            None => normal,
            Some(pose) => pose.apply_vector(&normal),
        })
    }
}

impl SplineSurface {
    /// Owns a copy of the immutable proof inputs, sharing the control net.
    pub fn from_offset(proof: OffsetSweepProof<'_>) -> Self {
        let (_, axis, _) = proof
            .definition()
            .original()
            .definition_frame()
            .expect("validated elliptic offset source is analytic");
        let axis = axis.scaled(&(Real::one() / axis.norm_squared().sqrt()));
        let axial = proof.axial().map(Clone::clone);
        let height = &axial[1] - &axial[0];
        let carrier = proof.surface().clone();
        let edge = BSplineCurve::boundary(&carrier, SplineBoundary::VMin);
        Self {
            data: Arc::new(Data {
                carrier,
                definition: proof.definition().clone(),
                axial,
                height,
                axis,
                edge,
                bound: proof.into_bound(),
                transport: None,
            }),
        }
    }
    /// Represented geometry; its finite domain is [0,4] x [0,1].
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Retained original exact signed normal-offset target.
    pub fn definition(&self) -> &EllipticOffsetDefinition {
        &self.data.definition
    }
    /// Exact lower and upper source-axis coordinates.
    pub fn axial(&self) -> [&Real; 2] {
        self.data.axial.each_ref()
    }
    /// Exact positive axial sweep length.
    pub fn height(&self) -> &Real {
        &self.data.height
    }
    /// Exact unit sweep direction; positive v follows this vector.
    pub fn axis(&self) -> &Vector3 {
        &self.data.axis
    }
    /// The shared lower edge, parameterised in the surface's u direction.
    pub fn lower_edge(&self) -> &BSplineCurve {
        &self.data.edge
    }
    /// Proved two-way representation-to-definition position bound in millimetres.
    pub fn bound(&self) -> &EvidenceBound {
        &self.data.bound
    }
    /// Positive identity of the immutable checked definition, not geometric apartness.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
    /// Source-law frame. The reference direction need not be a unit vector.
    pub fn frame(&self) -> (&Point3, &Vector3, &Vector3) {
        let (origin, _, reference) = self
            .definition()
            .original()
            .definition_frame()
            .expect("validated elliptic offset source is analytic");
        (origin, self.axis(), reference)
    }
    fn source_frame(&self) -> (&Self, Option<&Pose>) {
        match &self.data.transport {
            None => (self, None),
            Some(transport) => (&transport.source, Some(&transport.pose)),
        }
    }
    /// Rigid transport preserves regularity, material orientation and the metric
    /// bound. It does not reinterpret fitted controls as a new exact target.
    pub fn transformed(&self, pose: &Pose) -> Self {
        let (source, pose) = match &self.data.transport {
            None => (self, pose.clone()),
            Some(transport) => (&transport.source, transport.pose.then(pose)),
        };
        let inverse = pose.inverse();
        let carrier = source.carrier().transformed(&pose);
        let edge = BSplineCurve::boundary(&carrier, SplineBoundary::VMin);
        Self {
            data: Arc::new(Data {
                carrier,
                definition: source.definition().transformed(&pose),
                axial: source.data.axial.clone(),
                height: source.height().clone(),
                axis: pose.apply_vector(source.axis()),
                edge,
                bound: source.bound().clone(),
                transport: Some(Transport {
                    source: source.clone(),
                    pose,
                    inverse,
                }),
            }),
        }
    }
    /// Exact axial coordinate, relative to the retained source origin.
    pub fn axial_coordinate(&self, point: &Point3) -> Real {
        point.minus(self.frame().0).dot(self.axis())
    }
    fn project_base(&self, point: &Point3, axial: &Real) -> Point3 {
        point.plus(&self.axis().scaled(&(-(axial - &self.data.axial[0]))))
    }
    /// Exact finite incidence. A point outside either chart coordinate returns None;
    /// failed algebraic proof propagates rather than becoming non-incidence.
    pub fn parameters_at<'a>(
        &'a self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<SplinePointParameters<'a>>, Refusal> {
        dec.check_deadline()?;
        if let Some(transport) = &self.data.transport {
            return transport
                .source
                .parameters_at(&transport.inverse.apply_point(point), dec)
                .map(|parameters| {
                    parameters.map(|mut parameters| {
                        // The span witness remains in the immutable canonical frame.
                        parameters.surface = self;
                        parameters
                    })
                });
        }
        math::validate_point(point, dec)?;
        let axial = self.axial_coordinate(point);
        let v = (&axial - &self.data.axial[0]) / self.data.height.clone();
        if math::sign(&v, dec)? == Sign::Negative
            || math::sign(&(&v - &Real::one()), dec)? == Sign::Positive
        {
            return Ok(None);
        }
        let projected = self.project_base(point, &axial);
        Ok(self
            .lower_edge()
            .parameter_of(&projected, dec)?
            .map(|u| SplinePointParameters {
                surface: self,
                u,
                v,
            }))
    }

    /// Oriented side in the finite axial strip. Outside that strip there is
    /// no surface-side certificate; this is not solid containment.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        dec.check_deadline()?;
        if let Some(transport) = &self.data.transport {
            return transport
                .source
                .side_at(&transport.inverse.apply_point(point), dec);
        }
        math::validate_point(point, dec)?;
        let axial = self.axial_coordinate(point);
        if math::sign(&(&axial - &self.data.axial[0]), dec)? == Sign::Negative
            || math::sign(&(&axial - &self.data.axial[1]), dec)? == Sign::Positive
        {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "side point leaves the finite axial chart",
            ));
        }
        self.projected_side(&self.project_base(point, &axial), dec)
    }

    /// Encloses every chart preimage of a spatial box on one connected branch.
    /// Existence still requires a separate incidence proof.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        dec.check_deadline()?;
        let precision = point
            .iter()
            .map(Ball::precision)
            .max()
            .expect("three coordinates");
        if precision == 0 || precision > dec.policy().max_precision_bits {
            return Err(precision_refusal());
        }
        let bits = precision
            .max(64)
            .saturating_add(32)
            .min(dec.policy().max_precision_bits);
        if let Some(transport) = &self.data.transport {
            return transport
                .source
                .parameters_of_box(&transform_box(&transport.inverse, point, bits, dec)?, dec)
                .map(|parameters| {
                    parameters.map(|p| p.map(|value| value.with_precision(precision)))
                });
        }
        let axis = vector_box(self.axis(), bits, dec)?;
        let origin = math::point(self.frame().0, bits, dec)?;
        let displacement = std::array::from_fn(|i| &point[i] - &origin[i]);
        let axial = ball_dot(&displacement, &axis);
        let lower = dec.enclose_at(&self.data.axial[0], bits)?;
        let height = dec.enclose_at(self.height(), bits)?;
        let v = (&axial - &lower).checked_div(&height).ok_or_else(|| {
            Refusal::new(
                "numbers.divisor-not-separated",
                "spline sweep height is not separated from zero",
            )
        })?;
        let v_lo = v.lo().max(Rational::zero());
        let v_hi = v.hi().min(Rational::one());
        if v_lo > v_hi {
            return Ok(None);
        }
        let shift = &axial - &lower;
        let base: [Ball; 3] = std::array::from_fn(|i| &point[i] - &(&axis[i] * &shift));
        let curve = self.lower_edge();
        let mut row = Vec::with_capacity(4);
        let mut cover: Option<[Rational; 2]> = None;
        for (span, interval) in curve.spans() {
            dec.check_deadline()?;
            let first = span - 3;
            let controls = [
                math::point(curve.control(first), bits, dec)?,
                math::point(curve.control(first + 1), bits, dec)?,
                math::point(curve.control(first + 2), bits, dec)?,
                math::point(curve.control(first + 3), bits, dec)?,
            ];
            if (0..3).any(|i| {
                let hull = controls[0][i]
                    .hull(&controls[1][i])
                    .hull(&controls[2][i])
                    .hull(&controls[3][i]);
                base[i].hi() < hull.lo() || base[i].lo() > hull.hi()
            }) {
                continue;
            }
            let differences: [Vector3; 3] = std::array::from_fn(|i| {
                curve.control(first + i + 1).minus(curve.control(first + i))
            });
            // The opaque sweep proof establishes one open derivative half-plane
            // about the sum of these endpoint tangents on every Bezier span.
            let reference = differences[0].plus(&differences[2]);
            let reference_box = vector_box(&reference, bits, dec)?;
            let step = Real::from_rational(Rational::from_i64(3) / (interval[1] - interval[0]));
            let slopes = [
                dec.enclose_at(&(differences[0].dot(&reference) * step.clone()), bits)?,
                dec.enclose_at(&(differences[1].dot(&reference) * step.clone()), bits)?,
                dec.enclose_at(&(differences[2].dot(&reference) * step), bits)?,
            ];
            let slope = slopes[0].hull(&slopes[1]).hull(&slopes[2]);
            if slope.lo() <= Rational::zero() {
                return Err(Refusal::new(
                    "numbers.divisor-not-separated",
                    "spline inverse monotone projection is not separated",
                ));
            }
            let target = ball_dot(
                &std::array::from_fn(|i| &base[i] - &controls[0][i]),
                &reference_box,
            );
            let mut lo = interval[0].clone();
            let mut hi = interval[1].clone();
            let mut excluded = false;
            for _ in 0..bits.ilog2() + 4 {
                dec.check_deadline()?;
                let [position, tangent] =
                    curve.first_box_reusing([&lo, &hi], bits, dec, &mut row)?;
                if (0..3)
                    .any(|i| base[i].hi() < position[i].lo() || base[i].lo() > position[i].hi())
                {
                    excluded = true;
                    break;
                }
                let derivative = ball_dot(&tangent, &reference_box);
                let d_lo = derivative.lo().max(slope.lo());
                let d_hi = derivative.hi().min(slope.hi());
                if d_lo > d_hi {
                    return Err(Refusal::new(
                        "shapes.spline-inverse",
                        "inconsistent certified derivative covers",
                    ));
                }
                let middle = (&lo + &hi).mul_pow2(-1);
                let centre = curve.evaluate_box([&middle, &middle], bits, dec)?;
                let value = ball_dot(
                    &std::array::from_fn(|i| &centre[i] - &controls[0][i]),
                    &reference_box,
                );
                let correction = (&target - &value)
                    .checked_div(&ball_bounds(&d_lo, &d_hi, bits))
                    .ok_or_else(|| {
                        Refusal::new(
                            "numbers.divisor-not-separated",
                            "spline inverse derivative is not separated",
                        )
                    })?;
                let next_lo = lo.clone().max(&middle + &correction.lo());
                let next_hi = hi.clone().min(&middle + &correction.hi());
                if next_lo > next_hi {
                    excluded = true;
                    break;
                }
                if next_lo == lo && next_hi == hi {
                    break;
                }
                lo = next_lo;
                hi = next_hi;
            }
            if excluded {
                continue;
            }
            match &mut cover {
                None => cover = Some([lo, hi]),
                Some(previous) => {
                    if lo > previous[1] {
                        return Ok(None);
                    }
                    previous[1] = previous[1].clone().max(hi);
                }
            }
        }
        Ok(cover.map(|u| {
            [
                ball_bounds(&u[0], &u[1], precision),
                ball_bounds(&v_lo, &v_hi, precision),
            ]
        }))
    }
    /// Oriented side of the represented convex cross-section, after orthogonal
    /// projection. Negative is inside, zero is on the actual spline, positive
    /// is outside. This deliberately makes no finite axial-membership claim.
    /// Every negative result covers all curve spans by outward tangent supports.
    pub fn cross_section_side(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        dec.check_deadline()?;
        if let Some(transport) = &self.data.transport {
            return transport
                .source
                .cross_section_side(&transport.inverse.apply_point(point), dec);
        }
        math::validate_point(point, dec)?;
        let projected = self.project_base(point, &self.axial_coordinate(point));
        self.projected_side(&projected, dec)
    }

    fn projected_side(&self, projected: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let max_bits = dec.policy().max_precision_bits;
        if max_bits == 0 {
            return Err(precision_refusal());
        }
        if let Some(side) = filter_convex_side(self, projected, max_bits.min(128), dec)? {
            return Ok(side);
        }
        if self.lower_edge().parameter_of(projected, dec)?.is_some() {
            return Ok(Sign::Zero);
        }
        let mut queries = BTreeMap::new();
        let mut row = Vec::with_capacity(self.lower_edge().degree() + 1);
        let mut pending = Vec::new();
        let zero = Rational::zero();
        let two = Rational::from_i64(2);
        for (_, interval) in self.lower_edge().spans() {
            pending.push((interval[0].clone(), interval[1].clone(), 0u32));
            while let Some((lo, hi, depth)) = pending.pop() {
                dec.check_deadline()?;
                if depth >= max_bits {
                    return Err(precision_refusal());
                }
                let bits = (128u32.saturating_mul(1 + depth / 64)).min(max_bits);
                let query = match queries.entry(bits) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => entry.insert(HalfPlaneQuery {
                        point: math::point(projected, bits, dec)?,
                        axis: [
                            dec.enclose_at(&self.axis().x, bits)?,
                            dec.enclose_at(&self.axis().y, bits)?,
                            dec.enclose_at(&self.axis().z, bits)?,
                        ],
                    }),
                };
                let [position, tangent] =
                    self.lower_edge()
                        .first_box_reusing([&lo, &hi], bits, dec, &mut row)?;
                let normal: [Ball; 3] = std::array::from_fn(|i| {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    &tangent[j] * &query.axis[k] - &tangent[k] * &query.axis[j]
                });
                let terms: [Ball; 3] =
                    std::array::from_fn(|i| &normal[i] * &(&query.point[i] - &position[i]));
                let value = &(&terms[0] + &terms[1]) + &terms[2];
                if value.hi() < zero {
                    continue;
                }
                if value.lo() > zero {
                    return Ok(Sign::Positive);
                }
                let middle = (&lo + &hi) / two.clone();
                pending.push((middle.clone(), hi, depth + 1));
                pending.push((lo, middle, depth + 1));
            }
        }
        Ok(Sign::Negative)
    }
}
struct HalfPlaneQuery {
    point: [Ball; 3],
    axis: [Ball; 3],
}
fn precision_refusal() -> Refusal {
    Refusal::new(
        "numbers.precision-limit",
        "finite spline side could not be separated within the caller's precision policy",
    )
}

// The checked carrier is a closed convex composite cubic Bezier curve with
// unit weights. For a span C and unit sweep axis a, h=(C' cross a).(p-C)
// is its outward tangent half-plane polynomial. Strictly negative h on the
// complete loop proves inside; one positive sampled value proves outside.
// Unlike separate interval boxes for C and C', Bernstein coefficients retain
// their common parameter. An uncertain filter never proves coincidence.
fn filter_convex_side(
    surface: &SplineSurface,
    point: &Point3,
    bits: u32,
    dec: &Decider,
) -> Result<Option<Sign>, Refusal> {
    let point = math::point(point, bits, dec)?;
    let axis = [
        dec.enclose_at(&surface.axis().x, bits)?,
        dec.enclose_at(&surface.axis().y, bits)?,
        dec.enclose_at(&surface.axis().z, bits)?,
    ];
    let b2 = [1, 2, 1];
    let b3 = [1, 3, 3, 1];
    let b5 = [1, 5, 10, 10, 5, 1];
    let weights: [[Ball; 4]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let value = Rational::from_i64(b2[i] * b3[j]) / Rational::from_i64(b5[i + j]);
            Ball::from_rational(&value, bits)
        })
    });
    let half = Ball::from_rational(&(Rational::one() / Rational::from_i64(2)), bits);
    let zero = Rational::zero();
    // This is a bounded proof filter, not a modelling or contact tolerance.
    // The original exact path handles anything not settled by these cuts.
    let depth_limit = 24.min(bits / 2);
    let mut cells = 0;
    let mut pending = Vec::with_capacity(25);
    let curve = surface.lower_edge();
    for (span, _) in curve.spans() {
        dec.check_deadline()?;
        let first = span - 3;
        let controls = [
            math::point(curve.control(first), bits, dec)?,
            math::point(curve.control(first + 1), bits, dec)?,
            math::point(curve.control(first + 2), bits, dec)?,
            math::point(curve.control(first + 3), bits, dec)?,
        ];
        // Omit the common positive factor 3/(span_end-span_start).
        let normal: [[Ball; 3]; 3] = std::array::from_fn(|i| {
            let tangent: [Ball; 3] = std::array::from_fn(|k| &controls[i + 1][k] - &controls[i][k]);
            std::array::from_fn(|k| {
                let j = (k + 1) % 3;
                let l = (k + 2) % 3;
                &tangent[j] * &axis[l] - &tangent[l] * &axis[j]
            })
        });
        let delta: [[Ball; 3]; 4] =
            std::array::from_fn(|i| std::array::from_fn(|k| &point[k] - &controls[i][k]));
        let mut coefficients: [Ball; 6] = std::array::from_fn(|_| Ball::from_rational(&zero, bits));
        for i in 0..3 {
            for (j, delta) in delta.iter().enumerate() {
                let dot = &(&normal[i][0] * &delta[0])
                    + &(&normal[i][1] * &delta[1])
                    + &normal[i][2] * &delta[2];
                coefficients[i + j] = &coefficients[i + j] + &(&dot * &weights[i][j]);
            }
        }
        pending.push((coefficients, 0));
        while let Some((coefficients, depth)) = pending.pop() {
            dec.check_deadline()?;
            cells += 1;
            if cells > 1024 {
                return Ok(None);
            }
            if coefficients[0].lo() > zero || coefficients[5].lo() > zero {
                return Ok(Some(Sign::Positive));
            }
            if coefficients.iter().all(|value| value.hi() < zero) {
                continue;
            }
            if depth >= depth_limit {
                return Ok(None);
            }
            let [left, right] = split_half_planes(coefficients, &half);
            pending.push((right, depth + 1));
            pending.push((left, depth + 1));
        }
    }
    Ok(Some(Sign::Negative))
}

fn split_half_planes(mut row: [Ball; 6], half: &Ball) -> [[Ball; 6]; 2] {
    let mut left: [Option<Ball>; 6] = [const { None }; 6];
    let mut right: [Option<Ball>; 6] = [const { None }; 6];
    left[0] = Some(row[0].clone());
    right[5] = Some(row[5].clone());
    for level in 1..6 {
        for i in 0..6 - level {
            row[i] = &(&row[i] + &row[i + 1]) * half;
        }
        left[level] = Some(row[0].clone());
        right[5 - level] = Some(row[5 - level].clone());
    }
    [
        left.map(|value| value.expect("all left Bernstein coefficients assigned")),
        right.map(|value| value.expect("all right Bernstein coefficients assigned")),
    ]
}

fn ball_bounds(lo: &Rational, hi: &Rational, precision: u32) -> Ball {
    Ball::from_rational(lo, precision).hull(&Ball::from_rational(hi, precision))
}
fn vector_box(vector: &Vector3, precision: u32, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_at(&vector.x, precision)?,
        dec.enclose_at(&vector.y, precision)?,
        dec.enclose_at(&vector.z, precision)?,
    ])
}
fn ball_dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&a[0] * &b[0]) + &(&a[1] * &b[1]) + &a[2] * &b[2]
}
fn transform_box(
    pose: &Pose,
    point: &[Ball; 3],
    precision: u32,
    dec: &Decider,
) -> Result<[Ball; 3], Refusal> {
    let row = |i: usize| -> Result<Ball, Refusal> {
        let matrix = pose.rotation.matrix();
        let coefficients = [
            dec.enclose_at(&matrix[i][0], precision)?,
            dec.enclose_at(&matrix[i][1], precision)?,
            dec.enclose_at(&matrix[i][2], precision)?,
        ];
        Ok(ball_dot(&coefficients, point)
            + dec.enclose_at(
                [
                    &pose.translation.x,
                    &pose.translation.y,
                    &pose.translation.z,
                ][i],
                precision,
            )?)
    };
    Ok([row(0)?, row(1)?, row(2)?])
}
