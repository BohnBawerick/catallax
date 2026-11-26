//! Whole-edge distance certificates, independent of sampled admission.

#[cfg(test)]
mod tests;

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};

use catallax_numbers::{Ball, Decider, FactKind, Rational, Refusal};
use catallax_shapes::{BSplineCurve, Curve, GeomKey, ParamKind, Surface};

use super::edge_gap::{face_copy, traversal};
use super::gap::FaceCopyRef;
use super::{BITS, NumericalBody, NumericalEdge, NumericalFace, NumericalTrim};
use crate::{BodyId, EdgeId, Loop, Model, Tolerance};

/// The edge's identifier in the part passed to the proof request.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdgeProofKey {
    /// Handle in the exact part's model.
    Exact(EdgeId),
    /// Stable edge key in the numerical definition.
    Numerical(String),
}

/// Why an edge did not receive a whole-edge certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeProofFailure {
    /// A certified lower bound at some parameter exceeded the requested limit.
    AboveLimit,
    /// The edge used its share of the requested work cap.
    WorkCap,
    /// The caller cancelled before this edge was proved.
    Cancelled,
    /// There is no supported face copy or exact incidence proof for a use.
    Unsupported,
    /// A range could not be isolated at the working precision and subdivision depth.
    Unresolved,
}

/// Result for an entire edge, including all its face uses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdgeProofOutcome {
    /// A whole-edge certificate within the requested gap limit.
    Proved {
        /// Outward upper bound on the maximum gap, in millimetres.
        upper_mm: Rational,
    },
    /// No bound within the requested limit was proved for the whole edge.
    CouldNotProve(EdgeProofFailure),
}

/// A proof result tied to one edge of the input part.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EdgeProof {
    /// Edge identifier in this part.
    pub edge: EdgeProofKey,
    /// Result after checking the edge's face uses.
    pub outcome: EdgeProofOutcome,
}

/// Request-local results. No part, tolerance or admission evidence is modified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EdgeProofReport {
    /// One row per edge in the exact model or numerical definition order.
    pub edges: Vec<EdgeProof>,
    /// Number of work units spent across the part.
    pub work: u64,
}

#[derive(Clone, Copy)]
struct Work<'a> {
    left: u64,
    cancel: &'a AtomicBool,
}

impl Work<'_> {
    fn check_cancel(&self) -> Result<(), EdgeProofFailure> {
        if self.cancel.load(Ordering::Relaxed) {
            return Err(EdgeProofFailure::Cancelled);
        }
        Ok(())
    }

    fn spend(&mut self) -> Result<(), EdgeProofFailure> {
        self.check_cancel()?;
        self.left = self.left.checked_sub(1).ok_or(EdgeProofFailure::WorkCap)?;
        Ok(())
    }

    // A hard use cannot consume the whole part's allowance before other edges are tried.
    fn edge(
        &mut self,
        remaining_edges: usize,
        prove: impl FnOnce(&mut Self) -> Result<Rational, EdgeProofFailure>,
    ) -> EdgeProofOutcome {
        let slice = self.left.div_ceil(remaining_edges as u64);
        let mut work = Self {
            left: slice,
            ..*self
        };
        let result = work.check_cancel().and_then(|()| prove(&mut work));
        self.left -= slice - work.left;
        match work.check_cancel().and(result) {
            Ok(upper_mm) => EdgeProofOutcome::Proved { upper_mm },
            Err(reason) => EdgeProofOutcome::CouldNotProve(reason),
        }
    }
}

/// Proves zero gap from whole-carrier incidence for the edges of an exact part.
/// One work unit checks one edge use. A stored tolerance is not a proof.
pub fn prove_exact_edges(
    model: &Model,
    body: BodyId,
    work_cap: u64,
    cancel: &AtomicBool,
) -> EdgeProofReport {
    let edges = model.edges_of_body(body);
    let mut uses = Vec::new();
    for face in model.faces_of_body(body) {
        for boundary in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*boundary) {
                uses.extend(coedges.iter().map(|id| (model.coedge(*id).edge, face)));
            }
        }
    }
    uses.sort_unstable();
    let mut work = Work {
        left: work_cap,
        cancel,
    };
    let mut report = Vec::with_capacity(edges.len());
    for (i, &id) in edges.iter().enumerate() {
        let begin = uses.partition_point(|&(edge, _)| edge < id);
        let end = uses.partition_point(|&(edge, _)| edge <= id);
        let adjacent = &uses[begin..end];
        let outcome = work.edge(edges.len() - i, |work| {
            let edge = model.edge(id);
            if adjacent.is_empty() || edge.tolerance.is_some() {
                return Err(EdgeProofFailure::Unsupported);
            }
            for &(_, face) in adjacent {
                work.spend()?;
                let face = model.face(face);
                if face.tolerance.is_some()
                    || !model.geom().facts().holds(
                        FactKind::Incident,
                        GeomKey::Curve(edge.curve),
                        GeomKey::Surface(face.surface),
                    )
                {
                    return Err(EdgeProofFailure::Unsupported);
                }
            }
            Ok(Rational::zero())
        });
        report.push(EdgeProof {
            edge: EdgeProofKey::Exact(id),
            outcome,
        });
    }
    EdgeProofReport {
        edges: report,
        work: work_cap - work.left,
    }
}

/// Proves a distance upper bound no greater than `limit` for each numerical edge.
/// One work unit visits one parameter interval or projected trim endpoint.
/// Cancellation is observed between visits. Spline range evaluations within a
/// visit are synchronous. The result never relies on admission samples or tags.
pub fn prove_numerical_edges(
    body: &NumericalBody,
    limit: &Tolerance,
    work_cap: u64,
    cancel: &AtomicBool,
) -> EdgeProofReport {
    let definition = body.definition();
    let native: BTreeMap<_, _> = body
        .native_trims()
        .iter()
        .map(|trim| {
            (
                (trim.face.as_str(), trim.boundary, trim.use_index),
                &trim.curve,
            )
        })
        .collect();
    let mut uses = Vec::new();
    for (face_index, face) in definition.faces.iter().enumerate() {
        for (boundary, trims) in face.loops.iter().enumerate() {
            for (slot, trim) in trims.iter().enumerate() {
                uses.push((trim.edge, face_index, boundary, slot));
            }
        }
    }
    uses.sort_unstable();
    let mut uses = uses.as_slice();
    let dec = Decider::default();
    let limit_sq = limit.mm() * limit.mm();
    let mut work = Work {
        left: work_cap,
        cancel,
    };
    let mut report = Vec::with_capacity(definition.edges.len());
    for (i, edge) in definition.edges.iter().enumerate() {
        let (adjacent, rest) = uses.split_at(uses.partition_point(|&(edge, ..)| edge <= i));
        uses = rest;
        let outcome = work.edge(definition.edges.len() - i, |work| {
            if adjacent.is_empty() {
                return Err(EdgeProofFailure::Unsupported);
            }
            let mut upper = Rational::zero();
            for &(_, face, boundary, slot) in adjacent {
                work.check_cancel()?;
                let face = &definition.faces[face];
                let trim = &face.loops[boundary][slot];
                let native = native.get(&(face.key.as_str(), boundary, slot)).copied();
                upper = upper.max(prove_use(edge, face, trim, native, &limit_sq, work, &dec)?);
            }
            let rounded = ball(&upper).sqrt().expect("nonnegative squared gap").hi();
            Ok(rounded.min(limit.mm().clone()))
        });
        report.push(EdgeProof {
            edge: EdgeProofKey::Numerical(edge.key.clone()),
            outcome,
        });
    }
    EdgeProofReport {
        edges: report,
        work: work_cap - work.left,
    }
}

fn ball(q: &Rational) -> Ball {
    Ball::from_rational(q, BITS)
}

fn range(bounds: [&Rational; 2]) -> Ball {
    ball(bounds[0]).hull(&ball(bounds[1]))
}

fn vector(values: [&catallax_numbers::Real; 3], dec: &Decider) -> Result<[Ball; 3], Refusal> {
    let [x, y, z] = values;
    Ok([
        dec.enclose_at(x, BITS)?,
        dec.enclose_at(y, BITS)?,
        dec.enclose_at(z, BITS)?,
    ])
}

fn delta(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}

fn norm_squared(a: &[Ball; 3]) -> Ball {
    a[0].square() + a[1].square() + a[2].square()
}

fn position(curve: &Curve, t: [&Rational; 2], dec: &Decider) -> Result<[Ball; 3], Refusal> {
    match curve {
        Curve::BSpline(spline) => spline.evaluate_box(t, BITS, dec),
        _ => curve.evaluate(&range(t), dec)?.ok_or_else(|| {
            Refusal::new("topology.edge-proof-curve", "curve has no range evaluator")
        }),
    }
}

type UvBox = [[Rational; 2]; 2];

fn chart(surface: &Surface, uv: &UvBox, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    if let Surface::BSpline(s) = surface {
        return s
            .carrier()
            .evaluate_box(uv[0].each_ref(), uv[1].each_ref(), BITS, dec);
    }
    check_chart_domain(surface, uv)?;
    surface.evaluate(&range(uv[0].each_ref()), &range(uv[1].each_ref()), dec)
}

fn check_chart_domain(surface: &Surface, uv: &UvBox) -> Result<(), Refusal> {
    for ([lo, hi], domain) in uv.iter().zip(surface.domain()) {
        let inside = match domain {
            ParamKind::Ray => lo >= &Rational::zero(),
            ParamKind::Latitude => {
                let pole = &Ball::pi(BITS).lo() / &Rational::from(2);
                lo >= &(-&pole) && hi <= &pole
            }
            _ => true,
        };
        if !inside {
            return Err(Refusal::new(
                "topology.edge-proof-domain",
                "box may leave the chart",
            ));
        }
    }
    Ok(())
}

enum ParameterMap<'a> {
    Linear {
        start: &'a Rational,
        uv: &'a [Rational; 2],
        slope: [Rational; 2],
    },
    Native(&'a BSplineCurve),
    Projected {
        origin: [Ball; 3],
        normal: [Ball; 3],
    },
}

impl ParameterMap<'_> {
    fn uv(&self, t: [&Rational; 2], dec: &Decider) -> Result<UvBox, Refusal> {
        Ok(match self {
            Self::Linear { start, uv, slope } => {
                let x = t.map(|x| x - *start);
                std::array::from_fn(|i| {
                    let mut limits = x.each_ref().map(|x| &uv[i] + &(&slope[i] * x));
                    limits.sort();
                    limits
                })
            }
            Self::Native(spline) => {
                let [u, v, _] = spline.evaluate_box(t, BITS, dec)?;
                [[u.lo(), u.hi()], [v.lo(), v.hi()]]
            }
            Self::Projected { .. } => unreachable!(),
        })
    }

    fn slope(&self, t: [&Rational; 2], dec: &Decider) -> Result<[Ball; 2], Refusal> {
        Ok(match self {
            Self::Linear { slope, .. } => slope.each_ref().map(ball),
            Self::Native(spline) => {
                let [u, v, _] = spline.jet_box(t, BITS, dec)?.first;
                [u, v]
            }
            Self::Projected { .. } => unreachable!(),
        })
    }
}

struct Copies<'a> {
    curve: &'a Curve,
    surface: &'a Surface,
    map: ParameterMap<'a>,
    dec: &'a Decider,
}

fn normal_residual(vector: &[Ball; 3], normal: &[Ball; 3]) -> [Ball; 3] {
    let dot = &vector[0] * &normal[0] + &vector[1] * &normal[1] + &vector[2] * &normal[2];
    let zero = ball(&Rational::zero());
    [dot, zero.clone(), zero]
}

impl Copies<'_> {
    fn residual(&self, t: [&Rational; 2]) -> Result<[Ball; 3], Refusal> {
        let c = position(self.curve, t, self.dec)?;
        if let ParameterMap::Projected { origin, normal } = &self.map {
            // For the plane's orthogonal projection the residual is the signed normal distance.
            return Ok(normal_residual(&delta(&c, origin), normal));
        }
        let uv = self.map.uv(t, self.dec)?;
        Ok(delta(&c, &chart(self.surface, &uv, self.dec)?))
    }

    fn derivative(&self, t: [&Rational; 2]) -> Result<Option<[Ball; 3]>, Refusal> {
        let tangent = match self.curve {
            Curve::Line(line) => vector([&line.dir.x, &line.dir.y, &line.dir.z], self.dec)?,
            Curve::BSpline(spline) => spline.jet_box(t, BITS, self.dec)?.first,
            _ => return Ok(None),
        };
        if let ParameterMap::Projected { normal, .. } = &self.map {
            return Ok(Some(normal_residual(&tangent, normal)));
        }
        let uv = self.map.uv(t, self.dec)?;
        let [u, v] = self.map.slope(t, self.dec)?;
        let [su, sv] = match self.surface {
            Surface::BSpline(s) => {
                s.carrier()
                    .jet_box(uv[0].each_ref(), uv[1].each_ref(), BITS, self.dec)?
                    .first
            }
            s => {
                check_chart_domain(s, &uv)?;
                s.derivatives(&range(uv[0].each_ref()), &range(uv[1].each_ref()), self.dec)?
            }
        };
        Ok(Some(std::array::from_fn(|i| {
            &tangent[i] - &(&su[i] * &u + &sv[i] * &v)
        })))
    }
}

fn prove_use(
    edge: &NumericalEdge,
    face: &NumericalFace,
    trim: &NumericalTrim,
    native: Option<&BSplineCurve>,
    limit_sq: &Rational,
    work: &mut Work<'_>,
    dec: &Decider,
) -> Result<Rational, EdgeProofFailure> {
    if !matches!(
        edge.curve,
        Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) | Curve::BSpline(_)
    ) {
        return Err(EdgeProofFailure::Unsupported);
    }
    let [first, last] = traversal(trim, edge);
    let mut upper = Rational::zero();
    let map = match face_copy(face, trim, edge, native) {
        FaceCopyRef::Linear(uv) => {
            let span = last - first;
            if span == Rational::zero() {
                return Err(EdgeProofFailure::Unsupported);
            }
            ParameterMap::Linear {
                start: first,
                uv: &uv[0],
                slope: std::array::from_fn(|i| &(&uv[1][i] - &uv[0][i]) / &span),
            }
        }
        FaceCopyRef::Native(curve) => ParameterMap::Native(curve),
        FaceCopyRef::Projected => {
            let Surface::Plane(plane) = &face.surface else {
                return Err(EdgeProofFailure::Unsupported);
            };
            // The projected trim retains two declared endpoint copies as well.
            for (t, uv) in [first, last].into_iter().zip(&trim.uv) {
                work.spend()?;
                let residual = position(&edge.curve, [t, t], dec)
                    .and_then(|c| {
                        let uv = uv.each_ref().map(|x| [x.clone(), x.clone()]);
                        Ok(delta(&c, &chart(&face.surface, &uv, dec)?))
                    })
                    .map_err(|_| EdgeProofFailure::Unresolved)?;
                let squared = norm_squared(&residual);
                if squared.lo() > *limit_sq {
                    return Err(EdgeProofFailure::AboveLimit);
                }
                let bound = squared.hi();
                if bound > *limit_sq {
                    return Err(EdgeProofFailure::Unresolved);
                }
                upper = upper.max(bound);
            }
            let normal = &plane.normal;
            let origin = &plane.origin;
            let components = vector([&normal.x, &normal.y, &normal.z], dec)
                .map_err(|_| EdgeProofFailure::Unresolved)?;
            let length = dec
                .enclose_at(&normal.norm_squared().sqrt(), BITS)
                .ok()
                .and_then(|b| b.recip())
                .ok_or(EdgeProofFailure::Unresolved)?;
            ParameterMap::Projected {
                origin: vector([&origin.x, &origin.y, &origin.z], dec)
                    .map_err(|_| EdgeProofFailure::Unresolved)?,
                normal: components.map(|c| &c * &length),
            }
        }
    };
    let copies = Copies {
        curve: &edge.curve,
        surface: &face.surface,
        map,
        dec,
    };
    let (lo, hi) = (first.min(last), first.max(last));
    let bound = bisect(&copies, lo, hi, limit_sq, work, 0)?;
    Ok(upper.max(bound))
}

// Each accepted leaf covers a closed interval. The maximum of all accepted
// leaves bounds the whole use; no point test can accept a non-point interval.
fn bisect(
    copies: &Copies<'_>,
    lo: &Rational,
    hi: &Rational,
    limit_sq: &Rational,
    work: &mut Work<'_>,
    depth: u8,
) -> Result<Rational, EdgeProofFailure> {
    work.spend()?;
    if let Ok(residual) = copies.residual([lo, hi]) {
        let bound = norm_squared(&residual).hi();
        if bound <= *limit_sq {
            return Ok(bound);
        }
    }
    let mid = &(lo + hi) / &Rational::from(2);
    if let Ok(center) = copies.residual([&mid, &mid]) {
        if norm_squared(&center).lo() > *limit_sq {
            return Err(EdgeProofFailure::AboveLimit);
        }
        // The residual is continuous and piecewise smooth. Spline jet boxes
        // include both knot sides, so the mean-value enclosure holds across knots.
        if let Ok(Some(derivative)) = copies.derivative([lo, hi]) {
            let offset = range([&(lo - &mid), &(hi - &mid)]);
            let enclosure = std::array::from_fn(|i| &center[i] + &(&offset * &derivative[i]));
            let bound = norm_squared(&enclosure).hi();
            if bound <= *limit_sq {
                return Ok(bound);
            }
        }
    }
    if depth == 32 || lo == hi {
        return Err(EdgeProofFailure::Unresolved);
    }
    let left = bisect(copies, lo, &mid, limit_sq, work, depth + 1)?;
    let right = bisect(copies, &mid, hi, limit_sq, work, depth + 1)?;
    Ok(left.max(right))
}
