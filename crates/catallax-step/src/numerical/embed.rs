//! Numerical embedding of a retained exact B-rep whose geometry has no exact
//! STEP form.
//!
//! The topology, planes and native spline carriers are kept. Each exact
//! patch-section or patch-pair edge becomes a fitted cubic spline, and every
//! edge use on a spline face receives a fitted native UV curve at that edge's
//! own parameter. Binary64 arithmetic only proposes these fits. Numerical
//! admission rechecks every edge, trim and loop at 128 bits under the stated
//! tolerance, and the numerical writer rechecks its decimal payload. The
//! result is approximate; no exact answer or complete body bound follows.
//!
//! This follows OCCT's face/face section, which approximates each section as
//! one 3D B-spline and one 2D B-spline on each face (`IntTools_FaceFace` with
//! curve and both pcurve approximations enabled) before the Boolean stores
//! them on its section edges.

use super::*;
use catallax_numbers::Ball;
use catallax_shapes::{BSplineCurve, Plane, SplineSupport};
use catallax_topology::numerical::{
    NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalOptions,
    NumericalTrim, NumericalVertex,
};
use catallax_topology::{BodyId, EdgeId, Loop, Model, VertexId};
use std::collections::{BTreeMap, BTreeSet};

mod nurbs;
use nurbs::{Curve64, Surface64, concatenate, distance, dot, interpolate, lerp, sub};

const BITS: u32 = 128;
const SAMPLE_BITS: u32 = 64;
// Samples per traced piece or per native edge, doubled until the binary64
// fit meets a quarter of the tolerance. Admission samples every span at its
// eighths, so the fit is checked at the same places.
const DENSITIES: [usize; 6] = [4, 8, 16, 32, 64, 128];
// Admission checks every span of every loop pairwise, so fits stay small.
const MAX_SPANS: usize = 256;
const CHECK_SHARES: [f64; 7] = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875];

fn unsupported(detail: &str) -> Refusal {
    bad("unsupported", detail)
}

fn fit_failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("step.numerical-fit", detail)
}

fn f64_point(p: &Point3) -> [f64; 3] {
    [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| fit_failure("a fitted value is not finite"))
}

fn scalar(value: &Real, dec: &Decider) -> Result<Rational, Refusal> {
    Ok(dec.enclose_at(value, BITS)?.mid())
}

/// Whether every retained carrier of `body` has a numerical embedding here:
/// one closed shell of planes and retained rational patches or atlases, at
/// most one hole per face, and finite line, native spline, patch-section or
/// patch-pair edges. This is a cheap structural precheck, not admission.
pub fn embeds(model: &Model, body: BodyId) -> bool {
    let counts = model.counts(body);
    if counts.lumps != 1 || counts.shells != 1 {
        return false;
    }
    for face in model.faces_of_body(body) {
        let face = model.face(face);
        if !matches!(
            model.geom().surface(face.surface),
            Surface::Plane(_)
                | Surface::BSpline(
                    SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_)
                )
        ) || face.loops.is_empty()
            || face.loops.len() > 2
        {
            return false;
        }
        let mut seen = BTreeSet::new();
        for lp in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                return false;
            };
            for coedge in uses {
                if !seen.insert(model.coedge(*coedge).edge) {
                    return false;
                }
            }
        }
    }
    model.edges_of_body(body).into_iter().all(|id| {
        let edge = model.edge(id);
        edge.start.is_some()
            && edge.end.is_some()
            && matches!(
                model.geom().curve(edge.curve),
                Curve::Line(_) | Curve::BSpline(_) | Curve::PatchSection(_) | Curve::PatchPair(_)
            )
    })
}

struct SplineFace {
    surface: Surface64,
    // Binary64 copy of `bounds`, where inversion searches and snaps.
    region: [[f64; 2]; 2],
    // Exact parameter box that admitted trims stay inside: the patch domain,
    // or the box round an atlas face's charts. Its sides need not be binary64
    // values, so trims on a side take the exact side, never its binary64 copy.
    bounds: [[Rational; 2]; 2],
}

impl SplineFace {
    fn new(surface: Surface64, bounds: [[Rational; 2]; 2]) -> Self {
        Self {
            surface,
            region: bounds
                .each_ref()
                .map(|side| side.each_ref().map(Rational::to_f64_nearest)),
            bounds,
        }
    }

    // The exact value of a proposed parameter: a side within binary64 noise
    // is that exact side, and any other value is kept inside the box.
    fn exact(&self, axis: usize, value: f64) -> Result<Rational, Refusal> {
        let [lo, hi] = self.region[axis];
        let noise = 1e-12 * (hi - lo);
        if (value - lo).abs() <= noise {
            return Ok(self.bounds[axis][0].clone());
        }
        if (value - hi).abs() <= noise {
            return Ok(self.bounds[axis][1].clone());
        }
        Ok(exact(value)?
            .max(self.bounds[axis][0].clone())
            .min(self.bounds[axis][1].clone()))
    }
}

enum Carrier {
    Plane { origin: [f64; 3], normal: [f64; 3] },
    Spline(Box<SplineFace>),
}

fn carrier(surface: &Surface, dec: &mut Decider) -> Result<Carrier, Refusal> {
    // The exact box round a nonempty set of rectangles.
    let bounds = |domains: &mut dyn Iterator<Item = &[[Rational; 2]; 2]>| {
        let mut bounds: Option<[[Rational; 2]; 2]> = None;
        for domain in domains {
            bounds = Some(match bounds {
                None => domain.clone(),
                Some(box_) => std::array::from_fn(|axis| {
                    [
                        box_[axis][0].clone().min(domain[axis][0].clone()),
                        box_[axis][1].clone().max(domain[axis][1].clone()),
                    ]
                }),
            });
        }
        bounds.expect("a face has at least one chart")
    };
    Ok(match surface {
        Surface::Plane(Plane { origin, normal, .. }) => {
            let n = f64_point(&Point3::new(
                normal.x.clone(),
                normal.y.clone(),
                normal.z.clone(),
            ));
            let length = dot(n, n).sqrt();
            Carrier::Plane {
                origin: f64_point(origin),
                normal: n.map(|c| c / length),
            }
        }
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
            Carrier::Spline(Box::new(SplineFace::new(
                Surface64::new(patch.carrier(), dec)?,
                patch.domain().clone(),
            )))
        }
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
            dec.charge_work(atlas.charts().len() as u64)?;
            // A logical face covers only its own charts of the source; the
            // other rows of a closed source are a different face, and
            // admission keeps every trim on these charts.
            let bounds = bounds(&mut atlas.charts().iter().map(|chart| chart.patch().domain()));
            Carrier::Spline(Box::new(SplineFace::new(
                Surface64::new(atlas.source(), dec)?,
                bounds,
            )))
        }
        _ => {
            return Err(unsupported(
                "numerical embedding needs planes or rational spline faces",
            ));
        }
    })
}

#[derive(Clone)]
struct Uv {
    exact: [Rational; 2],
    approx: [f64; 2],
}

enum Native {
    Line { origin: [f64; 3], dir: [f64; 3] },
    Spline(Curve64),
}

impl Native {
    fn point(&self, t: f64) -> [f64; 3] {
        match self {
            Self::Line { origin, dir } => std::array::from_fn(|c| origin[c] + t * dir[c]),
            Self::Spline(curve) => curve.point(t),
        }
    }
}

struct Tolerances {
    sample: f64,
    fit: f64,
}

fn rational_curve<const D: usize>(
    knots: &[f64],
    controls: &[[f64; D]],
    domain: [&Rational; 2],
    ends: [[Rational; 3]; 2],
    face: Option<&SplineFace>,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    dec.charge_work((knots.len() + D * controls.len()) as u64)?;
    let mut exact_knots = Vec::with_capacity(knots.len());
    for (i, knot) in knots.iter().enumerate() {
        exact_knots.push(if i < 4 {
            domain[0].clone()
        } else if i + 4 >= knots.len() {
            domain[1].clone()
        } else {
            exact(*knot)?
        });
    }
    for pair in exact_knots.windows(2) {
        if pair[0] > pair[1] {
            return Err(fit_failure("fitted knots are not ordered"));
        }
    }
    let last = controls.len() - 1;
    let mut points = Vec::with_capacity(controls.len());
    for (i, control) in controls.iter().enumerate() {
        let xyz = if i == 0 {
            ends[0].clone()
        } else if i == last {
            ends[1].clone()
        } else {
            let mut xyz = [Rational::zero(), Rational::zero(), Rational::zero()];
            for c in 0..D {
                // A UV control a rounding step outside the face's exact box
                // would move the trim off its domain; the curve moves by that
                // step at most.
                xyz[c] = match face {
                    Some(face) => face.exact(c, control[c])?,
                    None => exact(control[c])?,
                };
            }
            xyz
        };
        let [x, y, z] = xyz;
        points.push(Point3::new(x.into(), y.into(), z.into()));
    }
    let weights = vec![Real::one(); points.len()];
    BSplineCurve::new(3, exact_knots, points, weights, dec)
}

fn uv_end(uv: &Uv) -> [Rational; 3] {
    [uv.exact[0].clone(), uv.exact[1].clone(), Rational::zero()]
}

// A native UV curve for a line or native spline edge on a spline face, at
// that edge's own parameter over `interval`, ascending.
fn native_uv(
    face: &SplineFace,
    edge: &Native,
    interval: [&Rational; 2],
    ends: [&Uv; 2],
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let range = interval.map(Rational::to_f64_nearest);
    for density in DENSITIES {
        let parameters: Vec<f64> = (0..=density)
            .map(|j| lerp(range, j as f64 / density as f64))
            .collect();
        let mut uvs = Vec::with_capacity(parameters.len());
        uvs.push(ends[0].approx);
        let mut sampled = true;
        for t in &parameters[1..density] {
            let point = edge.point(*t);
            let uv = face
                .surface
                .invert(point, face.region, uvs.last().copied(), dec)?;
            if distance(face.surface.point(uv), point) > tolerances.sample {
                sampled = false;
                break;
            }
            uvs.push(uv);
        }
        if !sampled {
            continue;
        }
        uvs.push(ends[1].approx);
        let (knots, controls) = interpolate(&parameters, &uvs, dec)?;
        let mut fits = true;
        'spans: for pair in parameters.windows(2) {
            for share in CHECK_SHARES {
                dec.charge_work(1)?;
                let t = lerp([pair[0], pair[1]], share);
                let uv = nurbs::evaluate(&knots, &controls, t);
                if distance(face.surface.point(uv), edge.point(t)) > tolerances.fit {
                    fits = false;
                    break 'spans;
                }
            }
        }
        if fits {
            return rational_curve(
                &knots,
                &controls,
                interval,
                ends.map(uv_end),
                Some(face),
                dec,
            );
        }
    }
    Err(fit_failure(
        "a native edge has no fitted UV curve within the tolerance",
    ))
}

// Interpolates each piece between consecutive `breaks` and joins them.
fn piecewise<const D: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    breaks: &[usize],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    let mut pieces = Vec::with_capacity(breaks.len() - 1);
    for pair in breaks.windows(2) {
        pieces.push(interpolate(
            &parameters[pair[0]..=pair[1]],
            &values[pair[0]..=pair[1]],
            dec,
        )?);
    }
    Ok(concatenate(pieces))
}

struct FittedUse<'a> {
    face: &'a SplineFace,
    ends: [&'a Uv; 2],
}

// One fitted cubic for an exact traced edge, from its start to its end vertex
// over [0, 1], and one native UV curve per spline face at the same parameter.
#[allow(clippy::too_many_arguments)]
fn fitted_edge(
    model: &Model,
    id: EdgeId,
    uses: &[FittedUse<'_>],
    planes: &[([f64; 3], [f64; 3])],
    ends: [(&[Rational; 3], [f64; 3]); 2],
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<(BSplineCurve, Vec<BSplineCurve>), Refusal> {
    let edge = model.edge(id);
    let start = model.vertex(edge.start.expect("finite edge")).point;
    let end = model.vertex(edge.end.expect("finite edge")).point;
    let forward = edge.sense == Sense::Same;
    let (from, to) = if forward { (start, end) } else { (end, start) };
    let trace =
        catallax_intersect::trace(model.geom(), edge.curve, Some(from), Some(to), dec)?.value;
    let mut pieces: Vec<_> = trace.pieces().iter().collect();
    if !forward {
        pieces.reverse();
    }
    let at = |share: f64| if forward { share } else { 1.0 - share };
    let evaluate = |piece: &catallax_intersect::TracePiece,
                    share: f64,
                    dec: &mut Decider|
     -> Result<[f64; 3], Refusal> {
        dec.charge_work(1)?;
        let t = Rational::from_f64_exact(at(share)).expect("finite share");
        let value = piece.evaluate(&Ball::from_rational(&t, SAMPLE_BITS), dec)?;
        Ok(value.map(|c| c.mid().to_f64_nearest()))
    };
    let first = evaluate(pieces[0], 0.0, dec)?;
    if distance(first, ends[0].1) > tolerances.sample {
        return Err(fit_failure("traced edge does not start at its vertex"));
    }
    let mut corners = vec![ends[0].1];
    for (index, piece) in pieces.iter().enumerate() {
        let corner = evaluate(piece, 1.0, dec)?;
        if index + 1 == pieces.len() {
            if distance(corner, ends[1].1) > tolerances.sample {
                return Err(fit_failure("traced edge does not end at its vertex"));
            }
            corners.push(ends[1].1);
        } else {
            corners.push(corner);
        }
    }
    let axis = monotone_axis(&corners);
    let mut best = f64::INFINITY;
    for density in DENSITIES {
        // A world coordinate that is monotone along the edge is a smooth
        // parameter across every traced piece. Without one, each piece keeps
        // its own affine trace parameter and pieces join with C0 continuity.
        let samples = match axis {
            Some(axis) => {
                coordinate_samples(&pieces, &corners, ends, axis, density, &evaluate, dec)?
            }
            None => None,
        };
        let (parameters, points, breaks) = match samples {
            Some(samples) => samples,
            None if density * pieces.len() > MAX_SPANS => break,
            None => {
                let mut points = vec![ends[0].1];
                let mut parameters = vec![0.0];
                let mut breaks = vec![0];
                for (index, piece) in pieces.iter().enumerate() {
                    for j in 1..=density {
                        points.push(if j == density {
                            corners[index + 1]
                        } else {
                            evaluate(piece, j as f64 / density as f64, dec)?
                        });
                        parameters
                            .push((index as f64 + j as f64 / density as f64) / pieces.len() as f64);
                    }
                    breaks.push(points.len() - 1);
                }
                *parameters.last_mut().unwrap() = 1.0;
                (parameters, points, breaks)
            }
        };
        let n = parameters.len() - 1;
        let mut uv_rows = Vec::with_capacity(uses.len());
        let mut sampled = true;
        'faces: for face_use in uses {
            let mut uvs = vec![face_use.ends[0].approx];
            for point in &points[1..n] {
                let uv = face_use.face.surface.invert(
                    *point,
                    face_use.face.region,
                    uvs.last().copied(),
                    dec,
                )?;
                if distance(face_use.face.surface.point(uv), *point) > tolerances.sample {
                    sampled = false;
                    break 'faces;
                }
                uvs.push(uv);
            }
            uvs.push(face_use.ends[1].approx);
            uv_rows.push(uvs);
        }
        if !sampled {
            continue;
        }
        let (knots, controls) = piecewise(&parameters, &points, &breaks, dec)?;
        let mut uv_fits = Vec::with_capacity(uv_rows.len());
        for uvs in &uv_rows {
            uv_fits.push(piecewise(&parameters, uvs, &breaks, dec)?.1);
        }
        let mut worst = 0.0f64;
        for pair in parameters.windows(2) {
            for share in CHECK_SHARES {
                dec.charge_work(1)?;
                let s = lerp([pair[0], pair[1]], share);
                let point = nurbs::evaluate(&knots, &controls, s);
                for (face_use, uv_controls) in uses.iter().zip(&uv_fits) {
                    let uv = nurbs::evaluate(&knots, uv_controls, s);
                    worst = worst.max(distance(face_use.face.surface.point(uv), point));
                }
                for (origin, normal) in planes {
                    worst = worst.max(dot(sub(point, *origin), *normal).abs());
                }
            }
        }
        best = best.min(worst);
        if worst.is_nan() || worst > tolerances.fit {
            continue;
        }
        let domain = [&Rational::zero(), &Rational::one()];
        let curve = rational_curve(
            &knots,
            &controls,
            domain,
            [ends[0].0.clone(), ends[1].0.clone()],
            None,
            dec,
        )?;
        let mut uv_curves = Vec::with_capacity(uses.len());
        for (face_use, uv_controls) in uses.iter().zip(&uv_fits) {
            uv_curves.push(rational_curve(
                &knots,
                uv_controls,
                domain,
                face_use.ends.map(uv_end),
                Some(face_use.face),
                dec,
            )?);
        }
        return Ok((curve, uv_curves));
    }
    Err(fit_failure(format!(
        "a traced edge has no fitted curve and pcurves within the tolerance; closest sampled residual {best:e} mm"
    )))
}

// The world axis along which the traced piece ends are most steadily
// monotone, when every step moves at least a tenth of its length along it.
fn monotone_axis(corners: &[[f64; 3]]) -> Option<usize> {
    let mut best = None;
    let mut best_share = 0.1;
    for axis in 0..3 {
        let sign = (corners[corners.len() - 1][axis] - corners[0][axis]).signum();
        let mut share = f64::INFINITY;
        for pair in corners.windows(2) {
            let length = distance(pair[0], pair[1]);
            let step = (pair[1][axis] - pair[0][axis]) * sign;
            share = share.min(if length > 0.0 { step / length } else { -1.0 });
        }
        if share > best_share {
            best_share = share;
            best = Some(axis);
        }
    }
    best
}

type Samples = (Vec<f64>, Vec<[f64; 3]>, Vec<usize>);
type Evaluate<'a> =
    dyn Fn(&catallax_intersect::TracePiece, f64, &mut Decider) -> Result<[f64; 3], Refusal> + 'a;

// Samples nearly uniform in the monotone coordinate. The parameter is the
// normalised coordinate of each actual sample, so the fit needs no exact
// placement. `None` when the samples are not strictly increasing in it.
fn coordinate_samples(
    pieces: &[&catallax_intersect::TracePiece],
    corners: &[[f64; 3]],
    ends: [(&[Rational; 3], [f64; 3]); 2],
    axis: usize,
    density: usize,
    evaluate: &Evaluate<'_>,
    dec: &mut Decider,
) -> Result<Option<Samples>, Refusal> {
    let lo = corners[0][axis];
    let width = corners[corners.len() - 1][axis] - lo;
    let mut points = vec![ends[0].1];
    let mut parameters = vec![0.0];
    let mut piece = 0;
    for j in 1..density {
        let target = lo + width * (j as f64 / density as f64);
        while piece + 1 < pieces.len() && (corners[piece + 1][axis] - target) * width < 0.0 {
            piece += 1;
        }
        let span = corners[piece + 1][axis] - corners[piece][axis];
        let share = ((target - corners[piece][axis]) / span).clamp(0.0, 1.0);
        let point = evaluate(pieces[piece], share, dec)?;
        let parameter = (point[axis] - lo) / width;
        if parameter.is_nan() || parameter <= *parameters.last().unwrap() {
            return Ok(None);
        }
        points.push(point);
        parameters.push(parameter);
    }
    if *parameters.last().unwrap() >= 1.0 {
        return Ok(None);
    }
    points.push(ends[1].1);
    parameters.push(1.0);
    let breaks = vec![0, points.len() - 1];
    Ok(Some((parameters, points, breaks)))
}

/// Generates numerical input and native UV trims from a retained exact body.
/// Vertices and native parameters are 128-bit midpoints, carriers are the
/// original planes and rational splines, and traced intersection edges are
/// fitted cubics. Fits aim at a quarter of `tolerance_mm`; admission by
/// [`catallax_topology::numerical::build_with_native_trims`] is still required.
pub fn definition_from_exact_model(
    model: &Model,
    body: BodyId,
    name: String,
    tolerance_mm: &Rational,
    dec: &mut Decider,
) -> Result<(NumericalDefinition, Vec<NumericalNativeTrim>), Refusal> {
    if !embeds(model, body) {
        return Err(unsupported(
            "numerical embedding needs one closed shell of supported carriers",
        ));
    }
    let tolerance = tolerance_mm.to_f64_nearest();
    let tolerances = Tolerances {
        sample: tolerance / 64.0,
        fit: tolerance / 4.0,
    };
    let vertex_ids = model.vertices_of_body(body);
    let vertex_index: BTreeMap<VertexId, usize> = vertex_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (*id, i))
        .collect();
    let mut vertices = Vec::with_capacity(vertex_ids.len());
    let mut approx = Vec::with_capacity(vertex_ids.len());
    for id in &vertex_ids {
        dec.charge_work(1)?;
        let p = model.geom().position(model.vertex(*id).point);
        let position = [scalar(&p.x, dec)?, scalar(&p.y, dec)?, scalar(&p.z, dec)?];
        approx.push(position.each_ref().map(Rational::to_f64_nearest));
        vertices.push(NumericalVertex {
            key: format!("{id:?}"),
            position,
        });
    }
    let face_ids = model.faces_of_body(body);
    let mut carriers = Vec::with_capacity(face_ids.len());
    let mut edge_faces: BTreeMap<EdgeId, Vec<usize>> = BTreeMap::new();
    let mut uvs: BTreeMap<(usize, usize), Uv> = BTreeMap::new();
    for (index, id) in face_ids.iter().enumerate() {
        dec.charge_work(1)?;
        let face = model.face(*id);
        let surface = model.geom().surface(face.surface);
        carriers.push(carrier(surface, dec)?);
        for lp in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                unreachable!("checked by embeds")
            };
            for coedge in uses {
                let edge = model.coedge(*coedge).edge;
                edge_faces.entry(edge).or_default().push(index);
                let edge = model.edge(edge);
                for vertex in [edge.start, edge.end].map(|v| vertex_index[&v.expect("finite")]) {
                    if uvs.contains_key(&(index, vertex)) {
                        continue;
                    }
                    let uv = match &carriers[index] {
                        Carrier::Plane { .. } => {
                            let p = vertices[vertex]
                                .position
                                .each_ref()
                                .map(|q| Ball::from_rational(q, BITS));
                            let at = surface
                                .parameters_of(&p, dec)?
                                .ok_or_else(|| bad("gap", "a vertex has no plane parameter"))?;
                            let exact = at.each_ref().map(Ball::mid);
                            Uv {
                                approx: exact.each_ref().map(Rational::to_f64_nearest),
                                exact,
                            }
                        }
                        Carrier::Spline(spline) => {
                            let target = approx[vertex];
                            let at = spline.surface.invert(target, spline.region, None, dec)?;
                            if distance(spline.surface.point(at), target) > tolerances.sample {
                                return Err(fit_failure("a vertex has no spline face parameter"));
                            }
                            Uv {
                                exact: [spline.exact(0, at[0])?, spline.exact(1, at[1])?],
                                approx: at,
                            }
                        }
                    };
                    uvs.insert((index, vertex), uv);
                }
            }
        }
    }
    let edge_ids = model.edges_of_body(body);
    let edge_index: BTreeMap<EdgeId, usize> = edge_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (*id, i))
        .collect();
    let mut edges = Vec::with_capacity(edge_ids.len());
    let mut uv_curves: BTreeMap<(usize, usize), BSplineCurve> = BTreeMap::new();
    for (index, id) in edge_ids.iter().enumerate() {
        dec.charge_work(1)?;
        let edge = model.edge(*id);
        let ends = [edge.start, edge.end].map(|v| vertex_index[&v.expect("finite")]);
        let faces = edge_faces.get(id).cloned().unwrap_or_default();
        let curve = model.geom().curve(edge.curve);
        match curve {
            Curve::Line(_) | Curve::BSpline(_) => {
                let mut parameters = Vec::with_capacity(2);
                for vertex in ends {
                    let p = model
                        .geom()
                        .position(model.vertex(vertex_ids[vertex]).point);
                    parameters.push(match curve {
                        Curve::Line(line) => scalar(
                            &(p.minus(&line.origin).dot(&line.dir) / line.dir.norm_squared()),
                            dec,
                        )?,
                        Curve::BSpline(spline) => {
                            let at = spline.unique_parameter_of(&p, dec)?.ok_or_else(|| {
                                bad("gap", "a vertex has no finite spline parameter")
                            })?;
                            let [lo, hi] = spline.domain();
                            scalar(at.value(), dec)?.max(lo.clone()).min(hi.clone())
                        }
                        _ => unreachable!(),
                    });
                }
                let parameters: [Rational; 2] = parameters.try_into().unwrap();
                let native = match curve {
                    Curve::Line(line) => Native::Line {
                        origin: f64_point(&line.origin),
                        dir: f64_point(&Point3::new(
                            line.dir.x.clone(),
                            line.dir.y.clone(),
                            line.dir.z.clone(),
                        )),
                    },
                    Curve::BSpline(spline) => Native::Spline(Curve64::new(spline, dec)?),
                    _ => unreachable!(),
                };
                let ascending = parameters[0] < parameters[1];
                let interval = if ascending {
                    [&parameters[0], &parameters[1]]
                } else {
                    [&parameters[1], &parameters[0]]
                };
                for face in &faces {
                    let Carrier::Spline(spline) = &carriers[*face] else {
                        continue;
                    };
                    let [a, b] = ends.map(|vertex| &uvs[&(*face, vertex)]);
                    let ends = if ascending { [a, b] } else { [b, a] };
                    uv_curves.insert(
                        (*face, index),
                        native_uv(spline, &native, interval, ends, &tolerances, dec)?,
                    );
                }
                edges.push(NumericalEdge {
                    key: format!("{id:?}"),
                    curve: curve.clone(),
                    vertices: ends,
                    parameters,
                });
            }
            Curve::PatchSection(_) | Curve::PatchPair(_) => {
                let mut uses = Vec::new();
                let mut spline_faces = Vec::new();
                let mut planes = Vec::new();
                for face in &faces {
                    match &carriers[*face] {
                        Carrier::Spline(spline) => {
                            spline_faces.push(*face);
                            uses.push(FittedUse {
                                face: spline,
                                ends: ends.map(|vertex| &uvs[&(*face, vertex)]),
                            });
                        }
                        Carrier::Plane { origin, normal } => planes.push((*origin, *normal)),
                    }
                }
                let (fit, fitted_uvs) = fitted_edge(
                    model,
                    *id,
                    &uses,
                    &planes,
                    ends.map(|vertex| (&vertices[vertex].position, approx[vertex])),
                    &tolerances,
                    dec,
                )?;
                for (face, uv) in spline_faces.into_iter().zip(fitted_uvs) {
                    uv_curves.insert((face, index), uv);
                }
                edges.push(NumericalEdge {
                    key: format!("{id:?}"),
                    curve: Curve::BSpline(fit),
                    vertices: ends,
                    parameters: [Rational::zero(), Rational::one()],
                });
            }
            _ => unreachable!("checked by embeds"),
        }
    }
    let mut faces = Vec::with_capacity(face_ids.len());
    let mut native_trims = Vec::new();
    for (index, id) in face_ids.iter().enumerate() {
        dec.charge_work(1)?;
        let face = model.face(*id);
        let key = format!("{id:?}");
        let mut loops = Vec::with_capacity(face.loops.len());
        for (boundary, lp) in face.loops.iter().enumerate() {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                unreachable!("checked by embeds")
            };
            let mut trims = Vec::with_capacity(coedges.len());
            for (use_index, coedge) in coedges.iter().enumerate() {
                let coedge = model.coedge(*coedge);
                let edge = edge_index[&coedge.edge];
                let mut walked = edges[edge].vertices;
                if coedge.sense == Sense::Reversed {
                    walked.reverse();
                }
                trims.push(NumericalTrim {
                    edge,
                    sense: coedge.sense,
                    uv: walked.map(|vertex| uvs[&(index, vertex)].exact.clone()),
                });
                if let Some(curve) = uv_curves.get(&(index, edge)) {
                    native_trims.push(NumericalNativeTrim {
                        face: key.clone(),
                        boundary,
                        use_index,
                        edge: edges[edge].key.clone(),
                        curve: curve.clone(),
                    });
                }
            }
            loops.push(trims);
        }
        faces.push(NumericalFace {
            key,
            surface: model.geom().surface(face.surface).clone(),
            sense: face.sense,
            loops,
        });
    }
    Ok((
        NumericalDefinition {
            name,
            vertices,
            edges,
            faces,
        },
        native_trims,
    ))
}

/// Embeds, admits and writes a retained exact body as an approximate STEP
/// artifact under `options`, stated before any work. Admission and the
/// writer's decimal payload are each checked at that tolerance.
pub fn build_from_exact_model(
    model: &Model,
    body: BodyId,
    name: String,
    options: NumericalOptions,
    step: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<(NumericalBody, NumericalStep), Refusal> {
    guarded(|| {
        let (definition, native_trims) =
            definition_from_exact_model(model, body, name, &options.tolerance_mm, dec)?;
        let body = catallax_topology::numerical::build_with_native_trims(
            definition,
            native_trims,
            options,
            dec,
        )?;
        let artifact = build(&body, step, dec)?;
        Ok((body, artifact))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::BSplineSurface;

    // A bilinear patch over [0, 1/3] x [1/10, 2/3]; no side is a binary64 value.
    fn face() -> SplineFace {
        let third = Rational::from_ratio(1, 3).unwrap();
        let tenth = Rational::from_ratio(1, 10).unwrap();
        let two_thirds = Rational::from_ratio(2, 3).unwrap();
        let mut dec = Decider::default();
        let source = BSplineSurface::new(
            [1, 1],
            [
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    third.clone(),
                    third.clone(),
                ],
                vec![
                    tenth.clone(),
                    tenth.clone(),
                    two_thirds.clone(),
                    two_thirds.clone(),
                ],
            ],
            [2, 2],
            [[0i64, 0], [0, 1], [1, 0], [1, 1]]
                .map(|[x, y]| Point3::new(Real::from(x), Real::from(y), Real::zero()))
                .to_vec(),
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap();
        SplineFace::new(
            Surface64::new(&source, &mut dec).unwrap(),
            [[Rational::zero(), third], [tenth, two_thirds]],
        )
    }

    #[test]
    fn trims_on_a_non_binary64_side_keep_the_exact_side() {
        let face = face();
        let third = Rational::from_ratio(1, 3).unwrap();
        let tenth = Rational::from_ratio(1, 10).unwrap();
        let two_thirds = Rational::from_ratio(2, 3).unwrap();
        // The binary64 copies of 1/3, 1/10 and 2/3 differ from the sides.
        for side in [&third, &tenth, &two_thirds] {
            assert_ne!(&exact(side.to_f64_nearest()).unwrap(), side);
        }
        assert_eq!(face.exact(0, third.to_f64_nearest()).unwrap(), third);
        assert_eq!(face.exact(1, tenth.to_f64_nearest()).unwrap(), tenth);
        assert_eq!(
            face.exact(1, two_thirds.to_f64_nearest() + 1e-16).unwrap(),
            two_thirds
        );
        // Outside the box clamps to the side; inside keeps the proposed value.
        assert_eq!(face.exact(0, 0.4).unwrap(), third);
        assert_eq!(face.exact(1, 0.0).unwrap(), tenth);
        assert_eq!(
            face.exact(0, 0.25).unwrap(),
            Rational::from_ratio(1, 4).unwrap()
        );
        // A fitted UV curve along the side u = 1/3 keeps every control on it.
        let curve = rational_curve(
            &[0.1, 0.1, 0.1, 0.1, 0.6, 0.6, 0.6, 0.6],
            &[[third.to_f64_nearest(), 0.1]; 4],
            [&tenth, &Rational::from_ratio(3, 5).unwrap()],
            [
                [third.clone(), tenth.clone(), Rational::zero()],
                [
                    third.clone(),
                    Rational::from_ratio(3, 5).unwrap(),
                    Rational::zero(),
                ],
            ],
            Some(&face),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(
            curve
                .controls()
                .all(|p| p.x.as_rational() == Some(third.clone()))
        );
    }
}
