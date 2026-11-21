//! Numerical embedding of retained exact topology and carriers.
//!
//! Planes and native spline carriers are kept. Each traced patch-section or
//! patch-pair edge becomes a fitted cubic, and every spline-face edge use gets
//! a fitted native UV curve at that edge's parameter. Binary64 arithmetic only
//! proposes fits. Numerical admission must recheck every edge, trim and loop
//! at 128 bits under the declared tolerance. This returns no admitted body,
//! exact answer, complete body bound, or STEP authority.
//!
//! Sines, cosines, arctangents and hypotenuses in binary64 go through the
//! `libm` crate (the musl library's routines in Rust), not the methods of
//! `f64`, which call the platform's library and differ in the last bit from
//! one machine to the next. The work that follows from these values is
//! counted, and the predictor's done-when lines ask for the same work on two
//! machines (Checkpoint 2 decisions, section 4.5).
//!
//! This follows OCCT's face/face section, which approximates each section as
//! one 3D B-spline and one 2D B-spline on each face (`IntTools_FaceFace` with
//! curve and both pcurve approximations enabled) before the Boolean stores
//! them on its section edges.

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, Curve, Plane, Point3, SplineSupport, Surface};
use catallax_topology::numerical::{
    NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalTrim,
    NumericalVertex,
};
use catallax_topology::{BodyId, EdgeId, Loop, Model, Sense, VertexId};
use std::collections::{BTreeMap, BTreeSet};

mod charts;
mod guides;
pub(crate) mod nurbs;
mod planar;
mod rings;
use nurbs::{Curve64, Surface64, concatenate, distance, dot, interpolate, lerp, sub};

const BITS: u32 = 128;
const SAMPLE_BITS: u32 = 64;
// Samples per traced piece or per native edge, doubled until the binary64
// fit meets a quarter of the tolerance. Admission samples every span at its
// eighths, so the fit is checked at the same places.
const DENSITIES: [usize; 7] = [4, 8, 16, 32, 64, 128, 256];
// Admission checks every span of every loop pairwise, so fits stay small.
const MAX_SPANS: usize = 256;
const CHECK_SHARES: [f64; 7] = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875];

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("features.numerical-embedding-unsupported", detail)
}

fn fit_failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("features.numerical-embedding-fit", detail)
}

fn f64_point(p: &Point3, dec: &Decider) -> Result<[f64; 3], Refusal> {
    Ok([
        scalar(&p.x, dec)?.to_f64_nearest(),
        scalar(&p.y, dec)?.to_f64_nearest(),
        scalar(&p.z, dec)?.to_f64_nearest(),
    ])
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| fit_failure("a fitted value is not finite"))
}

fn scalar(value: &Real, dec: &Decider) -> Result<Rational, Refusal> {
    Ok(dec.enclose_at(value, BITS)?.mid())
}

// Inverse results select a chart; they are not incidence certificates. At a
// seam, a small displaced proposal chooses a neighbouring representative.
// Callers recheck its physical point before retaining a trim.
fn inverse_seed(
    surface: &Surface,
    point: &mut [Ball; 3],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<([Rational; 2], bool), Refusal> {
    if let Some(uv) = surface.parameters_of(point, dec)? {
        return Ok((uv.each_ref().map(Ball::mid), false));
    }
    let displacement = Ball::from_rational(&tolerance.mul_pow2(-6), BITS);
    for axis in 0..3 {
        let shifted = &point[axis] + &displacement;
        let original = std::mem::replace(&mut point[axis], shifted);
        let proposal = surface.parameters_of(point, dec)?;
        point[axis] = original;
        if let Some(uv) = proposal {
            return Ok((uv.each_ref().map(Ball::mid), true));
        }
    }
    Err(fit_failure(
        "the point has no regular surface-chart inverse seed",
    ))
}

pub(crate) fn surface_parameters(
    surface: &Surface,
    point: &Point3,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    inverse_seed(
        surface,
        &mut [
            dec.enclose_budgeted_at(&point.x, BITS)?,
            dec.enclose_budgeted_at(&point.y, BITS)?,
            dec.enclose_budgeted_at(&point.z, BITS)?,
        ],
        tolerance,
        dec,
    )
    .map(|(uv, _)| uv)
}

/// Structural precheck for one shell with finite coedges or circular rings.
/// Native carriers and explicit line/circle/spline edges are retained; other
/// curves require a retained trace. This does not qualify charts or admission.
pub fn embeds(model: &Model, body: BodyId) -> bool {
    let counts = model.counts(body);
    if counts.lumps != 1 || counts.shells != 1 {
        return false;
    }
    for face in model.faces_of_body(body) {
        let face = model.face(face);
        if face.loops.is_empty() {
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
        (edge.start.is_some() && edge.end.is_some()) || rings::is_ring(model, id)
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

enum Carrier<'a> {
    Plane { origin: [f64; 3], normal: [f64; 3] },
    Spline(Box<SplineFace>),
    Native(&'a Surface),
}

impl Carrier<'_> {
    fn parametric(&self) -> Option<ParametricFace<'_>> {
        match self {
            Self::Plane { .. } => None,
            Self::Spline(face) => Some(ParametricFace::Spline(face)),
            Self::Native(surface) => Some(ParametricFace::Native(surface)),
        }
    }
}

fn carrier<'a>(surface: &'a Surface, dec: &mut Decider) -> Result<Carrier<'a>, Refusal> {
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
            let n = [
                scalar(&normal.x, dec)?.to_f64_nearest(),
                scalar(&normal.y, dec)?.to_f64_nearest(),
                scalar(&normal.z, dec)?.to_f64_nearest(),
            ];
            let length = dot(n, n).sqrt();
            Carrier::Plane {
                origin: f64_point(origin, dec)?,
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
        Surface::BSpline(source) => {
            let source = source.carrier();
            Carrier::Spline(Box::new(SplineFace::new(
                Surface64::new(source, dec)?,
                source.domain().map(|range| range.map(Clone::clone)),
            )))
        }
        _ => Carrier::Native(surface),
    })
}

#[derive(Clone, Copy)]
enum ParametricFace<'a> {
    Spline(&'a SplineFace),
    Native(&'a Surface),
}

impl<'a> ParametricFace<'a> {
    fn spline(self) -> Option<&'a SplineFace> {
        match self {
            Self::Spline(face) => Some(face),
            Self::Native(_) => None,
        }
    }

    fn point(self, uv: [f64; 2], dec: &mut Decider) -> Result<[f64; 3], Refusal> {
        let point = match self {
            Self::Spline(face) => face.surface.point(uv),
            Self::Native(surface) => {
                let at = [exact(uv[0])?, exact(uv[1])?]
                    .each_ref()
                    .map(|q| Ball::from_rational(q, BITS));
                surface
                    .evaluate(&at[0], &at[1], dec)?
                    .each_ref()
                    .map(|v| v.mid().to_f64_nearest())
            }
        };
        if !point.iter().all(|v| v.is_finite()) {
            return Err(fit_failure(
                "a native surface point is outside binary64 range",
            ));
        }
        Ok(point)
    }

    fn invert(
        self,
        point: [f64; 3],
        previous: Option<[f64; 2]>,
        tolerances: &Tolerances,
        dec: &mut Decider,
    ) -> Result<[f64; 2], Refusal> {
        match self {
            Self::Spline(face) => face.surface.invert(point, face.region, previous, dec),
            Self::Native(surface) => {
                let mut p = [exact(point[0])?, exact(point[1])?, exact(point[2])?]
                    .each_ref()
                    .map(|q| Ball::from_rational(q, BITS));
                let (seed, displaced) = inverse_seed(surface, &mut p, &tolerances.projection, dec)?;
                let mut uv = seed.each_ref().map(Rational::to_f64_nearest);
                if let Some(previous) = previous {
                    for axis in 0..2 {
                        if surface.domain()[axis] == catallax_shapes::ParamKind::Periodic {
                            let turns =
                                ((previous[axis] - uv[axis]) / std::f64::consts::TAU).round();
                            uv[axis] += turns * std::f64::consts::TAU;
                        }
                    }
                }
                if !uv.iter().all(|v| v.is_finite()) {
                    return Err(fit_failure(
                        "a native inverse parameter is outside binary64 range",
                    ));
                }
                if !displaced {
                    return Ok(uv);
                }
                // The seam displacement chooses a branch, not a different
                // physical endpoint. Refine its UV against the original point.
                for _ in 0..8 {
                    dec.charge_work(1)?;
                    let at = [exact(uv[0])?, exact(uv[1])?]
                        .each_ref()
                        .map(|q| Ball::from_rational(q, BITS));
                    let position = surface.evaluate(&at[0], &at[1], dec)?;
                    let error: [f64; 3] =
                        std::array::from_fn(|i| (&position[i] - &p[i]).mid().to_f64_nearest());
                    if dot(error, error).sqrt() <= tolerances.sample / 8.0 {
                        return Ok(uv);
                    }
                    let [du, dv] = surface
                        .derivatives(&at[0], &at[1], dec)?
                        .map(|d| d.each_ref().map(|q| q.mid().to_f64_nearest()));
                    let a = dot(du, du);
                    let b = dot(du, dv);
                    let c = dot(dv, dv);
                    let d = dot(du, error);
                    let e = dot(dv, error);
                    let determinant = a * c - b * b;
                    if determinant <= 0.0 || !determinant.is_finite() {
                        return Err(fit_failure(
                            "the native seam inverse tangent system is singular",
                        ));
                    }
                    uv = [
                        uv[0] - (c * d - b * e) / determinant,
                        uv[1] - (a * e - b * d) / determinant,
                    ];
                }
                Err(fit_failure(
                    "the native seam inverse did not refine its original point",
                ))
            }
        }
    }
}

#[derive(Clone)]
struct Uv {
    exact: [Rational; 2],
    approx: [f64; 2],
}

struct Circle64 {
    center: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    radius: f64,
}

impl Circle64 {
    fn new(circle: &catallax_shapes::Circle, dec: &mut Decider) -> Result<Self, Refusal> {
        let plane = Surface::Plane(Plane {
            origin: circle.center.clone(),
            normal: circle.normal.clone(),
            ref_dir: circle.ref_dir.clone(),
        });
        let zero = Ball::from_rational(&Rational::zero(), BITS);
        let [x, y] = plane.derivatives(&zero, &zero, dec)?;
        let result = Self {
            center: f64_point(&circle.center, dec)?,
            x: x.each_ref().map(|v| v.mid().to_f64_nearest()),
            y: y.each_ref().map(|v| v.mid().to_f64_nearest()),
            radius: scalar(&circle.radius_sq.sqrt(), dec)?.to_f64_nearest(),
        };
        if result.radius <= 0.0
            || !result.radius.is_finite()
            || !result
                .center
                .iter()
                .chain(&result.x)
                .chain(&result.y)
                .all(|v| v.is_finite())
        {
            return Err(fit_failure(
                "a retained circle has no finite binary64 proposal frame",
            ));
        }
        Ok(result)
    }

    fn point(&self, t: f64) -> [f64; 3] {
        let (sin, cos) = libm::sincos(t);
        std::array::from_fn(|i| self.center[i] + self.radius * (self.x[i] * cos + self.y[i] * sin))
    }

    fn phase(&self, point: [f64; 3]) -> f64 {
        let delta = sub(point, self.center);
        libm::atan2(dot(delta, self.y), dot(delta, self.x))
    }
}

enum Native {
    Line { origin: [f64; 3], dir: [f64; 3] },
    Circle(Circle64),
    Spline(Curve64),
}

impl Native {
    fn point(&self, t: f64) -> [f64; 3] {
        match self {
            Self::Line { origin, dir } => std::array::from_fn(|c| origin[c] + t * dir[c]),
            Self::Circle(circle) => circle.point(t),
            Self::Spline(curve) => curve.point(t),
        }
    }
}

struct Tolerances {
    sample: f64,
    fit: f64,
    projection: Rational,
}

fn rational_curve_definition<const D: usize>(
    knots: &[f64],
    controls: &[[f64; D]],
    domain: [&Rational; 2],
    ends: [[Rational; 3]; 2],
    face: Option<&SplineFace>,
    dec: &mut Decider,
) -> Result<(usize, Vec<Rational>, Vec<Point3>), Refusal> {
    dec.charge_work((knots.len() + D * controls.len()) as u64)?;
    let degree = knots.len() - controls.len() - 1;
    let mut exact_knots = Vec::with_capacity(knots.len());
    for (i, knot) in knots.iter().enumerate() {
        exact_knots.push(if i <= degree {
            domain[0].clone()
        } else if i + degree + 1 >= knots.len() {
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
    Ok((degree, exact_knots, points))
}

fn rational_curve<const D: usize>(
    knots: &[f64],
    controls: &[[f64; D]],
    domain: [&Rational; 2],
    ends: [[Rational; 3]; 2],
    face: Option<&SplineFace>,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let (degree, knots, points) =
        rational_curve_definition(knots, controls, domain, ends, face, dec)?;
    let weights = vec![Real::one(); points.len()];
    BSplineCurve::new(degree, knots, points, weights, dec)
}

fn uv_end(uv: &Uv) -> [Rational; 3] {
    [uv.exact[0].clone(), uv.exact[1].clone(), Rational::zero()]
}

// A native UV curve for a retained line, circle or spline at that edge's
// own parameter over `interval`, ascending.
fn native_uv(
    face: ParametricFace<'_>,
    edge: &Native,
    interval: [&Rational; 2],
    ends: [&Uv; 2],
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let range = interval.map(Rational::to_f64_nearest);
    let mut failure = None;
    let mut cuts = None;
    for density in DENSITIES {
        for partitioned in [false, true] {
            let sections = if partitioned {
                if cuts.is_none() {
                    cuts = Some(charts::parameter_cuts(
                        face, edge, range, ends, tolerances, dec,
                    )?);
                }
                cuts.as_ref().unwrap().as_slice()
            } else {
                &range
            };
            let pieces = sections.len() - 1;
            if partitioned && DENSITIES[0] * pieces > density {
                continue;
            }
            let mut parameters = Vec::with_capacity(density + 1);
            let mut breaks = Vec::new();
            parameters.push(range[0]);
            if partitioned {
                breaks.reserve(pieces + 1);
                breaks.push(0);
            }
            for (index, section) in sections.windows(2).enumerate() {
                let count = density / pieces + usize::from(index < density % pieces);
                for j in 1..=count {
                    parameters.push(if j == count {
                        section[1]
                    } else {
                        lerp([section[0], section[1]], j as f64 / count as f64)
                    });
                }
                if partitioned {
                    breaks.push(parameters.len() - 1);
                }
            }
            let mut uvs = Vec::with_capacity(parameters.len());
            uvs.push(ends[0].approx);
            let mut sampled = true;
            for t in &parameters[1..density] {
                let point = edge.point(*t);
                let uv = face.invert(point, uvs.last().copied(), tolerances, dec)?;
                let error = distance(face.point(uv, dec)?, point);
                if error > tolerances.sample {
                    failure = Some(("inverse", density, *t, error, tolerances.sample));
                    sampled = false;
                    break;
                }
                uvs.push(uv);
            }
            if !sampled {
                continue;
            }
            uvs.push(ends[1].approx);
            let (knots, controls) = if partitioned {
                piecewise_quartic(&parameters, &uvs, &breaks, dec)?
            } else {
                interpolate(&parameters, &uvs, dec)?
            };
            let mut fits = true;
            'spans: for pair in parameters.windows(2) {
                for share in CHECK_SHARES {
                    dec.charge_work(1)?;
                    let t = lerp([pair[0], pair[1]], share);
                    let uv = nurbs::evaluate(&knots, &controls, t);
                    let error = distance(face.point(uv, dec)?, edge.point(t));
                    if error > tolerances.fit {
                        failure = Some(("fit", density, t, error, tolerances.fit));
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
                    face.spline(),
                    dec,
                );
            }
        }
    }
    let (stage, density, parameter, error, limit) = failure.expect("all densities failed");
    let source_ends = [edge.point(range[0]), edge.point(range[1])];
    let surface_ends = [
        face.point(ends[0].approx, dec)?,
        face.point(ends[1].approx, dec)?,
    ];
    let endpoint_errors =
        std::array::from_fn::<_, 2, _>(|i| distance(source_ends[i], surface_ends[i]));
    let chart_knots = face.spline().map(|face| {
        [0, 1].map(|axis| {
            let knots = face.surface.internal_knots(axis);
            (
                knots.len(),
                knots
                    .iter()
                    .enumerate()
                    // Repeated knots are stored as identical values.
                    .filter(|(i, knot)| *i == 0 || knot.to_bits() != knots[*i - 1].to_bits())
                    .count(),
            )
        })
    });
    let chart_jets = face
        .spline()
        .map(|face| ends.map(|uv| face.surface.jet(uv.approx)));
    Err(fit_failure(format!(
        "native UV {stage} residual {error:e} mm exceeds {limit:e} mm at parameter {parameter:e} with density {density}; endpoint residuals {endpoint_errors:?}, curve ends {source_ends:?}, surface ends {surface_ends:?}, chart knot counts (total, unique) {chart_knots:?}, endpoint jets {chart_jets:?}"
    )))
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

fn piecewise_quartic<const D: usize>(
    parameters: &[f64],
    values: &[[f64; D]],
    breaks: &[usize],
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<[f64; D]>), Refusal> {
    let mut pieces = Vec::new();
    for pair in breaks.windows(2) {
        let mut from = pair[0];
        while from < pair[1] {
            let to = if pair[1] - from < 8 {
                pair[1]
            } else {
                from + 4
            };
            let indices = std::array::from_fn::<_, 5, _>(|i| from + (i * (to - from) + 2) / 4);
            pieces.push(nurbs::interpolate_quartic(
                &indices.map(|i| parameters[i]),
                &indices.map(|i| values[i]),
                dec,
            )?);
            from = to;
        }
    }
    Ok(concatenate(pieces))
}

struct FittedUse<'a> {
    face: ParametricFace<'a>,
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
                let uv = face_use
                    .face
                    .invert(*point, uvs.last().copied(), tolerances, dec)?;
                if distance(face_use.face.point(uv, dec)?, *point) > tolerances.sample {
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
                    worst = worst.max(distance(face_use.face.point(uv, dec)?, point));
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
                face_use.face.spline(),
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
/// Vertices are 128-bit midpoints. Native carriers and curves are retained;
/// angular ranges and fitted UV curves are numerical proposals, and other
/// traced intersection edges become fitted cubics. Fits aim at a quarter of
/// `tolerance_mm`; [`catallax_topology::numerical::build_with_native_trims`]
/// must still admit the resulting body.
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
        projection: tolerance_mm.clone(),
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
    let ring_vertices = rings::vertices(model, body, &mut vertices, &mut approx, dec)?;
    let edge_ends: BTreeMap<EdgeId, [usize; 2]> = model
        .edges_of_body(body)
        .into_iter()
        .map(|id| {
            let edge = model.edge(id);
            let ends = match [edge.start, edge.end] {
                [Some(a), Some(b)] => [vertex_index[&a], vertex_index[&b]],
                _ => [ring_vertices[&id]; 2],
            };
            (id, ends)
        })
        .collect();
    let face_ids = model.faces_of_body(body);
    let mut carriers = Vec::with_capacity(face_ids.len());
    let mut edge_faces: BTreeMap<EdgeId, Vec<usize>> = BTreeMap::new();
    let mut uvs: BTreeMap<(usize, usize), Uv> = BTreeMap::new();
    for (index, id) in face_ids.iter().enumerate() {
        dec.charge_work(1)?;
        let face = model.face(*id);
        let surface = model.geom().surface(face.surface);
        carriers.push(carrier(surface, dec)?);
        let mut previous_native_uv = None;
        for lp in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                unreachable!("checked by embeds")
            };
            for coedge in uses {
                let edge = model.coedge(*coedge).edge;
                edge_faces.entry(edge).or_default().push(index);
                for vertex in edge_ends[&edge] {
                    if uvs.contains_key(&(index, vertex)) {
                        continue;
                    }
                    let uv = match &carriers[index] {
                        Carrier::Plane { .. } => {
                            let p = vertices[vertex]
                                .position
                                .each_ref()
                                .map(|q| Ball::from_rational(q, BITS));
                            let at = surface.parameters_of(&p, dec)?.ok_or_else(|| {
                                Refusal::new(
                                    "features.numerical-embedding-gap",
                                    "a vertex has no plane parameter",
                                )
                            })?;
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
                        Carrier::Native(surface) => {
                            let face = ParametricFace::Native(surface);
                            let target = approx[vertex];
                            let at = face.invert(target, previous_native_uv, &tolerances, dec)?;
                            let error = distance(face.point(at, dec)?, target);
                            if error > tolerances.sample {
                                return Err(fit_failure(format!(
                                    "vertex {} on {id:?} has native UV {at:?} with residual {error:e} mm above {:e} mm",
                                    vertices[vertex].key, tolerances.sample
                                )));
                            }
                            previous_native_uv = Some(at);
                            Uv {
                                exact: [exact(at[0])?, exact(at[1])?],
                                approx: at,
                            }
                        }
                    };
                    uvs.insert((index, vertex), uv);
                }
            }
        }
    }
    let mut endpoint_uvs = BTreeMap::new();
    charts::lift(
        model,
        body,
        &carriers,
        &edge_ends,
        &approx,
        &mut uvs,
        &mut endpoint_uvs,
        &tolerances,
        dec,
    )?;
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
        let ends = edge_ends[id];
        let faces = edge_faces.get(id).cloned().unwrap_or_default();
        let curve = model.geom().curve(edge.curve);
        match curve {
            Curve::Line(_) | Curve::Circle(_) | Curve::BSpline(_) => {
                let native = match curve {
                    Curve::Line(line) => Native::Line {
                        origin: f64_point(&line.origin, dec)?,
                        dir: [
                            scalar(&line.dir.x, dec)?.to_f64_nearest(),
                            scalar(&line.dir.y, dec)?.to_f64_nearest(),
                            scalar(&line.dir.z, dec)?.to_f64_nearest(),
                        ],
                    },
                    Curve::Circle(circle) => Native::Circle(Circle64::new(circle, dec)?),
                    Curve::BSpline(spline) => Native::Spline(Curve64::new(spline, dec)?),
                    _ => unreachable!("retained primitive curve"),
                };
                let parameter = |vertex: usize, dec: &mut Decider| -> Result<Rational, Refusal> {
                    Ok(match curve {
                        Curve::Line(line) => {
                            let p = model
                                .geom()
                                .position(model.vertex(vertex_ids[vertex]).point);
                            scalar(
                                &(p.minus(&line.origin).dot(&line.dir) / line.dir.norm_squared()),
                                dec,
                            )?
                        }
                        Curve::Circle(_) => {
                            let Native::Circle(circle) = &native else {
                                unreachable!("retained circle proposal")
                            };
                            exact(circle.phase(approx[vertex]))?
                        }
                        Curve::BSpline(spline) => {
                            let p = model
                                .geom()
                                .position(model.vertex(vertex_ids[vertex]).point);
                            let at = spline.unique_parameter_of(&p, dec)?.ok_or_else(|| {
                                Refusal::new(
                                    "features.numerical-embedding-gap",
                                    "a vertex has no finite spline parameter",
                                )
                            })?;
                            let [lo, hi] = spline.domain();
                            scalar(at.value(), dec)?.max(lo.clone()).min(hi.clone())
                        }
                        _ => unreachable!("retained primitive curve"),
                    })
                };
                let mut parameters = [parameter(ends[0], dec)?, parameter(ends[1], dec)?];
                if matches!(curve, Curve::Circle(_)) {
                    let period = Ball::pi(BITS).mid().mul_pow2(1);
                    if edge.sense == Sense::Same && parameters[1] <= parameters[0] {
                        parameters[1] = &parameters[1] + &period;
                    } else if edge.sense == Sense::Reversed && parameters[1] >= parameters[0] {
                        parameters[1] = &parameters[1] - &period;
                    }
                }
                let ascending = parameters[0] < parameters[1];
                let interval = if ascending {
                    [&parameters[0], &parameters[1]]
                } else {
                    [&parameters[1], &parameters[0]]
                };
                for face in &faces {
                    let Some(parametric) = carriers[*face].parametric() else {
                        continue;
                    };
                    if ends[0] == ends[1]
                        && matches!(curve, Curve::Circle(_))
                        && matches!(parametric, ParametricFace::Native(_))
                    {
                        continue;
                    }
                    let [a, b] = ends.map(|vertex| {
                        endpoint_uvs
                            .get(&(*face, *id, vertex))
                            .unwrap_or(&uvs[&(*face, vertex)])
                    });
                    let ends = if ascending { [a, b] } else { [b, a] };
                    uv_curves.insert(
                        (*face, index),
                        native_uv(parametric, &native, interval, ends, &tolerances, dec).map_err(
                            |mut refusal| {
                                refusal.detail = format!(
                                    "native trim of {id:?} on {:?}, parameter range {:?}, UV ends {:?}: {}",
                                    face_ids[*face],
                                    interval.map(Rational::to_f64_nearest),
                                    ends.map(|uv| uv.approx),
                                    refusal.detail
                                );
                                refusal
                            },
                        )?,
                    );
                }
                edges.push(NumericalEdge {
                    key: format!("{id:?}"),
                    curve: curve.clone(),
                    vertices: ends,
                    parameters,
                });
            }
            _ => {
                let mut uses = Vec::new();
                let mut parametric_faces = Vec::new();
                let mut planes = Vec::new();
                for face in &faces {
                    match &carriers[*face] {
                        Carrier::Spline(_) | Carrier::Native(_) => {
                            parametric_faces.push(*face);
                            uses.push(FittedUse {
                                face: carriers[*face].parametric().expect("parametric carrier"),
                                ends: ends.map(|vertex| {
                                    endpoint_uvs
                                        .get(&(*face, *id, vertex))
                                        .unwrap_or(&uvs[&(*face, vertex)])
                                }),
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
                for (face, uv) in parametric_faces.into_iter().zip(fitted_uvs) {
                    uv_curves.insert((face, index), uv);
                }
                edges.push(NumericalEdge {
                    key: format!("{id:?}"),
                    curve: Curve::BSpline(fit),
                    vertices: ends,
                    parameters: [Rational::zero(), Rational::one()],
                });
            }
        }
    }
    let mut faces = Vec::with_capacity(face_ids.len());
    let mut native_trims = Vec::new();
    for (index, id) in face_ids.iter().enumerate() {
        dec.charge_work(1)?;
        let face = model.face(*id);
        let key = format!("{id:?}");
        let mut loops = Vec::with_capacity(face.loops.len());
        for lp in &face.loops {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                unreachable!("checked by embeds")
            };
            let mut trims = Vec::with_capacity(coedges.len());
            for coedge in coedges {
                let coedge = model.coedge(*coedge);
                let edge = edge_index[&coedge.edge];
                let mut walked = edges[edge].vertices;
                if coedge.sense == Sense::Reversed {
                    walked.reverse();
                }
                trims.push(NumericalTrim {
                    tag: catallax_topology::numerical::CopyTag::Independent,
                    edge,
                    sense: coedge.sense,
                    uv: walked.map(|vertex| {
                        endpoint_uvs
                            .get(&(index, coedge.edge, vertex))
                            .unwrap_or(&uvs[&(index, vertex)])
                            .exact
                            .clone()
                    }),
                });
            }
            loops.push(trims);
        }
        let mut numerical_face = NumericalFace {
            key,
            surface: model.geom().surface(face.surface).clone(),
            sense: face.sense,
            loops,
        };
        planar::outer_first(&mut numerical_face, &edges, dec)?;
        rings::close_band(&mut numerical_face, &vertices, &mut edges, dec)?;
        for (boundary, uses) in numerical_face.loops.iter().enumerate() {
            for (use_index, trim) in uses.iter().enumerate() {
                if let Some(curve) = uv_curves.get(&(index, trim.edge)) {
                    native_trims.push(NumericalNativeTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        face: numerical_face.key.clone(),
                        boundary,
                        use_index,
                        edge: edges[trim.edge].key.clone(),
                        curve: curve.clone(),
                    });
                }
            }
        }
        faces.push(numerical_face);
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

pub(crate) struct FilletGuide {
    pub curve: BSplineCurve,
    pub ends: [Point3; 2],
    pub check: catallax_numbers::NumericalCheck,
}

/// Fits a section-plane guide from the retained curve trace, not its carrier
/// classification. This includes Boolean intersections whose line/circle
/// identity is no longer explicit. The checks below bound sampled deviation,
/// not the deviation between samples or the completeness of a branch.
pub(crate) fn fillet_guide(
    model: &Model,
    edge_id: EdgeId,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<FilletGuide, Refusal> {
    let edge = model.edge(edge_id);
    let mut ends = [edge.start, edge.end].map(|v| v.map(|v| model.vertex(v).point));
    let reversed = edge.sense == Sense::Reversed;
    if reversed {
        ends.reverse();
    }
    let trace = catallax_intersect::trace(model.geom(), edge.curve, ends[0], ends[1], dec)?.value;
    let pieces = trace.pieces();
    if pieces.is_empty() {
        return Err(fit_failure(
            "a selected fillet edge has no retained trace pieces",
        ));
    }
    let at = |index: usize| {
        &pieces[if reversed {
            pieces.len() - 1 - index
        } else {
            index
        }]
    };
    let point = |piece: &catallax_intersect::TracePiece, share: f64, dec: &Decider| {
        let parameter = exact(if reversed { 1.0 - share } else { share })?;
        piece.evaluate(&Ball::from_rational(&parameter, BITS), dec)
    };
    let first = at(0).evaluate_exact(
        &if reversed {
            Rational::one()
        } else {
            Rational::zero()
        },
        dec,
    )?;
    let last = at(pieces.len() - 1).evaluate_exact(
        &if reversed {
            Rational::zero()
        } else {
            Rational::one()
        },
        dec,
    )?;
    let bound = tolerance.mul_pow2(-2);
    let proposal_bound = tolerance.to_f64_nearest() / 8.0;
    let mut fits = Vec::with_capacity(pieces.len());
    let mut samples = Vec::new();
    if DENSITIES[0] * pieces.len() > MAX_SPANS {
        return Err(fit_failure(
            "a fillet guide exceeds the trace fit's span limit",
        ));
    }
    // Keep the cheap cubic path. A smooth conic can exhaust its density
    // allowance, so also try quintics with the same samples and span cap.
    // Each attempt uses one common degree across all retained trace pieces.
    for degree in [3, 5] {
        fits.clear();
        samples.clear();
        'pieces: for index in 0..pieces.len() {
            let piece = at(index);
            let mut fitted = None;
            for density in DENSITIES {
                if density * pieces.len() > MAX_SPANS {
                    break;
                }
                if density < degree {
                    continue;
                }
                let mut parameters = Vec::with_capacity(density + 1);
                let mut points = Vec::with_capacity(density + 1);
                for j in 0..=density {
                    dec.charge_work(1)?;
                    let share = j as f64 / density as f64;
                    parameters.push((index as f64 + share) / pieces.len() as f64);
                    points.push(
                        point(piece, share, dec)?
                            .each_ref()
                            .map(|v| v.mid().to_f64_nearest()),
                    );
                }
                let (knots, controls) = if degree == 3 {
                    interpolate(&parameters, &points, dec)?
                } else {
                    nurbs::interpolate_quintic(&parameters, &points, dec)?
                };
                let mut acceptable = true;
                'checks: for j in 0..density {
                    for share in CHECK_SHARES {
                        dec.charge_work(1)?;
                        let local = (j as f64 + share) / density as f64;
                        let parameter = (index as f64 + local) / pieces.len() as f64;
                        let reference = point(piece, local, dec)?
                            .each_ref()
                            .map(|v| v.mid().to_f64_nearest());
                        if distance(nurbs::evaluate(&knots, &controls, parameter), reference)
                            > proposal_bound
                        {
                            acceptable = false;
                            break 'checks;
                        }
                    }
                }
                if acceptable {
                    fitted = Some((knots, controls, density));
                    break;
                }
            }
            let Some((knots, controls, density)) = fitted else {
                break 'pieces;
            };
            fits.push((knots, controls));
            for j in 0..density {
                for share in CHECK_SHARES {
                    let local = exact((j as f64 + share) / density as f64)?;
                    samples.push((index, local));
                }
            }
        }
        if fits.len() == pieces.len() {
            break;
        }
    }
    if fits.len() != pieces.len() {
        return Err(fit_failure(
            "a fillet guide cannot fit its retained trace within the declared tolerance and span allowance",
        ));
    }
    let coordinates = |p: &Point3| -> Result<[Rational; 3], Refusal> {
        Ok([scalar(&p.x, dec)?, scalar(&p.y, dec)?, scalar(&p.z, dec)?])
    };
    let (curve, clock) = guides::join(fits, [coordinates(&first)?, coordinates(&last)?], dec)?;
    let mut checked = Vec::with_capacity(samples.len());
    for (index, local) in samples {
        dec.charge_work(1)?;
        let width = &clock[index + 1] - &clock[index];
        let parameter = &clock[index] + &(&width * &local);
        let native = if reversed {
            &Rational::one() - &local
        } else {
            local
        };
        let reference = at(index).evaluate(&Ball::from_rational(&native, BITS), dec)?;
        let actual = curve.evaluate_box([&parameter, &parameter], BITS, dec)?;
        let delta: [Ball; 3] = std::array::from_fn(|axis| &actual[axis] - &reference[axis]);
        let error = (delta[0].square() + delta[1].square() + delta[2].square())
            .sqrt()
            .ok_or_else(|| fit_failure("fillet guide deviation cannot be enclosed"))?
            .hi();
        if error > bound {
            return Err(fit_failure(
                "the written fillet guide misses its retained trace at a checked sample",
            ));
        }
        checked.push(catallax_numbers::NumericalSample {
            parameters: vec![parameter],
            distance_upper_mm: error,
        });
    }
    Ok(FilletGuide {
        curve,
        ends: [first, last],
        check: catallax_numbers::NumericalCheck::Samples {
            name: "fillet.source-guide-sampled-deviation".into(),
            scope: format!("{edge_id:?}"),
            samples: checked,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::BSplineSurface;

    #[test]
    fn native_uv_follows_piecewise_carrier_parameters_without_extra_density() {
        let mut dec = Decider::default();
        let root_two = Real::from_i64(2).sqrt();
        let middle = &root_two - &Real::one();
        let half = Real::from(Rational::from_ratio(1, 2).unwrap());
        let weight = &(&Real::from_i64(2) + &root_two).sqrt() * &half;
        let diagonal = &Real::one() / &root_two;
        let radius = Real::from_i64(25);
        let rows = [
            [Real::one(), Real::zero()],
            [Real::one(), middle.clone()],
            [diagonal.clone(), diagonal],
            [middle, Real::one()],
            [Real::zero(), Real::one()],
        ];
        let third = Rational::from_ratio(1, 3).unwrap();
        let controls = rows
            .iter()
            .flat_map(|[x, y]| {
                [0, 10].map(|z| Point3::new(x * &radius, y * &radius, Real::from_i64(z)))
            })
            .collect();
        let source = BSplineSurface::new(
            [2, 1],
            [
                vec![Rational::zero(); 3]
                    .into_iter()
                    .chain([third.clone(), third])
                    .chain(vec![Rational::one(); 3])
                    .collect(),
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [5, 2],
            controls,
            [
                Real::one(),
                weight.clone(),
                Real::one(),
                weight,
                Real::one(),
            ]
            .into_iter()
            .flat_map(|w| [w.clone(), w])
            .collect(),
            &mut dec,
        )
        .unwrap();
        let face = SplineFace::new(
            Surface64::new(&source, &mut dec).unwrap(),
            std::array::from_fn(|_| [Rational::zero(), Rational::one()]),
        );
        let ends = [0, 1].map(|u| Uv {
            exact: [Rational::from_i64(u), Rational::zero()],
            approx: [u as f64, 0.0],
        });
        let end = Ball::pi(BITS).mid().mul_pow2(-1);
        let tolerance = Rational::parse("0.0000001").unwrap();
        let uv = native_uv(
            ParametricFace::Spline(&face),
            &Native::Circle(Circle64 {
                center: [0.0; 3],
                x: [1.0, 0.0, 0.0],
                y: [0.0, 1.0, 0.0],
                radius: 25.0,
            }),
            [&Rational::zero(), &end],
            [&ends[0], &ends[1]],
            &Tolerances {
                sample: 1e-7 / 64.0,
                fit: 2.5e-8,
                projection: tolerance,
            },
            &mut dec,
        )
        .unwrap();
        let mut maximum = 0.0_f64;
        for i in 0..=96 {
            let at = &end * &Rational::from_ratio(i, 96).unwrap();
            let parameters = uv.point_at(&Real::from(at.clone()), &mut dec).unwrap();
            let u = dec.enclose_at(&parameters.x, BITS).unwrap();
            let v = dec.enclose_at(&parameters.y, BITS).unwrap();
            let point = source
                .evaluate(&u, &v, &dec)
                .unwrap()
                .each_ref()
                .map(|v| v.mid().to_f64_nearest());
            let angle = at.to_f64_nearest();
            maximum = maximum.max(distance(
                point,
                [25.0 * angle.cos(), 25.0 * angle.sin(), 0.0],
            ));
        }
        assert!(
            maximum <= 1e-7,
            "carrier-chart physical residual={maximum:e} mm"
        );
    }

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

    #[test]
    fn half_circle_guides_keep_the_retained_branch_and_direction_at_fixed_density() {
        use catallax_shapes::{Circle, PointDef, RPoint3, RVector3};
        use catallax_topology::{Edge, Vertex};
        let (_, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            for lower in [false, true] {
                for sense in [Sense::Same, Sense::Reversed] {
                    let q = Rational::from_i64;
                    let mut model = Model::new();
                    let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                        center: RPoint3::new(q(0), q(0), q(0)).to_real(),
                        normal: RVector3::new(q(0), q(0), q(1)).to_real(),
                        ref_dir: RVector3::new(q(1), q(0), q(0)).to_real(),
                        radius_sq: q(400).into(),
                    }));
                    let vertices = [20, -20].map(|x| {
                        let point = model.geom_mut().add_point(PointDef::Explicit(
                            RPoint3::new(q(x), q(0), q(0)).to_real(),
                        ));
                        model.add_vertex(Vertex {
                            point,
                            tolerance: None,
                        })
                    });
                    let mut ends = if lower {
                        [vertices[1], vertices[0]]
                    } else {
                        vertices
                    };
                    if sense == Sense::Reversed {
                        ends.reverse();
                    }
                    let edge = model.add_edge(Edge {
                        curve,
                        sense,
                        start: Some(ends[0]),
                        end: Some(ends[1]),
                        tolerance: None,
                    });
                    let mut dec = Decider::default();
                    let tolerance = Rational::parse("0.0000001").unwrap();
                    let before = catallax_numbers::current_operation_work().used;
                    let guide = fillet_guide(&model, edge, &tolerance, &mut dec).unwrap();
                    let fitting_work = catallax_numbers::current_operation_work().used - before;
                    let mut worst: f64 = 0.0;
                    for i in 0..=96 {
                        let t = Rational::from_ratio(i, 96).unwrap();
                        let p = guide
                            .curve
                            .evaluate_box([&t, &t], BITS, &dec)
                            .unwrap()
                            .each_ref()
                            .map(|v| v.mid().to_f64_nearest());
                        worst = worst.max((p[0].hypot(p[1]) - 20.0).abs()).max(p[2].abs());
                        assert!(if lower { p[1] < 2.5e-8 } else { p[1] > -2.5e-8 });
                        if i == 0 || i == 96 {
                            let vertex = ends[usize::from(i == 96)];
                            let expected = model.geom().position(model.vertex(vertex).point);
                            assert!(distance(p, f64_point(&expected, &dec).unwrap()) <= 2.5e-8);
                        }
                    }
                    assert!(worst <= 2.5e-8, "radial/plane residual {worst}");
                    println!(
                        "half-guide lower={lower} sense={sense:?} work={fitting_work} radial_error_mm={worst}"
                    );
                }
            }
        });
        assert!(!work.exhausted);
    }
}
