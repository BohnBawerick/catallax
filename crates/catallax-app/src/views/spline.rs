//! Contour generators of spline faces.
//!
//! A face on an extruded Bezier profile gets exact lines (`extrusion.rs`).
//! A face on a finite rational patch, or a lofted face covered by such
//! patches, gets traced generators: `trace.rs` finds the zero set of the
//! contour function on each patch carrier the way SMLib's
//! `SmSurfaceSilhouette` does, and each traced generator is fitted with
//! `fit.rs` and marked approximate (decision Q14). OpenCASCADE has no closed
//! form here either: `Contap_Contour` marches such faces and
//! `HLRTopoBRep_DSFiller` approximates the marched points with
//! `BRepApprox_Approx`.
//!
//! The fit needs the true generator at any parameter, not only at traced
//! points. Along a stretch where the generator's surface parameter `u`
//! (or `v`) is monotone, the generator is the graph `v = w(u)` of the
//! implicit function `F(u, w(u)) = 0` (the implicit function theorem; Newton's
//! method solves it from the traced guess), with `w' = -F_u / F_v`. A traced
//! generator is cut where the other parameter takes over, and the fitted
//! pieces join in one cubic B-spline with a triple knot at each cut, which
//! keeps the joint point and lets the parameter speed change there (Piegl and
//! Tiller, "The NURBS Book", 2nd ed., section 5.2).
//!
//! Faces are trimmed in their own terms: a rational patch face must be the
//! complete patch, checked by `RationalPatchFace`, so the generators are
//! traced over the carrier's whole domain and end on its boundary rows. Other
//! spline faces refuse with `app.contour-trim` rather than give a partial
//! outline.

use std::cell::Cell;

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Curve, Point3, PointDef, SplineSupport, Surface, Vector3,
};
use catallax_topology::numerical::{gap, tolerance};
use catallax_topology::{FaceId, Model, RationalPatchFace};

use super::fit::{self, Fitted, Jet};
use super::float::{above, double, rational, unit_direction};
use super::{ContourArc, ContourFit};

mod bernstein;
mod extrusion;
mod trace;

use bernstein::Carrier;
use trace::{Branch, Tracer};

/// Source: judgement. A piece keeps its parameter while that parameter's
/// share of the tangent is at least half the other's, so the solve's
/// derivative `F_w` stays above `1 / sqrt(5)` of the gradient, and switches
/// only past that. Any factor above one avoids switching back and forth at
/// the diagonal; two leaves room for the [`trace`] turn of ten degrees.
const SWITCH_RATIO: f64 = 2.0;

/// Source: judgement. Newton's method on the traced guess converges
/// quadratically; sixteen steps only end one that stalls, and the fit then
/// measures whatever point it reached.
const SOLVE_STEPS: usize = 16;

fn trim(detail: &str) -> Refusal {
    Refusal::new("app.contour-trim", detail.to_owned())
}

/// The binary64 frame a carrier is traced and fitted in, with rule T1's
/// floor and rule T2's target.
struct Frame {
    /// Exact origin, the middle of the control box rounded to binary64.
    origin: [Rational; 3],
    floor: Rational,
    floor_mm: f64,
    target_mm: f64,
    /// Diagonal of the control box.
    diagonal: f64,
}

impl Frame {
    /// A positive-weight spline lies in the convex hull of its controls
    /// (Piegl and Tiller, section 4.2), so their box bounds the coordinate
    /// size rule T1 reads.
    fn new(surface: &BSplineSurface, dec: &mut Decider) -> Result<Self, Refusal> {
        let coarse = tolerance::floor_minimum_mm();
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        for point in surface.controls() {
            for (c, x) in [&point.x, &point.y, &point.z].into_iter().enumerate() {
                let x = double(x, &coarse, dec)?;
                low[c] = low[c].min(x);
                high[c] = high[c].max(x);
            }
        }
        let middle: [f64; 3] = std::array::from_fn(|c| low[c] + (high[c] - low[c]) / 2.0);
        let size = (0..3)
            .map(|c| (high[c] - middle[c]).max(middle[c] - low[c]))
            .fold(0.0_f64, f64::max)
            .next_up();
        let diagonal = (0..3)
            .map(|c| (high[c] - low[c]) * (high[c] - low[c]))
            .sum::<f64>()
            .sqrt();
        let floor = tolerance::floor_mm(&rational(size)?);
        let target = tolerance::own_fitted(&floor, None).target_mm;
        Ok(Self {
            origin: [
                rational(middle[0])?,
                rational(middle[1])?,
                rational(middle[2])?,
            ],
            floor_mm: floor.to_f64_nearest().next_up(),
            floor,
            target_mm: target.to_f64_nearest(),
            diagonal,
        })
    }

    fn point(&self, local: [f64; 3]) -> Result<Point3, Refusal> {
        Ok(Point3::new(
            Real::from(&self.origin[0] + &rational(local[0])?),
            Real::from(&self.origin[1] + &rational(local[1])?),
            Real::from(&self.origin[2] + &rational(local[2])?),
        ))
    }
}

/// A stretch of a traced generator monotone in one surface parameter.
struct Piece {
    /// The parameter followed: 0 for `u`.
    axis: usize,
    /// `+1` when it increases along the trace.
    sense: f64,
    from: f64,
    to: f64,
    /// Fit parameter at `from`; it then runs at unit speed in `|axis|`.
    start: f64,
    /// Traced points as `(followed, other)`, ascending in the first.
    guide: Vec<[f64; 2]>,
}

impl Piece {
    fn end(&self) -> f64 {
        self.start + (self.to - self.from).abs()
    }

    /// The generator point at fit parameter `t`, as `(u, v)`.
    fn solve(&self, carrier: &Carrier, t: f64, resolution: f64, count: &Cell<u64>) -> [f64; 2] {
        let c = self.from + self.sense * (t - self.start);
        let g = &self.guide;
        let i = g.partition_point(|x| x[0] <= c).clamp(1, g.len() - 1);
        let (a, b) = (g[i - 1], g[i]);
        let guess = if b[0] > a[0] {
            a[1] + (c - a[0]) / (b[0] - a[0]) * (b[1] - a[1])
        } else {
            a[1]
        };
        let other = 1 - self.axis;
        let mut uv = [0.0; 2];
        uv[self.axis] = c;
        uv[other] = guess;
        for _ in 0..SOLVE_STEPS {
            count.set(count.get() + 1);
            let level = carrier.level(uv[0], uv[1]);
            let slope = level.grad[other];
            if !above(slope.abs(), 0.0) {
                break;
            }
            let step = level.value / slope;
            uv[other] -= step;
            if !above(step.abs(), resolution) {
                break;
            }
        }
        uv
    }

    /// Point and derivative in the fit parameter.
    fn jet(&self, carrier: &Carrier, t: f64, resolution: f64, count: &Cell<u64>) -> Jet {
        let uv = self.solve(carrier, t, resolution, count);
        let jet = carrier.jet(uv[0], uv[1]);
        let level = carrier.level(uv[0], uv[1]);
        let other = 1 - self.axis;
        let slope = -level.grad[self.axis] / level.grad[other];
        let (along, across) = if self.axis == 0 {
            (jet.du, jet.dv)
        } else {
            (jet.dv, jet.du)
        };
        (
            jet.point,
            std::array::from_fn(|k| self.sense * (along[k] + across[k] * slope)),
        )
    }
}

/// Cuts a traced generator where its followed parameter changes.
fn pieces(branch: &Branch) -> Result<Vec<Piece>, Refusal> {
    let v = &branch.vertices;
    let dominant = |t: [f64; 2]| usize::from(t[0].abs() < t[1].abs());
    let mut cuts = vec![(dominant(v[0].tangent), 0)];
    for (i, vertex) in v.iter().enumerate().take(v.len() - 1).skip(1) {
        let axis = cuts[cuts.len() - 1].0;
        let t = vertex.tangent;
        if SWITCH_RATIO * t[axis].abs() < t[1 - axis].abs() {
            cuts.push((1 - axis, i));
        }
    }
    let mut out = Vec::with_capacity(cuts.len());
    let mut start = 0.0;
    for (n, &(axis, first)) in cuts.iter().enumerate() {
        let last = cuts.get(n + 1).map_or(v.len() - 1, |c| c.1);
        let (from, to) = (v[first].uv[axis], v[last].uv[axis]);
        let sense = if to > from { 1.0 } else { -1.0 };
        let mut guide: Vec<[f64; 2]> = v[first..=last]
            .iter()
            .map(|x| [x.uv[axis], x.uv[1 - axis]])
            .collect();
        if sense < 0.0 {
            guide.reverse();
        }
        if guide.windows(2).any(|w| !above(w[1][0], w[0][0])) {
            return Err(Refusal::new(
                "app.contour-trace",
                "a traced contour generator is not monotone where it was cut",
            ));
        }
        let piece = Piece {
            axis,
            sense,
            from,
            to,
            start,
            guide,
        };
        start = piece.end();
        out.push(piece);
    }
    Ok(out)
}

/// Fits one traced generator and records its arc.
#[allow(clippy::too_many_arguments)]
fn push_traced(
    scratch: &mut Model,
    face: FaceId,
    carrier: &Carrier,
    frame: &Frame,
    branch: &Branch,
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let pieces = pieces(branch)?;
    let count = Cell::new(0_u64);
    let mut fitted: Vec<Fitted> = Vec::with_capacity(pieces.len());
    let mut first: Option<[f64; 3]> = None;
    for (n, piece) in pieces.iter().enumerate() {
        dec.check_deadline()?;
        let resolution = trace::resolution(carrier, 1 - piece.axis);
        let generator = |t: f64| piece.jet(carrier, t, resolution, &count);
        let mut from = generator(piece.start);
        if let Some(previous) = fitted.last() {
            from.0 = previous.controls[previous.controls.len() - 1];
        }
        first.get_or_insert(from.0);
        let mut to = generator(piece.end());
        if branch.closed && n + 1 == pieces.len() {
            to.0 = first.expect("a first piece");
        }
        let fit = fit::fit(
            generator,
            (piece.start, from),
            (piece.end(), to),
            frame.target_mm,
            frame.floor_mm,
            gap::EDGE_WORK,
            dec,
        )?;
        fitted.push(fit);
    }
    dec.charge_work(count.get())?;
    // One cubic B-spline: double knots inside each piece, triple at cuts.
    let mut knots = Vec::new();
    let mut controls = Vec::new();
    for (n, fit) in fitted.iter().enumerate() {
        let breaks = &fit.breaks;
        if n == 0 {
            knots.extend(std::iter::repeat_n(rational(breaks[0])?, 4));
            controls.extend_from_slice(&fit.controls[..fit.controls.len() - 1]);
        } else {
            knots.extend(std::iter::repeat_n(rational(breaks[0])?, 3));
            let previous = &fitted[n - 1].controls;
            controls.push(previous[previous.len() - 1]);
            controls.extend_from_slice(&fit.controls[1..fit.controls.len() - 1]);
        }
        for &b in &breaks[1..breaks.len() - 1] {
            knots.extend(std::iter::repeat_n(rational(b)?, 2));
        }
    }
    let last = &fitted[fitted.len() - 1];
    knots.extend(std::iter::repeat_n(
        rational(last.breaks[last.breaks.len() - 1])?,
        4,
    ));
    controls.push(last.controls[last.controls.len() - 1]);
    let points = controls
        .iter()
        .map(|&local| frame.point(local))
        .collect::<Result<Vec<_>, _>>()?;
    let ends = if branch.closed {
        (None, None)
    } else {
        let geom = scratch.geom_mut();
        (
            Some(geom.add_point(PointDef::Explicit(points[0].clone()))),
            Some(geom.add_point(PointDef::Explicit(points[points.len() - 1].clone()))),
        )
    };
    let weights = vec![Real::one(); points.len()];
    let curve = BSplineCurve::new(3, knots, points, weights, dec)?;
    let curve = scratch.geom_mut().add_curve(Curve::BSpline(curve));
    let deviation = fitted
        .iter()
        .map(|f| f.deviation_mm)
        .fold(0.0_f64, f64::max);
    arcs.push(ContourArc {
        face,
        curve,
        start: ends.0,
        end: ends.1,
        on_faces: vec![face],
        fit: Some(ContourFit {
            deviation_mm: deviation,
            tolerance_mm: (gap::MARGIN * deviation).next_up(),
            target_mm: frame.target_mm,
            spans: fitted.iter().map(|f| f.breaks.len() - 1).sum(),
        }),
    });
    Ok(())
}

/// The traced generators of one complete patch carrier, for `face`.
fn traced(
    scratch: &mut Model,
    face: FaceId,
    surface: &BSplineSurface,
    view: &Vector3,
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let frame = Frame::new(surface, dec)?;
    let d = unit_direction(view, &frame.floor, dec)?;
    let origin = [&frame.origin[0], &frame.origin[1], &frame.origin[2]];
    let (carrier, work) = Carrier::new(surface, origin, d, &frame.floor, dec)?;
    dec.charge_work(work)?;
    let cap = gap::EDGE_WORK;
    let (mut tracer, work) = Tracer::new(&carrier, frame.diagonal, frame.floor_mm, cap, dec)?;
    let (branches, _) = tracer.branches(cap.saturating_sub(work), dec)?;
    for branch in &branches {
        push_traced(scratch, face, &carrier, &frame, branch, arcs, dec)?;
    }
    Ok(())
}

/// The generators of one spline face for the view direction `view`.
pub(super) fn generators(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    view: &Vector3,
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let Surface::BSpline(support) = scratch.geom().surface(scratch.face(face).surface) else {
        unreachable!("only spline faces are traced");
    };
    match support {
        SplineSupport::BezierExtrusion(_) => {
            extrusion::generators(scratch, body_faces, face, view, arcs, dec)
        }
        SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_) => {
            let cells = scratch
                .face_chart_cover(face)
                .map_or_else(|| vec![face], <[FaceId]>::to_vec);
            for cell in cells {
                let carrier = RationalPatchFace::of_chart_face(scratch, cell, dec)?
                    .ok_or_else(|| {
                        trim(
                            "contour generators of a rational spline face need complete \
                             patch charts",
                        )
                    })?
                    .support()
                    .carrier()
                    .clone();
                traced(scratch, face, &carrier, view, arcs, dec)?;
            }
            Ok(())
        }
        _ => Err(trim(
            "contour generators of fillet, draft and other spline bands need that band's \
             own trim; unfinished",
        )),
    }
}
