//! Rule T5's pair rule, and the small features rules T2 and T8 read.
//!
//! Two entities that share no topology are a pair: two different vertices, a
//! vertex and an edge that does not end at it, two edges with no vertex in
//! common, and a vertex or an edge and a face whose boundary it does not
//! touch. Their boxes go in a bounding-box tree, each widened by its reach,
//! and only the pairs whose widened boxes overlap are measured. On a part
//! with no small features the tree returns no pair and nothing is measured.
//! Every box test, sample, evaluation and foot search is charged as work.
//!
//! A face of a numerical body stores no tolerance (its surface is the
//! authoritative copy), so a pair with a face has `t_face = 0`. Its distance
//! is read to the face's interior: the foot of a point on the surface,
//! where that foot lies inside the face's loops. Nearer the boundary, the
//! face's own edges and vertices are the pair. Two faces both have `t = 0`,
//! so rule T5 cannot fire on them, and they are not paired.
//!
//! A face's loops are read in its parameters as chains of nodes. A use that
//! is straight in u/v is one node and a straight run to the next. A use with
//! a retained u/v curve, and a circle or spline on a plane (admission's
//! projected copy), is walked: at nine points, and at every place between
//! them where its u or its v turns, so that both run one way between two
//! nodes. A point is inside by the even-odd count of the runs a line of
//! constant v meets beyond it. A curved run is counted by its two ends where
//! the point lies outside their box, and on the curve itself where it lies
//! inside, so the zone between a chord and its arc is read on the arc.
//! OpenCASCADE's `BRepTopAdaptor_FClass2d` takes the same two steps: a
//! sampled polygon with the chord sag as its tolerance, then the true curves
//! (`BRepClass_FaceClassifier`) for a point inside that tolerance. SMLib's
//! `SmLoop::ContainsUVPoint` marches the true trim curve near the boundary.
//!
//! Distances are read in binary64 in the body's frame (rule T1), by the
//! sampling rule of the edge-gap checker (section 4.2 step 4) and Brent's
//! method. Two curves are then polished together by the Gauss-Newton of the
//! surface foot, and a curve that changes side of a face between two samples
//! is followed to the crossing by halving. They are searches for the least
//! distance, as the checker's are searches for the largest gap: the found
//! value is at least the true one.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Rational, Refusal};

use super::{
    EntityKind, EntityTolerance, PAIR_INFLATION, PairVerdict, fitted_feature_reach_mm, own_fitted,
    pair_verdict,
};
use crate::Sense;
use crate::numerical::gap::{
    FramedCurve, FramedSurface, FramedUvCurve, ROOT_DEPTH, brent_maximum, chebyshev_shares,
    nearest_on_curves, try_brent_maximum,
};
use crate::numerical::{
    NumericalBody, NumericalDefinition, NumericalNativeTrim, analytic, native::Index,
};

/// Two entities that share no topology, closer than a reach.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct FeaturePair {
    /// The first entity: its kind and key.
    pub a: (EntityKind, String),
    /// The second entity.
    pub b: (EntityKind, String),
    /// The least distance found between them, in millimetres.
    pub distance_mm: Rational,
    /// The sum of their tolerances, in millimetres.
    pub tolerance_sum_mm: Rational,
}

/// Source: Ericson, "Real-Time Collision Detection", 2005, section 6.2: a
/// leaf holds a few primitives so the tree stays shallow. Four, as the mesh
/// tree of `catallax-import` uses.
const LEAF: usize = 4;

/// Source: the edge-gap checker's sample count, section 4.2 step 4 of the
/// Checkpoint 2 decisions: `max(17, 3 x degree)` points on each piece.
const LEAST_POINTS: usize = 17;
/// See [`LEAST_POINTS`].
const POINTS_PER_DEGREE: usize = 3;

/// Source: the 9 by 9 chart grid that admission samples on a numerical
/// carrier (`NumericalBody`'s documentation). A face's u/v box is sampled on
/// it to seed the foot of a point on the face.
const FACE_GRID: usize = 9;

/// Source: the same allowance, nine points, admission gives the u/v walk of
/// each edge use on a numerical carrier. A curved copy of an edge use is
/// read at these points when a face's loops are written as nodes.
const LOOP_POINTS: usize = 9;

/// Source: judgement. Each round of alternating Brent searches on two
/// curves moves both feet nearer; three rounds settle a transversal pair.
const ALTERNATIONS: usize = 3;

type Bounds = [[f64; 2]; 3];

/// Two entity indices and the distance found between them.
type ClosePair = (usize, usize, f64);

fn widened(bounds: &Bounds, by: f64) -> Bounds {
    bounds.map(|[lo, hi]| [lo - by, hi + by])
}

fn overlap(a: &Bounds, b: &Bounds) -> bool {
    (0..3).all(|axis| a[axis][0] <= b[axis][1] && b[axis][0] <= a[axis][1])
}

fn box_distance(a: &Bounds, b: &Bounds) -> f64 {
    let mut sum = 0.0;
    for axis in 0..3 {
        let gap = (b[axis][0] - a[axis][1])
            .max(a[axis][0] - b[axis][1])
            .max(0.0);
        sum += gap * gap;
    }
    sum.sqrt()
}

fn union(a: &Bounds, b: &Bounds) -> Bounds {
    std::array::from_fn(|axis| [a[axis][0].min(b[axis][0]), a[axis][1].max(b[axis][1])])
}

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Whether `t` lies strictly between two parameters, given in either order.
/// A halving stops where it does not: no binary64 number is left between.
fn between(t: f64, a: f64, b: f64) -> bool {
    a.min(b) < t && t < a.max(b)
}

/// A binary bounding-box tree over widened entity boxes, split at the median
/// of the box centres along the longest axis (Ericson 2005, section 6.1).
enum Node {
    Leaf(Bounds, Vec<usize>),
    Inner(Bounds, Box<Node>, Box<Node>),
}

impl Node {
    fn bounds(&self) -> &Bounds {
        match self {
            Self::Leaf(bounds, _) | Self::Inner(bounds, ..) => bounds,
        }
    }

    fn build(items: &mut [usize], boxes: &[Bounds], dec: &mut Decider) -> Result<Self, Refusal> {
        dec.charge_work(items.len() as u64)?;
        let bounds = items
            .iter()
            .skip(1)
            .fold(boxes[items[0]], |acc, i| union(&acc, &boxes[*i]));
        if items.len() <= LEAF {
            return Ok(Self::Leaf(bounds, items.to_vec()));
        }
        let axis = (0..3)
            .max_by(|a, b| {
                (bounds[*a][1] - bounds[*a][0]).total_cmp(&(bounds[*b][1] - bounds[*b][0]))
            })
            .expect("three axes");
        let centre = |i: &usize| boxes[*i][axis][0] + boxes[*i][axis][1];
        items.sort_by(|a, b| centre(a).total_cmp(&centre(b)).then(a.cmp(b)));
        let (left, right) = items.split_at_mut(items.len() / 2);
        Ok(Self::Inner(
            bounds,
            Box::new(Self::build(left, boxes, dec)?),
            Box::new(Self::build(right, boxes, dec)?),
        ))
    }

    fn self_pairs(
        &self,
        boxes: &[Bounds],
        out: &mut Vec<(usize, usize)>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        match self {
            Self::Leaf(_, items) => {
                for (k, i) in items.iter().enumerate() {
                    for j in &items[k + 1..] {
                        dec.charge_work(1)?;
                        if overlap(&boxes[*i], &boxes[*j]) {
                            out.push(((*i).min(*j), (*i).max(*j)));
                        }
                    }
                }
                Ok(())
            }
            Self::Inner(_, left, right) => {
                left.self_pairs(boxes, out, dec)?;
                right.self_pairs(boxes, out, dec)?;
                left.cross(right, boxes, out, dec)
            }
        }
    }

    fn cross(
        &self,
        other: &Self,
        boxes: &[Bounds],
        out: &mut Vec<(usize, usize)>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        dec.charge_work(1)?;
        if !overlap(self.bounds(), other.bounds()) {
            return Ok(());
        }
        match (self, other) {
            (Self::Leaf(_, a), Self::Leaf(_, b)) => {
                for i in a {
                    for j in b {
                        dec.charge_work(1)?;
                        if overlap(&boxes[*i], &boxes[*j]) {
                            out.push(((*i).min(*j), (*i).max(*j)));
                        }
                    }
                }
                Ok(())
            }
            (Self::Inner(_, left, right), _) => {
                left.cross(other, boxes, out, dec)?;
                right.cross(other, boxes, out, dec)
            }
            (Self::Leaf(..), Self::Inner(_, left, right)) => {
                self.cross(left, boxes, out, dec)?;
                self.cross(right, boxes, out, dec)
            }
        }
    }
}

/// The pieces of a framed curve between its breaks, each with its sample
/// points and a box.
struct Piece {
    samples: Vec<(f64, [f64; 3])>,
    bounds: Bounds,
}

enum Shape {
    Point([f64; 3]),
    /// An edge, by its index in the edge table. Its pieces are sampled only
    /// when a pair needs them, so a part with no candidate pair samples
    /// nothing.
    Curve(usize),
    /// A face, by its index in the face table.
    Face(usize),
}

struct Entity {
    kind: EntityKind,
    key: String,
    tolerance: f64,
    shape: Shape,
    /// Vertex indices this entity is, ends at, or is bounded by, sorted.
    vertices: Vec<usize>,
    bounds: Bounds,
}

/// A curved copy of an edge use in its face's parameters.
enum Walk {
    /// A retained u/v curve.
    Retained(FramedUvCurve),
    /// The 3D curve of a circle or spline on a plane, read through the
    /// plane's chart: admission's projected copy.
    Projected(FramedCurve),
}

impl Walk {
    /// The copy's u/v at the edge parameter `t`.
    fn at(&mut self, surface: &FramedSurface, t: f64) -> [f64; 2] {
        match self {
            Self::Retained(curve) => curve.point(t),
            Self::Projected(curve) => surface
                .plane_parameters(&curve.point(t))
                .expect("a projected copy lies on a plane"),
        }
    }
}

/// One point of a face's loop in u/v, and how the loop runs on to the next.
struct LoopNode {
    uv: [f64; 2],
    /// Where the run to the next node follows a curved copy: the copy's
    /// index among the face's walks and its parameters at the two nodes.
    /// `None` for a straight run.
    along: Option<(usize, [f64; 2])>,
}

/// A face as the pair rule reads it: its surface, its u/v box and its loops
/// as chains of nodes in u/v, made when a pair first needs them.
struct FaceData {
    surface: FramedSurface,
    rect: [[f64; 2]; 2],
    loops: Vec<Vec<LoopNode>>,
    walks: Vec<Walk>,
}

impl FaceData {
    /// Whether `uv` lies inside the face's loops, by the even-odd rule, with
    /// u shifted by a period on a surface of revolution.
    ///
    /// A run counts where it meets the line of constant v beyond the point.
    /// A straight run is its chord. A curved run goes one way in u and in v
    /// between its nodes, so it lies in the box of its two ends: a point
    /// short of that box has the meeting beyond it, a point past the box has
    /// not, and for a point inside the box the meeting is found on the copy
    /// itself.
    fn inside(&mut self, uv: [f64; 2], dec: &mut Decider) -> Result<bool, Refusal> {
        let Self {
            surface,
            loops,
            walks,
            ..
        } = self;
        let shifts: &[f64] = match surface.period_u {
            Some(p) => &[0.0, -p, p],
            None => &[0.0],
        };
        for shift in shifts {
            let (u, v) = (uv[0] + shift, uv[1]);
            let mut inside = false;
            for nodes in loops.iter() {
                dec.charge_work(nodes.len() as u64)?;
                for (k, node) in nodes.iter().enumerate() {
                    let (a, b) = (node.uv, nodes[(k + 1) % nodes.len()].uv);
                    if (a[1] > v) == (b[1] > v) {
                        continue;
                    }
                    let beyond = match node.along {
                        None => u < a[0] + (v - a[1]) * (b[0] - a[0]) / (b[1] - a[1]),
                        Some(_) if u < a[0].min(b[0]) => true,
                        Some(_) if u >= a[0].max(b[0]) => false,
                        Some((walk, ends)) => {
                            u < meeting(&mut walks[walk], surface, ends, a[1] > v, v, dec)?
                        }
                    };
                    if beyond {
                        inside = !inside;
                    }
                }
            }
            if inside {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// The u at which a curved run meets the line of constant `v`. The run
/// goes one way in v from the parameter `from`, where it is `above` the
/// line or not, to `to`, where it is on the other side, so the meeting is
/// its one change of side. It is found by halving, [`ROOT_DEPTH`] times.
fn meeting(
    walk: &mut Walk,
    surface: &FramedSurface,
    [mut from, mut to]: [f64; 2],
    above: bool,
    v: f64,
    dec: &mut Decider,
) -> Result<f64, Refusal> {
    let mut u = None;
    for _ in 0..ROOT_DEPTH {
        let middle = from + (to - from) / 2.0;
        if !between(middle, from, to) {
            break;
        }
        dec.charge_work(1)?;
        let at = walk.at(surface, middle);
        if (at[1] > v) == above {
            from = middle;
        } else {
            to = middle;
        }
        u = Some(at[0]);
    }
    match u {
        Some(u) => Ok(u),
        None => {
            dec.charge_work(1)?;
            Ok(walk.at(surface, from)[0])
        }
    }
}

fn pieces(curve: &mut FramedCurve, dec: &mut Decider) -> Result<Vec<Piece>, Refusal> {
    let count = LEAST_POINTS.max(POINTS_PER_DEGREE * curve.degree);
    let count = count + count % 2;
    let shares = chebyshev_shares(count);
    let whole = curve.bounding_box();
    let breaks = curve.breaks.clone();
    let mut out = Vec::with_capacity(breaks.len() - 1);
    for window in breaks.windows(2) {
        let [a, b] = [window[0], window[1]];
        let mut samples = Vec::with_capacity(count + 2);
        let mut ts = vec![a];
        ts.extend(shares.iter().map(|s| a + s * (b - a)));
        ts.push(b);
        let mut bounds: Bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
        for t in ts {
            dec.charge_work(1)?;
            let p = curve.point(t);
            for axis in 0..3 {
                bounds[axis][0] = bounds[axis][0].min(p[axis]);
                bounds[axis][1] = bounds[axis][1].max(p[axis]);
            }
            samples.push((t, p));
        }
        // The samples' box can miss a bulge between them; the whole curve's
        // box cannot. Each piece keeps the samples' box widened by the
        // largest sample spacing, clipped to the whole box. The prune that
        // reads it only skips piece pairs, never entity pairs.
        let spacing = samples
            .windows(2)
            .map(|w| distance(&w[0].1, &w[1].1))
            .fold(0.0, f64::max);
        let bounds: Bounds = std::array::from_fn(|axis| {
            [
                (bounds[axis][0] - spacing).max(whole[axis][0]),
                (bounds[axis][1] + spacing).min(whole[axis][1]),
            ]
        });
        out.push(Piece { samples, bounds });
    }
    Ok(out)
}

/// The least distance from `point` to a curve's pieces within `reach`, by
/// the samples and a Brent polish of every sampled local minimum.
fn point_to_curve(
    point: &[f64; 3],
    curve: &mut FramedCurve,
    pieces: &[Piece],
    reach: f64,
    dec: &mut Decider,
) -> Result<Option<(f64, f64)>, Refusal> {
    let target = [[point[0]; 2], [point[1]; 2], [point[2]; 2]];
    let mut best: Option<(f64, f64)> = None;
    for piece in pieces {
        if box_distance(&piece.bounds, &target) > reach {
            continue;
        }
        let values: Vec<f64> = piece
            .samples
            .iter()
            .map(|(_, p)| distance(p, point))
            .collect();
        for (i, value) in values.iter().enumerate() {
            let left = i.checked_sub(1).map(|k| values[k]);
            let right = values.get(i + 1).copied();
            if left.is_some_and(|l| l < *value) || right.is_some_and(|r| r < *value) {
                continue;
            }
            let mut found = (*value, piece.samples[i].0);
            if let (Some(_), Some(_)) = (left, right) {
                let [a, b] = [piece.samples[i - 1].0, piece.samples[i + 1].0];
                let (share, negative, evaluations) =
                    brent_maximum(&mut |s| -distance(&curve.point(a + s * (b - a)), point));
                dec.charge_work(evaluations)?;
                if -negative < found.0 {
                    found = (-negative, a + share * (b - a));
                }
            }
            if best.is_none_or(|(d, _)| found.0 < d) {
                best = Some(found);
            }
        }
    }
    Ok(best)
}

/// The foot of `point` on a face's surface, from `seed` or, without one,
/// from the nearest point of a grid on the face's u/v box. `None` where the
/// search does not end at a finite u/v. The second value seeds the next
/// point along a curve: the foot, or the seed where there is no foot.
fn foot_on(
    point: &[f64; 3],
    face: &mut FaceData,
    seed: Option<[f64; 2]>,
    dec: &mut Decider,
) -> Result<(Option<[f64; 2]>, [f64; 2]), Refusal> {
    let seed = match seed {
        Some(seed) => seed,
        None => {
            let [[u0, u1], [v0, v1]] = face.rect;
            let mut best = (f64::INFINITY, [u0, v0]);
            dec.charge_work((FACE_GRID * FACE_GRID) as u64)?;
            for i in 0..FACE_GRID {
                for j in 0..FACE_GRID {
                    let uv = [
                        u0 + (u1 - u0) * i as f64 / (FACE_GRID - 1) as f64,
                        v0 + (v1 - v0) * j as f64 / (FACE_GRID - 1) as f64,
                    ];
                    let d = distance(&face.surface.point(uv), point);
                    if d < best.0 {
                        best = (d, uv);
                    }
                }
            }
            best.1
        }
    };
    // One more than the search spent: the surface's point at the foot, which
    // every caller reads next.
    let (foot, evaluations) = face.surface.foot(point, seed);
    dec.charge_work(evaluations + 1)?;
    if !foot.iter().all(|x| x.is_finite()) {
        return Ok((None, seed));
    }
    Ok((Some(foot), foot))
}

/// The distance from `point` to the interior of a face: its foot on the
/// surface, where that foot lies inside the face. `None` where it does not:
/// the face's boundary is read through its edges and vertices.
fn point_to_face(
    point: &[f64; 3],
    face: &mut FaceData,
    dec: &mut Decider,
) -> Result<Option<f64>, Refusal> {
    let (Some(foot), _) = foot_on(point, face, None, dec)? else {
        return Ok(None);
    };
    if !face.inside(foot, dec)? {
        return Ok(None);
    }
    Ok(Some(distance(&face.surface.point(foot), point)))
}

/// The least distance from a curve's pieces within `reach` of a face's box
/// to the face's interior.
///
/// A sample's foot is sought where a run of the curve next to it can come
/// within reach of the face's box: the box of the run's two samples, widened
/// by their spacing as a piece's box is. The nearest sample whose foot lies
/// inside the face is polished by Brent's method between its neighbours.
///
/// Two consecutive samples whose offsets from their feet point to opposite
/// sides of the surface have a crossing between them, however far from the
/// surface both are. It is followed by halving the run, [`ROOT_DEPTH`]
/// times, and it counts where its foot lies inside the face. OpenCASCADE's
/// `IntCurveSurface_Inter` brackets a crossing the same way, by a segment of
/// the curve's polygon that passes through the surface's polyhedron, and
/// then solves `S(u, v) = C(t)` by Newton's method, as SMLib's
/// `SmSurface::LocalCurveIntersect` does. Halving keeps the bracket and
/// needs no derivative of the curve.
fn curve_to_face(
    curve: &mut FramedCurve,
    pieces: &[Piece],
    face: &mut FaceData,
    face_box: &Bounds,
    reach: f64,
    dec: &mut Decider,
) -> Result<Option<f64>, Refusal> {
    let mut best: Option<(f64, [f64; 2], [f64; 2])> = None;
    // The runs that change side: their two parameters, the offset at the
    // first and the foot at the second.
    let mut crossings: Vec<([f64; 2], [f64; 3], [f64; 2])> = Vec::new();
    let mut seed = None;
    for piece in pieces {
        if box_distance(&piece.bounds, face_box) > reach {
            continue;
        }
        let near: Vec<bool> = piece
            .samples
            .windows(2)
            .map(|run| {
                let (p, q) = (&run[0].1, &run[1].1);
                let held: Bounds =
                    std::array::from_fn(|axis| [p[axis].min(q[axis]), p[axis].max(q[axis])]);
                box_distance(&widened(&held, distance(p, q)), face_box) <= reach
            })
            .collect();
        let mut before: Option<(f64, [f64; 3])> = None;
        for (k, (t, p)) in piece.samples.iter().enumerate() {
            // A sample with no run next to it that can come within reach of
            // the face's box has no foot sought.
            dec.charge_work(1)?;
            if !((k > 0 && near[k - 1]) || near.get(k).is_some_and(|near| *near)) {
                before = None;
                continue;
            }
            let (foot, next) = foot_on(p, face, seed, dec)?;
            seed = Some(next);
            let Some(foot) = foot else {
                before = None;
                continue;
            };
            let on = face.surface.point(foot);
            let offset: [f64; 3] = std::array::from_fn(|axis| p[axis] - on[axis]);
            if let Some((t_before, side)) = before
                && dot(&side, &offset) < 0.0
            {
                crossings.push(([t_before, *t], side, foot));
            }
            before = Some((*t, offset));
            if !face.inside(foot, dec)? {
                continue;
            }
            let d = distance(&on, p);
            let bracket = [
                piece.samples[k.saturating_sub(1)].0,
                piece.samples[(k + 1).min(piece.samples.len() - 1)].0,
            ];
            if best.is_none_or(|(b, ..)| d < b) {
                best = Some((d, foot, bracket));
            }
        }
    }
    let mut least: Option<f64> = None;
    if let Some((value, foot, [a, b])) = best {
        let (_, negative, _) = try_brent_maximum(|x| -> Result<f64, Refusal> {
            dec.charge_work(3)?;
            let p = curve.point(a + x * (b - a));
            let (uv, n) = face.surface.foot(&p, foot);
            dec.charge_work(n)?;
            if face.inside(uv, dec)? {
                Ok(-distance(&face.surface.point(uv), &p))
            } else {
                Ok(-value)
            }
        })?;
        least = Some(value.min(-negative));
    }
    for ([mut from, mut to], side, mut foot) in crossings {
        let mut found = None;
        for _ in 0..ROOT_DEPTH {
            let middle = from + (to - from) / 2.0;
            if !between(middle, from, to) {
                break;
            }
            let p = curve.point(middle);
            let (uv, evaluations) = face.surface.foot(&p, foot);
            // The curve's point, the search for its foot and the surface's
            // point there.
            dec.charge_work(evaluations + 2)?;
            if !uv.iter().all(|x| x.is_finite()) {
                break;
            }
            foot = uv;
            let on = face.surface.point(uv);
            let offset: [f64; 3] = std::array::from_fn(|axis| p[axis] - on[axis]);
            if dot(&side, &offset) > 0.0 {
                from = middle;
            } else {
                to = middle;
            }
            found = Some((uv, distance(&on, &p)));
        }
        if let Some((uv, d)) = found
            && face.inside(uv, dec)?
        {
            least = Some(least.map_or(d, |known| known.min(d)));
        }
    }
    Ok(least)
}

/// The least distance between two curves within `reach`: the nearest sample
/// pair of every pair of pieces whose boxes are within reach, polished by
/// alternating Brent searches on each curve between that pair's neighbours,
/// and then by [`nearest_on_curves`] over the whole of both edges. The
/// alternation nears a crossing by a constant factor a round and cannot
/// leave its first brackets; the second polish moves both parameters at
/// once and settles on a crossing to rounding. Where it leaves the range of
/// a closed edge, one that starts and ends at the same vertex, it has
/// passed that vertex: the parameter is taken round to the other end and
/// the polish runs once more from there.
fn curve_to_curve(
    a: (&mut FramedCurve, &[Piece], bool),
    b: (&mut FramedCurve, &[Piece], bool),
    reach: f64,
    dec: &mut Decider,
) -> Result<Option<f64>, Refusal> {
    let (curve_a, pieces_a, closed_a) = a;
    let (curve_b, pieces_b, closed_b) = b;
    let mut best: Option<f64> = None;
    for piece_a in pieces_a {
        for piece_b in pieces_b {
            dec.charge_work(1)?;
            if box_distance(&piece_a.bounds, &piece_b.bounds) > reach {
                continue;
            }
            let mut nearest = (f64::INFINITY, 0, 0);
            for (i, (_, p)) in piece_a.samples.iter().enumerate() {
                dec.charge_work(piece_b.samples.len() as u64)?;
                for (j, (_, q)) in piece_b.samples.iter().enumerate() {
                    let d = distance(p, q);
                    if d < nearest.0 {
                        nearest = (d, i, j);
                    }
                }
            }
            let bracket = |piece: &Piece, i: usize| {
                [
                    piece.samples[i.saturating_sub(1)].0,
                    piece.samples[(i + 1).min(piece.samples.len() - 1)].0,
                ]
            };
            let [a0, a1] = bracket(piece_a, nearest.1);
            let [b0, b1] = bracket(piece_b, nearest.2);
            let (mut s, mut t) = (piece_a.samples[nearest.1].0, piece_b.samples[nearest.2].0);
            let mut value = nearest.0;
            for _ in 0..ALTERNATIONS {
                let p = curve_a.point(s);
                let (share, negative, evaluations) =
                    brent_maximum(&mut |x| -distance(&curve_b.point(b0 + x * (b1 - b0)), &p));
                dec.charge_work(evaluations)?;
                if -negative <= value {
                    (t, value) = (b0 + share * (b1 - b0), -negative);
                }
                let q = curve_b.point(t);
                let (share, negative, evaluations) =
                    brent_maximum(&mut |x| -distance(&curve_a.point(a0 + x * (a1 - a0)), &q));
                dec.charge_work(evaluations)?;
                if -negative <= value {
                    (s, value) = (a0 + share * (a1 - a0), -negative);
                }
            }
            let mut at = [s, t];
            // The polish, and once more where a closed edge is re-entered.
            for _ in 0..2 {
                let (found, evaluations) = nearest_on_curves(curve_a, curve_b, at);
                dec.charge_work(evaluations)?;
                let round = [
                    reentered(curve_a, closed_a, found[0]),
                    reentered(curve_b, closed_b, found[1]),
                ];
                at = [round[0].unwrap_or(found[0]), round[1].unwrap_or(found[1])];
                if round.iter().all(Option::is_none) {
                    break;
                }
            }
            // The two points the distance is read at.
            dec.charge_work(2)?;
            if on_edge(curve_a, at[0]) && on_edge(curve_b, at[1]) {
                value = value.min(distance(&curve_a.point(at[0]), &curve_b.point(at[1])));
            }
            if best.is_none_or(|d| value < d) {
                best = Some(value);
            }
        }
    }
    Ok(best)
}

fn on_edge(curve: &FramedCurve, t: f64) -> bool {
    (curve.range[0]..=curve.range[1]).contains(&t)
}

/// A parameter beyond an end of a closed edge, taken round to the other
/// end: the edge's two ends are one point, so the curve goes on there.
/// `None` for a parameter that needs no such move.
fn reentered(curve: &FramedCurve, closed: bool, t: f64) -> Option<f64> {
    let [low, high] = curve.range;
    if closed && t < low {
        Some(t + (high - low))
    } else if closed && t > high {
        Some(t - (high - low))
    } else {
        None
    }
}

/// An entity's tolerance, found by its place: the
/// report lists vertices, then edges, in definition order.
fn tolerance_at(tolerances: &[EntityTolerance], index: usize, kind: EntityKind, key: &str) -> f64 {
    tolerances
        .get(index)
        .filter(|entity| entity.kind == kind && entity.key == key)
        .map_or(0.0, |entity| entity.tolerance_mm.to_f64_nearest())
}

/// The u/v box of a face: its loops' u/v ends and the controls of its
/// retained u/v curves, which hold those curves.
fn face_rect(
    definition: &NumericalDefinition,
    face: usize,
    native: &Index<'_>,
    dec: &mut Decider,
) -> Result<[[f64; 2]; 2], Refusal> {
    let mut rect = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
    let mut take = |uv: [f64; 2]| {
        for axis in 0..2 {
            rect[axis][0] = rect[axis][0].min(uv[axis]);
            rect[axis][1] = rect[axis][1].max(uv[axis]);
        }
    };
    for trim in definition.faces[face].loops.iter().flatten() {
        dec.charge_work(2)?;
        for uv in &trim.uv {
            take([uv[0].to_f64_nearest(), uv[1].to_f64_nearest()]);
        }
    }
    for curve in native.curves_for_face(face) {
        let curve = FramedUvCurve::new(curve, dec)?;
        for control in curve.controls() {
            take(control);
        }
    }
    Ok(rect)
}

/// The parameters of the nodes of a curved copy from `start` to `end`, in
/// the order of the walk: [`LOOP_POINTS`] evenly spaced points from the
/// start, every place between them where the copy's u or v turns, and the
/// end.
///
/// A turn shows where a coordinate that rose from one point to the next
/// falls at a later step, or the reverse, with no step between that moves it
/// at all: two points either side of a turn can round to the same value. It
/// is polished by Brent's method between the start of the first of those
/// two steps and the end of the second. Next to each end one more point is
/// read, a forward-difference step inside it (Dennis and Schnabel 1983,
/// section 5.4: the square root of the machine epsilon, here as a share of
/// the use), so a turn between an end and the first point after it shows as
/// well. The walk takes it that each coordinate turns at most once between
/// two of its points.
fn walk_parameters(
    walk: &mut Walk,
    surface: &FramedSurface,
    [start, end]: [f64; 2],
    dec: &mut Decider,
) -> Result<Vec<f64>, Refusal> {
    let span = end - start;
    let even = |k: usize| start + span * k as f64 / LOOP_POINTS as f64;
    let inset = f64::EPSILON.sqrt() * span;
    let mut read = vec![start, start + inset];
    read.extend((1..LOOP_POINTS).map(even));
    read.extend([end - inset, end]);
    dec.charge_work(read.len() as u64)?;
    let points: Vec<[f64; 2]> = read.iter().map(|t| walk.at(surface, *t)).collect();
    let mut nodes: Vec<f64> = (0..LOOP_POINTS).map(even).collect();
    nodes.push(end);
    for axis in 0..2 {
        // The last step that moved this coordinate: the point it began at,
        // and whether it rose.
        let mut last: Option<(usize, bool)> = None;
        for k in 1..read.len() {
            let step = points[k][axis] - points[k - 1][axis];
            if step == 0.0 || step.is_nan() {
                continue;
            }
            let rising = step > 0.0;
            if let Some((began, rose)) = last
                && rose != rising
            {
                let (from, to) = (read[began], read[k]);
                let sign = if rose { 1.0 } else { -1.0 };
                let (share, _, evaluations) =
                    brent_maximum(&mut |x| sign * walk.at(surface, from + x * (to - from))[axis]);
                dec.charge_work(evaluations)?;
                nodes.push(from + share * (to - from));
            }
            last = Some((k - 1, rising));
        }
    }
    nodes.sort_by(|a, b| {
        if span < 0.0 {
            b.total_cmp(a)
        } else {
            a.total_cmp(b)
        }
    });
    nodes.dedup();
    Ok(nodes)
}

/// The nodes of one curved copy from `start` to `end`, the face's walk
/// number `index`: each runs along the copy to the next, and the last one
/// to whatever node follows it in its loop.
fn walked(
    walk: &mut Walk,
    index: usize,
    surface: &FramedSurface,
    ends: [f64; 2],
    dec: &mut Decider,
) -> Result<Vec<LoopNode>, Refusal> {
    let parameters = walk_parameters(walk, surface, ends, dec)?;
    dec.charge_work(parameters.len() as u64)?;
    Ok(parameters
        .windows(2)
        .map(|run| LoopNode {
            uv: walk.at(surface, run[0]),
            along: Some((index, [run[0], run[1]])),
        })
        .collect())
}

/// A face's loops as chains of nodes in u/v, and the curved copies they
/// follow. A use that is straight in u/v is its start. A use with a
/// retained u/v curve, and a circle or spline on a plane, is walked
/// ([`walk_parameters`]): each of its nodes runs along the copy to the
/// next, the last one to the start of the next use.
fn face_loops(
    definition: &NumericalDefinition,
    face: usize,
    native: &Index<'_>,
    surface: &FramedSurface,
    centre: &[Rational; 3],
    dec: &mut Decider,
) -> Result<(Vec<Vec<LoopNode>>, Vec<Walk>), Refusal> {
    let source = &definition.faces[face];
    let mut loops = Vec::new();
    let mut walks = Vec::new();
    for (boundary, row) in source.loops.iter().enumerate() {
        let mut nodes = Vec::new();
        for (use_index, trim) in row.iter().enumerate() {
            dec.charge_work(1)?;
            let edge = &definition.edges[trim.edge];
            let [a, b] = edge.parameters.each_ref().map(Rational::to_f64_nearest);
            let [start, end] = if trim.sense == Sense::Same {
                [a, b]
            } else {
                [b, a]
            };
            let walk = match native.get(face, boundary, use_index) {
                _ if end - start == 0.0 => None,
                Some(curve) => Some(Walk::Retained(FramedUvCurve::new(curve, dec)?)),
                None if analytic::projected_edge(source, edge) => Some(Walk::Projected(
                    FramedCurve::new(&edge.curve, edge.parameters.each_ref(), centre, dec)?,
                )),
                None => None,
            };
            let Some(mut walk) = walk else {
                nodes.push(LoopNode {
                    uv: [
                        trim.uv[0][0].to_f64_nearest(),
                        trim.uv[0][1].to_f64_nearest(),
                    ],
                    along: None,
                });
                continue;
            };
            nodes.extend(walked(&mut walk, walks.len(), surface, [start, end], dec)?);
            walks.push(walk);
        }
        loops.push(nodes);
    }
    Ok((loops, walks))
}

struct Prepared {
    entities: Vec<Entity>,
    curves: Vec<FramedCurve>,
    surfaces: Vec<Option<FramedSurface>>,
    rects: Vec<[[f64; 2]; 2]>,
}

fn entities(
    definition: &NumericalDefinition,
    native: &Index<'_>,
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Prepared, Refusal> {
    let mut out = Vec::new();
    let mut curves = Vec::new();
    for (index, vertex) in definition.vertices.iter().enumerate() {
        dec.charge_work(1)?;
        let p: [f64; 3] =
            std::array::from_fn(|axis| (&vertex.position[axis] - &centre[axis]).to_f64_nearest());
        out.push(Entity {
            kind: EntityKind::Vertex,
            key: vertex.key.clone(),
            tolerance: tolerance_at(tolerances, index, EntityKind::Vertex, &vertex.key),
            shape: Shape::Point(p),
            vertices: vec![index],
            bounds: [[p[0]; 2], [p[1]; 2], [p[2]; 2]],
        });
    }
    let first_edge = definition.vertices.len();
    for (index, edge) in definition.edges.iter().enumerate() {
        dec.charge_work(1)?;
        let mut curve = FramedCurve::new(&edge.curve, edge.parameters.each_ref(), centre, dec)?;
        dec.charge_work(2)?;
        let bounds = curve.bounding_box();
        let mut vertices = edge.vertices.to_vec();
        vertices.sort_unstable();
        vertices.dedup();
        out.push(Entity {
            kind: EntityKind::Edge,
            key: edge.key.clone(),
            tolerance: tolerance_at(tolerances, first_edge + index, EntityKind::Edge, &edge.key),
            shape: Shape::Curve(index),
            vertices,
            bounds,
        });
        curves.push(curve);
    }
    let mut surfaces = Vec::new();
    let mut rects = Vec::new();
    for (index, face) in definition.faces.iter().enumerate() {
        dec.charge_work(1)?;
        let rect = face_rect(definition, index, native, dec)?;
        let surface = FramedSurface::new(&face.surface, centre, dec)?;
        let mut vertices = Vec::new();
        let mut bounds: Option<Bounds> = surface.bounding_box(rect);
        for trim in face.loops.iter().flatten() {
            dec.charge_work(1)?;
            let edge = definition.edges.get(trim.edge).ok_or_else(|| {
                Refusal::new(
                    "numerical.invalid-definition",
                    "face trim names an unknown edge",
                )
            })?;
            vertices.extend(edge.vertices);
            if matches!(face.surface, catallax_shapes::Surface::Plane(_)) {
                let edge_box = out[first_edge + trim.edge].bounds;
                bounds = Some(bounds.map_or(edge_box, |b| union(&b, &edge_box)));
            }
        }
        vertices.sort_unstable();
        vertices.dedup();
        out.push(Entity {
            kind: EntityKind::Face,
            key: face.key.clone(),
            tolerance: 0.0,
            shape: Shape::Face(index),
            vertices,
            bounds: bounds.unwrap_or([[f64::NEG_INFINITY, f64::INFINITY]; 3]),
        });
        surfaces.push(Some(surface));
        rects.push(rect);
    }
    Ok(Prepared {
        entities: out,
        curves,
        surfaces,
        rects,
    })
}

/// Whether two entities share topology: one vertex is the other, or an end
/// of it, or on the boundary of it; two edges or an edge and a face have a
/// vertex in common.
fn adjacent(a: &Entity, b: &Entity) -> bool {
    a.vertices
        .iter()
        .any(|v| b.vertices.binary_search(v).is_ok())
}

/// Every pair of entities that share no topology and are closer than their
/// pair's reach, with the distance found.
fn close_pairs(
    definition: &NumericalDefinition,
    native: &Index<'_>,
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    reach_of: impl Fn(f64) -> f64,
    pair_reach: impl Fn(f64, f64) -> f64,
    dec: &mut Decider,
) -> Result<(Vec<Entity>, Vec<ClosePair>), Refusal> {
    let Prepared {
        entities,
        mut curves,
        mut surfaces,
        rects,
    } = entities(definition, native, tolerances, centre, dec)?;
    if entities.len() < 2 {
        return Ok((entities, Vec::new()));
    }
    let boxes: Vec<Bounds> = entities
        .iter()
        .map(|entity| widened(&entity.bounds, reach_of(entity.tolerance)))
        .collect();
    let mut items: Vec<usize> = (0..entities.len()).collect();
    let tree = Node::build(&mut items, &boxes, dec)?;
    let mut candidates = Vec::new();
    tree.self_pairs(&boxes, &mut candidates, dec)?;
    candidates.sort_unstable();
    candidates.dedup();
    let mut sampled: Vec<Option<Vec<Piece>>> = (0..curves.len()).map(|_| None).collect();
    let mut faces: Vec<Option<FaceData>> = (0..surfaces.len()).map(|_| None).collect();
    let mut found = Vec::new();
    for (i, j) in candidates {
        let (a, b) = (&entities[i], &entities[j]);
        dec.charge_work(1)?;
        if matches!((&a.shape, &b.shape), (Shape::Face(_), Shape::Face(_))) || adjacent(a, b) {
            continue;
        }
        let reach = pair_reach(a.tolerance, b.tolerance);
        for shape in [&a.shape, &b.shape] {
            match shape {
                Shape::Curve(k) if sampled[*k].is_none() => {
                    sampled[*k] = Some(pieces(&mut curves[*k], dec)?);
                }
                Shape::Face(k) if faces[*k].is_none() => {
                    let surface = surfaces[*k].take().expect("prepared");
                    let (loops, walks) = face_loops(definition, *k, native, &surface, centre, dec)?;
                    faces[*k] = Some(FaceData {
                        surface,
                        rect: rects[*k],
                        loops,
                        walks,
                    });
                }
                _ => {}
            }
        }
        let distance = match (&a.shape, &b.shape) {
            (Shape::Point(p), Shape::Point(q)) => {
                dec.charge_work(1)?;
                Some(distance(p, q))
            }
            (Shape::Point(p), Shape::Curve(k)) | (Shape::Curve(k), Shape::Point(p)) => {
                let pieces = sampled[*k].as_deref().expect("sampled above");
                point_to_curve(p, &mut curves[*k], pieces, reach, dec)?.map(|f| f.0)
            }
            (Shape::Curve(ka), Shape::Curve(kb)) => {
                // Edges come after vertices in index order, so `ka < kb`.
                let (low, high) = curves.split_at_mut(*kb);
                // An edge whose two ends are one vertex is closed.
                curve_to_curve(
                    (
                        &mut low[*ka],
                        sampled[*ka].as_deref().expect("sampled above"),
                        a.vertices.len() == 1,
                    ),
                    (
                        &mut high[0],
                        sampled[*kb].as_deref().expect("sampled above"),
                        b.vertices.len() == 1,
                    ),
                    reach,
                    dec,
                )?
            }
            (Shape::Point(p), Shape::Face(k)) | (Shape::Face(k), Shape::Point(p)) => {
                let face = faces[*k].as_mut().expect("prepared above");
                point_to_face(p, face, dec)?
            }
            (Shape::Curve(c), Shape::Face(k)) | (Shape::Face(k), Shape::Curve(c)) => {
                let face_box = entities[if matches!(a.shape, Shape::Face(_)) {
                    i
                } else {
                    j
                }]
                .bounds;
                let face = faces[*k].as_mut().expect("prepared above");
                let pieces = sampled[*c].as_deref().expect("sampled above");
                curve_to_face(&mut curves[*c], pieces, face, &face_box, reach, dec)?
            }
            (Shape::Face(_), Shape::Face(_)) => unreachable!("skipped above"),
        };
        if let Some(d) = distance.filter(|d| *d < reach) {
            found.push((i, j, d));
        }
    }
    Ok((entities, found))
}

fn exact(value: f64) -> Rational {
    Rational::from_f64_exact(value).unwrap_or_else(Rational::zero)
}

fn feature(entities: &[Entity], i: usize, j: usize, d: f64) -> FeaturePair {
    let (a, b) = (&entities[i], &entities[j]);
    FeaturePair {
        a: (a.kind, a.key.clone()),
        b: (b.kind, b.key.clone()),
        distance_mm: exact(d),
        tolerance_sum_mm: exact(a.tolerance + b.tolerance),
    }
}

/// Rule T5 on an admitted definition: refuses with
/// `numerical.tolerance-pair` the first pair whose tolerances add up to more
/// than half their distance, and returns the pairs over a fifth as warnings.
pub(in crate::numerical) fn pair_rule_check(
    definition: &NumericalDefinition,
    native: &Index<'_>,
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Vec<FeaturePair>, Refusal> {
    let (entities, rows) = close_pairs(
        definition,
        native,
        tolerances,
        centre,
        |t| PAIR_INFLATION * t,
        |a, b| PAIR_INFLATION * (a + b),
        dec,
    )?;
    let mut warnings = Vec::new();
    for (i, j, d) in &rows {
        let pair = feature(&entities, *i, *j, *d);
        match pair_verdict(entities[*i].tolerance, entities[*j].tolerance, *d) {
            PairVerdict::Clear => {}
            PairVerdict::Warn => warnings.push(pair),
            PairVerdict::Error => {
                return Err(Refusal::new(
                    "numerical.tolerance-pair",
                    format!(
                        "{} {:?} and {} {:?} share no topology and are {} mm apart, but their tolerances add up to {} mm, more than half that (rule T5)",
                        pair.a.0.name(),
                        pair.a.1,
                        pair.b.0.name(),
                        pair.b.1,
                        d,
                        entities[*i].tolerance + entities[*j].tolerance
                    ),
                ));
            }
        }
    }
    Ok(warnings)
}

/// Every disjoint pair of vertices, edges and face interiors of `body` within
/// `reach_mm`, excluding face-to-face pairs: the small features of rules T2
/// (`L`) and T8.
/// A short edge shows as its two end vertices; a sliver face as two of its
/// edges.
pub fn features_within(
    body: &NumericalBody,
    reach_mm: &Rational,
    dec: &mut Decider,
) -> Result<Vec<FeaturePair>, Refusal> {
    let reach = reach_mm.to_f64_nearest();
    let report = body.tolerances();
    let frame = super::BodyFrame::of(body.definition(), dec)?;
    let sampled_carriers = BTreeSet::new();
    let native = Index::new(
        body.definition(),
        body.native_trims(),
        &sampled_carriers,
        dec,
    )?;
    let (entities, rows) = close_pairs(
        body.definition(),
        &native,
        &report.entities,
        &frame.centre,
        |_| reach / 2.0,
        |_, _| reach,
        dec,
    )?;
    Ok(rows
        .iter()
        .map(|(i, j, d)| feature(&entities, *i, *j, *d))
        .collect())
}

/// Rule T2's fit target for a provisional embedding. Only entity distances
/// matter here. The producer can tighten the fit before T5 checks the
/// final body's stored tolerances.
pub fn fitted_target_for_definition(
    definition: &NumericalDefinition,
    native_trims: &[NumericalNativeTrim],
    eps: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let reach = fitted_feature_reach_mm(eps).to_f64_nearest();
    let frame = super::BodyFrame::of(definition, dec)?;
    let sampled_carriers = BTreeSet::new();
    let native = Index::new(definition, native_trims, &sampled_carriers, dec)?;
    let (_, pairs) = close_pairs(
        definition,
        &native,
        &[],
        &frame.centre,
        |_| reach / 2.0,
        |_, _| reach,
        dec,
    )?;
    let smallest = pairs.into_iter().map(|(_, _, d)| d).reduce(f64::min);
    Ok(own_fitted(eps, smallest.map(exact).as_ref()).target_mm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Circle, Curve, Plane, RPoint3, RVector3, Surface};

    fn q(text: &str) -> Rational {
        Rational::parse(text).expect("number")
    }

    fn origin() -> [Rational; 3] {
        [q("0"), q("0"), q("0")]
    }

    /// The plane z = 0, with u along x and v along y.
    fn plane(dec: &mut Decider) -> FramedSurface {
        let surface = Surface::Plane(Plane {
            origin: RPoint3::new(q("0"), q("0"), q("0")).to_real(),
            normal: RVector3::new(q("0"), q("0"), q("1")).to_real(),
            ref_dir: RVector3::new(q("1"), q("0"), q("0")).to_real(),
        });
        FramedSurface::new(&surface, &origin(), dec).expect("a plane")
    }

    /// The circle of radius 10 about (3, -2, 0) in that plane between two
    /// parameters, as the projected copy of an edge use.
    fn arc(ends: [f64; 2], dec: &mut Decider) -> Walk {
        let curve = Curve::Circle(Circle {
            center: RPoint3::new(q("3"), q("-2"), q("0")).to_real(),
            normal: RVector3::new(q("0"), q("0"), q("1")).to_real(),
            radius_sq: q("100").into(),
            ref_dir: RVector3::new(q("1"), q("0"), q("0")).to_real(),
        });
        let ends = ends.map(|t| Rational::from_f64_exact(t).expect("finite"));
        Walk::Projected(
            FramedCurve::new(&curve, [&ends[0], &ends[1]], &origin(), dec).expect("a circle"),
        )
    }

    /// An arc from 1.5 to 3.2 radians: its v turns at pi/2, between the
    /// start and the first of the nine points, and its u turns at pi,
    /// between the last of them and the end. Both turns are nodes, in
    /// either direction of the walk.
    #[test]
    fn a_walk_has_a_node_at_a_turn_next_to_either_end() {
        let mut dec = Decider::default();
        let surface = plane(&mut dec);
        let turns = [std::f64::consts::FRAC_PI_2, std::f64::consts::PI];
        for ends in [[1.5, 3.2], [3.2, 1.5]] {
            let mut walk = arc(ends, &mut dec);
            let nodes = walk_parameters(&mut walk, &surface, ends, &mut dec).unwrap();
            assert_eq!(nodes.len(), LOOP_POINTS + 1 + turns.len(), "{nodes:?}");
            assert_eq!(nodes.first(), Some(&ends[0]));
            assert_eq!(nodes.last(), Some(&ends[1]));
            assert!(
                nodes
                    .windows(2)
                    .all(|run| (run[1] - run[0]) * (ends[1] - ends[0]) > 0.0),
                "{nodes:?}"
            );
            for turn in turns {
                assert!(
                    nodes.iter().any(|t| (t - turn).abs() < 1e-6),
                    "{turn}: {nodes:?}"
                );
            }
        }
    }

    /// A face bounded by one complete circle of radius 10 is read on the
    /// circle: all the way round, a point a billionth of the radius inside
    /// it is inside the face and one a billionth outside is not. A polygon
    /// through the walk's points is up to 0.6 mm inside the circle.
    #[test]
    fn a_disc_is_read_on_its_circle_and_not_on_chords() {
        let mut dec = Decider::default();
        let surface = plane(&mut dec);
        let ends = [0.0, std::f64::consts::TAU];
        let mut walk = arc(ends, &mut dec);
        let nodes = walked(&mut walk, 0, &surface, ends, &mut dec).unwrap();
        let mut face = FaceData {
            surface,
            rect: [[-7.0, 13.0], [-12.0, 8.0]],
            loops: vec![nodes],
            walks: vec![walk],
        };
        for k in 0..360 {
            let angle = 0.001 + f64::from(k) * std::f64::consts::TAU / 360.0;
            let at = |radius: f64| [3.0 + radius * angle.cos(), -2.0 + radius * angle.sin()];
            assert!(face.inside(at(10.0 - 1e-8), &mut dec).unwrap(), "{angle}");
            assert!(!face.inside(at(10.0 + 1e-8), &mut dec).unwrap(), "{angle}");
        }
    }
}
