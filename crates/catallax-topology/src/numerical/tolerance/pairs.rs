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
//! Distances are read in binary64 in the body's frame (rule T1), by the
//! sampling rule of the edge-gap checker (section 4.2 step 4) and Brent's
//! method. They are searches for the least distance, as the checker's are
//! searches for the largest gap: the found value is at least the true one.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Rational, Refusal};

use super::{
    EntityKind, EntityTolerance, PAIR_INFLATION, PairVerdict, fitted_feature_reach_mm, own_fitted,
    pair_verdict,
};
use crate::Sense;
use crate::numerical::gap::{
    FramedCurve, FramedSurface, FramedUvCurve, brent_maximum, chebyshev_shares, try_brent_maximum,
};
use crate::numerical::{NumericalBody, NumericalDefinition, NumericalNativeTrim, native::Index};

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
/// each edge use on a numerical carrier. A retained u/v curve is read at
/// these points when a face's loops are written as polygons.
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

/// A face as the pair rule reads it: its surface, its u/v box and its loops
/// as u/v polygons, made when a pair first needs them.
struct FaceData {
    surface: FramedSurface,
    rect: [[f64; 2]; 2],
    polygons: Vec<Vec<[f64; 2]>>,
}

impl FaceData {
    /// Whether `uv` lies inside the face's loops, by the even-odd rule, with
    /// u shifted by a period on a surface of revolution.
    fn inside(&self, uv: [f64; 2], dec: &mut Decider) -> Result<bool, Refusal> {
        let shifts: &[f64] = match self.surface.period_u {
            Some(p) => &[0.0, -p, p],
            None => &[0.0],
        };
        for shift in shifts {
            let (u, v) = (uv[0] + shift, uv[1]);
            let mut inside = false;
            for polygon in &self.polygons {
                dec.charge_work(polygon.len() as u64)?;
                for (k, a) in polygon.iter().enumerate() {
                    let b = polygon[(k + 1) % polygon.len()];
                    if (a[1] > v) != (b[1] > v)
                        && u < a[0] + (v - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                    {
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

/// The distance from `point` to the interior of a face: the foot of the
/// point on the surface, seeded by the nearest point of a grid on the face's
/// u/v box, where that foot lies inside the face. `None` where it does not:
/// the face's boundary is read through its edges and vertices. The foot is
/// returned either way, to seed the next point along a curve.
fn point_to_face(
    point: &[f64; 3],
    face: &mut FaceData,
    seed: Option<[f64; 2]>,
    dec: &mut Decider,
) -> Result<(Option<f64>, [f64; 2]), Refusal> {
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
    let (foot, evaluations) = face.surface.foot(point, seed);
    dec.charge_work(evaluations + 1)?;
    if !foot.iter().all(|x| x.is_finite()) {
        return Ok((None, seed));
    }
    if !face.inside(foot, dec)? {
        return Ok((None, foot));
    }
    Ok((Some(distance(&face.surface.point(foot), point)), foot))
}

/// The least distance from a curve's pieces within `reach` of a face's box
/// to the face's interior: every sample's foot, then a Brent polish around
/// the nearest.
fn curve_to_face(
    curve: &mut FramedCurve,
    pieces: &[Piece],
    face: &mut FaceData,
    face_box: &Bounds,
    reach: f64,
    dec: &mut Decider,
) -> Result<Option<f64>, Refusal> {
    let mut best: Option<(f64, f64, [f64; 2], [f64; 2])> = None;
    let mut seed = None;
    for piece in pieces {
        if box_distance(&piece.bounds, face_box) > reach {
            continue;
        }
        for (k, (t, p)) in piece.samples.iter().enumerate() {
            // A sample farther from the face's box than the reach cannot
            // come within it; no foot is sought for it.
            dec.charge_work(1)?;
            if box_distance(&[[p[0]; 2], [p[1]; 2], [p[2]; 2]], face_box) > reach {
                continue;
            }
            let (found, foot) = point_to_face(p, face, seed, dec)?;
            seed = Some(foot);
            let Some(d) = found else {
                continue;
            };
            let bracket = [
                piece.samples[k.saturating_sub(1)].0,
                piece.samples[(k + 1).min(piece.samples.len() - 1)].0,
            ];
            if best.is_none_or(|(b, ..)| d < b) {
                best = Some((d, *t, foot, bracket));
            }
        }
    }
    let Some((mut value, _, foot, [a, b])) = best else {
        return Ok(None);
    };
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
    value = value.min(-negative);
    Ok(Some(value))
}

/// The least distance between two curves within `reach`: the nearest sample
/// pair of every pair of pieces whose boxes are within reach, polished by
/// alternating Brent searches on each curve.
fn curve_to_curve(
    a: (&mut FramedCurve, &[Piece]),
    b: (&mut FramedCurve, &[Piece]),
    reach: f64,
    dec: &mut Decider,
) -> Result<Option<f64>, Refusal> {
    let (curve_a, pieces_a) = a;
    let (curve_b, pieces_b) = b;
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
            if best.is_none_or(|d| value < d) {
                best = Some(value);
            }
        }
    }
    Ok(best)
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

/// A face's loops as u/v polygons: each use's start, and, along a retained
/// u/v curve, the curve at evenly spaced parameters.
fn face_polygons(
    definition: &NumericalDefinition,
    face: usize,
    native: &Index<'_>,
    dec: &mut Decider,
) -> Result<Vec<Vec<[f64; 2]>>, Refusal> {
    let mut polygons = Vec::new();
    for (boundary, row) in definition.faces[face].loops.iter().enumerate() {
        let mut polygon = Vec::new();
        for (use_index, trim) in row.iter().enumerate() {
            dec.charge_work(1)?;
            match native.get(face, boundary, use_index) {
                Some(curve) => {
                    let edge = &definition.edges[trim.edge];
                    let [a, b] = edge.parameters.each_ref().map(Rational::to_f64_nearest);
                    let [start, end] = if trim.sense == Sense::Same {
                        [a, b]
                    } else {
                        [b, a]
                    };
                    let mut curve = FramedUvCurve::new(curve, dec)?;
                    dec.charge_work(LOOP_POINTS as u64)?;
                    for k in 0..LOOP_POINTS {
                        polygon.push(
                            curve.point(start + (end - start) * k as f64 / LOOP_POINTS as f64),
                        );
                    }
                }
                None => polygon.push([
                    trim.uv[0][0].to_f64_nearest(),
                    trim.uv[0][1].to_f64_nearest(),
                ]),
            }
        }
        polygons.push(polygon);
    }
    Ok(polygons)
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
                    faces[*k] = Some(FaceData {
                        surface: surfaces[*k].take().expect("prepared"),
                        rect: rects[*k],
                        polygons: face_polygons(definition, *k, native, dec)?,
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
                curve_to_curve(
                    (
                        &mut low[*ka],
                        sampled[*ka].as_deref().expect("sampled above"),
                    ),
                    (
                        &mut high[0],
                        sampled[*kb].as_deref().expect("sampled above"),
                    ),
                    reach,
                    dec,
                )?
            }
            (Shape::Point(p), Shape::Face(k)) | (Shape::Face(k), Shape::Point(p)) => {
                let face = faces[*k].as_mut().expect("prepared above");
                point_to_face(p, face, None, dec)?.0
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
