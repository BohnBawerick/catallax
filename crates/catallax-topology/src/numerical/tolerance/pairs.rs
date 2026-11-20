//! Rule T5's pair rule, and the small features rules T2 and T8 read.
//!
//! Two entities that share no topology are a pair: two different vertices, a
//! vertex and an edge that does not end at it, and two edges with no vertex in
//! common. Their boxes go in a bounding-box tree, each widened by its reach,
//! and only the pairs whose widened boxes overlap are measured. On a part
//! with no small features the tree returns no pair and nothing is measured.
//!
//! A face of a numerical body stores no tolerance (its surface is the
//! authoritative copy), so a pair with a face has `t_face = 0` and is read
//! through the face's edges and vertices. See the milestone note for what
//! that leaves out.
//!
//! Distances are read in binary64 in the body's frame (rule T1), by the
//! sampling rule of the edge-gap checker (section 4.2 step 4) and Brent's
//! method. They are searches for the least distance, as the checker's are
//! searches for the largest gap: the found value is at least the true one.

use catallax_numbers::{Decider, Rational, Refusal};

use super::{EntityKind, EntityTolerance, PAIR_INFLATION, PairVerdict, pair_verdict};
use crate::numerical::gap::{FramedCurve, brent_maximum, chebyshev_shares};
use crate::numerical::{NumericalBody, NumericalDefinition};

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
}

struct Entity {
    kind: EntityKind,
    key: String,
    tolerance: f64,
    shape: Shape,
    /// Vertex indices this entity is or ends at.
    vertices: [usize; 2],
    bounds: Bounds,
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

fn entities(
    definition: &NumericalDefinition,
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    dec: &mut Decider,
) -> Result<(Vec<Entity>, Vec<FramedCurve>), Refusal> {
    let mut out = Vec::new();
    let mut curves = Vec::new();
    let tolerance_of = |kind: EntityKind, key: &str| {
        tolerances
            .iter()
            .find(|entity| entity.kind == kind && entity.key == key)
            .map_or(0.0, |entity| entity.tolerance_mm.to_f64_nearest())
    };
    for (index, vertex) in definition.vertices.iter().enumerate() {
        dec.charge_work(1)?;
        let p: [f64; 3] =
            std::array::from_fn(|axis| (&vertex.position[axis] - &centre[axis]).to_f64_nearest());
        out.push(Entity {
            kind: EntityKind::Vertex,
            key: vertex.key.clone(),
            tolerance: tolerance_of(EntityKind::Vertex, &vertex.key),
            shape: Shape::Point(p),
            vertices: [index, index],
            bounds: [[p[0]; 2], [p[1]; 2], [p[2]; 2]],
        });
    }
    for (index, edge) in definition.edges.iter().enumerate() {
        dec.charge_work(1)?;
        let mut curve = FramedCurve::new(&edge.curve, edge.parameters.each_ref(), centre, dec)?;
        let bounds = curve.bounding_box();
        out.push(Entity {
            kind: EntityKind::Edge,
            key: edge.key.clone(),
            tolerance: tolerance_of(EntityKind::Edge, &edge.key),
            shape: Shape::Curve(index),
            vertices: edge.vertices,
            bounds,
        });
        curves.push(curve);
    }
    Ok((out, curves))
}

/// Whether two entities share topology: one vertex is the other, or an end
/// of it, or the two edges have an end in common.
fn adjacent(a: &Entity, b: &Entity) -> bool {
    a.vertices.iter().any(|v| b.vertices.contains(v))
}

/// Every pair of entities that share no topology and are closer than their
/// pair's reach, with the distance found.
fn close_pairs(
    definition: &NumericalDefinition,
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    reach_of: impl Fn(f64) -> f64,
    pair_reach: impl Fn(f64, f64) -> f64,
    dec: &mut Decider,
) -> Result<(Vec<Entity>, Vec<ClosePair>), Refusal> {
    let (entities, mut curves) = entities(definition, tolerances, centre, dec)?;
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
    let mut found = Vec::new();
    for (i, j) in candidates {
        let (a, b) = (&entities[i], &entities[j]);
        if adjacent(a, b) {
            continue;
        }
        let reach = pair_reach(a.tolerance, b.tolerance);
        for shape in [&a.shape, &b.shape] {
            if let Shape::Curve(k) = shape
                && sampled[*k].is_none()
            {
                sampled[*k] = Some(pieces(&mut curves[*k], dec)?);
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
    tolerances: &[EntityTolerance],
    centre: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Vec<FeaturePair>, Refusal> {
    let (entities, rows) = close_pairs(
        definition,
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

/// Every pair of vertices and edges of `body` that share no topology and are
/// closer than `reach_mm`: the small features of rules T2 (`L`) and T8.
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
    let (entities, rows) = close_pairs(
        body.definition(),
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
