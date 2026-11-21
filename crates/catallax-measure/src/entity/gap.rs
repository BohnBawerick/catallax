//! Exact minimum distance between two trimmed entities by strata.
//!
//! A closed face is the union of its open interior, its open boundary edges
//! and its vertices; an edge is its open interior and its vertices. The
//! minimum over two entities is attained on some pair of strata, at a pair of
//! points where the distance restricted to those strata is stationary. Each
//! stratum lies on an analytic carrier at a constant offset from a core (see
//! [`super::critical`]), so every stationary pair is a core stationary pair
//! moved by the offsets along a common normal. Every such pair admitted by
//! both closed trims is a witness; the least admitted value over all strata is
//! the minimum.
//!
//! Strata pairs whose box separation cannot beat the best witness are skipped.
//! When the best witness attains the carriers' own lower bound
//! `max(0, core distance - offsets)`, every other stratum is skipped. A stratum
//! pair without a closed form that cannot be skipped refuses
//! `measure.distance-uncovered`.
//!
//! One-parameter sets of stationary pairs share one value and are probed at
//! exact dyadic positions. A missed admitted member is harmless: the admitted
//! part of such a set is closed and, unless it is the whole set, ends on a
//! lower stratum, where the same value is stationary again.
//!
//! Between two whole solids the face pairs are visited in order of their
//! boxes' separation, and one best distance prunes every later pair (see
//! [`solids`]).

use super::critical::Core;
use super::critical::{CorePair, Critical, Dirs, Family, critical};
use super::{EdgeCurve, PairEntity, edge_curve, enclose, face_surface, root, sign, unit, zero};
use crate::{Accuracy, Enclosure, distance::PreparedFace};
use catallax_numbers::{Answer, Answered, Ball, Decider, MarkLog, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, Surface, Vector3, distance::Finite};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, VertexId};
use std::collections::BTreeMap;

mod stationary;

/// Dyadic probe levels along a slid family: at most 31 positions.
const SLIDE_LEVELS: u32 = 4;
/// Dyadic probe levels per quadrant of a turned family: at most 32 positions.
const TURN_LEVELS: u32 = 3;
/// Levels of each factor when a family both slides and turns.
const PRODUCT_LEVELS: u32 = 2;
const BOX_BITS: u32 = 64;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Face(FaceId),
    Edge(EdgeId),
    Vertex(VertexId),
}

/// Exact membership of a carrier point in a closed stratum.
enum Member {
    /// A closed surface face, a full circle edge or a vertex.
    Whole,
    Finite(Finite),
    Trim(FaceId),
    Segment(Point3, Point3),
    /// The counter-clockwise arc about `normal` from `from` to `to`.
    Arc {
        centre: Point3,
        normal: Vector3,
        from: Point3,
        to: Point3,
    },
}

struct Stratum {
    kind: Kind,
    /// `None` when the curve has no closed form here.
    core: Option<Core>,
    offset: Option<Real>,
    member: Member,
    extent: Option<[Ball; 3]>,
}

struct Scratch<'a> {
    model: &'a Model,
    copy: Option<Model>,
}

impl Scratch<'_> {
    fn model(&mut self) -> &mut Model {
        self.copy.get_or_insert_with(|| self.model.clone())
    }

    fn contains(&mut self, face: FaceId, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let model = self.copy.get_or_insert_with(|| self.model.clone());
        catallax_boolean::face_contains(model, face, p, dec)
    }
}

/// Half-plane class of `v` relative to the arc start `a`: 0 at angle zero,
/// 1 on the counter-clockwise side, 2 opposite, 3 on the clockwise side.
fn class(a: &Vector3, v: &Vector3, normal: &Vector3, dec: &mut Decider) -> Result<u8, Refusal> {
    Ok(match sign(&normal.dot(&a.cross(v)), dec)? {
        Sign::Positive => 1,
        Sign::Negative => 3,
        Sign::Zero => {
            if sign(&a.dot(v), dec)? == Sign::Positive {
                0
            } else {
                2
            }
        }
    })
}

impl Member {
    fn contains(
        &self,
        p: &Point3,
        scratch: &mut Scratch,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        match self {
            Member::Whole => Ok(true),
            Member::Finite(finite) => finite.contains_carrier_point(p, dec),
            Member::Trim(face) => scratch.contains(*face, p, dec),
            Member::Segment(a, b) => {
                let d = b.minus(a);
                let t = p.minus(a).dot(&d);
                Ok(sign(&t, dec)? != Sign::Negative
                    && sign(&(&t - &d.norm_squared()), dec)? != Sign::Positive)
            }
            Member::Arc {
                centre,
                normal,
                from,
                to,
            } => {
                let a = from.minus(centre);
                let x = p.minus(centre);
                let b = to.minus(centre);
                let (cx, cb) = (class(&a, &x, normal, dec)?, class(&a, &b, normal, dec)?);
                Ok(match cx.cmp(&cb) {
                    std::cmp::Ordering::Less => true,
                    std::cmp::Ordering::Greater => false,
                    std::cmp::Ordering::Equal => {
                        cx % 2 == 0 || sign(&normal.dot(&x.cross(&b)), dec)? != Sign::Negative
                    }
                })
            }
        }
    }
}

fn point_extent(p: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&p.x, BOX_BITS)?,
        dec.enclose_budgeted_at(&p.y, BOX_BITS)?,
        dec.enclose_budgeted_at(&p.z, BOX_BITS)?,
    ])
}

fn face_stratum(model: &Model, face: FaceId, dec: &mut Decider) -> Result<Stratum, Refusal> {
    let (core, offset) = match face_surface(model, face)? {
        Surface::Plane(s) => (
            Core::Plane {
                origin: s.origin.clone(),
                normal: unit(&s.normal),
            },
            None,
        ),
        Surface::Cylinder(s) => (
            Core::Line {
                origin: s.origin.clone(),
                dir: unit(&s.axis),
            },
            Some(root(&s.radius_sq)),
        ),
        Surface::Cone(s) => (
            Core::Cone {
                apex: s.apex.clone(),
                axis: unit(&s.axis),
                slope: s.slope.clone(),
            },
            None,
        ),
        Surface::Sphere(s) => (Core::Point(s.center.clone()), Some(s.radius.clone())),
        Surface::Torus(s) => (
            Core::Circle {
                centre: s.center.clone(),
                normal: unit(&s.axis),
                radius: s.major.clone(),
                reference: unit(&s.ref_dir),
            },
            Some(s.minor.clone()),
        ),
        _ => unreachable!("face_surface admits only analytic carriers"),
    };
    let member = if model.face(face).loops.is_empty() {
        Member::Whole
    } else if let Some(finite) = PreparedFace::recognize(model, face, dec)? {
        Member::Finite(finite)
    } else {
        Member::Trim(face)
    };
    let bounds = crate::face_bounding_box(
        model,
        face,
        &Accuracy {
            relative: Rational::zero(),
            absolute: Rational::from_ratio(1, 1_000_000).unwrap(),
        },
        dec,
    )?;
    let extent = bounds.marks.is_empty().then(|| {
        std::array::from_fn(|i| {
            Ball::from_rational(&bounds.value.min[i], BOX_BITS)
                .hull(&Ball::from_rational(&bounds.value.max[i], BOX_BITS))
        })
    });
    Ok(Stratum {
        kind: Kind::Face(face),
        core: Some(core),
        offset,
        member,
        extent,
    })
}

fn edge_stratum(model: &Model, edge: EdgeId, dec: &mut Decider) -> Result<Stratum, Refusal> {
    let unknown = Stratum {
        kind: Kind::Edge(edge),
        core: None,
        offset: None,
        member: Member::Whole,
        extent: None,
    };
    let Ok(curve) = edge_curve(model, edge, dec) else {
        return Ok(unknown);
    };
    Ok(match curve {
        EdgeCurve::Segment(a, b) => {
            let (ea, eb) = (point_extent(&a, dec)?, point_extent(&b, dec)?);
            Stratum {
                kind: Kind::Edge(edge),
                core: Some(Core::Line {
                    origin: a.clone(),
                    dir: unit(&b.minus(&a)),
                }),
                offset: None,
                member: Member::Segment(a, b),
                extent: Some(std::array::from_fn(|i| ea[i].hull(&eb[i]))),
            }
        }
        EdgeCurve::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } => {
            let normal = unit(&normal);
            let radius = root(&radius_sq);
            let e = model.edge(edge);
            let member = match (e.start, e.end) {
                (Some(s), Some(t)) if s != t => {
                    let at = |v: VertexId| model.geom().position(model.vertex(v).point);
                    let (from, to) = if e.sense == Sense::Same {
                        (at(s), at(t))
                    } else {
                        (at(t), at(s))
                    };
                    Member::Arc {
                        centre: center.clone(),
                        normal: normal.clone(),
                        from,
                        to,
                    }
                }
                _ => Member::Whole,
            };
            let centre = point_extent(&center, dec)?;
            let n = [&normal.x, &normal.y, &normal.z];
            let mut extent = centre;
            for i in 0..3 {
                let half = &radius * &(Real::one() - n[i].square()).sqrt();
                let half = dec.enclose_budgeted_at(&half, BOX_BITS)?;
                extent[i] = (&extent[i] - &half).hull(&(&extent[i] + &half));
            }
            Stratum {
                kind: Kind::Edge(edge),
                core: Some(Core::Circle {
                    centre: center.clone(),
                    normal,
                    radius,
                    reference: unit(&ref_dir),
                }),
                offset: None,
                member,
                extent: Some(extent),
            }
        }
    })
}

fn vertex_stratum(model: &Model, vertex: VertexId, dec: &mut Decider) -> Result<Stratum, Refusal> {
    let p = model.geom().position(model.vertex(vertex).point);
    Ok(Stratum {
        kind: Kind::Vertex(vertex),
        extent: Some(point_extent(&p, dec)?),
        core: Some(Core::Point(p)),
        offset: None,
        member: Member::Whole,
    })
}

/// The entity's own stratum first, then its edges and vertices.
fn strata(model: &Model, entity: PairEntity, dec: &mut Decider) -> Result<Vec<Stratum>, Refusal> {
    let mut out = Vec::new();
    let mut edges = Vec::new();
    let mut vertices = Vec::new();
    let add_edge = |e: EdgeId, edges: &mut Vec<EdgeId>, vertices: &mut Vec<VertexId>| {
        if !edges.contains(&e) {
            edges.push(e);
        }
        for v in [model.edge(e).start, model.edge(e).end]
            .into_iter()
            .flatten()
        {
            if !vertices.contains(&v) {
                vertices.push(v);
            }
        }
    };
    match entity {
        PairEntity::Face(face) => {
            out.push(face_stratum(model, face, dec)?);
            for lp in &model.face(face).loops {
                match model.loop_(*lp) {
                    Loop::Vertex(v) => {
                        if !vertices.contains(v) {
                            vertices.push(*v);
                        }
                    }
                    Loop::Coedges(cs) => {
                        for c in cs {
                            add_edge(model.coedge(*c).edge, &mut edges, &mut vertices);
                        }
                    }
                }
            }
        }
        PairEntity::Edge(edge) => {
            edge_curve(model, edge, dec)?;
            add_edge(edge, &mut edges, &mut vertices);
        }
    }
    for e in edges {
        out.push(edge_stratum(model, e, dec)?);
    }
    for v in vertices {
        out.push(vertex_stratum(model, v, dec)?);
    }
    Ok(out)
}

fn uncovered(detail: String) -> Refusal {
    Refusal::new("measure.distance-uncovered", detail)
}

/// A lower bound on the distance between two boxes.
fn separation(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    let mut sum = Rational::zero();
    for i in 0..3 {
        let gap = (a[i].lo() - b[i].hi())
            .max(b[i].lo() - a[i].hi())
            .max(Rational::zero());
        sum = sum + &gap * &gap;
    }
    Ball::from_rational(&sum, BOX_BITS)
        .sqrt()
        .map_or(Rational::zero(), |r| r.lo().max(Rational::zero()))
}

/// One way to offset a core pair onto both carriers.
#[derive(Clone)]
struct Combo {
    /// Explicit unit direction for meeting cores; `None` derives it from the pair.
    u: Option<Vector3>,
    sx: i64,
    sy: i64,
    value: Real,
    lo: Rational,
}

struct Best {
    value: Real,
    hi: Rational,
    lo: Rational,
}

struct Search<'s, 'a> {
    scratch: &'s mut Scratch<'a>,
    best: Option<Best>,
    /// An upper bound already attained by another pair of entities. A pair
    /// that cannot come below it cannot lower the minimum.
    cutoff: Option<Rational>,
}

fn signs(offset: &Option<Real>) -> &'static [i64] {
    if offset.is_some() { &[1, -1] } else { &[0] }
}

fn scaled_offset(offset: &Option<Real>, s: i64) -> Real {
    offset
        .as_ref()
        .map_or_else(Real::zero, |r| r * &Real::from_i64(s))
}

/// Rotates `v` by the exact angle `(c, s)` about the unit `axis`.
fn rotate(v: &Vector3, axis: &Vector3, c: &Rational, s: &Rational) -> Vector3 {
    let (c, s) = (Real::from(c.clone()), Real::from(s.clone()));
    let along = axis.scaled(&(&axis.dot(v) * &(Real::one() - c.clone())));
    v.scaled(&c).plus(&axis.cross(v).scaled(&s)).plus(&along)
}

/// Exact rational points of the unit circle, coarsest first, all quadrants.
fn turns(levels: u32) -> Vec<(Rational, Rational)> {
    let mut ts = vec![Rational::zero()];
    for level in 1..=levels {
        let count = 1i64 << level;
        for k in (1..count).step_by(2) {
            ts.push(Rational::from_ratio(k, count).unwrap());
        }
    }
    let mut out = Vec::new();
    for t in ts {
        let d = Rational::one() + &t * &t;
        let c = &(Rational::one() - &t * &t) / &d;
        let s = &(&t * &Rational::from_i64(2)) / &d;
        out.push((c.clone(), s.clone()));
        out.push((-&s, c.clone()));
        out.push((-&c, -&s));
        out.push((s, -c));
    }
    out
}

/// Dyadic midpoints of `[lo, hi]`, coarsest first.
fn slides(lo: &Rational, hi: &Rational, levels: u32) -> Vec<Rational> {
    let mut out = Vec::new();
    for level in 0..=levels {
        let count = 1i64 << level;
        for m in 0..count {
            out.push(lo + &(&(hi - lo) * &Rational::from_ratio(2 * m + 1, 2 * count).unwrap()));
        }
    }
    out
}

/// Parameters `t` with `p + t d` possibly inside `extent`, conservatively.
fn slide_range(
    p: &Point3,
    d: &Vector3,
    extent: &[Ball; 3],
    dec: &mut Decider,
) -> Result<Option<(Rational, Rational)>, Refusal> {
    let p = point_extent(p, dec)?;
    let d = [
        dec.enclose_budgeted_at(&d.x, BOX_BITS)?,
        dec.enclose_budgeted_at(&d.y, BOX_BITS)?,
        dec.enclose_budgeted_at(&d.z, BOX_BITS)?,
    ];
    let mut range: Option<(Rational, Rational)> = None;
    for i in 0..3 {
        if d[i].contains_zero() {
            continue;
        }
        let a = (&Ball::from_rational(&extent[i].lo(), BOX_BITS) - &p[i]).checked_div(&d[i]);
        let b = (&Ball::from_rational(&extent[i].hi(), BOX_BITS) - &p[i]).checked_div(&d[i]);
        let (Some(a), Some(b)) = (a, b) else { continue };
        let lo = a.lo().min(b.lo());
        let hi = a.hi().max(b.hi());
        range = Some(match range {
            None => (lo, hi),
            Some((l, h)) => (l.max(lo), h.min(hi)),
        });
    }
    Ok(range.filter(|(lo, hi)| lo <= hi))
}

impl Search<'_, '_> {
    /// The least upper bound attained so far, here or by another entity pair.
    fn bound(&self) -> Option<&Rational> {
        match (self.best.as_ref().map(|b| &b.hi), self.cutoff.as_ref()) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        }
    }

    /// Visits every pair of strata and returns the index pairs without a
    /// closed form.
    fn visit(
        &mut self,
        sa: &[Stratum],
        sb: &[Stratum],
        dec: &mut Decider,
    ) -> Result<Vec<(usize, usize)>, Refusal> {
        let mut deferred = Vec::new();
        'strata: for (i, x) in sa.iter().enumerate() {
            for (j, y) in sb.iter().enumerate() {
                dec.check_deadline()?;
                let top = i == 0 && j == 0;
                if !top && skip(self, x, y) {
                    continue;
                }
                let found = match (&x.core, &y.core) {
                    (Some(cx), Some(cy)) => critical(cx, cy, dec)?,
                    _ => None,
                };
                let Some(found) = found else {
                    if !parallel_faces(x, y, dec)? {
                        deferred.push((i, j));
                    }
                    continue;
                };
                if !self.strata_pair(x, y, &found, dec)? {
                    deferred.push((i, j));
                    continue;
                }
                // Every stratum lies on the carriers, so attaining their lower
                // bound ends the search.
                if top
                    && let Some(best) = &self.best
                    && zero(&(&best.value - &floor(x, y, &found, dec)?), dec)?
                {
                    break 'strata;
                }
            }
        }
        Ok(deferred)
    }

    /// The ways to offset `pair` onto both carriers, or `None` when the common
    /// normals form a two-dimensional set.
    fn combos(
        &self,
        x: &Stratum,
        y: &Stratum,
        pair: &CorePair,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Combo>>, Refusal> {
        let make = |u, sx, sy, dec: &mut Decider| -> Result<Combo, Refusal> {
            let v = &(&pair.dist - &scaled_offset(&x.offset, sx)) - &scaled_offset(&y.offset, sy);
            let value = if sign(&v, dec)? == Sign::Negative {
                -v
            } else {
                v
            };
            let lo = dec.enclose_budgeted_at(&value, BOX_BITS)?.lo();
            Ok(Combo {
                u,
                sx,
                sy,
                value,
                lo,
            })
        };
        if x.offset.is_none() && y.offset.is_none() {
            return Ok(Some(vec![make(None, 0, 0, dec)?]));
        }
        let us: Vec<Option<Vector3>> = match &pair.dirs {
            Dirs::Along => vec![None],
            Dirs::Fixed(list) => list.iter().cloned().map(Some).collect(),
            Dirs::All => return Ok(None),
        };
        let mut out = Vec::new();
        for u in us {
            for &sx in signs(&x.offset) {
                for &sy in signs(&y.offset) {
                    out.push(make(u.clone(), sx, sy, dec)?);
                }
            }
        }
        out.sort_by(|a, b| a.lo.cmp(&b.lo));
        Ok(Some(out))
    }

    fn may_improve(&self, combo: &Combo) -> bool {
        self.bound().is_none_or(|b| combo.lo < *b)
    }

    /// Offers the offset points of `pair`; true when both trims admit them.
    fn offer(
        &mut self,
        x: &Stratum,
        y: &Stratum,
        pair: &CorePair,
        combo: &Combo,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let u = match &combo.u {
            Some(u) => u.clone(),
            None if x.offset.is_none() && y.offset.is_none() => Vector3::zero(),
            None => pair
                .b
                .minus(&pair.a)
                .scaled(&(Real::one() / pair.dist.clone())),
        };
        let p = pair.a.plus(&u.scaled(&scaled_offset(&x.offset, combo.sx)));
        let q = pair.b.plus(&u.scaled(&scaled_offset(&y.offset, -combo.sy)));
        if !x.member.contains(&p, self.scratch, dec)?
            || !y.member.contains(&q, self.scratch, dec)?
        {
            return Ok(false);
        }
        self.consider(&combo.value, dec)?;
        Ok(true)
    }

    fn consider(&mut self, value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
        let ball = dec.enclose_budgeted_at(value, 128)?;
        let better = match &self.best {
            None => true,
            Some(best) if ball.hi() < best.lo => true,
            Some(best) if ball.lo() > best.hi => false,
            Some(best) => sign(&(value - &best.value), dec)? == Sign::Negative,
        };
        if better {
            self.best = Some(Best {
                value: value.clone(),
                lo: ball.lo(),
                hi: ball.hi(),
            });
        }
        Ok(())
    }

    /// Offers every stationary pair of two strata. False when no closed form exists.
    fn strata_pair(
        &mut self,
        x: &Stratum,
        y: &Stratum,
        critical_pairs: &Critical,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for pair in &critical_pairs.pairs {
            let Some(combos) = self.combos(x, y, pair, dec)? else {
                return Ok(false);
            };
            for combo in combos {
                if self.may_improve(&combo) {
                    self.offer(x, y, pair, &combo, dec)?;
                }
            }
        }
        for family in &critical_pairs.families {
            if !self.probe(x, y, family, dec)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn probe(
        &mut self,
        x: &Stratum,
        y: &Stratum,
        family: &Family,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Some(mut combos) = self.combos(x, y, &family.pair, dec)? else {
            return Ok(false);
        };
        combos.retain(|c| self.may_improve(c));
        if combos.is_empty() {
            return Ok(true);
        }
        let levels = |own| {
            if family.slide.is_some() && family.turn.is_some() {
                PRODUCT_LEVELS
            } else {
                own
            }
        };
        let offsets = match &family.slide {
            None => vec![None],
            Some(d) => {
                let (Some(ex), Some(ey)) = (&x.extent, &y.extent) else {
                    return Ok(true);
                };
                let (Some(rx), Some(ry)) = (
                    slide_range(&family.pair.a, d, ex, dec)?,
                    slide_range(&family.pair.b, d, ey, dec)?,
                ) else {
                    return Ok(true);
                };
                let (lo, hi) = (rx.0.max(ry.0), rx.1.min(ry.1));
                if lo > hi {
                    return Ok(true);
                }
                slides(&lo, &hi, levels(SLIDE_LEVELS))
                    .into_iter()
                    .map(Some)
                    .collect()
            }
        };
        let angles = match &family.turn {
            None => vec![None],
            Some(_) => turns(levels(TURN_LEVELS)).into_iter().map(Some).collect(),
        };
        for t in &offsets {
            for angle in &angles {
                dec.check_deadline()?;
                let mut pair = family.pair.clone();
                let mut combos_here = combos.clone();
                if let (Some(t), Some(d)) = (t, &family.slide) {
                    let shift = d.scaled(&Real::from(t.clone()));
                    pair.a = pair.a.plus(&shift);
                    pair.b = pair.b.plus(&shift);
                }
                if let (Some((c, s)), Some((centre, axis))) = (angle, &family.turn) {
                    let turn = |p: &Point3| centre.plus(&rotate(&p.minus(centre), axis, c, s));
                    pair.a = turn(&pair.a);
                    pair.b = turn(&pair.b);
                    for combo in &mut combos_here {
                        if let Some(u) = &combo.u {
                            combo.u = Some(rotate(u, axis, c, s));
                        }
                    }
                }
                let mut admitted = Vec::new();
                for (i, combo) in combos_here.iter().enumerate() {
                    if self.may_improve(combo) && self.offer(x, y, &pair, combo, dec)? {
                        admitted.push(i);
                    }
                }
                let mut i = 0;
                combos.retain(|_| {
                    let keep = !admitted.contains(&i);
                    i += 1;
                    keep
                });
                combos.retain(|c| self.may_improve(c));
                if combos.is_empty() {
                    return Ok(true);
                }
            }
        }
        Ok(true)
    }
}

/// The carriers' lower bound `max(0, least core distance - offsets)`, or zero
/// when either core is a cone.
fn floor(
    x: &Stratum,
    y: &Stratum,
    critical_pairs: &Critical,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let (Some(cx), Some(cy)) = (&x.core, &y.core) else {
        return Ok(Real::zero());
    };
    if cx.is_cone() || cy.is_cone() {
        return Ok(Real::zero());
    }
    let dists = critical_pairs
        .pairs
        .iter()
        .chain(critical_pairs.families.iter().map(|f| &f.pair))
        .map(|p| p.dist.clone());
    let mut least: Option<(Real, Rational)> = None;
    for d in dists {
        let lo = dec.enclose_budgeted_at(&d, BOX_BITS)?.lo();
        least = Some(match least {
            None => (d, lo),
            Some((m, mlo)) => {
                if sign(&(&d - &m), dec)? == Sign::Negative {
                    (d, lo)
                } else {
                    (m, mlo)
                }
            }
        });
    }
    let Some((least, _)) = least else {
        return Ok(Real::zero());
    };
    let value = &(&least - &scaled_offset(&x.offset, 1)) - &scaled_offset(&y.offset, 1);
    Ok(if sign(&value, dec)? == Sign::Positive {
        value
    } else {
        Real::zero()
    })
}

/// Parallel planar faces keep the planar witness search. Faces without exact
/// rational frames, or without a witness, are left to the strata search, where
/// their edges and vertices attain the same minimum (see [`parallel_faces`]).
/// With `general`, a witness search that finds nothing first falls back to the
/// general trimmed search, as single entity pairs always have.
fn parallel_planes(
    model: &Model,
    a: FaceId,
    b: FaceId,
    general: bool,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<Option<Answered<Enclosure>>, Refusal> {
    let (Surface::Plane(p), Surface::Plane(q)) = (face_surface(model, a)?, face_surface(model, b)?)
    else {
        return Ok(None);
    };
    if !zero(&p.normal.cross(&q.normal).norm_squared(), dec)? {
        return Ok(None);
    }
    if [&p.origin, &q.origin]
        .iter()
        .any(|o| o.as_rational().is_none())
        || [&p.normal, &q.normal]
            .iter()
            .any(|n| n.as_rational().is_none())
    {
        return Ok(None);
    }
    let width = if accuracy.absolute.sign() == Sign::Positive {
        accuracy.absolute.clone()
    } else {
        Rational::parse("1e-30").unwrap()
    };
    let options = catallax_shapes::distance::Options {
        width,
        ..Default::default()
    };
    let mut copy = model.clone();
    let pa = PreparedFace::new(&copy, a, dec)?;
    let pb = PreparedFace::new(&copy, b, dec)?;
    let bounds = if general {
        pa.parallel_face_to_face(&pb, &mut copy, &options, dec)?
    } else {
        match pa.parallel_witness(&pb, &mut copy, &options, dec)? {
            Some(bounds) => bounds,
            None => return Ok(None),
        }
    };
    let Some(upper) = bounds.upper else {
        return Ok(None);
    };
    let enclosure = Enclosure {
        lo: bounds.lower,
        hi: upper,
    };
    let mut marks = MarkLog::new();
    accuracy.mark(&mut marks, &[&enclosure], "minimum distance");
    Ok(Some(Ok(Answer::with_marks(enclosure, marks))))
}

pub(super) fn minimum(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let sa = strata(model, a, dec)?;
    let sb = strata(model, b, dec)?;
    // Entities sharing an edge or vertex touch there.
    if sa.iter().any(|x| sb.iter().any(|y| x.kind == y.kind)) {
        return Ok(Answer::exact(Enclosure::exact(Rational::zero())));
    }
    if let (PairEntity::Face(fa), PairEntity::Face(fb)) = (a, b)
        && let Some(answer) = parallel_planes(model, fa, fb, true, accuracy, dec)?
    {
        return answer;
    }
    let mut scratch = Scratch { model, copy: None };
    let mut search = Search {
        scratch: &mut scratch,
        best: None,
        cutoff: None,
    };
    let deferred = search.visit(&sa, &sb, dec)?;
    if let Some((x, y)) = deferred
        .into_iter()
        .map(|(i, j)| (&sa[i], &sb[j]))
        .find(|(x, y)| !skip(&search, x, y))
    {
        return Err(uncovered(format!(
            "no closed form between {:?} and {:?} on {a:?} and {b:?}",
            x.kind, y.kind
        )));
    }
    let Some(best) = search.best else {
        return Err(uncovered(format!(
            "no stationary pair of {a:?} and {b:?} lies on both"
        )));
    };
    enclose(&best.value, accuracy, "minimum distance", dec)
}

/// Whether the strata's boxes are at least as far apart as the best witness.
fn skip(search: &Search, x: &Stratum, y: &Stratum) -> bool {
    let (Some(bound), Some(ex), Some(ey)) = (search.bound(), &x.extent, &y.extent) else {
        return false;
    };
    separation(ex, ey) >= *bound
}

/// Whether `x` and `y` are two parallel planar faces. Their stationary pairs
/// form a two-parameter family with no closed form, but none is needed: the
/// pairs at the plane gap form a closed region of the overlap of the two
/// projected faces, and a boundary point of that region lies on an edge or a
/// vertex of one face. The same minimum is then stationary on that lower
/// stratum, which the search visits with the faces' own edges and vertices.
fn parallel_faces(x: &Stratum, y: &Stratum, dec: &mut Decider) -> Result<bool, Refusal> {
    match (x.kind, y.kind, &x.core, &y.core) {
        (
            Kind::Face(_),
            Kind::Face(_),
            Some(Core::Plane { normal: n1, .. }),
            Some(Core::Plane { normal: n2, .. }),
        ) => zero(&n1.cross(n2).norm_squared(), dec),
        _ => Ok(false),
    }
}

/// One face's strata, or the refusal of an unsupported carrier, and its box.
struct Prepared {
    face: FaceId,
    strata: Result<Vec<Stratum>, Refusal>,
    extent: Option<[Ball; 3]>,
}

fn prepare(model: &Model, face: FaceId, dec: &mut Decider) -> Result<Prepared, Refusal> {
    match strata(model, PairEntity::Face(face), dec) {
        Ok(strata) => Ok(Prepared {
            face,
            extent: strata[0].extent.clone(),
            strata: Ok(strata),
        }),
        Err(refusal) if refusal.code == "measure.unsupported-surface" => {
            // A carrier without closed forms still has a proven box, which is
            // enough to set the face aside when it is far from the minimum.
            let bounds = crate::face_bounding_box(
                model,
                face,
                &Accuracy {
                    relative: Rational::zero(),
                    absolute: Rational::from_ratio(1, 1_000_000).unwrap(),
                },
                dec,
            );
            dec.check_deadline()?;
            let extent = bounds.ok().filter(|b| b.marks.is_empty()).map(|b| {
                std::array::from_fn(|i| {
                    Ball::from_rational(&b.value.min[i], BOX_BITS)
                        .hull(&Ball::from_rational(&b.value.max[i], BOX_BITS))
                })
            });
            Ok(Prepared {
                face,
                strata: Err(refusal),
                extent,
            })
        }
        Err(refusal) => Err(refusal),
    }
}

fn lower(x: &Option<[Ball; 3]>, y: &Option<[Ball; 3]>) -> Rational {
    match (x, y) {
        (Some(x), Some(y)) => separation(x, y),
        _ => Rational::zero(),
    }
}

/// An open pair: strata without a closed form, or a face without a supported
/// carrier, whose boxes' separation is `lower`.
struct Open {
    lower: Rational,
    faces: (usize, usize),
    /// Indices into the two faces' strata, when both carriers are supported.
    strata: Option<(usize, usize)>,
    what: String,
}

/// Minimum distance between the closed boundaries of two bodies whose faces do
/// not cross, as `BRepExtrema_DistShapeShape` in OpenCASCADE prunes it: face
/// pairs are visited in order of their boxes' separation, and a pair whose
/// boxes are at least as far apart as the best distance found is skipped. The
/// strata search of each visited pair also skips anything that cannot come
/// below that distance. Strata pairs without a closed form whose boxes are
/// nearer than the final minimum are then settled face pair by face pair with
/// the certified interval search, which proves the two closed faces no nearer
/// than that minimum ([`PreparedFace::face_to_face_below`]). Any pair left
/// unsettled refuses, stating the enclosure reached.
pub(super) fn solids(
    model: &Model,
    a: BodyId,
    b: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let start = catallax_numbers::current_operation_work().used;
    let mut sides = [Vec::new(), Vec::new()];
    for (side, body) in sides.iter_mut().zip([a, b]) {
        for face in model.faces_of_body(body) {
            side.push(prepare(model, face, dec)?);
        }
    }
    let [pa, pb] = sides;
    let mut order = Vec::with_capacity(pa.len() * pb.len());
    for (i, x) in pa.iter().enumerate() {
        for (j, y) in pb.iter().enumerate() {
            order.push((lower(&x.extent, &y.extent), i, j));
        }
    }
    order.sort();
    let mut scratch = Scratch { model, copy: None };
    let mut found: Option<Enclosure> = None;
    let mut marks = MarkLog::new();
    let mut open: Vec<Open> = Vec::new();
    for (separated, i, j) in order {
        dec.check_deadline()?;
        let cutoff = found.as_ref().map(|f| f.hi.clone());
        if cutoff.as_ref().is_some_and(|c| separated >= *c) {
            break;
        }
        let (x, y) = (&pa[i], &pb[j]);
        let (Ok(sx), Ok(sy)) = (&x.strata, &y.strata) else {
            open.push(Open {
                lower: separated,
                faces: (i, j),
                strata: None,
                what: format!("{:?} and {:?}, unsupported carriers", x.face, y.face),
            });
            continue;
        };
        let answer = if sx.iter().any(|s| sy.iter().any(|t| s.kind == t.kind)) {
            Some(Answer::exact(Enclosure::exact(Rational::zero())))
        } else if let Some(answer) = parallel_planes(model, x.face, y.face, false, accuracy, dec)? {
            Some(answer?)
        } else {
            let mut search = Search {
                scratch: &mut scratch,
                best: None,
                cutoff,
            };
            for (si, ti) in search.visit(sx, sy, dec)? {
                let (s, t) = (&sx[si], &sy[ti]);
                open.push(Open {
                    lower: lower(&s.extent, &t.extent),
                    faces: (i, j),
                    strata: Some((si, ti)),
                    what: format!(
                        "{:?} and {:?} on {:?} and {:?}",
                        s.kind, t.kind, x.face, y.face
                    ),
                });
            }
            match search.best {
                Some(best) => Some(enclose(&best.value, accuracy, "minimum distance", dec)?),
                None => None,
            }
        };
        if let Some(answer) = answer {
            marks.extend(&answer.marks);
            found = Some(match found {
                None => answer.value,
                Some(f) => Enclosure {
                    lo: f.lo.min(answer.value.lo),
                    hi: f.hi.min(answer.value.hi),
                },
            });
        }
    }
    // Every open pair lies in the closures of its two faces, so a proof that
    // those faces are no nearer than the minimum settles all of its pairs.
    // Without any witness yet, a box diagonal of both solids bounds every
    // distance between them from above and serves as the ceiling instead.
    let ceiling = match &found {
        Some(f) => f.hi.clone(),
        None => (0..3).fold(Rational::zero(), |sum, k| {
            let hull = pa
                .iter()
                .chain(&pb)
                .filter_map(|p| p.extent.as_ref())
                .map(|e| (e[k].lo(), e[k].hi()))
                .reduce(|(l, h), (a, b)| (l.min(a), h.max(b)));
            sum + hull.map_or_else(Rational::zero, |(l, h)| h - l)
        }),
    };
    let mut pairs: Vec<(Rational, usize, usize)> = Vec::new();
    for o in open.iter().filter(|o| o.lower < ceiling) {
        match pairs.iter_mut().find(|p| (p.1, p.2) == o.faces) {
            Some(p) => p.0 = p.0.clone().min(o.lower.clone()),
            None => pairs.push((o.lower.clone(), o.faces.0, o.faces.1)),
        }
    }
    pairs.sort();
    // With a witness found, settling only tightens a bound the closed forms
    // already proved, so it may spend at most as much work again as the search
    // that found them: a pair whose faces make every certified cell costly,
    // such as planes bounded by spline edges, then keeps its box separation
    // instead of exhausting the request. Without a witness, settling is the
    // only route to an answer and runs within the request's own allowance, as
    // it does outside one. One cell budget, the certified search's own default
    // for one face pair, serves every open face pair in order.
    let spent = catallax_numbers::current_operation_work().used - start;
    let mut settle = || {
        let mut cells = catallax_shapes::distance::Options::default().max_cells;
        let mut prepared: BTreeMap<FaceId, Option<PreparedFace>> = BTreeMap::new();
        let mut unsettled: Vec<(Rational, usize, usize, String)> = Vec::new();
        let mut witness: Option<Rational> = None;
        let mut stopped = false;
        for (least, i, j) in &pairs {
            let (fa, fb) = (pa[*i].face, pb[*j].face);
            let what = open
                .iter()
                .find(|o| o.faces == (*i, *j))
                .map(|o| o.what.clone())
                .unwrap_or_default();
            if stopped || cells == 0 {
                unsettled.push((least.clone(), *i, *j, what));
                continue;
            }
            for face in [fa, fb] {
                if let std::collections::btree_map::Entry::Vacant(slot) = prepared.entry(face) {
                    slot.insert(PreparedFace::new(model, face, dec).ok());
                }
            }
            let (Some(Some(x)), Some(Some(y))) = (prepared.get(&fa), prepared.get(&fb)) else {
                stopped = dec.check_deadline().is_err();
                unsettled.push((least.clone(), *i, *j, what));
                continue;
            };
            match x.face_to_face_below(y, scratch.model(), &ceiling, cells, dec) {
                Ok(bounds) => {
                    cells = cells.saturating_sub(bounds.cells);
                    if bounds.witnesses.is_some()
                        && let Some(upper) = bounds.upper
                    {
                        witness = Some(witness.map_or(upper.clone(), |w| w.min(upper)));
                    }
                    if bounds.lower < ceiling {
                        unsettled.push((bounds.lower.max(least.clone()), *i, *j, what));
                    }
                }
                Err(_) => {
                    stopped = dec.check_deadline().is_err();
                    unsettled.push((least.clone(), *i, *j, what));
                }
            }
        }
        (unsettled, witness)
    };
    let (unsettled, witness) = if pairs.is_empty() {
        (Vec::new(), None)
    } else if found.is_none() || !catallax_numbers::operation_budget_active() {
        settle()
    } else if let Ok(allowance) = catallax_numbers::OperationAllowance::reserve(spent) {
        allowance.run(settle).0
    } else {
        let unsettled = pairs
            .iter()
            .map(|(least, i, j)| {
                let what = open
                    .iter()
                    .find(|o| o.faces == (*i, *j))
                    .map(|o| o.what.clone())
                    .unwrap_or_default();
                (least.clone(), *i, *j, what)
            })
            .collect();
        (unsettled, None)
    };
    dec.check_deadline()?;
    if !unsettled.is_empty() {
        return unsettled_answer(
            Unsettled {
                pa: &pa,
                pb: &pb,
                open: &open,
                unsettled: &unsettled,
                found,
                witness,
                ceiling,
                marks,
            },
            &mut scratch,
            (a, b),
            accuracy,
            dec,
        );
    }
    let Some(found) = found else {
        return Err(uncovered(format!(
            "no stationary pair of {a:?} and {b:?} lies on both boundaries"
        )));
    };
    Ok(Answer::with_marks(found, marks))
}

/// What the closed forms and the face-level settling left for
/// [`unsettled_answer`].
struct Unsettled<'u> {
    pa: &'u [Prepared],
    pb: &'u [Prepared],
    open: &'u [Open],
    /// Face pairs the certified interval search did not settle, with the lower
    /// bound it proved for each and the first open strata pair's description.
    unsettled: &'u [(Rational, usize, usize, String)],
    found: Option<Enclosure>,
    witness: Option<Rational>,
    ceiling: Rational,
    marks: MarkLog,
}

/// Settles the open strata pairs of unsettled face pairs by their certified
/// stationary pairs ([`stationary`]), then states the enclosure. A face pair
/// whose every open strata pair is settled no longer widens the answer with
/// the bounded mark; what stays open keeps the lower bound proved for it. The
/// pass runs as an attempt whose refusal returns its work: the answer is then
/// the one the face-level search gave.
fn unsettled_answer(
    u: Unsettled,
    scratch: &mut Scratch,
    (a, b): (BodyId, BodyId),
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let Unsettled {
        pa,
        pb,
        open,
        unsettled,
        found,
        witness,
        ceiling,
        mut marks,
    } = u;
    let before = found
        .as_ref()
        .map(|f| f.hi.clone())
        .into_iter()
        .chain(witness.clone())
        .min()
        .unwrap_or(ceiling);
    // The open strata pairs of every unsettled face pair, nearest first.
    let mut todo: Vec<usize> = (0..open.len())
        .filter(|&k| unsettled.iter().any(|u| (u.1, u.2) == open[k].faces))
        .collect();
    todo.sort_by(|&p, &q| open[p].lower.cmp(&open[q].lower).then(p.cmp(&q)));
    let certified = catallax_numbers::refund_on_refusal(|| {
        // One budget, the certified search's own default for one face pair,
        // serves every open strata pair in order.
        let cells = catallax_shapes::distance::Options::default().max_cells;
        let mut pass = stationary::Pass::new(scratch, Some(before), cells);
        let mut settled = BTreeMap::new();
        for &k in &todo {
            let o = &open[k];
            let (Some((s, t)), Ok(sx), Ok(sy)) =
                (o.strata, &pa[o.faces.0].strata, &pb[o.faces.1].strata)
            else {
                continue;
            };
            if let Some(result) = pass.settle(&sx[s], &sy[t], accuracy, dec)? {
                settled.insert(k, result);
            }
        }
        Ok::<_, Refusal>((settled, pass.witness, pass.slack))
    });
    let (settled, extra, slack) =
        certified.unwrap_or_else(|_| (BTreeMap::new(), None, Rational::zero()));
    // Each face pair's lower bound over its open strata pairs, and whether
    // every one of them is settled.
    let mut parts: Vec<(Option<Rational>, bool, &str)> = Vec::new();
    for (proved, i, j, what) in unsettled {
        let mut least: Option<Rational> = None;
        let mut settled_all = true;
        for (k, o) in open.iter().enumerate().filter(|(_, o)| o.faces == (*i, *j)) {
            let value = match settled.get(&k) {
                Some(r) => {
                    settled_all &= !r.exhausted;
                    r.low.as_ref().map(|l| l.clone().max(o.lower.clone()))
                }
                None => {
                    settled_all = false;
                    Some(o.lower.clone())
                }
            };
            if let Some(v) = value {
                least = Some(least.map_or(v.clone(), |l| l.min(v)));
            }
        }
        parts.push((
            least.map(|l| l.max(proved.clone())),
            settled_all,
            what.as_str(),
        ));
    }
    let hi = found
        .as_ref()
        .map(|f| f.hi.clone())
        .into_iter()
        .chain(witness)
        .chain(extra)
        .min();
    let lo = found
        .as_ref()
        .map(|f| f.lo.clone())
        .into_iter()
        .chain(parts.iter().filter_map(|p| p.0.clone()))
        .min();
    // The first open face pair with the least lower bound names the mark.
    let open_part = parts
        .iter()
        .filter(|p| !p.1)
        .min_by(|p, q| match (&p.0, &q.0) {
            (Some(x), Some(y)) => x.cmp(y),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
    let Some(hi) = hi else {
        return Err(match open_part {
            Some(p) => uncovered(format!(
                "no closed form between {} and no witness pair; minimum gap is at least {} mm",
                p.2,
                lo.unwrap_or_else(Rational::zero).to_decimal_bounds(6).0
            )),
            None => uncovered(format!(
                "no stationary pair of {a:?} and {b:?} lies on both boundaries"
            )),
        });
    };
    // Without any value below the bound in force, nothing lies below the
    // least upper bound attained. Dropped solutions widen the lower end by
    // their summed diagonals (see [`stationary`]).
    let lo = lo
        .map_or_else(|| hi.clone(), |l| (&l - &slack).max(Rational::zero()))
        .min(hi.clone());
    let enclosure = Enclosure { lo, hi };
    match open_part {
        // A proven enclosure, wider than closed forms would give: the bounded
        // answer kind, marked when it misses the requested accuracy.
        Some(p) => accuracy.mark(
            &mut marks,
            &[&enclosure],
            &format!("minimum distance, with no closed form between {},", p.2),
        ),
        None => accuracy.mark(&mut marks, &[&enclosure], "minimum distance"),
    }
    Ok(Answer::with_marks(enclosure, marks))
}
