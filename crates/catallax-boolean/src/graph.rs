//! Boundary intersections trimmed to both operands, with canonical exact
//! points and curves. Temporary carrier cuts are removed from the answer.

use crate::geometry::{canonical_point, coincident, on_arc, refusal};
use crate::points::curve_sample;
use crate::region::{Location, RegionCache};
use crate::stages::{IntersectionGraph, SectionEdge};
use catallax_intersect::{
    Component, CurveSurface, Cyclic, SurfacePair, curve_is_closed, cyclic_order,
    intersect_curve_surface, intersect_surfaces, on_curve, order_on_curve, same_curve,
};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, MarkLog, Outcome, Refusal, StrictOrder,
    ZeroProof,
};
use catallax_shapes::{CurveId, GeomKey, PointId, SurfaceId};
use catallax_topology::{BodyId, BodyKind, FaceId, Loop, Model};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) fn sort_points(
    model: &Model,
    curve: CurveId,
    points: &mut [PointId],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if points.len() < 2 {
        return Ok(());
    }
    let open = !curve_is_closed(model.geom(), curve, dec)?;
    let anchor = points[0];
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 {
            let before = if open {
                match order_on_curve(model.geom(), curve, points[j], points[j - 1], dec) {
                    Outcome::Proven(StrictOrder::Less) => true,
                    Outcome::Proven(StrictOrder::Greater) => false,
                    Outcome::ProvenCoincident(_) => {
                        return Err(refusal(
                            "boolean.imprint-ambiguous",
                            "canonical curve cuts were not distinct",
                        ));
                    }
                    Outcome::NotExact(n) => {
                        return Err(refusal("boolean.imprint-ambiguous", n.detail));
                    }
                    Outcome::Refused(e) => return Err(e),
                }
            } else if j == 1 {
                false
            } else {
                match cyclic_order(model.geom(), curve, [anchor, points[j], points[j - 1]], dec) {
                    Outcome::Proven(Cyclic::Forward) => true,
                    Outcome::Proven(Cyclic::Backward) => false,
                    Outcome::ProvenCoincident(_) => {
                        return Err(refusal(
                            "boolean.imprint-ambiguous",
                            "cyclic cuts were not distinct",
                        ));
                    }
                    Outcome::NotExact(n) => {
                        return Err(refusal("boolean.imprint-ambiguous", n.detail));
                    }
                    Outcome::Refused(e) => return Err(e),
                }
            };
            if !before {
                break;
            }
            points.swap(j, j - 1);
            j -= 1;
        }
    }
    Ok(())
}

fn face_edges(model: &Model, face: FaceId) -> Vec<CurveId> {
    let mut result = Vec::new();
    for &l in &model.face(face).loops {
        if let Loop::Coedges(cs) = model.loop_(l) {
            for &c in cs {
                let curve = model.edge(model.coedge(c).edge).curve;
                if !result.contains(&curve) {
                    result.push(curve);
                }
            }
        }
    }
    result
}

fn add_curve(
    model: &mut Model,
    pool: &mut Vec<CurveId>,
    curve: CurveId,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    for (index, &other) in pool.iter().enumerate() {
        if coincident(
            same_curve(model.geom(), curve, other, dec),
            "boolean.imprint-ambiguous",
        )? {
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Coincident,
                a: GeomKey::Curve(curve),
                b: GeomKey::Curve(other),
                at: None,
                origin: FactOrigin::Proof(ZeroProof::Algebraic),
            });
            return Ok(index);
        }
    }
    let index = pool.len();
    pool.push(curve);
    Ok(index)
}

/// Canonicalizes the hits of `curve` on every operand surface in query order.
/// Returns the raw hits, their marks and the surfaces whose whole meeting with
/// `curve` is now canonical.
fn carrier_hits(
    model: &mut Model,
    curve: CurveId,
    surfaces: &[SurfaceId],
    points: &mut Vec<PointId>,
    dec: &mut Decider,
) -> Result<(Vec<PointId>, MarkLog, BTreeSet<SurfaceId>), Refusal> {
    let mut found = Vec::new();
    let mut marks = MarkLog::new();
    let mut complete = BTreeSet::new();
    for &surface in surfaces {
        let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
        if !matches!(answer.value, CurveSurface::Contained) {
            complete.insert(surface);
        }
        if let CurveSurface::Meets(hits) = answer.value {
            for hit in hits {
                if !model.geom().facts().holds(
                    FactKind::Incident,
                    GeomKey::Point(hit.point),
                    GeomKey::Curve(curve),
                ) {
                    return Err(refusal(
                        "boolean.imprint-ambiguous",
                        "intersection branch membership lacks an incidence proof",
                    ));
                }
                canonical_point(model, points, hit.point, dec, "boolean.imprint-ambiguous")?;
                found.push(hit.point);
            }
        }
        marks.extend(&answer.marks);
    }
    Ok((found, marks, complete))
}

/// Borrowed graph state shared by every carrier's pieces.
struct Context<'a> {
    first: &'a [FaceId],
    second: &'a [FaceId],
    mandatory: &'a BTreeSet<PointId>,
    point_faces: &'a BTreeMap<PointId, SharedFaces>,
}

/// The section pieces of `curve` between consecutive sorted `cuts`.
fn carrier_pieces(
    model: &mut Model,
    curve: CurveId,
    cuts: &[PointId],
    closed: bool,
    context: &Context<'_>,
    regions: &mut RegionCache,
    dec: &mut Decider,
) -> Result<Vec<SectionEdge>, Refusal> {
    let arcs = if !closed {
        cuts.windows(2)
            .map(|w| (Some(w[0]), Some(w[1])))
            .collect::<Vec<_>>()
    } else if cuts.is_empty() {
        vec![(None, None)]
    } else {
        (0..cuts.len())
            .map(|i| (Some(cuts[i]), Some(cuts[(i + 1) % cuts.len()])))
            .collect()
    };
    let mut pieces: Vec<SectionEdge> = Vec::new();
    for (start, end) in arcs {
        let sample = curve_sample(model, curve, start, end, dec)?;
        let Some((faces_of_first, faces_of_second)) =
            shared_faces(model, context.first, context.second, sample, regions, dec)?
        else {
            continue;
        };
        let edge = SectionEdge {
            curve,
            start,
            end,
            faces_of_first,
            faces_of_second,
        };
        if let Some(last) = pieces.last_mut()
            && last.end == edge.start
            && same_faces(last, &edge)
            && edge.start.is_some_and(|p| {
                !context.mandatory.contains(&p)
                    && context.point_faces.get(&p).is_some_and(|(first, second)| {
                        first == &edge.faces_of_first && second == &edge.faces_of_second
                    })
            })
        {
            last.end = edge.end;
            continue;
        }
        pieces.push(edge);
    }
    if closed && pieces.len() > 1 {
        let n = pieces.len();
        if pieces[n - 1].end == pieces[0].start
            && same_faces(&pieces[n - 1], &pieces[0])
            && pieces[0].start.is_some_and(|p| {
                !context.mandatory.contains(&p)
                    && context.point_faces.get(&p).is_some_and(|(first, second)| {
                        first == &pieces[0].faces_of_first && second == &pieces[0].faces_of_second
                    })
            })
        {
            let last = pieces.pop().expect("nonempty sections");
            pieces[0].start = last.start;
        }
    }
    if closed
        && pieces.len() == 1
        && pieces[0].start == pieces[0].end
        && pieces[0].start.is_some_and(|p| {
            !context.mandatory.contains(&p)
                && context.point_faces.get(&p).is_some_and(|(first, second)| {
                    first == &pieces[0].faces_of_first && second == &pieces[0].faces_of_second
                })
        })
    {
        pieces[0].start = None;
        pieces[0].end = None;
    }
    Ok(pieces)
}

/// Discovery order of the complete curve/surface pass, recovered after
/// skipping inactive carriers. The skipped carriers are queried only when
/// their hits could come first.
#[derive(Default)]
struct Order {
    /// The active carrier whose hits first found each point; seeds are 0.
    finder: BTreeMap<PointId, usize>,
    /// Skipped carriers that cross each point before its finder.
    crossers: BTreeMap<PointId, BTreeSet<usize>>,
    /// Position in the pruned discovery order.
    index: BTreeMap<PointId, usize>,
    /// Hits of each queried skipped carrier, in query order.
    hits: BTreeMap<usize, Vec<PointId>>,
}

impl Order {
    /// The complete pass's discovery key of `point`: its carrier and position.
    fn key(
        &mut self,
        model: &mut Model,
        pool: &[CurveId],
        surfaces: &[SurfaceId],
        point: PointId,
        dec: &mut Decider,
    ) -> Result<(usize, usize), Refusal> {
        for skipped in self.crossers.get(&point).cloned().unwrap_or_default() {
            if let std::collections::btree_map::Entry::Vacant(entry) = self.hits.entry(skipped) {
                let (found, ..) =
                    carrier_hits(model, pool[skipped], surfaces, &mut Vec::new(), dec)?;
                entry.insert(found);
            }
            for (rank, &hit) in self.hits[&skipped].iter().enumerate() {
                let mut known = vec![point];
                if canonical_point(model, &mut known, hit, dec, "boolean.imprint-ambiguous")?
                    == point
                {
                    return Ok((skipped, rank));
                }
            }
        }
        let finder = self.finder.get(&point).copied().unwrap_or(usize::MAX);
        Ok((finder, self.index[&point]))
    }

    /// The position in sorted `cuts` of the cut the complete pass finds first.
    fn first(
        &mut self,
        model: &mut Model,
        pool: &[CurveId],
        surfaces: &[SurfaceId],
        cuts: &[PointId],
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut candidates: Vec<_> = (0..cuts.len())
            .filter_map(|i| Some((*self.crossers.get(&cuts[i])?.first()?, i)))
            .collect();
        if candidates.is_empty() {
            return Ok(0);
        }
        candidates.sort_unstable();
        let mut best = (usize::MAX, usize::MAX);
        let mut start = 0;
        for (i, &cut) in cuts.iter().enumerate() {
            if !self.crossers.contains_key(&cut) {
                let key = (
                    self.finder.get(&cut).copied().unwrap_or(usize::MAX),
                    self.index[&cut],
                );
                if key < best {
                    (best, start) = (key, i);
                }
            }
        }
        // A skipped carrier's hit comes first only from an earlier carrier.
        for (lower, i) in candidates {
            if lower > best.0 {
                break;
            }
            let key = self.key(model, pool, surfaces, cuts[i], dec)?;
            if key < best {
                (best, start) = (key, i);
            }
        }
        Ok(start)
    }
}

pub(crate) fn containing_faces(
    model: &mut Model,
    faces: &[FaceId],
    point: PointId,
    regions: &mut RegionCache,
    dec: &mut Decider,
) -> Result<Vec<FaceId>, Refusal> {
    let mut result = Vec::new();
    for &face in faces {
        if coincident(
            model.geom().side(model.face(face).surface, point, dec),
            "boolean.imprint-ambiguous",
        )? && regions.locate(model, face, point, dec)? != Location::Outside
        {
            result.push(face);
        }
    }
    Ok(result)
}

type SharedFaces = (Vec<FaceId>, Vec<FaceId>);

// Every graph consumer requires membership on both operands. Exact exclusion
// from the smaller face set therefore settles the conjunction without testing
// the other set; incidence alone never discharges either finite trim test.
fn shared_faces(
    model: &mut Model,
    first: &[FaceId],
    second: &[FaceId],
    point: PointId,
    regions: &mut RegionCache,
    dec: &mut Decider,
) -> Result<Option<SharedFaces>, Refusal> {
    let reversed = first.len() > second.len();
    let (smaller, larger) = if reversed {
        (second, first)
    } else {
        (first, second)
    };
    let small_faces = containing_faces(model, smaller, point, regions, dec)?;
    if small_faces.is_empty() {
        return Ok(None);
    }
    let large_faces = containing_faces(model, larger, point, regions, dec)?;
    if large_faces.is_empty() {
        return Ok(None);
    }
    Ok(Some(if reversed {
        (large_faces, small_faces)
    } else {
        (small_faces, large_faces)
    }))
}

fn same_faces(a: &SectionEdge, b: &SectionEdge) -> bool {
    a.faces_of_first == b.faces_of_first && a.faces_of_second == b.faces_of_second
}

// A failed pruning probe must not change the carrier/point definitions chosen
// by the ordinary graph. The append-only model snapshot shares its history;
// all optional witnesses and support meetings stay private to this probe.
fn boundaries_apart(model: &Model, first: &[FaceId], second: &[FaceId], dec: &mut Decider) -> bool {
    let mut probe = model.clone();
    let mut bounds = crate::bounds::FaceBounds::new(false);
    let mut pairs = BTreeMap::new();
    for &a in first {
        for &b in second {
            if bounds.excludes(&mut probe, a, b, dec) {
                continue;
            }
            let surfaces = (probe.face(a).surface, probe.face(b).surface);
            let disjoint = *pairs.entry(surfaces).or_insert_with(|| {
                intersect_surfaces(probe.geom_mut(), surfaces.0, surfaces.1, dec).is_ok_and(
                    |answer| {
                        answer.marks.is_empty() && matches!(answer.value, SurfacePair::Disjoint)
                    },
                )
            });
            if !disjoint {
                return false;
            }
        }
    }
    true
}

pub(crate) fn intersect(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    dec: &mut Decider,
) -> Answered<IntersectionGraph> {
    section(model, a, b, true, dec)
}

/// The section graph. `prune` lets finite-face exclusions skip curve/surface
/// queries; the graph is the same either way.
pub(crate) fn section(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    prune: bool,
    dec: &mut Decider,
) -> Answered<IntersectionGraph> {
    assert_ne!(a, b, "boolean operands must be distinct bodies");
    for body in [a, b] {
        if !matches!(model.body(body).kind, BodyKind::Solid | BodyKind::Empty) {
            return Err(refusal(
                "boolean.not-solid",
                "operand is not solid or empty",
            ));
        }
    }
    let first = model.faces_of_body(a);
    let second = model.faces_of_body(b);
    let algebraic = first.iter().chain(&second).any(|face| {
        !model
            .geom()
            .surface(model.face(*face).surface)
            .has_rational_definition()
    });
    if !algebraic && boundaries_apart(model, &first, &second, dec) {
        return Ok(Answer::exact(IntersectionGraph {
            a,
            b,
            edges: Vec::new(),
            points: Vec::new(),
            singular: BTreeSet::new(),
            marks: MarkLog::new(),
        }));
    }
    let mut bounds = crate::bounds::FaceBounds::new(algebraic);
    let mut regions = RegionCache::default();
    let mut pool = Vec::new();
    let mut points = Vec::new();
    let mut mandatory = BTreeSet::new();
    let mut singular = BTreeSet::new();
    for body in [a, b] {
        for vertex in model.vertices_of_body(body) {
            let point = canonical_point(
                model,
                &mut points,
                model.vertex(vertex).point,
                dec,
                "boolean.imprint-ambiguous",
            )?;
            mandatory.insert(point);
        }
    }
    let mut pairs: BTreeMap<(SurfaceId, SurfaceId), SurfacePair> = BTreeMap::new();
    // Preserve every canonical carrier and its insertion order. Optional bound
    // witnesses stay on a fork, outside the graph's canonical point pool.
    // Surfaces proven to contain each carrier: both supports of a section
    // component, or the face surface carrying a coincident face's edge.
    let mut supports: Vec<BTreeSet<SurfaceId>> = Vec::new();
    let mut pruning_model = (prune && !algebraic).then(|| model.clone());
    let mut finite_bounds = crate::bounds::FaceBounds::new(false);
    let mut active = Vec::new();
    for &fa in &first {
        for &fb in &second {
            let sa = model.face(fa).surface;
            let sb = model.face(fb).surface;
            if algebraic && bounds.disjoint(model, fa, fb, dec)? {
                continue;
            }
            let pair = if let Some(value) = pairs.get(&(sa, sb)) {
                value.clone()
            } else {
                let answer = match intersect_surfaces(model.geom_mut(), sa, sb, dec) {
                    Ok(answer) => answer,
                    Err(refusal) => {
                        if refusal.code == "intersect.traced-pair-required"
                            && bounds.disjoint(model, fa, fb, dec)?
                        {
                            // This proves only this trimmed pair disjoint. Never
                            // cache it as a statement about the whole supports.
                            continue;
                        }
                        return Err(refusal);
                    }
                };
                assert!(
                    answer.marks.is_empty(),
                    "surface configuration is rational and exact"
                );
                pairs.insert((sa, sb), answer.value.clone());
                answer.value
            };
            // Disjoint and point-only supports have no carrier activity to
            // prove. Cache the optional finite bound only when adding a curve.
            let mut relevance = None;
            let mut relevant = |dec: &mut Decider| {
                *relevance.get_or_insert_with(|| match &mut pruning_model {
                    Some(probe) => !finite_bounds.excludes(probe, fa, fb, dec),
                    None => true,
                })
            };
            match pair {
                SurfacePair::Disjoint => {}
                SurfacePair::Coincident { .. } => {
                    for face in [fa, fb] {
                        let surface = model.face(face).surface;
                        for curve in face_edges(model, face) {
                            let relevant = relevant(dec);
                            let index = add_curve(model, &mut pool, curve, dec)?;
                            active.resize(pool.len(), false);
                            active[index] |= relevant;
                            supports.resize(pool.len(), BTreeSet::new());
                            supports[index].insert(surface);
                        }
                    }
                }
                SurfacePair::Meets(parts) => {
                    for part in parts {
                        match part {
                            Component::Point { point } => {
                                canonical_point(
                                    model,
                                    &mut points,
                                    point,
                                    dec,
                                    "boolean.imprint-ambiguous",
                                )?;
                            }
                            Component::Curve {
                                curve,
                                singular: nodes,
                                ..
                            } => {
                                let relevant = relevant(dec);
                                let index = add_curve(model, &mut pool, curve, dec)?;
                                active.resize(pool.len(), false);
                                active[index] |= relevant;
                                supports.resize(pool.len(), BTreeSet::new());
                                supports[index].extend([sa, sb]);
                                for node in nodes {
                                    let point = canonical_point(
                                        model,
                                        &mut points,
                                        node,
                                        dec,
                                        "boolean.imprint-ambiguous",
                                    )?;
                                    mandatory.insert(point);
                                    singular.insert(point);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let surfaces: Vec<SurfaceId> = first
        .iter()
        .chain(&second)
        .map(|f| model.face(*f).surface)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut marks = MarkLog::new();
    // Surfaces whose whole meeting with each carrier is now canonical.
    let mut complete: Vec<BTreeSet<SurfaceId>> = vec![BTreeSet::new(); pool.len()];
    let mut order = Order {
        finder: points.iter().map(|&p| (p, 0)).collect(),
        ..Order::default()
    };
    for (index, &curve) in pool.iter().enumerate() {
        // An inactive carrier has no potentially meeting generating face pair.
        // The complete curve/pair pass below still includes its crossings.
        if !active[index] {
            continue;
        }
        let known = points.len();
        let (_, found_marks, found_complete) =
            carrier_hits(model, curve, &surfaces, &mut points, dec)?;
        marks.extend(&found_marks);
        complete[index] = found_complete;
        for &point in &points[known..] {
            order.finder.insert(point, index);
        }
    }
    for i in 0..pool.len() {
        for j in i + 1..pool.len() {
            // Two carriers meet only on each other's supports. When one
            // carrier's complete meeting with the other's support is already
            // canonical, the pair can add no point; equality is never inferred.
            // A skipped carrier has no complete meetings of its own. Its pair
            // with an earlier active carrier that covers it can still be
            // skipped: every crossing was then found by that carrier or before
            // it, so the complete pass would not find it first through the
            // skipped one.
            let covered = |c: usize, d: usize| !supports[c].is_disjoint(&complete[d]);
            let skip = match (active[i], active[j]) {
                (true, true) => covered(i, j) || covered(j, i),
                (true, false) => covered(j, i),
                (false, _) => false,
            };
            if skip {
                continue;
            }
            let crate::curve_pair::Meetings::Points(meetings) =
                crate::curve_pair::meetings(model, pool[i], pool[j], dec)?
            else {
                continue;
            };
            for point in meetings {
                let canonical =
                    canonical_point(model, &mut points, point, dec, "boolean.imprint-ambiguous")?;
                // The complete pass could find this crossing first through
                // the skipped carrier.
                let skipped = if active[i] { j } else { i };
                if active[i] != active[j]
                    && order
                        .finder
                        .get(&canonical)
                        .is_none_or(|&found| skipped < found)
                {
                    order.crossers.entry(canonical).or_default().insert(skipped);
                }
            }
        }
    }
    order.index = points.iter().enumerate().map(|(i, &p)| (p, i)).collect();
    let mut point_faces = BTreeMap::new();
    for &point in &points {
        if let Some(faces) = shared_faces(model, &first, &second, point, &mut regions, dec)? {
            point_faces.insert(point, faces);
        }
    }
    let mut edges = Vec::new();
    for &curve in &pool {
        let mut cuts = Vec::new();
        for &point in &points {
            if coincident(
                on_curve(model.geom(), curve, point, dec),
                "boolean.imprint-ambiguous",
            )? {
                cuts.push(point);
            }
        }
        sort_points(model, curve, &mut cuts, dec)?;
        let closed = curve_is_closed(model.geom(), curve, dec)?;
        let context = Context {
            first: &first,
            second: &second,
            mandatory: &mandatory,
            point_faces: &point_faces,
        };
        let mut pieces = carrier_pieces(model, curve, &cuts, closed, &context, &mut regions, dec)?;
        // A closed carrier starts at the cut the pass found first. Restore the
        // complete pass's choice when a skipped carrier could change it.
        if closed && !pieces.is_empty() {
            let start = order.first(model, &pool, &surfaces, &cuts, dec)?;
            if start != 0 {
                cuts.rotate_left(start);
                pieces = carrier_pieces(model, curve, &cuts, closed, &context, &mut regions, dec)?;
            }
        }
        edges.extend(pieces);
    }
    let mut isolated = Vec::new();
    for point in points {
        if !point_faces.contains_key(&point) {
            continue;
        }
        let mut on_section = false;
        for edge in &edges {
            if on_arc(
                model,
                edge.curve,
                edge.start,
                edge.end,
                point,
                dec,
                "boolean.imprint-ambiguous",
            )? {
                on_section = true;
                break;
            }
        }
        if !on_section {
            isolated.push(point);
        }
    }
    if isolated.iter().any(|p| order.crossers.contains_key(p)) {
        let mut keyed = Vec::with_capacity(isolated.len());
        for point in isolated {
            keyed.push((order.key(model, &pool, &surfaces, point, dec)?, point));
        }
        keyed.sort_unstable_by_key(|&(key, _)| key);
        isolated = keyed.into_iter().map(|(_, point)| point).collect();
    }
    let graph = IntersectionGraph {
        a,
        b,
        edges,
        points: isolated,
        singular,
        marks: marks.clone(),
    };
    Ok(Answer::with_marks(graph, marks))
}
