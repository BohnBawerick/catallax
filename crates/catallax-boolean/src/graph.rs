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
    // Surfaces proven to contain each carrier: both supports of a section
    // component, or the face surface carrying a coincident face's edge.
    let mut supports: Vec<BTreeSet<SurfaceId>> = Vec::new();
    // Preserve every canonical carrier and its insertion order. Optional bound
    // witnesses stay on a fork, outside the graph's canonical point pool.
    let mut pruning_model = (!algebraic).then(|| model.clone());
    let mut finite_bounds = crate::bounds::FaceBounds::new(false);
    let mut active = Vec::new();
    // A coincident face pair that may overlap carries a two-dimensional common
    // region, where a carrier can lie inside both faces without generating it.
    let mut overlapping_coincidence = false;
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
                    overlapping_coincidence |= relevant(dec);
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
    for (index, &curve) in pool.iter().enumerate() {
        // An inactive carrier has no potentially meeting generating face pair.
        // The complete curve/pair pass below still includes its crossings.
        if !active[index] {
            continue;
        }
        for &surface in &surfaces {
            let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
            if !matches!(answer.value, CurveSurface::Contained) {
                complete[index].insert(surface);
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
                    canonical_point(
                        model,
                        &mut points,
                        hit.point,
                        dec,
                        "boolean.imprint-ambiguous",
                    )?;
                }
            }
            marks.extend(&answer.marks);
        }
    }
    for i in 0..pool.len() {
        for j in i + 1..pool.len() {
            // Two carriers meet only on each other's supports. When one
            // carrier's complete meeting with the other's support is already
            // canonical, the pair can add no point; equality is never inferred.
            let covered = |c: usize, d: usize| !supports[c].is_disjoint(&complete[d]);
            if covered(i, j) || covered(j, i) {
                continue;
            }
            let crate::curve_pair::Meetings::Points(meetings) =
                crate::curve_pair::meetings(model, pool[i], pool[j], dec)?
            else {
                continue;
            };
            for point in meetings {
                canonical_point(model, &mut points, point, dec, "boolean.imprint-ambiguous")?;
            }
        }
    }
    let mut point_faces = BTreeMap::new();
    for &point in &points {
        if let Some(faces) = shared_faces(model, &first, &second, point, &mut regions, dec)? {
            point_faces.insert(point, faces);
        }
    }
    let mut edges = Vec::new();
    for (index, &curve) in pool.iter().enumerate() {
        // A section arc lies in both supports of an overlapping face pair, so
        // its carrier was generated by that pair. An inactive carrier without
        // overlapping coincidence has no arc on both boundaries.
        if !active[index] && !overlapping_coincidence {
            continue;
        }
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
                shared_faces(model, &first, &second, sample, &mut regions, dec)?
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
                    !mandatory.contains(&p)
                        && point_faces.get(&p).is_some_and(|(first, second)| {
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
                    !mandatory.contains(&p)
                        && point_faces.get(&p).is_some_and(|(first, second)| {
                            first == &pieces[0].faces_of_first
                                && second == &pieces[0].faces_of_second
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
                !mandatory.contains(&p)
                    && point_faces.get(&p).is_some_and(|(first, second)| {
                        first == &pieces[0].faces_of_first && second == &pieces[0].faces_of_second
                    })
            })
        {
            pieces[0].start = None;
            pieces[0].end = None;
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
