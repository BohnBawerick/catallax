//! Boundary intersections trimmed to both operands, with canonical exact
//! points and curves. Temporary carrier cuts are removed from the answer.

use crate::geometry::{canonical_point, coincident, on_arc, refusal};
use crate::points::curve_sample;
use crate::region::{Location, locate};
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
) -> Result<(), Refusal> {
    for &other in pool.iter() {
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
            return Ok(());
        }
    }
    pool.push(curve);
    Ok(())
}

pub(crate) fn containing_faces(
    model: &mut Model,
    faces: &[FaceId],
    point: PointId,
    dec: &mut Decider,
) -> Result<Vec<FaceId>, Refusal> {
    let mut result = Vec::new();
    for &face in faces {
        if coincident(
            model.geom().side(model.face(face).surface, point, dec),
            "boolean.imprint-ambiguous",
        )? && locate(model, face, point, dec)? != Location::Outside
        {
            result.push(face);
        }
    }
    Ok(result)
}

fn same_faces(a: &SectionEdge, b: &SectionEdge) -> bool {
    a.faces_of_first == b.faces_of_first && a.faces_of_second == b.faces_of_second
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
    let mut bounds = crate::bounds::FaceBounds::new(algebraic);
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
    for &fa in &first {
        for &fb in &second {
            let sa = model.face(fa).surface;
            let sb = model.face(fb).surface;
            // Algebraic supports can meet far outside their finite trims. Use
            // the certified section bounds before constructing those branches;
            // an overlapping box leaves the ordinary intersection path intact.
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
            match pair {
                SurfacePair::Disjoint => {}
                SurfacePair::Coincident { .. } => {
                    for face in [fa, fb] {
                        for curve in face_edges(model, face) {
                            add_curve(model, &mut pool, curve, dec)?;
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
                                add_curve(model, &mut pool, curve, dec)?;
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
    for &curve in &pool {
        for &surface in &surfaces {
            let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
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
        point_faces.insert(
            point,
            (
                containing_faces(model, &first, point, dec)?,
                containing_faces(model, &second, point, dec)?,
            ),
        );
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
            let faces_of_first = containing_faces(model, &first, sample, dec)?;
            let faces_of_second = containing_faces(model, &second, sample, dec)?;
            if faces_of_first.is_empty() || faces_of_second.is_empty() {
                continue;
            }
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
                        && point_faces[&p].0 == edge.faces_of_first
                        && point_faces[&p].1 == edge.faces_of_second
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
                        && point_faces[&p].0 == pieces[0].faces_of_first
                        && point_faces[&p].1 == pieces[0].faces_of_second
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
                    && point_faces[&p].0 == pieces[0].faces_of_first
                    && point_faces[&p].1 == pieces[0].faces_of_second
            })
        {
            pieces[0].start = None;
            pieces[0].end = None;
        }
        edges.extend(pieces);
    }
    let mut isolated = Vec::new();
    for point in points {
        if point_faces[&point].0.is_empty() || point_faces[&point].1.is_empty() {
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
