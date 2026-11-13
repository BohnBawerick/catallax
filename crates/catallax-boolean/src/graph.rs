//! Boundary intersections trimmed to both operands, with canonical exact
//! points and curves. Temporary carrier cuts are removed from the answer.

use crate::geometry::{canonical_point, coincident, on_arc, refusal, undecided};
use crate::points::curve_sample;
use crate::region::{Location, RegionCache};
use crate::stages::{IntersectionGraph, SectionEdge};
use catallax_intersect::{
    Component, CurveSurface, Cyclic, SurfacePair, curve_is_closed, cyclic_order,
    intersect_curve_surface, intersect_curve_surface_with_line_parameter, intersect_surfaces,
    on_curve, order_on_curve, same_curve,
};
use catallax_numbers::{
    Answer, Answered, Ball, Decider, Fact, FactKind, FactOrigin, MarkLog, Outcome, Refusal,
    StrictOrder, ZeroProof,
};
use catallax_shapes::{CurveId, GeomKey, PointId, SurfaceId};
use catallax_topology::{BodyId, BodyKind, FaceId, Loop, Model, VertexId};
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
                        return Err(undecided("boolean.imprint-ambiguous", n));
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
                        return Err(undecided("boolean.imprint-ambiguous", n));
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
    // The first equal carrier ends the search; workers compare against later
    // carriers ahead. Merging a carrier into the shared pool is an ordered
    // append.
    #[cfg(feature = "work-span")]
    let _append = catallax_numbers::work_span::append();
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::search("carrier merge");
    let mut ahead = catallax_numbers::run_ahead!(dec, pool.len());
    for (index, &other) in pool.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, pool) = (model.clone(), pool.clone());
            move |index: usize, dec: &mut Decider| {
                let _ = same_curve(model.geom(), curve, pool[index], dec);
            }
        });
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
/// Returns the raw hits, their marks and, for each surface whose whole meeting
/// with `curve` is now canonical, that meeting's canonical points.
fn carrier_hits(
    model: &mut Model,
    curve: CurveId,
    surfaces: &[SurfaceId],
    points: &mut Vec<PointId>,
    exclusions: bool,
    trimmed: &mut crate::trim_bounds::TrimBounds,
    dec: &mut Decider,
) -> Result<(Vec<PointId>, MarkLog, Complete), Refusal> {
    let query = if trimmed.has_fitted_band() {
        intersect_curve_surface_with_line_parameter
    } else {
        intersect_curve_surface
    };
    let mut found = Vec::new();
    let mut marks = MarkLog::new();
    let mut complete = Complete::new();
    let enclosure = exclusions
        .then(|| crate::carrier_bounds::enclosure(model.geom(), curve))
        .flatten();
    let mut validated = false;
    for &surface in surfaces {
        if validated
            && let Some(bounds) = &enclosure
            && crate::carrier_bounds::excludes(bounds, model.geom().surface(surface))
        {
            // A strict exclusion proves the complete meeting empty, including
            // all tangencies and isolated points. Retain its coverage proof.
            complete.insert(surface, Vec::new());
            continue;
        }
        let answer = match query(model.geom_mut(), curve, surface, dec) {
            Ok(answer) => answer,
            Err(refusal) => {
                // A refused query on a bounded support is empty when the
                // complete carrier's box misses the complete support's box.
                if matches!(
                    refusal.code,
                    "intersect.unsupported-curve"
                        | "intersect.unsupported-pair"
                        | "intersect.traced-pair-required"
                        | "intersect.parametric-surface"
                ) && let Some(bounds) = crate::carrier_bounds::enclosure(model.geom(), curve)
                    && trimmed.carrier_misses_support(model, curve, &bounds, surface, dec)
                {
                    complete.insert(surface, Vec::new());
                    continue;
                }
                return Err(named_curve(model, curve, surface, refusal));
            }
        };
        // Keep the first ordinary query's validation/refusal path unchanged.
        validated |= !matches!(
            model.geom().surface(surface),
            catallax_shapes::Surface::TorusSheet(_)
        );
        if !matches!(answer.value, CurveSurface::Contained) {
            complete.insert(surface, Vec::new());
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
                let canonical =
                    canonical_point(model, points, hit.point, dec, "boolean.imprint-ambiguous")?;
                complete.get_mut(&surface).unwrap().push(canonical);
                found.push(hit.point);
            }
        }
        marks.extend(&answer.marks);
    }
    Ok((found, marks, complete))
}

/// A refused surface pair, named by the two support kinds it needs.
fn named_pair(model: &Model, a: SurfaceId, b: SurfaceId, mut refusal: Refusal) -> Refusal {
    let [a, b] = [a, b].map(|surface| model.geom().surface(surface));
    // Analytic refusals already name their pair. Spline pairs did not.
    if !matches!(a, catallax_shapes::Surface::BSpline(_))
        && !matches!(b, catallax_shapes::Surface::BSpline(_))
    {
        return refusal;
    }
    refusal.detail = format!(
        "{} with {}: {}",
        crate::trim_bounds::surface_name(a),
        crate::trim_bounds::surface_name(b),
        refusal.detail
    );
    refusal
}

/// A refused curve query on a spline support, named by the two kinds.
fn named_curve(model: &Model, curve: CurveId, surface: SurfaceId, mut refusal: Refusal) -> Refusal {
    let surface = model.geom().surface(surface);
    if !matches!(surface, catallax_shapes::Surface::BSpline(_)) {
        return refusal;
    }
    let curve = match model.geom().curve(curve) {
        catallax_shapes::Curve::Line(_) => "line",
        catallax_shapes::Curve::Circle(_) => "circle",
        catallax_shapes::Curve::Ellipse(_) => "ellipse",
        catallax_shapes::Curve::Intersection(_) => "traced section curve",
        catallax_shapes::Curve::BSpline(_) => "spline curve",
        _ => "curve",
    };
    refusal.detail = format!(
        "{curve} with {}: {}",
        crate::trim_bounds::surface_name(surface),
        refusal.detail
    );
    refusal
}

/// Each surface whose whole meeting with a carrier is canonical, with the
/// canonical points of that meeting.
type Complete = BTreeMap<SurfaceId, Vec<PointId>>;

/// Borrowed graph state shared by every carrier's pieces.
struct Context<'a> {
    first: &'a [FaceId],
    second: &'a [FaceId],
    multiway: bool,
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
        let Some((faces_of_first, faces_of_second)) = shared_faces(
            model,
            context.first,
            context.second,
            context.multiway,
            sample,
            regions,
            dec,
        )?
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
#[derive(Clone, Default)]
struct Order {
    /// The active carrier whose hits first found each point; seeds are 0.
    finder: BTreeMap<PointId, usize>,
    /// Skipped carriers that cross each point before its finder.
    crossers: BTreeMap<PointId, BTreeSet<usize>>,
    /// Position in the pruned discovery order.
    index: BTreeMap<PointId, usize>,
    /// Hits of each queried skipped carrier, in query order.
    hits: BTreeMap<usize, Vec<PointId>>,
    /// Whether complete-carrier bounds may prove support queries empty.
    carrier_exclusions: bool,
    /// Boxes that settle a refused support query, as in the complete pass.
    trimmed: crate::trim_bounds::TrimBounds,
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
                let (found, ..) = carrier_hits(
                    model,
                    pool[skipped],
                    surfaces,
                    &mut Vec::new(),
                    self.carrier_exclusions,
                    &mut self.trimmed,
                    dec,
                )?;
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
        if crate::geometry::on_surface(
            model,
            model.face(face).surface,
            point,
            dec,
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
    multiway: bool,
    point: PointId,
    regions: &mut RegionCache,
    dec: &mut Decider,
) -> Result<Option<SharedFaces>, Refusal> {
    if multiway {
        let faces = containing_faces(model, first, point, regions, dec)?;
        return Ok((faces.len() >= 2).then(|| (faces, Vec::new())));
    }
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
    section(model, a, b, true, true, dec)
}

struct Inputs<'a> {
    bodies: [BodyId; 2],
    first: &'a [FaceId],
    second: &'a [FaceId],
    vertices: Vec<VertexId>,
    multiway: bool,
    clip: Option<[Ball; 3]>,
    carrier_exclusions: bool,
}

pub(crate) fn arrangement(
    model: &mut Model,
    source: BodyId,
    faces: &[FaceId],
    dec: &mut Decider,
) -> Answered<IntersectionGraph> {
    let mut vertices = BTreeSet::new();
    for &face in faces {
        for lp in &model.face(face).loops {
            match model.loop_(*lp) {
                Loop::Vertex(vertex) => {
                    vertices.insert(*vertex);
                }
                Loop::Coedges(coedges) => {
                    for coedge in coedges {
                        let edge = model.edge(model.coedge(*coedge).edge);
                        vertices.extend([edge.start, edge.end].into_iter().flatten());
                    }
                }
            }
        }
    }
    let clip = crate::bounds::body_box(model, source, dec)?;
    section_inputs(
        model,
        Inputs {
            bodies: [source, source],
            first: faces,
            second: faces,
            vertices: vertices.into_iter().collect(),
            multiway: true,
            clip,
            carrier_exclusions: true,
        },
        true,
        dec,
    )
}

/// The section graph. `prune` lets finite-face exclusions skip curve/surface
/// queries; the graph is the same either way.
pub(crate) fn section(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    prune: bool,
    carrier_exclusions: bool,
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
    let vertices = [a, b]
        .into_iter()
        .flat_map(|body| model.vertices_of_body(body))
        .collect();
    section_inputs(
        model,
        Inputs {
            bodies: [a, b],
            first: &first,
            second: &second,
            vertices,
            multiway: false,
            clip: None,
            carrier_exclusions,
        },
        prune,
        dec,
    )
}

fn section_inputs(
    model: &mut Model,
    inputs: Inputs<'_>,
    prune: bool,
    dec: &mut Decider,
) -> Answered<IntersectionGraph> {
    let Inputs {
        bodies: [a, b],
        first,
        second,
        vertices,
        multiway,
        clip,
        carrier_exclusions,
    } = inputs;
    let algebraic = first.iter().chain(second).any(|face| {
        !model
            .geom()
            .surface(model.face(*face).surface)
            .has_rational_definition()
    });
    if !multiway && !algebraic && boundaries_apart(model, first, second, dec) {
        return Ok(Answer::exact(IntersectionGraph {
            a,
            b,
            edges: Vec::new(),
            points: Vec::new(),
            singular: BTreeSet::new(),
            marks: MarkLog::new(),
        }));
    }
    let mut pass = PairPass {
        bounds: crate::bounds::FaceBounds::new(algebraic),
        trimmed: crate::trim_bounds::TrimBounds::for_faces(model, first.iter().chain(second)),
        finite_bounds: crate::bounds::FaceBounds::new(false),
        pruning_model: None,
        pairs: BTreeMap::new(),
        #[cfg(feature = "work-span")]
        pairs_span: catallax_numbers::work_span::Instance::new(),
        pool: Vec::new(),
        points: Vec::new(),
        mandatory: BTreeSet::new(),
        singular: BTreeSet::new(),
        supports: Vec::new(),
        active: Vec::new(),
        overlapping_coincidence: false,
    };
    let mut regions = RegionCache::default();
    for vertex in vertices {
        let point = canonical_point(
            model,
            &mut pass.points,
            model.vertex(vertex).point,
            dec,
            "boolean.imprint-ambiguous",
        )?;
        pass.mandatory.insert(point);
    }
    pass.pruning_model = (prune && !algebraic).then(|| model.clone());
    let face_pairs: Vec<(FaceId, FaceId)> = first
        .iter()
        .flat_map(|&fa| second.iter().map(move |&fb| (fa, fb)))
        .filter(|&(fa, fb)| !(multiway && fb <= fa))
        .collect();
    // Each pair's support intersection is its own proof. Merging its carriers
    // and points into the pass stays in pair order on this thread; workers
    // intersect later pairs ahead against the pass as it stood.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("face pairs");
    let mut ahead = catallax_numbers::run_ahead!(dec, face_pairs.len());
    for (index, &pair) in face_pairs.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, pass, face_pairs, clip) = (
                model.clone(),
                pass.clone(),
                face_pairs.clone(),
                clip.clone(),
            );
            move |index: usize, dec: &mut Decider| {
                let _ = face_pair(
                    &mut model.clone(),
                    &mut pass.clone(),
                    face_pairs[index],
                    algebraic,
                    clip.as_ref(),
                    dec,
                );
            }
        });
        face_pair(model, &mut pass, pair, algebraic, clip.as_ref(), dec)?;
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let PairPass {
        mut trimmed,
        pool,
        mut points,
        mandatory,
        singular,
        supports,
        mut active,
        overlapping_coincidence,
        ..
    } = pass;
    // Otherwise every section piece lies on a carrier of a potentially meeting
    // pair, so a skipped carrier keeps no piece and a point found only on
    // skipped carriers is neither a kept cut nor an isolated point.
    if overlapping_coincidence {
        active.fill(true);
    }
    let surfaces: Vec<SurfaceId> = first
        .iter()
        .chain(second)
        .map(|f| model.face(*f).surface)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut marks = MarkLog::new();
    // Surfaces whose whole meeting with each carrier is now canonical.
    let mut complete: Vec<Complete> = vec![Complete::new(); pool.len()];
    let mut order = Order {
        finder: points.iter().map(|&p| (p, 0)).collect(),
        carrier_exclusions,
        trimmed: crate::trim_bounds::TrimBounds::for_faces(model, first.iter().chain(second)),
        ..Order::default()
    };
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("carrier hits");
    let mut ahead = catallax_numbers::run_ahead!(dec, pool.len());
    for (index, &curve) in pool.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        // An inactive carrier has no potentially meeting generating face pair.
        // The complete curve/pair pass below still includes its crossings.
        if !active[index] {
            continue;
        }
        ahead.at(index, || {
            let (model, pool, active, surfaces, points, trimmed) = (
                model.clone(),
                pool.clone(),
                active.clone(),
                surfaces.clone(),
                points.clone(),
                trimmed.clone(),
            );
            move |index: usize, dec: &mut Decider| {
                if active[index] {
                    let _ = carrier_hits(
                        &mut model.clone(),
                        pool[index],
                        &surfaces,
                        &mut points.clone(),
                        carrier_exclusions,
                        &mut trimmed.clone(),
                        dec,
                    );
                }
            }
        });
        let known = points.len();
        let (_, found_marks, found_complete) = carrier_hits(
            model,
            curve,
            &surfaces,
            &mut points,
            carrier_exclusions,
            &mut trimmed,
            dec,
        )?;
        marks.extend(&found_marks);
        complete[index] = found_complete;
        for &point in &points[known..] {
            order.finder.insert(point, index);
        }
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let carrier_pairs: Vec<(usize, usize)> = (0..pool.len())
        .flat_map(|i| (i + 1..pool.len()).map(move |j| (i, j)))
        .collect();
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("carrier pairs");
    let mut ahead = catallax_numbers::run_ahead!(dec, carrier_pairs.len());
    for (index, &pair) in carrier_pairs.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, pairs, pool, active, supports, complete, points, order) = (
                model.clone(),
                carrier_pairs.clone(),
                pool.clone(),
                active.clone(),
                supports.clone(),
                complete.clone(),
                points.clone(),
                order.clone(),
            );
            move |index: usize, dec: &mut Decider| {
                let carriers = Carriers {
                    pool: &pool,
                    active: &active,
                    supports: &supports,
                    complete: &complete,
                };
                // Most pairs are settled by their supports alone; copy the
                // pass only for a pair that needs a proof.
                if carriers.settled_without_proof(pairs[index]) {
                    return;
                }
                let _ = carrier_pair(
                    &mut model.clone(),
                    pairs[index],
                    carriers,
                    &mut points.clone(),
                    &mut order.clone(),
                    dec,
                );
            }
        });
        let carriers = Carriers {
            pool: &pool,
            active: &active,
            supports: &supports,
            complete: &complete,
        };
        carrier_pair(model, pair, carriers, &mut points, &mut order, dec)?;
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    order.index = points.iter().enumerate().map(|(i, &p)| (p, i)).collect();
    let mut point_faces = BTreeMap::new();
    // The faces that hold each point are found for each point on its own.
    // Other threads find them for later points whole, each on a copy of the
    // model and with the region cache as they were when the loop was handed
    // over. This thread takes such a result when the point would have been
    // placed the same on the model as it is now, and the other thread's copy
    // of the cache held all the point asked for.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("point faces");
    let mut ahead = catallax_numbers::values_ahead!(dec, points.len());
    for (index, &point) in points.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        let found = ahead.take_if(
            index,
            dec,
            || {
                let (base, first, second, points, regions) = (
                    std::sync::Arc::new(model.clone()),
                    first.to_vec(),
                    second.to_vec(),
                    points.clone(),
                    regions.clone(),
                );
                move |index: usize, dec: &mut Decider| {
                    regions.ahead(&base, |model, regions| {
                        shared_faces(
                            model,
                            &first,
                            &second,
                            multiway,
                            points[index],
                            regions,
                            dec,
                        )
                    })
                }
            },
            |ahead| ahead.holds(model),
        );
        let faces = match found {
            Some(ahead) => ahead.commit(model),
            None => shared_faces(model, first, second, multiway, point, &mut regions, dec),
        }?;
        if let Some(faces) = faces {
            point_faces.insert(point, faces);
        }
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let mut edges = Vec::new();
    // Each carrier is cut into its section pieces on its own, and handed
    // over like the points above. A carrier that has to ask where a skipped
    // carrier crosses adds to the discovery order's own lists, which the
    // carriers after it read; another thread's result for such a carrier is
    // not taken, and this thread cuts it itself.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("carrier pieces");
    let mut ahead = catallax_numbers::values_ahead!(dec, pool.len());
    for (index, &curve) in pool.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        let found = ahead.take_if(
            index,
            dec,
            || {
                let (base, pool, surfaces, points, first, second, mandatory, point_faces) = (
                    std::sync::Arc::new(model.clone()),
                    pool.clone(),
                    surfaces.clone(),
                    points.clone(),
                    first.to_vec(),
                    second.to_vec(),
                    mandatory.clone(),
                    point_faces.clone(),
                );
                let (regions, order) = (regions.clone(), order.clone());
                move |index: usize, dec: &mut Decider| {
                    let context = Context {
                        first: &first,
                        second: &second,
                        multiway,
                        mandatory: &mandatory,
                        point_faces: &point_faces,
                    };
                    let mut order = order.clone();
                    let asked = order.hits.len();
                    let ahead = regions.ahead(&base, |model, regions| {
                        carrier_sections(
                            model,
                            pool[index],
                            (&pool, &surfaces, &points),
                            &context,
                            regions,
                            &mut order,
                            dec,
                        )
                    });
                    (ahead, order.hits.len() != asked)
                }
            },
            |(ahead, asked)| !*asked && ahead.holds(model),
        );
        let context = Context {
            first,
            second,
            multiway,
            mandatory: &mandatory,
            point_faces: &point_faces,
        };
        edges.extend(match found {
            Some((ahead, _)) => ahead.commit(model),
            None => carrier_sections(
                model,
                curve,
                (&pool, &surfaces, &points),
                &context,
                &mut regions,
                &mut order,
                dec,
            ),
        }?);
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let mut isolated = Vec::new();
    // Whether a point lies on a section is a value of what every thread
    // shares, so other threads test later points whole.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("isolated points");
    let mut ahead = catallax_numbers::values_ahead!(dec, points.len());
    for (index, &point) in points.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        if !point_faces.contains_key(&point) {
            continue;
        }
        let adopted = ahead.take(index, dec, || {
            let (model, points, edges) = (model.clone(), points.clone(), edges.clone());
            let kept: BTreeSet<PointId> = point_faces.keys().copied().collect();
            move |index: usize, dec: &mut Decider| {
                // The loop skips a point with no faces; its value is not read.
                if !kept.contains(&points[index]) {
                    return Ok(true);
                }
                on_section(&model, &edges, points[index], dec)
            }
        });
        let on = match adopted {
            Some(on) => on,
            None => on_section(model, &edges, point, dec),
        }?;
        if !on {
            isolated.push(point);
        }
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
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

/// What the face-pair pass builds: the canonical carriers and points, in the
/// order the pairs produce them.
#[derive(Clone)]
struct PairPass {
    bounds: crate::bounds::FaceBounds,
    trimmed: crate::trim_bounds::TrimBounds,
    finite_bounds: crate::bounds::FaceBounds,
    // Optional bound witnesses stay on this fork, outside the graph's
    // canonical point pool.
    pruning_model: Option<Model>,
    pairs: BTreeMap<(SurfaceId, SurfaceId), SurfacePair>,
    // Work-span observation only: names the pair cache's entries.
    #[cfg(feature = "work-span")]
    pairs_span: catallax_numbers::work_span::Instance,
    // Every canonical carrier in insertion order.
    pool: Vec<CurveId>,
    points: Vec<PointId>,
    mandatory: BTreeSet<PointId>,
    singular: BTreeSet<PointId>,
    // Surfaces proven to contain each carrier: both supports of a section
    // component, or the face surface carrying a coincident face's edge.
    supports: Vec<BTreeSet<SurfaceId>>,
    active: Vec<bool>,
    // A coincident face pair that may overlap has a two-dimensional common
    // region, where a carrier can keep a section piece without a potentially
    // meeting generating pair.
    overlapping_coincidence: bool,
}

/// Intersects one face pair's supports and merges what they share into `pass`.
fn face_pair(
    model: &mut Model,
    pass: &mut PairPass,
    (fa, fb): (FaceId, FaceId),
    algebraic: bool,
    clip: Option<&[Ball; 3]>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if let Some(clip) = clip
        && pass
            .finite_bounds
            .excludes_clipped(model, fa, fb, clip, dec)
    {
        return Ok(());
    }
    let sa = model.face(fa).surface;
    let sb = model.face(fb).surface;
    if algebraic && pass.bounds.disjoint(model, fa, fb, dec)? {
        return Ok(());
    }
    #[cfg(feature = "work-span")]
    catallax_numbers::work_span::read_at(
        "surface pair cache",
        pass.pairs_span.id(),
        (u64::from(sa.index()) << 32) | u64::from(sb.index()),
    );
    let pair = if let Some(value) = pass.pairs.get(&(sa, sb)) {
        value.clone()
    } else {
        let answer = match intersect_surfaces(model.geom_mut(), sa, sb, dec) {
            Ok(answer) => answer,
            Err(refusal) => {
                if refusal.code == "intersect.traced-pair-required"
                    && pass.bounds.disjoint(model, fa, fb, dec)?
                {
                    // This proves only this trimmed pair disjoint. Never
                    // cache it as a statement about the whole supports.
                    return Ok(());
                }
                // Only a refused pair next to a fitted band reaches
                // this wider box proof, so every other pair keeps its
                // path and its work.
                if matches!(
                    refusal.code,
                    "intersect.traced-pair-required" | "intersect.unsupported-pair"
                ) && pass.trimmed.separated(model, fa, fb, dec)
                {
                    return Ok(());
                }
                if !matches!(
                    refusal.code,
                    "intersect.traced-pair-required" | "intersect.unsupported-pair"
                ) {
                    return Err(named_pair(model, sa, sb, refusal));
                }
                match catallax_intersect::intersect_surfaces_traced(model.geom_mut(), sa, sb, dec) {
                    Ok(answer) => answer,
                    Err(refusal) => return Err(named_pair(model, sa, sb, refusal)),
                }
            }
        };
        assert!(
            answer.marks.is_empty(),
            "surface configuration is rational and exact"
        );
        pass.pairs.insert((sa, sb), answer.value.clone());
        #[cfg(feature = "work-span")]
        catallax_numbers::work_span::write_at(
            "surface pair cache",
            pass.pairs_span.id(),
            (u64::from(sa.index()) << 32) | u64::from(sb.index()),
        );
        answer.value
    };
    let PairPass {
        finite_bounds,
        pruning_model,
        pool,
        points,
        mandatory,
        singular,
        supports,
        active,
        overlapping_coincidence,
        ..
    } = pass;
    // Disjoint and point-only supports have no carrier activity to
    // prove. Cache the optional finite bound only when adding a curve.
    let mut relevance = None;
    let mut relevant = |dec: &mut Decider| {
        *relevance.get_or_insert_with(|| match pruning_model {
            Some(probe) => !finite_bounds.excludes(probe, fa, fb, dec),
            None => true,
        })
    };
    match pair {
        SurfacePair::Disjoint => {}
        SurfacePair::Coincident { .. } => {
            *overlapping_coincidence |= relevant(dec);
            for face in [fa, fb] {
                let surface = model.face(face).surface;
                for curve in face_edges(model, face) {
                    let relevant = relevant(dec);
                    let index = add_curve(model, pool, curve, dec)?;
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
                        canonical_point(model, points, point, dec, "boolean.imprint-ambiguous")?;
                    }
                    Component::Curve {
                        curve,
                        singular: nodes,
                        ..
                    } => {
                        let relevant = relevant(dec);
                        let index = add_curve(model, pool, curve, dec)?;
                        active.resize(pool.len(), false);
                        active[index] |= relevant;
                        supports.resize(pool.len(), BTreeSet::new());
                        supports[index].extend([sa, sb]);
                        for node in nodes {
                            let point = canonical_point(
                                model,
                                points,
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
    Ok(())
}

/// The carriers of one pass, as the carrier-pair pass reads them.
#[derive(Clone, Copy)]
struct Carriers<'a> {
    pool: &'a [CurveId],
    active: &'a [bool],
    supports: &'a [BTreeSet<SurfaceId>],
    complete: &'a [Complete],
}

impl Carriers<'_> {
    /// Whether `carrier_pair` returns for this pair before any decision.
    fn settled_without_proof(&self, (i, j): (usize, usize)) -> bool {
        let covered = |c: usize, d: usize| {
            self.supports[c]
                .iter()
                .any(|s| self.complete[d].contains_key(s))
        };
        match (self.active[i], self.active[j]) {
            (true, true) => covered(i, j) || covered(j, i),
            (true, false) => covered(j, i),
            (false, false) => true,
            (false, true) => false,
        }
    }
}

/// Canonicalizes where carriers `i < j` cross, unless their supports already
/// account for every crossing.
fn carrier_pair(
    model: &mut Model,
    (i, j): (usize, usize),
    carriers: Carriers<'_>,
    points: &mut Vec<PointId>,
    order: &mut Order,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let Carriers {
        pool,
        active,
        supports,
        complete,
    } = carriers;
    // Two carriers meet only on each other's supports. When one
    // carrier's complete meeting with the other's support is already
    // canonical, the pair can add no point; equality is never inferred.
    let covered = |c: usize, d: usize| supports[c].iter().any(|s| complete[d].contains_key(s));
    let skip = match (active[i], active[j]) {
        (true, true) => covered(i, j) || covered(j, i),
        // A skipped carrier has no complete meetings of its own. A
        // later one's crossings were all found by the earlier active
        // carrier or before it, so its hits could not come first.
        (true, false) => covered(j, i),
        // An earlier skipped carrier could find a crossing first. Its
        // crossings with a covering active carrier are that carrier's
        // canonical hits on the skipped carrier's supports.
        (false, true) if covered(i, j) => {
            for surface in &supports[i] {
                for &point in complete[j].get(surface).into_iter().flatten() {
                    if order.finder.get(&point).is_some_and(|&found| i >= found) {
                        continue;
                    }
                    if coincident(
                        on_curve(model.geom(), pool[i], point, dec),
                        "boolean.imprint-ambiguous",
                    )? {
                        order.crossers.entry(point).or_default().insert(i);
                    }
                }
            }
            true
        }
        // Crossings of two skipped carriers lie on no kept piece.
        (false, false) => true,
        (false, true) => false,
    };
    if skip {
        return Ok(());
    }
    let crate::curve_pair::Meetings::Points(meetings) =
        crate::curve_pair::meetings(model, pool[i], pool[j], dec)?
    else {
        return Ok(());
    };
    for point in meetings {
        let canonical = canonical_point(model, points, point, dec, "boolean.imprint-ambiguous")?;
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
    Ok(())
}

/// The section pieces of one carrier between the canonical points on it.
fn carrier_sections(
    model: &mut Model,
    curve: CurveId,
    (pool, surfaces, points): (&[CurveId], &[SurfaceId], &[PointId]),
    context: &Context<'_>,
    regions: &mut RegionCache,
    order: &mut Order,
    dec: &mut Decider,
) -> Result<Vec<SectionEdge>, Refusal> {
    let mut cuts = Vec::new();
    for &point in points {
        if coincident(
            on_curve(model.geom(), curve, point, dec),
            "boolean.imprint-ambiguous",
        )? {
            cuts.push(point);
        }
    }
    sort_points(model, curve, &mut cuts, dec)?;
    let closed = curve_is_closed(model.geom(), curve, dec)?;
    let mut pieces = carrier_pieces(model, curve, &cuts, closed, context, regions, dec)?;
    // A closed carrier starts at the cut the pass found first. Restore the
    // complete pass's choice when a skipped carrier could change it.
    if closed && !pieces.is_empty() {
        let start = order.first(model, pool, surfaces, &cuts, dec)?;
        if start != 0 {
            cuts.rotate_left(start);
            pieces = carrier_pieces(model, curve, &cuts, closed, context, regions, dec)?;
        }
    }
    Ok(pieces)
}

/// Whether `point` lies on a section piece.
fn on_section(
    model: &Model,
    edges: &[SectionEdge],
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for edge in edges {
        if on_arc(
            model,
            edge.curve,
            edge.start,
            edge.end,
            point,
            dec,
            "boolean.imprint-ambiguous",
        )? {
            return Ok(true);
        }
    }
    Ok(false)
}
