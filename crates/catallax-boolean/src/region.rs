//! Face membership by certified boundary homotopies in rational charts.
//! Spherical winding is relative to a fixed rational reference point, so
//! changing a projection pole cannot change a cell's winding signature.
//! Circular torus cuts use seam-independent cells of the two angular circles.

use crate::geometry::{natural_ends, on_arc, product, refusal, sign};
use crate::project::{
    Box2, Point2, Projection, align, ball, common_point, separated, shift, winding,
};
use catallax_intersect::{CurveTrace, TracePiece, trace};
use catallax_numbers::{Ball, Decider, FactKind, Rational, Refusal, Sign};
use catallax_shapes::{GeomKey, Point3, PointId, SplineSupport, Surface};
use catallax_topology::{
    EdgeId, FaceId, Loop, LoopId, Model, Sense, SplinePrism, envelope_retains_chart, pipe_rows,
    saddle_rows,
};
use std::collections::{HashMap, hash_map::Entry};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Location {
    Inside,
    Outside,
    Boundary,
}
struct Segment {
    entire: Box2,
    start: Box2,
    end: Box2,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct ProjectionKey {
    surface: catallax_shapes::SurfaceId,
    choice: i64,
    edge: EdgeId,
    piece: usize,
    bits: u32,
    from: Rational,
    to: Rational,
}
type TraceEnclosures = std::collections::BTreeMap<ProjectionKey, Option<Box2>>;
const MAX_TRACE_BOXES: usize = 4096;

#[derive(Clone, Default)]
struct TraceCache {
    curves: HashMap<EdgeId, CurveTrace>,
    boxes: TraceEnclosures,
    // How many values this cache and the coordinates beside it did not hold
    // and so computed, while a count is kept.
    missed: Option<usize>,
    // Work-span observation only: names this workspace's entries, and gives
    // each projection key a dense number.
    #[cfg(feature = "work-span")]
    span: catallax_numbers::work_span::Instance,
    #[cfg(feature = "work-span")]
    span_ids: std::collections::BTreeMap<ProjectionKey, u64>,
}

/// Proofs prepared during a batch of queries against one append-only model.
/// The owner must not reuse this workspace with another model or decider.
#[derive(Clone, Default)]
pub(crate) struct RegionCache {
    traces: TraceCache,
    coordinates: HashMap<(PointId, u32), [Ball; 3]>,
}

// Shared across pieces, projection choices and precision retries. Precision
// controls arithmetic enclosures, never the amount of subdivision work.
const MAX_SAMPLE_CELLS: usize = 4096;

struct SamplingWork<'a> {
    traces: &'a mut TraceCache,
    remaining: usize,
}

// Counts a value a cache did not hold, while a count is kept.
fn count(missed: &mut Option<usize>) {
    if let Some(missed) = missed {
        *missed += 1;
    }
}

/// What another thread hands back for a loop item that adds to the model and
/// uses the loop's region cache: the item as [`crate::geometry::plan`] ran
/// it, and whether its copy of the cache held everything the item asked for.
pub(crate) struct Ahead<T> {
    planned: crate::geometry::Planned<T>,
    missed: bool,
}

impl<T, E> Ahead<Result<T, E>>
where
    Result<T, E>: catallax_topology::Relocate,
{
    /// Whether the writer may take the item, with `model` as it is now.
    ///
    /// Not when the item refused: a refusal is the writer's to give, from
    /// its own run. Not when the item would have read something else on
    /// `model`. And not when the item's copy of the region cache missed a
    /// value. The writer's cache holds everything the copy began with, so an
    /// item that only found values finds them here too and charges the same.
    /// An item that missed one may find it here, filled by an earlier item,
    /// and would then have charged less.
    pub(crate) fn holds(&self, model: &Model) -> bool {
        self.planned.value().is_ok() && !self.missed && self.planned.holds(model)
    }

    /// Does to `model` what the item did to its copy and returns the item's
    /// value. Only for an item that [`Ahead::holds`] on `model`. The writer's
    /// region cache needs nothing: the item computed nothing it lacks.
    pub(crate) fn commit(self, model: &mut Model) -> Result<T, E> {
        self.planned.commit(model)
    }
}

struct Sampling<'a, const CACHE: bool> {
    piece: &'a TracePiece,
    edge: EdgeId,
    index: usize,
    enclosures: &'a mut TraceEnclosures,
    missed: &'a mut Option<usize>,
    projection: &'a Projection,
    chart: Option<(catallax_shapes::SurfaceId, i64)>,
    queries: &'a [Box2],
    bits: u32,
    remaining: &'a mut usize,
    dec: &'a mut Decider,
    #[cfg(feature = "work-span")]
    span: u64,
    #[cfg(feature = "work-span")]
    span_ids: &'a mut std::collections::BTreeMap<ProjectionKey, u64>,
}

impl<const CACHE: bool> Sampling<'_, CACHE> {
    fn project(
        &mut self,
        parameter: Ball,
        from: &Rational,
        to: &Rational,
    ) -> Result<Option<Box2>, Refusal> {
        if !CACHE {
            return self
                .projection
                .project_certified(&self.piece.enclose(&parameter, self.dec)?, self.dec);
        }
        self.dec.check_deadline()?;
        let room = self.enclosures.len() < MAX_TRACE_BOXES;
        let (surface, choice) = self
            .chart
            .expect("prepared projections retain their chart identity");
        let key = ProjectionKey {
            surface,
            choice,
            edge: self.edge,
            piece: self.index,
            bits: self.bits,
            from: from.clone(),
            to: to.clone(),
        };
        #[cfg(feature = "work-span")]
        let span_key = {
            let next = self.span_ids.len() as u64;
            *self.span_ids.entry(key.clone()).or_insert(next)
        };
        #[cfg(feature = "work-span")]
        catallax_numbers::work_span::read_at("trace projection cache", self.span, span_key);
        match self.enclosures.entry(key) {
            std::collections::btree_map::Entry::Occupied(entry) => Ok(entry.get().clone()),
            std::collections::btree_map::Entry::Vacant(entry) => {
                count(self.missed);
                let value = self
                    .projection
                    .project_certified(&self.piece.enclose(&parameter, self.dec)?, self.dec)?;
                if room {
                    #[cfg(feature = "work-span")]
                    catallax_numbers::work_span::write_at(
                        "trace projection cache",
                        self.span,
                        span_key,
                    );
                    Ok(entry.insert(value).clone())
                } else {
                    Ok(value)
                }
            }
        }
    }
    fn sample(
        &mut self,
        from: Rational,
        to: Rational,
        depth: u32,
        out: &mut Vec<Segment>,
    ) -> Result<(), Refusal> {
        let mut pending = vec![(from, to, depth)];
        while let Some((from, to, depth)) = pending.pop() {
            if let Some(mid) = self.cell(&from, &to, depth, out)? {
                // Pop the left half first, preserving the recursive traversal.
                pending.push((mid.clone(), to, depth + 1));
                pending.push((from, mid, depth + 1));
            }
        }
        Ok(())
    }

    fn cell(
        &mut self,
        from: &Rational,
        to: &Rational,
        depth: u32,
        out: &mut Vec<Segment>,
    ) -> Result<Option<Rational>, Refusal> {
        if *self.remaining == 0 {
            return Err(refusal(
                "boolean.sampling-budget",
                "face membership exhausted its deterministic sampling cell allowance",
            ));
        }
        *self.remaining -= 1;
        let projection = self.projection;
        let queries = self.queries;
        let bits = self.bits;
        let parameter = ball(from, to, bits);
        let projected = self.project(parameter, from, to)?;
        let mid = (from + to).mul_pow2(-1);
        if let Some(mut entire) = projected {
            let start = self.project(Ball::from_rational(from, bits), from, from)?;
            let end = self.project(Ball::from_rational(to, bits), to, to)?;
            if let (Some(mut start), Some(mut end)) = (start, end) {
                if projection.periodic() {
                    align(&mut start, &entire);
                    align(&mut end, &entire);
                }
                for i in 0..2 {
                    entire[i] = entire[i].hull(&start[i]).hull(&end[i]);
                }
                let mut clear = true;
                for query in queries {
                    let disjoint = if projection.periodic() {
                        let mut q = query.clone();
                        align(&mut q, &entire);
                        separated(&entire, &q, self.dec)?
                    } else {
                        separated(&entire, query, self.dec)?
                    };
                    if !disjoint {
                        clear = false;
                        break;
                    }
                }
                if clear && depth >= 2 {
                    out.push(Segment { entire, start, end });
                    return Ok(None);
                }
                let Some(mut middle) = self.project(Ball::from_rational(&mid, bits), &mid, &mid)?
                else {
                    return Err(refusal(
                        "boolean.classify-ambiguous",
                        "a trace midpoint has no separated chart at this precision",
                    ));
                };
                if projection.periodic() {
                    align(&mut middle, &entire);
                }
                if (0..2)
                    .all(|i| entire[i].lo() == middle[i].lo() && entire[i].hi() == middle[i].hi())
                {
                    return Err(refusal(
                        "boolean.classify-ambiguous",
                        "trace variation is below the coefficient enclosure precision",
                    ));
                }
            }
        } else if self
            .project(Ball::from_rational(&mid, bits), &mid, &mid)?
            .is_none()
        {
            return Err(refusal(
                "boolean.classify-ambiguous",
                "a trace midpoint has no separated chart at this precision",
            ));
        }
        if depth >= bits {
            return Err(refusal(
                "boolean.classify-ambiguous",
                "trace boxes do not separate the boundary from the query at this precision",
            ));
        }
        Ok(Some(mid))
    }
}

fn polygon(segments: &mut [Segment], periodic: bool) -> Result<Vec<(Point2, Point2)>, Refusal> {
    if segments.is_empty() {
        return Ok(Vec::new());
    }
    if periodic {
        for i in 1..segments.len() {
            let reference = segments[i - 1].end.clone();
            let amount = align(&mut segments[i].start, &reference);
            shift(&mut segments[i].entire, &amount);
            shift(&mut segments[i].end, &amount);
        }
    }
    let n = segments.len();
    let mut closing = segments[n - 1].end.clone();
    let closing_shift = if periodic {
        align(&mut closing, &segments[0].start)
    } else {
        Rational::zero()
    };
    let mut last_box = segments[n - 1].entire.clone();
    if periodic {
        shift(&mut last_box, &closing_shift);
    }
    let first = common_point(&[&closing, &segments[0].start, &last_box, &segments[0].entire])
        .ok_or_else(|| {
            refusal(
                "boolean.classify-ambiguous",
                "trace endpoints have no common certified chart box",
            )
        })?;
    let mut vertices = vec![first.clone()];
    for i in 1..n {
        vertices.push(
            common_point(&[
                &segments[i - 1].end,
                &segments[i].start,
                &segments[i - 1].entire,
                &segments[i].entire,
            ])
            .ok_or_else(|| {
                refusal(
                    "boolean.classify-ambiguous",
                    "adjacent trace boxes did not certify a common endpoint",
                )
            })?,
        );
    }
    let mut last = first;
    last[0] = &last[0] - &closing_shift;
    vertices.push(last);
    Ok(vertices
        .windows(2)
        .map(|w| (w[0].clone(), w[1].clone()))
        .collect())
}

fn build<const CACHE: bool>(
    model: &Model,
    loops: &[LoopId],
    projection: (&Projection, Option<(catallax_shapes::SurfaceId, i64)>),
    queries: &[Box2],
    bits: u32,
    work: &mut SamplingWork<'_>,
    dec: &mut Decider,
) -> Result<Vec<Vec<(Point2, Point2)>>, Refusal> {
    let (projection, chart) = projection;
    let mut polygons = Vec::new();
    for &loop_id in loops {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            continue;
        };
        let mut segments = Vec::new();
        for &coedge_id in coedges {
            let coedge = model.coedge(coedge_id);
            let edge = model.edge(coedge.edge);
            #[cfg(feature = "work-span")]
            catallax_numbers::work_span::read_at(
                "edge trace cache",
                work.traces.span.id(),
                u64::from(coedge.edge.index()),
            );
            let traced = match work.traces.curves.entry(coedge.edge) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    count(&mut work.traces.missed);
                    let (from, to) = natural_ends(model, coedge.edge);
                    let answer = trace(model.geom(), edge.curve, from, to, dec)?;
                    if !answer.marks.is_empty() {
                        return Err(refusal(
                            "boolean.classify-ambiguous",
                            "a marked trace cannot certify face membership",
                        ));
                    }
                    #[cfg(feature = "work-span")]
                    catallax_numbers::work_span::write_at(
                        "edge trace cache",
                        work.traces.span.id(),
                        u64::from(coedge.edge.index()),
                    );
                    entry.insert(answer.value)
                }
            };
            let mut edge_segments = Vec::new();
            for (index, piece) in traced.pieces().iter().enumerate() {
                Sampling::<CACHE> {
                    #[cfg(feature = "work-span")]
                    span: work.traces.span.id(),
                    #[cfg(feature = "work-span")]
                    span_ids: &mut work.traces.span_ids,
                    piece,
                    edge: coedge.edge,
                    index,
                    enclosures: &mut work.traces.boxes,
                    missed: &mut work.traces.missed,
                    projection,
                    chart,
                    queries,
                    bits,
                    remaining: &mut work.remaining,
                    dec,
                }
                .sample(
                    Rational::zero(),
                    Rational::one(),
                    0,
                    &mut edge_segments,
                )?;
            }
            if product(coedge.sense, edge.sense) == Sense::Reversed {
                edge_segments.reverse();
                for segment in &mut edge_segments {
                    std::mem::swap(&mut segment.start, &mut segment.end);
                }
            }
            segments.extend(edge_segments);
        }
        polygons.push(polygon(&mut segments, projection.periodic())?);
    }
    Ok(polygons)
}

fn off_carriers(
    model: &Model,
    loops: &[LoopId],
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let coordinates = [point.x.clone(), point.y.clone(), point.z.clone()];
    for &l in loops {
        if let Loop::Coedges(cs) = model.loop_(l) {
            for &c in cs {
                let curve = model.edge(model.coedge(c).edge).curve;
                let mut apart = false;
                for equation in model.geom().curve(curve).equations(model.geom()).ok_or_else(|| catallax_numbers::Refusal::new("boolean.parametric-query", "curve equation unavailable; this boolean path needs a finite parametric solver"))? {
                    if sign(
                        &equation.eval_real(&coordinates),
                        dec,
                        "boolean.classify-ambiguous",
                    )? != Sign::Zero
                    {
                        apart = true;
                        break;
                    }
                }
                if !apart {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

fn reference(
    model: &Model,
    surface: &Surface,
    loops: &[LoopId],
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    if !matches!(surface, Surface::Sphere(_)) {
        return Ok(None);
    }
    for k in 0..crate::geometry::DIRECTION_ATTEMPTS {
        let projection = Projection::new(surface, i64::from(k))?;
        let p = projection
            .pole_point()
            .expect("sphere projection has a pole");
        if off_carriers(model, loops, &p, dec)? {
            return Ok(Some(p));
        }
    }
    Err(refusal(
        "boolean.classify-ambiguous",
        "no rational reference point separated from the boundary carriers",
    ))
}

impl RegionCache {
    /// Runs a loop item on a thread that is ahead of the writer, on a copy of
    /// `base` and with a copy of this cache, which is the loop's cache as it
    /// was when the loop was handed over.
    pub(crate) fn ahead<T>(
        &self,
        base: &std::sync::Arc<Model>,
        work: impl FnOnce(&mut Model, &mut RegionCache) -> T,
    ) -> Ahead<T> {
        let mut regions = self.clone();
        regions.traces.missed = Some(0);
        let planned = crate::geometry::plan(base, |model| work(model, &mut regions));
        Ahead {
            planned,
            missed: regions.traces.missed != Some(0),
        }
    }

    pub(crate) fn signatures(
        &mut self,
        model: &Model,
        surface: &Surface,
        loop_groups: &[Vec<LoopId>],
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Vec<i64>, Refusal> {
        self.signatures_with::<false>(model, surface, loop_groups, point, dec, None)
    }

    pub(crate) fn signatures_prepared(
        &mut self,
        model: &Model,
        surface: catallax_shapes::SurfaceId,
        loop_groups: &[Vec<LoopId>],
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Vec<i64>, Refusal> {
        self.signatures_with::<true>(
            model,
            model.geom().surface(surface),
            loop_groups,
            point,
            dec,
            Some(surface),
        )
    }

    fn signatures_with<const CACHE: bool>(
        &mut self,
        model: &Model,
        surface: &Surface,
        loop_groups: &[Vec<LoopId>],
        point: PointId,
        dec: &mut Decider,
        surface_id: Option<catallax_shapes::SurfaceId>,
    ) -> Result<Vec<i64>, Refusal> {
        dec.check_deadline()?;
        let surface = &*crate::geometry::chart(surface);
        if let Surface::Torus(torus) = surface
            && let Some(signature) =
                crate::torus_region::signature(model, torus, loop_groups, point, dec)?
        {
            return Ok(signature);
        }
        if let Some(direct) =
            crate::fast_region::signatures(model, surface, loop_groups, point, dec)?
        {
            return Ok(direct);
        }
        let loops: Vec<_> = loop_groups.iter().flatten().copied().collect();
        let reference = reference(model, surface, &loops, dec)?;
        #[cfg(feature = "work-span")]
        let coordinates_span = self.traces.span.id();
        let mut work = SamplingWork {
            traces: &mut self.traces,
            remaining: MAX_SAMPLE_CELLS,
        };
        let limit = dec.policy().max_precision_bits.max(1);
        let mut bits = limit.min(128);
        loop {
            let width = Rational::one().mul_pow2(-i64::from(bits / 2));
            #[cfg(feature = "work-span")]
            catallax_numbers::work_span::read_at(
                "point coordinate cache",
                coordinates_span,
                (u64::from(point.index()) << 32) | u64::from(bits),
            );
            let coordinates = match self.coordinates.entry((point, bits)) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    count(&mut work.traces.missed);
                    #[cfg(feature = "work-span")]
                    let _written = catallax_numbers::work_span::write_on_drop_at(
                        "point coordinate cache",
                        coordinates_span,
                        (u64::from(point.index()) << 32) | u64::from(bits),
                    );
                    // Exact coordinates may already meet the width at 64 bits.
                    // Promote their working precision too, or chart arithmetic
                    // keeps rounding at 64 bits through every precision retry.
                    entry.insert(
                        model
                            .geom()
                            .enclose_point(point, &width, dec)?
                            .map(|value| value.with_precision(bits)),
                    )
                }
            };
            let choices = if matches!(surface, Surface::Sphere(_)) {
                crate::geometry::DIRECTION_ATTEMPTS
            } else {
                1
            };
            for choice in 0..choices {
                let projection = Projection::certified(surface, i64::from(choice), dec)?;
                if let Some(pole) = projection.pole_point()
                    && !off_carriers(model, &loops, &pole, dec)?
                {
                    continue;
                }
                let Some(query) = projection.project_certified(coordinates, dec)? else {
                    continue;
                };
                let mut queries = vec![query];
                if let Some(p) = &reference {
                    let coordinates = [
                        dec.enclose_at(&p.x, bits)?,
                        dec.enclose_at(&p.y, bits)?,
                        dec.enclose_at(&p.z, bits)?,
                    ];
                    let Some(projected) = projection.project_certified(&coordinates, dec)? else {
                        continue;
                    };
                    queries.push(projected);
                }
                let polygons = match build::<CACHE>(
                    model,
                    &loops,
                    (&projection, surface_id.map(|id| (id, i64::from(choice)))),
                    &queries,
                    bits,
                    &mut work,
                    dec,
                ) {
                    Ok(p) => p,
                    Err(e) if e.code == "boolean.classify-ambiguous" => break,
                    Err(e) => return Err(e),
                };
                let query = queries[0].clone().map(|b| b.mid());
                let reference = queries.get(1).map(|q| q.clone().map(|b| b.mid()));
                let mut offset = 0;
                let mut result = Vec::new();
                for group in loop_groups {
                    let end = offset
                        + group
                            .iter()
                            .filter(|l| matches!(model.loop_(**l), Loop::Coedges(_)))
                            .count();
                    let mut value = 0;
                    for polygon in &polygons[offset..end] {
                        value += winding(polygon, &query, projection.periodic(), dec)?;
                        if let Some(reference) = &reference {
                            value -= winding(polygon, reference, false, dec)?;
                        }
                    }
                    result.push(value);
                    offset = end;
                }
                return Ok(result);
            }
            if bits == limit {
                break;
            }
            bits = bits.saturating_mul(2).min(limit);
        }
        Err(refusal(
            "boolean.classify-ambiguous",
            "no certified chart separated the boundary from the query",
        ))
    }

    pub(crate) fn locate(
        &mut self,
        model: &mut Model,
        face: FaceId,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        self.locate_with::<true>(model, face, point, dec, None)
    }

    pub(crate) fn locate_prepared(
        &mut self,
        model: &mut Model,
        face: FaceId,
        point: PointId,
        dec: &mut Decider,
        targets: &mut std::collections::BTreeMap<FaceId, Vec<i64>>,
    ) -> Result<Location, Refusal> {
        self.locate_with::<false>(model, face, point, dec, Some(targets))
    }

    fn locate_with<const PROBE: bool>(
        &mut self,
        model: &mut Model,
        face: FaceId,
        point: PointId,
        dec: &mut Decider,
        targets: Option<&mut std::collections::BTreeMap<FaceId, Vec<i64>>>,
    ) -> Result<Location, Refusal> {
        dec.check_deadline()?;
        let f = model.face(face).clone();
        if f.loops.is_empty() {
            return Ok(Location::Inside);
        }
        if let Some(location) =
            crate::spline_fillet::locate(model, face, &model.geom().position(point), dec)?
        {
            return Ok(location);
        }
        if let Some(location) = whole_chart(model, face, point, dec)? {
            return Ok(location);
        }
        if let Some(location) = crate::fast_region::locate(model, face, point, dec)? {
            return Ok(location);
        }
        if let Some(location) = band_caps(model, face, point, dec)? {
            return Ok(location);
        }
        for &l in &f.loops {
            match model.loop_(l) {
                Loop::Vertex(v) => {
                    if crate::geometry::coincident(
                        model.geom().same_point(model.vertex(*v).point, point, dec),
                        "boolean.classify-ambiguous",
                    )? {
                        return Ok(Location::Boundary);
                    }
                }
                Loop::Coedges(cs) => {
                    for &c in cs {
                        let e = model.coedge(c).edge;
                        let (a, b) = natural_ends(model, e);
                        if on_arc(
                            model,
                            model.edge(e).curve,
                            a,
                            b,
                            point,
                            dec,
                            "boolean.classify-ambiguous",
                        )? {
                            return Ok(Location::Boundary);
                        }
                    }
                }
            }
        }
        // A proved injective spline projection supplies the complete winding
        // query directly. Auxiliary sections need no separate curve-equivalence
        // solver when this face already has such a chart.
        let projected_spline =
            PROBE && Projection::spline(model.geom().surface(f.surface), dec)?.is_some();
        if PROBE
            && !projected_spline
            && !matches!(
                model.geom().surface(f.surface),
                Surface::Torus(_) | Surface::TorusSheet(_)
            )
            && let Some(location) = crate::probed_region::locate(model, face, point, dec)?
        {
            return Ok(location);
        }
        let surface = crate::geometry::chart(model.geom().surface(f.surface)).into_owned();
        let groups = vec![f.loops.clone()];
        let value = if PROBE {
            self.signatures(model, &surface, &groups, point, dec)?
        } else {
            self.signatures_prepared(model, f.surface, &groups, point, dec)?
        };
        if matches!(surface, Surface::Sphere(_) | Surface::Torus(_)) {
            let barriers: Vec<EdgeId> = f
                .loops
                .iter()
                .flat_map(|l| match model.loop_(*l) {
                    Loop::Coedges(cs) => cs.iter().map(|c| model.coedge(*c).edge).collect(),
                    Loop::Vertex(_) => Vec::new(),
                })
                .collect();
            let mut compute = || {
                crate::points::left_probe(model, face, &barriers, dec, |model, p, _, dec| {
                    Ok(Some(if PROBE {
                        self.signatures(model, &surface, &groups, p, dec)?
                    } else {
                        self.signatures_prepared(model, f.surface, &groups, p, dec)?
                    }))
                })
            };
            let fresh;
            let target = if let Some(targets) = targets {
                match targets.entry(face) {
                    std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                    std::collections::btree_map::Entry::Vacant(entry) => entry.insert(compute()?),
                }
            } else {
                fresh = compute()?;
                &fresh
            };
            if &value == target {
                return Ok(Location::Inside);
            }
            if matches!(surface, Surface::Torus(_)) || (value[0] - target[0]).abs() == 1 {
                return Ok(Location::Outside);
            }
        } else {
            let value = if f.sense == Sense::Same {
                value[0]
            } else {
                -value[0]
            };
            if value == 1 {
                return Ok(Location::Inside);
            }
            if value == 0 {
                return Ok(Location::Outside);
            }
        }
        Err(refusal(
            "boolean.classify-ambiguous",
            "boundary winding does not certify a unit face region",
        ))
    }
}

// A sphere face bounded only by whole sphere band rails. Each rail's angle
// about its band's turn axis rises along its parameter, so like a circle
// walked round its normal, a rail walked forward keeps the face in the cap
// that axis points to. The face holds a point exactly when the point is in the face's
// cap of each rail. `None` leaves other faces, and points that no rail places,
// to the general test.
fn band_caps(
    model: &Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let f = model.face(face);
    let Surface::Sphere(sphere) = model.geom().surface(f.surface) else {
        return Ok(None);
    };
    let query = model.geom().position(point);
    let mut inside = true;
    let mut any = false;
    for &l in &f.loops {
        let Loop::Coedges(coedges) = model.loop_(l) else {
            return Ok(None);
        };
        for &c in coedges {
            let coedge = model.coedge(c);
            let Some(above) =
                catallax_topology::band_cap(model, coedge.edge, &sphere.center, &query, dec)?
            else {
                return Ok(None);
            };
            let forward = ((model.edge(coedge.edge).sense == Sense::Same)
                == (coedge.sense == Sense::Same))
                == (f.sense == Sense::Same);
            inside &= above == forward;
            any = true;
        }
    }
    Ok(any.then_some(if inside {
        Location::Inside
    } else {
        Location::Outside
    }))
}

// A point recorded on a finite envelope lies in its closed chart. Only a face
// that keeps that whole chart turns this into membership; any other face stays
// on its trimmed-domain path, so support incidence never decides it.
fn whole_chart(
    model: &Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let surface = model.face(face).surface;
    let recorded = model.geom().facts().holds(
        FactKind::Incident,
        GeomKey::Point(point),
        GeomKey::Surface(surface),
    );
    if !recorded {
        return Ok(None);
    }
    let position = model.geom().position(point);
    let on_boundary = match model.geom().surface(surface) {
        Surface::BSpline(SplineSupport::RationalPatch(_)) => {
            if crate::spline_fillet::trimmed_patch(model, face) {
                return Ok(None);
            }
            let Some(checked) = catallax_topology::RationalPatchFace::of_face(model, face, dec)?
            else {
                return Ok(None);
            };
            return Ok(Some(match checked.locate(&position, dec)? {
                Sign::Negative => Location::Inside,
                Sign::Zero => Location::Boundary,
                Sign::Positive => Location::Outside,
            }));
        }
        Surface::BSpline(SplineSupport::BezierExtrusion(_)) => {
            let Some(band) = catallax_topology::BezierExtrusionBand::of_face(model, face, dec)?
            else {
                return Ok(None);
            };
            return Ok(Some(match band.locate(&position, dec)? {
                Sign::Negative => Location::Inside,
                Sign::Zero => Location::Boundary,
                Sign::Positive => Location::Outside,
            }));
        }
        Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
            if !envelope_retains_chart(model, face, dec)? {
                return Ok(None);
            }
            envelope.on_boundary(&position, dec)?
        }
        Surface::BSpline(SplineSupport::RollingPipe(pipe)) => {
            if pipe_rows(model, face, dec)?.is_none() {
                return Ok(None);
            }
            pipe.on_boundary(&position, dec)?
        }
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
            if saddle_rows(model, face).is_none() {
                return Ok(None);
            }
            saddle.on_boundary(&position, dec)?
        }
        Surface::BSpline(SplineSupport::LinearFillet(band)) => {
            if catallax_topology::LinearFilletBand::of_face(model, face, dec)?.is_none() {
                return Ok(None);
            }
            band.on_boundary(&position, dec)?
        }
        Surface::BSpline(SplineSupport::RollingMitre(band)) => {
            if catallax_topology::RollingMitreBand::of_face(model, face, dec)?.is_none() {
                return Ok(None);
            }
            band.on_boundary(&position, dec)?
        }
        Surface::BSpline(SplineSupport::RollingSphere(band)) => {
            if saddle_rows(model, face).is_none() {
                return Ok(None);
            }
            band.on_boundary(&position, dec)?
        }
        _ => return Ok(None),
    };
    Ok(Some(if on_boundary {
        Location::Boundary
    } else {
        Location::Inside
    }))
}

// A complete checked prism bounds each face by its shared rings alone, so the
// section and the axial interval decide closed membership with no boundary
// tracing. This answers for the one face asked about. It stays out of `locate`,
// which also serves whole-body parity: a body short of a cap must not borrow a
// face from a complete prism and count a crossing through it.
pub(crate) fn prism_face(
    model: &Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    if let Some(band) = catallax_topology::BezierExtrusionBand::of_face(model, face, dec)? {
        return Ok(Some(
            band.locate(&model.geom().position(point), dec)? != Sign::Positive,
        ));
    }
    if let Some(cell) = catallax_topology::DraftedCylinder::of_face(model, face, dec)? {
        let position = model.geom().position(point);
        return Ok(Some(
            cell.locate_on_face(face, &position, dec)? != Sign::Positive,
        ));
    }
    let Some(prism) = SplinePrism::of_face(model, face, dec)? else {
        return Ok(None);
    };
    let position = model.geom().position(point);
    Ok(Some(
        prism.locate_on_face(face, &position, dec)? != Sign::Positive,
    ))
}

pub(crate) fn locate(
    model: &mut Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    RegionCache::default().locate(model, face, point, dec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn periodic_polygon_keeps_exact_period_shifts_when_midpoints_round() {
        let r = |s| Rational::parse(s).unwrap();
        let vertices = [[r("1/3"), r("0")], [r("4/3"), r("1")], [r("5/3"), r("-1")]];
        for bits in [8, 64, 128] {
            let mut segments = Vec::new();
            for i in 0..3 {
                let period = if i == 0 { r("0") } else { r("-4") };
                let endpoint = |v: &[Rational; 2]| {
                    [
                        Ball::from_rational(&(&v[0] + &period), bits),
                        Ball::from_rational(&v[1], bits),
                    ]
                };
                let start = endpoint(&vertices[i]);
                let end = endpoint(&vertices[(i + 1) % 3]);
                let entire = std::array::from_fn(|axis| start[axis].hull(&end[axis]));
                segments.push(Segment { entire, start, end });
            }
            let polygon = polygon(&mut segments, true).unwrap();
            assert_eq!(
                winding(&polygon, &[r("1"), r("0")], false, &mut Decider::default()).unwrap(),
                -1,
            );
        }
    }

    #[test]
    fn prepared_carrier_queries_keep_distinct_windings_and_honor_cancellation() {
        use catallax_numbers::Real;
        use catallax_shapes::{PointDef, Pose, RVector3, Rotation};
        use catallax_topology::cylinder;

        let mut model = Model::new();
        let mut dec = Decider::default();
        let first = cylinder(
            &mut model,
            Rational::from_i64(3),
            Rational::from_i64(5),
            &Pose::identity(),
        )
        .unwrap();
        let carrier = model
            .faces_of_body(first)
            .into_iter()
            .find_map(|face| {
                let id = model.face(face).surface;
                matches!(model.geom().surface(id), Surface::Cylinder(_)).then_some(id)
            })
            .unwrap();
        let pose = Pose::new(
            Rotation::identity(),
            RVector3::new(Rational::from_i64(2), Rational::zero(), Rational::zero()).into(),
            &mut dec,
        )
        .unwrap();
        let second = cylinder(
            &mut model,
            Rational::from_i64(3),
            Rational::from_i64(5),
            &pose,
        )
        .unwrap();
        let result = crate::boolean(&mut model, first, second, crate::BoolOp::Subtract, &mut dec)
            .unwrap()
            .body;
        let face = model
            .faces_of_body(result)
            .into_iter()
            .find(|&face| model.face(face).surface == carrier)
            .unwrap();
        let groups = [model.face(face).loops.clone()];
        let surface = model.geom().surface(carrier).clone();
        let mut prepared = RegionCache::default();
        for ([x, y, z], expected) in [
            ([-3, 0, 2], 1),
            ([3, 0, 2], 0),
            ([-3, 0, -1], 0),
            ([0, 3, 2], 1),
        ] {
            let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
                Real::from_i64(x),
                Real::from_i64(y),
                Real::from_i64(z),
            )));
            assert_eq!(
                prepared
                    .signatures(&model, &surface, &groups, point, &mut dec)
                    .unwrap(),
                vec![expected],
                "query {x},{y},{z}",
            );
        }
        dec.set_deadline(Some(std::time::Instant::now()));
        let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
            Real::from_i64(-3),
            Real::zero(),
            Real::from_i64(2),
        )));
        assert_eq!(
            prepared
                .signatures(&model, &surface, &groups, point, &mut dec)
                .unwrap_err()
                .code,
            "numbers.deadline",
        );
    }

    // A face on a cylinder that another cylinder cut, five points on the
    // carrier to ask about, and a second point at each of the same places.
    // The face's winding needs traces of its edges.
    fn cut_cylinder() -> (
        Model,
        catallax_shapes::SurfaceId,
        Vec<Vec<LoopId>>,
        [Vec<PointId>; 2],
    ) {
        use catallax_numbers::Real;
        use catallax_shapes::{PointDef, Pose, RVector3, Rotation};
        use catallax_topology::cylinder;

        let mut model = Model::new();
        let mut dec = Decider::default();
        let radius = Rational::from_i64(3);
        let height = Rational::from_i64(5);
        let first = cylinder(
            &mut model,
            radius.clone(),
            height.clone(),
            &Pose::identity(),
        )
        .unwrap();
        let carrier = model
            .faces_of_body(first)
            .into_iter()
            .find_map(|face| {
                let id = model.face(face).surface;
                matches!(model.geom().surface(id), Surface::Cylinder(_)).then_some(id)
            })
            .unwrap();
        let pose = Pose::new(
            Rotation::identity(),
            RVector3::new(Rational::from_i64(2), Rational::zero(), Rational::zero()).into(),
            &mut dec,
        )
        .unwrap();
        let second = cylinder(&mut model, radius, height, &pose).unwrap();
        let result = crate::boolean(&mut model, first, second, crate::BoolOp::Subtract, &mut dec)
            .unwrap()
            .body;
        let face = model
            .faces_of_body(result)
            .into_iter()
            .find(|&face| model.face(face).surface == carrier)
            .unwrap();
        let groups = vec![model.face(face).loops.clone()];
        let points = [(); 2].map(|()| {
            [[-3, 0, 2], [3, 0, 2], [-3, 0, -1], [0, 3, 2], [0, -3, 1]]
                .map(|[x, y, z]| {
                    model.geom_mut().add_point(PointDef::Explicit(Point3::new(
                        Real::from_i64(x),
                        Real::from_i64(y),
                        Real::from_i64(z),
                    )))
                })
                .to_vec()
        });
        (model, carrier, groups, points)
    }

    // The answers with their whole refusals, the decisions made, the work
    // used and whether the allowance ran out, and what the cache holds.
    type Asked = (
        Vec<Result<Vec<i64>, Refusal>>,
        catallax_numbers::DecisionStats,
        (u64, bool),
        Held,
    );

    // What a region cache holds, in a fixed order: each traced edge with the
    // number of its pieces, every projected box, and every point enclosure.
    #[derive(PartialEq, Debug)]
    struct Held {
        traces: Vec<(EdgeId, usize)>,
        boxes: String,
        coordinates: Vec<String>,
    }

    fn held(regions: &RegionCache) -> Held {
        let mut traces: Vec<_> = regions
            .traces
            .curves
            .iter()
            .map(|(edge, trace)| (*edge, trace.pieces().len()))
            .collect();
        traces.sort_unstable();
        let mut coordinates: Vec<_> = regions.coordinates.iter().collect();
        coordinates.sort_unstable_by_key(|(key, _)| **key);
        Held {
            traces,
            boxes: format!("{:?}", regions.traces.boxes),
            coordinates: coordinates
                .into_iter()
                .map(|entry| format!("{entry:?}"))
                .collect(),
        }
    }

    // Asks for windings in three loops that share one cache, as the loops of
    // a Boolean share theirs. With `threads` above one, other threads answer
    // later items whole, each with the cache as it was when its loop began.
    //
    // The first loop asks about each point with the cache empty. The second
    // asks about the same points again and finds everything. The third asks
    // twice about each of the second points at the same places: the traces
    // are found, the point's enclosure is missed the first time and found
    // the second.
    fn ask(threads: usize, limit: u64) -> Asked {
        use catallax_numbers::{AheadSite, ValuesAhead};
        // The face is built on one thread, so that only the loops below take
        // anything from another.
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(1));
        let (mut model, carrier, groups, [points, twins]) = cut_cylinder();
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        catallax_numbers::set_run_ahead_wait(threads > 1);
        let twice: Vec<PointId> = twins.iter().flat_map(|&twin| [twin, twin]).collect();
        let ((answers, stats, cache), work) =
            catallax_numbers::with_operation_budget(limit, || {
                let mut dec = Decider::default();
                let mut regions = RegionCache::default();
                static SITE: AheadSite = AheadSite::new();
                let mut answers = Vec::new();
                for asked in [&points, &points, &twice] {
                    let mut ahead = ValuesAhead::new(&SITE, &dec, asked.len());
                    for (index, &point) in asked.iter().enumerate() {
                        let found = ahead.take_if(
                            index,
                            &mut dec,
                            || {
                                let (base, regions, groups, asked) = (
                                    std::sync::Arc::new(model.clone()),
                                    regions.clone(),
                                    groups.clone(),
                                    asked.clone(),
                                );
                                move |index: usize, dec: &mut Decider| {
                                    regions.ahead(&base, |model, regions| {
                                        regions.signatures_prepared(
                                            model,
                                            carrier,
                                            &groups,
                                            asked[index],
                                            dec,
                                        )
                                    })
                                }
                            },
                            |ahead| ahead.holds(&model),
                        );
                        let answer = match found {
                            Some(ahead) => ahead.commit(&mut model),
                            None => regions
                                .signatures_prepared(&model, carrier, &groups, point, &mut dec),
                        };
                        answers.push(answer);
                    }
                }
                (answers, dec.stats().clone(), held(&regions))
            });
        catallax_numbers::set_threads(None);
        catallax_numbers::set_run_ahead_wait(false);
        (answers, stats, (work.used, work.exhausted), cache)
    }

    // Another thread's item is taken only when its copy of the loop's cache
    // held all it asked for. An item that missed a value the writer holds by
    // then would have charged the writer for work it never does.
    #[test]
    fn an_item_that_missed_the_loops_cache_is_the_writers_to_run() {
        let one = ask(1, u64::MAX);
        let first = [vec![1], vec![0], vec![0], vec![1], vec![1]];
        assert_eq!(one.0[..5], first.clone().map(Ok::<_, Refusal>));
        assert_eq!(one.0[5..10], first.map(Ok::<_, Refusal>));
        assert!(!one.2.1, "the full allowance is enough");
        assert!(
            !one.3.traces.is_empty() && one.3.boxes.len() > 2,
            "traces and their boxes are kept: {:?}",
            one.3.traces
        );
        assert_eq!(one.3.coordinates.len(), 10, "one enclosure for each point");
        let before = catallax_numbers::run_ahead_adopted();
        assert_eq!(ask(2, u64::MAX), one, "two threads");
        // Only the second loop's later items found everything: four of them.
        // In the first loop every item missed the traces. In the third every
        // item missed its point's enclosure on the other thread, and for
        // every second one the writer held it by then.
        assert_eq!(catallax_numbers::run_ahead_adopted() - before, 4);
        assert_eq!(ask(4, u64::MAX), one, "four threads");
        // At allowances that run out part-way the same refusal falls on the
        // same item, with the same work, and the cache stops where it did.
        let total = one.2.0;
        for limit in [total - 1, total / 2, total / 3, total / 5] {
            let one = ask(1, limit);
            assert!(one.2.1, "limit {limit} of {total} runs out");
            assert!(one.0.iter().any(Result::is_err), "limit {limit} of {total}");
            assert_eq!(ask(3, limit), one, "limit {limit}");
        }
    }
}
