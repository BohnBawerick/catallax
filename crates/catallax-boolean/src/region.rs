//! Face membership by certified boundary homotopies in rational charts.
//! Spherical winding is relative to a fixed rational reference point, so
//! changing a projection pole cannot change a cell's winding signature.
//! Circular torus cuts use seam-independent cells of the two angular circles.

use crate::geometry::{natural_ends, on_arc, product, refusal, sign};
use crate::project::{
    Box2, Point2, Projection, align, ball, common_point, separated, shift, winding,
};
use catallax_intersect::{CurveTrace, TracePiece, trace};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};
use catallax_shapes::{Point3, PointId, Surface};
use catallax_topology::{EdgeId, FaceId, Loop, LoopId, Model, Sense};
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

#[derive(PartialEq, Eq, PartialOrd, Ord)]
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

#[derive(Default)]
struct TraceCache {
    curves: HashMap<EdgeId, CurveTrace>,
    boxes: TraceEnclosures,
}

/// Proofs prepared during a batch of queries against one append-only model.
/// The owner must not reuse this workspace with another model or decider.
#[derive(Default)]
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

struct Sampling<'a, const CACHE: bool> {
    piece: &'a TracePiece,
    edge: EdgeId,
    index: usize,
    enclosures: &'a mut TraceEnclosures,
    projection: &'a Projection,
    chart: Option<(catallax_shapes::SurfaceId, i64)>,
    queries: &'a [Box2],
    bits: u32,
    remaining: &'a mut usize,
    dec: &'a mut Decider,
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
                .project_certified(&self.piece.evaluate(&parameter, self.dec)?, self.dec);
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
        match self.enclosures.entry(key) {
            std::collections::btree_map::Entry::Occupied(entry) => Ok(entry.get().clone()),
            std::collections::btree_map::Entry::Vacant(entry) => {
                let value = self
                    .projection
                    .project_certified(&self.piece.evaluate(&parameter, self.dec)?, self.dec)?;
                if room {
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
            let traced = match work.traces.curves.entry(coedge.edge) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    let (from, to) = natural_ends(model, coedge.edge);
                    let answer = trace(model.geom(), edge.curve, from, to, dec)?;
                    if !answer.marks.is_empty() {
                        return Err(refusal(
                            "boolean.classify-ambiguous",
                            "a marked trace cannot certify face membership",
                        ));
                    }
                    entry.insert(answer.value)
                }
            };
            let mut edge_segments = Vec::new();
            for (index, piece) in traced.pieces().iter().enumerate() {
                Sampling::<CACHE> {
                    piece,
                    edge: coedge.edge,
                    index,
                    enclosures: &mut work.traces.boxes,
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
        let mut work = SamplingWork {
            traces: &mut self.traces,
            remaining: MAX_SAMPLE_CELLS,
        };
        let limit = dec.policy().max_precision_bits.max(1);
        let mut bits = limit.min(128);
        loop {
            let width = Rational::one().mul_pow2(-i64::from(bits / 2));
            let coordinates = match self.coordinates.entry((point, bits)) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    entry.insert(model.geom().enclose_point(point, &width, dec)?)
                }
            };
            let choices = if matches!(surface, Surface::Sphere(_)) {
                crate::geometry::DIRECTION_ATTEMPTS
            } else {
                1
            };
            for choice in 0..choices {
                let projection = Projection::new(surface, i64::from(choice))?;
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
        if let Some(location) = crate::fast_region::locate(model, face, point, dec)? {
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
        if PROBE
            && !matches!(model.geom().surface(f.surface), Surface::Torus(_))
            && let Some(location) = crate::probed_region::locate(model, face, point, dec)?
        {
            return Ok(location);
        }
        let surface = model.geom().surface(f.surface).clone();
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

pub(crate) fn signatures(
    model: &Model,
    surface: &Surface,
    loop_groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Vec<i64>, Refusal> {
    RegionCache::default().signatures(model, surface, loop_groups, point, dec)
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
}
