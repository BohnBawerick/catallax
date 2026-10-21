//! Face membership by certified boundary homotopies in rational charts.
//! Spherical winding is relative to a fixed rational reference point, so
//! changing a projection pole cannot change a cell's winding signature.

use crate::geometry::{natural_ends, on_arc, product, refusal, sign};
use crate::project::{
    Box2, Point2, Projection, align, ball, common_point, separated, shift, winding,
};
use catallax_intersect::{CurveTrace, TracePiece, trace};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{PointId, RPoint3, Surface};
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

struct Sampling<'a> {
    piece: &'a TracePiece,
    projection: &'a Projection,
    queries: &'a [Box2],
    bits: u32,
    dec: &'a mut Decider,
}

impl Sampling<'_> {
    fn sample(
        &mut self,
        from: Rational,
        to: Rational,
        depth: u32,
        out: &mut Vec<Segment>,
    ) -> Result<(), Refusal> {
        let piece = self.piece;
        let projection = self.projection;
        let queries = self.queries;
        let bits = self.bits;
        let dec = &mut *self.dec;
        let parameter = ball(&from, &to, bits);
        let projected = projection.project_certified(&piece.evaluate(&parameter, dec)?, dec)?;
        let mid = (&from + &to).mul_pow2(-1);
        if let Some(mut entire) = projected {
            let start = projection.project_certified(
                &piece.evaluate(&Ball::from_rational(&from, bits), dec)?,
                dec,
            )?;
            let end = projection
                .project_certified(&piece.evaluate(&Ball::from_rational(&to, bits), dec)?, dec)?;
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
                        separated(&entire, &q, dec)?
                    } else {
                        separated(&entire, query, dec)?
                    };
                    if !disjoint {
                        clear = false;
                        break;
                    }
                }
                if clear && depth >= 2 {
                    out.push(Segment { entire, start, end });
                    return Ok(());
                }
                let Some(mut middle) = projection.project_certified(
                    &piece.evaluate(&Ball::from_rational(&mid, bits), dec)?,
                    dec,
                )?
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
        } else if projection
            .project_certified(&piece.evaluate(&Ball::from_rational(&mid, bits), dec)?, dec)?
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
        self.sample(from, mid.clone(), depth + 1, out)?;
        self.sample(mid, to, depth + 1, out)
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

fn build(
    model: &Model,
    loops: &[LoopId],
    projection: &Projection,
    queries: &[Box2],
    bits: u32,
    traces: &mut HashMap<EdgeId, CurveTrace>,
    dec: &mut Decider,
) -> Result<Vec<Vec<(Point2, Point2)>>, Refusal> {
    let mut polygons = Vec::new();
    for &loop_id in loops {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            continue;
        };
        let mut segments = Vec::new();
        for &coedge_id in coedges {
            let coedge = model.coedge(coedge_id);
            let edge = model.edge(coedge.edge);
            let traced = match traces.entry(coedge.edge) {
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
            for piece in traced.pieces() {
                Sampling {
                    piece,
                    projection,
                    queries,
                    bits,
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
    point: &RPoint3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let coordinates = [point.x.clone(), point.y.clone(), point.z.clone()];
    for &l in loops {
        if let Loop::Coedges(cs) = model.loop_(l) {
            for &c in cs {
                let curve = model.edge(model.coedge(c).edge).curve;
                let mut apart = false;
                for equation in model.geom().curve(curve).equations(model.geom()) {
                    if sign(
                        &Real::from(equation.eval(&coordinates)),
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
) -> Result<Option<RPoint3>, Refusal> {
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

pub(crate) fn signatures(
    model: &Model,
    surface: &Surface,
    loop_groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Vec<i64>, Refusal> {
    let mut direct = Vec::with_capacity(loop_groups.len());
    for group in loop_groups {
        if let Some(value) = crate::fast_region::signature(model, surface, group, point, dec)? {
            direct.push(value);
        } else {
            break;
        }
    }
    if direct.len() == loop_groups.len() {
        return Ok(direct);
    }
    let loops: Vec<_> = loop_groups.iter().flatten().copied().collect();
    let reference = reference(model, surface, &loops, dec)?;
    let mut traces = HashMap::new();
    let limit = dec.policy().max_precision_bits.max(1);
    let mut bits = limit.min(128);
    loop {
        let width = Rational::one().mul_pow2(-i64::from(bits / 2));
        let coordinates = model.geom().enclose_point(point, &width, dec)?;
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
            let Some(query) = projection.project_certified(&coordinates, dec)? else {
                continue;
            };
            let mut queries = vec![query];
            if let Some(p) = &reference {
                let coordinates = [&p.x, &p.y, &p.z].map(|x| Ball::from_rational(x, bits));
                let Some(projected) = projection.project_certified(&coordinates, dec)? else {
                    continue;
                };
                queries.push(projected);
            }
            let polygons = match build(model, &loops, &projection, &queries, bits, &mut traces, dec)
            {
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
    model: &mut Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
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
    if let Some(location) = crate::probed_region::locate(model, face, point, dec)? {
        return Ok(location);
    }
    let surface = model.geom().surface(f.surface).clone();
    let groups = vec![f.loops.clone()];
    let value = signatures(model, &surface, &groups, point, dec)?[0];
    if matches!(surface, Surface::Sphere(_)) {
        let barriers: Vec<EdgeId> = f
            .loops
            .iter()
            .flat_map(|l| match model.loop_(*l) {
                Loop::Coedges(cs) => cs.iter().map(|c| model.coedge(*c).edge).collect(),
                Loop::Vertex(_) => Vec::new(),
            })
            .collect();
        let target = crate::points::left_probe(model, face, &barriers, dec, |model, p, _, dec| {
            Ok(Some(signatures(model, &surface, &groups, p, dec)?[0]))
        })?;
        if value == target {
            return Ok(Location::Inside);
        }
        if (value - target).abs() == 1 {
            return Ok(Location::Outside);
        }
    } else {
        let value = if f.sense == Sense::Same {
            value
        } else {
            -value
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
}
