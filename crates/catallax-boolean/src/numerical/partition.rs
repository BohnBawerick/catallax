//! Kept cone loops from finite native arcs and split source boundaries.

use super::*;
use catallax_topology::numerical::{NumericalEdge, NumericalFace, NumericalTrim};
use std::collections::{BTreeMap, BTreeSet};

pub(super) struct Input<'a> {
    pub key: &'a str,
    pub sense: Sense,
    pub tool: &'a Cylinder,
    pub intersection: NumericalNativeIntersection,
    pub rings: Vec<NumericalConeRing>,
    pub cone_quarters: Vec<Surface>,
    pub cylinder_source: NumericalCylinderSource,
    pub options: NumericalOptions,
    pub checks: Vec<NumericalCheck>,
}

struct Topology<'a> {
    key: &'a str,
    definition: NumericalDefinition,
    uv: Vec<[Rational; 2]>,
    native: BTreeMap<usize, BSplineCurve>,
    ring_links: Vec<NumericalRingLink>,
    boundary_edges: BTreeMap<(usize, Rational, usize, usize), usize>,
}
fn flip(sense: Sense) -> Sense {
    if sense == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}
fn reverse(mut boundary: Vec<NumericalTrim>) -> Vec<NumericalTrim> {
    boundary.reverse();
    for trim in &mut boundary {
        trim.sense = flip(trim.sense);
        trim.uv.swap(0, 1)
    }
    boundary
}

fn in_quarter(uv: &[Rational; 2], quarter: usize) -> bool {
    uv[0] >= q(quarter as i64) && uv[0] <= q(quarter as i64 + 1) && uv[1] >= q(0) && uv[1] <= q(1)
}
fn at_quarter(uv: &[Rational; 2], quarter: usize) -> [Rational; 2] {
    let mut result = uv.clone();
    if quarter == 3 && result[0] == q(0) {
        result[0] = q(4)
    }
    result
}

impl Input<'_> {
    pub fn build(self, dec: &mut Decider) -> Result<NumericalFacePartition, Refusal> {
        let Self {
            key,
            sense,
            tool,
            intersection,
            rings,
            cone_quarters,
            cylinder_source,
            options,
            mut checks,
        } = self;
        let source = &intersection.windows().sources()[0];
        let events::Events {
            vertices,
            uv,
            ends,
            mut records,
            repairs,
            checks: event_checks,
        } = events::collect(key, &intersection, &options.tolerance_mm, dec)?;
        checks.extend(event_checks);
        let mut topology = Topology {
            key,
            definition: NumericalDefinition {
                name: key.into(),
                vertices,
                edges: Vec::new(),
                faces: Vec::new(),
            },
            uv,
            native: BTreeMap::new(),
            ring_links: Vec::new(),
            boundary_edges: BTreeMap::new(),
        };
        let material = material::CylinderSide::new(tool, dec)?;
        material.covers(source, intersection.windows().cylinder_parameters(), dec)?;
        checks.push(NumericalCheck::Passed {
            name: "finite.tool-band-covers-cone-control-hull".into(),
            scope: key.into(),
        });
        let mut cylinder_trims = Vec::new();
        let mut arc_uses = BTreeMap::<usize, Vec<NumericalTrim>>::new();
        for (index, piece) in intersection.pieces().iter().enumerate() {
            if events::point_piece(piece, dec)? {
                continue;
            }
            let quarter = (0..4)
                .find(|i| piece.source_domains()[0][0][0] == q(*i as i64))
                .ok_or_else(|| bad("arc lost cone quarter"))?;
            let a = *ends
                .get(&(index, 0))
                .ok_or_else(|| bad("missing owned arc start"))?;
            let b = *ends
                .get(&(index, 1))
                .ok_or_else(|| bad("missing owned arc end"))?;
            let edge = topology.definition.edges.len();
            topology.definition.edges.push(NumericalEdge {
                key: format!("{key}/fit-{index}"),
                curve: Curve::BSpline(piece.curve().clone()),
                vertices: [a, b],
                parameters: piece.curve().domain().map(Clone::clone),
            });
            topology.native.insert(edge, piece.uv()[0].clone());
            let forward = material.outside_left(source, &piece.uv()[0], dec)?;
            let mut uv = [
                at_quarter(&topology.uv[a], quarter),
                at_quarter(&topology.uv[b], quarter),
            ];
            if !forward {
                uv.swap(0, 1)
            }
            arc_uses.entry(quarter).or_default().push(NumericalTrim {
                tag: catallax_topology::numerical::CopyTag::Independent,
                edge,
                sense: if forward {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                uv,
            });
            let cylinder_quarter = (0..4)
                .find(|i| piece.source_domains()[1][0][0] == q(*i as i64))
                .ok_or_else(|| bad("arc lost cylinder quarter"))?;
            cylinder_trims.push(NumericalCylinderTrim {
                piece: index,
                edge,
                quarter: cylinder_quarter,
                curve: piece.uv()[1].clone(),
            });
            checks.push(NumericalCheck::Passed {
                name: "finite.sampled-native-outside-left".into(),
                scope: format!("{key}/fit-{index}"),
            });
        }
        for (quarter, surface) in cone_quarters.into_iter().enumerate() {
            let mut uses = arc_uses.remove(&quarter).unwrap_or_default();
            for side in 0..4 {
                let fixed_axis = if side % 2 == 0 { 1 } else { 0 };
                let fixed = match side {
                    0 => q(0),
                    1 => q(quarter as i64 + 1),
                    2 => q(1),
                    _ => q(quarter as i64),
                };
                let free = 1 - fixed_axis;
                charge_visits(topology.uv.len(), dec)?;
                let mut points: Vec<_> = topology
                    .uv
                    .iter()
                    .enumerate()
                    .filter_map(|(id, uv)| {
                        let at = at_quarter(uv, quarter);
                        (in_quarter(&at, quarter) && at[fixed_axis] == fixed).then_some((
                            at[free].clone(),
                            id,
                            at,
                        ))
                    })
                    .collect();
                sort_by(&mut points, |a, b| a.0.cmp(&b.0), dec)?;
                if side >= 2 {
                    charge_visits(points.len(), dec)?;
                    points.reverse()
                }
                for pair in points.windows(2) {
                    if pair[0].0 == pair[1].0 {
                        return Err(bad(
                            "distinct owned events occupy the same source boundary parameter",
                        ));
                    }
                    let midpoint = std::array::from_fn(|i| (&pair[0].2[i] + &pair[1].2[i]) / q(2));
                    if !material.outside(source, &midpoint, dec)? {
                        continue;
                    }
                    uses.push(topology.boundary(
                        source,
                        &rings,
                        &intersection,
                        (fixed_axis, &fixed),
                        [
                            (pair[0].1, pair[0].2.clone()),
                            (pair[1].1, pair[1].2.clone()),
                        ],
                        dec,
                    )?);
                }
            }
            let loops = topology.walk(uses, dec)?;
            let mut outers = Vec::new();
            let mut holes = Vec::new();
            for boundary in loops {
                let polygon = topology.polygon(&boundary, dec)?;
                let area = area(&polygon, dec)?;
                match area.sign() {
                    Sign::Positive => outers.push((boundary, polygon)),
                    Sign::Negative => holes.push((boundary, polygon)),
                    Sign::Zero => {
                        return Err(bad(
                            "sampled finite loop has zero area; no loop was deleted",
                        ));
                    }
                }
            }
            let mut ownership = vec![Vec::new(); outers.len()];
            for (boundary, polygon) in holes {
                let mut owners = Vec::new();
                for (i, (_, outer)) in outers.iter().enumerate() {
                    dec.charge_work(1)?;
                    if inside(&polygon[0], outer, dec)? {
                        owners.push(i)
                    }
                }
                let [owner] = owners.as_slice() else {
                    return Err(bad("finite hole has no unique enclosing retained region"));
                };
                ownership[*owner].push(boundary);
            }
            for (region, ((outer, _), holes)) in outers.into_iter().zip(ownership).enumerate() {
                let loops = std::iter::once(outer)
                    .chain(holes)
                    .map(|row| {
                        if sense == Sense::Same {
                            row
                        } else {
                            reverse(row)
                        }
                    })
                    .collect();
                topology.definition.faces.push(NumericalFace {
                    key: format!("{key}/quarter-{quarter}/region-{region}"),
                    surface: surface.clone(),
                    sense,
                    loops,
                });
            }
        }
        let mut native_trims = Vec::new();
        for face in &topology.definition.faces {
            for (boundary, row) in face.loops.iter().enumerate() {
                for (use_index, trim) in row.iter().enumerate() {
                    if let Some(curve) = topology.native.get(&trim.edge) {
                        native_trims.push(NumericalNativeTrim {
                            tag: catallax_topology::numerical::CopyTag::Independent,
                            face: face.key.clone(),
                            boundary,
                            use_index,
                            edge: topology.definition.edges[trim.edge].key.clone(),
                            curve: curve.clone(),
                        })
                    }
                }
            }
        }
        compact(
            &mut topology.definition,
            &mut records,
            &mut topology.ring_links,
            &mut cylinder_trims,
            dec,
        )?;
        checks.extend(super::checks::open(
            &topology.definition,
            &native_trims,
            dec,
        )?);
        checks.push(NumericalCheck::Passed {
            name: "finite.owned-closed-native-loop-walks".into(),
            scope: key.into(),
        });
        dec.charge_work(checks.len() as u64 + repairs.len() as u64)?;
        let evidence = NumericalEvidence::new(options.tolerance_mm, checks, repairs)?;
        let work = NumericalFaceWork {
            boundaries: 0,
            intersection: intersection.work_stages(),
            supports: 0,
            topology: 0,
        };
        Ok(NumericalFacePartition {
            definition: topology.definition,
            native_trims,
            intersection,
            rings,
            ring_links: topology.ring_links,
            events: records,
            cylinder_source,
            cylinder_trims,
            evidence,
            work,
        })
    }
}

impl Topology<'_> {
    fn boundary(
        &mut self,
        source: &BSplineSurface,
        rings: &[NumericalConeRing],
        intersection: &NumericalNativeIntersection,
        coordinate: (usize, &Rational),
        ends: [(usize, [Rational; 2]); 2],
        dec: &mut Decider,
    ) -> Result<NumericalTrim, Refusal> {
        let (fixed_axis, fixed) = coordinate;
        let free = 1 - fixed_axis;
        let forward = ends[0].1[free] < ends[1].1[free];
        let sorted = if forward {
            [&ends[0], &ends[1]]
        } else {
            [&ends[1], &ends[0]]
        };
        let canonical = if fixed_axis == 0 && fixed == &q(4) {
            q(0)
        } else {
            fixed.clone()
        };
        let index = (fixed_axis, canonical.clone(), sorted[0].0, sorted[1].0);
        let edge = if let Some(index) = self.boundary_edges.get(&index) {
            *index
        } else {
            let count = source.control_count();
            let (degree, knots, controls, weights) = if fixed_axis == 1 {
                let row = usize::from(fixed == &q(1));
                (
                    source.degrees()[0],
                    source.knots()[0].to_vec(),
                    (0..count[0])
                        .map(|i| source.controls()[i * count[1] + row].clone())
                        .collect(),
                    (0..count[0])
                        .map(|i| source.weights()[i * count[1] + row].clone())
                        .collect(),
                )
            } else {
                let pole = (0..4)
                    .find(|i| canonical == q(*i as i64))
                    .ok_or_else(|| bad("source seam is not a native quarter boundary"))?
                    * 2;
                (
                    source.degrees()[1],
                    source.knots()[1].to_vec(),
                    (0..count[1])
                        .map(|i| source.controls()[pole * count[1] + i].clone())
                        .collect(),
                    (0..count[1])
                        .map(|i| source.weights()[pole * count[1] + i].clone())
                        .collect(),
                )
            };
            let curve = BSplineCurve::new(degree, knots, controls, weights, dec)?;
            let edge = self.definition.edges.len();
            self.definition.edges.push(NumericalEdge {
                key: format!(
                    "{}/source-{fixed_axis}-{canonical}-{}-{}",
                    self.key, sorted[0].0, sorted[1].0
                ),
                curve: Curve::BSpline(curve),
                vertices: [sorted[0].0, sorted[1].0],
                parameters: [sorted[0].1[free].clone(), sorted[1].1[free].clone()],
            });
            self.boundary_edges.insert(index, edge);
            if fixed_axis == 1 {
                let ring = usize::from(fixed == &q(1));
                let Surface::Cone(cone) = &intersection.windows().originals()[0] else {
                    unreachable!()
                };
                let native = [sorted[0].1[0].clone(), sorted[1].1[0].clone()];
                let angles = [
                    band::angle(&rings[ring], cone, &native[0], dec)?,
                    band::angle(&rings[ring], cone, &native[1], dec)?,
                ];
                self.ring_links.push(NumericalRingLink {
                    edge,
                    ring,
                    source_u: native,
                    original_parameters_rad: angles,
                });
            }
            edge
        };
        Ok(NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            edge,
            sense: if forward {
                Sense::Same
            } else {
                Sense::Reversed
            },
            uv: [ends[0].1.clone(), ends[1].1.clone()],
        })
    }
    fn vertex_ends(&self, trim: &NumericalTrim) -> [usize; 2] {
        let mut v = self.definition.edges[trim.edge].vertices;
        if trim.sense == Sense::Reversed {
            v.swap(0, 1)
        }
        v
    }
    fn walk(
        &self,
        uses: Vec<NumericalTrim>,
        dec: &mut Decider,
    ) -> Result<Vec<Vec<NumericalTrim>>, Refusal> {
        let mut outgoing = BTreeMap::<usize, Vec<usize>>::new();
        for (i, trim) in uses.iter().enumerate() {
            dec.charge_work(1)?;
            outgoing
                .entry(self.vertex_ends(trim)[0])
                .or_default()
                .push(i)
        }
        let mut passed = BTreeSet::new();
        let mut loops = Vec::new();
        for start in 0..uses.len() {
            dec.charge_work(1)?;
            if passed.contains(&start) {
                continue;
            }
            let mut row = Vec::new();
            let mut current = start;
            loop {
                dec.charge_work(1)?;
                if !passed.insert(current) {
                    return Err(bad("finite loop reused an already walked directed edge"));
                }
                row.push(uses[current].clone());
                let end = self.vertex_ends(&uses[current])[1];
                if end == self.vertex_ends(&uses[start])[0] {
                    break;
                }
                charge_visits(outgoing.get(&end).map_or(0, Vec::len), dec)?;
                let candidates: Vec<_> = outgoing
                    .get(&end)
                    .into_iter()
                    .flatten()
                    .filter(|i| !passed.contains(i))
                    .copied()
                    .collect();
                if candidates.is_empty() {
                    return Err(bad(format!(
                        "finite retained boundary is open at vertex {end}, uv={:?}",
                        self.uv[end].each_ref().map(Rational::to_f64_nearest),
                    )));
                }
                current = if candidates.len() == 1 {
                    candidates[0]
                } else {
                    self.next(&uses[current], &uses, &candidates, dec)?
                };
            }
            loops.push(row);
        }
        Ok(loops)
    }
    fn direction(
        &self,
        trim: &NumericalTrim,
        at_end: bool,
        divisor: i64,
        dec: &mut Decider,
    ) -> Result<[Rational; 2], Refusal> {
        if let Some(curve) = self.native.get(&trim.edge) {
            let edge = &self.definition.edges[trim.edge];
            let mut p = edge.parameters.clone();
            if trim.sense == Sense::Reversed {
                p.swap(0, 1)
            }
            if at_end {
                p.swap(0, 1)
            }
            let near = &p[0] + &((&p[1] - &p[0]) / q(divisor));
            let a = events::rational_point(curve.point_at(&Real::from(p[0].clone()), dec)?)?;
            let b = events::rational_point(curve.point_at(&Real::from(near), dec)?)?;
            Ok([&b[0] - &a[0], &b[1] - &a[1]])
        } else {
            let mut p = trim.uv.clone();
            if at_end {
                p.swap(0, 1)
            }
            Ok([&p[1][0] - &p[0][0], &p[1][1] - &p[0][1]])
        }
    }
    fn next(
        &self,
        incoming: &NumericalTrim,
        uses: &[NumericalTrim],
        candidates: &[usize],
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut previous = None;
        for divisor in [1024, 2048, 4096, 8192] {
            let a = self.direction(incoming, true, divisor, dec)?;
            let mut options = candidates
                .iter()
                .map(|i| {
                    dec.charge_work(1)?;
                    let b = self.direction(&uses[*i], false, divisor, dec)?;
                    Ok((
                        *i,
                        [&a[0] * &b[0] + &a[1] * &b[1], &a[1] * &b[0] - &a[0] * &b[1]],
                    ))
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            sort_by(
                &mut options,
                |a, b| angular(&a.1, &b.1).then_with(|| a.0.cmp(&b.0)),
                dec,
            )?;
            if options.len() > 1
                && angular(&options[0].1, &options[1].1) == std::cmp::Ordering::Equal
            {
                continue;
            }
            if previous == Some(options[0].0) {
                return Ok(options[0].0);
            }
            previous = Some(options[0].0);
        }
        Err(bad(
            "finite native inward directions need further tangent refinement",
        ))
    }
    fn polygon(
        &self,
        boundary: &[NumericalTrim],
        dec: &mut Decider,
    ) -> Result<Vec<[Rational; 2]>, Refusal> {
        let mut result = Vec::new();
        for trim in boundary {
            dec.charge_work(1)?;
            result.push(trim.uv[0].clone());
            if let Some(curve) = self.native.get(&trim.edge) {
                let edge = &self.definition.edges[trim.edge];
                let mut p = edge.parameters.clone();
                if trim.sense == Sense::Reversed {
                    p.swap(0, 1)
                }
                charge_visits(
                    curve
                        .knots()
                        .len()
                        .checked_mul(3)
                        .ok_or_else(|| bad("finite knot scan bound overflowed"))?,
                    dec,
                )?;
                let mut breaks: BTreeSet<_> = curve
                    .knots()
                    .iter()
                    .filter(|t| {
                        **t > p[0].clone().min(p[1].clone()) && **t < p[0].clone().max(p[1].clone())
                    })
                    .cloned()
                    .collect();
                breaks.extend(p.clone());
                let mut ordered: Vec<_> = breaks.into_iter().collect();
                if p[0] > p[1] {
                    ordered.reverse()
                }
                for span in ordered.windows(2) {
                    for i in 0..8 {
                        if span[0] == p[0] && i == 0 {
                            continue;
                        }
                        let t = &span[0]
                            + &((&span[1] - &span[0]) * Rational::from_ratio(i, 8).unwrap());
                        let at = events::rational_point(curve.point_at(&Real::from(t), dec)?)?;
                        result.push([at[0].clone(), at[1].clone()]);
                    }
                }
            }
        }
        Ok(result)
    }
}
fn angular(a: &[Rational; 2], b: &[Rational; 2]) -> std::cmp::Ordering {
    let half = |p: &[Rational; 2]| {
        usize::from(
            p[1].sign() == Sign::Negative
                || (p[1].sign() == Sign::Zero && p[0].sign() == Sign::Negative),
        )
    };
    half(a)
        .cmp(&half(b))
        .then_with(|| match (&a[0] * &b[1] - &a[1] * &b[0]).sign() {
            Sign::Positive => std::cmp::Ordering::Less,
            Sign::Negative => std::cmp::Ordering::Greater,
            Sign::Zero => std::cmp::Ordering::Equal,
        })
}
fn area(polygon: &[[Rational; 2]], dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut sum = q(0);
    for (a, b) in polygon
        .iter()
        .zip(polygon.iter().cycle().skip(1))
        .take(polygon.len())
    {
        dec.charge_work(1)?;
        sum = sum + (&a[0] * &b[1] - &a[1] * &b[0]);
    }
    Ok(sum)
}
fn inside(
    p: &[Rational; 2],
    polygon: &[[Rational; 2]],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut result = false;
    for (a, b) in polygon
        .iter()
        .zip(polygon.iter().cycle().skip(1))
        .take(polygon.len())
    {
        dec.charge_work(1)?;
        if (a[1] > p[1]) != (b[1] > p[1])
            && ((&b[0] - &a[0]) * (&p[1] - &a[1]) - (&b[1] - &a[1]) * (&p[0] - &a[0])).sign()
                == if b[1] > a[1] {
                    Sign::Positive
                } else {
                    Sign::Negative
                }
        {
            result = !result
        }
    }
    Ok(result)
}
fn compact(
    definition: &mut NumericalDefinition,
    events: &mut [NumericalFiniteEvent],
    rings: &mut Vec<NumericalRingLink>,
    cylinders: &mut Vec<NumericalCylinderTrim>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut visits = 0usize;
    for face in &definition.faces {
        dec.charge_work(1)?;
        visits = visits
            .checked_add(2)
            .ok_or_else(|| bad("finite face work bound overflowed"))?;
        for row in &face.loops {
            dec.charge_work(1)?;
            visits = visits
                .checked_add(2)
                .ok_or_else(|| bad("finite loop work bound overflowed"))?;
            visits = visits
                .checked_add(
                    row.len()
                        .checked_mul(2)
                        .ok_or_else(|| bad("finite trim work bound overflowed"))?,
                )
                .ok_or_else(|| bad("finite trim work bound overflowed"))?;
        }
    }
    for (count, repetitions) in [
        (definition.edges.len(), 8),
        (events.len(), 1),
        (rings.len(), 2),
        (cylinders.len(), 1),
    ] {
        visits = visits
            .checked_add(
                count
                    .checked_mul(repetitions)
                    .ok_or_else(|| bad("finite compaction work bound overflowed"))?,
            )
            .ok_or_else(|| bad("finite compaction work bound overflowed"))?;
    }
    charge_visits(visits, dec)?;
    let used_edges: BTreeSet<_> = definition
        .faces
        .iter()
        .flat_map(|f| f.loops.iter().flatten().map(|t| t.edge))
        .collect();
    let used_vertices: BTreeSet<_> = used_edges
        .iter()
        .flat_map(|i| definition.edges[*i].vertices)
        .collect();
    let vertices: BTreeMap<_, _> = used_vertices
        .iter()
        .enumerate()
        .map(|(new, old)| (*old, new))
        .collect();
    let edges: BTreeMap<_, _> = used_edges
        .iter()
        .enumerate()
        .map(|(new, old)| (*old, new))
        .collect();
    definition.vertices = used_vertices
        .into_iter()
        .map(|i| definition.vertices[i].clone())
        .collect();
    definition.edges = used_edges
        .into_iter()
        .map(|i| {
            let mut e = definition.edges[i].clone();
            e.vertices = e.vertices.map(|i| vertices[&i]);
            e
        })
        .collect();
    for trim in definition
        .faces
        .iter_mut()
        .flat_map(|f| f.loops.iter_mut().flatten())
    {
        trim.edge = edges[&trim.edge]
    }
    for event in events {
        event.vertex = event.vertex.and_then(|i| vertices.get(&i).copied())
    }
    rings.retain(|r| edges.contains_key(&r.edge));
    for r in rings {
        r.edge = edges[&r.edge]
    }
    for c in cylinders {
        c.edge = *edges
            .get(&c.edge)
            .ok_or_else(|| bad("a retained fitted arc was omitted from finite material loops"))?
    }
    dec.check_deadline()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hole_containment_charges_each_outer_polygon_visit() {
        let polygon: Vec<_> = (0..16)
            .flat_map(|_| [[q(0), q(0)], [q(1), q(0)], [q(1), q(1)], [q(0), q(1)]])
            .collect();
        let (answer, work) = catallax_numbers::with_operation_budget(8, || {
            inside(&[q(0), q(0)], &polygon, &mut Decider::default())
        });
        assert!(answer.is_err());
        assert!(
            work.exhausted,
            "outer polygon traversal escaped its allowance: {work:?}"
        );
        let (answer, work) = catallax_numbers::with_operation_budget(1024, || {
            inside(&[q(0), q(0)], &polygon, &mut Decider::default())
        });
        assert!(answer.is_ok());
        assert!(work.used >= 64);
        assert!(!work.exhausted);
    }

    #[test]
    fn polygon_visits_and_sort_honor_an_expired_deadline() {
        let polygon = [[q(0), q(0)], [q(1), q(0)], [q(1), q(1)], [q(0), q(1)]];
        let mut dec = Decider::default();
        dec.set_deadline(Some(
            std::time::Instant::now() - std::time::Duration::from_secs(1),
        ));
        assert_eq!(
            inside(&[q(0), q(0)], &polygon, &mut dec).unwrap_err().code,
            "numbers.deadline"
        );
        assert_eq!(
            area(&polygon, &mut dec).unwrap_err().code,
            "numbers.deadline"
        );
        let mut values = [3, 1, 2];
        assert_eq!(
            sort_by(&mut values, Ord::cmp, &mut dec).unwrap_err().code,
            "numbers.deadline"
        );
        assert_eq!(values, [3, 1, 2]);
    }

    #[test]
    fn sort_reserves_its_checked_bound_before_comparisons() {
        let mut values = [3, 1, 2];
        let (answer, work) = catallax_numbers::with_operation_budget(8, || {
            sort_by(&mut values, Ord::cmp, &mut Decider::default())
        });
        assert!(answer.is_err());
        assert!(work.exhausted);
        assert_eq!(values, [3, 1, 2]);
    }
}
