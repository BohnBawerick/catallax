//! Exact removal of numerical chart seams from the published topology.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{BSplineCurve, Curve, RationalAtlas, SplineSupport, Surface};

use crate::*;

/// A complete logical atlas, with one input face per chart in chart order.
/// All input faces must use those exact checked chart definition nodes.
pub struct AtlasFaceGroup {
    /// Connected regular exact atlas, including its G1 seam identifications.
    pub atlas: RationalAtlas,
    /// Input chart faces, in `atlas.charts()` order.
    pub faces: Vec<FaceId>,
}

/// Exact cell topology retained for numerical integration and exchange. It is
/// not selectable model topology. The constructor proves that coalescing only
/// removes chart seams and joins exact boundary curves; it proves no new global
/// nonintersection claim about an arbitrary caller-assembled input body.
#[derive(Clone, Debug)]
pub struct ChartExpansion {
    pub(crate) body: BodyId,
    pub(crate) faces: BTreeMap<FaceId, FaceId>,
    pub(crate) edges: BTreeMap<EdgeId, Vec<(EdgeId, bool)>>,
    pub(crate) needs_cells: bool,
}

impl ChartExpansion {
    /// Exact sewn chart body, with no approximation or fitted replacement.
    pub fn body(&self) -> BodyId {
        self.body
    }
    /// Whether a nonanalytic logical spline face still requires its cell cover.
    /// Exact plane/line reductions use the ordinary analytic consumers.
    pub fn needs_cells(&self) -> bool {
        self.needs_cells
    }
    /// Logical picking face corresponding to one expanded face.
    pub fn face(&self, expanded: FaceId) -> FaceId {
        self.faces[&expanded]
    }
    /// Logical edges and their expanded pieces in logical start-to-end order.
    /// A Boolean says whether that piece's start-to-end order is retained.
    pub fn edges(&self) -> &BTreeMap<EdgeId, Vec<(EdgeId, bool)>> {
        &self.edges
    }
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("topology.rational-atlas", detail)
}
fn reverse(sense: Sense) -> Sense {
    if sense == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}

struct Group {
    surface: catallax_shapes::SurfaceId,
    sense: Sense,
    origins: Vec<catallax_recipe::Origin>,
}

type Use = (usize, Sense);
type Chain = Vec<(EdgeId, bool)>;

fn chains(
    model: &Model,
    edges: &[EdgeId],
    degree: &BTreeMap<VertexId, usize>,
) -> Result<Vec<Chain>, Refusal> {
    let mut adjacent: BTreeMap<VertexId, Vec<EdgeId>> = BTreeMap::new();
    let mut unused: BTreeSet<_> = edges.iter().copied().collect();
    let mut result = Vec::new();
    for edge in edges {
        let value = model.edge(*edge);
        match (value.start, value.end) {
            (Some(a), Some(b)) => {
                adjacent.entry(a).or_default().push(*edge);
                adjacent.entry(b).or_default().push(*edge);
            }
            (None, None) => {
                result.push(vec![(*edge, true)]);
                unused.remove(edge);
            }
            _ => return Err(bad("a physical boundary edge has only one endpoint")),
        }
    }
    let kept = |vertex: &VertexId| degree[vertex] != 2 || adjacent[vertex].len() != 2;
    while let Some(first) = unused.first().copied() {
        let start = unused
            .iter()
            .find_map(|edge| {
                let value = model.edge(*edge);
                [value.start, value.end]
                    .into_iter()
                    .flatten()
                    .find(|v| kept(v))
            })
            .unwrap_or_else(|| model.edge(first).start.expect("non-ring boundary"));
        let mut at = start;
        let mut chain = Vec::new();
        loop {
            let edge = adjacent[&at]
                .iter()
                .find(|edge| unused.contains(edge))
                .copied()
                .ok_or_else(|| bad("a boundary chain is not sewn"))?;
            unused.remove(&edge);
            let value = model.edge(edge);
            let forward = value.start == Some(at);
            at = if forward { value.end } else { value.start }.expect("non-ring boundary");
            chain.push((edge, forward));
            if at == start || kept(&at) {
                break;
            }
        }
        result.push(chain);
    }
    Ok(result)
}

fn loops(model: &mut Model, uses: BTreeMap<EdgeId, Sense>) -> Result<Vec<LoopId>, Refusal> {
    let mut outgoing: BTreeMap<VertexId, (EdgeId, Sense)> = BTreeMap::new();
    let mut rings = Vec::new();
    for (edge, sense) in uses {
        let value = model.edge(edge);
        let start = if sense == Sense::Same {
            value.start
        } else {
            value.end
        };
        if let Some(start) = start {
            if outgoing.insert(start, (edge, sense)).is_some() {
                return Err(bad("a logical face boundary branches"));
            }
        } else {
            rings.push((edge, sense));
        }
    }
    let mut result = Vec::new();
    for (edge, sense) in rings {
        let coedge = model.add_coedge(Coedge { edge, sense });
        result.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    while let Some(start) = outgoing.first_key_value().map(|(v, _)| *v) {
        let mut at = start;
        let mut coedges = Vec::new();
        loop {
            let (edge, sense) = outgoing
                .remove(&at)
                .ok_or_else(|| bad("a logical face boundary is open"))?;
            let value = model.edge(edge);
            at = if sense == Sense::Same {
                value.end
            } else {
                value.start
            }
            .ok_or_else(|| bad("a finite logical edge has no endpoint"))?;
            coedges.push(model.add_coedge(Coedge { edge, sense }));
            if at == start {
                break;
            }
        }
        result.push(model.add_loop(Loop::Coedges(coedges)));
    }
    Ok(result)
}

/// Publishes maximal logical atlas faces and only physical boundary edges.
/// The exact input cell body remains an immutable consumer expansion. Every
/// chart face and incidence is checked before any seam can disappear.
pub fn coalesce_rational_charts(
    model: &mut Model,
    body: BodyId,
    atlases: Vec<AtlasFaceGroup>,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    if !model.check(body).is_empty() {
        return Err(bad("the input cell body is not a sewn solid"));
    }
    let all: BTreeSet<_> = model.faces_of_body(body).into_iter().collect();
    let mut owner: BTreeMap<FaceId, usize> = BTreeMap::new();
    let mut groups = Vec::new();
    for atlas in atlases {
        if atlas.faces.len() != atlas.atlas.charts().len() || atlas.faces.is_empty() {
            return Err(bad("the atlas needs exactly one face per chart"));
        }
        let sense = model.face(atlas.faces[0]).sense;
        let mut origins = Vec::new();
        let mut chart_rows = Vec::with_capacity(atlas.faces.len());
        for (face, chart) in atlas.faces.iter().zip(atlas.atlas.charts()) {
            if !all.contains(face) || owner.insert(*face, groups.len()).is_some() {
                return Err(bad("a chart face is absent or claimed twice"));
            }
            let Some(checked) = RationalPatchFace::of_face(model, *face, dec)? else {
                return Err(bad("an atlas cell is not a complete rational patch face"));
            };
            if model.face(*face).sense != sense
                || !checked.support().shares_definition_nodes(chart.patch())
            {
                return Err(bad(
                    "an atlas cell changes its checked definition or orientation",
                ));
            }
            chart_rows.push(checked.rows().map(|row| row.0));
            origins.extend_from_slice(model.face_origins(*face));
        }
        let mut row_counts = BTreeMap::new();
        for edge in chart_rows.iter().flatten() {
            *row_counts.entry(*edge).or_insert(0usize) += 1;
        }
        let sides = [
            catallax_shapes::SplineBoundary::VMin,
            catallax_shapes::SplineBoundary::VMax,
            catallax_shapes::SplineBoundary::UMin,
            catallax_shapes::SplineBoundary::UMax,
        ];
        for (index, chart) in atlas.atlas.charts().iter().enumerate() {
            for (row, side) in sides.iter().enumerate() {
                if let Some(neighbor) = chart.neighbor(*side) {
                    if chart_rows[index][row] != chart_rows[neighbor][row ^ 1] {
                        return Err(bad(
                            "the sewn cell edges do not follow the proved atlas transitions",
                        ));
                    }
                } else if row_counts[&chart_rows[index][row]] > 1 {
                    return Err(bad("a logical component meets itself at an unjoined row"));
                }
            }
        }
        let surface = match atlas.atlas.rational_plane() {
            Some(plane) => Surface::Plane(plane),
            None => Surface::BSpline(atlas.atlas.into()),
        };
        let surface = model.geom_mut().add_surface(surface);
        groups.push(Group {
            surface,
            sense,
            origins,
        });
    }
    for face in &all {
        if owner.contains_key(face) {
            continue;
        }
        let value = model.face(*face);
        if matches!(
            model.geom().surface(value.surface),
            Surface::BSpline(SplineSupport::RationalAtlas(_))
        ) {
            return Err(bad("nested chart expansions must be flattened first"));
        }
        owner.insert(*face, groups.len());
        groups.push(Group {
            surface: value.surface,
            sense: value.sense,
            origins: model.face_origins(*face).to_vec(),
        });
    }
    let mut incidences: BTreeMap<EdgeId, Vec<Use>> = BTreeMap::new();
    for face in &all {
        for lp in &model.face(*face).loops {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                return Err(bad("an atlas complex cannot contain apex-only loops"));
            };
            for coedge in coedges {
                let c = model.coedge(*coedge);
                incidences
                    .entry(c.edge)
                    .or_default()
                    .push((owner[face], c.sense));
            }
        }
    }
    let mut pairs: BTreeMap<[usize; 2], Vec<EdgeId>> = BTreeMap::new();
    let mut degree: BTreeMap<VertexId, usize> = BTreeMap::new();
    for (edge, uses) in &incidences {
        let [(a, sa), (b, sb)] = uses.as_slice() else {
            return Err(bad("an input edge lacks two face uses"));
        };
        if sa == sb {
            return Err(bad("an input edge has inconsistent material traversal"));
        }
        if a == b {
            continue;
        }
        let mut pair = [*a, *b];
        pair.sort_unstable();
        pairs.entry(pair).or_default().push(*edge);
        let value = model.edge(*edge);
        for v in [value.start, value.end].into_iter().flatten() {
            *degree.entry(v).or_default() += 1;
        }
    }
    let mut logical_uses = vec![BTreeMap::new(); groups.len()];
    let mut edge_map = BTreeMap::new();
    for edges in pairs.values() {
        for mut chain in chains(model, edges, &degree)? {
            let logical = if chain.len() == 1 {
                chain[0].1 = true;
                chain[0].0
            } else {
                let first = model.edge(chain[0].0);
                let last = model.edge(chain.last().expect("nonempty chain").0);
                let start = if chain[0].1 { first.start } else { first.end };
                let end = if chain.last().expect("nonempty chain").1 {
                    last.end
                } else {
                    last.start
                };
                let ring = start == end;
                let pieces = chain
                    .iter()
                    .map(|(edge, forward)| {
                        let value = model.edge(*edge);
                        let Curve::BSpline(curve) = model.geom().curve(value.curve) else {
                            return Err(bad(
                                "joining analytic boundary pieces needs an exact curve reduction",
                            ));
                        };
                        Ok((curve.clone(), *forward == (value.sense == Sense::Same)))
                    })
                    .collect::<Result<Vec<_>, Refusal>>()?;
                let curve = BSplineCurve::joined(&pieces, dec)?;
                let curve = model.geom_mut().add_curve(Curve::BSpline(curve));
                model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: (!ring).then_some(start).flatten(),
                    end: (!ring).then_some(end).flatten(),
                })
            };
            let logical = if let Some(line) = rational_line(model, logical) {
                let edge = model.edge(logical).clone();
                let curve = model.geom_mut().add_curve(Curve::Line(line));
                model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: edge.start,
                    end: edge.end,
                })
            } else {
                logical
            };
            let constructions = chain
                .iter()
                .flat_map(|(edge, _)| model.edge_constructions(*edge).iter().cloned())
                .collect();
            model.record_edge_constructions(logical, constructions);
            for (edge, forward) in &chain {
                for (group, sense) in &incidences[edge] {
                    let sense = if *forward { *sense } else { reverse(*sense) };
                    if logical_uses[*group]
                        .insert(logical, sense)
                        .is_some_and(|old| old != sense)
                    {
                        return Err(bad("a joined boundary reverses traversal within one face"));
                    }
                }
            }
            edge_map.insert(logical, chain);
        }
    }
    let mut logical_faces = Vec::with_capacity(groups.len());
    for (group, uses) in groups.into_iter().zip(logical_uses) {
        let loops = loops(model, uses)?;
        record_incidence(model, &loops, group.surface);
        logical_faces.push(model.add_named_face(
            Face {
                surface: group.surface,
                sense: group.sense,
                loops,
            },
            group.origins,
        ));
    }
    let face_map: BTreeMap<_, _> = owner
        .iter()
        .map(|(old, group)| (*old, logical_faces[*group]))
        .collect();
    let mut covers = vec![Vec::new(); logical_faces.len()];
    for (old, group) in &owner {
        covers[*group].push(*old);
    }
    for (face, cover) in logical_faces.iter().zip(covers) {
        model.face_expansions.insert(*face, Arc::from(cover));
    }
    let source = model.body(body).clone();
    let mut lumps = Vec::new();
    let mut seen_groups = BTreeSet::new();
    for lump in source.lumps {
        let mut shells = Vec::new();
        for shell in model.lump(lump).shells.clone() {
            let ids: BTreeSet<_> = model
                .shell(shell)
                .faces
                .iter()
                .map(|f| face_map[f])
                .collect();
            if ids.iter().any(|f| !seen_groups.insert(*f)) {
                return Err(bad("a logical face spans different input shells"));
            }
            shells.push(model.add_shell(Shell {
                faces: ids.into_iter().collect(),
            }));
        }
        lumps.push(model.add_lump(Lump { shells }));
    }
    let logical = model.add_body(Body {
        kind: source.kind,
        lumps,
        marks: source.marks,
    });
    if !model.check(logical).is_empty() {
        return Err(bad("the coalesced logical body is not sewn"));
    }
    let needs_cells = logical_faces.iter().any(|face| {
        matches!(
            model.geom().surface(model.face(*face).surface),
            Surface::BSpline(SplineSupport::RationalAtlas(_))
        )
    });
    model.chart_expansions.insert(
        logical,
        Arc::new(ChartExpansion {
            body,
            faces: face_map,
            edges: edge_map,
            needs_cells,
        }),
    );
    Ok(logical)
}

fn record_incidence(model: &mut Model, loops: &[LoopId], surface: catallax_shapes::SurfaceId) {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::GeomKey;
    for lp in loops {
        let count = match model.loop_(*lp) {
            Loop::Coedges(cs) => cs.len(),
            Loop::Vertex(_) => 0,
        };
        for index in 0..count {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                unreachable!()
            };
            let edge = model.edge(model.coedge(cs[index]).edge).clone();
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(edge.curve),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("topology.rational-atlas"),
            });
            for vertex in [edge.start, edge.end].into_iter().flatten() {
                let point = model.vertex(vertex).point;
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Point(point),
                    b: GeomKey::Curve(edge.curve),
                    at: None,
                    origin: FactOrigin::Construction("topology.rational-atlas"),
                });
            }
        }
    }
}

fn rational_line(model: &Model, edge: EdgeId) -> Option<catallax_shapes::Line> {
    let edge = model.edge(edge);
    let Curve::BSpline(curve) = model.geom().curve(edge.curve) else {
        return None;
    };
    let origin = model
        .geom()
        .position(model.vertex(edge.start?).point)
        .as_rational()?;
    let end = model
        .geom()
        .position(model.vertex(edge.end?).point)
        .as_rational()?;
    let direction = end.minus(&origin);
    let length = direction.norm_squared();
    if length.is_zero() {
        return None;
    }
    for point in curve.controls() {
        let relative = point.as_rational()?.minus(&origin);
        let along = relative.dot(&direction);
        if !relative.cross(&direction).is_zero()
            || along.sign() == catallax_numbers::Sign::Negative
            || along > length
        {
            return None;
        }
    }
    // Positive weights keep the complete curve in this segment. Continuity
    // and the two proved endpoints give complete coverage, not a sampled fit.
    Some(catallax_shapes::Line {
        origin: origin.to_real(),
        dir: direction.to_real(),
    })
}
