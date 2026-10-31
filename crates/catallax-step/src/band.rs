//! Export-only cuts for complete cylindrical, conical and toroidal bands.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, FactKind, Outcome, Real, Refusal, StrictOrder, StrictSign};
use catallax_shapes::{CurveId, GeomKey, Point3, Surface, SurfaceId, Torus, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, VertexId};

use crate::number;

pub(crate) enum Parameter {
    Exact(Real),
    /// A multiple of a quarter turn, on a torus's tube angle.
    QuarterTurns(i64),
}

pub(crate) struct Association {
    pub surface: SurfaceId,
    pub v: Parameter,
    pub v_scale: Real,
    pub increasing: bool,
}

pub(crate) struct Ring {
    pub origin: Point3,
    pub radius: Real,
    pub associations: Vec<Association>,
}

pub(crate) struct Band {
    pub surface: SurfaceId,
    /// Ring edges from low to high. A band ending at its cone's apex has
    /// only the high ring.
    pub edges: Vec<EdgeId>,
    pub forward: Vec<bool>,
    pub apex: Option<VertexId>,
    pub low_v: Real,
    pub v_scale: Real,
    pub radius_bound: Real,
    pub seam_direction: Vector3,
    /// A closed chain of trimmed edges in place of one complete ring.
    pub chain: Option<Chain>,
}

/// Coedges of one closed boundary, walked from the seam vertex in the chart's
/// increasing direction at the low end and its decreasing one at the high end.
pub(crate) struct Chain {
    pub uses: Vec<(EdgeId, bool)>,
    pub vertex: VertexId,
    pub at_low: bool,
}

#[derive(Default)]
pub(crate) struct Plan {
    pub faces: BTreeMap<FaceId, Rc<Band>>,
    pub rings: BTreeMap<EdgeId, Rc<Ring>>,
    pub curves: BTreeMap<CurveId, CurveForm>,
    pub references: BTreeMap<SurfaceId, Vector3>,
    pub envelope_edges: BTreeMap<EdgeId, Rc<crate::junction::EnvelopeEdge>>,
    /// Circular arcs of a junction, with one exact line per curved face.
    pub arcs: BTreeMap<EdgeId, Rc<Ring>>,
}

pub(crate) struct SourceRing {
    pub(crate) edge: EdgeId,
    pub(crate) curve: CurveId,
    pub(crate) center: Point3,
    pub(crate) normal: Vector3,
    pub(crate) radius_sq: Real,
    pub(crate) radius: Real,
    pub(crate) v: Real,
    /// The tube angle on a torus, which replaces `v`.
    pub(crate) turns: Option<i64>,
    pub(crate) increasing: bool,
}

pub(crate) struct SourceBand {
    pub(crate) face: FaceId,
    pub(crate) surface: SurfaceId,
    pub(crate) axis: Vector3,
    pub(crate) reference: Vector3,
    pub(crate) v_scale: Real,
    pub(crate) rings: Vec<SourceRing>,
    pub(crate) apex: Option<VertexId>,
    pub(crate) chain: Option<(Chain, Real, Real)>,
}

impl Plan {
    pub(crate) fn prepare(
        model: &Model,
        bodies: impl Iterator<Item = BodyId> + Clone,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if crate::junction::present(model, bodies.clone()) {
            return crate::junction::plan(model, bodies, dec);
        }
        if crate::junction::tangent_present(model, bodies.clone()) {
            return crate::junction::tangent_plan(model, bodies, dec);
        }
        let mut has_band = false;
        let mut has_cone = false;
        let mut has_torus = false;
        for body in bodies.clone() {
            for lump in &model.body(body).lumps {
                for shell in &model.lump(*lump).shells {
                    for &face in &model.shell(*shell).faces {
                        let definition = model.face(face);
                        match model.geom().surface(definition.surface) {
                            Surface::Plane(_) => {}
                            Surface::Sphere(_) | Surface::Torus(_)
                                if definition.loops.is_empty() => {}
                            Surface::Cylinder(_) | Surface::Cone(_)
                                if band_loops(model, face).is_some() =>
                            {
                                has_band = true;
                                has_cone |= matches!(
                                    model.geom().surface(definition.surface),
                                    Surface::Cone(_)
                                );
                            }
                            Surface::Torus(_) if band_loops(model, face).is_some() => {
                                has_torus = true
                            }
                            _ => return Ok(Self::default()),
                        }
                    }
                }
            }
        }
        // A torus band beside a cone has no chart route, so only that scene
        // comes here. Every other torus keeps its charted bytes.
        if !has_band || (has_torus && !has_cone) {
            return Ok(Self::default());
        }
        let mut sources = Vec::new();
        let mut seen = BTreeSet::new();
        for body in bodies {
            for lump in &model.body(body).lumps {
                for shell in &model.lump(*lump).shells {
                    for &face in &model.shell(*shell).faces {
                        if !seen.insert(face) {
                            continue;
                        }
                        let definition = model.face(face);
                        match model.geom().surface(definition.surface) {
                            Surface::Plane(_) => {}
                            Surface::Sphere(_) | Surface::Torus(_)
                                if definition.loops.is_empty() => {}
                            Surface::Cylinder(_) | Surface::Cone(_) => {
                                let Some(band) = source_band(model, face, dec)? else {
                                    // Do not partly change a file that still needs the
                                    // arbitrary-trim or singular-apex export contract.
                                    return Ok(Self::default());
                                };
                                sources.push(band);
                            }
                            Surface::Torus(_) => {
                                let Some(band) = torus_band(model, face, dec)? else {
                                    return Ok(Self::default());
                                };
                                sources.push(band);
                            }
                            _ => return Ok(Self::default()),
                        }
                    }
                }
            }
        }
        Ok(assemble(model, sources))
    }
}

// Connected bands must cut shared rings at the same physical meridian. A
// reference rotation changes the parametrisation, not the carrier. A chain
// fixes its component's reference, because its seam starts at a vertex.
pub(crate) fn assemble(model: &Model, sources: Vec<SourceBand>) -> Plan {
    let mut parents: Vec<_> = (0..sources.len()).collect();
    let mut edge_owners = BTreeMap::new();
    let mut curve_owners = BTreeMap::new();
    let mut surface_owners = BTreeMap::new();
    for (index, band) in sources.iter().enumerate() {
        if let Some(previous) = surface_owners.insert(band.surface, index) {
            join(&mut parents, index, previous);
        }
        for ring in &band.rings {
            if let Some(previous) = edge_owners.insert(ring.edge, index) {
                join(&mut parents, index, previous);
            }
            if let Some(previous) = curve_owners.insert(ring.curve, index) {
                join(&mut parents, index, previous);
            }
        }
    }
    let mut plan = Plan::default();
    let mut rings: BTreeMap<EdgeId, Ring> = BTreeMap::new();
    let mut chains = BTreeMap::new();
    for (index, source) in sources.iter().enumerate() {
        if source.chain.is_some() {
            chains.insert(root(&mut parents, index), &source.reference);
        }
    }
    for (index, source) in sources.iter().enumerate() {
        let parent = root(&mut parents, index);
        let reference = chains
            .get(&parent)
            .copied()
            .unwrap_or(&sources[parent].reference);
        plan.references.insert(source.surface, reference.clone());
        // Use the certified carrier frame, not a normalized difference
        // of algebraic endpoint heights that must be re-proved unit.
        let seam_direction = match model.geom().surface(source.surface) {
            Surface::Cylinder(_) => source.axis.clone(),
            Surface::Cone(cone) => source.axis.plus(&reference.scaled(&cone.slope)),
            // The seam is a tube circle, which turns about this axis's normal.
            Surface::Torus(_) => source.axis.clone(),
            _ => unreachable!(),
        };
        plan.faces.insert(
            source.face,
            Rc::new(Band {
                surface: source.surface,
                edges: source.rings.iter().map(|r| r.edge).collect(),
                forward: source
                    .rings
                    .iter()
                    .map(|r| r.increasing == (model.edge(r.edge).sense == Sense::Same))
                    .collect(),
                apex: source.apex,
                // The written cone's apex is its v origin.
                low_v: match &source.chain {
                    Some((chain, v, _)) if chain.at_low => v.clone(),
                    _ if source.apex.is_some() => Real::zero(),
                    // The seam circle's own parameter is the tube angle.
                    _ if source.rings[0].turns.is_some() => Real::zero(),
                    _ => source.rings[0].v.clone(),
                },
                v_scale: source.v_scale.clone(),
                radius_bound: source.rings.iter().fold(
                    source
                        .chain
                        .as_ref()
                        .map_or_else(Real::zero, |(_, _, radius)| radius.clone()),
                    |sum, r| &sum + &r.radius,
                ),
                seam_direction,
                chain: source.chain.as_ref().map(|(chain, _, _)| Chain {
                    uses: chain.uses.clone(),
                    vertex: chain.vertex,
                    at_low: chain.at_low,
                }),
            }),
        );
        for ring in &source.rings {
            let origin = number::point_plus(&ring.center, reference, &ring.radius);
            rings
                .entry(ring.edge)
                .or_insert_with(|| Ring {
                    origin,
                    radius: ring.radius.clone(),
                    associations: Vec::with_capacity(2),
                })
                .associations
                .push(Association {
                    surface: source.surface,
                    v: ring
                        .turns
                        .map_or_else(|| Parameter::Exact(ring.v.clone()), Parameter::QuarterTurns),
                    v_scale: source.v_scale.clone(),
                    increasing: ring.increasing,
                });
            plan.curves
                .entry(ring.curve)
                .or_insert_with(|| CurveForm::Circle {
                    center: ring.center.clone(),
                    normal: ring.normal.clone(),
                    radius_sq: ring.radius_sq.clone(),
                    ref_dir: reference.clone(),
                });
        }
    }
    plan.rings = rings
        .into_iter()
        .map(|(edge, ring)| (edge, Rc::new(ring)))
        .collect();
    plan
}

// Two complete rings, or for a cone one complete ring and a vertex loop.
fn band_loops(model: &Model, face: FaceId) -> Option<(Vec<EdgeId>, Option<VertexId>)> {
    let definition = model.face(face);
    let [a, b] = definition.loops.as_slice() else {
        return None;
    };
    let edge = |loop_id| {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            return None;
        };
        let [coedge] = coedges.as_slice() else {
            return None;
        };
        let edge = model.coedge(*coedge).edge;
        let definition = model.edge(edge);
        (definition.start.is_none() && definition.end.is_none()).then_some(edge)
    };
    let cone = matches!(model.geom().surface(definition.surface), Surface::Cone(_));
    match (model.loop_(*a), model.loop_(*b)) {
        (Loop::Vertex(v), _) if cone => Some((vec![edge(*b)?], Some(*v))),
        (_, Loop::Vertex(v)) if cone => Some((vec![edge(*a)?], Some(*v))),
        _ => Some((vec![edge(*a)?, edge(*b)?], None)),
    }
}

pub(crate) fn source_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let face_definition = model.face(face);
    let Some((edges, apex)) = band_loops(model, face) else {
        return Ok(None);
    };
    let surface = model.geom().surface(face_definition.surface);
    let Some((origin, axis, reference, slant)) = carrier_frame(surface) else {
        return Ok(None);
    };
    if let Some(vertex) = apex {
        let point = model.geom().position(model.vertex(vertex).point);
        let constructed = [
            (&point.x, &origin.x),
            (&point.y, &origin.y),
            (&point.z, &origin.z),
        ]
        .into_iter()
        .all(|(a, b)| a.same_definition(b));
        if !constructed && !zero(&point.minus(origin).norm_squared(), dec)? {
            return Ok(None);
        }
    }
    let mut rings = Vec::with_capacity(edges.len());
    for edge in edges {
        let Some(ring) = source_ring(model, face, edge, dec)? else {
            return Ok(None);
        };
        rings.push(ring);
    }
    if rings.len() == 2 {
        match dec.compare(&rings[0].v, &rings[1].v) {
            Outcome::Proven(StrictOrder::Greater) => rings.swap(0, 1),
            Outcome::Proven(StrictOrder::Less) => {}
            Outcome::Refused(refusal) => return Err(refusal),
            _ => return Ok(None),
        }
    }
    Ok(Some(SourceBand {
        face,
        surface: face_definition.surface,
        axis,
        reference,
        v_scale: slant,
        rings,
        apex,
        chain: None,
    }))
}

pub(crate) fn carrier_frame(surface: &Surface) -> Option<(&Point3, Vector3, Vector3, Real)> {
    match surface {
        Surface::Cylinder(s) => Some((
            &s.origin,
            number::unit(&s.axis),
            number::unit(&s.ref_dir),
            Real::one(),
        )),
        Surface::Cone(s) => Some((
            &s.apex,
            number::unit(&s.axis),
            number::unit(&s.ref_dir),
            number::sqrt_real(&(Real::one() + &s.slope * &s.slope)),
        )),
        _ => None,
    }
}

pub(crate) fn source_ring(
    model: &Model,
    face: FaceId,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<Option<SourceRing>, Refusal> {
    let face_definition = model.face(face);
    let surface = model.geom().surface(face_definition.surface);
    let Some((origin, axis, ..)) = carrier_frame(surface) else {
        return Ok(None);
    };
    {
        let definition = model.edge(edge);
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), definition.curve, dec)?
        else {
            return Ok(None);
        };
        let offset = center.minus(origin);
        let height = offset.dot(&axis);
        let curve_key = GeomKey::Curve(definition.curve);
        let surface_key = GeomKey::Surface(face_definition.surface);
        if model
            .geom()
            .facts()
            .holds(FactKind::Coaxial, curve_key, surface_key)
        {
            dec.note_fact_decision();
        } else if !zero(&normal.cross(&axis).norm_squared(), dec)?
            || !zero(&offset.cross(&axis).norm_squared(), dec)?
        {
            return Ok(None);
        }
        if matches!(surface, Surface::Cone(_))
            && !matches!(dec.sign(&height), Outcome::Proven(StrictSign::Positive))
        {
            return Ok(None);
        }
        // Whole-carrier incidence already proves the circular section's
        // radius equation. Re-eliminating an algebraic drill angle here can
        // exhaust the request after construction has certified that identity.
        if model
            .geom()
            .facts()
            .holds(FactKind::Incident, curve_key, surface_key)
        {
            dec.note_fact_decision();
        } else {
            let expected_radius_sq = match surface {
                Surface::Cylinder(s) => s.radius_sq.clone(),
                Surface::Cone(s) => (&s.slope * &height).square(),
                _ => unreachable!(),
            };
            if !zero(&(&radius_sq - &expected_radius_sq), dec)? {
                return Ok(None);
            }
        }
        let increasing = match dec.sign(&normal.dot(&axis)) {
            Outcome::Proven(sign) => sign == StrictSign::Positive,
            Outcome::Refused(refusal) => return Err(refusal),
            _ => return Ok(None),
        };
        Ok(Some(SourceRing {
            edge,
            curve: definition.curve,
            center,
            normal,
            radius: number::sqrt_real(&radius_sq),
            radius_sq,
            v: height,
            turns: None,
            increasing,
        }))
    }
}

// The band of a torus between a tube-side parallel at a quarter turn and its
// innermost parallel, as a concave ring fillet leaves it. The tube angle of a
// ring is an exact quarter turn; any other parallel keeps the unseamed encoder.
pub(crate) fn torus_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let definition = model.face(face);
    let Surface::Torus(torus) = model.geom().surface(definition.surface) else {
        return Ok(None);
    };
    let Some((edges, None)) = band_loops(model, face) else {
        return Ok(None);
    };
    let axis = number::unit(&torus.axis);
    let mut rings = Vec::with_capacity(2);
    for edge in edges {
        let Some(ring) = torus_ring(model, torus, &axis, edge, dec)? else {
            return Ok(None);
        };
        rings.push(ring);
    }
    rings.sort_by_key(|ring| ring.turns);
    if rings[0].turns == rings[1].turns {
        return Ok(None);
    }
    // Two parallels bound two bands. The low ring runs with the chart's u, read
    // on the face's side, only on the quarter between them.
    let low = &rings[0];
    let along = model.edge(low.edge).sense == Sense::Same;
    let uses = ring_use(model, face, low.edge)?;
    if (low.increasing == (along == uses)) != (definition.sense == Sense::Same) {
        return Ok(None);
    }
    Ok(Some(SourceBand {
        face,
        surface: definition.surface,
        axis,
        reference: number::unit(&torus.ref_dir),
        v_scale: torus.minor.clone(),
        rings,
        apex: None,
        chain: None,
    }))
}

// A torus parallel at an exact quarter turn of the tube, or `None`.
pub(crate) fn torus_ring(
    model: &Model,
    torus: &Torus,
    axis: &Vector3,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<Option<SourceRing>, Refusal> {
    let source = model.edge(edge);
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = curve_form(model.geom(), source.curve, dec)?
    else {
        return Ok(None);
    };
    let offset = center.minus(&torus.center);
    if !zero(&normal.cross(axis).norm_squared(), dec)?
        || !zero(&offset.cross(axis).norm_squared(), dec)?
    {
        return Ok(None);
    }
    let height = offset.dot(axis);
    let radius = number::sqrt_real(&radius_sq);
    let turns = if zero(&height, dec)? {
        if !zero(&(&radius - &(&torus.major - &torus.minor)), dec)? {
            return Ok(None);
        }
        2
    } else {
        if !zero(&(&radius - &torus.major), dec)? {
            return Ok(None);
        }
        if zero(&(&height - &torus.minor), dec)? {
            1
        } else if zero(&(&height + &torus.minor), dec)? {
            3
        } else {
            return Ok(None);
        }
    };
    let increasing = match dec.sign(&normal.dot(axis)) {
        Outcome::Proven(sign) => sign == StrictSign::Positive,
        Outcome::Refused(refusal) => return Err(refusal),
        _ => return Ok(None),
    };
    Ok(Some(SourceRing {
        edge,
        curve: source.curve,
        center,
        normal,
        radius,
        radius_sq,
        v: Real::zero(),
        turns: Some(turns),
        increasing,
    }))
}

// Whether a ring's coedge in this face follows its edge.
pub(crate) fn ring_use(model: &Model, face: FaceId, edge: EdgeId) -> Result<bool, Refusal> {
    for lp in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(*lp)
            && let [coedge] = coedges.as_slice()
            && model.coedge(*coedge).edge == edge
        {
            return Ok(model.coedge(*coedge).sense == Sense::Same);
        }
    }
    Err(Refusal::new(
        "step.not-solid",
        "a band ring is not a loop of its face",
    ))
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Refused(refusal) => Err(refusal),
        _ => Ok(false),
    }
}

fn root(parents: &mut [usize], mut index: usize) -> usize {
    while parents[index] != index {
        parents[index] = parents[parents[index]];
        index = parents[index];
    }
    index
}

fn join(parents: &mut [usize], a: usize, b: usize) {
    let a = root(parents, a);
    let b = root(parents, b);
    parents[b] = a;
}
