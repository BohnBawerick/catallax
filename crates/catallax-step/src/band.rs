//! Export-only cuts for complete cylindrical, conical and toroidal bands.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, FactKind, Outcome, Real, Refusal, StrictOrder, StrictSign};
use catallax_shapes::{CurveId, GeomKey, Point3, Surface, SurfaceId, Torus, TorusSheet, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, VertexId};

use crate::number;

pub(crate) enum Parameter {
    Exact(Real),
    /// A multiple of a quarter turn, on a torus's tube angle.
    QuarterTurns(i64),
    /// Any other tube angle on a torus.
    Tube(Tube),
}

/// The tube angle `atan2(height, across)` in `[0, 2 pi)`, one turn more when
/// `wraps`: `height` runs along the torus axis and `across` away from it,
/// both from the tube's centre circle.
#[derive(Clone)]
pub(crate) struct Tube {
    pub height: Real,
    pub across: Real,
    pub wraps: bool,
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
    /// The low ring's tube angle on a torus band whose rings are not quarter
    /// turns. Its seam then starts on that ring, at this angle.
    pub tube_start: Option<Tube>,
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
    /// A tube angle that is not a quarter turn, which replaces `v`.
    pub(crate) tube: Option<Tube>,
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
        let mut has_blend_torus = false;
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
                            Surface::Torus(_) | Surface::TorusSheet(_)
                                if band_loops(model, face).is_some() =>
                            {
                                has_band = true;
                                has_torus = true;
                                has_blend_torus |= model
                                    .face_origins(face)
                                    .iter()
                                    .any(|origin| origin.role.name() == "blend");
                            }
                            _ => return Ok(Self::default()),
                        }
                    }
                }
            }
        }
        // Unblended torus scenes keep their charted bytes; complete fillet bands may use this route.
        if !has_band || (has_torus && !has_cone && !has_blend_torus) {
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
                            Surface::Torus(_) | Surface::TorusSheet(_) => {
                                let band = match torus_band(model, face, dec)? {
                                    Some(band) => Some(band),
                                    None => tube_band(model, face, dec)?,
                                };
                                let Some(band) = band else {
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
            Surface::Torus(_) | Surface::TorusSheet(_) => source.axis.clone(),
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
                tube_start: source.rings.first().and_then(|ring| ring.tube.clone()),
                v_scale: source.v_scale.clone(),
                radius_bound: match model.geom().surface(source.surface) {
                    Surface::TorusSheet(sheet) => &sheet.major + &sheet.minor,
                    _ => source.rings.iter().fold(
                        source
                            .chain
                            .as_ref()
                            .map_or_else(Real::zero, |(_, _, radius)| radius.clone()),
                        |sum, r| &sum + &r.radius,
                    ),
                },
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
                    v: match (&ring.tube, ring.turns) {
                        (Some(tube), _) => Parameter::Tube(tube.clone()),
                        (None, Some(turns)) => Parameter::QuarterTurns(turns),
                        (None, None) => Parameter::Exact(ring.v.clone()),
                    },
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
        match decided(dec.compare(&rings[0].v, &rings[1].v))? {
            Outcome::Proven(StrictOrder::Greater) => rings.swap(0, 1),
            Outcome::Proven(StrictOrder::Less) => {}
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
            && !matches!(
                decided(dec.sign(&height))?,
                Outcome::Proven(StrictSign::Positive)
            )
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
        let increasing = match decided(dec.sign(&normal.dot(&axis)))? {
            Outcome::Proven(sign) => sign == StrictSign::Positive,
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
            tube: None,
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

/// The two coaxial parallels bounding a torus band face, low ring first, and
/// the low ring's loop index. The band runs up the tube from its low ring, so
/// that ring runs with the chart's u read on the face's side.
pub(crate) fn tube_rings(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<(Vector3, Vec<SourceRing>, usize)>, Refusal> {
    let definition = model.face(face);
    let surface = model.geom().surface(definition.surface);
    let axis = match surface {
        Surface::Torus(torus) => &torus.axis,
        Surface::TorusSheet(torus) => &torus.axis,
        _ => return Ok(None),
    };
    let Some((edges, None)) = band_loops(model, face) else {
        return Ok(None);
    };
    let axis = number::unit(axis);
    let mut rings = Vec::with_capacity(2);
    for edge in edges {
        let Some(ring) = tube_ring(model, surface, &axis, edge, dec)? else {
            return Ok(None);
        };
        rings.push(ring);
    }
    let mut low = None;
    for (index, ring) in rings.iter().enumerate() {
        let along = model.edge(ring.edge).sense == Sense::Same;
        let uses = ring_use(model, face, ring.edge)?;
        if (ring.increasing == (along == uses)) == (definition.sense == Sense::Same)
            && low.replace(index).is_some()
        {
            return Ok(None);
        }
    }
    let Some(low) = low else {
        return Ok(None);
    };
    if low == 1 {
        rings.swap(0, 1);
    }
    Ok(Some((axis, rings, low)))
}

// The band of a torus between two parallels at any tube angles, beside a
// cone. Its high ring is a turn higher when its angle is below the low
// ring's, as on a band across the outer equator.
fn tube_band(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<SourceBand>, Refusal> {
    let Some((axis, mut rings, _)) = tube_rings(model, face, dec)? else {
        return Ok(None);
    };
    let [first, second] = [&rings[0], &rings[1]].map(|ring| ring.tube.as_ref().unwrap());
    let Some(below) = tube_below(second, first, dec)? else {
        return Ok(None);
    };
    let definition = model.face(face);
    let surface = model.geom().surface(definition.surface);
    if matches!(surface, Surface::TorusSheet(_)) {
        for ring in &rings {
            match decided(dec.sign(&ring.radius))? {
                Outcome::Proven(StrictSign::Positive) => {}
                _ => return Ok(None),
            }
        }
        if avoids_inner_equator(first, second, below, dec)? != Some(true) {
            return Ok(None);
        }
    }
    rings[1].tube.as_mut().unwrap().wraps = below;
    let (Surface::Torus(Torus { minor, ref_dir, .. })
    | Surface::TorusSheet(TorusSheet { minor, ref_dir, .. })) = surface
    else {
        unreachable!("tube rings bound a toroidal face")
    };
    Ok(Some(SourceBand {
        face,
        surface: definition.surface,
        axis,
        reference: number::unit(ref_dir),
        v_scale: minor.clone(),
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
    let increasing = match decided(dec.sign(&normal.dot(axis)))? {
        Outcome::Proven(sign) => sign == StrictSign::Positive,
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
        tube: None,
        increasing,
    }))
}

// A coaxial torus parallel at any tube angle, or `None`. Only bands beside a
// cone come here: other torus bands keep their charted bytes.
fn tube_ring(
    model: &Model,
    surface: &Surface,
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
    let (Surface::Torus(Torus {
        center: origin,
        major,
        minor,
        ..
    })
    | Surface::TorusSheet(TorusSheet {
        center: origin,
        major,
        minor,
        ..
    })) = surface
    else {
        return Ok(None);
    };
    let offset = center.minus(origin);
    if !zero(&normal.cross(axis).norm_squared(), dec)?
        || !zero(&offset.cross(axis).norm_squared(), dec)?
    {
        return Ok(None);
    }
    let height = offset.dot(axis);
    let radius = number::sqrt_real(&radius_sq);
    let across = &radius - major;
    if !zero(
        &(&(&height * &height) + &(&across * &across) - minor * minor),
        dec,
    )? {
        return Ok(None);
    }
    let increasing = match decided(dec.sign(&normal.dot(axis)))? {
        Outcome::Proven(sign) => sign == StrictSign::Positive,
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
        turns: None,
        tube: Some(Tube {
            height,
            across,
            wraps: false,
        }),
        increasing,
    }))
}

// Positive endpoint radii stay positive throughout a sheet band only if its
// tube interval avoids the inner equator, where its radius is minimal.
fn avoids_inner_equator(
    first: &Tube,
    second: &Tube,
    wraps: bool,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let inner = Tube {
        height: Real::zero(),
        across: Real::from_i64(-1),
        wraps: false,
    };
    let (Some(before), Some(after)) = (
        tube_below(first, &inner, dec)?,
        tube_below(&inner, second, dec)?,
    ) else {
        return Ok(None);
    };
    let crosses = if wraps {
        before || after
    } else {
        before && after
    };
    Ok(Some(!crosses))
}

// Whether the tube angle of `first` in `[0, 2 pi)` is below that of `second`.
pub(crate) fn tube_below(
    first: &Tube,
    second: &Tube,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let mut half = |tube: &Tube| -> Result<Option<bool>, Refusal> {
        // The upper half turn holds the outer equator but not the inner one.
        Ok(match decided(dec.sign(&tube.height))? {
            Outcome::Proven(sign) => Some(sign == StrictSign::Negative),
            _ => match decided(dec.sign(&tube.across))? {
                Outcome::Proven(sign) => Some(sign == StrictSign::Negative),
                _ => None,
            },
        })
    };
    let (Some(a), Some(b)) = (half(first)?, half(second)?) else {
        return Ok(None);
    };
    if a != b {
        return Ok(Some(b));
    }
    let turn = &(&first.across * &second.height) - &(&first.height * &second.across);
    Ok(match decided(dec.sign(&turn))? {
        Outcome::Proven(sign) => Some(sign == StrictSign::Positive),
        _ => None,
    })
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
    Ok(matches!(
        decided(dec.sign(value))?,
        Outcome::ProvenCoincident(_)
    ))
}

// Band recognition declines only on a proven answer. An undecided one refuses
// with its `numbers.*` code: declining would write the legacy layout, so the
// bytes would depend on the work allowance.
fn decided<T>(outcome: Outcome<T>) -> Result<Outcome<T>, Refusal> {
    match outcome {
        Outcome::NotExact(decision) => Err(Refusal::from(decision)),
        Outcome::Refused(refusal) => Err(refusal),
        decided => Ok(decided),
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Pose, RVector3, Rotation};

    #[test]
    fn sheet_band_intervals_reject_the_inner_equator_and_complement() {
        let tube = |(across, height): (i64, i64)| Tube {
            across: Real::from_i64(across),
            height: Real::from_i64(height),
            wraps: false,
        };
        for (low, high, wraps, expected) in [
            ((1, 1), (0, 1), false, Some(true)),
            ((0, 1), (1, 1), true, Some(false)),
            ((0, -1), (0, 1), true, Some(true)),
            ((0, 1), (0, -1), false, Some(false)),
            ((1, -1), (1, 1), true, Some(true)),
            ((1, 1), (1, -1), false, Some(false)),
            ((-1, 1), (-1, -1), false, Some(false)),
            ((-1, -1), (-1, 1), true, Some(true)),
            ((-1, 0), (0, 1), true, None),
            ((0, 1), (-1, 0), false, None),
        ] {
            assert_eq!(
                avoids_inner_equator(&tube(low), &tube(high), wraps, &mut Decider::default())
                    .unwrap(),
                expected,
                "{low:?} -> {high:?}, wraps={wraps}"
            );
        }
    }

    #[test]
    fn a_filleted_torus_cannot_bypass_a_trimmed_spheres_chart() {
        std::thread::Builder::new()
            .stack_size(8 * 1024 * 1024)
            .spawn(|| {
                let mut case = include_str!("../../../parts/fillets/cases.jsonl")
                    .lines()
                    .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
                    .find(|case| case["id"] == "torus-slice")
                    .unwrap();
                let recipe = catallax_recipe::Recipe::parse(&case["recipe"].to_string()).unwrap();
                let mut dec = Decider::default();
                let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
                case["recipe"]["result"] = "j".into();
                let unblended =
                    catallax_recipe::Recipe::parse(&case["recipe"].to_string()).unwrap();
                let unblended = catallax_features::evaluate(&unblended, &mut dec).unwrap();
                assert!(
                    Plan::prepare(&unblended.model, [unblended.result].into_iter(), &mut dec)
                        .unwrap()
                        .faces
                        .is_empty()
                );
                assert!(
                    !Plan::prepare(&built.model, [built.result].into_iter(), &mut dec)
                        .unwrap()
                        .faces
                        .is_empty()
                );
                let sphere =
                    catallax_topology::sphere(&mut built.model, 5.into(), &Pose::identity())
                        .unwrap();
                let pose = Pose::from_rational(
                    Rotation::identity(),
                    RVector3::new((-10).into(), (-10).into(), 3.into()),
                );
                let block = catallax_topology::block(
                    &mut built.model,
                    [20.into(), 20.into(), 10.into()],
                    &pose,
                )
                .unwrap();
                let cap = catallax_boolean::boolean(
                    &mut built.model,
                    sphere,
                    block,
                    catallax_boolean::BoolOp::Intersect,
                    &mut dec,
                )
                .unwrap()
                .body;
                assert!(built.model.check(cap).is_empty());
                for bodies in [[built.result, cap], [cap, built.result]] {
                    let plan = Plan::prepare(&built.model, bodies.into_iter(), &mut dec).unwrap();
                    assert!(plan.faces.is_empty());
                    assert!(plan.rings.is_empty());
                    assert!(plan.references.is_empty());
                }
            })
            .unwrap()
            .join()
            .unwrap();
    }

    // A smaller work or precision allowance may refuse, but must never write
    // another layout: an undecided band test used to fall back to the legacy
    // charted export. With less work the file is byte-identical. With less
    // precision the last digit of a written decimal may differ, so the
    // entities and their order are compared.
    #[test]
    fn a_smaller_allowance_refuses_or_writes_the_same_layout() {
        fn layout(file: &[u8]) -> Vec<String> {
            String::from_utf8_lossy(file)
                .lines()
                .map(|line| line.split('(').next().unwrap_or(line).to_string())
                .collect()
        }
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                let cases: Vec<serde_json::Value> =
                    include_str!("../../../parts/fillets/cases.jsonl")
                        .lines()
                        .map(|line| serde_json::from_str(line).unwrap())
                        .collect();
                for id in ["cone-top", "torus-slice", "cylinder-rim"] {
                    let case = cases.iter().find(|case| case["id"] == id).unwrap();
                    let recipe =
                        catallax_recipe::Recipe::parse(&case["recipe"].to_string()).unwrap();
                    let built =
                        catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
                    let options = crate::StepOptions::default();
                    let write = |dec: &mut Decider| {
                        let mut bytes = Vec::new();
                        crate::write_step(&built.model, &[built.result], &options, dec, &mut bytes)
                            .map(|_| bytes)
                    };
                    let full = write(&mut Decider::default()).unwrap();
                    let mut refused = 0;
                    for (max_precision_bits, max_work) in [
                        (16_384, 0),
                        (16_384, 10),
                        (16_384, 100),
                        (16_384, 1_000),
                        (16_384, 10_000),
                        (64, 1_000_000_000),
                        (128, 1_000_000_000),
                    ] {
                        let mut dec = Decider::new(catallax_numbers::BudgetPolicy {
                            max_precision_bits,
                            max_work,
                        });
                        let at = format!("{id}: {max_precision_bits} bits, {max_work} work");
                        match write(&mut dec) {
                            Ok(bytes) if max_precision_bits == 16_384 => {
                                assert!(bytes == full, "{at} wrote another file");
                            }
                            Ok(bytes) => {
                                assert_eq!(layout(&bytes), layout(&full), "{at}");
                            }
                            Err(_) => refused += 1,
                        }
                    }
                    assert!(
                        refused > 0,
                        "{id}: the sweep must reach a refusing allowance"
                    );
                }
            })
            .unwrap()
            .join()
            .unwrap();
    }
}
