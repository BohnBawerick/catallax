//! Exact export-only topology over regular periodic charts. Spherical callers
//! certify that both chart poles are outside the requested trimmed domain.

use crate::{band, chart, cut, number, trim, wire};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Shell};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) struct Prepared {
    pub model: Model,
    pub bodies: Vec<BodyId>,
    pub plan: chart::Plan,
    face_map: BTreeMap<FaceId, FaceId>,
    edge_map: BTreeMap<EdgeId, Vec<EdgeId>>,
}

/// Preserve the legacy encoder for scenes outside the regular-chart contract,
/// and for complete bands already handled without changing their STEP bytes.
pub(crate) fn prepare(
    model: &Model,
    bodies: impl Iterator<Item = BodyId> + Clone,
    options: &crate::StepOptions,
    dec: &mut Decider,
) -> Result<Option<Prepared>, Refusal> {
    let faces = || {
        bodies
            .clone()
            .flat_map(|body| model.body(body).lumps.iter())
            .flat_map(|lump| model.lump(*lump).shells.iter())
            .flat_map(|shell| model.shell(*shell).faces.iter().copied())
    };
    if !faces().any(|face| {
        let definition = model.face(face);
        !definition.loops.is_empty()
            && matches!(
                model.geom().surface(definition.surface),
                Surface::Cylinder(_) | Surface::Sphere(_)
            )
    }) {
        return Ok(None);
    }
    // Classify structure before doing exact membership work. In particular,
    // this is not an arbitrary cone, torus, or elliptic-cylinder cutover.
    for face in faces() {
        let definition = model.face(face);
        match model.geom().surface(definition.surface) {
            Surface::Plane(_) => {}
            Surface::Sphere(_) | Surface::Torus(_) if definition.loops.is_empty() => {}
            Surface::Cylinder(_) | Surface::Sphere(_)
                if definition
                    .loops
                    .iter()
                    .all(|id| matches!(model.loop_(*id), Loop::Coedges(_))) => {}
            Surface::Cone(_) if band::is_complete(model, face, dec)? => {}
            _ => return Ok(None),
        }
    }
    let mut geometry = trim::Geometry::new(model);
    let mut selected = BTreeSet::new();
    let mut poles = BTreeMap::new();
    for face in faces() {
        if selected.contains(&face) {
            continue;
        }
        let definition = model.face(face);
        match model.geom().surface(definition.surface) {
            Surface::Cylinder(_) | Surface::Cone(_) => {}
            Surface::Sphere(sphere) if !definition.loops.is_empty() => {
                let points = *poles.entry(definition.surface).or_insert_with(|| {
                    let radius = number::unit(&sphere.axis).scaled(&sphere.radius);
                    [
                        geometry.point(sphere.center.plus(&radius.scaled(&(-Real::one())))),
                        geometry.point(sphere.center.plus(&radius)),
                    ]
                });
                for point in points {
                    if geometry.sphere_contains(model, face, point, dec)? != trim::Location::Outside
                    {
                        return Ok(None);
                    }
                }
            }
            _ => continue,
        }
        selected.insert(face);
    }
    if selected.is_empty() {
        return Ok(None);
    }
    let bodies: Vec<_> = bodies.collect();
    let mut prepared = prepare_regular(model, &bodies, selected, geometry, dec)?;
    // The raw translator cannot resolve separate vertices below its 1e-7
    // model-unit precision. Do not turn a rounded-away sliver into full
    // periodic coverage. Keep the legacy bytes rather than regularizing it.
    let resolution = options
        .curve_tolerance
        .abs()
        .max(Rational::from_ratio(1, 10_000_000).unwrap());
    for &edge in prepared.plan.uses.keys() {
        let definition = prepared.model.edge(edge);
        let (Some(start), Some(end)) = (definition.start, definition.end) else {
            unreachable!("chart edge without endpoints")
        };
        for vertex in [start, end] {
            if let std::collections::btree_map::Entry::Vacant(entry) =
                prepared.plan.written_vertices.entry(vertex)
            {
                let point = number::vertex(&prepared.model, vertex, options, dec)?;
                if point.error() > options.curve_tolerance.abs() {
                    return Ok(None);
                }
                entry.insert(point);
            }
        }
        if start == end {
            continue;
        }
        let a = &prepared.plan.written_vertices[&start];
        let b = &prepared.plan.written_vertices[&end];
        let distance_squared = a
            .xyz
            .iter()
            .zip(&b.xyz)
            .fold(Rational::zero(), |sum, (a, b)| {
                let delta = &a.value - &b.value;
                sum + (&delta * &delta)
            });
        let margin = resolution.clone() + a.error() + b.error();
        if distance_squared <= &margin * &margin {
            return Ok(None);
        }
    }
    Ok(Some(prepared))
}

pub(crate) fn prepare_regular(
    model: &Model,
    bodies: &[BodyId],
    faces: BTreeSet<FaceId>,
    mut geometry: trim::Geometry,
    dec: &mut Decider,
) -> Result<Prepared, Refusal> {
    assert!(!faces.is_empty(), "empty periodic preparation");
    let mut cuts = cut::Cuts::default();
    let mut meridians = BTreeMap::new();
    for face in faces {
        let mut meridian = chart::meridian(&mut geometry, model, face, dec)?;
        cuts.add(&geometry, &mut meridian, dec)?;
        meridians.insert(face, meridian);
    }
    cuts.rings(&mut geometry, model, meridians.keys().copied());
    let splits = cuts.prepare(&mut geometry, model, dec)?;
    let mut prepared = model.clone();
    *prepared.geom_mut() = geometry.store;
    let splits = splits.apply(&mut prepared);
    let mut added_edges = splits.added_edges;
    let mut face_map = BTreeMap::new();
    let mut charts = BTreeMap::new();
    for (face, meridian) in meridians {
        let chart = wire::cut_face(&mut prepared, face, meridian, &splits)?;
        added_edges += chart.seams.len();
        let (replacement, chart) = chart.store(&mut prepared, face);
        face_map.insert(face, replacement);
        charts.insert(replacement, chart);
    }
    let mut seen = BTreeSet::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            if !seen.insert(face) || face_map.contains_key(&face) {
                continue;
            }
            let source = model.face(face);
            let changed = source.loops.iter().any(|id| match model.loop_(*id) {
                Loop::Coedges(coedges) => coedges
                    .iter()
                    .any(|c| splits.replacements.contains_key(&model.coedge(*c).edge)),
                Loop::Vertex(_) => false,
            });
            if !changed {
                continue;
            }
            let mut face_definition = source.clone();
            let mut loops = Vec::with_capacity(source.loops.len());
            for id in &source.loops {
                loops.push(match model.loop_(*id) {
                    Loop::Coedges(coedges) => {
                        let expanded = splits.expand(&prepared, coedges);
                        let coedges = expanded
                            .into_iter()
                            .map(|coedge| prepared.add_coedge(coedge))
                            .collect();
                        prepared.add_loop(Loop::Coedges(coedges))
                    }
                    Loop::Vertex(_) => *id,
                });
            }
            face_definition.loops = loops;
            face_map.insert(face, prepared.add_face(face_definition));
        }
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let faces = model
                    .shell(*id)
                    .faces
                    .iter()
                    .map(|face| face_map.get(face).copied().unwrap_or(*face))
                    .collect();
                *id = prepared.add_shell(Shell { faces });
            }
            *id = prepared.add_lump(lump);
        }
        replacements.push(prepared.add_body(body));
    }
    let plan = chart::Plan::new(&prepared, charts, added_edges, splits.added_vertices);
    let edge_map = splits
        .replacements
        .into_iter()
        .map(|(edge, replacements)| (edge, replacements.into_iter().map(|c| c.edge).collect()))
        .collect();
    Ok(Prepared {
        model: prepared,
        bodies: replacements,
        plan,
        face_map,
        edge_map,
    })
}

impl Prepared {
    pub fn errors(
        &self,
        faces: BTreeMap<FaceId, Rational>,
        edges: BTreeMap<EdgeId, Rational>,
    ) -> (BTreeMap<FaceId, Rational>, BTreeMap<EdgeId, Rational>) {
        let faces = faces
            .into_iter()
            .map(|(face, error)| (self.face_map.get(&face).copied().unwrap_or(face), error))
            .collect();
        let mut remapped = BTreeMap::new();
        for (edge, error) in edges {
            if let Some(replacements) = self.edge_map.get(&edge) {
                for &replacement in replacements {
                    remapped.insert(replacement, error.clone());
                }
            } else {
                remapped.insert(edge, error);
            }
        }
        (faces, remapped)
    }
}
