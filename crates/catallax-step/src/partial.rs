//! Exact export-only topology over regular periodic charts. Spherical callers
//! certify that both chart poles are outside the requested trimmed domain.
//! An interior pole first moves onto its own export-only cap.

use crate::{chart, cut, frames, number, poles, trim, wire};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{Curve, Surface};
use catallax_topology::{BodyId, Coedge, EdgeId, FaceId, Loop, Model, Sense, Shell};
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
/// With `reframe`, every trimmed spherical face is first moved onto its
/// sphere's perpendicular frame.
/// `retry_frame` is set only for a certified interior pole with no clear cap.
pub(crate) fn prepare(
    model: &Model,
    bodies: impl Iterator<Item = BodyId> + Clone,
    options: &crate::StepOptions,
    reframe: bool,
    retry_frame: &mut bool,
    dec: &mut Decider,
) -> Result<Option<Prepared>, Refusal> {
    if !reframe {
        return prepare_charts(model, bodies, options, true, true, retry_frame, dec);
    }
    let bodies: Vec<_> = bodies.collect();
    let reframed = frames::perpendicular(model, &bodies);
    let Some(mut prepared) = prepare_charts(
        &reframed.model,
        reframed.bodies.iter().copied(),
        options,
        true,
        true,
        retry_frame,
        dec,
    )?
    else {
        return Ok(None);
    };
    for (source, moved) in reframed.faces {
        let charted = prepared.face_map.get(&moved).copied().unwrap_or(moved);
        prepared.face_map.insert(source, charted);
    }
    Ok(Some(prepared))
}

fn prepare_charts(
    model: &Model,
    bodies: impl Iterator<Item = BodyId> + Clone,
    options: &crate::StepOptions,
    split_poles: bool,
    boundary_frames: bool,
    retry_frame: &mut bool,
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
                Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Torus(_)
            )
    }) {
        return Ok(None);
    }
    // Classify structure before doing exact membership work. In particular,
    // this is not an arbitrary cone, elliptic-cylinder or torus-trim cutover.
    // Trimmed tori need explicit circular trims and a proven exterior v seam.
    for face in faces() {
        let definition = model.face(face);
        match model.geom().surface(definition.surface) {
            Surface::Plane(_) if unchanged_planar_boundary(model, face) => {}
            Surface::Sphere(_) | Surface::Torus(_) if definition.loops.is_empty() => {}
            Surface::Torus(_) if circular_boundary(model, face) => {}
            Surface::Cylinder(_) | Surface::Sphere(_)
                if definition
                    .loops
                    .iter()
                    .all(|id| matches!(model.loop_(*id), Loop::Coedges(_))) => {}
            _ => return Ok(None),
        }
    }
    let mut geometry = trim::Geometry::new(model);
    let mut selected = BTreeSet::new();
    let mut poles = BTreeMap::new();
    let mut interior = Vec::new();
    // A pole on a trim fits neither chart nor cap in the stored frame.
    let mut boundary_poles = BTreeSet::new();
    for face in faces() {
        if selected.contains(&face) {
            continue;
        }
        let definition = model.face(face);
        match model.geom().surface(definition.surface) {
            Surface::Cylinder(_) => {}
            Surface::Torus(_) if !definition.loops.is_empty() => {}
            Surface::Sphere(sphere) if !definition.loops.is_empty() => {
                let points = *poles.entry(definition.surface).or_insert_with(|| {
                    let radius = number::unit(&sphere.axis).scaled(&sphere.radius);
                    [
                        geometry.point(sphere.center.plus(&radius.scaled(&(-Real::one())))),
                        geometry.point(sphere.center.plus(&radius)),
                    ]
                });
                let mut inside = [false; 2];
                for (index, point) in points.into_iter().enumerate() {
                    match geometry.sphere_contains(model, face, point, dec)? {
                        Some(trim::Location::Outside) => {}
                        Some(trim::Location::Inside) if split_poles => inside[index] = true,
                        Some(trim::Location::Boundary) if boundary_frames => {
                            boundary_poles.insert(definition.surface);
                        }
                        _ => return Ok(None),
                    }
                }
                if boundary_poles.contains(&definition.surface) {
                    continue;
                }
                if inside.contains(&true) {
                    interior.push((face, inside));
                    continue;
                }
            }
            _ => continue,
        }
        selected.insert(face);
    }
    let bodies: Vec<_> = bodies.collect();
    if !boundary_poles.is_empty() {
        let Some(reframed) = frames::reframe(
            model,
            &bodies,
            &boundary_poles,
            split_poles,
            &mut geometry,
            dec,
        )?
        else {
            return Ok(None);
        };
        let Some(mut prepared) = prepare_charts(
            &reframed.model,
            reframed.bodies.iter().copied(),
            options,
            split_poles,
            false,
            retry_frame,
            dec,
        )?
        else {
            return Ok(None);
        };
        for (source, framed) in reframed.faces {
            let charted = prepared.face_map.get(&framed).copied().unwrap_or(framed);
            prepared.face_map.insert(source, charted);
        }
        return Ok(Some(prepared));
    }
    if !interior.is_empty() {
        let Some(split) = poles::split(model, &bodies, &interior, &mut geometry, dec)? else {
            *retry_frame = true;
            return Ok(None);
        };
        let Some(mut prepared) = prepare_charts(
            &split.model,
            split.bodies.iter().copied(),
            options,
            false,
            false,
            retry_frame,
            dec,
        )?
        else {
            return Ok(None);
        };
        prepared.plan.added_edges += split.rings;
        for (source, holed) in split.faces {
            let charted = prepared.face_map.get(&holed).copied().unwrap_or(holed);
            prepared.face_map.insert(source, charted);
        }
        return Ok(Some(prepared));
    }
    if selected.is_empty() {
        return Ok(None);
    }
    let Some(mut prepared) = prepare_regular(model, &bodies, selected, geometry, dec)? else {
        return Ok(None);
    };
    // A touching spherical domain needs simple chart cycles and coupled
    // pcurve fits. Other repeated corners still keep the legacy encoding.
    for (&face, chart) in &prepared.plan.faces {
        if !matches!(
            prepared
                .model
                .geom()
                .surface(prepared.model.face(face).surface),
            Surface::Sphere(_)
        ) {
            continue;
        }
        let mut corners = BTreeSet::new();
        for id in &prepared.model.face(face).loops {
            if chart.touching {
                corners.clear();
            }
            let Loop::Coedges(coedges) = prepared.model.loop_(*id) else {
                continue;
            };
            for id in coedges {
                let coedge = prepared.model.coedge(*id);
                let edge = prepared.model.edge(coedge.edge);
                // Sides follow the curve's natural direction, as in the walks.
                let sides = chart.sides[id];
                let (a, b) = if edge.sense == Sense::Same {
                    (edge.start, edge.end)
                } else {
                    (edge.end, edge.start)
                };
                let (vertex, side) = if coedge.sense == edge.sense {
                    (a, sides[0])
                } else {
                    (b, sides[1])
                };
                if let Some(vertex) = vertex
                    && !corners.insert((vertex, side))
                {
                    return Ok(None);
                }
            }
        }
    }
    // The raw translator cannot resolve separate vertices below its 1e-7
    // model-unit precision. Do not turn a rounded-away sliver into full
    // periodic coverage. Keep the legacy bytes rather than regularizing it.
    // Check every written edge: an uncharted sliver leaves the scene unresolved.
    let resolution = options
        .curve_tolerance
        .abs()
        .max(Rational::from_ratio(1, 10_000_000).unwrap());
    let mut edges = BTreeSet::new();
    for &body in &prepared.bodies {
        for face in prepared.model.faces_of_body(body) {
            for id in &prepared.model.face(face).loops {
                if let Loop::Coedges(coedges) = prepared.model.loop_(*id) {
                    edges.extend(coedges.iter().map(|c| prepared.model.coedge(*c).edge));
                }
            }
        }
    }
    for edge in edges {
        let definition = prepared.model.edge(edge);
        let (Some(start), Some(end)) = (definition.start, definition.end) else {
            assert!(
                !prepared.plan.uses.contains_key(&edge),
                "chart edge without endpoints"
            );
            continue;
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

/// A split edge can leave a loop passing one vertex twice, where a closed
/// trim touches the rest of its boundary. On a reversed planar face the raw
/// translator's wire check rejects such a loop unless it starts at that
/// vertex. Other faces keep their loop start, and with it their bytes.
fn start_at_repeated_vertex(model: &Model, wire: &mut [Coedge]) {
    let starts: Vec<_> = wire
        .iter()
        .map(|coedge| {
            let edge = model.edge(coedge.edge);
            let (a, b) = if edge.sense == Sense::Same {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            if coedge.sense == edge.sense { a } else { b }
        })
        .collect();
    let mut seen = BTreeSet::new();
    let Some(repeated) = starts
        .iter()
        .flatten()
        .find(|vertex| !seen.insert(**vertex))
    else {
        return;
    };
    let first = starts
        .iter()
        .position(|vertex| *vertex == Some(*repeated))
        .expect("a repeated vertex starts a coedge");
    wire.rotate_left(first);
}

// Uncharted caps may gain exact subdivisions, but must keep the same native
// analytic space curves. An opaque intersection could otherwise receive the
// chart fitter's different 3D approximation without a proven chart for this face.
fn unchanged_planar_boundary(model: &Model, face: FaceId) -> bool {
    let definition = model.face(face);
    !definition.loops.is_empty()
        && definition.loops.iter().all(|id| {
            let Loop::Coedges(coedges) = model.loop_(*id) else {
                return false;
            };
            !coedges.is_empty()
                && coedges.iter().all(|id| {
                    match model.geom().curve(model.edge(model.coedge(*id).edge).curve) {
                        Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_) => true,
                        Curve::Intersection(curve) => {
                            curve
                                .surfaces
                                .iter()
                                .any(|id| matches!(model.geom().surface(*id), Surface::Plane(_)))
                                && curve.surfaces.iter().all(|id| {
                                    matches!(
                                        model.geom().surface(*id),
                                        Surface::Plane(_)
                                            | Surface::Cylinder(_)
                                            | Surface::Sphere(_)
                                    )
                                })
                        }
                        _ => false,
                    }
                })
        })
}

// Opaque torus intersections are quartic. Their membership queries are not
// part of this cutover and can exhaust an operation allowance that the legacy
// encoder never spends, so they keep the legacy bytes.
fn circular_boundary(model: &Model, face: FaceId) -> bool {
    let definition = model.face(face);
    !definition.loops.is_empty()
        && definition.loops.iter().all(|id| {
            let Loop::Coedges(coedges) = model.loop_(*id) else {
                return false;
            };
            !coedges.is_empty()
                && coedges.iter().all(|id| {
                    matches!(
                        model.geom().curve(model.edge(model.coedge(*id).edge).curve),
                        Curve::Circle(_)
                    )
                })
        })
}

fn touching_spherical_trims(model: &Model, face: FaceId) -> bool {
    if !matches!(
        model.geom().surface(model.face(face).surface),
        Surface::Sphere(_)
    ) {
        return false;
    }
    let mut starts = BTreeSet::new();
    for id in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*id) else {
            continue;
        };
        for id in coedges {
            let coedge = model.coedge(*id);
            let edge = model.edge(coedge.edge);
            let vertex = if coedge.sense == Sense::Same {
                edge.start
            } else {
                edge.end
            };
            if let Some(vertex) = vertex
                && !starts.insert(vertex)
            {
                return true;
            }
        }
    }
    false
}

pub(crate) fn prepare_regular(
    model: &Model,
    bodies: &[BodyId],
    faces: BTreeSet<FaceId>,
    mut geometry: trim::Geometry,
    dec: &mut Decider,
) -> Result<Option<Prepared>, Refusal> {
    assert!(!faces.is_empty(), "empty periodic preparation");
    let mut cuts = cut::Cuts::default();
    let mut meridians = BTreeMap::new();
    let mut touching = BTreeSet::new();
    for face in faces {
        let coupled = touching_spherical_trims(model, face);
        if coupled {
            touching.insert(face);
        }
        let mut meridian = chart::meridian(&mut geometry, model, face, coupled, dec)?;
        if matches!(
            geometry.store.surface(model.face(face).surface),
            Surface::Torus(_)
        ) && !chart::torus_seam_outside(&mut geometry, model, face, &meridian, dec)?
        {
            return Ok(None);
        }
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
        let mut chart = wire::cut_face(&mut prepared, face, meridian, &splits)?;
        if touching.contains(&face) {
            chart = chart.split_touching(&prepared)?;
        }
        added_edges += chart.seams.len();
        let (replacement, mut chart) = chart.store(&mut prepared, face);
        chart.touching = touching.contains(&face);
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
                        let mut expanded = splits.expand(&prepared, coedges);
                        if source.sense == Sense::Reversed
                            && matches!(model.geom().surface(source.surface), Surface::Plane(_))
                        {
                            start_at_repeated_vertex(&prepared, &mut expanded);
                        }
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
    Ok(Some(Prepared {
        model: prepared,
        bodies: replacements,
        plan,
        face_map,
        edge_map,
    }))
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
