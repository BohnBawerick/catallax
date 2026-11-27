//! Checks of datum targets and position controls against the part's exact
//! geometry, made before they are written and available to a caller that maps
//! imported items onto an imported model.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Plane, Point3, Surface, Vector3};
use catallax_topology::{FaceId, Loop, Model, Sense};

use super::measure::{
    Shape, face_edges, feature_shape, is_positive, is_zero, is_zero_vector, oriented_plane,
};
use super::target::{self, DatumTarget};
use super::{Feature, ItemKind, Pmi, PositionControl, PositionFrame};
use crate::graph::{Field, Graph, Id};

fn refuse(item: &str, what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-invalid", format!("{item}: {what}"))
}

/// The complete trim of a planar face from its loops, or `None` when the face
/// is not planar or an edge is neither a line nor a circle.
pub(super) fn face_boundary(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<PlanarBoundary>, Refusal> {
    let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
        return Ok(None);
    };
    let mut vertices = Vec::new();
    let mut indices = BTreeMap::new();
    let mut edges = Vec::new();
    for boundary in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*boundary) else {
            continue;
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            let ends = match (edge.start, edge.end) {
                (None, None) => None,
                (Some(a), Some(b)) => {
                    let mut pair = [a, b].map(|vertex| {
                        *indices.entry(vertex).or_insert_with(|| {
                            vertices.push(model.geom().position(model.vertex(vertex).point));
                            vertices.len() - 1
                        })
                    });
                    if edge.sense == Sense::Reversed {
                        pair.swap(0, 1);
                    }
                    Some(pair)
                }
                _ => return Ok(None),
            };
            match curve_form(model.geom(), edge.curve, dec)? {
                CurveForm::Line { .. } => {
                    let Some(pair) = ends else {
                        return Ok(None);
                    };
                    edges.push(PlanarEdge::Segment(pair));
                }
                CurveForm::Circle {
                    center,
                    normal,
                    radius_sq,
                    ref_dir,
                } => edges.push(PlanarEdge::Arc {
                    circle: Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    },
                    // A circle closed at one vertex is the whole circle.
                    ends: ends.filter(|[a, b]| a != b),
                }),
                _ => return Ok(None),
            }
        }
    }
    Ok((!edges.is_empty()).then(|| PlanarBoundary {
        plane: plane.clone(),
        vertices,
        edges,
    }))
}

/// A planar face and its normal pointing out of the material.
pub(super) fn outward_plane(model: &Model, face: FaceId) -> Option<(&Plane, Vector3)> {
    let face = model.face(face);
    let Surface::Plane(plane) = model.geom().surface(face.surface) else {
        return None;
    };
    let normal = if face.sense == Sense::Same {
        plane.normal.clone()
    } else {
        plane.normal.scaled(&Real::from_i64(-1))
    };
    Some((plane, normal))
}

fn check_target(
    model: &Model,
    target: &DatumTarget<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let boundary = face_boundary(model, target.face, dec)?.ok_or_else(|| {
        refuse(
            item,
            "a datum target needs a planar face bounded by lines and circles",
        )
    })?;
    match target::off_face(&boundary, target, dec)? {
        Some(reason) => Err(refuse(item, reason)),
        None => Ok(()),
    }
}

fn fields<'a>(graph: &'a Graph, field: &Field, name: &str) -> Option<&'a [Field]> {
    let Field::Ref(id) = field else { return None };
    let (kind, fields) = graph.record(*id);
    (kind == name).then_some(fields)
}

fn written_real(value: &str) -> Real {
    let value = value.replace(".E", ".0E");
    let value = if value.ends_with('.') {
        format!("{value}0")
    } else {
        value
    };
    Real::from_rational(crate::number::formatted_value(&value))
}

fn values(fields: &[Field]) -> Option<[Real; 3]> {
    let values = match fields.get(1)? {
        Field::Text(values) => values
            .strip_prefix('(')?
            .strip_suffix(')')?
            .split(',')
            .collect::<Vec<_>>(),
        Field::List(values) => values
            .iter()
            .map(|v| match v {
                Field::Text(text) => Some(text.as_ref()),
                _ => None,
            })
            .collect::<Option<Vec<_>>>()?,
        _ => return None,
    };
    let values: [&str; 3] = values.try_into().ok()?;
    Some(values.map(written_real))
}

fn written_plane(graph: &Graph, placement: &Field) -> Option<Plane> {
    let frame = fields(graph, placement, "AXIS2_PLACEMENT_3D")?;
    let [x, y, z] = values(fields(graph, frame.get(1)?, "CARTESIAN_POINT")?)?;
    let origin = Point3::new(x, y, z);
    let [x, y, z] = values(fields(graph, frame.get(2)?, "DIRECTION")?)?;
    let normal = Vector3::new(x, y, z);
    let [x, y, z] = values(fields(graph, frame.get(3)?, "DIRECTION")?)?;
    let reference = Vector3::new(x, y, z);
    let ref_dir = normal.cross(&reference.cross(&normal));
    Some(Plane {
        origin,
        normal,
        ref_dir,
    })
}

fn written_boundary(graph: &Graph, face: Id) -> Option<(PlanarBoundary, bool)> {
    let ("ADVANCED_FACE", [_, Field::List(bounds), surface, Field::Text(sense)]) =
        graph.record(face)
    else {
        return None;
    };
    let plane = written_plane(graph, fields(graph, surface, "PLANE")?.get(1)?)?;
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for bound in bounds {
        let cycle = fields(graph, bound, "FACE_BOUND")
            .or_else(|| fields(graph, bound, "FACE_OUTER_BOUND"))?
            .get(1)?;
        if fields(graph, cycle, "VERTEX_LOOP").is_some() {
            continue;
        }
        let Field::List(oriented) = fields(graph, cycle, "EDGE_LOOP")?.get(1)? else {
            return None;
        };
        for edge in oriented {
            let edge = fields(
                graph,
                fields(graph, edge, "ORIENTED_EDGE")?.get(3)?,
                "EDGE_CURVE",
            )?;
            let [
                _,
                Field::Ref(start),
                Field::Ref(end),
                curve,
                Field::Text(same),
            ] = edge
            else {
                return None;
            };
            let ends = if start == end {
                None
            } else {
                let mut ends = [0, 0];
                for (slot, vertex) in ends.iter_mut().zip([*start, *end]) {
                    let (_, fields) = graph.record(vertex);
                    let point = fields.get(1)?;
                    let [x, y, z] = values(self::fields(graph, point, "CARTESIAN_POINT")?)?;
                    *slot = vertices.len();
                    vertices.push(Point3::new(x, y, z));
                }
                if same.as_ref() == ".F." {
                    ends.swap(0, 1);
                }
                Some(ends)
            };
            let curve = fields(graph, curve, "SEAM_CURVE")
                .or_else(|| fields(graph, curve, "SURFACE_CURVE"))
                .and_then(|fields| fields.get(1))
                .unwrap_or(curve);
            if fields(graph, curve, "LINE").is_some() {
                edges.push(PlanarEdge::Segment(ends?));
            } else {
                let [_, placement, Field::Text(radius)] = fields(graph, curve, "CIRCLE")? else {
                    return None;
                };
                let circle = written_plane(graph, placement)?;
                edges.push(PlanarEdge::Arc {
                    circle: Circle {
                        center: circle.origin,
                        normal: circle.normal,
                        radius_sq: written_real(radius).square(),
                        ref_dir: circle.ref_dir,
                    },
                    ends,
                });
            }
        }
    }
    (!edges.is_empty()).then_some((
        PlanarBoundary {
            plane,
            vertices,
            edges,
        },
        sense.as_ref() == ".T.",
    ))
}

pub(super) fn check_written_target(
    graph: &Graph,
    faces: &[Id],
    target: &DatumTarget<FaceId>,
    normal: &Vector3,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let [face] = faces else {
        return Err(refuse(item, "the datum target's face was split on export"));
    };
    let (boundary, same) = written_boundary(graph, *face)
        .ok_or_else(|| refuse(item, "the datum target's written trim cannot be checked"))?;
    let plane = &boundary.plane;
    if is_zero_vector(&plane.ref_dir, dec)?
        || !is_zero_vector(&normal.cross(&plane.normal), dec)?
        || is_positive(&normal.dot(&plane.normal), dec)? != same
    {
        return Err(refuse(
            item,
            "the datum target's direction differs from the written face",
        ));
    }
    for vertex in &boundary.vertices {
        if !is_zero(&vertex.minus(&plane.origin).dot(&plane.normal), dec)? {
            return Err(refuse(
                item,
                "rounding changed the datum target's written trim plane",
            ));
        }
    }
    for edge in &boundary.edges {
        if let PlanarEdge::Arc { circle, ends } = edge {
            if is_zero_vector(&circle.ref_dir, dec)?
                || !is_zero_vector(&circle.normal.cross(&plane.normal), dec)?
                || !is_zero(&circle.center.minus(&plane.origin).dot(&plane.normal), dec)?
            {
                return Err(refuse(
                    item,
                    "rounding changed the datum target's written circular trim",
                ));
            }
            if let Some(ends) = ends {
                for end in ends {
                    let d = boundary.vertices[*end].minus(&circle.center).norm_squared();
                    if !is_zero(&(&d - &circle.radius_sq), dec)? {
                        return Err(refuse(
                            item,
                            "rounding changed the datum target's written circular trim",
                        ));
                    }
                }
            }
        }
    }
    if let Some(reason) = target::off_face(&boundary, target, dec)? {
        return Err(refuse(item, format!("on the written face, {reason}")));
    }
    Ok(())
}

/// Whether a datum's feature is a feature of size: one cylinder, or two
/// opposed parallel planes.
fn feature_of_size(
    model: &Model,
    feature: &Feature<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if feature.iter().all(|face| {
        matches!(
            model.geom().surface(model.face(*face).surface),
            Surface::Cylinder(_)
        )
    }) {
        return Ok(matches!(
            feature_shape(model, feature, item, dec),
            Ok(Shape::Axis { .. })
        ));
    }
    let [a, b] = feature.as_slice() else {
        return Ok(false);
    };
    let (Ok(a), Ok(b)) = (
        oriented_plane(model, &[*a], item, dec),
        oriented_plane(model, &[*b], item, dec),
    ) else {
        return Ok(false);
    };
    Ok(is_zero_vector(&a.normal.cross(&b.normal), dec)?
        && is_positive(&(-a.normal.dot(&b.normal)), dec)?
        && !is_zero(&b.point.minus(&a.point).dot(&a.normal), dec)?)
}

fn check_open_entry(
    model: &Model,
    entry: FaceId,
    at: &Point3,
    normal: &Vector3,
    written: Option<&BTreeMap<FaceId, Vec<Id>>>,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let bodies = written.is_none().then(|| model.bodies_with_face(entry));
    let mut pending = vec![entry];
    let mut visited = BTreeSet::new();
    while let Some(face) = pending.pop() {
        if !visited.insert(face)
            || written.is_some_and(|faces| !faces.contains_key(&face))
            || bodies.as_ref().is_some_and(|bodies| {
                face != entry
                    && !model
                        .bodies_with_face(face)
                        .iter()
                        .any(|body| bodies.contains(body))
            })
        {
            continue;
        }
        let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
            continue;
        };
        if !is_zero_vector(&plane.normal.cross(normal), dec)?
            || !is_zero(&at.minus(&plane.origin).dot(normal), dec)?
        {
            continue;
        }
        let boundary = face_boundary(model, face, dec)?
            .ok_or_else(|| refuse(item, "the projected zone entry trim cannot be checked"))?;
        if boundary.contains(at, dec)? {
            return Err(refuse(item, "a projected zone's axis meets a closed end"));
        }
        pending.extend(
            face_edges(model, face)
                .flat_map(|edge| model.coedges_of_edge(edge))
                .filter_map(|coedge| model.owner_of_coedge(coedge).map(|(_, face)| face)),
        );
    }
    Ok(())
}

fn check_inward_trim(
    model: &Model,
    edges: &BTreeSet<catallax_topology::EdgeId>,
    origin: &Point3,
    normal: &Vector3,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let unknown = || refuse(item, "the cylinder's axial trim cannot be checked");
    let mut inward = false;
    let mut check = |point: &Point3, dec: &mut Decider| -> Result<(), Refusal> {
        let distance = point.minus(origin).dot(normal);
        if is_positive(&distance, dec)? {
            return Err(refuse(
                item,
                "the projected zone points back into its controlled cylinder",
            ));
        }
        inward |= is_positive(&(-distance), dec)?;
        Ok(())
    };
    for edge in edges {
        let edge = model.edge(*edge);
        match curve_form(model.geom(), edge.curve, dec)? {
            CurveForm::Line { .. } => {
                for vertex in [edge.start, edge.end] {
                    let vertex = vertex.ok_or_else(unknown)?;
                    check(&model.geom().position(model.vertex(vertex).point), dec)?;
                }
            }
            CurveForm::Circle {
                center,
                normal: axis,
                ..
            } if is_zero_vector(&axis.cross(normal), dec)? => {
                check(&center, dec)?;
            }
            _ => return Err(unknown()),
        }
    }
    if !inward {
        return Err(unknown());
    }
    Ok(())
}

fn check_position(
    model: &Model,
    control: &PositionControl<FaceId>,
    datums: &BTreeMap<&str, &Feature<FaceId>>,
    targeted: &BTreeSet<&str>,
    item: &str,
    written: Option<&BTreeMap<FaceId, Vec<Id>>>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut axes = Vec::with_capacity(control.features.len());
    for feature in &control.features {
        match feature_shape(model, feature, item, dec) {
            Ok(Shape::Axis {
                point, direction, ..
            }) => axes.push((point, direction)),
            Ok(Shape::Plane { .. }) | Err(_) => {
                return Err(refuse(
                    item,
                    "a position control's features are each one cylinder",
                ));
            }
        }
    }
    for segment in control.frame.segments() {
        for datum in &segment.datums {
            if datum.material.is_none() {
                continue;
            }
            let letter = datum.letter.as_str();
            if targeted.contains(letter) {
                return Err(refuse(
                    item,
                    format!(
                        "datum {letter} is set up by datum targets and has no material boundary"
                    ),
                ));
            }
            let eligible = match datums.get(letter) {
                Some(feature) => feature_of_size(model, feature, item, dec)?,
                None => false,
            };
            if !eligible {
                return Err(refuse(
                    item,
                    format!("datum {letter} is not a feature of size and has no material boundary"),
                ));
            }
        }
    }
    if let PositionFrame::Single {
        projected_zone: Some(zone),
        ..
    } = &control.frame
    {
        let direction = target::vector(&zone.direction);
        if zone.entry_faces.len() != control.features.len() {
            return Err(refuse(
                item,
                "a projected zone needs one entry face per feature",
            ));
        }
        let entries: BTreeMap<_, BTreeSet<_>> = zone
            .entry_faces
            .iter()
            .map(|entry| (*entry, face_edges(model, *entry).collect()))
            .collect();
        for ((feature, (point, axis)), entry) in
            control.features.iter().zip(&axes).zip(&zone.entry_faces)
        {
            let (plane, normal) = outward_plane(model, *entry)
                .ok_or_else(|| refuse(item, "a projected zone's entry face is planar"))?;
            if !is_zero_vector(&axis.cross(&normal), dec)? {
                return Err(refuse(
                    item,
                    "a controlled axis is not square to its entry face",
                ));
            }
            if !is_zero_vector(&direction.cross(&normal), dec)?
                || !is_positive(&direction.dot(&normal), dec)?
            {
                return Err(refuse(
                    item,
                    "the projected zone does not point out of the part from its entry face",
                ));
            }
            let edges: BTreeSet<_> = feature
                .iter()
                .flat_map(|face| face_edges(model, *face))
                .collect();
            let mut meeting = entries.iter().filter(|(_, own)| !own.is_disjoint(&edges));
            if meeting.next().map(|(face, _)| face) != Some(entry) || meeting.next().is_some() {
                return Err(refuse(
                    item,
                    "each controlled feature needs exactly one projected zone entry face it meets",
                ));
            }
            let offset = plane.origin.minus(point).dot(&normal) / axis.dot(&normal);
            check_open_entry(
                model,
                *entry,
                &point.plus(&axis.scaled(&offset)),
                &normal,
                written,
                item,
                dec,
            )?;
            check_inward_trim(model, &edges, &plane.origin, &normal, item, dec)?;
        }
    }
    Ok(())
}

/// Checks the datum targets and position controls of `pmi` against `model`:
/// each target lies on its trimmed planar face, each controlled feature is a
/// cylinder, a material boundary is only on a datum feature of size that no
/// target sets up, and a projected zone starts on a face that meets its
/// feature square to its axis and points out of the part. Refuses
/// `step.pmi-invalid` with the item's name. Other items are not checked here.
pub fn check_on_model(model: &Model, pmi: &Pmi<FaceId>, dec: &mut Decider) -> Result<(), Refusal> {
    check_on_model_for_export(model, pmi, None, dec)
}

pub(super) fn check_on_model_for_export(
    model: &Model,
    pmi: &Pmi<FaceId>,
    written: Option<&BTreeMap<FaceId, Vec<Id>>>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut datums = BTreeMap::new();
    let mut targeted = BTreeSet::new();
    for item in &pmi.items {
        if item
            .kind
            .faces()
            .into_iter()
            .any(|face| !model.has_face(*face))
        {
            return Err(refuse(
                &item.id,
                "a face it is attached to is not in the model",
            ));
        }
        match &item.kind {
            ItemKind::Datum(datum) => {
                datums.insert(datum.letter.as_str(), &datum.feature);
            }
            ItemKind::DatumTarget(target) => {
                targeted.insert(target.datum.as_str());
            }
            _ => {}
        }
    }
    for item in &pmi.items {
        match &item.kind {
            ItemKind::DatumTarget(target) => check_target(model, target, &item.id, dec)?,
            ItemKind::Position(control) => {
                check_position(model, control, &datums, &targeted, &item.id, written, dec)?;
            }
            _ => {}
        }
    }
    Ok(())
}
