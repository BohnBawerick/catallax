//! RWStepShape and RWStepRepr fields, retaining inherited components.
use super::decode::{boolean, component, reference, references};
use super::*;

fn fields<'a>(
    e: &'a Entity,
    kind: &str,
    simple_offset: usize,
    id: EntityId,
) -> Result<(&'a [crate::part21::Param], usize), ImportRefusal> {
    Ok((
        component(e, kind, id)?,
        if e.parts.len() == 1 { simple_offset } else { 0 },
    ))
}

pub(super) fn read(e: &Entity, id: EntityId) -> Result<Topology, ImportRefusal> {
    let kinds = [
        "BREP_WITH_VOIDS",
        "FACETED_BREP",
        "MANIFOLD_SOLID_BREP",
        "SHELL_BASED_SURFACE_MODEL",
        "FACE_BASED_SURFACE_MODEL",
        "ORIENTED_CLOSED_SHELL",
        "ORIENTED_OPEN_SHELL",
        "CLOSED_SHELL",
        "OPEN_SHELL",
        "CONNECTED_FACE_SET",
        "ORIENTED_FACE",
        "ADVANCED_FACE",
        "FACE_SURFACE",
        "FACE",
        "FACE_OUTER_BOUND",
        "FACE_BOUND",
        "EDGE_LOOP",
        "POLY_LOOP",
        "VERTEX_LOOP",
        "ORIENTED_EDGE",
        "EDGE_CURVE",
        "EDGE",
        "VERTEX_POINT",
    ];
    let kind = kinds
        .into_iter()
        .find(|kind| e.part(kind).is_some())
        .ok_or_else(|| ImportRefusal::Unsupported {
            entity: format!("#{id}"),
            type_name: e.type_name().into(),
        })?;
    let (a, n) = fields(e, kind, 1, id)?;
    Ok(match kind {
        "BREP_WITH_VOIDS" => {
            let outer = if e.parts.len() == 1 {
                reference(a, 1, id)?
            } else {
                reference(component(e, "MANIFOLD_SOLID_BREP", id)?, 0, id)?
            };
            Topology::Solid {
                outer,
                voids: references(a, if e.parts.len() == 1 { 2 } else { 0 }, id)?,
            }
        }
        "FACETED_BREP" | "MANIFOLD_SOLID_BREP" => {
            let (a, n) = if e.parts.len() == 1 {
                (a, n)
            } else {
                fields(e, "MANIFOLD_SOLID_BREP", 1, id)?
            };
            Topology::Solid {
                outer: reference(a, n, id)?,
                voids: Vec::new(),
            }
        }
        "SHELL_BASED_SURFACE_MODEL" => Topology::SurfaceModel(references(a, n, id)?),
        "FACE_BASED_SURFACE_MODEL" => Topology::FaceModel(references(a, n, id)?),
        "ORIENTED_CLOSED_SHELL" | "ORIENTED_OPEN_SHELL" => Topology::OrientedShell(
            reference(a, if e.parts.len() == 1 { 2 } else { 0 }, id)?,
            boolean(a, if e.parts.len() == 1 { 3 } else { 1 }, id)?,
        ),
        "CLOSED_SHELL" | "OPEN_SHELL" | "CONNECTED_FACE_SET" => {
            let (a, n) = if e.parts.len() == 1 || kind == "CONNECTED_FACE_SET" {
                (a, n)
            } else {
                fields(e, "CONNECTED_FACE_SET", 1, id)?
            };
            Topology::Shell {
                faces: references(a, n, id)?,
                closed: match kind {
                    "CLOSED_SHELL" => Some(true),
                    "OPEN_SHELL" => Some(false),
                    _ => None,
                },
            }
        }
        "ORIENTED_FACE" => Topology::OrientedFace(
            reference(a, if e.parts.len() == 1 { 2 } else { 0 }, id)?,
            boolean(a, if e.parts.len() == 1 { 3 } else { 1 }, id)?,
        ),
        "ADVANCED_FACE" | "FACE_SURFACE" | "FACE" => {
            let bounds = if e.parts.len() == 1 || kind == "FACE" {
                references(a, n, id)?
            } else {
                references(component(e, "FACE", id)?, 0, id)?
            };
            let support = if kind == "FACE" {
                None
            } else {
                let (a, n) = if e.parts.len() == 1 {
                    (a, 2)
                } else {
                    fields(e, "FACE_SURFACE", 2, id)?
                };
                Some((reference(a, n, id)?, boolean(a, n + 1, id)?))
            };
            Topology::Face {
                bounds,
                surface: support.map(|s| s.0),
                same_sense: support.map(|s| s.1),
            }
        }
        "FACE_BOUND" | "FACE_OUTER_BOUND" => {
            let (a, n) = if e.parts.len() == 1 || kind == "FACE_BOUND" {
                (a, n)
            } else {
                fields(e, "FACE_BOUND", 1, id)?
            };
            Topology::FaceBound {
                loop_id: reference(a, n, id)?,
                orientation: boolean(a, n + 1, id)?,
                outer: kind == "FACE_OUTER_BOUND",
            }
        }
        "EDGE_LOOP" => Topology::EdgeLoop(references(a, n, id)?),
        "POLY_LOOP" => Topology::PolyLoop(references(a, n, id)?),
        "VERTEX_LOOP" => Topology::VertexLoop(reference(a, n, id)?),
        "ORIENTED_EDGE" => Topology::OrientedEdge(
            reference(a, if e.parts.len() == 1 { 3 } else { 0 }, id)?,
            boolean(a, if e.parts.len() == 1 { 4 } else { 1 }, id)?,
        ),
        "EDGE_CURVE" | "EDGE" => {
            let (ends, offset) = if e.parts.len() == 1 || kind == "EDGE" {
                (a, n)
            } else {
                fields(e, "EDGE", 1, id)?
            };
            let support = if kind == "EDGE" {
                None
            } else {
                let offset = if e.parts.len() == 1 { 3 } else { 0 };
                Some((reference(a, offset, id)?, boolean(a, offset + 1, id)?))
            };
            Topology::Edge {
                vertices: [
                    reference(ends, offset, id)?,
                    reference(ends, offset + 1, id)?,
                ],
                curve: support.map(|s| s.0),
                same_sense: support.map(|s| s.1),
            }
        }
        "VERTEX_POINT" => Topology::Vertex(reference(a, n, id)?),
        _ => unreachable!("selected supported topology"),
    })
}

pub(super) fn placement_link(e: &Entity, id: EntityId) -> Result<PlacementLink, ImportRefusal> {
    for kind in [
        "REPRESENTATION_MAP",
        "MAPPED_ITEM",
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
        "ITEM_DEFINED_TRANSFORMATION",
    ] {
        let Some(a) = e.part(kind) else { continue };
        return Ok(match kind {
            "REPRESENTATION_MAP" => {
                PlacementLink::RepresentationMap(reference(a, 0, id)?, reference(a, 1, id)?)
            }
            "MAPPED_ITEM" => {
                let n = usize::from(e.parts.len() == 1);
                PlacementLink::MappedItem(reference(a, n, id)?, reference(a, n + 1, id)?)
            }
            "ITEM_DEFINED_TRANSFORMATION" => PlacementLink::ItemDefinedTransformation([
                reference(a, 2, id)?,
                reference(a, 3, id)?,
            ]),
            _ => {
                let (base, n, transform) = if e.parts.len() == 1 {
                    (a, 2, 4)
                } else {
                    (component(e, "REPRESENTATION_RELATIONSHIP", id)?, 2, 0)
                };
                PlacementLink::Relationship(
                    [reference(base, n, id)?, reference(base, n + 1, id)?],
                    reference(a, transform, id)?,
                )
            }
        });
    }
    Err(ImportRefusal::Unsupported {
        entity: format!("#{id}"),
        type_name: e.type_name().into(),
    })
}
