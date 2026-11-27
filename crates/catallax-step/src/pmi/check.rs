//! Checks of datum targets and position controls against the part's exact
//! geometry, made before they are written and available to a caller that maps
//! imported items onto an imported model.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Surface, Vector3};
use catallax_topology::{FaceId, Loop, Model, Sense};

use super::measure::{
    Shape, face_edges, feature_shape, is_positive, is_zero, is_zero_vector, oriented_plane,
};
use super::target::{self, DatumTarget};
use super::{Feature, ItemKind, Pmi, PositionControl, PositionFrame};

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

/// The normal of a planar face pointing out of the material.
pub(super) fn outward_normal(model: &Model, face: FaceId) -> Option<Vector3> {
    let face = model.face(face);
    let Surface::Plane(plane) = model.geom().surface(face.surface) else {
        return None;
    };
    Some(if face.sense == Sense::Same {
        plane.normal.clone()
    } else {
        plane.normal.scaled(&Real::from_i64(-1))
    })
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

fn check_position(
    model: &Model,
    control: &PositionControl<FaceId>,
    datums: &BTreeMap<&str, &Feature<FaceId>>,
    targeted: &BTreeSet<&str>,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut axes = Vec::with_capacity(control.features.len());
    for feature in &control.features {
        match feature_shape(model, feature, item, dec) {
            Ok(Shape::Axis { direction, .. }) => axes.push(direction),
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
        for ((feature, axis), entry) in control.features.iter().zip(&axes).zip(&zone.entry_faces) {
            let normal = outward_normal(model, *entry)
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
            let edges: BTreeSet<_> = face_edges(model, *entry).collect();
            if !feature
                .iter()
                .any(|face| face_edges(model, *face).any(|edge| edges.contains(&edge)))
            {
                return Err(refuse(
                    item,
                    "a projected zone's entry face does not meet its controlled feature",
                ));
            }
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
    let mut datums = BTreeMap::new();
    let mut targeted = BTreeSet::new();
    for item in &pmi.items {
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
                check_position(model, control, &datums, &targeted, &item.id, dec)?;
            }
            _ => {}
        }
    }
    Ok(())
}
