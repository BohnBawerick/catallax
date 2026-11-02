//! Complete line/circular-arc boundary definitions from physical face loops.

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{
    Circle, Surface,
    distance::{PlanarBoundary, PlanarEdge},
};
use catallax_topology::{FaceId, Loop, Model, Sense};
use std::collections::BTreeMap;

pub(super) fn boundary(
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
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            continue;
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            let endpoints = match (edge.start, edge.end) {
                (None, None) => None,
                (Some(a), Some(b)) => {
                    let mut pair = [a, b].map(|vertex| {
                        *indices.entry(vertex).or_insert_with(|| {
                            let index = vertices.len();
                            vertices.push(model.geom().position(model.vertex(vertex).point));
                            index
                        })
                    });
                    if edge.sense == Sense::Reversed {
                        pair.swap(0, 1);
                    }
                    Some(pair)
                }
                _ => {
                    return Err(Refusal::new(
                        "measure.distance-planar-trim",
                        "a finite planar trim edge has only one endpoint",
                    ));
                }
            };
            match curve_form(model.geom(), edge.curve, dec)? {
                CurveForm::Line { .. } => {
                    let Some(pair) = endpoints else {
                        return Err(Refusal::new(
                            "measure.distance-planar-trim",
                            "a planar line trim needs finite endpoints",
                        ));
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
                    ends: endpoints,
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
