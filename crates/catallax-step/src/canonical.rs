//! Export-only analytic carriers for the rings of a torus face. A ring that a
//! torus meets a coaxial cylinder or plane in is stored as an intersection
//! component, a quartic to every membership query, though it is a circle. Its
//! chart cuts need that circle, which `curve_form` already writes. The circle
//! replaces the component on a copy of the model, with the same incidences.

use std::collections::BTreeMap;

use catallax_intersect::CurveForm;
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Refusal};
use catallax_shapes::{Circle, Curve, GeomKey, Surface};
use catallax_topology::{BodyId, Coedge, Edge, EdgeId, Face, FaceId, Loop, Model, Shell};

use crate::poles::Split;

/// The model with every vertexless torus ring that is an intersection
/// component and a circle moved onto an explicit circle, or `None` when no
/// ring is.
pub(crate) fn circles(
    model: &Model,
    bodies: &[BodyId],
    dec: &mut Decider,
) -> Result<Option<Split>, Refusal> {
    let mut seen = std::collections::BTreeSet::new();
    let mut rings: Vec<EdgeId> = Vec::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            if !seen.insert(face) {
                continue;
            }
            let definition = model.face(face);
            if !matches!(
                model.geom().surface(definition.surface),
                Surface::Torus(_) | Surface::TorusSheet(_)
            ) {
                continue;
            }
            for loop_id in &definition.loops {
                let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                    continue;
                };
                for coedge in coedges {
                    let edge = model.coedge(*coedge).edge;
                    let source = model.edge(edge);
                    if source.start.is_none()
                        && source.end.is_none()
                        && matches!(model.geom().curve(source.curve), Curve::Intersection(_))
                        && !rings.contains(&edge)
                    {
                        rings.push(edge);
                    }
                }
            }
        }
    }
    let mut replaced = model.clone();
    let mut edges = BTreeMap::<EdgeId, EdgeId>::new();
    for edge in rings {
        let source = model.edge(edge).clone();
        let Curve::Intersection(component) = model.geom().curve(source.curve).clone() else {
            unreachable!("a ring was listed as an intersection component")
        };
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } = catallax_intersect::curve_form(model.geom(), source.curve, dec)?
        else {
            continue;
        };
        let circle = replaced.geom_mut().add_curve(Curve::Circle(Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        }));
        for surface in component.surfaces {
            replaced.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(circle),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("step.chart-ring-circle"),
            });
        }
        edges.insert(
            edge,
            replaced.add_edge(Edge {
                curve: circle,
                ..source
            }),
        );
    }
    if edges.is_empty() {
        return Ok(None);
    }
    let mut faces = BTreeMap::<FaceId, FaceId>::new();
    seen.clear();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            if !seen.insert(face) {
                continue;
            }
            let definition = model.face(face);
            let mut changed = false;
            let mut loops = Vec::with_capacity(definition.loops.len());
            for loop_id in &definition.loops {
                let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                    loops.push(*loop_id);
                    continue;
                };
                if !coedges
                    .iter()
                    .any(|coedge| edges.contains_key(&model.coedge(*coedge).edge))
                {
                    loops.push(*loop_id);
                    continue;
                }
                changed = true;
                let coedges = coedges
                    .iter()
                    .map(|coedge| {
                        let source = model.coedge(*coedge);
                        replaced.add_coedge(Coedge {
                            edge: edges.get(&source.edge).copied().unwrap_or(source.edge),
                            sense: source.sense,
                        })
                    })
                    .collect();
                loops.push(replaced.add_loop(Loop::Coedges(coedges)));
            }
            if changed {
                faces.insert(
                    face,
                    replaced.add_face(Face {
                        loops,
                        ..definition.clone()
                    }),
                );
            }
        }
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let shell = model
                    .shell(*id)
                    .faces
                    .iter()
                    .map(|face| faces.get(face).copied().unwrap_or(*face))
                    .collect();
                *id = replaced.add_shell(Shell { faces: shell });
            }
            *id = replaced.add_lump(lump);
        }
        replacements.push(replaced.add_body(body));
    }
    Ok(Some(Split {
        model: replaced,
        bodies: replacements,
        faces,
        rings: 0,
    }))
}
