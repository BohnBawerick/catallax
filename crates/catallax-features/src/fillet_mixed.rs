//! Fillets whose selected convex and concave edges meet at a vertex.
//!
//! Section 4.5.1 forms the concave blends of such an operation first and
//! then rolls the convex sheets over them. At a mixed vertex the concave
//! blend leaves a sharp boundary arc against the vertex's other faces. That
//! arc is tangent to the convex edges it joins, so the convex pass continues
//! along the chain they form and closes the junction with a patch that
//! touches the concave blend. That patch is the vertex's junction sector and
//! takes its `corner.vertex` name.
//!
//! Both passes reuse the ordinary strip construction with its own proofs.
//! Only vertices where every incident edge is selected take this route; a
//! mixed vertex that leaves an edge out keeps the strip construction's
//! refusal.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Real, Refusal};
use catallax_recipe::{EdgeName, EdgeSelection, Origin, Role, StepId, VertexName};
use catallax_topology::{BodyId, EdgeId, Model};

use super::{convex, select, strips_of, uses};
use crate::{ConsumedFace, fillet::FilletDefinition};

/// Whether an origin is a blend this step made.
fn own_blend(origin: &Origin, step: &StepId) -> bool {
    &origin.step == step && matches!(origin.role, Role::Blend { .. })
}

/// Builds the operation in a concave and then a convex pass when some
/// selected vertex mixes the two, or `None` when none does or a mixed vertex
/// leaves an incident edge unselected.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    selected: &[EdgeId],
    radius: &Real,
    step: &StepId,
    earlier: &[FilletDefinition],
    consumed: &mut Vec<ConsumedFace>,
    dec: &mut Decider,
) -> Result<Option<(BodyId, Vec<FilletDefinition>)>, Refusal> {
    let input_uses = uses(model, body);
    let mut convexity = BTreeMap::new();
    for edge in selected {
        let Some(incident) = input_uses.get(edge) else {
            return Ok(None);
        };
        // An edge the strip construction cannot classify keeps its refusal.
        let Ok(kind) = convex(model, *edge, incident, dec) else {
            return Ok(None);
        };
        convexity.insert(*edge, kind);
    }
    // Each concave edge's mixed endpoints, with the names of the faces there
    // besides its own two.
    let mut ends: Vec<(EdgeId, VertexName, BTreeSet<Vec<Origin>>)> = Vec::new();
    for vertex in model.vertices_of_body(body) {
        let incident: Vec<_> = input_uses
            .keys()
            .copied()
            .filter(|edge| {
                let e = model.edge(*edge);
                e.start == Some(vertex) || e.end == Some(vertex)
            })
            .collect();
        let kinds: BTreeSet<_> = incident
            .iter()
            .filter_map(|edge| convexity.get(edge))
            .collect();
        if kinds.len() != 2 {
            continue;
        }
        if incident.iter().any(|edge| !convexity.contains_key(edge)) {
            return Ok(None);
        }
        let Some(name) = model.vertex_name(body, vertex) else {
            return Ok(None);
        };
        for edge in &incident {
            if convexity[edge] {
                continue;
            }
            let own: Vec<_> = input_uses[edge].iter().map(|(f, _)| *f).collect();
            let others = incident
                .iter()
                .flat_map(|e| input_uses[e].iter().map(|(f, _)| *f))
                .filter(|f| !own.contains(f))
                .map(|f| model.face_origins(f).to_vec())
                .collect();
            ends.push((*edge, name.clone(), others));
        }
    }
    if ends.is_empty() {
        return Ok(None);
    }
    let blend_names: BTreeMap<_, _> = selected
        .iter()
        .map(|edge| (*edge, model.edge_name(body, *edge)))
        .collect();
    let concave: Vec<_> = selected
        .iter()
        .copied()
        .filter(|edge| !convexity[edge])
        .collect();
    let convex_names: Vec<_> = selected
        .iter()
        .filter(|edge| convexity[*edge])
        .filter_map(|edge| blend_names[edge].clone())
        .collect();
    let (rounded, mut definitions) = strips_of(
        model, body, concave, radius, true, step, earlier, consumed, dec,
    )?;
    // The boundary arc each concave blend leaves at a mixed vertex: between
    // its blend face and one of the vertex's other faces.
    let mut arcs: Vec<(EdgeName, VertexName)> = Vec::new();
    let mut known = earlier.to_vec();
    known.extend(definitions.iter().cloned());
    let candidates = select(model, rounded, &EdgeSelection::All, true, &known, dec)?;
    let rounded_uses = uses(model, rounded);
    let mut second = Vec::new();
    for edge in candidates {
        let incident = &rounded_uses[&edge];
        if !convex(model, edge, incident, dec)? {
            continue;
        }
        let Some(name) = model.edge_name(rounded, edge) else {
            continue;
        };
        if convex_names.contains(&name) {
            second.push(edge);
            continue;
        }
        let faces: Vec<_> = incident
            .iter()
            .map(|(f, _)| model.face_origins(*f).to_vec())
            .collect();
        let found = ends.iter().find(|(concave, _, others)| {
            let Some(Some(blend)) = blend_names.get(concave) else {
                return false;
            };
            let blend = Origin {
                step: step.clone(),
                role: Role::Blend {
                    edge: blend.clone(),
                },
            };
            faces.iter().any(|f| f == std::slice::from_ref(&blend))
                && faces.iter().any(|f| others.contains(f))
        });
        if let Some((_, vertex, _)) = found {
            arcs.push((name, vertex.clone()));
            second.push(edge);
        }
    }
    let (result, later) = strips_of(
        model, rounded, second, radius, true, step, &known, consumed, dec,
    )?;
    let corner = |origin: &Origin| -> Option<Origin> {
        let Role::Blend { edge } = &origin.role else {
            return None;
        };
        if &origin.step != step || !edge.between().iter().any(|o| own_blend(o, step)) {
            return None;
        }
        let vertex = arcs.iter().find(|(arc, _)| arc == edge)?.1.clone();
        Some(Origin {
            step: step.clone(),
            role: Role::Corner { vertex },
        })
    };
    model.map_face_origins(result, |origin| {
        corner(origin).unwrap_or_else(|| origin.clone())
    });
    for mut definition in later {
        if let Some(renamed) = corner(&definition.origin) {
            definition.origin = renamed;
        }
        definitions.push(definition);
    }
    Ok(Some((result, definitions)))
}
