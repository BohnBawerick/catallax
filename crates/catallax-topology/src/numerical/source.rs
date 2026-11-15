//! Explicit source conversion. Admission and numerical evidence are still required.

use super::*;
use crate::VertexId;
use std::collections::BTreeMap;

fn scalar(value: &catallax_numbers::Real, dec: &Decider) -> Result<Rational, Refusal> {
    Ok(dec.enclose_at(value, BITS)?.mid())
}

/// Generates finite numerical input from a retained line/patch model. Vertex and
/// parameter values are rounded to 128-bit midpoints; the original carriers are
/// retained. This returns no admitted body or incidence authority. Call [`build`]
/// to check the generated topology, finite trims and residuals under a tolerance.
pub fn definition_from_model(
    model: &Model,
    body: BodyId,
    name: String,
    dec: &mut Decider,
) -> Result<NumericalDefinition, Refusal> {
    let vertices = model.vertices_of_body(body);
    let vertex_map: BTreeMap<_, _> = vertices
        .iter()
        .enumerate()
        .map(|(i, id)| (*id, i))
        .collect();
    let points: BTreeMap<_, _> = vertices
        .iter()
        .map(|id| (*id, model.geom().position(model.vertex(*id).point)))
        .collect();
    let mut result = NumericalDefinition {
        name,
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    for id in &vertices {
        dec.charge_work(1)?;
        let p = &points[id];
        result.vertices.push(NumericalVertex {
            key: format!("{id:?}"),
            position: [scalar(&p.x, dec)?, scalar(&p.y, dec)?, scalar(&p.z, dec)?],
        });
    }
    let edges = model.edges_of_body(body);
    let edge_map: BTreeMap<_, _> = edges.iter().enumerate().map(|(i, id)| (*id, i)).collect();
    for id in &edges {
        dec.charge_work(1)?;
        let edge = model.edge(*id);
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(bad(
                "unsupported",
                "numerical source conversion requires finite edge ends",
            ));
        };
        let curve = model.geom().curve(edge.curve).clone();
        let mut parameters = Vec::new();
        for vertex in [a, b] {
            let p = &points[&vertex];
            let at = match &curve {
                Curve::Line(line) => scalar(
                    &(p.minus(&line.origin).dot(&line.dir) / line.dir.norm_squared()),
                    dec,
                )?,
                Curve::BSpline(curve) => {
                    let at = curve.unique_parameter_of(p, dec)?.ok_or_else(|| {
                        bad("definition", "source vertex has no finite spline parameter")
                    })?;
                    let [lo, hi] = curve.domain();
                    scalar(at.value(), dec)?.max(lo.clone()).min(hi.clone())
                }
                _ => {
                    return Err(bad(
                        "unsupported",
                        "numerical source conversion requires line or native spline edges",
                    ));
                }
            };
            parameters.push(at);
        }
        result.edges.push(NumericalEdge {
            key: format!("{id:?}"),
            curve,
            vertices: [vertex_map[&a], vertex_map[&b]],
            parameters: parameters.try_into().unwrap(),
        });
    }
    for id in model.faces_of_body(body) {
        dec.charge_work(1)?;
        let face = model.face(id);
        let surface = match model.geom().surface(face.surface) {
            Surface::Plane(plane) => Surface::Plane(plane.clone()),
            Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
                Surface::BSpline(SplineSupport::RationalPatch(patch.clone()))
            }
            Surface::BSpline(SplineSupport::RationalAtlas(atlas)) if atlas.charts().len() == 1 => {
                Surface::BSpline(SplineSupport::RationalPatch(
                    atlas.charts()[0].patch().clone(),
                ))
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical source conversion requires planes or one retained rational patch per face",
                ));
            }
        };
        let mut parameters = BTreeMap::<VertexId, [Rational; 2]>::new();
        let mut loops = Vec::new();
        for lp in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                return Err(bad(
                    "unsupported",
                    "numerical source conversion requires explicit finite boundary loops",
                ));
            };
            let mut boundary = Vec::new();
            for use_id in uses {
                let coedge = model.coedge(*use_id);
                let edge = model.edge(coedge.edge);
                let mut ends = [edge.start.unwrap(), edge.end.unwrap()];
                if coedge.sense == Sense::Reversed {
                    ends.reverse();
                }
                let mut uv = Vec::new();
                for vertex in ends {
                    if let Some(at) = parameters.get(&vertex) {
                        uv.push(at.clone());
                        continue;
                    }
                    let point = &points[&vertex];
                    let at = match &surface {
                        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
                            let at = patch.parameters_at(point, dec)?.ok_or_else(|| {
                                bad("trim", "source vertex has no retained patch parameter")
                            })?;
                            [scalar(&at[0], dec)?, scalar(&at[1], dec)?]
                        }
                        Surface::Plane(_) => {
                            let p = [
                                dec.enclose_at(&point.x, BITS)?,
                                dec.enclose_at(&point.y, BITS)?,
                                dec.enclose_at(&point.z, BITS)?,
                            ];
                            curved::parameters(&surface, &p, dec)?
                        }
                        _ => unreachable!(),
                    };
                    parameters.insert(vertex, at.clone());
                    uv.push(at);
                }
                boundary.push(NumericalTrim {
                    edge: edge_map[&coedge.edge],
                    sense: coedge.sense,
                    uv: uv.try_into().unwrap(),
                });
            }
            loops.push(boundary);
        }
        result.faces.push(NumericalFace {
            key: format!("{id:?}"),
            surface,
            sense: face.sense,
            loops,
        });
    }
    Ok(result)
}
