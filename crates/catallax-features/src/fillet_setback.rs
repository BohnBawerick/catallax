//! Exact station-clipped affected regions for straight plane/plane fillets.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeSelection, FilletSetback, Origin, Role, StepId};
use catallax_shapes::{Circle, Curve, Line, Plane, Surface};
use catallax_topology::{Body, BodyId, BodyKind, Edge, EdgeId, Lump, Model, Sense, Shell};

use super::{convex, resolve, spine, uses, vertex};
use crate::{
    extrude::{face, reverse},
    fillet::{FilletDefinition, hull_clear},
    profile::sign,
};

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-overflow", detail)
}

fn incompatible(detail: &str) -> Refusal {
    Refusal::new("fillet.incompatible-setbacks", detail)
}

fn line_edge(
    model: &mut Model,
    from: catallax_topology::VertexId,
    to: catallax_topology::VertexId,
) -> EdgeId {
    let a = model.geom().position(model.vertex(from).point);
    let b = model.geom().position(model.vertex(to).point);
    let curve = model.geom_mut().add_curve(Curve::Line(Line {
        dir: b.minus(&a),
        origin: a,
    }));
    model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(from),
        end: Some(to),
    })
}

fn reversed(path: Vec<(EdgeId, Sense)>) -> Vec<(EdgeId, Sense)> {
    path.into_iter()
        .rev()
        .map(|(edge, sense)| (edge, reverse(sense)))
        .collect()
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    selection: &EdgeSelection,
    radius: &Real,
    limits: &[FilletSetback],
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    if sign(radius, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "fillet.non-positive-radius",
            "fillet radius must be positive",
        ));
    }
    let input_uses = uses(model, body);
    let selected = resolve(model, selection, &input_uses, dec)?;
    let [edge] = selected.as_slice() else {
        return Err(unsupported(
            "station clipping currently requires one selected straight edge",
        ));
    };
    let edge = *edge;
    let input = model.edge(edge).clone();
    let (Curve::Line(_), Some(start), Some(end)) =
        (model.geom().curve(input.curve), input.start, input.end)
    else {
        return Err(unsupported(
            "station clipping requires a smooth straight edge with two endpoints",
        ));
    };
    let endpoints = [start, end];
    let names = endpoints.map(|v| model.vertex_name(body, v));
    let mut setbacks: [Option<Real>; 2] = [None, None];
    for limit in limits {
        let distance: Real = limit.setback.value().clone().into();
        let resolved = resolve(
            model,
            &EdgeSelection::Refs(vec![limit.edge.clone()]),
            &input_uses,
            dec,
        )?;
        if resolved.as_slice() != [edge] {
            return Err(incompatible(
                "a setback must resolve to exactly the selected input edge",
            ));
        }
        let index = match names
            .each_ref()
            .map(|name| name.as_ref() == Some(&limit.at))
        {
            [true, false] => 0,
            [false, true] => 1,
            [false, false] => {
                return Err(incompatible(
                    "a setback vertex is not an endpoint of its selected input edge",
                ));
            }
            [true, true] => {
                return Err(incompatible(
                    "the setback vertex name matches both input endpoints",
                ));
            }
        };
        if sign(&distance, dec)? != Sign::Positive {
            return Err(incompatible("a setback distance must be positive"));
        }
        if let Some(previous) = &setbacks[index]
            && sign(&(&distance - previous), dec)? != Sign::Zero
        {
            return Err(incompatible(
                "conflicting distances constrain the same endpoint station",
            ));
        }
        setbacks[index] = Some(distance);
    }
    let [Some(first), Some(last)] = setbacks.each_ref() else {
        return Err(unsupported(
            "this station construction needs explicit limits at both endpoints",
        ));
    };
    let positions = endpoints.map(|v| model.geom().position(model.vertex(v).point));
    let along = positions[1].minus(&positions[0]);
    let squared = along.norm_squared();
    if sign(&squared, dec)? != Sign::Positive {
        return Err(incompatible(
            "the selected edge has no positive station interval",
        ));
    }
    let length = super::fillet_envelope::real_sqrt(&squared);
    if sign(&(&length - &(first + last)), dec)? != Sign::Positive {
        return Err(incompatible(
            "the two station planes leave no positive edge run",
        ));
    }
    let direction = along.scaled(&(&Real::one() / &length));
    let stations = [
        positions[0].plus(&direction.scaled(first)),
        positions[1].plus(&direction.scaled(&(-last))),
    ];
    let incident = &input_uses[&edge];
    let [(a, _), (b, _)] = incident.as_slice() else {
        return Err(unsupported(
            "the selected edge must separate two input faces",
        ));
    };
    let faces = [*a, *b];
    let supports = faces.map(|f| model.geom().surface(model.face(f).surface).clone());
    if !supports
        .iter()
        .all(|support| matches!(support, Surface::Plane(_)))
    {
        return Err(unsupported(
            "station-clipped curved contact regions need their own construction",
        ));
    }
    let origin = Origin {
        step: step.clone(),
        role: Role::Blend {
            edge: model.edge_name(body, edge).ok_or_else(|| {
                unsupported("the selected edge needs a persistent origin-pair name")
            })?,
        },
    };
    let definition = FilletDefinition::new(
        origin,
        supports,
        faces.map(|f| model.face(f).sense == Sense::Same),
        model.geom().curve(input.curve).clone(),
        radius.clone(),
        convex(model, edge, incident, dec)?,
        dec,
    )?;
    let hull = definition
        .cut_hull([(&stations[0], &along), (&stations[1], &along)], dec)?
        .ok_or_else(|| unsupported("the contact region has no straight station hull"))?;
    // The exact tangent triangle encloses every normal section, not only its rails.
    // Strict separation preserves every other input face's entire trimmed point set.
    for other in model.faces_of_body(body) {
        if !faces.contains(&other)
            && !hull_clear(&hull, model.geom().surface(model.face(other).surface), dec)?
        {
            return Err(unsupported(
                "the station region is not proved clear of an unselected neighbour",
            ));
        }
    }
    for section in [0, 3] {
        for (offset, face) in faces.iter().enumerate() {
            if !catallax_boolean::face_contains(model, *face, &hull[section + offset + 1], dec)? {
                return Err(incompatible(
                    "a station contact lies outside its selected input face",
                ));
            }
        }
    }
    let turn = sign(
        &hull[1]
            .minus(&hull[0])
            .cross(&hull[2].minus(&hull[0]))
            .dot(&along),
        dec,
    )?;
    if turn == Sign::Zero {
        return Err(unsupported("the station cross-section is degenerate"));
    }
    let vertices = hull.each_ref().map(|p| vertex(model, p.clone()));
    let mut sections = Vec::with_capacity(2);
    for offset in [0, 3] {
        let centre = spine(&definition, &hull[offset])?;
        let first = hull[offset + 1].minus(&centre);
        let last = hull[offset + 2].minus(&centre);
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: centre,
            normal: first.cross(&last),
            radius_sq: radius.square(),
            ref_dir: first,
        }));
        let arc = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(vertices[offset + 1]),
            end: Some(vertices[offset + 2]),
        });
        sections.push([
            line_edge(model, vertices[offset], vertices[offset + 1]),
            arc,
            line_edge(model, vertices[offset + 2], vertices[offset]),
        ]);
    }
    let rails = [0, 1, 2].map(|i| line_edge(model, vertices[i], vertices[i + 3]));
    let mut patches = Vec::with_capacity(5);
    for station in 0..2 {
        let mut path: Vec<_> = sections[station]
            .iter()
            .map(|edge| (*edge, Sense::Same))
            .collect();
        if (station == 0) == (turn == Sign::Positive) {
            path = reversed(path);
        }
        let normal = if station == 0 {
            along.scaled(&Real::from_i64(-1))
        } else {
            along.clone()
        };
        patches.push(face(
            model,
            Surface::Plane(Plane {
                origin: stations[station].clone(),
                normal,
                ref_dir: hull[3 * station + 1].minus(&hull[3 * station]),
            }),
            Sense::Same,
            vec![path],
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::Corner {
                    vertex: names[station].clone().ok_or_else(|| {
                        unsupported("the station endpoint needs a persistent vertex name")
                    })?,
                },
            }],
        ));
    }
    for side in 0..3 {
        let next = (side + 1) % 3;
        let mut path = vec![
            (sections[0][side], Sense::Same),
            (rails[next], Sense::Same),
            (sections[1][side], Sense::Reversed),
            (rails[side], Sense::Reversed),
        ];
        if turn == Sign::Negative {
            path = reversed(path);
        }
        let (surface, sense, origins) = if side == 1 {
            (
                definition.surface.clone(),
                Sense::Reversed,
                vec![definition.origin.clone()],
            )
        } else {
            let source = faces[usize::from(side == 2)];
            (
                model.geom().surface(model.face(source).surface).clone(),
                if definition.convex {
                    model.face(source).sense
                } else {
                    reverse(model.face(source).sense)
                },
                model.face_origins(source).to_vec(),
            )
        };
        patches.push(face(model, surface, sense, vec![path], None, origins));
    }
    let shell = model.add_shell(Shell { faces: patches });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let region = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    let faults = model.check(region);
    if !faults.is_empty() {
        return Err(Refusal::new(
            "fillet.invalid-station-region",
            format!("{faults:?}"),
        ));
    }
    let operation = if definition.convex {
        catallax_boolean::BoolOp::Subtract
    } else {
        catallax_boolean::BoolOp::Unite
    };
    let result = catallax_boolean::boolean(model, body, region, operation, dec)
        .map_err(|failure| failure.refusal)?
        .body;
    Ok((result, vec![definition]))
}
