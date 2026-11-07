//! Shared-rail construction of affine-radius fillets on convex planar skins.

use super::*;
use catallax_recipe::{RadiusLaw, VariableFilletEdge};
use catallax_shapes::{LinearFillet, Plane};
use catallax_topology::{BodyKind, LinearFilletBand};

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-variable-contact", detail)
}

struct Input {
    edge: EdgeId,
    ends: [VertexId; 2],
    faces: [FaceId; 2],
    senses: [Sense; 2],
    radii: [Real; 2],
}
struct Band {
    input: Input,
    caps: [FaceId; 2],
    support: LinearFillet,
    origin: Origin,
}

// A valid closed, outward-oriented polyhedral skin whose faces all support
// every input vertex is its convex half-space intersection. Convex planar
// loops exclude hidden holes and self-crossing face domains. No tessellation
// or sampled material label supplies any part of this proof.
fn convex_skin(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Vec<(FaceId, Plane, Vector3)>, Refusal> {
    let input = model.body(body);
    let [lump] = input.lumps.as_slice() else {
        return Err(unsupported(
            "variable planar fillets need one convex material lump",
        ));
    };
    if input.kind != BodyKind::Solid
        || model.lump(*lump).shells.len() != 1
        || !model.check(body).is_empty()
    {
        return Err(unsupported(
            "variable planar fillets need one valid closed convex shell",
        ));
    }
    let vertices = model
        .vertices_of_body(body)
        .into_iter()
        .map(|v| model.geom().position(model.vertex(v).point))
        .collect::<Vec<_>>();
    let mut planes = Vec::new();
    for face in model.faces_of_body(body) {
        dec.check_deadline()?;
        let f = model.face(face);
        let Surface::Plane(plane) = model.geom().surface(f.surface) else {
            return Err(unsupported(
                "this variable-radius constructor needs planar input faces",
            ));
        };
        let [lp] = f.loops.as_slice() else {
            return Err(unsupported("convex input faces have one polygon loop"));
        };
        let Loop::Coedges(cs) = model.loop_(*lp) else {
            return Err(unsupported("a convex source face needs a polygon boundary"));
        };
        if cs.len() < 3 {
            return Err(unsupported(
                "a convex source face needs at least three sides",
            ));
        }
        let normal = outward(model, face, &plane.origin, dec)?;
        for point in &vertices {
            if sign(&point.minus(&plane.origin).dot(&normal), dec)? == Sign::Positive {
                return Err(unsupported(
                    "input planes do not prove a convex material shell",
                ));
            }
        }
        let mut polygon = Vec::with_capacity(cs.len());
        for c in cs {
            let c = model.coedge(*c);
            let e = model.edge(c.edge);
            let Curve::Line(line) = model.geom().curve(e.curve) else {
                return Err(unsupported(
                    "a convex planar source needs straight polygon edges",
                ));
            };
            let (Some(a), Some(b)) = edge_blend::ends(model, c.edge, c.sense) else {
                return Err(unsupported("polygon edges need two vertices"));
            };
            let a = model.geom().position(model.vertex(a).point);
            let b = model.geom().position(model.vertex(b).point);
            if sign(&line.dir.norm_squared(), dec)? != Sign::Positive
                || sign(&line.dir.dot(&normal), dec)? != Sign::Zero
                || sign(&line.origin.minus(&plane.origin).dot(&normal), dec)? != Sign::Zero
                || sign(&a.minus(&plane.origin).dot(&normal), dec)? != Sign::Zero
                || sign(&b.minus(&plane.origin).dot(&normal), dec)? != Sign::Zero
                || sign(&a.minus(&line.origin).cross(&line.dir).norm_squared(), dec)? != Sign::Zero
                || sign(&b.minus(&line.origin).cross(&line.dir).norm_squared(), dec)? != Sign::Zero
            {
                return Err(unsupported("source polygon incidence is not exact"));
            }
            polygon.push(a);
        }
        let mut area = Real::zero();
        for i in 0..polygon.len() {
            let a = &polygon[i];
            let b = &polygon[(i + 1) % polygon.len()];
            let along = b.minus(a);
            if sign(&along.norm_squared(), dec)? != Sign::Positive {
                return Err(unsupported("source polygon contains a collapsed edge"));
            }
            let into = normal.cross(&along);
            for p in &polygon {
                if sign(&p.minus(a).dot(&into), dec)? == Sign::Negative {
                    return Err(unsupported(
                        "source face has no proved convex polygon domain",
                    ));
                }
            }
            area = &area
                + &a.minus(&plane.origin)
                    .cross(&b.minus(&plane.origin))
                    .dot(&normal);
        }
        if sign(&area, dec)? != Sign::Positive {
            return Err(unsupported(
                "source face loop has the wrong material orientation",
            ));
        }
        planes.push((face, plane.clone(), normal));
    }
    Ok(planes)
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    laws: &[VariableFilletEdge],
    step: &StepId,
    earlier: &[FilletDefinition],
    consumed: &mut Vec<crate::ConsumedFace>,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    let uses = uses(model, body);
    let mut selected = BTreeSet::new();
    let mut endpoints = BTreeMap::new();
    let mut inputs = Vec::with_capacity(laws.len());
    let mut constant = None;
    let mut same_constant = true;
    for law in laws {
        let found = resolve_refs(model, std::slice::from_ref(&law.edge), &uses, dec)?;
        let [edge] = found.as_slice() else {
            return Err(Refusal::new(
                "fillet.ambiguous-edge",
                "each radius law must resolve to one input edge",
            ));
        };
        if !selected.insert(*edge) {
            return Err(Refusal::new(
                "fillet.conflicting-radius-laws",
                "two radius laws resolve to the same input edge",
            ));
        }
        let e = model.edge(*edge);
        let (Some(a), Some(b)) = (e.start, e.end) else {
            return Err(unsupported(
                "an endpoint law needs an open edge with two endpoints",
            ));
        };
        let matches = [a, b].map(|v| {
            model
                .vertex_name(body, v)
                .is_some_and(|name| name == law.from)
        });
        let ends = match matches {
            [true, false] => [a, b],
            [false, true] => [b, a],
            _ => {
                return Err(Refusal::new(
                    "fillet.invalid-law-endpoint",
                    "from must name exactly one endpoint of its selected edge",
                ));
            }
        };
        let RadiusLaw::Linear { start, end } = &law.radius else {
            return Err(unsupported("the radius law has no construction definition"));
        };
        for (vertex, radius) in ends.iter().copied().zip([start.value(), end.value()]) {
            if let Some(previous) = endpoints.insert(vertex, radius)
                && previous != radius
            {
                return Err(Refusal::new(
                    "fillet.incompatible-radius-laws",
                    "selected laws assign different radii to their shared input vertex",
                ));
            }
        }
        let radii = [
            Real::from(start.value().clone()),
            Real::from(end.value().clone()),
        ];
        if start != end {
            same_constant = false;
        }
        match &constant {
            None => constant = Some(start.clone()),
            Some(radius) if radius != start => same_constant = false,
            _ => (),
        }
        let [(a, sa), (b, sb)] = uses[edge].as_slice() else {
            return Err(bad("a variable fillet edge needs two incident faces"));
        };
        inputs.push(Input {
            edge: *edge,
            ends,
            faces: [*a, *b],
            senses: [*sa, *sb],
            radii,
        });
    }
    if inputs.is_empty() {
        return Err(unsupported(
            "a variable fillet needs a nonempty law selection",
        ));
    }
    if same_constant {
        let refs = EdgeSelection::Refs(laws.iter().map(|law| law.edge.clone()).collect());
        let radius = Real::from(constant.unwrap().value().clone());
        return super::build(
            model, body, &refs, &radius, true, step, earlier, consumed, dec,
        );
    }
    let planes = convex_skin(model, body, dec)?;
    let mut bands = Vec::<Band>::with_capacity(inputs.len());
    for input in inputs {
        if !convex(model, input.edge, &uses[&input.edge], dec)? {
            return Err(unsupported(
                "concave variable bands require a nonconvex material-region proof",
            ));
        }
        let points = input
            .ends
            .map(|v| model.geom().position(model.vertex(v).point));
        let Curve::Line(line) = model.geom().curve(model.edge(input.edge).curve) else {
            return Err(unsupported(
                "linear station laws currently need a straight input edge",
            ));
        };
        if sign(
            &points[1].minus(&points[0]).cross(&line.dir).norm_squared(),
            dec,
        )? != Sign::Zero
        {
            return Err(bad(
                "selected straight edge has inconsistent endpoint incidence",
            ));
        }
        let supports = input.faces.map(|face| {
            planes
                .iter()
                .find(|(f, _, _)| *f == face)
                .unwrap()
                .1
                .clone()
        });
        let same_sense = input
            .faces
            .map(|face| model.face(face).sense == Sense::Same);
        let support =
            LinearFillet::new(supports, same_sense, points, input.radii.clone(), true, dec)?;
        let mut caps = [input.faces[0]; 2];
        for (end, old) in input.ends.iter().enumerate() {
            let incident: BTreeSet<_> = uses
                .iter()
                .filter(|(edge, _)| {
                    let e = model.edge(**edge);
                    e.start == Some(*old) || e.end == Some(*old)
                })
                .flat_map(|(_, faces)| faces.iter().map(|(face, _)| *face))
                .collect();
            if incident.len() != 3 {
                return Err(Refusal::new(
                    UNSUPPORTED_END,
                    "variable band ends require a regular three-face input vertex",
                ));
            }
            caps[end] = *incident
                .iter()
                .find(|face| !input.faces.contains(face))
                .unwrap();
            let plane = &planes
                .iter()
                .find(|(face, _, _)| *face == caps[end])
                .unwrap()
                .1;
            if sign(&plane.normal.cross(support.along()).norm_squared(), dec)? != Sign::Zero {
                return Err(Refusal::new(
                    UNSUPPORTED_END,
                    "an oblique variable-band end needs its exact end-section trimming construction",
                ));
            }
        }
        // Positive weights place every represented point in this exact hull.
        // Strict clearance from nonincident source faces excludes every
        // unhandled roll-over event over the whole band and affected region.
        for (face, plane, normal) in &planes {
            for point in support.carrier().controls() {
                let side = sign(&point.minus(&plane.origin).dot(normal), dec)?;
                if side == Sign::Positive
                    || (!input.faces.contains(face)
                        && !caps.contains(face)
                        && side != Sign::Negative)
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "variable contact region reaches a neighbouring face and needs continuation",
                    ));
                }
            }
        }
        for previous in &bands {
            if !hulls_apart(
                support.carrier().controls(),
                previous.support.carrier().controls(),
                dec,
            )? {
                return Err(Refusal::new(
                    UNSUPPORTED_CORNER,
                    "variable contact regions meet and need a compatible junction construction",
                ));
            }
        }
        let name = model
            .edge_name(body, input.edge)
            .ok_or_else(|| bad("variable fillet source edge has no unambiguous persistent name"))?;
        bands.push(Band {
            input,
            caps,
            support,
            origin: Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            },
        });
    }
    let mut edit = Replacement::default();
    let mut end_changes = BTreeMap::new();
    let mut made = Vec::with_capacity(bands.len());
    for band in bands {
        let input = model.edge(band.input.edge).clone();
        let mut feet = [[band.input.ends[0]; 2]; 2];
        for (end, old) in band.input.ends.iter().copied().enumerate() {
            for (side, rail) in feet.iter_mut().enumerate() {
                rail[end] = if band.support.zero_end(end) {
                    old
                } else {
                    vertex(model, band.support.foot(side, end).clone())
                };
            }
        }
        let named_forward = input.start == Some(band.input.ends[0]);
        let mut path = Vec::with_capacity(4);
        for (side, rail) in feet.iter().enumerate() {
            let a = band.support.foot(side, 0).clone();
            let b = band.support.foot(side, 1).clone();
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: a.clone(),
                dir: b.minus(&a),
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: if named_forward {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                start: Some(rail[usize::from(!named_forward)]),
                end: Some(rail[usize::from(named_forward)]),
            });
            edit.edges
                .insert((band.input.faces[side], band.input.edge), edge);
            path.push((edge, reverse(band.input.senses[side])));
        }
        for (end, old) in band.input.ends.iter().copied().enumerate() {
            if let Some(circle) = band.support.section(end) {
                let curve = model.geom_mut().add_curve(Curve::Circle(circle));
                let cap = model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: Some(feet[0][end]),
                    end: Some(feet[1][end]),
                });
                let first_finishes = if reverse(band.input.senses[0]) == Sense::Same {
                    input.end == Some(old)
                } else {
                    input.start == Some(old)
                };
                path.push((
                    cap,
                    if first_finishes {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                ));
                edit.gaps.insert((band.caps[end], old), Gap::Edge(cap));
            }
            for (edge, adjacent) in &uses {
                if *edge == band.input.edge {
                    continue;
                }
                let e = model.edge(*edge);
                if e.start != Some(old) && e.end != Some(old) {
                    continue;
                }
                let side = if adjacent
                    .iter()
                    .any(|(face, _)| *face == band.input.faces[0])
                {
                    0
                } else {
                    1
                };
                let replacement = feet[side][end];
                if replacement != old {
                    end_changes.insert((*edge, old), replacement);
                }
            }
        }
        let path = close(model, path)?;
        let sense = if band.support.orientation() > 0 {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let face = crate::extrude::face(
            model,
            Surface::BSpline(SplineSupport::LinearFillet(band.support)),
            sense,
            vec![path],
            None,
            vec![band.origin],
        );
        edit.faces.push((band.input.faces[0], face));
        made.push(face);
    }
    for (old, incident) in &uses {
        let e = model.edge(*old).clone();
        let start = e
            .start
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        let end = e
            .end
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        if start == e.start && end == e.end {
            continue;
        }
        let edge = model.add_edge(Edge { start, end, ..e });
        for (face, _) in incident {
            edit.edges.insert((*face, *old), edge);
        }
    }
    let result = edge_blend::apply(model, body, &edit)?;
    for face in made {
        LinearFilletBand::of_face(model, face, dec)?
            .ok_or_else(|| bad("new variable band lost its exact carrier"))?;
    }
    Ok((result, Vec::new()))
}
