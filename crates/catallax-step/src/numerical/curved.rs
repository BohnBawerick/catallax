//! Retained curved carriers and all physical loops in an approximate artifact.

use super::*;
use catallax_shapes::SplineSupport;

mod native;
mod payload;

fn direction(w: &mut Writer<'_>, values: &[Rational; 3]) -> Result<Id, Refusal> {
    let tokens = values
        .iter()
        .map(|q| native::token(w, q))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(w.graph
        .add("DIRECTION", vec![text("''"), Field::List(tokens)]))
}

fn plane(w: &mut Writer<'_>, values: &payload::PlaneValues, name: &str) -> Result<Id, Refusal> {
    let point = w.point(&values.origin)?;
    let normal = direction(w, &values.normal)?;
    let reference = direction(w, &values.reference)?;
    let axis = w.graph.add(
        "AXIS2_PLACEMENT_3D",
        vec![text("''"), point.into(), normal.into(), reference.into()],
    );
    Ok(w.graph.add("PLANE", vec![text(string(name)), axis.into()]))
}

fn line(w: &mut Writer<'_>, values: &payload::LineValues, name: &str) -> Result<Id, Refusal> {
    let point = w.point(&values.origin)?;
    let direction = direction(w, &values.direction)?;
    let magnitude = w.scalar(&values.magnitude)?;
    let vector = w
        .graph
        .add("VECTOR", vec![text("''"), direction.into(), magnitude]);
    Ok(w.graph.add(
        "LINE",
        vec![text(string(name)), point.into(), vector.into()],
    ))
}

pub(super) fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    let payload = payload::build(body, options, dec)?;
    let definition = payload.body.definition();
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: payload.scalar_error,
        direction_error: payload.direction_error,
    };
    let mut vertices = Vec::new();
    for vertex in &definition.vertices {
        let point = w.point(&vertex.position)?;
        vertices.push(w.graph.add(
            "VERTEX_POINT",
            vec![text(string(&vertex.key)), point.into()],
        ));
    }
    let mut surfaces = Vec::new();
    for (index, face) in definition.faces.iter().enumerate() {
        surfaces.push(match &face.surface {
            Surface::Plane(_) => plane(
                &mut w,
                payload.planes[index].as_ref().expect("prepared plane"),
                &face.key,
            )?,
            Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
                native::surface(&mut w, patch.carrier())?
            }
            _ => return Err(bad("unsupported", "unprepared numerical STEP surface")),
        });
    }
    let context = w.graph.add(
        "GEOMETRIC_REPRESENTATION_CONTEXT",
        vec![text("''"), text("''"), text("2")],
    );
    let mut edges = Vec::new();
    for (index, edge) in definition.edges.iter().enumerate() {
        let curve = match &edge.curve {
            Curve::Line(_) => line(
                &mut w,
                payload.lines[index].as_ref().expect("prepared line"),
                &edge.key,
            )?,
            Curve::BSpline(spline) => native::curve(&mut w, spline)?,
            _ => return Err(bad("unsupported", "unprepared numerical STEP edge")),
        };
        let mut pcurves = Vec::new();
        for (face_index, face) in definition.faces.iter().enumerate() {
            if matches!(face.surface, Surface::Plane(_)) {
                continue;
            }
            for trim in face
                .loops
                .iter()
                .flatten()
                .filter(|trim| trim.edge == index)
            {
                let mut uv = trim.uv.clone();
                if trim.sense == Sense::Reversed {
                    uv.reverse()
                }
                pcurves.push(native::pcurve(
                    &mut w,
                    surfaces[face_index],
                    context,
                    &edge.parameters,
                    uv,
                )?);
            }
        }
        let geometry = if pcurves.is_empty() {
            curve
        } else {
            w.graph.add(
                "SURFACE_CURVE",
                vec![text("''"), curve.into(), refs(pcurves), text(".CURVE_3D.")],
            )
        };
        edges.push(w.graph.add(
            "EDGE_CURVE",
            vec![
                text(string(&edge.key)),
                vertices[edge.vertices[0]].into(),
                vertices[edge.vertices[1]].into(),
                geometry.into(),
                boolean(edge.parameters[0] < edge.parameters[1]),
            ],
        ));
    }
    let mut faces = Vec::new();
    for (index, face) in definition.faces.iter().enumerate() {
        let mut bounds = Vec::new();
        for (number, boundary) in face.loops.iter().enumerate() {
            let uses = boundary
                .iter()
                .map(|trim| {
                    w.graph.add(
                        "ORIENTED_EDGE",
                        vec![
                            text("''"),
                            text("*"),
                            text("*"),
                            edges[trim.edge].into(),
                            boolean(trim.sense == Sense::Same),
                        ],
                    )
                })
                .collect::<Vec<_>>();
            let loop_ = w.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]);
            bounds.push(w.graph.add(
                if number == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                },
                vec![text("''"), loop_.into(), boolean(true)],
            ));
        }
        faces.push(w.graph.add(
            "ADVANCED_FACE",
            vec![
                text(string(&face.key)),
                refs(bounds),
                surfaces[index].into(),
                boolean(face.sense == Sense::Same),
            ],
        ));
    }
    let shell = w.graph.add("CLOSED_SHELL", vec![text("''"), refs(faces)]);
    let solid = w.graph.add(
        "MANIFOLD_SOLID_BREP",
        vec![text("'numerical approximate'"), shell.into()],
    );
    w.finish(
        solid,
        body,
        payload.body.evidence().clone(),
        Some(payload.spline_rounding),
        vec![
            "source.retained-native-carriers-and-physical-loops".into(),
            "file.decimal-model-numerical-admission".into(),
            "file.canonical-byte-identity".into(),
            "file.curved-seam-contract".into(),
            "file.approximate-category-no-global-bound".into(),
        ],
    )
}
