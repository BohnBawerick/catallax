//! Retained curved carriers and all physical loops in an approximate artifact.

use super::*;
use catallax_shapes::SplineSupport;
use catallax_topology::numerical::{NumericalFace, NumericalTrim};

mod native;
mod payload;

pub(super) type PcurveUse<'a> = (usize, &'a NumericalTrim);

pub(super) fn pcurve_uses<'a>(
    faces: &'a [NumericalFace],
    edge_count: usize,
    dec: &mut Decider,
) -> Result<Vec<Vec<PcurveUse<'a>>>, Refusal> {
    dec.charge_work(edge_count as u64)?;
    let mut uses = vec![Vec::new(); edge_count];
    for (index, face) in faces.iter().enumerate() {
        dec.charge_work(1)?;
        for boundary in &face.loops {
            dec.charge_work(1)?;
            for trim in boundary {
                dec.charge_work(1)?;
                if !matches!(face.surface, Surface::Plane(_)) {
                    uses[trim.edge].push((index, trim));
                }
            }
        }
    }
    Ok(uses)
}

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
    let uses = pcurve_uses(&definition.faces, definition.edges.len(), w.dec)?;
    let mut edges = Vec::new();
    for (index, edge) in definition.edges.iter().enumerate() {
        let curve = match &edge.curve {
            Curve::Line(_) => line(
                &mut w,
                payload.lines[index].as_ref().expect("prepared line"),
                &edge.key,
            )?,
            Curve::BSpline(spline) => native::finite_curve(&mut w, spline, &edge.parameters)?,
            _ => return Err(bad("unsupported", "unprepared numerical STEP edge")),
        };
        let mut pcurves = Vec::new();
        for (face_index, trim) in &uses[index] {
            w.dec.charge_work(1)?;
            let mut uv = trim.uv.clone();
            if trim.sense == Sense::Reversed {
                uv.reverse()
            }
            pcurves.push(native::pcurve(
                &mut w,
                surfaces[*face_index],
                context,
                &edge.parameters,
                uv,
            )?);
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::Plane;

    #[test]
    fn pcurve_index_charges_skipped_planar_uses_once() {
        let trim = NumericalTrim {
            edge: 0,
            sense: Sense::Same,
            uv: std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero())),
        };
        let face = NumericalFace {
            key: "planar".into(),
            surface: Surface::Plane(Plane {
                origin: Point3::new(0.into(), 0.into(), 0.into()),
                normal: Vector3::new(0.into(), 0.into(), 1.into()),
                ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
            }),
            sense: Sense::Same,
            loops: vec![vec![trim; 16]],
        };
        let faces = vec![face; 32];
        // Each explicit reservation also pays Decider's cancellation checkpoint.
        let expected = 32 + 1 + 2 * (32 + 32 + 32 * 16);
        let (result, work) = with_operation_budget(expected, || {
            pcurve_uses(&faces, 32, &mut Decider::default())
        });
        assert!(result.unwrap().iter().all(Vec::is_empty));
        assert_eq!(work.used, expected);
        assert!(!work.exhausted);
        let (result, work) = with_operation_budget(expected - 1, || {
            pcurve_uses(&faces, 32, &mut Decider::default())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
