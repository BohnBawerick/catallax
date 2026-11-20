//! Retained curved carriers and all physical loops in an approximate artifact.

use super::*;
use catallax_shapes::{BSplineCurve, SplineSupport};
use catallax_topology::numerical::{NumericalFace, NumericalNativeTrim, NumericalTrim};
use std::collections::BTreeMap;

mod native;
mod payload;

pub(super) type PcurveUse<'a> = (usize, &'a NumericalTrim, Option<&'a BSplineCurve>);

pub(super) fn pcurve_uses<'a>(
    faces: &'a [NumericalFace],
    native_trims: &'a [NumericalNativeTrim],
    edge_count: usize,
    dec: &mut Decider,
) -> Result<Vec<Vec<PcurveUse<'a>>>, Refusal> {
    dec.charge_work(edge_count as u64)?;
    let mut uses = vec![Vec::new(); edge_count];
    let mut native = BTreeMap::new();
    for trim in native_trims {
        dec.charge_work(1)?;
        native.insert(
            (trim.face.as_str(), trim.boundary, trim.use_index),
            &trim.curve,
        );
    }
    for (index, face) in faces.iter().enumerate() {
        dec.charge_work(1)?;
        for (boundary_index, boundary) in face.loops.iter().enumerate() {
            dec.charge_work(1)?;
            for (use_index, trim) in boundary.iter().enumerate() {
                dec.charge_work(1)?;
                let curve = native
                    .get(&(face.key.as_str(), boundary_index, use_index))
                    .copied();
                if !matches!(face.surface, Surface::Plane(_)) || curve.is_some() {
                    uses[trim.edge].push((index, trim, curve));
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
    audit: bool,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: Rational::zero(),
        direction_error: Rational::zero(),
        magnitude_error: Rational::zero(),
        audit,
    };
    let payload = payload::build(body, &mut w)?;
    let definition = &payload.definition;
    w.error = w.error.max(payload.scalar_error);
    w.direction_error = w.direction_error.max(payload.direction_error);
    w.magnitude_error = w.magnitude_error.max(payload.magnitude_error);
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
        surfaces.push(match &payload.surfaces[index] {
            payload::SurfaceValues::Plane(values) => plane(&mut w, values, &face.key)?,
            payload::SurfaceValues::Native(entity) => *entity,
            payload::SurfaceValues::Spline => match &face.surface {
                Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
                    native::surface(&mut w, patch.carrier())?
                }
                Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
                    native::surface(&mut w, atlas.source())?
                }
                _ => return Err(bad("unsupported", "unprepared numerical STEP surface")),
            },
        });
    }
    let context = w.graph.add(
        "GEOMETRIC_REPRESENTATION_CONTEXT",
        vec![text("''"), text("''"), text("2")],
    );
    let uses = pcurve_uses(
        &definition.faces,
        &payload.native_trims,
        definition.edges.len(),
        w.dec,
    )?;
    let mut edges = Vec::new();
    let mut parameter = super::analytic::ParameterRounding {
        error: Rational::zero(),
    };
    let mut parameter_checks = Vec::new();
    for (index, edge) in definition.edges.iter().enumerate() {
        let curve = match &payload.curves[index] {
            payload::CurveValues::Line(values) => line(&mut w, values, &edge.key)?,
            payload::CurveValues::Native(entity) => *entity,
            payload::CurveValues::Spline => match &edge.curve {
                Curve::BSpline(spline) => native::finite_curve(&mut w, spline, &edge.parameters)?,
                _ => return Err(bad("unsupported", "unprepared numerical STEP edge")),
            },
        };
        let mut pcurves = Vec::new();
        for (face_index, trim, native_curve) in &uses[index] {
            w.dec.charge_work(1)?;
            if let Some(curve) = native_curve {
                pcurves.push(native::native_pcurve(
                    &mut w,
                    surfaces[*face_index],
                    context,
                    &edge.parameters,
                    curve,
                )?);
                continue;
            }
            let mut uv = trim.uv.clone();
            if trim.sense == Sense::Reversed {
                uv.reverse()
            }
            if matches!(
                definition.faces[*face_index].surface,
                Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Torus(_)
            ) {
                let (curve, actual) = super::analytic::pcurve(
                    &mut w,
                    &mut parameter,
                    surfaces[*face_index],
                    &edge.parameters,
                    &uv,
                )?;
                // An audit measures the u/v line as it was written.
                if w.audit {
                    parameter_checks.push(super::analytic::pcurve_checks(
                        edge,
                        &definition.faces[*face_index],
                        &actual,
                        None,
                        w.dec,
                    )?);
                }
                pcurves.push(curve);
            } else {
                pcurves.push(native::pcurve(
                    &mut w,
                    surfaces[*face_index],
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
                if uses[index].len() == 2 && uses[index][0].0 == uses[index][1].0 {
                    "SEAM_CURVE"
                } else {
                    "SURFACE_CURVE"
                },
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
            let start = if matches!(
                face.surface,
                Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Torus(_)
            ) {
                super::analytic::seam_loop_start(boundary, &definition.edges, face.sense, w.dec)?
            } else {
                0
            };
            let uses = boundary
                .iter()
                .cycle()
                .skip(start)
                .take(boundary.len())
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
    let solid = w.representation(faces, body)?;
    let emitted_native_trims = payload.native_trims.clone();
    let mut checks = payload.phase_checks;
    checks.extend(parameter_checks);
    let mut spline_rounding = payload.spline_rounding;
    spline_rounding.max_parameter_error = spline_rounding.max_parameter_error.max(parameter.error);
    let mut artifact = w.finish(
        solid,
        body,
        Payload::Curved {
            definition: payload.definition,
            native_trims: payload.native_trims,
            checks,
        },
        Some(spline_rounding),
        vec![
            "source.retained-native-carriers-and-physical-loops".into(),
            "file.rounding-error-added-to-tolerance".into(),
            "file.canonical-byte-identity".into(),
            "file.curved-seam-contract".into(),
            "file.approximate-category-no-global-bound".into(),
        ],
    )?;
    artifact.report.emitted_native_trims = emitted_native_trims;
    Ok(artifact)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::Plane;

    #[test]
    fn pcurve_index_charges_skipped_planar_uses_once() {
        let trim = NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
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
            pcurve_uses(&faces, &[], 32, &mut Decider::default())
        });
        assert!(result.unwrap().iter().all(Vec::is_empty));
        assert_eq!(work.used, expected);
        assert!(!work.exhausted);
        let (result, work) = with_operation_budget(expected - 1, || {
            pcurve_uses(&faces, &[], 32, &mut Decider::default())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
