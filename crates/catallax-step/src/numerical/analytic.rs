//! Approximate finite analytic exchange with explicit cylindrical pcurves.
//! Decimal carrier and pcurve values are rechecked as numerical geometry.

use super::*;
use catallax_numbers::{Ball, NumericalCheck, NumericalSample};
use catallax_shapes::{Circle, Cylinder, Line, Plane};
use catallax_topology::numerical::{self as n, NumericalDefinition};

const BITS: u32 = 128;

pub(super) fn applies(body: &NumericalBody) -> bool {
    let d = body.definition();
    d.edges
        .iter()
        .all(|edge| matches!(edge.curve, Curve::Line(_) | Curve::Circle(_)))
        && d.faces
            .iter()
            .all(|face| matches!(face.surface, Surface::Plane(_) | Surface::Cylinder(_)))
        && (d
            .edges
            .iter()
            .any(|edge| matches!(edge.curve, Curve::Circle(_)))
            || d.faces
                .iter()
                .any(|face| matches!(face.surface, Surface::Cylinder(_))))
}

fn round3(p: &[Rational; 3], options: &NumericalStepOptions) -> Result<[Rational; 3], Refusal> {
    Ok([
        written(&p[0], options)?,
        written(&p[1], options)?,
        written(&p[2], options)?,
    ])
}
fn vector_from(p: &[Rational; 3]) -> Vector3 {
    Vector3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}
struct WrittenFrame {
    entity: Id,
    origin: [Rational; 3],
    normal: [Rational; 3],
    reference: [Rational; 3],
}
fn placement(
    w: &mut Writer<'_>,
    origin: &Point3,
    normal: &Vector3,
    reference: &Vector3,
) -> Result<WrittenFrame, Refusal> {
    let o = rational_point(origin)?;
    let n = rational_vector(normal)?;
    let x = rational_vector(reference)?;
    let point = w.point(&o)?;
    let axis = w.direction(&n)?;
    let reference = w.direction(&x)?;
    let id = w.graph.add(
        "AXIS2_PLACEMENT_3D",
        vec![text("''"), point.into(), axis.into(), reference.into()],
    );
    let n = direction_values(&n, w.options, w.dec)?.values;
    let x = direction_values(&x, w.options, w.dec)?.values;
    // gp_Ax2 projects its reference direction before normalization.
    let y = charged_cross(&n, &x, w.dec)?;
    let x = charged_cross(&y, &n, w.dec)?;
    Ok(WrittenFrame {
        entity: id,
        origin: round3(&o, w.options)?,
        normal: n,
        reference: x,
    })
}

fn near(mut value: Rational, reference: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    let half = Ball::pi(BITS).mid();
    let period = &half * &Rational::from_i64(2);
    for _ in 0..4 {
        let delta = arithmetic(&value, reference, '-', dec)?;
        if delta > half {
            value = arithmetic(&value, &period, '-', dec)?;
        } else if delta < -&half {
            value = arithmetic(&value, &period, '+', dec)?;
        } else {
            return Ok(value);
        }
    }
    Err(bad(
        "unsupported",
        "analytic parameter requires a distant periodic chart",
    ))
}
fn angle(x: &Ball, y: &Ball, reference: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    if y.contains_zero() && x.hi().sign() == Sign::Negative {
        // The known source chart selects the periodic branch. A midpoint here
        // remains numerical and is checked against the actual written carrier.
        if y.lo().is_zero() && y.hi().is_zero() {
            return near(Ball::pi(BITS).mid(), reference, dec);
        }
    }
    let value = Ball::atan2(y, x)
        .ok_or_else(|| bad("rounding-topology", "written inverse reaches its axis"))?;
    let at = if value.hi() - value.lo() > Rational::one() {
        if x.hi().sign() != Sign::Negative || !y.contains_zero() {
            return Err(bad("rounding-topology", "written angle is unresolved"));
        }
        Ball::pi(BITS).mid()
    } else {
        value.mid()
    };
    near(at, reference, dec)
}
fn local(
    origin: &Point3,
    axis: &Vector3,
    reference: &Vector3,
    p: &[Rational; 3],
) -> Result<[Ball; 3], Refusal> {
    Ok(crate::pcurve::Placement::new(
        rational_point(origin)?,
        rational_vector(axis)?,
        rational_vector(reference)?,
    )?
    .point(p, BITS))
}
fn uv(
    surface: &Surface,
    p: &[Rational; 3],
    old: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    match surface {
        Surface::Plane(plane) => {
            let [x, y, _] = local(&plane.origin, &plane.normal, &plane.ref_dir, p)?;
            Ok([x.mid(), y.mid()])
        }
        Surface::Cylinder(cylinder) => {
            let [x, y, z] = local(&cylinder.origin, &cylinder.axis, &cylinder.ref_dir, p)?;
            Ok([angle(&x, &y, &old[0], dec)?, z.mid()])
        }
        _ => Err(bad(
            "unsupported",
            "analytic writer requires planes or cylinders",
        )),
    }
}

fn rounded_surfaces(w: &mut Writer<'_>, d: &mut NumericalDefinition) -> Result<Vec<Id>, Refusal> {
    let mut result = Vec::new();
    for face in &mut d.faces {
        w.dec.charge_work(1)?;
        let (id, surface) = match &face.surface {
            Surface::Plane(plane) => {
                let WrittenFrame {
                    entity: frame,
                    origin: o,
                    normal: n,
                    reference: x,
                } = placement(w, &plane.origin, &plane.normal, &plane.ref_dir)?;
                (
                    w.graph
                        .add("PLANE", vec![text(string(&face.key)), frame.into()]),
                    Surface::Plane(Plane {
                        origin: point_from(&o),
                        normal: vector_from(&n),
                        ref_dir: vector_from(&x),
                    }),
                )
            }
            Surface::Cylinder(cylinder) => {
                let WrittenFrame {
                    entity: frame,
                    origin: o,
                    normal: n,
                    reference: x,
                } = placement(w, &cylinder.origin, &cylinder.axis, &cylinder.ref_dir)?;
                let radius = crate::number::sqrt_real(&cylinder.radius_sq)
                    .as_rational()
                    .ok_or_else(|| {
                        bad("unsupported", "analytic writer requires a rational radius")
                    })?;
                let value = w.scalar(&radius)?;
                let radius = written(&radius, w.options)?;
                (
                    w.graph.add(
                        "CYLINDRICAL_SURFACE",
                        vec![text(string(&face.key)), frame.into(), value],
                    ),
                    Surface::Cylinder(Cylinder {
                        origin: point_from(&o),
                        axis: vector_from(&n),
                        radius_sq: Real::from(&radius * &radius),
                        ref_dir: vector_from(&x),
                    }),
                )
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "analytic writer requires planes or cylinders",
                ));
            }
        };
        face.surface = surface;
        result.push(id);
    }
    Ok(result)
}

fn rounded_curves(w: &mut Writer<'_>, d: &mut NumericalDefinition) -> Result<Vec<Id>, Refusal> {
    let mut result = Vec::new();
    for edge in &mut d.edges {
        w.dec.charge_work(1)?;
        let (id, curve, parameters) = match &edge.curve {
            Curve::Line(line) => {
                let origin = rational_point(&line.origin)?;
                let direction = rational_vector(&line.dir)?;
                let point = w.point(&origin)?;
                let axis = w.direction(&direction)?;
                let vector = w
                    .graph
                    .add("VECTOR", vec![text("''"), axis.into(), text("1.")]);
                let id = w.graph.add(
                    "LINE",
                    vec![text(string(&edge.key)), point.into(), vector.into()],
                );
                let o = round3(&origin, w.options)?;
                let v = direction_values(&direction, w.options, w.dec)?.values;
                let norm = crate::number::sqrt(&dot(&v, &v, w.dec)?);
                let dir = Vector3::new(
                    &Real::from(v[0].clone()) / &norm,
                    &Real::from(v[1].clone()) / &norm,
                    &Real::from(v[2].clone()) / &norm,
                );
                let mut ts = Vec::new();
                for index in edge.vertices {
                    let delta = difference(&d.vertices[index].position, &o, w.dec)?;
                    let value = &Real::from(dot(&delta, &v, w.dec)?) / &norm;
                    ts.push(w.dec.enclose_budgeted_at(&value, BITS)?.mid());
                }
                (
                    id,
                    Curve::Line(Line {
                        origin: point_from(&o),
                        dir,
                    }),
                    ts.try_into().unwrap(),
                )
            }
            Curve::Circle(circle) => {
                // A full circle's 3D carrier has no distinguished reference
                // phase. Put its retained endpoint at phase zero so a STEP
                // reader reconstructs the same natural range as its pcurves.
                // The center, axis, radius and source vertex remain unchanged.
                let reference = if edge.vertices[0] == edge.vertices[1] {
                    let center = rational_point(&circle.center)?;
                    vector_from(&difference(
                        &d.vertices[edge.vertices[0]].position,
                        &center,
                        w.dec,
                    )?)
                } else {
                    circle.ref_dir.clone()
                };
                let WrittenFrame {
                    entity: frame,
                    origin: o,
                    normal: n,
                    reference: x,
                } = placement(w, &circle.center, &circle.normal, &reference)?;
                let radius = crate::number::sqrt_real(&circle.radius_sq)
                    .as_rational()
                    .ok_or_else(|| {
                        bad("unsupported", "analytic writer requires a rational radius")
                    })?;
                let value = w.scalar(&radius)?;
                let radius = written(&radius, w.options)?;
                let id = w
                    .graph
                    .add("CIRCLE", vec![text(string(&edge.key)), frame.into(), value]);
                let circle = Circle {
                    center: point_from(&o),
                    normal: vector_from(&n),
                    radius_sq: Real::from(&radius * &radius),
                    ref_dir: vector_from(&x),
                };
                let mut ts = Vec::new();
                for (i, index) in edge.vertices.into_iter().enumerate() {
                    let [a, b, _] = local(
                        &circle.center,
                        &circle.normal,
                        &circle.ref_dir,
                        &d.vertices[index].position,
                    )?;
                    ts.push(angle(&a, &b, &edge.parameters[i], w.dec)?);
                }
                // Keep the edge's actual oriented finite arc, including a seam
                // crossing; changing a periodic branch does not change geometry.
                let period = Ball::pi(BITS).mid() * Rational::from_i64(2);
                if edge.parameters[0] < edge.parameters[1] && ts[1] <= ts[0] {
                    ts[1] = &ts[1] + &period;
                }
                if edge.parameters[0] > edge.parameters[1] && ts[1] >= ts[0] {
                    ts[1] = &ts[1] - &period;
                }
                (id, Curve::Circle(circle), ts.try_into().unwrap())
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "analytic writer requires lines or circular arcs",
                ));
            }
        };
        edge.curve = curve;
        edge.parameters = parameters;
        result.push(id);
    }
    Ok(result)
}

struct ParameterRounding {
    // Lexical selector error across parameter coefficients with different units.
    // It is no common-unit geometric bound; actual 3D distance samples, rather
    // than this diagnostic maximum, enforce the source tolerance.
    error: Rational,
}

fn closed_circle_parameters(
    definition: &mut NumericalDefinition,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let period = Ball::pi(BITS).mid() * Rational::from_i64(2);
    let mut reversed = Vec::new();
    for (index, edge) in definition.edges.iter_mut().enumerate() {
        dec.charge_work(1)?;
        if edge.vertices[0] != edge.vertices[1] || !matches!(edge.curve, Curve::Circle(_)) {
            continue;
        }
        if (&edge.parameters[1] - &edge.parameters[0]).abs() != period {
            return Err(bad(
                "rounding-topology",
                "closed circle must retain one complete turn",
            ));
        }
        // A STEP EDGE_CURVE does not carry this parameter range. OCCT projects
        // the shared endpoint and reconstructs the ascending natural full turn
        // of the anchored carrier above. Its pcurve must use that same range.
        // Changing an edge's direction also changes every coedge sense,
        // retaining each physical boundary walk.
        if edge.parameters[1] < edge.parameters[0] {
            reversed.push(index);
        }
        edge.parameters = [Rational::zero(), period.clone()];
    }
    for face in &mut definition.faces {
        dec.charge_work(1)?;
        for trim in face.loops.iter_mut().flatten() {
            dec.charge_work(1)?;
            if reversed.binary_search(&trim.edge).is_ok() {
                trim.sense = match trim.sense {
                    Sense::Same => Sense::Reversed,
                    Sense::Reversed => Sense::Same,
                };
            }
        }
    }
    Ok(())
}
impl ParameterRounding {
    fn scalar(&mut self, q: &Rational, w: &mut Writer<'_>) -> Result<(Field, Rational), Refusal> {
        w.dec.charge_work(1)?;
        let value = crate::number::rational(q, w.options.significant_digits);
        self.error = self.error.clone().max(value.error);
        Ok((text(value.text), written(q, w.options)?))
    }
}

fn pcurve(
    w: &mut Writer<'_>,
    rounding: &mut ParameterRounding,
    surface: Id,
    ts: &[Rational; 2],
    ends: &[[Rational; 2]; 2],
) -> Result<(Id, [[Ball; 2]; 2]), Refusal> {
    let dt = arithmetic(&ts[1], &ts[0], '-', w.dec)?;
    if dt.is_zero() {
        return Err(bad(
            "rounding-topology",
            "written edge has no parameter extent",
        ));
    }
    let mut derivative = Vec::new();
    let mut origin = Vec::new();
    for axis in 0..2 {
        let delta = arithmetic(&ends[1][axis], &ends[0][axis], '-', w.dec)?;
        let d = arithmetic(&delta, &dt, '/', w.dec)?;
        let displacement = arithmetic(&d, &ts[0], '*', w.dec)?;
        origin.push(arithmetic(&ends[0][axis], &displacement, '-', w.dec)?);
        derivative.push(d);
    }
    let squared = arithmetic(&derivative[0], &derivative[0], '*', w.dec)?
        + arithmetic(&derivative[1], &derivative[1], '*', w.dec)?;
    let length = Ball::from_rational(&squared, BITS)
        .sqrt()
        .ok_or_else(|| bad("rounding-topology", "pcurve direction has no norm"))?
        .mid();
    if length.is_zero() {
        return Err(bad("rounding-topology", "pcurve direction vanished"));
    }
    let scale = derivative[0].abs().max(derivative[1].abs());
    let mut point_fields = Vec::new();
    let mut o = Vec::new();
    let mut direction_fields = Vec::new();
    let mut v = Vec::new();
    for axis in 0..2 {
        let (field, value) = rounding.scalar(&origin[axis], w)?;
        point_fields.push(field);
        o.push(value);
        let ratio = arithmetic(&derivative[axis], &scale, '/', w.dec)?;
        let (field, value) = rounding.scalar(&ratio, w)?;
        direction_fields.push(field);
        v.push(value);
    }
    let (magnitude, m) = rounding.scalar(&length, w)?;
    let point = w.graph.add(
        "CARTESIAN_POINT",
        vec![text("''"), Field::List(point_fields)],
    );
    let direction = w
        .graph
        .add("DIRECTION", vec![text("''"), Field::List(direction_fields)]);
    let vector = w
        .graph
        .add("VECTOR", vec![text("''"), direction.into(), magnitude]);
    let line = w
        .graph
        .add("LINE", vec![text("''"), point.into(), vector.into()]);
    let context = w.graph.add_complex(vec![
        ("GEOMETRIC_REPRESENTATION_CONTEXT", vec![text("2")]),
        ("REPRESENTATION_CONTEXT", vec![text("''"), text("''")]),
    ]);
    let definition = w.graph.add(
        "DEFINITIONAL_REPRESENTATION",
        vec![text("''"), refs([line]), context.into()],
    );
    let id = w.graph.add(
        "PCURVE",
        vec![text("''"), surface.into(), definition.into()],
    );
    let norm = Ball::from_rational(&(&v[0] * &v[0] + &v[1] * &v[1]), BITS)
        .sqrt()
        .ok_or_else(|| bad("rounding-topology", "written pcurve direction is invalid"))?;
    let mut actual = std::array::from_fn(|_| {
        std::array::from_fn(|_| Ball::from_rational(&Rational::zero(), BITS))
    });
    for endpoint in 0..2 {
        for axis in 0..2 {
            let slope = Ball::from_rational(&v[axis], BITS)
                .checked_div(&norm)
                .ok_or_else(|| bad("rounding-topology", "written pcurve direction vanished"))?
                * Ball::from_rational(&m, BITS);
            actual[endpoint][axis] = Ball::from_rational(&o[axis], BITS)
                + slope * Ball::from_rational(&ts[endpoint], BITS);
        }
    }
    Ok((id, actual))
}

fn pcurve_checks(
    edge: &n::NumericalEdge,
    face: &n::NumericalFace,
    actual: &[[Ball; 2]; 2],
    dec: &mut Decider,
) -> Result<NumericalCheck, Refusal> {
    let mut samples = Vec::new();
    for i in 0..=8 {
        dec.charge_work(1)?;
        let share = Rational::from_ratio(i, 8).unwrap();
        let t = &edge.parameters[0] + &(&share * &(&edge.parameters[1] - &edge.parameters[0]));
        let s = Ball::from_rational(&share, BITS);
        let uv: [Ball; 2] = std::array::from_fn(|axis| {
            &actual[0][axis] + &(&s * &(&actual[1][axis] - &actual[0][axis]))
        });
        let p = edge
            .curve
            .evaluate(&Ball::from_rational(&t, BITS), dec)?
            .ok_or_else(|| bad("unsupported", "written analytic edge has no evaluator"))?;
        let q = face.surface.evaluate(&uv[0], &uv[1], dec)?;
        let delta: [Ball; 3] = std::array::from_fn(|axis| &p[axis] - &q[axis]);
        let distance = (delta[0].square() + delta[1].square() + delta[2].square())
            .sqrt()
            .ok_or_else(|| bad("rounding-topology", "written pcurve distance is unresolved"))?
            .hi();
        samples.push(NumericalSample {
            parameters: vec![t, uv[0].mid(), uv[1].mid()],
            distance_upper_mm: distance,
        });
    }
    Ok(NumericalCheck::Samples {
        name: "file.actual-decimal-pcurve-edge-surface-distance".into(),
        scope: format!("{} on {}", edge.key, face.key),
        samples,
    })
}

pub(super) fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    if options.product_name.trim().is_empty()
        || !(1..=17).contains(&options.significant_digits)
        || options.max_scalar_rounding_error_mm.sign() != Sign::Positive
    {
        return Err(bad(
            "options",
            "named product, 1..17 digits and positive scalar rounding limit required",
        ));
    }
    let count = body.counts();
    if count.shells != 1 || count.lumps != 1 {
        return Err(bad("unsupported", "one retained closed shell required"));
    }
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: Rational::zero(),
        direction_error: Rational::zero(),
    };
    let mut d = body.definition().clone();
    let mut vertices = Vec::new();
    for vertex in &mut d.vertices {
        let point = w.point(&vertex.position)?;
        vertices.push(w.graph.add(
            "VERTEX_POINT",
            vec![text(string(&vertex.key)), point.into()],
        ));
        vertex.position = round3(&vertex.position, options)?;
    }
    let surfaces = rounded_surfaces(&mut w, &mut d)?;
    let geometry = rounded_curves(&mut w, &mut d)?;
    closed_circle_parameters(&mut d, w.dec)?;
    for face in &mut d.faces {
        for trim in face.loops.iter_mut().flatten() {
            let edge = &d.edges[trim.edge];
            let ends = if trim.sense == Sense::Same {
                edge.vertices
            } else {
                [edge.vertices[1], edge.vertices[0]]
            };
            for (i, index) in ends.into_iter().enumerate() {
                trim.uv[i] = uv(
                    &face.surface,
                    &d.vertices[index].position,
                    &trim.uv[i],
                    w.dec,
                )?;
            }
        }
    }
    // Validate the actual rounded chart before classifying a repeated edge as
    // a seam. Admission requires a complete periodic cell for such an edge;
    // use count alone is not a geometric seam condition.
    let admitted = n::build(d, body.options().clone(), w.dec)?;
    let d = admitted.definition();
    let mut parameter = ParameterRounding {
        error: Rational::zero(),
    };
    let mut edges = Vec::new();
    let mut parameter_checks = Vec::new();
    let uses = curved::pcurve_uses(&d.faces, d.edges.len(), w.dec)?;
    for (index, edge) in d.edges.iter().enumerate() {
        let mut pcurves = Vec::new();
        for (f, trim) in &uses[index] {
            let face = &d.faces[*f];
            w.dec.charge_work(1)?;
            let ends = if trim.sense == Sense::Same {
                trim.uv.clone()
            } else {
                [trim.uv[1].clone(), trim.uv[0].clone()]
            };
            let (curve, actual) = pcurve(
                &mut w,
                &mut parameter,
                surfaces[*f],
                &edge.parameters,
                &ends,
            )?;
            parameter_checks.push(pcurve_checks(edge, face, &actual, w.dec)?);
            pcurves.push(curve);
        }
        let curve = if pcurves.is_empty() {
            geometry[index]
        } else {
            w.graph.add(
                if uses[index].len() == 2 && uses[index][0].0 == uses[index][1].0 {
                    "SEAM_CURVE"
                } else {
                    "SURFACE_CURVE"
                },
                vec![
                    text(string(&edge.key)),
                    geometry[index].into(),
                    refs(pcurves),
                    text(".CURVE_3D."),
                ],
            )
        };
        edges.push(w.graph.add(
            "EDGE_CURVE",
            vec![
                text(string(&edge.key)),
                vertices[edge.vertices[0]].into(),
                vertices[edge.vertices[1]].into(),
                curve.into(),
                boolean(edge.parameters[0] < edge.parameters[1]),
            ],
        ));
    }
    // The canonical projected vertex walk checks finite trim topology. Actual
    // independently rounded pcurves get their own distance enclosures, including
    // endpoints; no rounded parameter gap is silently snapped or called exact.
    let mut checks = admitted.evidence().checks().to_vec();
    checks.extend(parameter_checks);
    let emitted = NumericalEvidence::new(body.options().tolerance_mm.clone(), checks, Vec::new())?;
    let mut faces = Vec::new();
    for (index, face) in d.faces.iter().enumerate() {
        let mut bounds = Vec::new();
        for (i, boundary) in face.loops.iter().enumerate() {
            let mut uses = Vec::new();
            for trim in boundary {
                w.dec.charge_work(1)?;
                uses.push(w.graph.add(
                    "ORIENTED_EDGE",
                    vec![
                        text("''"),
                        text("*"),
                        text("*"),
                        edges[trim.edge].into(),
                        boolean(trim.sense == Sense::Same),
                    ],
                ));
            }
            let loop_ = w.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]);
            bounds.push(w.graph.add(
                if i == 0 {
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
    let description = format!(
        "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};complete_body_bound=unestablished;max_scalar_rounding_error_mm={};max_direction_ratio_rounding_error={};max_pcurve_parameter_scalar_rounding_error={}",
        body.definition().name,
        body.options().tolerance_mm,
        body.evidence().checks().len(),
        body.evidence().repairs().len(),
        w.error,
        w.direction_error,
        parameter.error
    );
    let mut file = w.graph.write(&[solid], |roots| {
        crate::export::approximate_header(&options.product_name, roots, &description)
    });
    file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    let artifact = NumericalStep {
        bytes: file.into_bytes(),
        report: NumericalStepReport {
            original_definition: body.source_definition().clone(),
            source_history: body.history().to_vec(),
            source_name: body.definition().name.clone(),
            source_definition: body.definition().clone(),
            source_evidence: body.evidence().clone(),
            emitted_evidence: emitted,
            scalar_rounding_limit_mm: options.max_scalar_rounding_error_mm.clone(),
            max_scalar_rounding_error_mm: w.error,
            max_direction_ratio_rounding_error: w.direction_error,
            counts: [count.faces, count.edges, count.vertices],
            spline_rounding: None,
            checks: vec![
                "source.retained-shared-analytic-boundaries".into(),
                "file.decimal-analytic-carriers-and-pcurves-numerical-admission".into(),
                "file.explicit-cylindrical-pcurves".into(),
                "file.canonical-byte-identity".into(),
                "file.approximate-category-no-global-bound".into(),
                format!(
                    "file.pcurve-parameter-scalar-rounding-max={}",
                    parameter.error
                ),
            ],
        },
    };
    artifact.check_written(artifact.bytes(), w.dec)?;
    Ok(artifact)
}
