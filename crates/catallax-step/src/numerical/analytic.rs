//! Approximate finite analytic exchange with explicit periodic pcurves.
//! Decimal carrier and pcurve values are rechecked as numerical geometry.

use super::*;
use catallax_numbers::{Ball, NumericalCheck, NumericalSample};
use catallax_shapes::{Circle, Cylinder, Line, Plane};
use catallax_topology::numerical::{self as n, NumericalDefinition};

const BITS: u32 = 128;

mod cones;

pub(super) fn applies(body: &NumericalBody) -> bool {
    let d = body.definition();
    body.native_trims().is_empty()
        && d.edges
            .iter()
            .all(|edge| matches!(edge.curve, Curve::Line(_) | Curve::Circle(_)))
        && d.faces.iter().all(|face| {
            matches!(
                face.surface,
                Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_)
            )
        })
        && (d
            .edges
            .iter()
            .any(|edge| matches!(edge.curve, Curve::Circle(_)))
            || d.faces
                .iter()
                .any(|face| matches!(face.surface, Surface::Cylinder(_) | Surface::Cone(_))))
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
        Surface::Cone(cone) => cones::uv(cone, p, old, dec),
        _ => Err(bad(
            "unsupported",
            "analytic writer requires planes, cylinders or cones",
        )),
    }
}

fn rounded_surfaces(
    w: &mut Writer<'_>,
    d: &mut NumericalDefinition,
    cones: &mut Vec<Option<cones::Parameters>>,
) -> Result<Vec<Id>, Refusal> {
    let mut result = Vec::new();
    for face in &mut d.faces {
        w.dec.charge_work(1)?;
        let mut cone_parameters = None;
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
            Surface::Cone(cone) => {
                let (id, surface, parameters) = cones::rounded(w, &face.key, cone)?;
                cone_parameters = Some(parameters);
                (id, surface)
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
        cones.push(cone_parameters);
    }
    Ok(result)
}

fn rounded_curves(
    w: &mut Writer<'_>,
    d: &mut NumericalDefinition,
) -> Result<(Vec<Id>, Vec<NumericalCheck>), Refusal> {
    let mut result = Vec::new();
    let mut phase_checks = Vec::new();
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
                // EDGE_CURVE omits its range. Anchor the natural lower end at
                // zero so neither a finite arc nor its pcurves cross the 3D
                // carrier's periodic origin when a reader projects vertices.
                // Only the reference phase changes: retain center, axis, radius,
                // endpoint identities, finite walk and every face's UV chart.
                let closed = edge.vertices[0] == edge.vertices[1];
                let anchor = usize::from(!closed && edge.parameters[1] < edge.parameters[0]);
                // Keep the existing bounded source-chart policy, even though
                // the emitted frame no longer depends on its absolute phase.
                if !closed {
                    let _ = near(Rational::zero(), &edge.parameters[anchor], w.dec)?;
                }
                let center = round3(&rational_point(&circle.center)?, w.options)?;
                let reference = vector_from(&difference(
                    &d.vertices[edge.vertices[anchor]].position,
                    &center,
                    w.dec,
                )?);
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
                phase_checks.push(circle_phase(
                    &circle,
                    &radius,
                    &d.vertices[edge.vertices[anchor]].position,
                    &edge.key,
                    closed,
                    w,
                )?);
                let mut ts = Vec::new();
                for (i, index) in edge.vertices.into_iter().enumerate() {
                    let [a, b, _] = local(
                        &circle.center,
                        &circle.normal,
                        &circle.ref_dir,
                        &d.vertices[index].position,
                    )?;
                    let reference = if closed {
                        edge.parameters[i].clone()
                    } else {
                        arithmetic(&edge.parameters[i], &edge.parameters[anchor], '-', w.dec)?
                    };
                    ts.push(angle(&a, &b, &reference, w.dec)?);
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
                if !closed {
                    ts[anchor] = Rational::zero();
                    let extent = &ts[1 - anchor];
                    if extent.sign() != Sign::Positive || extent >= &period {
                        return Err(bad(
                            "rounding-topology",
                            "written distinct endpoints do not select a positive arc below one turn",
                        ));
                    }
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
    Ok((result, phase_checks))
}

fn circle_phase(
    circle: &Circle,
    radius: &Rational,
    vertex: &[Rational; 3],
    key: &str,
    closed: bool,
    w: &mut Writer<'_>,
) -> Result<NumericalCheck, Refusal> {
    let kind = if closed { "closed-circle" } else { "arc" };
    let code = if closed {
        "step.numerical-closed-circle-phase"
    } else {
        "step.numerical-arc-phase"
    };
    w.dec.charge_work(1)?;
    let [x, y, _] = local(&circle.center, &circle.normal, &circle.ref_dir, vertex)?;
    if x.lo().sign() != Sign::Positive {
        return Err(Refusal::new(
            code,
            "written endpoint does not select the positive reference ray",
        ));
    }
    // For x > 0, |atan(y/x)| <= |y|/x. Check the actual written frame in
    // angular and physical units, independently of the looser source tolerance.
    // The angular cap follows OCCT Precision::Angular(); the physical selector
    // cap stays below the raw builder's 1e-7 mm precision. Neither is exactness.
    let tangent = y.lo().abs().max(y.hi().abs());
    let angle_upper = arithmetic(&tangent, &x.lo(), '/', w.dec)?;
    let displacement = arithmetic(radius, &angle_upper, '*', w.dec)?;
    let angular_limit = Rational::from_ratio(1, 1_000_000_000_000).unwrap();
    let distance_limit = w
        .options
        .max_scalar_rounding_error_mm
        .clone()
        .min(Rational::from_ratio(1, 1_000_000_000).unwrap());
    if angle_upper > angular_limit || displacement > distance_limit {
        return Err(Refusal::new(
            code,
            "decimal direction does not establish the forced zero phase within its angular and physical selector limits",
        ));
    }
    Ok(NumericalCheck::Samples {
        name: format!("file.{kind}-zero-phase-tangent-bound"),
        scope: format!("{key};angular_limit_rad=1e-12;distance_limit_mm={distance_limit}"),
        samples: vec![NumericalSample {
            parameters: vec![Rational::zero()],
            distance_upper_mm: displacement,
        }],
    })
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
    cone: Option<&cones::Parameters>,
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
        let q = if let Some(cone) = cone {
            cone.evaluate(&face.surface, &uv[0], &uv[1], dec)?
        } else {
            face.surface.evaluate(&uv[0], &uv[1], dec)?
        };
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

fn seam_loop_start(
    boundary: &[n::NumericalTrim],
    edges: &[n::NumericalEdge],
    face_sense: Sense,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    // OCCT TranslateEdgeLoop selects and installs both seam pcurves on each
    // occurrence. Its final choice depends on the cumulative edge/wire/face
    // orientation of the last occurrence. Give every retained seam a forward
    // final occurrence. Only rotate a closed walk; never reverse it, exchange
    // the seam's two charts, or alter geometry or parameter values.
    let forward = |trim: &n::NumericalTrim| {
        ((trim.sense == Sense::Same)
            == (edges[trim.edge].parameters[0] < edges[trim.edge].parameters[1]))
            == (face_sense == Sense::Same)
    };
    let mut last = std::collections::BTreeMap::new();
    for trim in boundary {
        dec.charge_work(1)?;
        let entry = last.entry(trim.edge).or_insert((0_usize, false));
        entry.0 += 1;
        entry.1 = forward(trim);
    }
    last.retain(|_, (count, _)| *count == 2);
    let mut wrong = last.values().filter(|(_, forward)| !forward).count();
    if wrong == 0 {
        return Ok(0);
    }
    // Moving the first use to the end changes only that edge's last use. This
    // linear walk tests every cyclic start without rebuilding the other uses.
    for (index, trim) in boundary.iter().enumerate() {
        dec.charge_work(1)?;
        if let Some((_, previous)) = last.get_mut(&trim.edge) {
            wrong -= usize::from(!*previous);
            *previous = forward(trim);
            wrong += usize::from(!*previous);
            if wrong == 0 {
                return Ok((index + 1) % boundary.len());
            }
        }
    }
    // A more complicated arrangement may have no such cyclic start. Keep its
    // admitted walk instead of inventing a different ordering or topology.
    Ok(0)
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
    let mut cone_parameters = Vec::with_capacity(d.faces.len());
    let surfaces = rounded_surfaces(&mut w, &mut d, &mut cone_parameters)?;
    let (geometry, phase_checks) = rounded_curves(&mut w, &mut d)?;
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
    let uses = curved::pcurve_uses(&d.faces, &[], d.edges.len(), w.dec)?;
    for (index, edge) in d.edges.iter().enumerate() {
        let mut pcurves = Vec::new();
        for (f, trim, _) in &uses[index] {
            let face = &d.faces[*f];
            w.dec.charge_work(1)?;
            let mut ends = if trim.sense == Sense::Same {
                trim.uv.clone()
            } else {
                [trim.uv[1].clone(), trim.uv[0].clone()]
            };
            if let (Surface::Cylinder(cylinder), Curve::Circle(circle)) =
                (&face.surface, &edge.curve)
                && edge.vertices[0] != edge.vertices[1]
                && ends[0][1] != ends[1][1]
            {
                // A circle normal parallel to the cylinder axis gives a
                // constant-height pcurve (OCCT ProjLib_Cylinder::Project).
                // This cheap check uses the written rational directions.
                // Keep existing flat maps, tilted carriers, the finite U chart
                // and the canonical shared-vertex walk used for admission.
                let normal = rational_vector(&circle.normal)?;
                let axis = rational_vector(&cylinder.axis)?;
                if charged_cross(&normal, &axis, w.dec)?
                    .iter()
                    .all(Rational::is_zero)
                {
                    let [_, _, height] = local(
                        &cylinder.origin,
                        &cylinder.axis,
                        &cylinder.ref_dir,
                        &rational_point(&circle.center)?,
                    )?;
                    let height = height.mid();
                    ends[0][1] = height.clone();
                    ends[1][1] = height;
                }
            }
            let (curve, actual) = pcurve(
                &mut w,
                &mut parameter,
                surfaces[*f],
                &edge.parameters,
                &ends,
            )?;
            parameter_checks.push(pcurve_checks(
                edge,
                face,
                &actual,
                cone_parameters[*f].as_ref(),
                w.dec,
            )?);
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
    checks.extend(phase_checks);
    checks.extend(parameter_checks);
    let emitted = NumericalEvidence::new(body.options().tolerance_mm.clone(), checks, Vec::new())?;
    let mut faces = Vec::new();
    for (index, face) in d.faces.iter().enumerate() {
        let mut bounds = Vec::new();
        for (i, boundary) in face.loops.iter().enumerate() {
            let mut uses = Vec::new();
            let start = if matches!(face.surface, Surface::Cylinder(_) | Surface::Cone(_)) {
                seam_loop_start(boundary, &d.edges, face.sense, w.dec)?
            } else {
                0
            };
            for trim in boundary.iter().cycle().skip(start).take(boundary.len()) {
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
    let mut description = format!(
        "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};complete_body_bound=unestablished;max_scalar_rounding_error_mm={};max_direction_ratio_rounding_error={};max_pcurve_parameter_scalar_rounding_error={}",
        body.definition().name,
        body.options().tolerance_mm,
        body.evidence().checks().len(),
        body.evidence().repairs().len(),
        w.error,
        w.direction_error,
        parameter.error
    );
    let has_cones = cone_parameters.iter().any(Option::is_some);
    let has_cylinders = d
        .faces
        .iter()
        .any(|face| matches!(face.surface, Surface::Cylinder(_)));
    let angle_error = cone_parameters
        .iter()
        .flatten()
        .map(|p| p.angle_rounding_error.clone())
        .max()
        .unwrap_or_else(Rational::zero);
    if has_cones {
        description.push_str(&format!(
            ";max_cone_semi_angle_rounding_error_rad={angle_error}"
        ));
    }
    let mut file_checks = vec![
        "source.retained-shared-analytic-boundaries".into(),
        "file.decimal-analytic-carriers-and-pcurves-numerical-admission".into(),
        if has_cones && !has_cylinders {
            "file.explicit-conical-pcurves"
        } else {
            "file.explicit-cylindrical-pcurves"
        }
        .into(),
        "file.canonical-byte-identity".into(),
        "file.approximate-category-no-global-bound".into(),
        format!(
            "file.pcurve-parameter-scalar-rounding-max={}",
            parameter.error
        ),
    ];
    if has_cones {
        if has_cylinders {
            file_checks.push("file.explicit-conical-pcurves".into());
        }
        file_checks.push("file.cone-step-v-axial-height".into());
        file_checks.push(format!(
            "file.cone-semi-angle-rounding-max-rad={angle_error}"
        ));
    }
    let mut file = w.graph.write(&[solid], |roots| {
        crate::export::approximate_header(&options.product_name, roots, &description)
    });
    file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    let artifact = NumericalStep {
        bytes: file.into_bytes(),
        report: NumericalStepReport {
            original_definition: body.source_definition().clone(),
            original_native_trims: body.source_native_trims().to_vec(),
            source_native_trims: body.native_trims().to_vec(),
            emitted_native_trims: Vec::new(),
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
            checks: file_checks,
        },
    };
    artifact.check_written(artifact.bytes(), w.dec)?;
    Ok(artifact)
}
