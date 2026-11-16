//! Reconstruct the actual decimal carriers before numerical admission.

use super::*;
use catallax_numbers::Ball;
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Line, Plane, RationalAtlas, RationalPatch, SplineSupport,
};

mod spline;

pub(super) struct LineValues {
    pub origin: [Rational; 3],
    pub direction: [Rational; 3],
    pub magnitude: Rational,
}

pub(super) struct PlaneValues {
    pub origin: [Rational; 3],
    pub normal: [Rational; 3],
    pub reference: [Rational; 3],
}

#[derive(Default)]
struct Errors {
    scalar: Rational,
    direction: Rational,
    spline: NumericalSplineRounding,
}

pub(super) struct Payload {
    pub body: NumericalBody,
    pub lines: Vec<Option<LineValues>>,
    pub planes: Vec<Option<PlaneValues>>,
    pub scalar_error: Rational,
    pub direction_error: Rational,
    pub spline_rounding: NumericalSplineRounding,
}

fn rounded(
    value: &Real,
    error: &mut Rational,
    options: &NumericalStepOptions,
    length: bool,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    dec.charge_work(1)?;
    let bits = 128.min(dec.policy().max_precision_bits);
    let bound = dec.enclose_budgeted_at(value, bits)?;
    let middle = bound.mid();
    let token = written(&middle, options)?;
    let actual = (&token - &bound.lo())
        .abs()
        .max((&bound.hi() - &token).abs());
    if length && actual > options.max_scalar_rounding_error_mm {
        return Err(bad(
            "rounding-limit",
            "written length exceeds its separate rounding limit",
        ));
    }
    *error = error.clone().max(actual);
    Ok(token)
}

fn point(
    source: &Point3,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    [&source.x, &source.y, &source.z]
        .map(|value| rounded(value, &mut errors.scalar, options, true, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unreachable!())
}

fn unit(vector: &Vector3) -> Vector3 {
    vector.scaled(&(&Real::one() / &vector.norm_squared().sqrt()))
}

fn direction(
    source: &Vector3,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<([Rational; 3], Vector3), Refusal> {
    let source = unit(source);
    let values: [Rational; 3] = [&source.x, &source.y, &source.z]
        .map(|value| rounded(value, &mut errors.direction, options, false, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .expect("fixed coordinate count");
    if values.iter().all(Rational::is_zero) {
        return Err(bad("rounding-topology", "written direction vanished"));
    }
    let actual = Vector3::new(
        values[0].clone().into(),
        values[1].clone().into(),
        values[2].clone().into(),
    );
    Ok((values, unit(&actual)))
}

fn parameter(
    value: &Rational,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    rounded(
        &Real::from(value.clone()),
        &mut errors.spline.max_parameter_error,
        options,
        false,
        dec,
    )
}

pub(super) fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Payload, Refusal> {
    let mut definition = body.definition().clone();
    definition.name.push_str(" STEP decimal payload");
    let mut errors = Errors::default();
    for vertex in &mut definition.vertices {
        vertex.position = point(&point_from(&vertex.position), &mut errors, options, dec)?;
    }
    let mut lines = Vec::new();
    for edge in &mut definition.edges {
        let original_parameters = edge.parameters.clone();
        edge.parameters = edge
            .parameters
            .each_ref()
            .map(|value| parameter(value, &mut errors, options, dec))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .expect("fixed coordinate count");
        if edge.parameters[0] == edge.parameters[1] {
            return Err(bad("rounding-topology", "finite edge parameters collapsed"));
        }
        match &edge.curve {
            Curve::Line(source) => {
                let origin = point(&source.origin, &mut errors, options, dec)?;
                let (direction, actual) = direction(&source.dir, &mut errors, options, dec)?;
                let magnitude = rounded(
                    &source.dir.norm_squared().sqrt(),
                    &mut errors.scalar,
                    options,
                    true,
                    dec,
                )?;
                if magnitude.sign() != Sign::Positive {
                    return Err(bad("rounding-topology", "line magnitude vanished"));
                }
                edge.curve = Curve::Line(Line {
                    origin: point_from(&origin),
                    dir: actual.scaled(&Real::from(magnitude.clone())),
                });
                lines.push(Some(LineValues {
                    origin,
                    direction,
                    magnitude,
                }));
            }
            Curve::BSpline(source) => {
                edge.curve = Curve::BSpline(spline::curve(
                    source,
                    &original_parameters,
                    &mut errors,
                    options,
                    dec,
                )?);
                lines.push(None);
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical curved STEP needs retained line or native spline edges",
                ));
            }
        }
    }
    let bits = 128.min(dec.policy().max_precision_bits);
    let mut planes = Vec::new();
    for face in &mut definition.faces {
        match &face.surface {
            Surface::Plane(source) => {
                let origin = point(&source.origin, &mut errors, options, dec)?;
                let (normal, n) = direction(&source.normal, &mut errors, options, dec)?;
                let (reference, r) = direction(&source.ref_dir, &mut errors, options, dec)?;
                let x = unit(&n.cross(&r).cross(&n));
                face.surface = Surface::Plane(Plane {
                    origin: point_from(&origin),
                    normal: n,
                    ref_dir: x,
                });
                for trim in face.loops.iter_mut().flatten() {
                    let edge = &definition.edges[trim.edge];
                    let vertices = if trim.sense == Sense::Same {
                        edge.vertices
                    } else {
                        [edge.vertices[1], edge.vertices[0]]
                    };
                    for (slot, vertex) in vertices.into_iter().enumerate() {
                        let world = definition.vertices[vertex]
                            .position
                            .each_ref()
                            .map(|q| Ball::from_rational(q, bits));
                        let uv = face.surface.parameters_of(&world, dec)?.ok_or_else(|| {
                            bad("rounding-topology", "written plane inverse is unresolved")
                        })?;
                        trim.uv[slot] = uv.each_ref().map(Ball::mid);
                    }
                }
                planes.push(Some(PlaneValues {
                    origin,
                    normal,
                    reference,
                }));
            }
            Surface::BSpline(SplineSupport::RationalPatch(source)) => {
                let carrier = spline::surface(source.carrier(), &mut errors, options, dec)?;
                face.surface = Surface::BSpline(RationalPatch::new(carrier, dec)?.into());
                for trim in face.loops.iter_mut().flatten() {
                    for uv in &mut trim.uv {
                        for value in uv {
                            *value = parameter(value, &mut errors, options, dec)?;
                        }
                    }
                }
                planes.push(None);
            }
            Surface::BSpline(SplineSupport::RationalAtlas(source)) => {
                face.surface =
                    Surface::BSpline(spline::atlas(source, &mut errors, options, dec)?.into());
                for trim in face.loops.iter_mut().flatten() {
                    for uv in &mut trim.uv {
                        for value in uv {
                            *value = parameter(value, &mut errors, options, dec)?;
                        }
                    }
                }
                planes.push(None);
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical curved STEP needs planes or retained rational patches",
                ));
            }
        }
    }
    let mut native_trims = body.native_trims().to_vec();
    if !native_trims.is_empty() {
        let mut intervals = std::collections::BTreeMap::new();
        for edge in &body.definition().edges {
            dec.charge_work(1)?;
            intervals.insert(edge.key.as_str(), &edge.parameters);
        }
        for trim in &mut native_trims {
            let parameters = intervals
                .get(trim.edge.as_str())
                .ok_or_else(|| bad("rounding-topology", "native trim lost its stable edge"))?;
            trim.curve = spline::uv_curve(&trim.curve, parameters, &mut errors, options, dec)?;
        }
    }
    let emitted = if body.kind() == BodyKind::Surface {
        catallax_topology::numerical::build_surface_with_native_trims(
            definition, native_trims, body.options().clone(), dec,
        )?
    } else {
        catallax_topology::numerical::build_with_sampled_carriers(
            definition, native_trims, body.sampled_carriers().clone(), body.options().clone(), dec,
        )?
    };
    Ok(Payload {
        body: emitted,
        lines,
        planes,
        scalar_error: errors.scalar,
        direction_error: errors.direction,
        spline_rounding: errors.spline,
    })
}
