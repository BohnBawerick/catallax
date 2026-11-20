//! Reconstruct the actual decimal carriers before numerical admission.

use super::*;
use catallax_numbers::{Ball, NumericalCheck};
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

// One value per edge or face of one body, built once and read once: boxing
// the common variant would only add an allocation each.
#[allow(clippy::large_enum_variant)]
pub(super) enum CurveValues {
    Line(LineValues),
    Native(Id),
    Spline,
}

#[allow(clippy::large_enum_variant)]
pub(super) enum SurfaceValues {
    Plane(PlaneValues),
    Native(Id),
    Spline,
}

#[derive(Default)]
struct Errors {
    scalar: Rational,
    direction: Rational,
    spline: NumericalSplineRounding,
}

pub(super) struct Payload {
    pub body: NumericalBody,
    pub curves: Vec<CurveValues>,
    pub surfaces: Vec<SurfaceValues>,
    pub scalar_error: Rational,
    pub direction_error: Rational,
    pub spline_rounding: NumericalSplineRounding,
    pub phase_checks: Vec<NumericalCheck>,
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

pub(super) fn build(body: &NumericalBody, w: &mut Writer<'_>) -> Result<Payload, Refusal> {
    let options = w.options;
    let mut definition = body.definition().clone();
    // STEP carries no copy tag and both copies are rounded below, so the
    // check of the written values measures every use and every retained
    // u/v curve.
    definition.reset_copy_tags();
    definition.name.push_str(" STEP decimal payload");
    let mut errors = Errors::default();
    for vertex in &mut definition.vertices {
        vertex.position = point(&point_from(&vertex.position), &mut errors, options, w.dec)?;
    }
    let mut curves = Vec::with_capacity(definition.edges.len());
    let mut phase_checks = Vec::new();
    for edge in &mut definition.edges {
        if let Curve::Circle(source) = &edge.curve {
            let written =
                super::super::analytic::rounded_circle(w, source, edge, &definition.vertices)?;
            edge.parameters = written.parameters;
            if edge.vertices[0] != edge.vertices[1] {
                for value in &mut edge.parameters {
                    *value = parameter(value, &mut errors, options, w.dec)?;
                }
            }
            edge.curve = Curve::Circle(written.circle);
            phase_checks.push(written.phase);
            curves.push(CurveValues::Native(written.entity));
            continue;
        }
        let original_parameters = edge.parameters.clone();
        for value in &mut edge.parameters {
            *value = parameter(value, &mut errors, options, w.dec)?;
        }
        if edge.parameters[0] == edge.parameters[1] {
            return Err(bad("rounding-topology", "finite edge parameters collapsed"));
        }
        match &edge.curve {
            Curve::Line(source) => {
                let origin = point(&source.origin, &mut errors, options, w.dec)?;
                let (direction, actual) = direction(&source.dir, &mut errors, options, w.dec)?;
                let magnitude = rounded(
                    &source.dir.norm_squared().sqrt(),
                    &mut errors.scalar,
                    options,
                    true,
                    w.dec,
                )?;
                if magnitude.sign() != Sign::Positive {
                    return Err(bad("rounding-topology", "line magnitude vanished"));
                }
                edge.curve = Curve::Line(Line {
                    origin: point_from(&origin),
                    dir: actual.scaled(&Real::from(magnitude.clone())),
                });
                curves.push(CurveValues::Line(LineValues {
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
                    w.dec,
                )?);
                curves.push(CurveValues::Spline);
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical curved STEP needs retained lines, circles or native spline edges",
                ));
            }
        }
    }
    super::super::analytic::closed_circle_parameters(&mut definition, w.dec)?;
    let bits = 128.min(w.dec.policy().max_precision_bits);
    let mut surfaces = Vec::with_capacity(definition.faces.len());
    for face in &mut definition.faces {
        match &face.surface {
            Surface::Plane(source) => {
                let origin = point(&source.origin, &mut errors, options, w.dec)?;
                let (normal, n) = direction(&source.normal, &mut errors, options, w.dec)?;
                let (reference, r) = direction(&source.ref_dir, &mut errors, options, w.dec)?;
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
                        let uv = face.surface.parameters_of(&world, w.dec)?.ok_or_else(|| {
                            bad("rounding-topology", "written plane inverse is unresolved")
                        })?;
                        trim.uv[slot] = uv.each_ref().map(Ball::mid);
                    }
                }
                surfaces.push(SurfaceValues::Plane(PlaneValues {
                    origin,
                    normal,
                    reference,
                }));
            }
            Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Torus(_) => {
                let (entity, surface) =
                    super::super::analytic::rounded_native_surface(w, &face.key, &face.surface)?;
                face.surface = surface;
                for trim in face.loops.iter_mut().flatten() {
                    let edge = &definition.edges[trim.edge];
                    let vertices = if trim.sense == Sense::Same {
                        edge.vertices
                    } else {
                        [edge.vertices[1], edge.vertices[0]]
                    };
                    for (slot, vertex) in vertices.into_iter().enumerate() {
                        trim.uv[slot] = super::super::analytic::uv(
                            &face.surface,
                            &definition.vertices[vertex].position,
                            &trim.uv[slot],
                            w.dec,
                        )?;
                    }
                }
                super::super::analytic::normalize_cylinder_cell(face, w.dec)?;
                surfaces.push(SurfaceValues::Native(entity));
            }
            Surface::BSpline(SplineSupport::RationalPatch(source)) => {
                let carrier = spline::surface(source.carrier(), &mut errors, options, w.dec)?;
                face.surface = Surface::BSpline(RationalPatch::new(carrier, w.dec)?.into());
                for trim in face.loops.iter_mut().flatten() {
                    for uv in &mut trim.uv {
                        for value in uv {
                            *value = parameter(value, &mut errors, options, w.dec)?;
                        }
                    }
                }
                surfaces.push(SurfaceValues::Spline);
            }
            Surface::BSpline(SplineSupport::RationalAtlas(source)) => {
                let written =
                    spline::atlas(source, &mut errors, options, w.dec).map_err(|error| {
                        Refusal::new(error.code, format!("face {:?}: {}", face.key, error.detail))
                    })?;
                face.surface = Surface::BSpline(written.into());
                for trim in face.loops.iter_mut().flatten() {
                    for uv in &mut trim.uv {
                        for value in uv {
                            *value = parameter(value, &mut errors, options, w.dec)?;
                        }
                    }
                }
                surfaces.push(SurfaceValues::Spline);
            }
            _ => {
                return Err(bad(
                    "unsupported",
                    "numerical curved STEP needs planes, cylinders, spheres, ring tori or retained rational patches",
                ));
            }
        }
    }
    let mut native_trims = body.native_trims().to_vec();
    for trim in &mut native_trims {
        trim.tag = catallax_topology::numerical::CopyTag::Independent;
    }
    if !native_trims.is_empty() {
        let mut intervals = std::collections::BTreeMap::new();
        for (source, edge) in body.definition().edges.iter().zip(&definition.edges) {
            w.dec.charge_work(1)?;
            intervals.insert(source.key.as_str(), (&source.parameters, edge));
        }
        for trim in &mut native_trims {
            let (parameters, edge) = intervals
                .get(trim.edge.as_str())
                .ok_or_else(|| bad("rounding-topology", "native trim lost its stable edge"))?;
            if matches!(edge.curve, Curve::Circle(_)) && edge.vertices[0] == edge.vertices[1] {
                return Err(bad(
                    "unsupported",
                    "a closed circle with a fitted pcurve needs a jointly representable full-period knot domain",
                ));
            }
            trim.curve = spline::uv_curve(
                &trim.curve,
                parameters,
                &edge.parameters,
                &mut errors,
                options,
                w.dec,
            )?;
        }
    }
    let build = if body.kind() == BodyKind::Surface {
        catallax_topology::numerical::build_surface_with_sampled_carriers
    } else {
        catallax_topology::numerical::build_with_sampled_carriers
    };
    let emitted = build(
        definition,
        native_trims,
        body.sampled_carriers().clone(),
        body.options().clone(),
        w.dec,
    )?;
    Ok(Payload {
        body: emitted,
        curves,
        surfaces,
        scalar_error: errors.scalar,
        direction_error: errors.direction,
        spline_rounding: errors.spline,
        phase_checks,
    })
}
