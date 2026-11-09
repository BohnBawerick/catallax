//! Finite curve queries scoped to checked polynomial extrusion sections.

use crate::{Contact, CurveHit, CurveSurface};
use catallax_numbers::{
    Decider, FactKind, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof,
};
use catallax_shapes::{
    BSplineCurve, BezierExtrusion, Curve, CurveId, GeomKey, GeomStore, Meet, MeetOf, Plane, Point3,
    PointDef, Position, SplineSupport, Surface, SurfaceId, Vector3,
};

mod pair;

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(p) if p.proof != ZeroProof::SeparationBound => Ok(Sign::Zero),
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a polynomial curve query needs an exact incidence decision",
        )),
        Outcome::NotExact(error) => Err(Refusal::new(error.code, error.detail)),
        Outcome::Refused(error) => Err(error),
    }
}
fn equal(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn polynomial(curve: &BSplineCurve, dec: &mut Decider) -> Result<bool, Refusal> {
    let n = curve.control_count();
    if !(3..=4).contains(&n) || curve.degree() + 1 != n || curve.knots().len() != 2 * n {
        return Ok(false);
    }
    if curve.knots()[..n].iter().any(|k| !k.is_zero())
        || curve.knots()[n..].iter().any(|k| *k != Rational::one())
    {
        return Ok(false);
    }
    for weight in curve.weights() {
        if sign(&(weight - &Real::one()), dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}
fn matches_section(
    curve: &BSplineCurve,
    support: &BezierExtrusion,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(support.section_parameter(curve, dec)?.is_some())
}

// Incidence narrows the search. Matching the retained controls supplies the
// complete finite-span proof; incidence alone does not authorize the query.
fn section_data(
    geom: &GeomStore,
    id: CurveId,
    dec: &mut Decider,
) -> Result<Option<(BSplineCurve, Vector3)>, Refusal> {
    let Curve::BSpline(curve) = geom.curve(id) else {
        return Ok(None);
    };
    for fact in geom.facts().iter() {
        if fact.kind != FactKind::Incident {
            continue;
        }
        let surface = match (fact.a, fact.b) {
            (GeomKey::Curve(c), GeomKey::Surface(s)) | (GeomKey::Surface(s), GeomKey::Curve(c))
                if c == id =>
            {
                s
            }
            _ => continue,
        };
        if let Surface::BSpline(SplineSupport::BezierExtrusion(support)) = geom.surface(surface)
            && matches_section(curve, support, dec)?
        {
            return Ok(Some((curve.clone(), support.direction().clone())));
        }
    }
    Ok(None)
}

pub(super) fn section(
    geom: &GeomStore,
    id: CurveId,
    dec: &mut Decider,
) -> Result<Option<BSplineCurve>, Refusal> {
    Ok(section_data(geom, id, dec)?.map(|(curve, _)| curve))
}

pub(super) use pair::meetings;

pub(super) fn same(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let first = section(geom, a, dec)?;
    let second = section(geom, b, dec)?;
    let (curve, other, both) = match (first, second) {
        (Some(a), b_section) => (a, geom.curve(b), b_section.is_some()),
        (None, Some(b)) => (b, geom.curve(a), false),
        (None, None) => return Ok(None),
    };
    let Curve::BSpline(other) = other else {
        // A complete embedded finite open arc is neither an unbounded line nor
        // a closed circle or ellipse.
        return Ok(
            matches!(other, Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_)).then_some(false),
        );
    };
    if polynomial(other, dec)? && curve.control_count() == other.control_count() {
        for reverse in [false, true] {
            let mut same = true;
            for i in 0..curve.control_count() {
                let j = if reverse {
                    curve.control_count() - 1 - i
                } else {
                    i
                };
                if !equal(curve.control(i), other.control(j), dec)? {
                    same = false;
                    break;
                }
            }
            if same {
                return Ok(Some(true));
            }
        }
    }
    if both {
        let a = [curve.control(0), curve.control(curve.control_count() - 1)];
        let b = [other.control(0), other.control(other.control_count() - 1)];
        let forward = equal(a[0], b[0], dec)? && equal(a[1], b[1], dec)?;
        let reverse = equal(a[0], b[1], dec)? && equal(a[1], b[0], dec)?;
        if !forward && !reverse {
            return Ok(Some(false));
        }
    }
    Ok(None)
}

fn plane_hits(
    geom: &mut GeomStore,
    id: CurveId,
    surface: SurfaceId,
    curve: &BSplineCurve,
    plane: &Plane,
    dec: &mut Decider,
) -> Result<CurveSurface, Refusal> {
    let coordinates = curve.query_span_polynomials(curve.degree(), dec)?;
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let mut equation = RealPoly::zero(1);
    for i in 0..3 {
        equation = equation + (&coordinates[i] - &coordinates[3].scale(origin[i])).scale(normal[i]);
    }
    let mut zero = true;
    for (_, coefficient) in equation.terms() {
        if sign(&coefficient, dec)? != Sign::Zero {
            zero = false;
            break;
        }
    }
    if zero {
        return Ok(CurveSurface::Contained);
    }
    let mut hits = Vec::new();
    for parameter in equation.real_roots(dec)? {
        if sign(&parameter, dec)? == Sign::Negative
            || sign(&(&Real::one() - &parameter), dec)? == Sign::Negative
        {
            continue;
        }
        let mut derivative = equation.clone();
        let mut multiplicity = None;
        for order in 1..=curve.degree() {
            derivative = derivative.partial(0);
            if sign(&derivative.eval_real(std::slice::from_ref(&parameter)), dec)? != Sign::Zero {
                multiplicity = Some(order);
                break;
            }
        }
        let order = multiplicity.ok_or_else(|| {
            Refusal::new(
                "intersect.bezier-root-multiplicity",
                "the finite section root multiplicity was not proved",
            )
        })?;
        let point = geom.add_point(PointDef::Meet(Meet {
            of: MeetOf::CurveSurface(id, surface),
            position: Position::Closed(curve.point_at(&parameter, dec)?),
        }));
        hits.push(CurveHit {
            point,
            contact: if order % 2 == 0 {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        });
    }
    Ok(if hits.is_empty() {
        CurveSurface::Disjoint
    } else {
        CurveSurface::Meets(hits)
    })
}

pub(super) fn hits(
    geom: &mut GeomStore,
    id: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<CurveSurface>, Refusal> {
    let Some((curve, normal)) = section_data(geom, id, dec)? else {
        return Ok(None);
    };
    match geom.surface(surface).clone() {
        Surface::Plane(plane) => plane_hits(geom, id, surface, &curve, &plane, dec).map(Some),
        Surface::BSpline(SplineSupport::BezierExtrusion(support)) => {
            if matches_section(&curve, &support, dec)? {
                return Ok(Some(CurveSurface::Contained));
            }
            if sign(&normal.cross(support.direction()).norm_squared(), dec)? != Sign::Zero {
                return Err(Refusal::new(
                    "intersect.bezier-curve-surface",
                    "finite section-to-wall queries require parallel extrusion directions",
                ));
            }
            let direction = support.direction();
            let height = curve
                .control(0)
                .minus(&support.controls()[0])
                .dot(direction)
                / direction.norm_squared();
            if sign(&height, dec)? == Sign::Negative
                || sign(&(&Real::one() - &height), dec)? == Sign::Negative
            {
                return Ok(Some(CurveSurface::Disjoint));
            }
            let lower = support.lower_edge();
            let offset = direction.scaled(&height);
            let other = BSplineCurve::new(
                lower.degree(),
                lower.knots().to_vec(),
                support.controls().iter().map(|p| p.plus(&offset)).collect(),
                vec![Real::one(); lower.control_count()],
                dec,
            )?;
            let answer = pair::prove(&curve, &other, &normal, direction, dec)?;
            match answer {
                crate::BezierSectionPair::Coincident => Ok(Some(CurveSurface::Contained)),
                crate::BezierSectionPair::Points(points) => {
                    let coordinates = curve.span_polynomials(curve.degree(), dec)?;
                    let derivatives: Vec<_> =
                        coordinates[..3].iter().map(|p| p.partial(0)).collect();
                    let mut hits = Vec::new();
                    for position in points {
                        let parameter =
                            curve.unique_parameter_of(&position, dec)?.ok_or_else(|| {
                                Refusal::new(
                                    "intersect.bezier-curve-pair",
                                    "a proved section meeting has no parameter",
                                )
                            })?;
                        let tangent = Vector3::new(
                            derivatives[0].eval_real(std::slice::from_ref(parameter.value())),
                            derivatives[1].eval_real(std::slice::from_ref(parameter.value())),
                            derivatives[2].eval_real(std::slice::from_ref(parameter.value())),
                        );
                        let normal = support.normal_at(&position, dec)?;
                        let contact = if sign(&tangent.dot(&normal), dec)? == Sign::Zero {
                            Contact::Tangent
                        } else {
                            Contact::Crossing
                        };
                        let point = geom.add_point(PointDef::Meet(Meet {
                            of: MeetOf::CurveSurface(id, surface),
                            position: Position::Closed(position),
                        }));
                        hits.push(CurveHit { point, contact });
                    }
                    Ok(Some(if hits.is_empty() {
                        CurveSurface::Disjoint
                    } else {
                        CurveSurface::Meets(hits)
                    }))
                }
            }
        }
        _ => Err(Refusal::new(
            "intersect.bezier-curve-surface",
            "this surface has no checked finite Bezier section solver",
        )),
    }
}

pub(super) fn line_points(
    geom: &GeomStore,
    id: CurveId,
    line: CurveId,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let Some(curve) = section(geom, id, dec)? else {
        return Ok(None);
    };
    let crate::CurveForm::Line { origin, dir } = crate::curve_form(geom, line, dec)? else {
        return Ok(None);
    };
    let coordinates = curve.query_span_polynomials(curve.degree(), dec)?;
    let origins = [&origin.x, &origin.y, &origin.z];
    let delta: Vec<_> = (0..3)
        .map(|i| &coordinates[i] - &coordinates[3].scale(origins[i]))
        .collect();
    let equations = [
        delta[1].scale(&dir.z) - delta[2].scale(&dir.y),
        delta[2].scale(&dir.x) - delta[0].scale(&dir.z),
        delta[0].scale(&dir.y) - delta[1].scale(&dir.x),
    ];
    let mut selected = None;
    for equation in &equations {
        for (_, coefficient) in equation.terms() {
            if sign(&coefficient, dec)? != Sign::Zero {
                selected = Some(equation);
                break;
            }
        }
        if selected.is_some() {
            break;
        }
    }
    let equation = selected.ok_or_else(|| {
        Refusal::new(
            "intersect.bezier-line-overlap",
            "a collinear finite section needs a bounded-overlap result",
        )
    })?;
    let mut result = Vec::new();
    for parameter in equation.real_roots(dec)? {
        if sign(&parameter, dec)? == Sign::Negative
            || sign(&(&Real::one() - &parameter), dec)? == Sign::Negative
        {
            continue;
        }
        let mut on_line = true;
        for equation in &equations {
            if sign(&equation.eval_real(std::slice::from_ref(&parameter)), dec)? != Sign::Zero {
                on_line = false;
                break;
            }
        }
        if on_line {
            result.push(curve.point_at(&parameter, dec)?);
        }
    }
    Ok(Some(result))
}

pub(super) fn normal_acceleration(
    geom: &GeomStore,
    id: CurveId,
    position: &Point3,
    tangent: &catallax_shapes::Vector3,
    dec: &mut Decider,
) -> Result<Option<catallax_shapes::Vector3>, Refusal> {
    let Some(curve) = section(geom, id, dec)? else {
        return Ok(None);
    };
    let parameter = curve.unique_parameter_of(position, dec)?.ok_or_else(|| {
        Refusal::new(
            "intersect.bezier-point-off-curve",
            "normal acceleration needs a point on the finite section",
        )
    })?;
    let first = parameter.first(dec)?;
    let speed_sq = first.norm_squared();
    let tangent_sq = tangent.norm_squared();
    if sign(&speed_sq, dec)? != Sign::Positive
        || sign(&tangent_sq, dec)? != Sign::Positive
        || sign(&first.cross(tangent).norm_squared(), dec)? != Sign::Zero
    {
        return Err(Refusal::new(
            "intersect.bezier-tangent-mismatch",
            "normal acceleration needs a nonzero tangent to the section",
        ));
    }
    let coordinates = curve.span_polynomials(parameter.span(), dec)?;
    let second = catallax_shapes::Vector3::new(
        coordinates[0]
            .partial(0)
            .partial(0)
            .eval_real(std::slice::from_ref(parameter.value())),
        coordinates[1]
            .partial(0)
            .partial(0)
            .eval_real(std::slice::from_ref(parameter.value())),
        coordinates[2]
            .partial(0)
            .partial(0)
            .eval_real(std::slice::from_ref(parameter.value())),
    );
    let projection = first.dot(&second) / speed_sq.clone();
    Ok(Some(
        second
            .minus(&first.scaled(&projection))
            .scaled(&(tangent_sq / speed_sq)),
    ))
}
