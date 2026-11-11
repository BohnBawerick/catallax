//! Complete finite meetings through a section's retained cutting surface.

use super::*;
use catallax_numbers::ZeroProof;
use catallax_shapes::{BSplineCurve, Plane, SectionHits, SectionSurface, Sphere};

mod spline;

fn unresolved(detail: &str) -> Refusal {
    refusal("boolean.patch-section-meetings", detail)
}

fn exact_incidence(
    model: &mut Model,
    curve: CurveId,
    position: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    // Query an explicit point before giving it a meeting definition.
    let point = model
        .geom_mut()
        .add_point(PointDef::Explicit(position.clone()));
    match on_curve(model.geom(), curve, point, dec) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(true),
        Outcome::ProvenCoincident(_) => {
            Err(unresolved("section meetings need exact curve incidence"))
        }
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

fn store(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    position: Point3,
    points: &mut Vec<PointId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let point = model.geom_mut().add_point(PointDef::Meet(Meet {
        of: MeetOf::Curves(a, b),
        position: Position::Closed(position),
    }));
    canonical_point(model, points, point, dec, "boolean.imprint-ambiguous")?;
    Ok(())
}

fn constraints(
    model: &Model,
    curve: CurveId,
    dec: &mut Decider,
) -> Result<Vec<SectionSurface>, Refusal> {
    Ok(match model.geom().curve(curve) {
        Curve::PatchSection(branch) => vec![branch.surface().clone()],
        Curve::Line(line) => {
            let mut planes = Vec::new();
            for axis in 0..3 {
                let normal = line.dir.cross(&crate::geometry::unit(axis));
                match dec.sign(&normal.norm_squared()) {
                    Outcome::Proven(_) => planes.push(SectionSurface::Plane(Plane {
                        origin: line.origin.clone(),
                        normal,
                        ref_dir: line.dir.clone(),
                    })),
                    Outcome::ProvenCoincident(proof)
                        if proof.proof != ZeroProof::SeparationBound => {}
                    Outcome::ProvenCoincident(_) => {
                        return Err(unresolved(
                            "a line constraint needs an exact nonzero normal",
                        ));
                    }
                    Outcome::NotExact(error) => return Err(Refusal::from(error)),
                    Outcome::Refused(error) => return Err(error),
                }
            }
            planes
        }
        Curve::Circle(circle) => vec![
            SectionSurface::Plane(Plane {
                origin: circle.center.clone(),
                normal: circle.normal.clone(),
                ref_dir: circle.ref_dir.clone(),
            }),
            SectionSurface::Sphere(Sphere {
                center: circle.center.clone(),
                radius: circle.radius_sq.sqrt(),
                axis: circle.normal.clone(),
                ref_dir: circle.ref_dir.clone(),
            }),
        ],
        Curve::Ellipse(ellipse) => vec![SectionSurface::Plane(Plane {
            origin: ellipse.center.clone(),
            normal: ellipse.normal.clone(),
            ref_dir: ellipse.major_dir.clone(),
        })],
        Curve::Intersection(branch) => branch
            .surfaces
            .iter()
            .filter_map(|surface| SectionSurface::of(model.geom().surface(*surface)))
            .collect(),
        _ => Vec::new(),
    })
}

pub(super) fn meetings(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<Meetings>, Refusal> {
    let (a, b) = if matches!(model.geom().curve(a), Curve::PatchSection(_)) {
        (a, b)
    } else if matches!(model.geom().curve(b), Curve::PatchSection(_)) {
        (b, a)
    } else {
        return Ok(None);
    };
    let Curve::PatchSection(branch) = model.geom().curve(a).clone() else {
        unreachable!()
    };
    if a == b
        || model.geom().curve(a).same_definition(model.geom().curve(b))
        || model
            .geom()
            .facts()
            .holds(FactKind::Coincident, GeomKey::Curve(a), GeomKey::Curve(b))
    {
        return Ok(Some(Meetings::Coincident));
    }
    if let Curve::PatchSection(other) = model.geom().curve(b)
        && let Some(same) = branch.same_component(other, dec)?
    {
        if same {
            return Ok(Some(Meetings::Coincident));
        }
        // Distinct maximal components of this same section are disjoint.
        if branch
            .surface()
            .surface()
            .same_definition(&other.surface().surface())
        {
            return Ok(Some(Meetings::Points(Vec::new())));
        }
    }
    let mut points = Vec::new();
    let surface = model.geom_mut().add_surface(branch.surface().surface());
    let answer = if let Curve::BSpline(curve) = model.geom().curve(b).clone() {
        let equation = model
            .geom()
            .implicit(surface)
            .ok_or_else(|| unresolved("the section cutter has no polynomial equation"))?;
        spline::hits(&curve, equation, dec)?
    } else {
        let answer = intersect_curve_surface(model.geom_mut(), b, surface, dec)?;
        if !answer.marks.is_empty() {
            return Err(unresolved(
                "a marked meeting cannot establish an exact section split",
            ));
        }
        match answer.value {
            CurveSurface::Contained => None,
            CurveSurface::Disjoint => Some(Vec::new()),
            CurveSurface::Meets(hits) => Some(
                hits.into_iter()
                    .map(|hit| model.geom().position(hit.point))
                    .collect(),
            ),
        }
    };
    if let Some(candidates) = answer {
        // b meets a only among its complete meetings with a's cutter.
        // Finite carrier and component membership still belong to a.
        for position in candidates {
            if branch.contains(&position, dec)? {
                store(model, a, b, position, &mut points, dec)?;
            }
        }
        return Ok(Some(Meetings::Points(points)));
    }
    for constraint in constraints(model, b, dec)? {
        let SectionHits::Points(candidates) = branch.surface_hits(&constraint, dec)? else {
            continue;
        };
        // Every meeting is among this complete finite set. The remaining
        // carrier equations and b's branch selector must still hold exactly.
        for position in candidates {
            if exact_incidence(model, b, &position, dec)? {
                store(model, a, b, position, &mut points, dec)?;
            }
        }
        return Ok(Some(Meetings::Points(points)));
    }
    Err(unresolved(&format!(
        "a patch section and {:?} share the cutting equation; a finite common-surface curve solver or component-overlap proof is missing",
        model.geom().curve(b).kind(),
    )))
}
