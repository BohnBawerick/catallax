//! Complete finite meetings through a section's retained cutting surface.

use super::*;
use catallax_numbers::{RealPoly, ZeroProof};
use catallax_shapes::{
    BSplineCurve, Plane, SectionHits, SectionSurface, Sphere, SplineSupport, Surface, SurfaceId,
};

mod finite_line;

pub(super) fn finite_line_cover(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<FiniteLineCover>, Refusal> {
    finite_line::cover(model, a, b, dec)
}
mod spline;

pub(super) fn native_line_meetings(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<Meetings>, Refusal> {
    finite_line::meetings(model, a, b, dec)
}

pub(crate) enum FiniteSurfacePoints {
    Contained,
    Points(Vec<Point3>),
}

/// Complete finite point sets for clipping, without two-sided contact labels.
pub(crate) fn finite_surface_points(
    model: &Model,
    curve: CurveId,
    surface: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<FiniteSurfacePoints>, Refusal> {
    if matches!(
        model.geom().curve(curve),
        Curve::Circle(_) | Curve::Ellipse(_)
    ) && let Surface::BSpline(support) = model.geom().surface(surface)
    {
        let equations: [RealPoly; 2] = model
            .geom()
            .curve(curve)
            .equations(model.geom())
            .expect("a conic has two world equations")
            .try_into()
            .expect("a conic has two world equations");
        let mut points: Vec<Point3> = Vec::new();
        let mut append = |patch: &catallax_shapes::RationalPatch| -> Result<(), Refusal> {
            'candidate: for at in patch.quadric_parameters(&equations, dec)? {
                let point = patch.point_at(&at, dec)?;
                for previous in &points {
                    if sign(
                        &point.minus(previous).norm_squared(),
                        dec,
                        "boolean.patch-section-meetings",
                    )? == Sign::Zero
                    {
                        continue 'candidate;
                    }
                }
                points.push(point);
            }
            Ok(())
        };
        match support {
            SplineSupport::RationalPatch(patch) => append(patch)?,
            SplineSupport::RationalAtlas(atlas) => {
                for chart in atlas.charts() {
                    append(chart.patch())?;
                }
            }
            _ => return Ok(None),
        }
        return Ok(Some(FiniteSurfacePoints::Points(points)));
    }
    if let (Curve::BSpline(spline), Surface::BSpline(support)) =
        (model.geom().curve(curve), model.geom().surface(surface))
    {
        let answer = match support {
            SplineSupport::RationalPatch(patch) => patch.spline_hits(spline, dec)?,
            SplineSupport::RationalAtlas(atlas) => atlas.spline_hits(spline, dec)?,
            _ => None,
        };
        return Ok(answer.map(|answer| match answer {
            SectionHits::Contained => FiniteSurfacePoints::Contained,
            SectionHits::Points(points) => FiniteSurfacePoints::Points(points),
        }));
    }
    if let (Curve::PatchPair(branch), Surface::Plane(plane)) =
        (model.geom().curve(curve), model.geom().surface(surface))
    {
        return Ok(Some(match branch.plane_hits(plane, dec)? {
            SectionHits::Contained => FiniteSurfacePoints::Contained,
            SectionHits::Points(points) => FiniteSurfacePoints::Points(points),
        }));
    }
    let (Curve::BSpline(spline), Surface::Plane(plane)) =
        (model.geom().curve(curve), model.geom().surface(surface))
    else {
        return Ok(None);
    };
    if sign(
        &plane.normal.norm_squared(),
        dec,
        "boolean.patch-section-meetings",
    )? == Sign::Zero
    {
        return Err(unresolved("a finite plane query needs a nonzero normal"));
    }
    // Retained extrusion sections have a complete source solver before span substitution.
    if catallax_intersect::is_bezier_extrusion_section(model.geom(), curve, dec)? {
        return Ok(None);
    }
    let equation = model
        .geom()
        .implicit(surface)
        .ok_or_else(|| unresolved("the finite plane has no polynomial equation"))?;
    Ok(Some(match spline::hits(spline, equation, dec)? {
        None => FiniteSurfacePoints::Contained,
        Some(points) => FiniteSurfacePoints::Points(points),
    }))
}

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
    if matches!(model.geom().curve(curve), Curve::Line(_))
        && let Some(equations) = model.geom().curve_equations(curve)
        && equations
            .iter()
            .any(|equation| equation.rational().is_none())
    {
        // The complete line constraints were validated before this candidate
        // check. Admit their affine coefficients before composing them with
        // an independent section root; all equations must still vanish exactly.
        let coordinates = [position.x.clone(), position.y.clone(), position.z.clone()];
        for equation in equations {
            let admitted;
            let equation = if equation.rational().is_some() {
                equation
            } else {
                let terms = equation
                    .terms()
                    .into_iter()
                    .map(|(monomial, mut coefficient)| {
                        if coefficient.rational_ref().is_none()
                            && let Some(rational) = dec.rational_value(&coefficient)?
                        {
                            coefficient = rational.into();
                        }
                        Ok((monomial.exponents().to_vec(), coefficient))
                    })
                    .collect::<Result<Vec<_>, Refusal>>()?;
                admitted = RealPoly::from_terms(equation.nvars(), terms);
                &admitted
            };
            match dec.sign(&equation.eval_real(&coordinates)) {
                Outcome::Proven(_) => return Ok(false),
                Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {}
                Outcome::ProvenCoincident(_) => {
                    return Err(unresolved("section meetings need exact curve incidence"));
                }
                Outcome::NotExact(error) => return Err(Refusal::from(error)),
                Outcome::Refused(error) => return Err(error),
            }
        }
        return Ok(true);
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, Real, with_operation_budget};
    use catallax_shapes::{Line, Vector3};

    #[test]
    fn line_incidence_admits_conjugate_coefficients_before_an_independent_cubic() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let quadratic = RealPoly::from_terms(
                1,
                vec![
                    (vec![0], Real::from_i64(5)),
                    (vec![1], Real::from_i64(-43)),
                    (vec![2], Real::from_i64(86)),
                ],
            );
            let [a, b]: [Real; 2] = quadratic.real_roots(&mut dec)?.try_into().unwrap();
            let sum = &a + &b;
            let cubic = RealPoly::from_terms(
                1,
                vec![(vec![3], Real::one()), (vec![0], Real::from_i64(-2))],
            )
            .real_roots(&mut dec)?
            .remove(0);
            let mut model = Model::new();
            // a+b=1/2, so this is the complete vertical line x=2, y=0.
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: Point3::new(&sum * &Real::from_i64(4), Real::zero(), Real::zero()),
                dir: Vector3::new(
                    &sum - &Real::from(Rational::from_ratio(1, 2).unwrap()),
                    Real::zero(),
                    Real::one(),
                ),
            }));
            let on = Point3::new(Real::from_i64(2), Real::zero(), cubic.clone());
            assert!(exact_incidence(&mut model, curve, &on, &mut dec)?);
            let off = Point3::new(
                Real::from(Rational::from_ratio(2049, 1024).unwrap()),
                Real::zero(),
                cubic,
            );
            assert!(!exact_incidence(&mut model, curve, &off, &mut dec)?);
            Ok(())
        });
        result.unwrap();
        assert!(!work.exhausted);
    }

    #[test]
    fn ellipse_clipping_keeps_both_saddle_hits_and_excludes_the_finite_outside() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let half = Real::from(Rational::from_ratio(1, 2).unwrap());
            let mut model = Model::new();
            let definition = catallax_shapes::Ellipse::new(
                Point3::new(Real::zero(), Real::zero(), half.clone()),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                Real::from_i64(4),
                Real::one(),
                &mut dec,
            )?;
            let ellipse = model.geom_mut().add_curve(Curve::Ellipse(definition));
            for (low, expected) in [
                (Rational::zero(), 2),
                (Rational::from_ratio(3, 4).unwrap(), 1),
            ] {
                let knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
                let controls = [low, Rational::from_i64(2)]
                    .into_iter()
                    .flat_map(|x| {
                        [Rational::zero(), Rational::from_i64(2)].map(|y| {
                            Point3::new(x.clone().into(), y.clone().into(), (&x * &y).into())
                        })
                    })
                    .collect();
                let carrier = catallax_shapes::BSplineSurface::new(
                    [1, 1],
                    [knots.clone(), knots],
                    [2, 2],
                    controls,
                    vec![Real::one(); 4],
                    &mut dec,
                )?;
                let patch = catallax_shapes::RationalPatch::new(carrier, &mut dec)?;
                let surface = model
                    .geom_mut()
                    .add_surface(Surface::BSpline(SplineSupport::RationalPatch(patch)));
                let Some(FiniteSurfacePoints::Points(points)) =
                    finite_surface_points(&model, ellipse, surface, &mut dec)?
                else {
                    panic!("isolated finite ellipse/saddle points");
                };
                assert_eq!(points.len(), expected);
                let mut below_one = 0;
                for point in points {
                    for residual in [
                        &point.z - &half,
                        &point.x * &point.y - half.clone(),
                        point.x.square() / Real::from_i64(4) + point.y.square() - Real::one(),
                    ] {
                        assert_eq!(
                            sign(&residual, &mut dec, "test.ellipse-saddle")?,
                            Sign::Zero
                        );
                    }
                    let x = dec.enclose_budgeted_at(&point.x, 64)?;
                    if x.hi() < Rational::one() {
                        below_one += 1;
                    } else {
                        assert!(x.lo() > Rational::one());
                    }
                }
                assert_eq!(below_one, expected - 1);
            }
            Ok(())
        });
        result.unwrap();
        assert!(!work.exhausted);
    }
}
