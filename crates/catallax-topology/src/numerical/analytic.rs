//! Finite analytic numerical admission. No exact incidence is inferred.

use super::*;

pub(super) fn validate_circle(
    circle: &catallax_shapes::Circle,
    parameters: &[Rational; 2],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    catallax_shapes::Circle::new(
        circle.center.clone(),
        circle.normal.clone(),
        circle.radius_sq.clone(),
        circle.ref_dir.clone(),
        dec,
    )?;
    let turn = Ball::pi(BITS).lo() * Rational::from_i64(2);
    if (&parameters[1] - &parameters[0]).abs() >= turn {
        return Err(bad(
            "unsupported",
            "a finite circular arc must span less than one turn",
        ));
    }
    Ok(())
}

pub(super) fn supports(surface: &Surface) -> bool {
    matches!(
        surface,
        Surface::Plane(_)
            | Surface::Cylinder(_)
            | Surface::BSpline(SplineSupport::RationalPatch(_))
    )
}

pub(super) fn validate_surface(surface: &Surface, dec: &mut Decider) -> Result<(), Refusal> {
    if let Surface::Cylinder(cylinder) = surface {
        catallax_shapes::Cylinder::from_radius_sq(
            cylinder.origin.clone(),
            cylinder.axis.clone(),
            cylinder.radius_sq.clone(),
            cylinder.ref_dir.clone(),
            dec,
        )?;
    }
    Ok(())
}

pub(super) fn projected(face: &NumericalFace, edges: &[NumericalEdge]) -> bool {
    curved::projected(face, edges)
        || (matches!(face.surface, Surface::Plane(_))
            && face.loops.iter().flatten().any(|trim| {
                edges
                    .get(trim.edge)
                    .is_some_and(|edge| matches!(edge.curve, Curve::Circle(_)))
            }))
}

pub(super) fn projected_edge(face: &NumericalFace, edge: &NumericalEdge) -> bool {
    matches!(face.surface, Surface::Plane(_))
        && matches!(edge.curve, Curve::Circle(_) | Curve::BSpline(_))
}

pub(super) fn check_loops(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if let Surface::Cylinder(_) = face.surface {
        // A finite chart may be translated across a seam, but must not cover a
        // full turn or overlap itself periodically. Noncontractible bands need
        // a separate seam construction and are deliberately not guessed here.
        let mut range: Option<[Rational; 2]> = None;
        for uv in face.loops.iter().flatten().flat_map(|trim| &trim.uv) {
            dec.charge_work(1)?;
            range = Some(match range {
                None => [uv[0].clone(), uv[0].clone()],
                Some([lo, hi]) => [lo.min(uv[0].clone()), hi.max(uv[0].clone())],
            });
        }
        if let Some([lo, hi]) = range
            && hi - lo >= Ball::pi(BITS).lo() * Rational::from_i64(2)
        {
            return Err(bad(
                "trim",
                "cylindrical trim covers a complete turn or more",
            ));
        }
        return super::check_loops(face, dec);
    }
    if !projected(face, edges)
        || !face.loops.iter().flatten().any(|trim| {
            edges
                .get(trim.edge)
                .is_some_and(|edge| matches!(edge.curve, Curve::Circle(_)))
        })
    {
        return curved::check_loops(face, edges, dec);
    }
    let mut sampled = face.clone();
    sampled.loops.clear();
    for boundary in &face.loops {
        let mut polygon = Vec::new();
        for trim in boundary {
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "trim names an unknown edge"))?;
            let [a, b] = if trim.sense == Sense::Same {
                edge.parameters.each_ref()
            } else {
                [&edge.parameters[1], &edge.parameters[0]]
            };
            let mut previous = trim.uv[0].clone();
            for i in 1..=SUBDIVISIONS {
                dec.charge_work(1)?;
                let uv = if i == SUBDIVISIONS {
                    trim.uv[1].clone()
                } else {
                    let s = Rational::from_ratio(i, SUBDIVISIONS).unwrap();
                    let t = a + &(&s * &(b - a));
                    let point = edge
                        .curve
                        .evaluate(&Ball::from_rational(&t, BITS), dec)?
                        .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
                    curved::parameters(&face.surface, &point, dec)?
                };
                polygon.push(NumericalTrim {
                    edge: trim.edge,
                    sense: trim.sense,
                    uv: [previous, uv.clone()],
                });
                previous = uv;
            }
        }
        sampled.loops.push(polygon);
    }
    super::check_loops(&sampled, dec)
}
