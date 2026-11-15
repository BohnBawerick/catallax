//! Finite analytic numerical admission. No exact incidence is inferred.

use super::*;

pub(super) fn validate_circle(
    circle: &catallax_shapes::Circle,
    parameters: &[Rational; 2],
    closed: bool,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    catallax_shapes::Circle::new(
        circle.center.clone(),
        circle.normal.clone(),
        circle.radius_sq.clone(),
        circle.ref_dir.clone(),
        dec,
    )?;
    let span = (&parameters[1] - &parameters[0]).abs();
    // A closed edge carries one explicit rational approximation of a full
    // period. Coincident positions with different vertex IDs do not close it.
    let valid = if closed {
        span == period()
    } else {
        span < Ball::pi(BITS).lo() * Rational::from_i64(2)
    };
    if !valid {
        return Err(bad(
            "unsupported",
            "a circle needs an open arc below one turn or an explicitly closed single period",
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
        // Ordinary contractible charts remain below one turn. A full band
        // must instead provide one shared seam with two translated UV uses.
        if !periodic_cell(face, edges, dec)? {
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
                    "a complete cylindrical band requires a paired seam",
                ));
            }
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

fn period() -> Rational {
    Ball::pi(BITS).mid() * Rational::from_i64(2)
}

/// Recognize a fundamental parameter cell, rather than admitting an arbitrary
/// wide polygon. The two seam branches may slope after decimal frame rounding,
/// but must be exact translations by one canonical period at both endpoints.
fn periodic_cell(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [boundary] = face.loops.as_slice() else {
        return Ok(false);
    };
    if boundary.len() != 4 {
        return Ok(false);
    }
    for offset in 0..2 {
        dec.charge_work(1)?;
        let [a, b, c, d] = std::array::from_fn::<_, 4, _>(|i| &boundary[(i + offset) % 4]);
        let Some(first) = edges.get(a.edge) else {
            return Ok(false);
        };
        let Some(second) = edges.get(c.edge) else {
            return Ok(false);
        };
        let Some(seam) = edges.get(b.edge) else {
            return Ok(false);
        };
        if a.edge == c.edge
            || b.edge != d.edge
            || b.sense == d.sense
            || !matches!(first.curve, Curve::Circle(_))
            || !matches!(second.curve, Curve::Circle(_))
            || first.vertices[0] != first.vertices[1]
            || second.vertices[0] != second.vertices[1]
            || !matches!(seam.curve, Curve::Line(_))
            || seam.vertices[0] == seam.vertices[1]
        {
            continue;
        }
        let p = period();
        if [&a.uv, &c.uv]
            .iter()
            .any(|uv| uv[0][1] != uv[1][1] || (&uv[1][0] - &uv[0][0]).abs() != p)
        {
            continue;
        }
        let translation = &b.uv[0][0] - &d.uv[1][0];
        if translation.abs() == p
            && &b.uv[1][0] - &d.uv[0][0] == translation
            && b.uv[0][1] == d.uv[1][1]
            && b.uv[1][1] == d.uv[0][1]
        {
            return Ok(true);
        }
    }
    Ok(false)
}
