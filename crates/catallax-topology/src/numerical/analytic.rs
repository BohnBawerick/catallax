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
            | Surface::Cone(_)
            | Surface::BSpline(SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_))
    )
}

pub(super) fn validate_surface(surface: &Surface, dec: &mut Decider) -> Result<(), Refusal> {
    if let Surface::Cone(cone) = surface {
        catallax_shapes::Cone::new(
            cone.apex.clone(),
            cone.axis.clone(),
            cone.slope.clone(),
            cone.ref_dir.clone(),
            dec,
        )?;
    }
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
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if matches!(face.surface, Surface::Cone(_)) {
        for uv in face.loops.iter().flatten().flat_map(|trim| &trim.uv) {
            dec.charge_work(1)?;
            if uv[1].sign() != Sign::Positive {
                return Err(bad(
                    "trim",
                    "a finite cone chart must remain strictly above its apex",
                ));
            }
        }
    }
    if matches!(face.surface, Surface::Cylinder(_) | Surface::Cone(_)) {
        // Ordinary contractible charts remain below one turn. A full band
        // must instead provide one shared seam with two translated UV uses.
        if !periodic_cell(face, edges, dec)? {
            let mut seen = BTreeSet::new();
            for trim in face.loops.iter().flatten() {
                dec.charge_work(1)?;
                if !seen.insert(trim.edge) {
                    return Err(bad(
                        "trim",
                        "a paired seam must establish one complete angular period",
                    ));
                }
            }
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
                    "a complete periodic band requires a paired seam",
                ));
            }
        }
        return super::check_loops(face, None, dec);
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
    // Chord polygons omit circular bulges. Complete rims additionally require
    // continuous projected-disk containment and separation under the stated cap.
    let mut complete = false;
    for boundary in &face.loops {
        for trim in boundary {
            dec.charge_work(1)?;
            if let Some(edge) = edges.get(trim.edge)
                && edge.vertices[0] == edge.vertices[1]
                && matches!(edge.curve, Curve::Circle(_))
            {
                if boundary.len() != 1 {
                    return Err(bad(
                        "trim",
                        "a complete planar circle must own its boundary",
                    ));
                }
                complete = true;
            }
        }
    }
    if complete && face.loops.len() > 1 {
        circle_clearance(face, edges, tolerance, dec)?;
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
                    let point = crate::numerical::curve_at(&edge.curve, &t, dec)?
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
    super::check_loops(&sampled, Some((face, edges)), dec)
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

/// A complete finite circular boundary, independently of its material role.
pub(super) fn complete_circle<'a>(
    boundary: &[NumericalTrim],
    edges: &'a [NumericalEdge],
) -> Option<&'a catallax_shapes::Circle> {
    let [trim] = boundary else { return None };
    let edge = edges.get(trim.edge)?;
    if edge.vertices[0] != edge.vertices[1] {
        return None;
    }
    match &edge.curve {
        Curve::Circle(circle) => Some(circle),
        _ => None,
    }
}

struct ProjectedCircle {
    center: [Ball; 2],
    enclosing_radius: Rational,
    interior_radius: Option<Rational>,
}

fn projected_circle(
    face: &NumericalFace,
    circle: &catallax_shapes::Circle,
    outer: bool,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<ProjectedCircle, Refusal> {
    let center = [
        dec.enclose_budgeted_at(&circle.center.x, BITS)?,
        dec.enclose_budgeted_at(&circle.center.y, BITS)?,
        dec.enclose_budgeted_at(&circle.center.z, BITS)?,
    ];
    let center = face
        .surface
        .parameters_of(&center, dec)?
        .ok_or_else(|| bad("trim", "circular boundary has no finite planar projection"))?;
    let radius = dec
        .enclose_budgeted_at(&circle.radius_sq, BITS)?
        .sqrt()
        .ok_or_else(|| bad("trim", "circular boundary radius is unresolved"))?;
    let interior_radius = if outer {
        let Surface::Plane(plane) = &face.surface else {
            return Err(bad("trim", "circular boundary needs a planar support"));
        };
        let normal = [
            dec.enclose_budgeted_at(&circle.normal.x, BITS)?,
            dec.enclose_budgeted_at(&circle.normal.y, BITS)?,
            dec.enclose_budgeted_at(&circle.normal.z, BITS)?,
        ];
        let plane_normal = [
            dec.enclose_budgeted_at(&plane.normal.x, BITS)?,
            dec.enclose_budgeted_at(&plane.normal.y, BITS)?,
            dec.enclose_budgeted_at(&plane.normal.z, BITS)?,
        ];
        let dot = &normal[0] * &plane_normal[0]
            + &normal[1] * &plane_normal[1]
            + &normal[2] * &plane_normal[2];
        let normal_sq = normal[0].square() + normal[1].square() + normal[2].square();
        let plane_sq =
            plane_normal[0].square() + plane_normal[1].square() + plane_normal[2].square();
        let cosine = dot
            .square()
            .checked_div(&(normal_sq * plane_sq))
            .and_then(|value| value.sqrt())
            .ok_or_else(|| bad("trim", "circular outer projection is unresolved"))?;
        // A tilted circle projects to an ellipse. Its inscribed disk uses
        // r*|cos(angle)|, not the enclosing radius used for holes.
        Some(&radius.lo() * &cosine.lo() - tolerance)
    } else {
        None
    };
    Ok(ProjectedCircle {
        center,
        enclosing_radius: &radius.hi() + tolerance,
        interior_radius,
    })
}

fn segment_distance_squared(
    center: &[Ball; 2],
    trim: &NumericalTrim,
    dec: &mut Decider,
) -> Result<Ball, Refusal> {
    dec.charge_work(1)?;
    let [a, b] = &trim.uv;
    let v = [&b[0] - &a[0], &b[1] - &a[1]];
    let length = &v[0] * &v[0] + &v[1] * &v[1];
    if length.is_zero() {
        return Err(bad("trim", "empty boundary segment"));
    }
    let delta: [Ball; 2] = std::array::from_fn(|i| &center[i] - &Ball::from_rational(&a[i], BITS));
    let parameter = (&delta[0] * &Ball::from_rational(&v[0], BITS)
        + &delta[1] * &Ball::from_rational(&v[1], BITS))
        .checked_div(&Ball::from_rational(&length, BITS))
        .ok_or_else(|| bad("trim", "boundary segment projection is unresolved"))?;
    let zero = Rational::zero();
    let one = Rational::one();
    let lo = Ball::from_rational(&parameter.lo().max(zero.clone()).min(one.clone()), BITS);
    let hi = Ball::from_rational(&parameter.hi().max(zero).min(one), BITS);
    let parameter = lo.hull(&hi);
    let residual: [Ball; 2] =
        std::array::from_fn(|i| &delta[i] - &(&parameter * &Ball::from_rational(&v[i], BITS)));
    Ok(residual[0].square() + residual[1].square())
}

fn circle_clearance(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut circles = Vec::with_capacity(face.loops.len());
    for (index, boundary) in face.loops.iter().enumerate() {
        circles.push(
            complete_circle(boundary, edges)
                .map(|circle| projected_circle(face, circle, index == 0, tolerance, dec))
                .transpose()?,
        );
    }
    for (i, circle) in circles.iter().enumerate() {
        let Some(circle) = circle else { continue };
        for (j, boundary) in face.loops.iter().enumerate() {
            if i == j {
                continue;
            }
            dec.charge_work(1)?;
            if let Some(other) = &circles[j] {
                if i > j {
                    continue;
                }
                let distance = ((&circle.center[0] - &other.center[0]).square()
                    + (&circle.center[1] - &other.center[1]).square())
                .sqrt()
                .ok_or_else(|| bad("trim", "circle separation is unresolved"))?;
                let clear = if let Some(radius) = &circle.interior_radius {
                    distance.hi() + &other.enclosing_radius < *radius
                } else {
                    distance.lo() > &circle.enclosing_radius + &other.enclosing_radius
                };
                if !clear {
                    return Err(bad("trim", "circular boundaries lack continuous clearance"));
                }
                continue;
            }
            for trim in boundary {
                if !edges
                    .get(trim.edge)
                    .is_some_and(|edge| matches!(edge.curve, Curve::Line(_)))
                {
                    return Err(bad(
                        "trim",
                        "complete circular boundaries need finite-line or complete-circle neighbours",
                    ));
                }
                if let Some(radius) = &circle.interior_radius {
                    dec.charge_work(1)?;
                    let delta: [Ball; 2] = std::array::from_fn(|axis| {
                        &circle.center[axis] - &Ball::from_rational(&trim.uv[0][axis], BITS)
                    });
                    if radius.sign() != Sign::Positive
                        || (delta[0].square() + delta[1].square()).hi() >= radius * radius
                    {
                        return Err(bad("trim", "hole leaves the circular outer boundary"));
                    }
                } else if segment_distance_squared(&circle.center, trim, dec)?.lo()
                    <= &circle.enclosing_radius * &circle.enclosing_radius
                {
                    return Err(bad(
                        "trim",
                        "complete circular hole lacks continuous clearance above the requested tolerance",
                    ));
                }
            }
        }
    }
    Ok(())
}
