use super::*;
use catallax_numbers::Sign;

struct BoundaryUpdate {
    vertices: BTreeMap<String, Vertex>,
    curves: BTreeMap<String, Curve>,
    intersection_error_upper_mm: Rational,
}

pub(super) fn move_face(
    faces: &mut [Face],
    key: &str,
    translation: &Vec3,
    dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let mut surface = faces[selected].surface.clone();
    let frame = match &mut surface {
        Surface::Plane(p) | Surface::Cylinder(p, _) => p,
        _ => {
            return Err(bad(
                "direct.support-family",
                "move requires a planar or cylindrical support",
            ));
        }
    };
    frame.origin = a::offset(&frame.origin, translation, &Rational::one(), '+', dec)?;
    replace(faces, key, &surface, dec)
}

pub(super) fn offset(
    faces: &mut [Face],
    key: &str,
    distance: &Rational,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let signed = if faces[selected].same_sense {
        distance.clone()
    } else {
        -distance
    };
    let mut surface = faces[selected].surface.clone();
    let offset_error_upper_mm = match &mut surface {
        Surface::Plane(p) => {
            let square = a::dot(&p.axis, &p.axis, dec)?;
            let [lo, hi] = length_bounds(&square, dec)?;
            let mid = a::arithmetic(
                &a::arithmetic(&lo, &hi, '+', dec)?,
                &Rational::from_i64(2),
                '/',
                dec,
            )?;
            let scale = a::arithmetic(&signed, &mid, '/', dec)?;
            let e0 =
                a::arithmetic(&a::arithmetic(&scale, &lo, '*', dec)?, &signed, '-', dec)?.abs();
            let e1 =
                a::arithmetic(&a::arithmetic(&scale, &hi, '*', dec)?, &signed, '-', dec)?.abs();
            p.origin = a::offset(&p.origin, &p.axis, &scale, '+', dec)?;
            a::maximum(e0, e1, dec)?
        }
        Surface::Cylinder(_, radius) => {
            *radius = a::arithmetic(radius, &signed, '+', dec)?;
            if radius.sign() != Sign::Positive {
                return Err(bad(
                    "direct.radius",
                    "offset would collapse or reverse the cylindrical support",
                ));
            }
            Rational::zero()
        }
        _ => {
            return Err(bad(
                "direct.support-family",
                "offset requires a planar or cylindrical support",
            ));
        }
    };
    let endpoints = replace(faces, key, &surface, dec)?;
    let intersection_error_upper_mm = match endpoints {
        Construction::Rational => Rational::zero(),
        Construction::Numerical {
            intersection_error_upper_mm,
            ..
        } => intersection_error_upper_mm,
    };
    if offset_error_upper_mm > options.tolerance_mm
        || intersection_error_upper_mm > options.tolerance_mm
    {
        return Err(bad(
            "direct.offset-bound",
            "offset or endpoint construction exceeds the stated tolerance",
        ));
    }
    Ok(accuracy(offset_error_upper_mm, intersection_error_upper_mm))
}

pub(super) fn replace(
    faces: &mut [Face],
    key: &str,
    surface: &Surface,
    dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let update = match (&faces[selected].surface, surface) {
        (Surface::Plane(old), Surface::Plane(new)) => {
            plane_vertices(faces, selected, old, new, dec)?
        }
        (Surface::Cylinder(old, r0), Surface::Cylinder(new, r1)) => {
            cylinder_vertices(faces, selected, old, r0, new, r1, dec)?
        }
        _ => {
            return Err(bad(
                "direct.support-family",
                "replacement preserves the planar or cylindrical support family",
            ));
        }
    };
    if let (Surface::Plane(old), Surface::Plane(new)) = (&faces[selected].surface, surface)
        && a::dot(&old.axis, &new.axis, dec)?.sign() == Sign::Negative
    {
        faces[selected].same_sense = !faces[selected].same_sense;
    }
    update_vertices(faces, &update.vertices, &update.curves, dec)?;
    faces[selected].surface = surface.clone();
    Ok(accuracy(
        Rational::zero(),
        update.intersection_error_upper_mm,
    ))
}

fn selected_vertices(face: &Face) -> BTreeMap<String, Vertex> {
    let mut result = BTreeMap::new();
    for boundary in &face.loops {
        for edge in &boundary.edges {
            result.insert(edge.start.name.clone(), edge.start.clone());
            result.insert(edge.end.name.clone(), edge.end.clone());
        }
        if let Some(v) = &boundary.vertex {
            result.insert(v.name.clone(), v.clone());
        }
    }
    result
}

fn project_along(
    point: &Vec3,
    direction: &Vec3,
    support: &Placement,
    dec: &mut Decider,
) -> Result<Vec3, ImportError> {
    let denominator = a::dot(direction, &support.axis, dec)?;
    if denominator.is_zero() {
        return Err(bad(
            "direct.intersection",
            "the incident direction is parallel to the replacement plane",
        ));
    }
    let delta = a::sub(&support.origin, point, dec)?;
    let scale = a::arithmetic(&a::dot(&delta, &support.axis, dec)?, &denominator, '/', dec)?;
    Ok(a::offset(point, direction, &scale, '+', dec)?)
}

fn plane_vertices(
    faces: &[Face],
    selected: usize,
    old: &Placement,
    new: &Placement,
    dec: &mut Decider,
) -> Result<BoundaryUpdate, ImportError> {
    let unchanged_normal = parallel(&old.axis, &new.axis, dec)?;
    if unchanged_normal
        && a::dot(&a::sub(&new.origin, &old.origin, dec)?, &new.axis, dec)?.is_zero()
    {
        return Err(bad(
            "direct.no-change",
            "the replacement plane is the same geometric support",
        ));
    }
    let incidence = vertex_owners(faces, dec)?;
    let mut vertices = selected_vertices(&faces[selected]);
    let mut error = Rational::zero();
    for (key, vertex) in &mut vertices {
        let adjacent: Vec<_> = incidence
            .get(key)
            .into_iter()
            .flatten()
            .copied()
            .filter(|i| *i != selected)
            .collect();
        if adjacent.len() == 2
            && adjacent
                .iter()
                .all(|i| matches!(faces[*i].surface, Surface::Plane(_)))
        {
            vertex.point = triple(
                [
                    new,
                    plane(&faces[adjacent[0]].surface)?,
                    plane(&faces[adjacent[1]].surface)?,
                ],
                dec,
            )?;
        } else if adjacent.len() == 1 && unchanged_normal {
            let Surface::Cylinder(cylinder, radius) = &faces[adjacent[0]].surface else {
                return Err(bad(
                    "direct.neighbourhood",
                    "the boundary lacks two independent adjacent planes",
                ));
            };
            if !parallel(&cylinder.axis, &new.axis, dec)? {
                return Err(bad(
                    "direct.neighbourhood",
                    "a cylindrical rim must stay perpendicular to its wall axis",
                ));
            }
            vertex.point = project_along(&vertex.point, &cylinder.axis, new, dec)?;
            error = a::maximum(
                error,
                circle_rim_error(&vertex.point, cylinder, radius, new, dec)?,
                dec,
            )?;
        } else {
            return Err(bad(
                "direct.neighbourhood",
                "the replacement has an ambiguous or unsupported vertex neighbourhood",
            ));
        }
    }
    let mut circles = BTreeMap::new();
    for boundary in &faces[selected].loops {
        for edge in &boundary.edges {
            match &edge.curve {
                Curve::Circle(frame, radius)
                    if unchanged_normal && parallel(&frame.axis, &new.axis, dec)? =>
                {
                    let mut next = frame.clone();
                    next.origin = project_along(&frame.origin, &old.axis, new, dec)?;
                    circles.insert(edge.edge.clone(), Curve::Circle(next, radius.clone()));
                }
                Curve::Line(..) => {}
                _ => {
                    return Err(bad(
                        "direct.boundary-curve",
                        "replacement needs straight edges or perpendicular circular rims",
                    ));
                }
            }
        }
    }
    Ok(BoundaryUpdate {
        vertices,
        curves: circles,
        intersection_error_upper_mm: error,
    })
}

fn cylinder_vertices(
    faces: &[Face],
    selected: usize,
    old: &Placement,
    old_radius: &Rational,
    new: &Placement,
    new_radius: &Rational,
    dec: &mut Decider,
) -> Result<BoundaryUpdate, ImportError> {
    if new_radius.sign() != Sign::Positive || !parallel(&old.axis, &new.axis, dec)? {
        return Err(bad(
            "direct.neighbourhood",
            "cylindrical replacement requires a positive radius and parallel axis",
        ));
    }
    let square = a::dot(&old.axis, &old.axis, dec)?;
    if square.is_zero() {
        return Err(bad("direct.normal", "cylinder axis is zero"));
    }
    let delta = a::sub(&new.origin, &old.origin, dec)?;
    let axial = a::arithmetic(&a::dot(&delta, &old.axis, dec)?, &square, '/', dec)?;
    let shift = a::offset(&delta, &old.axis, &axial, '-', dec)?;
    if crate::exact::is_zero(&shift) && old_radius == new_radius {
        return Err(bad(
            "direct.no-change",
            "the replacement cylinder is the same geometric support",
        ));
    }
    let incidence = vertex_owners(faces, dec)?;
    let ownership = owners(faces, dec)?;
    let mut vertices = selected_vertices(&faces[selected]);
    let mut error = Rational::zero();
    for (key, vertex) in &mut vertices {
        let mut cap = None;
        for adjacent in incidence[key].iter().copied().filter(|i| *i != selected) {
            let plane = plane(&faces[adjacent].surface)?;
            if !parallel(&plane.axis, &old.axis, dec)? {
                return Err(bad(
                    "direct.neighbourhood",
                    "cylinder edits require perpendicular planar end supports",
                ));
            }
            cap.get_or_insert(plane);
        }
        let cap = cap.ok_or_else(|| {
            bad(
                "direct.neighbourhood",
                "a cylindrical boundary endpoint requires a planar end support",
            )
        })?;
        let relative = a::sub(&vertex.point, &old.origin, dec)?;
        let axial = a::arithmetic(&a::dot(&relative, &old.axis, dec)?, &square, '/', dec)?;
        let centre = a::offset(&old.origin, &old.axis, &axial, '+', dec)?;
        let radial = a::sub(&vertex.point, &centre, dec)?;
        let [lo, hi] = length_bounds(&a::dot(&radial, &radial, dec)?, dec)?;
        let length = a::arithmetic(
            &a::arithmetic(&lo, &hi, '+', dec)?,
            &Rational::from_i64(2),
            '/',
            dec,
        )?;
        let scale = a::arithmetic(new_radius, &length, '/', dec)?;
        let centre = a::offset(&centre, &shift, &Rational::one(), '+', dec)?;
        let point = a::offset(&centre, &radial, &scale, '+', dec)?;
        vertex.point = project_along(&point, &old.axis, cap, dec)?;
        for adjacent in incidence[key].iter().copied().filter(|i| *i != selected) {
            let cap = plane(&faces[adjacent].surface)?;
            error = a::maximum(
                error,
                circle_rim_error(&vertex.point, new, new_radius, cap, dec)?,
                dec,
            )?;
        }
    }
    let mut curves = BTreeMap::new();
    for boundary in &faces[selected].loops {
        for edge in &boundary.edges {
            match &edge.curve {
                Curve::Circle(frame, _) if parallel(&frame.axis, &old.axis, dec)? => {
                    let other = ownership[&edge.edge]
                        .iter()
                        .copied()
                        .find(|i| *i != selected)
                        .ok_or_else(|| {
                            bad(
                                "direct.neighbourhood",
                                "a circular rim requires its retained cap",
                            )
                        })?;
                    let cap = plane(&faces[other].surface)?;
                    let mut next = frame.clone();
                    next.origin = project_along(&new.origin, &new.axis, cap, dec)?;
                    curves.insert(edge.edge.clone(), Curve::Circle(next, new_radius.clone()));
                }
                Curve::Line(..) => {}
                _ => {
                    return Err(bad(
                        "direct.boundary-curve",
                        "cylinder edits require straight seams and circular end rims",
                    ));
                }
            }
        }
    }
    Ok(BoundaryUpdate {
        vertices,
        curves,
        intersection_error_upper_mm: error,
    })
}

fn length_bounds(square: &Rational, dec: &mut Decider) -> Result<[Rational; 2], ImportError> {
    if square.sign() != Sign::Positive {
        return Err(bad(
            "direct.normal",
            "a radial or normal direction has zero length",
        ));
    }
    let bits = square
        .numerator()
        .bit_length()
        .saturating_add(square.denominator().bit_length());
    dec.charge_work(
        bits.saturating_mul(bits)
            .saturating_mul(1 + u64::from(bits.max(1).ilog2())),
    )?;
    if let Some(length) = crate::exact::rational_sqrt(square) {
        return Ok([length.clone(), length]);
    }
    let (lo, hi) = a::sqrt_enclosure(square, dec)?;
    if lo.sign() != Sign::Positive {
        return Err(bad(
            "direct.normal-bound",
            "direction length could not be separated from zero",
        ));
    }
    Ok([lo, hi])
}

/// Radial and axial directions are perpendicular. Their L1 error bounds the
/// distance to the exact cylinder/cap intersection at this radial phase.
fn circle_rim_error(
    point: &Vec3,
    cylinder: &Placement,
    radius: &Rational,
    cap: &Placement,
    dec: &mut Decider,
) -> Result<Rational, ImportError> {
    let delta = a::sub(point, &cylinder.origin, dec)?;
    let square = a::dot(&cylinder.axis, &cylinder.axis, dec)?;
    let along = a::arithmetic(&a::dot(&delta, &cylinder.axis, dec)?, &square, '/', dec)?;
    let radial = a::offset(&delta, &cylinder.axis, &along, '-', dec)?;
    let radial_square = a::dot(&radial, &radial, dec)?;
    let radius_square = a::arithmetic(radius, radius, '*', dec)?;
    let mut error = if radial_square == radius_square {
        Rational::zero()
    } else {
        let (lo, hi) = a::sqrt_enclosure(&radial_square, dec)?;
        a::maximum(
            a::arithmetic(&lo, radius, '-', dec)?.abs(),
            a::arithmetic(&hi, radius, '-', dec)?.abs(),
            dec,
        )?
    };
    let separation = a::dot(&a::sub(point, &cap.origin, dec)?, &cap.axis, dec)?.abs();
    if !separation.is_zero() {
        let [lo, _] = length_bounds(&a::dot(&cap.axis, &cap.axis, dec)?, dec)?;
        error = a::arithmetic(
            &error,
            &a::arithmetic(&separation, &lo, '/', dec)?,
            '+',
            dec,
        )?;
    }
    Ok(error)
}

fn accuracy(
    offset_error_upper_mm: Rational,
    intersection_error_upper_mm: Rational,
) -> Construction {
    if offset_error_upper_mm.is_zero() && intersection_error_upper_mm.is_zero() {
        Construction::Rational
    } else {
        Construction::Numerical {
            offset_error_upper_mm,
            intersection_error_upper_mm,
        }
    }
}
