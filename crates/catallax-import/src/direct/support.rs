use super::*;
use catallax_numbers::Sign;

pub(super) fn move_face(
    faces: &mut [Face], key: &str, translation: &Vec3, dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let mut surface = faces[selected].surface.clone();
    let frame = match &mut surface {
        Surface::Plane(p) | Surface::Cylinder(p, _) => p,
        _ => return Err(bad("direct.support-family", "move requires a planar or cylindrical support")),
    };
    frame.origin = a::offset(&frame.origin, translation, &Rational::one(), '+', dec)?;
    replace(faces, key, surface, dec)
}

pub(super) fn offset(
    faces: &mut [Face], key: &str, distance: &Rational,
    options: &n::NumericalOptions, dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let signed = if faces[selected].same_sense { distance.clone() } else { -distance };
    let mut surface = faces[selected].surface.clone();
    let construction = match &mut surface {
        Surface::Plane(p) => {
            let square = a::dot(&p.axis, &p.axis, dec)?;
            if square.sign() != Sign::Positive {
                return Err(bad("direct.normal", "the selected support has no normal"));
            }
            // Newton integer roots use a conservative reservation before the
            // existing exact helper; ordinary rational normals remain exact.
            let bits = square.numerator().bit_length().saturating_add(square.denominator().bit_length()) as u64;
            dec.charge_work(bits.saturating_mul(bits).saturating_mul(1 + u64::from(bits.max(1).ilog2())))?;
            let (scale, error) = if let Some(length) = crate::exact::rational_sqrt(&square) {
                (a::arithmetic(&signed, &length, '/', dec)?, Rational::zero())
            } else {
                let (lo, hi) = a::sqrt_enclosure(&square, dec)?;
                if lo.sign() != Sign::Positive {
                    return Err(bad("direct.normal-bound", "normal length could not be separated from zero"));
                }
                let mid = a::arithmetic(&a::arithmetic(&lo, &hi, '+', dec)?, &Rational::from_i64(2), '/', dec)?;
                let scale = a::arithmetic(&signed, &mid, '/', dec)?;
                let e0 = a::arithmetic(&a::arithmetic(&scale, &lo, '*', dec)?, &signed, '-', dec)?.abs();
                let e1 = a::arithmetic(&a::arithmetic(&scale, &hi, '*', dec)?, &signed, '-', dec)?.abs();
                (scale, a::maximum(e0, e1, dec)?)
            };
            if error > options.tolerance_mm {
                return Err(bad("direct.offset-bound", "parallel displacement enclosure exceeds the stated tolerance"));
            }
            p.origin = a::offset(&p.origin, &p.axis, &scale, '+', dec)?;
            if error.is_zero() { Construction::Rational }
            else { Construction::BoundedOffset { error_upper_mm: error } }
        }
        Surface::Cylinder(_, radius) => {
            *radius = a::arithmetic(radius, &signed, '+', dec)?;
            if radius.sign() != Sign::Positive {
                return Err(bad("direct.radius", "offset would collapse or reverse the cylindrical support"));
            }
            Construction::Rational
        }
        _ => return Err(bad("direct.support-family", "offset requires a planar or cylindrical support")),
    };
    replace(faces, key, surface, dec)?;
    Ok(construction)
}

pub(super) fn replace(
    faces: &mut [Face], key: &str, surface: Surface, dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let (vertices, circles) = match (&faces[selected].surface, &surface) {
        (Surface::Plane(old), Surface::Plane(new)) => plane_vertices(faces, selected, old, new, dec)?,
        (Surface::Cylinder(old, r0), Surface::Cylinder(new, r1)) => cylinder_vertices(faces, selected, old, r0, new, r1, dec)?,
        _ => return Err(bad("direct.support-family", "replacement preserves the planar or cylindrical support family")),
    };
    if let (Surface::Plane(old), Surface::Plane(new)) = (&faces[selected].surface, &surface) {
        if a::dot(&old.axis, &new.axis, dec)?.sign() == Sign::Negative {
            faces[selected].same_sense = !faces[selected].same_sense;
        }
    }
    update_vertices(faces, &vertices, &circles, dec)?;
    faces[selected].surface = surface;
    Ok(Construction::Rational)
}

fn selected_vertices(face: &Face) -> BTreeMap<String, Vertex> {
    let mut result = BTreeMap::new();
    for boundary in &face.loops {
        for edge in &boundary.edges {
            result.insert(edge.start.name.clone(), edge.start.clone());
            result.insert(edge.end.name.clone(), edge.end.clone());
        }
        if let Some(v) = &boundary.vertex { result.insert(v.name.clone(), v.clone()); }
    }
    result
}

fn project_along(point: &Vec3, direction: &Vec3, support: &Placement, dec: &mut Decider) -> Result<Vec3, ImportError> {
    let denominator = a::dot(direction, &support.axis, dec)?;
    if denominator.is_zero() {
        return Err(bad("direct.intersection", "the incident direction is parallel to the replacement plane"));
    }
    let delta = a::sub(&support.origin, point, dec)?;
    let scale = a::arithmetic(&a::dot(&delta, &support.axis, dec)?, &denominator, '/', dec)?;
    Ok(a::offset(point, direction, &scale, '+', dec)?)
}

fn plane_vertices(
    faces: &[Face], selected: usize, old: &Placement, new: &Placement, dec: &mut Decider,
) -> Result<(BTreeMap<String, Vertex>, BTreeMap<String, Curve>), ImportError> {
    let unchanged_normal = parallel(&old.axis, &new.axis, dec)?;
    if unchanged_normal && a::dot(&a::sub(&new.origin, &old.origin, dec)?, &new.axis, dec)?.is_zero() {
        return Err(bad("direct.no-change", "the replacement plane is the same geometric support"));
    }
    let incidence = vertex_owners(faces, dec)?;
    let mut vertices = selected_vertices(&faces[selected]);
    for (key, vertex) in &mut vertices {
        let adjacent: Vec<_> = incidence.get(key).into_iter().flatten().copied().filter(|i| *i != selected).collect();
        if adjacent.len() == 2 && adjacent.iter().all(|i| matches!(faces[*i].surface, Surface::Plane(_))) {
            vertex.point = triple([new, plane(&faces[adjacent[0]].surface)?, plane(&faces[adjacent[1]].surface)?], dec)?;
        } else if adjacent.len() == 1 && unchanged_normal {
            let Surface::Cylinder(cylinder, _) = &faces[adjacent[0]].surface else {
                return Err(bad("direct.neighbourhood", "the boundary lacks two independent adjacent planes"));
            };
            if !parallel(&cylinder.axis, &new.axis, dec)? {
                return Err(bad("direct.neighbourhood", "a cylindrical rim must stay perpendicular to its wall axis"));
            }
            vertex.point = project_along(&vertex.point, &cylinder.axis, new, dec)?;
        } else {
            return Err(bad("direct.neighbourhood", "the replacement has an ambiguous or unsupported vertex neighbourhood"));
        }
    }
    let mut circles = BTreeMap::new();
    for boundary in &faces[selected].loops {
        for edge in &boundary.edges {
            match &edge.curve {
                Curve::Circle(frame, radius) if unchanged_normal && parallel(&frame.axis, &new.axis, dec)? => {
                    let mut next = frame.clone();
                    next.origin = project_along(&frame.origin, &old.axis, new, dec)?;
                    circles.insert(edge.edge.clone(), Curve::Circle(next, radius.clone()));
                }
                Curve::Line(..) => {}
                _ => return Err(bad("direct.boundary-curve", "replacement needs straight edges or perpendicular circular rims")),
            }
        }
    }
    Ok((vertices, circles))
}

fn cylinder_vertices(
    faces: &[Face], selected: usize, old: &Placement, old_radius: &Rational,
    new: &Placement, new_radius: &Rational, dec: &mut Decider,
) -> Result<(BTreeMap<String, Vertex>, BTreeMap<String, Curve>), ImportError> {
    if new_radius.sign() != Sign::Positive || !parallel(&old.axis, &new.axis, dec)? {
        return Err(bad("direct.neighbourhood", "cylindrical replacement requires a positive radius and parallel axis"));
    }
    let square = a::dot(&old.axis, &old.axis, dec)?;
    if square.is_zero() { return Err(bad("direct.normal", "cylinder axis is zero")); }
    let delta = a::sub(&new.origin, &old.origin, dec)?;
    let axial = a::arithmetic(&a::dot(&delta, &old.axis, dec)?, &square, '/', dec)?;
    let shift = a::offset(&delta, &old.axis, &axial, '-', dec)?;
    if crate::exact::is_zero(&shift) && old_radius == new_radius {
        return Err(bad("direct.no-change", "the replacement cylinder is the same geometric support"));
    }
    let factor = a::arithmetic(new_radius, old_radius, '/', dec)?;
    let incidence = vertex_owners(faces, dec)?;
    let mut vertices = selected_vertices(&faces[selected]);
    for (key, vertex) in &mut vertices {
        for adjacent in incidence.get(key).into_iter().flatten().copied().filter(|i| *i != selected) {
            let cap = plane(&faces[adjacent].surface)?;
            if !parallel(&cap.axis, &old.axis, dec)? {
                return Err(bad("direct.neighbourhood", "cylinder edits require perpendicular planar end supports"));
            }
        }
        let relative = a::sub(&vertex.point, &old.origin, dec)?;
        let axial = a::arithmetic(&a::dot(&relative, &old.axis, dec)?, &square, '/', dec)?;
        let centre = a::offset(&old.origin, &old.axis, &axial, '+', dec)?;
        let radial = a::sub(&vertex.point, &centre, dec)?;
        let centre = a::offset(&centre, &shift, &Rational::one(), '+', dec)?;
        vertex.point = a::offset(&centre, &radial, &factor, '+', dec)?;
    }
    let mut circles = BTreeMap::new();
    for boundary in &faces[selected].loops {
        for edge in &boundary.edges {
            match &edge.curve {
                Curve::Circle(frame, radius) if parallel(&frame.axis, &old.axis, dec)? => {
                    let mut next = frame.clone();
                    next.origin = a::offset(&frame.origin, &shift, &Rational::one(), '+', dec)?;
                    circles.insert(edge.edge.clone(), Curve::Circle(next, a::arithmetic(radius, &factor, '*', dec)?));
                }
                Curve::Line(..) => {}
                _ => return Err(bad("direct.boundary-curve", "cylinder edits require straight seams and circular end rims")),
            }
        }
    }
    Ok((vertices, circles))
}
