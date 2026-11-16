use super::*;

/// The four-sided SMLib construction: two side supports and two end supports.
/// Unlike a low-level face deletion, this reconstructs every shared edge use.
pub(super) fn face(
    faces: &mut Vec<Face>,
    key: &str,
    fillet: bool,
    dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let removed = &faces[selected];
    if removed.loops.len() != 1
        || !removed.loops[0].outer
        || removed.loops[0].vertex.is_some()
        || removed.loops[0].edges.len() != 4
    {
        return Err(bad(
            "direct.four-sided",
            "face removal requires one four-sided outer boundary",
        ));
    }
    let boundary = &removed.loops[0].edges;
    let edge_keys: BTreeSet<_> = boundary.iter().map(|e| e.edge.as_str()).collect();
    let corner_keys: BTreeSet<_> = boundary
        .iter()
        .flat_map(|e| [&e.start.name, &e.end.name])
        .collect();
    if edge_keys.len() != 4 || corner_keys.len() != 4 {
        return Err(bad(
            "direct.four-sided",
            "periodic seams or repeated corners do not define a removable four-sided feature",
        ));
    }
    let ownership = owners(faces, dec)?;
    let mut adjacent = Vec::with_capacity(4);
    for edge in boundary {
        let owners = &ownership[&edge.edge];
        if owners.len() != 2 {
            return Err(bad(
                "direct.neighbourhood",
                "each removed boundary edge needs one distinct neighbouring face",
            ));
        }
        adjacent.push(*owners.iter().find(|i| **i != selected).expect("two owners"));
    }
    if adjacent.iter().copied().collect::<BTreeSet<_>>().len() != 4 {
        return Err(bad(
            "direct.neighbourhood",
            "four distinct side and end faces are required",
        ));
    }
    let incidence = vertex_owners(faces, dec)?;
    if corner_keys.iter().any(|key| incidence[*key].len() != 3) {
        return Err(bad(
            "direct.neighbourhood",
            "high-valence corner removal is not uniquely defined",
        ));
    }
    let mut partition = None;
    for start in 0..2 {
        let ends = [start, start + 2];
        let sides = [(start + 1) % 4, (start + 3) % 4];
        let e0 = plane(&faces[adjacent[ends[0]]].surface)?;
        let e1 = plane(&faces[adjacent[ends[1]]].surface)?;
        let s0 = plane(&faces[adjacent[sides[0]]].surface)?;
        let s1 = plane(&faces[adjacent[sides[1]]].surface)?;
        if parallel(&e0.axis, &e1.axis, dec)? && !parallel(&s0.axis, &s1.axis, dec)? {
            if partition.is_some() {
                return Err(bad(
                    "direct.neighbourhood",
                    "the side/end support assignment is ambiguous",
                ));
            }
            partition = Some((ends, sides));
        }
    }
    let (ends, sides) = partition.ok_or_else(|| {
        bad(
            "direct.neighbourhood",
            "removal needs parallel end planes and intersecting side planes",
        )
    })?;
    let side_planes = [
        plane(&faces[adjacent[sides[0]]].surface)?,
        plane(&faces[adjacent[sides[1]]].surface)?,
    ];
    if fillet {
        let Surface::Cylinder(cylinder, radius) = &removed.surface else {
            return Err(bad(
                "direct.fillet-family",
                "fillet removal requires a cylindrical blend",
            ));
        };
        for side in side_planes {
            if !a::dot(&side.axis, &cylinder.axis, dec)?.is_zero() {
                return Err(bad(
                    "direct.fillet-tangency",
                    "a fillet side plane must be parallel to the cylinder axis",
                ));
            }
            let delta = a::sub(&cylinder.origin, &side.origin, dec)?;
            let distance = a::dot(&delta, &side.axis, dec)?;
            let lhs = a::arithmetic(&distance, &distance, '*', dec)?;
            let rr = a::arithmetic(radius, radius, '*', dec)?;
            let rhs = a::arithmetic(&rr, &a::dot(&side.axis, &side.axis, dec)?, '*', dec)?;
            if lhs != rhs {
                return Err(bad(
                    "direct.fillet-tangency",
                    "retained side supports do not certify cylindrical tangency",
                ));
            }
        }
        for end in ends {
            if !parallel(
                &plane(&faces[adjacent[end]].surface)?.axis,
                &cylinder.axis,
                dec,
            )? {
                return Err(bad(
                    "direct.fillet-end",
                    "fillet end planes must be perpendicular to its axis",
                ));
            }
        }
    }
    let mut vertices = BTreeMap::new();
    for end in ends {
        let edge = &boundary[end];
        let point = triple(
            [
                side_planes[0],
                side_planes[1],
                plane(&faces[adjacent[end]].surface)?,
            ],
            dec,
        )?;
        let vertex = Vertex {
            name: edge.start.name.clone(),
            point,
        };
        vertices.insert(edge.start.name.clone(), vertex.clone());
        vertices.insert(edge.end.name.clone(), vertex);
    }
    let rails = [&boundary[sides[0]], &boundary[sides[1]]];
    if rails.iter().any(|e| !matches!(e.curve, Curve::Line(..))) {
        return Err(bad(
            "direct.side-curve",
            "side rails must be straight support intersections",
        ));
    }
    let canonical = rails[0].edge.clone();
    let rail_keys = [rails[0].edge.clone(), rails[1].edge.clone()];
    let start = vertices[&rails[0].start.name].clone();
    let end = vertices[&rails[0].end.name].clone();
    let direction = a::sub(&end.point, &start.point, dec)?;
    if crate::exact::is_zero(&direction) {
        return Err(bad(
            "direct.topology-change",
            "extended side supports would collapse the replacement edge",
        ));
    }
    let collapsed = [
        boundary[ends[0]].edge.clone(),
        boundary[ends[1]].edge.clone(),
    ];
    faces.remove(selected);
    for face in faces.iter_mut() {
        for boundary in &mut face.loops {
            boundary.edges.retain(|e| !collapsed.contains(&e.edge));
            for edge in &mut boundary.edges {
                if rail_keys.contains(&edge.edge) {
                    let previous_start = &vertices[&edge.start.name];
                    let previous_end = &vertices[&edge.end.name];
                    if previous_start.name != start.name {
                        if previous_start.name != end.name || previous_end.name != start.name {
                            return Err(bad(
                                "direct.edge-correspondence",
                                "side rails do not share the reconstructed end vertices",
                            ));
                        }
                        edge.forward = !edge.forward;
                    } else if previous_end.name != end.name {
                        return Err(bad(
                            "direct.edge-correspondence",
                            "side rails do not connect opposite end supports",
                        ));
                    }
                    edge.edge = canonical.clone();
                    edge.start = start.clone();
                    edge.end = end.clone();
                    edge.curve = Curve::Line(start.point.clone(), direction.clone());
                    edge.curve_same_sense = true;
                }
            }
        }
    }
    update_vertices(faces, &vertices, &BTreeMap::new(), dec)?;
    Ok(Construction::Rational)
}

fn coaxial(a0: &Placement, b: &Placement, dec: &mut Decider) -> Result<bool, ImportError> {
    Ok(parallel(&a0.axis, &b.axis, dec)?
        && parallel(&a0.axis, &a::sub(&b.origin, &a0.origin, dec)?, dec)?)
}

/// Connected inward cylindrical walls, including internal blind/step caps.
/// Only complete inner loops on retained exterior faces may disappear.
pub(super) fn hole(
    faces: &mut Vec<Face>,
    key: &str,
    dec: &mut Decider,
) -> Result<Construction, ImportError> {
    let selected = index(faces, key)?;
    let Surface::Cylinder(axis, _) = &faces[selected].surface else {
        return Err(bad(
            "direct.hole-wall",
            "select an inward cylindrical hole wall",
        ));
    };
    if faces[selected].same_sense {
        return Err(bad(
            "direct.hole-wall",
            "an outside cylinder is not a cavity wall",
        ));
    }
    let ownership = owners(faces, dec)?;
    let mut walls = BTreeSet::new();
    for (i, face) in faces.iter().enumerate() {
        if !face.same_sense {
            if let Surface::Cylinder(candidate, _) = &face.surface {
                if coaxial(axis, candidate, dec)? {
                    walls.insert(i);
                }
            }
        }
    }
    let mut removable = walls.clone();
    for (i, face) in faces.iter().enumerate() {
        let Surface::Plane(cap) = &face.surface else {
            continue;
        };
        if !parallel(&axis.axis, &cap.axis, dec)? {
            continue;
        }
        let mut neighbours = BTreeSet::new();
        let mut circles_only = true;
        for boundary in &face.loops {
            if boundary.vertex.is_some() || boundary.edges.is_empty() {
                circles_only = false;
            }
            for edge in &boundary.edges {
                if let Curve::Circle(frame, _) = &edge.curve {
                    circles_only &= coaxial(axis, frame, dec)?;
                } else {
                    circles_only = false;
                }
                neighbours.extend(ownership[&edge.edge].iter().copied().filter(|j| *j != i));
            }
        }
        if circles_only && !neighbours.is_empty() && neighbours.iter().all(|j| walls.contains(j)) {
            removable.insert(i);
        }
    }
    let mut removed = BTreeSet::from([selected]);
    let mut pending = vec![selected];
    while let Some(i) = pending.pop() {
        for boundary in &faces[i].loops {
            for edge in &boundary.edges {
                dec.charge_work(1)?;
                for other in &ownership[&edge.edge] {
                    if removable.contains(other) && removed.insert(*other) {
                        pending.push(*other);
                    }
                }
            }
        }
    }
    let mut mouths = 0;
    for (i, face) in faces.iter_mut().enumerate() {
        if removed.contains(&i) {
            continue;
        }
        let mut keep = Vec::with_capacity(face.loops.len());
        for boundary in face.loops.drain(..) {
            let touches: Vec<_> = boundary
                .edges
                .iter()
                .map(|e| ownership[&e.edge].iter().any(|j| removed.contains(j)))
                .collect();
            if touches.iter().any(|v| *v) {
                if boundary.outer || boundary.vertex.is_some() || touches.iter().any(|v| !*v) {
                    return Err(bad(
                        "direct.hole-mouth",
                        "hole removal requires complete inner loops, not an outside wall or a partial slot",
                    ));
                }
                mouths += 1;
            } else {
                keep.push(boundary);
            }
        }
        face.loops = keep;
    }
    if mouths == 0 {
        return Err(bad(
            "direct.hole-mouth",
            "no complete exterior hole mouth was found",
        ));
    }
    let mut i = 0;
    faces.retain(|_| {
        let keep = !removed.contains(&i);
        i += 1;
        keep
    });
    Ok(Construction::Rational)
}
