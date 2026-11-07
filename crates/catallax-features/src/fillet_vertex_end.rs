//! A selected planar strip ending at a convex four-face vertex.

use super::*;

/// A closed planar skin with straight edges whose vertices lie in every
/// outward halfspace is convex. The halfspaces then certify entire trimmed
/// arcs, not just their endpoints, without a nearest-face sampling decision.
fn convex_polyhedron(
    model: &Model,
    body: BodyId,
    faces: &[FaceId],
    uses: &Uses,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if uses
        .keys()
        .any(|edge| !matches!(model.geom().curve(model.edge(*edge).curve), Curve::Line(_)))
    {
        return Err(bad(
            "a multiface strip end needs a straight-edged convex polyhedron",
        ));
    }
    let vertices = model.vertices_of_body(body);
    for &face in faces {
        let p = plane(model, face)?;
        let outside = if model.face(face).sense == Sense::Same {
            Sign::Positive
        } else {
            Sign::Negative
        };
        for vertex in &vertices {
            let at = model.geom().position(model.vertex(*vertex).point);
            if sign(&at.minus(&p.origin).dot(&p.normal), dec)? == outside {
                return Err(bad(
                    "a multiface strip end needs a convex input material sector",
                ));
            }
        }
    }
    Ok(())
}

/// All plane crossings on the complete conic are enumerated. No crossing may
/// lie strictly inside the retained arc; its one-sided material sense is then
/// fixed by one interior point. Tangent interior events also decline here.
fn contained_arc(
    model: &mut Model,
    faces: &[FaceId],
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let e = model.edge(edge).clone();
    let [Some(start), Some(end)] = [e.start, e.end] else {
        return Err(bad("a multiface end curve needs two distinct endpoints"));
    };
    let [start, end] = [start, end].map(|v| model.vertex(v).point);
    let [from, to] = if e.sense == Sense::Same {
        [start, end]
    } else {
        [end, start]
    };
    let sample =
        catallax_intersect::arc_point(model.geom(), e.curve, Some(from), Some(to), dec)?.value;
    for &face in faces {
        let p = plane(model, face)?;
        let outside = if model.face(face).sense == Sense::Same {
            Sign::Positive
        } else {
            Sign::Negative
        };
        if sign(&sample.minus(&p.origin).dot(&p.normal), dec)? == outside {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a multiface end arc leaves the input solid",
            ));
        }
        let surface = model.face(face).surface;
        let meeting =
            catallax_intersect::intersect_curve_surface(model.geom_mut(), e.curve, surface, dec)?
                .value;
        if let CurveSurface::Meets(hits) = meeting {
            for hit in hits {
                let at = model.geom().position(hit.point);
                if same(&at, &model.geom().position(start), dec)?
                    || same(&at, &model.geom().position(end), dec)?
                {
                    continue;
                }
                if on_trimmed(model, edge, hit.point, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "an interior end-arc event needs further neighbour continuation",
                    ));
                }
            }
        }
    }
    Ok(())
}

/// At a convex four-face vertex, the two unselected end walls clip the strip
/// in two conic arcs. Their junction is the unique hit of their shared input
/// edge in the strip's shorter normal sector. Both walls are transverse to
/// the cylinder axis, so projecting these arcs onto its normal circle is
/// one-to-one. Excluding the other rail selects each shorter sector interval;
/// together the arcs cover the complete strip end. No new edge run starts.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    old: VertexId,
    run: &mut Run,
    incident: &BTreeSet<FaceId>,
    uses: &Uses,
    edit: &mut Replacement,
    end_changes: &mut BTreeMap<(EdgeId, VertexId), VertexId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if incident.len() != 4 || !run.definition.convex {
        return Err(bad("a two-wall strip end needs a convex four-face vertex"));
    }
    let Surface::Cylinder(cylinder) = &run.definition.surface else {
        return Err(bad("a two-wall strip end needs a planar cylindrical strip"));
    };
    let axis = cylinder.axis.clone();
    let all_faces = model.faces_of_body(body);
    convex_polyhedron(model, body, &all_faces, uses, dec)?;
    let mut sides = [None; 2];
    let mut junction = None;
    for (edge, adjacent) in uses {
        let e = model.edge(*edge);
        if *edge == run.edge || (e.start != Some(old) && e.end != Some(old)) {
            continue;
        }
        let [first, second] = adjacent.as_slice() else {
            return Err(bad("a multiface end edge needs two incident faces"));
        };
        let own = run
            .faces
            .map(|face| adjacent.iter().any(|(f, _)| *f == face));
        match own {
            [true, false] | [false, true] => {
                let i = usize::from(own[1]);
                let wall = if first.0 == run.faces[i] {
                    second.0
                } else {
                    first.0
                };
                if sides[i].replace((*edge, wall)).is_some() {
                    return Err(bad("a strip support has more than one end-wall edge"));
                }
            }
            [false, false] => {
                if junction.replace((*edge, [first.0, second.0])).is_some() {
                    return Err(bad("a two-wall strip end has more than one junction edge"));
                }
            }
            _ => return Err(bad("another input edge pairs both strip supports")),
        }
    }
    let [Some((a_edge, a_wall)), Some((b_edge, b_wall))] = sides else {
        return Err(bad("each strip support must meet one end wall"));
    };
    let (junction, walls) =
        junction.ok_or_else(|| bad("the two end walls need one shared edge"))?;
    if a_wall == b_wall || !walls.contains(&a_wall) || !walls.contains(&b_wall) {
        return Err(bad("the end-wall fan is not one four-face cycle"));
    }
    let input = model.edge(run.edge);
    let far = if input.start == Some(old) {
        input.end
    } else {
        input.start
    }
    .ok_or_else(|| bad("a multiface strip end needs a finite input edge"))?;
    let at = model.geom().position(model.vertex(old).point);
    let along = model.geom().position(model.vertex(far).point).minus(&at);
    for wall in [a_wall, b_wall] {
        let p = plane(model, wall)?;
        if sign(&axis.dot(&p.normal), dec)? == Sign::Zero
            || sign(&outward(model, wall, &at, dec)?.dot(&along), dec)? != Sign::Negative
        {
            return Err(bad(
                "the end walls must be transverse and face the remaining run",
            ));
        }
    }
    let a = wall_foot(model, uses, run, run.faces[0], a_wall, old, dec)?;
    let b = wall_foot(model, uses, run, run.faces[1], b_wall, old, dec)?;
    let surface = model.geom_mut().add_surface(run.definition.surface.clone());
    let curve = model.edge(junction).curve;
    let meeting =
        catallax_intersect::intersect_curve_surface(model.geom_mut(), curve, surface, dec)?.value;
    let CurveSurface::Meets(hits) = meeting else {
        return Err(bad("the end-wall junction does not cross the strip"));
    };
    let mut hub = None;
    for hit in hits {
        let point = model.geom().position(hit.point);
        if on_trimmed(model, junction, hit.point, dec)? && run.definition.in_sector(&point, dec)? {
            if same(&point, &a, dec)? || same(&point, &b, dec)? || hub.replace(hit.point).is_some()
            {
                return Err(bad(
                    "the end-wall junction needs one strict strip-sector hit",
                ));
            }
        }
    }
    let hub = hub.ok_or_else(|| bad("the end-wall junction has no hit in the strip sector"))?;
    let hub = model.add_vertex(Vertex { point: hub });
    let feet = [vertex(model, a), vertex(model, b)];
    let mut caps = Vec::with_capacity(2);
    for (wall, from, to, other) in [
        (a_wall, feet[0], hub, run.faces[1]),
        (b_wall, hub, feet[1], run.faces[0]),
    ] {
        let avoid = line_plane(rail(run, other)?, plane(model, wall)?, dec)?;
        let support = model.geom().surface(model.face(wall).surface).clone();
        let cap = connector(
            model,
            &run.definition.surface,
            &support,
            from,
            to,
            &[&run.definition],
            Some(&avoid),
            dec,
        )?;
        contained_arc(model, &all_faces, cap, dec)?;
        edit.gaps.insert((wall, old), cap);
        caps.push(cap);
    }
    for (face, foot) in run.faces.into_iter().zip(feet) {
        run.ends.insert((old, face), foot);
    }
    run.caps.insert(old, caps);
    end_changes.insert((a_edge, old), feet[0]);
    end_changes.insert((b_edge, old), feet[1]);
    end_changes.insert((junction, old), hub);
    Ok(())
}
