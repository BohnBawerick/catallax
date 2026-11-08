//! Planar strips at convex four-face ends and partial corners.

use super::*;
use catallax_shapes::{Plane, Sphere, SurfaceId};

/// Re-expresses one coordinate using a plane equality already proved by the
/// caller. This is the same point, not a projection or tolerance-based snap.
fn plane_coordinates(point: Point3, plane: &Plane, dec: &mut Decider) -> Result<Point3, Refusal> {
    let mut coordinates = [point.x, point.y, point.z];
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    for axis in 0..3 {
        if sign(normal[axis], dec)? == Sign::Zero {
            continue;
        }
        let mut offset = Real::zero();
        for other in 0..3 {
            if other != axis {
                offset = &offset + &(&(&coordinates[other] - origin[other]) * normal[other]);
            }
        }
        coordinates[axis] = origin[axis] - &(&offset / normal[axis]);
        let [x, y, z] = coordinates;
        return Ok(Point3::new(x, y, z));
    }
    Err(bad("a plane equality needs a nonzero normal"))
}

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

/// Complete-conic plane events fix the side of the whole retained arc.
/// All input and cone halfspaces share one interior-point certificate.
pub(super) fn contained_arc(
    model: &mut Model,
    faces: &[FaceId],
    extra: &[(SurfaceId, Sign)],
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
    for index in 0..faces.len() + extra.len() {
        let (surface, outside) = if index < faces.len() {
            let face = model.face(faces[index]);
            (
                face.surface,
                if face.sense == Sense::Same {
                    Sign::Positive
                } else {
                    Sign::Negative
                },
            )
        } else {
            extra[index - faces.len()]
        };
        let Surface::Plane(p) = model.geom().surface(surface) else {
            return Err(bad("a multiface arc certificate needs a plane"));
        };
        if sign(&sample.minus(&p.origin).dot(&p.normal), dec)? == outside {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a multiface arc leaves its prescribed halfspace",
            ));
        }
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
        if on_trimmed(model, junction, hit.point, dec)?
            && run.definition.in_sector(&point, dec)?
            && (same(&point, &a, dec)?
                || same(&point, &b, dec)?
                || hub.replace(hit.point).is_some())
        {
            return Err(bad(
                "the end-wall junction needs one strict strip-sector hit",
            ));
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
        contained_arc(model, &all_faces, &[], cap, dec)?;
        edit.gaps.insert((wall, old), Gap::Edge(cap));
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

pub(super) struct PlanarPair {
    pub(super) faces: Vec<FaceId>,
    pub(super) common: FaceId,
    pub(super) outer: [FaceId; 2],
    pub(super) wall: FaceId,
}

pub(super) fn planar_pair(
    model: &Model,
    body: BodyId,
    runs: [&Run; 2],
    incident: &BTreeSet<FaceId>,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<PlanarPair, Refusal> {
    let [first, second] = runs;
    if incident.len() != 4
        || !first.definition.convex
        || !second.definition.convex
        || sign(&(&first.definition.radius - &second.definition.radius), dec)? != Sign::Zero
    {
        return Err(bad(
            "a two-strip four-face junction needs equal convex planar strips",
        ));
    }
    let faces = model.faces_of_body(body);
    convex_polyhedron(model, body, &faces, uses, dec)?;
    let common = *first
        .faces
        .iter()
        .find(|face| second.faces.contains(face))
        .ok_or_else(|| bad("two selected strips must share one support"))?;
    let outer = [
        *first.faces.iter().find(|face| **face != common).unwrap(),
        *second.faces.iter().find(|face| **face != common).unwrap(),
    ];
    if outer[0] == outer[1] {
        return Err(bad("two selected strips need distinct outer supports"));
    }
    let wall = *incident
        .iter()
        .find(|face| **face != common && !outer.contains(face))
        .ok_or_else(|| bad("a two-strip four-face junction needs one end wall"))?;
    Ok(PlanarPair {
        faces,
        common,
        outer,
        wall,
    })
}

/// The selected A/B and B/C sheets supply the A/B/C sphere sector. Cliff
/// closure keeps its exposed A/C side planar and retains both lower strips;
/// smooth closure continues that side's full-radius contact sheet to the wall.
#[allow(clippy::too_many_arguments)]
pub(super) fn build_pair(
    model: &mut Model,
    body: BodyId,
    old: VertexId,
    runs: [&mut Run; 2],
    incident: &BTreeSet<FaceId>,
    uses: &Uses,
    edit: &mut Replacement,
    end_changes: &mut BTreeMap<(EdgeId, VertexId), VertexId>,
    step: &StepId,
    closure: catallax_recipe::FilletCorner,
    dec: &mut Decider,
) -> Result<Option<Continuation>, Refusal> {
    if closure == catallax_recipe::FilletCorner::Smooth {
        return super::fillet_vertex_smooth::build(
            model,
            body,
            old,
            runs,
            incident,
            uses,
            edit,
            end_changes,
            step,
            dec,
        )
        .map(Some);
    }
    let [first, second] = runs;
    let PlanarPair {
        faces: all_faces,
        common,
        outer,
        wall,
    } = planar_pair(model, body, [&*first, &*second], incident, uses, dec)?;
    let contact = line_line(rail(first, common)?, rail(second, common)?, dec)?;
    let normal = unit(&outward(model, common, &contact, dec)?);
    let radius = first.definition.radius.clone();
    let centre = contact.plus(&normal.scaled(&(-&radius)));
    let at = model.geom().position(model.vertex(old).point);
    let normals = [
        unit(&outward(model, outer[0], &at, dec)?),
        normal.clone(),
        unit(&outward(model, outer[1], &at, dec)?),
    ];
    let cross = normals[0].cross(&normals[2]);
    let turn = sign(&cross.dot(&normals[1]), dec)?;
    let cliff_normal = match turn {
        Sign::Positive => cross.clone(),
        Sign::Negative => cross.scaled(&-Real::one()),
        Sign::Zero => return Err(bad("the active corner normals must be independent")),
    };
    // A stationary cliff lets the final cell represent the whole homotopy.
    // At either other cone boundary a point outside the ball is already
    // consumed by that selected strip. A moving cliff needs its swept set.
    if sign(&at.minus(&centre).dot(&cliff_normal), dec)? != Sign::Zero {
        return Err(bad(
            "a moving normal-cone cliff needs its swept affected region",
        ));
    }
    let cliff = Plane {
        origin: at.clone(),
        normal: cliff_normal.clone(),
        ref_dir: normals[0].clone(),
    };
    let centre = plane_coordinates(centre, &cliff, dec)?;
    let antipode = centre.plus(&normal.scaled(&(-&radius)));
    let sum = normals.iter().fold(Vector3::zero(), |sum, n| sum.plus(n));
    for n in &normals {
        if sign(&sum.dot(n), dec)? != Sign::Positive {
            return Err(bad("the active corner normals need one open hemisphere"));
        }
    }
    for face in &all_faces {
        let p = plane(model, *face)?;
        let n = unit(&outward(model, *face, &at, dec)?);
        let depth = &centre.minus(&p.origin).dot(&n) + &radius;
        let side = sign(&depth, dec)?;
        if side == Sign::Positive
            || (first.faces.contains(face) || second.faces.contains(face)) && side != Sign::Zero
        {
            return Err(bad(
                "the active corner ball must lie inside every input support",
            ));
        }
    }
    let cliff_surface = model.geom_mut().add_surface(Surface::Plane(cliff.clone()));
    let sphere = Surface::Sphere(Sphere::new(
        centre.clone(),
        radius.clone(),
        normals[0].clone(),
        normals[1].minus(&normals[0].scaled(&normals[1].dot(&normals[0]))),
        dec,
    )?);
    for run in [&*first, &*second] {
        let Surface::Cylinder(cylinder) = &run.definition.surface else {
            return Err(bad("a planar junction strip must reduce to a cylinder"));
        };
        if sign(
            &centre
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
        {
            return Err(bad(
                "the two strip axes do not share the common-contact ball centre",
            ));
        }
        for face in &all_faces {
            if !run.faces.contains(face)
                && sign(&cylinder.axis.dot(&plane(model, *face)?.normal), dec)? == Sign::Zero
                && !run
                    .definition
                    .avoids_support(model.geom().surface(model.face(*face).surface), dec)?
            {
                return Err(bad(
                    "a parallel input support needs a lower-patch extremum proof",
                ));
            }
        }
        let input = model.edge(run.edge);
        let far = if input.start == Some(old) {
            input.end
        } else {
            input.start
        }
        .ok_or_else(|| bad("a four-face junction needs finite input edges"))?;
        let along = model.geom().position(model.vertex(far).point).minus(&at);
        if sign(&cylinder.axis.dot(&plane(model, wall)?.normal), dec)? == Sign::Zero
            || sign(&outward(model, wall, &at, dec)?.dot(&along), dec)? != Sign::Negative
        {
            return Err(bad(
                "the end wall must be transverse and face both remaining runs",
            ));
        }
    }
    let mut boundaries = [None; 2];
    for (edge, adjacent) in uses {
        let e = model.edge(*edge);
        if *edge == first.edge
            || *edge == second.edge
            || (e.start != Some(old) && e.end != Some(old))
        {
            continue;
        }
        if adjacent.len() != 2 || !adjacent.iter().any(|(face, _)| *face == wall) {
            return Err(bad("an unselected junction edge must meet the end wall"));
        }
        let own = outer.map(|face| adjacent.iter().any(|(f, _)| *f == face));
        let i = match own {
            [true, false] => 0,
            [false, true] => 1,
            _ => return Err(bad("the end-wall fan does not pair the two outer supports")),
        };
        if boundaries[i].replace(*edge).is_some() {
            return Err(bad("an outer support has more than one end-wall edge"));
        }
    }
    let [Some(a_edge), Some(b_edge)] = boundaries else {
        return Err(bad("both outer supports must meet the end wall"));
    };
    let a = wall_foot(model, uses, first, outer[0], wall, old, dec)?;
    let b = wall_foot(model, uses, second, outer[1], wall, old, dec)?;
    let surfaces = [
        first.definition.surface.clone(),
        second.definition.surface.clone(),
    ]
    .map(|surface| model.geom_mut().add_surface(surface));
    let meeting =
        catallax_intersect::intersect_surfaces(model.geom_mut(), surfaces[0], surfaces[1], dec)?
            .value;
    let SurfacePair::Meets(components) = meeting else {
        return Err(bad("the two strip cylinders have no common junction curve"));
    };
    let wall_surface = model.face(wall).surface;
    let contact_vertex = vertex(model, contact);
    let contact_point = model.vertex(contact_vertex).point;
    let mut junction = None;
    for component in components {
        let Component::Curve { curve, .. } = component else {
            continue;
        };
        if !matches!(model.geom().curve(curve), Curve::Ellipse(_)) {
            return Err(bad(
                "equal intersecting strip axes must reduce to ellipse junctions",
            ));
        }
        if decided(catallax_intersect::on_curve(
            model.geom(),
            curve,
            contact_point,
            dec,
        ))?
        .is_some()
        {
            continue;
        }
        let meeting = catallax_intersect::intersect_curve_surface(
            model.geom_mut(),
            curve,
            wall_surface,
            dec,
        )?
        .value;
        let CurveSurface::Meets(hits) = meeting else {
            continue;
        };
        for hit in hits {
            let p = model.geom().position(hit.point);
            if !first.definition.in_sector(&p, dec)? || !second.definition.in_sector(&p, dec)? {
                continue;
            }
            let mut inside = true;
            for face in &all_faces {
                let plane = plane(model, *face)?;
                let outside = if model.face(*face).sense == Sense::Same {
                    Sign::Positive
                } else {
                    Sign::Negative
                };
                if sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? == outside {
                    inside = false;
                    break;
                }
            }
            if !inside {
                continue;
            }
            if same(&p, &model.geom().position(contact_point), dec)?
                || same(&p, &a, dec)?
                || same(&p, &b, dec)?
                || junction.replace((curve, hit.point)).is_some()
            {
                return Err(bad(
                    "the two-strip end needs one strict junction in both sectors",
                ));
            }
        }
    }
    let (curve, point) =
        junction.ok_or_else(|| bad("no end-wall junction lies in both strip sectors"))?;
    let hub = model.add_vertex(Vertex { point });
    let avoid = model.geom_mut().add_point(PointDef::Explicit(antipode));
    let sense = match decided(catallax_intersect::cyclic_order(
        model.geom(),
        curve,
        [contact_point, avoid, point],
        dec,
    ))? {
        Some(Cyclic::Forward) => Sense::Reversed,
        Some(Cyclic::Backward) => Sense::Same,
        None => return Err(bad("the two-strip junction has a degenerate half-ellipse")),
    };
    let cosine = normals[0].dot(&normals[2]);
    let duals = [
        normals[0].minus(&normals[2].scaled(&cosine)),
        normals[2].minus(&normals[0].scaled(&cosine)),
    ];
    let meeting =
        catallax_intersect::intersect_curve_surface(model.geom_mut(), curve, cliff_surface, dec)?
            .value;
    let CurveSurface::Meets(hits) = meeting else {
        return Err(bad(
            "the strip junction does not cross the normal-cone cliff",
        ));
    };
    let mut crossing = None;
    for hit in hits {
        let p = model.geom().position(hit.point);
        if same(&p, &model.geom().position(contact_point), dec)?
            || same(&p, &model.geom().position(point), dec)?
        {
            continue;
        }
        let order = decided(catallax_intersect::cyclic_order(
            model.geom(),
            curve,
            [contact_point, hit.point, point],
            dec,
        ))?;
        if !matches!(
            (sense, order),
            (Sense::Same, Some(Cyclic::Forward)) | (Sense::Reversed, Some(Cyclic::Backward))
        ) {
            continue;
        }
        let delta = p.minus(&centre);
        if sign(&delta.dot(&duals[0]), dec)? != Sign::Positive
            || sign(&delta.dot(&duals[1]), dec)? != Sign::Positive
            || sign(&(&delta.norm_squared() - &radius.square()), dec)? != Sign::Positive
        {
            continue;
        }
        if crossing.replace(hit.point).is_some() {
            return Err(bad(
                "the normal-cone cliff needs one strict junction crossing",
            ));
        }
    }
    let q = model.add_vertex(Vertex {
        point: crossing.ok_or_else(|| bad("no exposed normal-cone cliff crossing"))?,
    });
    let shared = model.add_edge(Edge {
        curve,
        sense,
        start: Some(q),
        end: Some(hub),
    });
    contained_arc(
        model,
        &all_faces,
        &[(cliff_surface, Sign::Positive)],
        shared,
        dec,
    )?;
    let cone_sides = duals.map(|normal| {
        let surface = model.geom_mut().add_surface(Surface::Plane(Plane {
            origin: centre.clone(),
            normal,
            ref_dir: cliff_normal.clone(),
        }));
        (surface, Sign::Negative)
    });
    let feet = [vertex(model, a), vertex(model, b)];
    let sphere_feet = [normals[0].clone(), normals[2].clone()]
        .map(|n| vertex(model, centre.plus(&n.scaled(&radius))));
    let name = model
        .vertex_name(body, old)
        .ok_or_else(|| bad("a partial corner needs an unambiguous input vertex name"))?;
    let origin = Origin {
        step: step.clone(),
        role: Role::Corner { vertex: name },
    };
    let mut wall_caps = Vec::with_capacity(2);
    let mut cliff_caps = Vec::with_capacity(2);
    let mut sphere_caps = Vec::with_capacity(3);
    let mut shared_use = None;
    for (i, run) in [&mut *first, &mut *second].into_iter().enumerate() {
        let foot = feet[i];
        let sphere_foot = sphere_feet[i];
        let delta = model
            .geom()
            .position(model.vertex(foot).point)
            .minus(&centre);
        if sign(&delta.dot(&cliff_normal), dec)? != Sign::Negative {
            return Err(bad("the end wall must leave a strict lower cylinder patch"));
        }
        let avoid = line_plane(rail(run, common)?, plane(model, wall)?, dec)?;
        let support = model.geom().surface(model.face(wall).surface).clone();
        let (from, to) = if i == 0 { (foot, hub) } else { (hub, foot) };
        let wall_cap = connector(
            model,
            &run.definition.surface,
            &support,
            from,
            to,
            &[&run.definition],
            Some(&avoid),
            dec,
        )?;
        contained_arc(
            model,
            &all_faces,
            &[(cliff_surface, Sign::Positive)],
            wall_cap,
            dec,
        )?;
        wall_caps.push((wall_cap, Sense::Same));

        let avoid = line_plane(rail(run, common)?, &cliff, dec)?;
        let cliff_cap = connector(
            model,
            &run.definition.surface,
            &Surface::Plane(cliff.clone()),
            sphere_foot,
            q,
            &[&run.definition],
            Some(&avoid),
            dec,
        )?;
        contained_arc(model, &all_faces, &cone_sides, cliff_cap, dec)?;

        let input = model.edge(run.edge).clone();
        let index = run.faces.iter().position(|face| *face == outer[i]).unwrap();
        let curve = model
            .geom_mut()
            .add_curve(run.definition.contacts[index].clone());
        let (start, end) = if input.start == Some(old) {
            (foot, sphere_foot)
        } else {
            (sphere_foot, foot)
        };
        let rail = model.add_edge(Edge {
            curve,
            sense: input.sense,
            start: Some(start),
            end: Some(end),
        });
        let rail_sense = reverse(run.senses[index]);
        let from = edge_blend::ends(model, rail, rail_sense).1.unwrap();
        let rest = chain(model, &[cliff_cap, shared, wall_cap], from)?;
        let mut path = Vec::with_capacity(4);
        path.push((rail, rail_sense));
        for (edge, sense) in rest {
            if edge == cliff_cap {
                cliff_caps.push((edge, reverse(sense)));
            } else if edge == shared {
                if shared_use.is_some_and(|previous| previous == sense) {
                    return Err(bad(
                        "the lower patches must use their common ellipse oppositely",
                    ));
                }
                shared_use = Some(sense);
            }
            path.push((edge, sense));
        }
        let face = crate::extrude::face(
            model,
            run.definition.surface.clone(),
            Sense::Same,
            vec![path],
            None,
            vec![run.definition.origin.clone()],
        );
        edit.faces.push((outer[i], face));
        edit.gaps.insert((outer[i], old), Gap::Edge(rail));
        for face in run.faces {
            run.ends.insert(
                (old, face),
                if face == common {
                    contact_vertex
                } else {
                    sphere_foot
                },
            );
        }
        let n = &normals[if i == 0 { 0 } else { 2 }];
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: centre.clone(),
            normal: n.cross(&normals[1]),
            radius_sq: radius.square(),
            ref_dir: n.clone(),
        }));
        let cap = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(sphere_foot),
            end: Some(contact_vertex),
        });
        let at_first = if reverse(run.senses[0]) == Sense::Same {
            input.end
        } else {
            input.start
        };
        let from = run.ends[&(old, run.faces[if at_first == Some(old) { 0 } else { 1 }])];
        sphere_caps.push((
            cap,
            if from == sphere_foot {
                Sense::Reversed
            } else {
                Sense::Same
            },
        ));
        run.caps.insert(old, vec![cap]);
    }
    let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
        center: centre,
        normal: cross,
        radius_sq: radius.square(),
        ref_dir: normals[0].clone(),
    }));
    let arc = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(sphere_feet[0]),
        end: Some(sphere_feet[1]),
    });
    let arc_sense =
        if edge_blend::ends(model, sphere_caps[0].0, sphere_caps[0].1).0 == Some(sphere_feet[0]) {
            Sense::Reversed
        } else {
            Sense::Same
        };
    sphere_caps.push((arc, arc_sense));
    cliff_caps.push((arc, reverse(arc_sense)));
    let path = close(model, sphere_caps)?;
    let face = crate::extrude::face(
        model,
        sphere,
        Sense::Same,
        vec![path],
        None,
        vec![origin.clone()],
    );
    edit.faces.push((outer[0], face));
    let path = close(model, cliff_caps)?;
    let face = crate::extrude::face(
        model,
        Surface::Plane(cliff),
        Sense::Same,
        vec![path],
        None,
        vec![origin],
    );
    edit.faces.push((outer[0], face));
    edit.gaps.insert((wall, old), Gap::Path(wall_caps));
    end_changes.insert((a_edge, old), feet[0]);
    end_changes.insert((b_edge, old), feet[1]);
    Ok(None)
}
