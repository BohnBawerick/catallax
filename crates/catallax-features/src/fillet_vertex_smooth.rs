//! Full-radius contact-sheet continuation at a partial planar junction.

use super::fillet_vertex_end::{PlanarPair, contained_arc, planar_pair};
use super::*;
use catallax_shapes::{Plane, Sphere};

pub(super) fn gap_sense(
    model: &Model,
    uses: &Uses,
    edge: EdgeId,
    face: FaceId,
    old: VertexId,
) -> Result<Sense, Refusal> {
    let sense = uses[&edge]
        .iter()
        .find_map(|(own, sense)| (*own == face).then_some(*sense))
        .ok_or_else(|| bad("a continuation boundary has no source-face use"))?;
    let (from, to) = edge_blend::ends(model, edge, sense);
    if from == Some(old) {
        Ok(Sense::Same)
    } else if to == Some(old) {
        Ok(Sense::Reversed)
    } else {
        Err(bad(
            "a continuation boundary does not meet its input vertex",
        ))
    }
}

/// The open A/C side of the A/B/C sphere joins the same-radius A/C cylinder.
/// Its axis points away from B and its finite end is the actual D plane.
/// Complete conic/halfspace events certify that end, not just its endpoints.
/// Projection along the transverse axis is one-to-one; each ruling between
/// the sphere cap and the end cap is affine in every input halfspace. Thus
/// the caps certify the whole cylinder face, including its interior.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    old: VertexId,
    runs: [&mut Run; 2],
    incident: &BTreeSet<FaceId>,
    uses: &Uses,
    edit: &mut Replacement,
    end_changes: &mut BTreeMap<(EdgeId, VertexId), VertexId>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Continuation, Refusal> {
    let [first, second] = runs;
    let PlanarPair {
        faces,
        common,
        outer,
        wall,
    } = planar_pair(model, body, [&*first, &*second], incident, uses, dec)?;
    let at = model.geom().position(model.vertex(old).point);
    let radius = first.definition.radius.clone();
    let normals = [
        unit(&outward(model, outer[0], &at, dec)?),
        unit(&outward(model, common, &at, dec)?),
        unit(&outward(model, outer[1], &at, dec)?),
    ];
    // The equal-depth bisectors pass through the old vertex. Solving them
    // together retains exact stationary coordinates, rather than introducing
    // an equivalent but noncanonical intersection of offset contact rails.
    let bisector = normals[1]
        .minus(&normals[0])
        .cross(&normals[2].minus(&normals[0]));
    let divisor = bisector.dot(&normals[0]);
    if sign(&divisor, dec)? == Sign::Zero {
        return Err(bad("the active corner normals must be independent"));
    }
    let centre = at.plus(&bisector.scaled(&(-&radius / divisor)));
    for face in &faces {
        let support = plane(model, *face)?;
        let n = unit(&outward(model, *face, &at, dec)?);
        let side = sign(&(&centre.minus(&support.origin).dot(&n) + &radius), dec)?;
        if side == Sign::Positive
            || (first.faces.contains(face) || second.faces.contains(face)) && side != Sign::Zero
        {
            return Err(bad(
                "the active corner ball must lie inside every input support",
            ));
        }
    }
    let mut axis = normals[0].cross(&normals[2]);
    match sign(&axis.dot(&normals[1]), dec)? {
        Sign::Negative => {}
        Sign::Positive => axis = axis.scaled(&-Real::one()),
        Sign::Zero => return Err(bad("the open contact sheet needs an outgoing axis")),
    }
    if sign(&axis.dot(&outward(model, wall, &at, dec)?), dec)? != Sign::Positive {
        return Err(bad(
            "the end wall must transversely close the outgoing contact sheet",
        ));
    }
    let origin = Origin {
        step: step.clone(),
        role: Role::Corner {
            vertex: model
                .vertex_name(body, old)
                .ok_or_else(|| bad("a smooth corner needs an unambiguous input vertex name"))?,
        },
    };
    let definition = FilletDefinition::new(
        origin.clone(),
        outer.map(|face| model.geom().surface(model.face(face).surface).clone()),
        outer.map(|face| model.face(face).sense == Sense::Same),
        Curve::Line(Line {
            origin: at.clone(),
            dir: axis.clone(),
        }),
        radius.clone(),
        true,
        dec,
    )?;
    for surface in [
        &first.definition.surface,
        &second.definition.surface,
        &definition.surface,
    ] {
        let Surface::Cylinder(cylinder) = surface else {
            return Err(bad(
                "a planar junction contact sheet must reduce to a cylinder",
            ));
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
                "the contact sheet axes must share the active ball centre",
            ));
        }
    }
    let sphere = Surface::Sphere(Sphere::new(
        centre.clone(),
        radius.clone(),
        normals[0].clone(),
        normals[1].minus(&normals[0].scaled(&normals[1].dot(&normals[0]))),
        dec,
    )?);
    let feet = normals
        .clone()
        .map(|n| vertex(model, centre.plus(&n.scaled(&radius))));
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
        let index = match own {
            [true, false] => 0,
            [false, true] => 1,
            _ => return Err(bad("the end-wall fan does not pair the two outer supports")),
        };
        if boundaries[index].replace(*edge).is_some() {
            return Err(bad("an outer support has more than one end-wall edge"));
        }
    }
    let [Some(a_edge), Some(c_edge)] = boundaries else {
        return Err(bad("both outer supports must meet the end wall"));
    };
    let boundaries = [a_edge, c_edge];
    let mut end_feet = [None; 2];
    for (index, edge) in boundaries.into_iter().enumerate() {
        let Curve::Line(line) = &definition.contacts[index] else {
            return Err(bad("a planar continuation needs straight contact rails"));
        };
        let point = line_plane(line, plane(model, wall)?, dec)?;
        let foot = vertex(model, point);
        if !on_trimmed(model, edge, model.vertex(foot).point, dec)? {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a continuation rail leaves its end-wall edge",
            ));
        }
        end_feet[index] = Some(foot);
    }
    let [a_end, c_end] = end_feet.map(Option::unwrap);
    let start_cap = meridian_edge(
        model,
        Circle {
            center: centre.clone(),
            normal: normals[0].cross(&normals[2]),
            radius_sq: radius.square(),
            ref_dir: normals[0].clone(),
        },
        feet[0],
        feet[2],
        dec,
    )?;
    let wall_support = model.geom().surface(model.face(wall).surface).clone();
    let end_cap = connector(
        model,
        &definition.surface,
        &wall_support,
        a_end,
        c_end,
        &[&definition],
        None,
        None,
        dec,
    )?;
    let dot = normals[0].dot(&normals[2]);
    let extra = [
        (axis.clone(), normals[0].clone()),
        (normals[0].minus(&normals[2].scaled(&dot)), axis.clone()),
        (normals[2].minus(&normals[0].scaled(&dot)), axis.clone()),
    ]
    .map(|(normal, ref_dir)| {
        let surface = model.geom_mut().add_surface(Surface::Plane(Plane {
            origin: centre.clone(),
            ref_dir,
            normal,
        }));
        (surface, Sign::Negative)
    });
    contained_arc(model, &faces, &extra, end_cap, dec)?;
    let mut rails = [None; 2];
    for (index, (from, to)) in [(feet[0], a_end), (feet[2], c_end)].into_iter().enumerate() {
        let Curve::Line(line) = &definition.contacts[index] else {
            unreachable!()
        };
        let sense = match sign(
            &model
                .geom()
                .position(model.vertex(to).point)
                .minus(&model.geom().position(model.vertex(from).point))
                .dot(&line.dir),
            dec,
        )? {
            Sign::Positive => Sense::Same,
            Sign::Negative => Sense::Reversed,
            Sign::Zero => return Err(bad("a continued contact rail must have positive extent")),
        };
        let curve = model
            .geom_mut()
            .add_curve(definition.contacts[index].clone());
        rails[index] = Some(model.add_edge(Edge {
            curve,
            sense,
            start: Some(from),
            end: Some(to),
        }));
    }
    let [a_rail, c_rail] = rails.map(Option::unwrap);
    let rail_sense = reverse(gap_sense(model, uses, a_edge, outer[0], old)?);
    let cylinder_loop = match rail_sense {
        Sense::Same => vec![
            (a_rail, Sense::Same),
            (end_cap, Sense::Same),
            (c_rail, Sense::Reversed),
            (start_cap, Sense::Reversed),
        ],
        Sense::Reversed => vec![
            (a_rail, Sense::Reversed),
            (start_cap, Sense::Same),
            (c_rail, Sense::Same),
            (end_cap, Sense::Reversed),
        ],
    };
    let mut sphere_caps = Vec::with_capacity(3);
    for (edge, sense) in cylinder_loop.iter().copied().skip(1) {
        if edge == start_cap {
            sphere_caps.push((edge, reverse(sense)));
        } else if edge == c_rail && sense != reverse(gap_sense(model, uses, c_edge, outer[1], old)?)
        {
            return Err(bad("continued rails must oppose both source-face gaps"));
        } else if edge == end_cap && sense != gap_sense(model, uses, a_edge, wall, old)? {
            return Err(bad("the continuation cap must oppose its end-wall gap"));
        }
    }
    for (index, run) in [&mut *first, &mut *second].into_iter().enumerate() {
        let outer_foot = feet[if index == 0 { 0 } else { 2 }];
        for face in run.faces {
            run.ends.insert(
                (old, face),
                if face == common { feet[1] } else { outer_foot },
            );
        }
        let cap = meridian_edge(
            model,
            Circle {
                center: centre.clone(),
                normal: normals[if index == 0 { 0 } else { 2 }].cross(&normals[1]),
                radius_sq: radius.square(),
                ref_dir: normals[if index == 0 { 0 } else { 2 }].clone(),
            },
            outer_foot,
            feet[1],
            dec,
        )?;
        sphere_caps.push(corner_cap_use(model, run, old, cap)?);
        run.caps.insert(old, vec![cap]);
    }
    let sphere_loop = close(model, sphere_caps)?;
    let sphere_face = crate::extrude::face(
        model,
        sphere,
        Sense::Same,
        vec![sphere_loop],
        None,
        vec![origin.clone()],
    );
    let cylinder_face = crate::extrude::face(
        model,
        definition.surface.clone(),
        Sense::Same,
        vec![cylinder_loop],
        None,
        vec![origin],
    );
    edit.faces
        .extend([(outer[0], sphere_face), (outer[0], cylinder_face)]);
    edit.gaps.insert((outer[0], old), Gap::Edge(a_rail));
    edit.gaps.insert((outer[1], old), Gap::Edge(c_rail));
    edit.gaps.insert((wall, old), Gap::Edge(end_cap));
    end_changes.insert((a_edge, old), a_end);
    end_changes.insert((c_edge, old), c_end);
    let on_axis = at.plus(&axis.scaled(&(centre.minus(&at).dot(&axis) / axis.norm_squared())));
    let p = |vertex_id| model.geom().position(model.vertex(vertex_id).point);
    Ok(Continuation {
        definition,
        region: RegionBound::Sheet {
            hull: [on_axis, p(feet[0]), p(feet[2]), at, p(a_end), p(c_end)],
            centre,
            points: std::cell::OnceCell::new(),
        },
    })
}
