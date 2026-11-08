//! Exact spherical normal cells and their finite internal contact cylinder.

use super::*;
use crate::fillet_corner::{FanTopology, PlanarFan, SphereSector};

pub(super) struct Input<'a> {
    pub body: BodyId,
    pub old: VertexId,
    pub active: &'a [usize],
    pub faces: Vec<FaceId>,
    pub directions: Vec<Vector3>,
    pub uses: &'a Uses,
    pub step: &'a StepId,
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new(UNSUPPORTED_CORNER, detail)
}

pub(super) fn build(
    model: &mut Model,
    runs: &mut [Run],
    edit: &mut Replacement,
    input: Input<'_>,
    dec: &mut Decider,
) -> Result<Continuation, Refusal> {
    let Input {
        body,
        old,
        active,
        faces,
        directions,
        uses,
        step,
    } = input;
    let faces: [FaceId; 4] = faces
        .try_into()
        .map_err(|_| unsupported("a split contact fan needs four input supports"))?;
    let directions: [Vector3; 4] = directions
        .try_into()
        .map_err(|_| unsupported("a split contact fan needs four material normals"))?;
    let indices: [usize; 4] = active
        .try_into()
        .map_err(|_| unsupported("a split contact fan needs four selected runs"))?;
    let at = model.geom().position(model.vertex(old).point);
    let radius = runs[indices[0]].definition.radius.clone();
    let convex = runs[indices[0]].definition.convex;
    let mut pairs = [[0; 2]; 4];
    for (slot, index) in indices.into_iter().enumerate() {
        let run = &runs[index];
        if run.definition.convex != convex
            || sign(&(&run.definition.radius - &radius), dec)? != Sign::Zero
        {
            return Err(unsupported(
                "a split contact fan needs one radius and material side",
            ));
        }
        for (side, face) in run.faces.into_iter().enumerate() {
            pairs[slot][side] = faces
                .iter()
                .position(|f| *f == face)
                .ok_or_else(|| unsupported("a selected run leaves its input fan"))?;
        }
    }
    let fan = PlanarFan::new(at.clone(), directions, pairs, radius, dec)?;
    let FanTopology::Split {
        vertices,
        bridge,
        endpoints,
    } = &fan.topology
    else {
        return Err(unsupported(
            "the split-fan path requires two distinct offset vertices",
        ));
    };
    for (slot, index) in indices.into_iter().enumerate() {
        let run = &runs[index];
        let centre = &vertices[endpoints[slot]].centre;
        let (Curve::Line(edge), Surface::Cylinder(cylinder)) =
            (&run.definition.edge, &run.definition.surface)
        else {
            return Err(unsupported(
                "a planar contact fan needs straight cylindrical runs",
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
            return Err(unsupported(
                "a selected strip axis misses its offset endpoint",
            ));
        }
        let input = model.edge(run.edge);
        let far = if input.start == Some(old) {
            input.end
        } else {
            input.start
        }
        .ok_or_else(|| unsupported("a contact fan needs finite selected runs"))?;
        let far = model.geom().position(model.vertex(far).point);
        let length = far.minus(&at).dot(&edge.dir);
        let station = centre.minus(&at).dot(&edge.dir);
        let remaining = &length - &station;
        if sign(&(&station * &length), dec)? != Sign::Positive
            || sign(&(&remaining * &length), dec)? != Sign::Positive
        {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "an offset junction endpoint is not strictly inside its selected run",
            ));
        }
    }
    let origin = Origin {
        step: step.clone(),
        role: Role::Corner {
            vertex: model
                .vertex_name(body, old)
                .ok_or_else(|| unsupported("a contact network needs its original vertex name"))?,
        },
    };
    let axis = fan.directions[bridge[0]].cross(&fan.directions[bridge[1]]);
    let definition = FilletDefinition::new(
        origin.clone(),
        bridge.map(|i| model.geom().surface(model.face(faces[i]).surface).clone()),
        bridge.map(|i| model.face(faces[i]).sense == Sense::Same),
        Curve::Line(Line {
            origin: at.clone(),
            dir: axis,
        }),
        fan.radius.clone(),
        convex,
        dec,
    )?;
    let Surface::Cylinder(cylinder) = &definition.surface else {
        return Err(unsupported(
            "an internal planar contact sheet must be cylindrical",
        ));
    };
    for vertex in vertices {
        if sign(
            &vertex
                .centre
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
        {
            return Err(unsupported(
                "the internal sheet must pass through both ball centres",
            ));
        }
    }
    let continuation = Continuation {
        definition,
        region: RegionBound::Network {
            vertex: at,
            centres: vertices.each_ref().map(|v| v.centre.clone()),
            points: std::cell::OnceCell::new(),
        },
    };
    for face in model.faces_of_body(body) {
        if faces.contains(&face) {
            continue;
        }
        let support = model.geom().surface(model.face(face).surface);
        if !hull_clear(continuation.region_hull(), support, dec)? {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-incident face is not proved clear of the contact network",
            ));
        }
    }
    let sectors = vertices.each_ref().map(|v| {
        SphereSector::new(
            v.centre.clone(),
            fan.radius.clone(),
            v.contacts.map(|i| fan.directions[i].clone()),
            convex,
            dec,
        )
    });
    let [first, second] = sectors;
    let sectors = [first?, second?];
    let mut feet = [[None; 4]; 2];
    for (node, (v, sector)) in vertices.iter().zip(&sectors).enumerate() {
        for (face, point) in v.contacts.iter().zip(&sector.feet) {
            feet[node][*face] = Some(vertex(model, point.clone()));
        }
    }
    let mut sphere_caps = [Vec::with_capacity(3), Vec::with_capacity(3)];
    for (slot, index) in indices.into_iter().enumerate() {
        let node = endpoints[slot];
        let pair = fan.pairs[slot];
        let [a, b] = pair.map(|i| feet[node][i].unwrap());
        let run = &mut runs[index];
        for (face, foot) in run.faces.into_iter().zip([a, b]) {
            run.ends.insert((old, face), foot);
        }
        let cap = meridian_edge(
            model,
            Circle {
                center: vertices[node].centre.clone(),
                normal: fan.directions[pair[0]].cross(&fan.directions[pair[1]]),
                radius_sq: fan.radius.square(),
                ref_dir: fan.directions[pair[0]].clone(),
            },
            a,
            b,
            dec,
        )?;
        sphere_caps[node].push(corner_cap_use(model, run, old, cap)?);
        run.caps.insert(old, vec![cap]);
    }
    let mut caps = [None; 2];
    for node in 0..2 {
        caps[node] = Some(meridian_edge(
            model,
            Circle {
                center: vertices[node].centre.clone(),
                normal: fan.directions[bridge[0]].cross(&fan.directions[bridge[1]]),
                radius_sq: fan.radius.square(),
                ref_dir: fan.directions[bridge[0]].clone(),
            },
            feet[node][bridge[0]].unwrap(),
            feet[node][bridge[1]].unwrap(),
            dec,
        )?);
    }
    let caps = caps.map(Option::unwrap);
    let mut rails = [None; 2];
    let mut rail_senses = [Sense::Same; 2];
    for (side, support) in bridge.iter().enumerate() {
        let from = feet[0][*support].unwrap();
        let to = feet[1][*support].unwrap();
        let Curve::Line(line) = &continuation.definition.contacts[side] else {
            return Err(unsupported(
                "a planar internal sheet needs straight contact rails",
            ));
        };
        let along = model
            .geom()
            .position(model.vertex(to).point)
            .minus(&model.geom().position(model.vertex(from).point))
            .dot(&line.dir);
        let sense = match sign(&along, dec)? {
            Sign::Positive => Sense::Same,
            Sign::Negative => Sense::Reversed,
            Sign::Zero => return Err(unsupported("an internal contact rail has zero extent")),
        };
        let curve = model
            .geom_mut()
            .add_curve(continuation.definition.contacts[side].clone());
        rails[side] = Some(model.add_edge(Edge {
            curve,
            sense,
            start: Some(from),
            end: Some(to),
        }));
        let own = indices
            .iter()
            .enumerate()
            .find_map(|(slot, index)| {
                (endpoints[slot] == 0 && runs[*index].faces.contains(&faces[*support]))
                    .then_some(*index)
            })
            .ok_or_else(|| unsupported("an internal rail lacks its source-face gap"))?;
        rail_senses[side] =
            fillet_vertex_smooth::gap_sense(model, uses, runs[own].edge, faces[*support], old)?;
    }
    let [first_rail, second_rail] = rails.map(Option::unwrap);
    if rail_senses[1] != reverse(rail_senses[0]) {
        return Err(unsupported(
            "the internal rails do not oppose both source-face gaps",
        ));
    }
    let cylinder_loop = if rail_senses[0] == Sense::Same {
        vec![
            (first_rail, Sense::Same),
            (caps[1], Sense::Same),
            (second_rail, Sense::Reversed),
            (caps[0], Sense::Reversed),
        ]
    } else {
        vec![
            (first_rail, Sense::Reversed),
            (caps[0], Sense::Same),
            (second_rail, Sense::Same),
            (caps[1], Sense::Reversed),
        ]
    };
    for (edge, sense) in &cylinder_loop {
        for node in 0..2 {
            if *edge == caps[node] {
                sphere_caps[node].push((*edge, reverse(*sense)));
            }
        }
    }
    for (node, boundary) in sphere_caps.into_iter().enumerate() {
        let path = close(model, boundary)?;
        let face = crate::extrude::face(
            model,
            sectors[node].surface.clone(),
            if convex { Sense::Same } else { Sense::Reversed },
            vec![path],
            None,
            vec![origin.clone()],
        );
        edit.faces.push((faces[vertices[node].contacts[0]], face));
    }
    let face = crate::extrude::face(
        model,
        continuation.definition.surface.clone(),
        if convex { Sense::Same } else { Sense::Reversed },
        vec![cylinder_loop],
        None,
        vec![origin],
    );
    edit.faces.push((faces[bridge[0]], face));
    for (support, rail) in bridge.iter().zip([first_rail, second_rail]) {
        edit.gaps.insert((faces[*support], old), Gap::Edge(rail));
    }
    Ok(continuation)
}
