//! Planar four-contact normal fans compiled without B-rep topology.

use super::*;
use catallax_features::fillet_corner::{FanTopology, PlanarFan, SphereSector};
use catallax_shapes::{Line, Sphere};

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-corner", detail)
}

fn cone_contains(
    normals: &[Vector3],
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for a in 0..normals.len() {
        for b in a + 1..normals.len() {
            for c in b + 1..normals.len() {
                let determinant =
                    geometry::sign(&normals[a].dot(&normals[b].cross(&normals[c])), dec)?;
                if determinant == Sign::Zero {
                    continue;
                }
                let mut inside = true;
                for [first, second] in [[b, c], [c, a], [a, b]] {
                    let sign = geometry::sign(
                        &direction.dot(&normals[first].cross(&normals[second])),
                        dec,
                    )?;
                    if sign != Sign::Zero && sign != determinant {
                        inside = false;
                        break;
                    }
                }
                if inside {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

fn certify_cells(
    vertex: &Point3,
    centres: [Option<&Point3>; 2],
    normals: &[Vector3],
    input: &[Plane],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let half = Real::one() / Real::from(2);
    let balls = centres.map(|centre| {
        centre.map(|centre| {
            let h = vertex.minus(centre);
            (
                centre.plus(&h.scaled(&half)),
                &h.norm_squared() / &Real::from(4),
            )
        })
    });
    // For d in a contact cone, the source inequalities give
    // |d|² <= d·(vertex-centre). The whole affected cell therefore lies in
    // this diametral ball. Bridge cells lie in the convex hull of both balls.
    for plane in input {
        match geometry::sign(&vertex.minus(&plane.origin).dot(&plane.normal), dec)? {
            Sign::Positive => {
                return Err(unsupported(
                    "a fan endpoint lies outside its certified input",
                ));
            }
            Sign::Zero => {
                if !cone_contains(normals, &plane.normal, dec)? {
                    return Err(unsupported(
                        "an unselected incident support needs another fan cell",
                    ));
                }
                continue;
            }
            Sign::Negative => {}
        }
        let norm_squared = plane.normal.norm_squared();
        for (midpoint, radius_squared) in balls.iter().flatten() {
            let clearance = plane.origin.minus(midpoint).dot(&plane.normal);
            if geometry::sign(&clearance, dec)? != Sign::Positive
                || geometry::sign(
                    &(&(&clearance * &clearance) - &(radius_squared * &norm_squared)),
                    dec,
                )? != Sign::Positive
            {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "a non-incident input wall is not proved clear of the fan's affected cells",
                ));
            }
        }
    }
    Ok(())
}

fn cap(
    definition: &FilletDefinition,
    interval: &Option<[Rational; 2]>,
    vertex: &Point3,
    centre: &Point3,
    dec: &mut Decider,
) -> Result<Surface, Refusal> {
    let (Curve::Line(edge), Surface::Cylinder(cylinder), Some(interval)) =
        (&definition.edge, &definition.surface, interval)
    else {
        return Err(unsupported("a fan needs finite straight cylindrical runs"));
    };
    let ends = interval
        .each_ref()
        .map(|t| edge.origin.plus(&edge.dir.scaled(&Real::from(t.clone()))));
    let far = if geometry::sign(&vertex.minus(&ends[0]).norm_squared(), dec)? == Sign::Zero {
        &ends[1]
    } else if geometry::sign(&vertex.minus(&ends[1]).norm_squared(), dec)? == Sign::Zero {
        &ends[0]
    } else {
        return Err(unsupported(
            "a fan endpoint is not an endpoint of its selected run",
        ));
    };
    let along = far.minus(vertex);
    let station = centre.minus(vertex).dot(&along);
    if geometry::sign(&station, dec)? != Sign::Positive
        || geometry::sign(&(&along.norm_squared() - &station), dec)? != Sign::Positive
    {
        return Err(Refusal::new(
            "fillet.unsupported-overflow",
            "a fan endpoint consumes its selected run",
        ));
    }
    let normal = match geometry::sign(&cylinder.axis.dot(&along), dec)? {
        Sign::Positive => cylinder.axis.scaled(&-Real::one()),
        Sign::Negative => cylinder.axis.clone(),
        Sign::Zero => return Err(unsupported("a fan cap needs a nonzero axial interval")),
    };
    let Surface::Plane(support) = &definition.supports[0] else {
        return Err(unsupported("a planar fan needs planar supports"));
    };
    Ok(Surface::Plane(Plane {
        origin: centre.clone(),
        normal,
        ref_dir: support.normal.clone(),
    }))
}

fn sector_planes(sector: &SphereSector) -> Vec<Surface> {
    sector
        .dual_basis()
        .into_iter()
        .enumerate()
        .map(|(i, dual)| {
            Surface::Plane(Plane {
                origin: sector.centre.clone(),
                normal: dual.scaled(&-Real::one()),
                ref_dir: sector.directions[(i + 1) % 3].clone(),
            })
        })
        .collect()
}

impl Plan {
    pub(super) fn network_fillet(
        &mut self,
        root: usize,
        source_faces: &[usize],
        definitions: &[(FilletDefinition, Option<[Rational; 2]>)],
        endpoints: &[(Point3, Vec<usize>)],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        if !endpoints.iter().any(|(_, runs)| runs.len() == 4) {
            return Ok(None);
        }
        let input = self.convex_planar_input(root, dec)?;
        if definitions.iter().any(|(definition, _)| {
            !definition.convex
                || !matches!(definition.surface, Surface::Cylinder(_))
                || definition
                    .supports
                    .iter()
                    .any(|support| !matches!(support, Surface::Plane(_)))
        }) {
            return Err(unsupported(
                "field contact fans need convex planar cylindrical runs",
            ));
        }
        let mut caps: Vec<Vec<Surface>> = vec![Vec::new(); definitions.len()];
        let mut regions = Vec::new();
        let mut bridges = Vec::new();
        for (vertex, runs) in endpoints {
            if runs.len() == 1 {
                continue;
            }
            if runs.len() == 2 {
                return Err(unsupported(
                    "a contact network with a partial junction needs its own cell proof",
                ));
            }
            let mut supports: Vec<(&Surface, bool)> = Vec::with_capacity(runs.len());
            let mut pairs = Vec::with_capacity(runs.len());
            for index in runs {
                let definition = &definitions[*index].0;
                let mut pair = [0; 2];
                for (slot, support) in definition.supports.iter().enumerate() {
                    pair[slot] = if let Some(index) = supports
                        .iter()
                        .position(|(other, _)| support.same_definition(other))
                    {
                        if supports[index].1 != definition.same_sense[slot] {
                            return Err(unsupported(
                                "one source support has inconsistent material sides",
                            ));
                        }
                        index
                    } else {
                        supports.push((support, definition.same_sense[slot]));
                        supports.len() - 1
                    };
                }
                pairs.push(pair);
            }
            if supports.len() != runs.len() {
                return Err(unsupported(
                    "a fully selected fan needs one exposed run per incident support",
                ));
            }
            let mut names = Vec::with_capacity(supports.len());
            let mut directions = Vec::with_capacity(supports.len());
            for (support, same) in &supports {
                let face = source_faces
                    .iter()
                    .map(|i| &self.faces[*i])
                    .find(|face| face.surface.same_definition(support))
                    .ok_or_else(|| unsupported("a contact fan lost an input support name"))?;
                names.push(face.origin.clone());
                let Surface::Plane(plane) = support else {
                    unreachable!()
                };
                let normal = unit(&plane.normal);
                directions.push(if *same {
                    normal
                } else {
                    normal.scaled(&-Real::one())
                });
            }
            names.sort();
            names.dedup();
            let origin = Origin {
                step: step.clone(),
                role: Role::Corner {
                    vertex: VertexName::new(names).map_err(|error| bad(&error.to_string()))?,
                },
            };
            if let [a, b, c] = runs.as_slice() {
                let corner = CornerDefinition::new(
                    origin.clone(),
                    [*a, *b, *c].map(|i| &definitions[i].0),
                    dec,
                )?;
                certify_cells(
                    vertex,
                    [Some(corner.centre()), None],
                    &directions,
                    &input,
                    dec,
                )?;
                for index in runs {
                    caps[*index].push(cap(
                        &definitions[*index].0,
                        &definitions[*index].1,
                        vertex,
                        corner.centre(),
                        dec,
                    )?);
                }
                let sector = corner.sector;
                regions.push((sector_planes(&sector), sector.surface, origin));
                continue;
            }
            let indices: [usize; 4] = runs.as_slice().try_into().map_err(|_| {
                unsupported("a field contact fan needs three or four selected runs")
            })?;
            let directions: [Vector3; 4] = directions
                .try_into()
                .map_err(|_| unsupported("a field contact fan needs four supports"))?;
            let pairs: [[usize; 2]; 4] = pairs
                .try_into()
                .map_err(|_| unsupported("a field contact fan needs four exposed pairs"))?;
            let radius = definitions[indices[0]].0.radius.clone();
            for index in indices {
                if geometry::sign(&(&definitions[index].0.radius - &radius), dec)? != Sign::Zero {
                    return Err(unsupported("one contact network needs a common radius"));
                }
            }
            let fan = PlanarFan::new(vertex.clone(), directions, pairs, radius, dec)?;
            match &fan.topology {
                FanTopology::Single { centre } => {
                    certify_cells(vertex, [Some(centre), None], &fan.directions, &input, dec)?;
                    for index in indices {
                        caps[index].push(cap(
                            &definitions[index].0,
                            &definitions[index].1,
                            vertex,
                            centre,
                            dec,
                        )?);
                    }
                    let mut planes = Vec::with_capacity(4);
                    for [a, b] in fan.pairs {
                        let mut normal = fan.directions[a].cross(&fan.directions[b]);
                        let other = (0..4).find(|i| *i != a && *i != b).unwrap();
                        if geometry::sign(&normal.dot(&fan.directions[other]), dec)?
                            == Sign::Positive
                        {
                            normal = normal.scaled(&-Real::one());
                        }
                        planes.push(Surface::Plane(Plane {
                            origin: centre.clone(),
                            normal,
                            ref_dir: fan.directions[a].clone(),
                        }));
                    }
                    let reference = fan.directions[1].minus(
                        &fan.directions[0].scaled(&fan.directions[0].dot(&fan.directions[1])),
                    );
                    let sphere = Sphere::new(
                        centre.clone(),
                        fan.radius.clone(),
                        fan.directions[0].clone(),
                        reference,
                        dec,
                    )?;
                    regions.push((planes, Surface::Sphere(sphere), origin));
                }
                FanTopology::Split {
                    vertices,
                    bridge,
                    endpoints,
                } => {
                    certify_cells(
                        vertex,
                        vertices.each_ref().map(|node| Some(&node.centre)),
                        &fan.directions,
                        &input,
                        dec,
                    )?;
                    for (slot, index) in indices.into_iter().enumerate() {
                        caps[index].push(cap(
                            &definitions[index].0,
                            &definitions[index].1,
                            vertex,
                            &vertices[endpoints[slot]].centre,
                            dec,
                        )?);
                    }
                    for node in vertices {
                        let sector = SphereSector::new(
                            node.centre.clone(),
                            fan.radius.clone(),
                            node.contacts.map(|i| fan.directions[i].clone()),
                            true,
                            dec,
                        )?;
                        regions.push((sector_planes(&sector), sector.surface, origin.clone()));
                    }
                    let definition = FilletDefinition::new(
                        origin,
                        bridge.map(|i| supports[i].0.clone()),
                        bridge.map(|i| supports[i].1),
                        Curve::Line(Line {
                            origin: vertex.clone(),
                            dir: fan.directions[bridge[0]].cross(&fan.directions[bridge[1]]),
                        }),
                        fan.radius.clone(),
                        true,
                        dec,
                    )?;
                    let axis = vertices[1].centre.minus(&vertices[0].centre);
                    let caps = [0, 1].map(|i| {
                        Surface::Plane(Plane {
                            origin: vertices[i].centre.clone(),
                            normal: if i == 0 {
                                axis.scaled(&-Real::one())
                            } else {
                                axis.clone()
                            },
                            ref_dir: fan.directions[bridge[0]].clone(),
                        })
                    });
                    bridges.push((definition, caps));
                }
            }
        }
        let mut result = root;
        for ((definition, _), caps) in definitions.iter().zip(caps) {
            if let [Surface::Plane(first), Surface::Plane(second)] = caps.as_slice() {
                if geometry::sign(&first.normal.dot(&second.normal), dec)? != Sign::Negative
                    || geometry::sign(&second.origin.minus(&first.origin).dot(&first.normal), dec)?
                        != Sign::Negative
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "contact fan caps consume their shared strip",
                    ));
                }
            } else if caps.len() > 2 {
                return Err(unsupported(
                    "a straight selected run has more than two endpoint caps",
                ));
            }
            let node = self.analytic_strip(definition, caps, dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        for (definition, caps) in bridges {
            let node = self.analytic_strip(&definition, Vec::from(caps), dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        for (planes, surface, origin) in regions {
            let node = self.analytic_region(planes, surface, origin, dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        Ok(Some(result))
    }
}
