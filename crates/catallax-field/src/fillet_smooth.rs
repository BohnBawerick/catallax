//! Algebraic normal-cone cells for smooth planar four-face junctions.

use super::*;
use catallax_features::fillet_corner::SphereSector;
use catallax_shapes::Line;

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-corner", detail)
}

impl Plan {
    pub(super) fn convex_planar_input(
        &self,
        root: usize,
        dec: &mut Decider,
    ) -> Result<Vec<Plane>, Refusal> {
        match self.instructions[root] {
            Instruction::Leaf(index) => {
                let leaf = &self.leaves[index];
                if let geometry::Leaf::Feature(feature) = leaf {
                    return feature.convex_planes(dec)?.ok_or_else(|| {
                        unsupported("smooth field junctions need a certified convex planar prism")
                    });
                }
                if matches!(leaf, geometry::Leaf::Draft(_)) {
                    return Err(unsupported(
                        "a draft junction needs its retained outward halfspaces",
                    ));
                }
                leaf.certify_planar_convexity(dec).map_err(|error| {
                    if error.code == "draft.unsupported-arrangement" {
                        unsupported("smooth field junctions need a certified convex planar input")
                    } else {
                        error
                    }
                })?;
                leaf.constraints()
                    .iter()
                    .map(|constraint| {
                        let Surface::Plane(plane) = constraint.surface() else {
                            return Err(unsupported(
                                "smooth junction input constraints must be planar",
                            ));
                        };
                        Ok(plane.clone())
                    })
                    .collect()
            }
            Instruction::Boolean {
                kind: BooleanKind::Intersect,
                a,
                b,
            } => {
                let mut planes = self.convex_planar_input(a, dec)?;
                planes.extend(self.convex_planar_input(b, dec)?);
                Ok(planes)
            }
            _ => Err(unsupported(
                "smooth field junctions need a certified convex planar input",
            )),
        }
    }

    pub(super) fn analytic_region(
        &mut self,
        planes: Vec<Surface>,
        surface: Surface,
        origin: Origin,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        // These auxiliary cones and cylinders are unbounded. The original
        // finite source clips the result; do not assign a fictitious bound.
        let cone = self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(planes, None)?),
            Vec::new(),
            dec,
        )?;
        let face = geometry::face(surface.clone(), origin, dec)?;
        let interior = self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(vec![surface], None)?),
            vec![face],
            dec,
        )?;
        self.boolean(BooleanKind::Subtract, cone, interior, dec)
    }

    pub(super) fn smooth_fillet(
        &mut self,
        root: usize,
        source_faces: &[usize],
        definitions: &[(FilletDefinition, Option<[Rational; 2]>)],
        endpoints: &[(Point3, Vec<usize>)],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let mut caps: Vec<Vec<Surface>> = Vec::new();
        let mut continuations = Vec::new();
        let mut sectors = Vec::new();
        let mut input_planes = None;
        for (at, runs) in endpoints {
            let [first_index, second_index] = runs.as_slice() else {
                continue;
            };
            let mut incident = Vec::new();
            for &index in source_faces {
                let Surface::Plane(plane) = &self.faces[index].surface else {
                    continue;
                };
                if geometry::sign(&at.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Zero {
                    incident.push(index);
                }
            }
            if incident.len() != 4 {
                return Err(unsupported(
                    "the smooth junction lacks a certified terminating contact sheet",
                ));
            }
            if input_planes.is_none() {
                input_planes = Some(self.convex_planar_input(root, dec)?);
            }
            let planes = input_planes.as_ref().unwrap();
            let first = &definitions[*first_index].0;
            let second = &definitions[*second_index].0;
            if !first.convex
                || !second.convex
                || first.same_sense != [true; 2]
                || second.same_sense != [true; 2]
                || geometry::sign(&(&first.radius - &second.radius), dec)? != Sign::Zero
            {
                return Err(unsupported(
                    "smooth field junctions need equal convex planar contacts",
                ));
            }
            let shared = first
                .supports
                .iter()
                .position(|a| second.supports.iter().any(|b| a.same_definition(b)))
                .ok_or_else(|| unsupported("two smooth strips must share one support"))?;
            let second_shared = second
                .supports
                .iter()
                .position(|s| s.same_definition(&first.supports[shared]))
                .unwrap();
            let supports = [
                first.supports[1 - shared].clone(),
                first.supports[shared].clone(),
                second.supports[1 - second_shared].clone(),
            ];
            let [a, b, c] = supports.each_ref().map(|support| match support {
                Surface::Plane(plane) => Ok(unit(&plane.normal)),
                _ => Err(unsupported(
                    "smooth field junctions need three independent planes",
                )),
            });
            let normals = [a?, b?, c?];
            let bisector = normals[1]
                .minus(&normals[0])
                .cross(&normals[2].minus(&normals[0]));
            let divisor = bisector.dot(&normals[0]);
            if geometry::sign(&divisor, dec)? == Sign::Zero {
                return Err(unsupported(
                    "smooth field corner normals must be independent",
                ));
            }
            let centre = at.plus(&bisector.scaled(&(-&first.radius / divisor)));
            for plane in planes {
                if geometry::sign(
                    &(&centre.minus(&plane.origin).dot(&unit(&plane.normal)) + &first.radius),
                    dec,
                )? == Sign::Positive
                {
                    return Err(unsupported(
                        "the smooth corner ball leaves an input halfspace",
                    ));
                }
            }
            let vertex = VertexName::new(
                incident
                    .iter()
                    .map(|i| self.faces[*i].origin.clone())
                    .collect(),
            )
            .map_err(|e| unsupported(&e.to_string()))?;
            let origin = Origin {
                step: step.clone(),
                role: Role::Corner { vertex },
            };
            let mut axis = normals[0].cross(&normals[2]);
            match geometry::sign(&axis.dot(&normals[1]), dec)? {
                Sign::Negative => {}
                Sign::Positive => axis = axis.scaled(&-Real::one()),
                Sign::Zero => {
                    return Err(unsupported("the open contact sheet needs an outgoing axis"));
                }
            }
            let wall = incident
                .iter()
                .find_map(|i| {
                    let surface = &self.faces[*i].surface;
                    (!supports.iter().any(|s| s.same_definition(surface))).then_some(surface)
                })
                .ok_or_else(|| unsupported("a smooth field junction needs one end wall"))?;
            let Surface::Plane(wall) = wall else {
                unreachable!()
            };
            let mut outward_wall = None;
            for plane in planes {
                if geometry::sign(&plane.normal.cross(&wall.normal).norm_squared(), dec)?
                    == Sign::Zero
                    && geometry::sign(&plane.origin.minus(&wall.origin).dot(&wall.normal), dec)?
                        == Sign::Zero
                {
                    outward_wall = Some(&plane.normal);
                    break;
                }
            }
            let outward_wall = outward_wall
                .ok_or_else(|| unsupported("the end wall has no certified outward halfspace"))?;
            if geometry::sign(&axis.dot(outward_wall), dec)? != Sign::Positive {
                return Err(unsupported(
                    "the end wall must transversely close the outgoing sheet",
                ));
            }
            let continuation = FilletDefinition::new(
                origin.clone(),
                [supports[0].clone(), supports[2].clone()],
                [true; 2],
                Curve::Line(Line {
                    origin: at.clone(),
                    dir: axis.clone(),
                }),
                first.radius.clone(),
                true,
                dec,
            )?;
            for (index, other_normal) in [(*first_index, &normals[2]), (*second_index, &normals[0])]
            {
                let Surface::Cylinder(cylinder) = &definitions[index].0.surface else {
                    return Err(unsupported(
                        "a planar smooth strip must reduce to a cylinder",
                    ));
                };
                let mut outgoing = cylinder.axis.clone();
                match geometry::sign(&outgoing.dot(other_normal), dec)? {
                    Sign::Negative => {}
                    Sign::Positive => outgoing = outgoing.scaled(&-Real::one()),
                    Sign::Zero => {
                        return Err(unsupported(
                            "a selected strip needs a transverse corner cap",
                        ));
                    }
                }
                if caps.is_empty() {
                    caps.resize_with(definitions.len(), Vec::new);
                }
                if caps[index].is_empty() {
                    caps[index].reserve_exact(4);
                }
                caps[index].push(Surface::Plane(Plane {
                    origin: centre.clone(),
                    normal: outgoing.scaled(&-Real::one()),
                    ref_dir: normals[1].clone(),
                }));
            }
            continuations.push((
                continuation,
                Surface::Plane(Plane {
                    origin: centre.clone(),
                    normal: axis.scaled(&-Real::one()),
                    ref_dir: normals[0].clone(),
                }),
            ));
            sectors.push((
                SphereSector::new(centre, first.radius.clone(), normals, true, dec)?,
                origin,
            ));
        }
        if sectors.is_empty() {
            return Ok(None);
        }
        if endpoints.iter().any(|(_, runs)| runs.len() > 2) {
            return Err(unsupported(
                "mixed smooth and higher-contact junctions need complete field normal fans",
            ));
        }
        let mut result = root;
        for (index, (definition, _)) in definitions.iter().enumerate() {
            if caps[index].len() > 1 {
                let [Surface::Plane(first), Surface::Plane(second)] = caps[index].as_slice() else {
                    return Err(unsupported(
                        "a straight strip has more than two corner caps",
                    ));
                };
                if geometry::sign(&first.normal.dot(&second.normal), dec)? != Sign::Negative
                    || geometry::sign(&second.origin.minus(&first.origin).dot(&first.normal), dec)?
                        != Sign::Negative
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "smooth corner caps consume the selected strip",
                    ));
                }
            }
            let node = self.analytic_strip(definition, std::mem::take(&mut caps[index]), dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        for (definition, cap) in continuations {
            let node = self.analytic_strip(&definition, vec![cap], dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        for (sector, origin) in sectors {
            let planes = sector
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
                .collect();
            let node = self.analytic_region(planes, sector.surface, origin, dec)?;
            result = self.boolean(BooleanKind::Subtract, result, node, dec)?;
        }
        Ok(Some(result))
    }

    pub(super) fn analytic_strip(
        &mut self,
        definition: &FilletDefinition,
        mut caps: Vec<Surface>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let (Surface::Cylinder(cylinder), [Surface::Plane(a), Surface::Plane(b)]) =
            (&definition.surface, &definition.supports)
        else {
            return Err(unsupported(
                "smooth field strips need planar cylindrical contacts",
            ));
        };
        if !definition.convex {
            return Err(unsupported("analytic field strips need convex contacts"));
        }
        let [na, nb] = [0, 1].map(|i| {
            let normal = unit(&[a, b][i].normal);
            if definition.same_sense[i] {
                normal
            } else {
                normal.scaled(&-Real::one())
            }
        });
        let dot = na.dot(&nb);
        for normal in [na.minus(&nb.scaled(&dot)), nb.minus(&na.scaled(&dot))] {
            caps.push(Surface::Plane(Plane {
                origin: cylinder.origin.clone(),
                normal: normal.scaled(&-Real::one()),
                ref_dir: cylinder.axis.clone(),
            }));
        }
        self.analytic_region(
            caps,
            definition.surface.clone(),
            definition.origin.clone(),
            dec,
        )
    }
}
