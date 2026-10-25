//! Straight normal images of exactly one recovered finite polygonal patch.

use super::recovery;
use crate::{Face, Plan, geometry};
use catallax_features::draft::{
    DraftDefinition,
    arrangement::{Patch, Skin},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{
    BooleanKind, EdgeName, FaceRef, Number, Operation, Origin, Pose, Role, StepId, ThickenSide,
    Vector3 as RecipeVector,
};
use catallax_shapes::{Plane, Surface};
use std::collections::BTreeMap;

fn selection_error(error: Refusal) -> Refusal {
    match error.code {
        "recipe.missing-reference" => Refusal::new("reference.lost-face", error.detail),
        "recipe.ambiguous-reference" => Refusal::new("reference.ambiguous-face", error.detail),
        _ => error,
    }
}

/// Exact box reduction of a convex arrangement term. This preserves the existing
/// closed-box moment path for rectilinear polygons; other terms retain their
/// complete analytic halfspaces. It is never a bounding-box approximation.
fn box_leaf(
    surfaces: &[Surface],
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<(geometry::Leaf, Vec<Face>)>, Refusal> {
    let mut lower: [Option<Rational>; 3] = std::array::from_fn(|_| None);
    let mut upper = lower.clone();
    for surface in surfaces {
        let Surface::Plane(plane) = surface else {
            return Ok(None);
        };
        let (Some(normal), Some(origin)) = (plane.normal.as_rational(), plane.origin.as_rational())
        else {
            return Ok(None);
        };
        let n = [normal.x, normal.y, normal.z];
        let p = [origin.x, origin.y, origin.z];
        if n.iter().filter(|x| !x.is_zero()).count() != 1 {
            return Ok(None);
        }
        let axis = n
            .iter()
            .position(|x| !x.is_zero())
            .expect("one nonzero normal coordinate");
        let bound = if n[axis].sign() == Sign::Positive {
            &mut upper[axis]
        } else {
            &mut lower[axis]
        };
        *bound = Some(match bound.take() {
            Some(old) if n[axis].sign() == Sign::Positive => old.min(p[axis].clone()),
            Some(old) => old.max(p[axis].clone()),
            None => p[axis].clone(),
        });
    }
    let [Some(x), Some(y), Some(z)] = lower else {
        return Ok(None);
    };
    let [Some(a), Some(b), Some(c)] = upper else {
        return Ok(None);
    };
    let size = [&a - &x, &b - &y, &c - &z];
    if size.iter().any(|q| q.sign() != Sign::Positive) {
        return Err(recovery::bad("a swept convex term has no positive volume"));
    }
    let mut pose = Pose::identity();
    pose.translation = RecipeVector([x, y, z].map(Number::new));
    Ok(Some(geometry::compile(
        &Operation::Block {
            size: size.map(Number::new),
            pose,
        },
        step,
        dec,
    )?))
}

impl Plan {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn boolean_planar_thicken(
        &mut self,
        node: usize,
        refs: &[FaceRef],
        thickness: &Rational,
        side: ThickenSide,
        step: &StepId,
        fallback: Refusal,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let skin = match recovery::skin(self, node, dec) {
            Ok(skin) => skin,
            Err(e) if e.code.starts_with("draft.") => return Err(fallback),
            Err(e) => return Err(e),
        };
        let selected = skin.resolve_faces(refs, dec).map_err(selection_error)?;
        let mut selected = selected.into_iter();
        let (Some(own), None) = (selected.next(), selected.next()) else {
            return Err(fallback);
        };
        let patch = &skin.patches[own];
        let front = matches!(side, ThickenSide::Front);
        let distance = Real::from(if front { thickness.clone() } else { -thickness });
        let Some(Surface::Plane(far)) =
            Surface::Plane(patch.plane.clone()).parallel_support(&distance, dec)?
        else {
            unreachable!("a plane parallel never collapses")
        };
        let shift = far.origin.minus(&patch.plane.origin);
        let mut remap = BTreeMap::new();
        let mut vertices = Vec::new();
        for i in patch.loops.iter().flatten() {
            if let std::collections::btree_map::Entry::Vacant(entry) = remap.entry(*i) {
                entry.insert(vertices.len());
                vertices.push(skin.vertices[*i].clone());
            }
        }
        let count = vertices.len();
        let tops: Vec<_> = vertices.iter().map(|p| p.plus(&shift)).collect();
        vertices.extend(tops);
        let mut incidence = BTreeMap::<_, Vec<usize>>::new();
        for (index, patch) in skin.patches.iter().enumerate() {
            for cycle in &patch.loops {
                for i in 0..cycle.len() {
                    let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                    incidence
                        .entry((a.min(b), a.max(b)))
                        .or_default()
                        .push(index);
                }
            }
        }
        let mut patches = Vec::new();
        for cycle in &patch.loops {
            for i in 0..cycle.len() {
                let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                let uses = &incidence[&(a.min(b), a.max(b))];
                let others: Vec<_> = uses.iter().copied().filter(|f| *f != own).collect();
                let [other] = others.as_slice() else {
                    return Err(recovery::bad(
                        "each polygonal rim needs exactly one distinct physical neighbour",
                    ));
                };
                if uses.len() != 2 {
                    return Err(recovery::bad(
                        "a swept polygon needs a closed manifold source boundary",
                    ));
                }
                let mut origins = Vec::new();
                for of in &patch.origins {
                    for across in &skin.patches[*other].origins {
                        let edge = EdgeName::new(of.clone(), across.clone())
                            .map_err(|e| recovery::bad(&e.to_string()))?;
                        origins.push(Origin {
                            step: step.clone(),
                            role: Role::Rim { edge },
                        });
                    }
                }
                origins.sort();
                origins.dedup();
                if origins.is_empty() {
                    return Err(recovery::bad("a swept rim has no complete source origins"));
                }
                let mut wire = vec![remap[&a], remap[&b], remap[&b] + count, remap[&a] + count];
                if !front {
                    wire.reverse();
                }
                patches.push(Patch {
                    plane: Plane {
                        origin: skin.vertices[a].clone(),
                        normal: skin.vertices[b]
                            .minus(&skin.vertices[a])
                            .cross(&patch.plane.normal),
                        ref_dir: patch.plane.normal.clone(),
                    },
                    loops: vec![wire],
                    origins,
                });
            }
        }
        let loops = |offset: usize, flip: bool| {
            patch
                .loops
                .iter()
                .map(|cycle| {
                    let mut cycle: Vec<_> = cycle.iter().map(|i| remap[i] + offset).collect();
                    if flip {
                        cycle.reverse();
                    }
                    cycle
                })
                .collect()
        };
        let mut base = patch.plane.clone();
        let mut far = far;
        if front {
            base.normal = base.normal.scaled(&Real::from_i64(-1));
        } else {
            far.normal = far.normal.scaled(&Real::from_i64(-1));
        }
        patches.push(Patch {
            plane: base,
            loops: loops(0, front),
            origins: patch.origins.clone(),
        });
        patches.push(Patch {
            plane: far,
            loops: loops(count, !front),
            origins: patch
                .origins
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Offset {
                        of: Box::new(of.clone()),
                    },
                })
                .collect(),
        });
        let swept = Skin { vertices, patches };
        // Identity subdivision only certifies the closed normal-sweep definition;
        // it does not move a support. Publish the exact convex material terms
        // through existing analytic leaves. Only maximal physical boundary
        // supports get face candidates; internal partition walls get no origins.
        let identity = DraftDefinition::new(
            patch.plane.normal.clone(),
            patch.plane.clone(),
            &Rational::zero(),
            dec,
        )?;
        let components = swept.draft(&Default::default(), &identity, dec)?;
        if components.len() != 1 {
            return Err(recovery::bad(
                "one finite face must sweep to one connected solid",
            ));
        }
        let bound = &self.node_radius(node, dec)? + thickness;
        let mut result = 0;
        for component in components {
            for term in &component.terms {
                let surfaces: Vec<_> = term
                    .sides
                    .iter()
                    .map(|(index, same)| {
                        let mut plane = component.planes[*index].clone();
                        if !same {
                            plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                        }
                        Surface::Plane(plane)
                    })
                    .collect();
                let (leaf, native) = match box_leaf(&surfaces, step, dec)? {
                    Some((leaf, faces)) => (leaf, faces),
                    None => (
                        geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
                        Vec::new(),
                    ),
                };
                let mut faces = Vec::new();
                let mut roles = Vec::new();
                for face in &component.faces {
                    if term.sides.contains(&(face.support, face.same_sense)) {
                        for origin in &face.origins {
                            let face = geometry::face(
                                Surface::Plane(component.plane(face)),
                                origin.clone(),
                                dec,
                            )?;
                            let role = native
                                .iter()
                                .find(|n| n.support.same_definition(&face.support))
                                .map(|n| n.origin.role.clone());
                            if !native.is_empty() && role.is_none() {
                                return Err(recovery::bad(
                                    "physical boundary has no exact primitive chart",
                                ));
                            }
                            roles.push(role);
                            faces.push(face);
                        }
                    }
                }
                let first = self.faces.len();
                let term = self.leaf(leaf, faces, dec)?;
                for (offset, role) in roles.into_iter().enumerate() {
                    if let Some(role) = role {
                        self.definition_roles.insert(first + offset, role);
                    }
                }
                result = self.boolean(BooleanKind::Union, result, term, dec)?;
            }
        }
        Ok(result)
    }
}
