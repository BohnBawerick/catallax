//! The structural checks: everything about a recipe that needs no geometry.

use std::collections::{BTreeMap, BTreeSet};

use crate::ids::{ElementId, StepId};
use crate::names::{
    EdgeName, EdgeRef, EdgeSelection, FaceRef, HolePart, Origin, OriginPattern, Role, VertexName,
};
use crate::operation::{Grade, HoleDepth, HoleHead, HoleTip, Operation, ValueKind};
use crate::recipe::{FORMAT, FORMAT_V4, Recipe, RecipeError};

pub(crate) fn check(recipe: &Recipe) -> Result<(), RecipeError> {
    if recipe.format != FORMAT && recipe.format != FORMAT_V4 {
        return Err(RecipeError::Format {
            found: recipe.format.clone(),
        });
    }
    let mut index = BTreeMap::new();
    for (i, step) in recipe.steps.iter().enumerate() {
        if recipe.format == FORMAT
            && matches!(
                &step.op,
                Operation::ExtrudeThrough { side: Some(_), .. }
                    | Operation::ExtrudeTo { .. }
                    | Operation::ExtrudeNext { .. }
            )
        {
            return Err(RecipeError::Parse {
                message: "extrusion end conditions require catallax-recipe/4".into(),
            });
        }
        if index.contains_key(&step.id) {
            return Err(RecipeError::DuplicateStep {
                step: step.id.clone(),
            });
        }
        let ctx = Ctx {
            recipe,
            index: &index,
            step: &step.id,
        };
        ctx.inputs(&step.op)?;
        ctx.elements(&step.op)?;
        ctx.rules(&step.op)?;
        for pattern in patterns(&step.op) {
            ctx.pattern(pattern)?;
        }
        index.insert(&step.id, i);
    }
    match index.get(&recipe.result) {
        Some(&i) if recipe.steps[i].op.kind() == ValueKind::Body => Ok(()),
        _ => Err(RecipeError::Result {
            result: recipe.result.clone(),
        }),
    }
}

/// Every origin pattern a step's references use.
fn patterns(op: &Operation) -> Vec<&OriginPattern> {
    fn face(r: &FaceRef) -> impl Iterator<Item = &OriginPattern> {
        std::iter::once(&r.face).chain(&r.touching)
    }
    fn edge(r: &EdgeRef) -> impl Iterator<Item = &OriginPattern> {
        r.between.iter().chain(&r.touching)
    }
    fn grade(g: &Grade) -> Vec<&OriginPattern> {
        match g {
            Grade::Distance { faces, .. } => faces.iter().flat_map(face).collect(),
            Grade::Constant(_) | Grade::Linear { .. } => vec![],
        }
    }
    match op {
        Operation::ExtrudeTo { limit, .. } => face(limit).collect(),
        Operation::Fillet {
            edges: EdgeSelection::Refs(refs),
            ..
        }
        | Operation::Chamfer {
            edges: EdgeSelection::Refs(refs),
            ..
        } => refs.iter().flat_map(edge).collect(),
        Operation::Shell { open: faces, .. }
        | Operation::Thicken { faces, .. }
        | Operation::Draft { faces, .. } => faces.iter().flat_map(face).collect(),
        Operation::Lattice { radius: g, .. } | Operation::Gyroid { level: g, .. } => grade(g),
        _ => vec![],
    }
}

/// True when the operation makes faces of its own or carries faces of its inputs into
/// faces it names. A boolean and a draft keep their inputs' names and make none.
fn makes_faces(op: &Operation) -> bool {
    !matches!(
        op,
        Operation::Sketch { .. }
            | Operation::Path { .. }
            | Operation::Boolean { .. }
            | Operation::Draft { .. }
    )
}

struct Ctx<'r> {
    recipe: &'r Recipe,
    /// The steps before the one being checked, by id.
    index: &'r BTreeMap<&'r StepId, usize>,
    /// The step being checked.
    step: &'r StepId,
}

impl Ctx<'_> {
    fn inputs(&self, op: &Operation) -> Result<(), RecipeError> {
        for (input, expected) in op.inputs() {
            let Some(&i) = self.index.get(input) else {
                return Err(RecipeError::UnknownInput {
                    step: self.step.clone(),
                    input: input.clone(),
                });
            };
            let found = self.recipe.steps[i].op.kind();
            if found != expected {
                return Err(RecipeError::WrongInput {
                    step: self.step.clone(),
                    input: input.clone(),
                    expected,
                    found,
                });
            }
        }
        Ok(())
    }

    fn elements(&self, op: &Operation) -> Result<(), RecipeError> {
        let mut seen = BTreeSet::new();
        for element in op.elements() {
            if !seen.insert(element) {
                return Err(RecipeError::DuplicateElement {
                    step: self.step.clone(),
                    element: element.clone(),
                });
            }
        }
        Ok(())
    }

    /// Rules particular to one operation.
    fn rules(&self, op: &Operation) -> Result<(), RecipeError> {
        match op {
            Operation::Loft { sections, .. } => {
                if sections.len() < 2 {
                    return Err(RecipeError::LoftSections {
                        step: self.step.clone(),
                    });
                }
                for section in sections {
                    let Operation::Sketch { loops, .. } = self.op(&section.profile) else {
                        unreachable!("inputs() checked the section is a sketch");
                    };
                    let mut used = BTreeSet::new();
                    for start in &section.start {
                        self.element(&section.profile, start)?;
                        let lp = loops
                            .iter()
                            .position(|l| l.0.iter().any(|s| s.id() == start))
                            .expect("element() found the segment");
                        used.insert(lp);
                    }
                    if used.len() != section.start.len() || used.len() != loops.len() {
                        return Err(RecipeError::LoftStarts {
                            step: self.step.clone(),
                            profile: section.profile.clone(),
                        });
                    }
                }
                Ok(())
            }
            Operation::Chamfer {
                edges: EdgeSelection::All,
                other_distance: Some(_),
                ..
            } => Err(RecipeError::AsymmetricAll {
                step: self.step.clone(),
            }),
            Operation::ImportStep { file, .. } | Operation::ImportMesh { file, .. } => {
                let hex = |c: u8| c.is_ascii_digit() || (b'a'..=b'f').contains(&c);
                if file.sha256.len() == 64 && file.sha256.bytes().all(hex) {
                    Ok(())
                } else {
                    Err(RecipeError::FileHash {
                        step: self.step.clone(),
                    })
                }
            }
            _ => Ok(()),
        }
    }

    /// The operation of an earlier step that `inputs()` or `origin()` already found.
    fn op(&self, id: &StepId) -> &Operation {
        &self.recipe.steps[self.index[id]].op
    }

    /// Checks that the earlier step `of` has the element `element`.
    fn element(&self, of: &StepId, element: &ElementId) -> Result<(), RecipeError> {
        if self.op(of).elements().contains(&element) {
            Ok(())
        } else {
            Err(RecipeError::UnknownElement {
                step: self.step.clone(),
                of: of.clone(),
                element: element.clone(),
            })
        }
    }

    /// Finds the step a name mentions, which must come before step number `before`.
    fn earlier(&self, origin: &StepId, before: usize) -> Result<usize, RecipeError> {
        match self.index.get(origin) {
            Some(&i) if i < before => Ok(i),
            _ => Err(RecipeError::UnknownOrigin {
                step: self.step.clone(),
                origin: origin.clone(),
            }),
        }
    }

    fn pattern(&self, pattern: &OriginPattern) -> Result<(), RecipeError> {
        let i = self.earlier(&pattern.step, usize::MAX)?;
        if !makes_faces(&self.recipe.steps[i].op) {
            return Err(RecipeError::NoFaces {
                step: self.step.clone(),
                origin: pattern.step.clone(),
            });
        }
        match &pattern.role {
            Some(role) => self.role(i, role),
            None => Ok(()),
        }
    }

    /// Checks an origin made by a step before step number `before`.
    fn origin(&self, origin: &Origin, before: usize) -> Result<(), RecipeError> {
        let i = self.earlier(&origin.step, before)?;
        self.role(i, &origin.role)
    }

    fn edge(&self, edge: &EdgeName, before: usize) -> Result<(), RecipeError> {
        edge.between()
            .iter()
            .try_for_each(|o| self.origin(o, before))
    }

    fn vertex(&self, vertex: &VertexName, before: usize) -> Result<(), RecipeError> {
        vertex
            .faces()
            .iter()
            .try_for_each(|o| self.origin(o, before))
    }

    /// Checks that step number `i` makes faces with `role`. Names inside the role must
    /// come from steps before `i`, since they name faces of its input.
    fn role(&self, i: usize, role: &Role) -> Result<(), RecipeError> {
        let step = &self.recipe.steps[i];
        let wrong = |name: String| RecipeError::WrongRole {
            step: self.step.clone(),
            origin: step.id.clone(),
            role: name,
        };
        let no_copy = |instance: u32| RecipeError::NoSuchCopy {
            step: self.step.clone(),
            origin: step.id.clone(),
            instance,
        };
        use Operation as O;
        match (&step.op, role) {
            (O::Block { .. } | O::Cylinder { .. }, Role::Start | Role::End) => Ok(()),
            (O::Cone { base_radius, .. }, Role::Start)
                if base_radius.value().sign() == catallax_numbers::Sign::Positive =>
            {
                Ok(())
            }
            (O::Cone { top_radius, .. }, Role::End)
                if top_radius.value().sign() == catallax_numbers::Sign::Positive =>
            {
                Ok(())
            }
            (
                O::Block { .. },
                Role::Side {
                    segment,
                    path: None,
                },
            ) if matches!(segment.as_str(), "x_min" | "x_max" | "y_min" | "y_max") => Ok(()),
            (
                O::Cylinder { .. } | O::Sphere { .. } | O::Cone { .. } | O::Torus { .. },
                Role::Side {
                    segment,
                    path: None,
                },
            ) if segment.as_str() == "wall" => Ok(()),
            (
                O::Extrude { .. }
                | O::ExtrudeTo { .. }
                | O::ExtrudeNext { .. }
                | O::Revolve { .. }
                | O::Sweep { .. }
                | O::Loft { .. },
                Role::Start | Role::End,
            ) => Ok(()),
            (
                O::ExtrudeThrough {
                    side:
                        Some(
                            crate::operation::ThroughSide::Forward
                            | crate::operation::ThroughSide::Backward,
                        ),
                    ..
                },
                Role::Start,
            ) => Ok(()),
            (
                O::Extrude { profile, .. }
                | O::ExtrudeThrough { profile, .. }
                | O::ExtrudeTo { profile, .. }
                | O::ExtrudeNext { profile, .. }
                | O::Revolve { profile, .. },
                Role::Side {
                    segment,
                    path: None,
                },
            ) => self.element(profile, segment),
            (
                O::Sweep { profile, path, .. },
                Role::Side {
                    segment,
                    path: Some(along),
                },
            ) => {
                self.element(profile, segment)?;
                self.element(path, along)
            }
            (
                O::Loft { sections, .. },
                Role::Side {
                    segment,
                    path: None,
                },
            ) => self.element(&sections[0].profile, segment),
            (
                O::Hole {
                    at,
                    depth,
                    tip,
                    head,
                    ..
                },
                Role::Hole { at: id, part },
            ) => {
                if !at.iter().any(|p| &p.id == id) {
                    return Err(RecipeError::UnknownElement {
                        step: self.step.clone(),
                        of: step.id.clone(),
                        element: id.clone(),
                    });
                }
                let blind = matches!(depth, HoleDepth::Blind(_));
                let made = match part {
                    HolePart::Entry | HolePart::Wall => true,
                    HolePart::Floor => blind && (*tip == HoleTip::Flat || head.is_some()),
                    HolePart::Tip => blind && matches!(tip, HoleTip::Cone(_)),
                    HolePart::CounterboreWall | HolePart::CounterboreFloor => {
                        matches!(head, Some(HoleHead::Counterbore { .. }))
                    }
                    HolePart::Countersink => matches!(head, Some(HoleHead::Countersink { .. })),
                };
                if made {
                    Ok(())
                } else {
                    Err(wrong(format!("hole.{}", part.name())))
                }
            }
            (O::Fillet { .. }, Role::Blend { edge })
            | (O::Chamfer { .. }, Role::Bevel { edge })
            | (O::Shell { .. } | O::Offset { .. }, Role::Round { edge })
            | (O::Thicken { .. }, Role::Rim { edge }) => self.edge(edge, i),
            (O::Fillet { .. } | O::Chamfer { .. }, Role::Corner { vertex })
            | (O::Shell { .. } | O::Offset { .. }, Role::RoundCorner { vertex }) => {
                self.vertex(vertex, i)
            }
            (O::Shell { .. }, Role::Inner { of }) | (O::Thicken { .. }, Role::Offset { of }) => {
                self.origin(of, i)
            }
            (O::Pattern { layout, .. }, Role::Copy { instance, of }) => {
                if *instance == 0 || *instance >= layout.count() {
                    return Err(no_copy(*instance));
                }
                self.origin(of, i)
            }
            (O::Mirror { .. }, Role::Copy { instance, of }) => {
                if *instance != 1 {
                    return Err(no_copy(*instance));
                }
                self.origin(of, i)
            }
            (O::Lattice { .. }, Role::Beams)
            | (O::Gyroid { .. }, Role::Surface)
            | (O::SmoothBlend { .. }, Role::Smooth)
            | (O::ImportStep { .. }, Role::Imported { .. })
            | (O::ImportMesh { .. }, Role::Mesh) => Ok(()),
            _ => Err(wrong(role.name().into())),
        }
    }
}
