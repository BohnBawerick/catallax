//! One polygonal planar face swept from independent field definitions.

mod recovery;
mod sweep;

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{FaceRef, Origin, Role, StepId, ThickenSide};

use crate::{Instruction, Plan, geometry};

impl Plan {
    /// Sweeps the one start or end face a selection resolves to along its normal.
    /// Any source or selection outside that family reports `fallback`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn planar_cap_thicken(
        &mut self,
        node: usize,
        refs: &[FaceRef],
        t: &Rational,
        side: ThickenSide,
        step: &StepId,
        fallback: Refusal,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let Instruction::Leaf(leaf) = self.instructions[node] else {
            return self.boolean_planar_thicken(node, refs, t, side, step, fallback, dec);
        };
        let geometry::Leaf::Feature(feature) = &self.leaves[leaf] else {
            return Err(fallback);
        };
        let members: Vec<_> = (0..self.faces.len())
            .filter(|i| self.faces[*i].leaf == leaf)
            .collect();
        let mut cap = None;
        for reference in refs {
            let found: Vec<_> = members
                .iter()
                .copied()
                .filter(|i| reference.face.matches(&self.faces[*i].origin))
                .filter(|i| {
                    reference.touching.iter().all(|pattern| {
                        members.iter().any(|j| {
                            j != i
                                && self.definition_role(*j) != self.definition_role(*i)
                                && pattern.matches(&self.faces[*j].origin)
                        })
                    })
                })
                .collect();
            let Some(first) = found.first() else {
                return Err(Refusal::new(
                    "reference.lost-face",
                    format!("no input face resolves {reference:?}"),
                ));
            };
            let role = self.definition_role(*first).clone();
            if found.iter().any(|i| *self.definition_role(*i) != role)
                || !matches!(role, Role::Start | Role::End)
                || cap.as_ref().is_some_and(|c| *c != role)
            {
                return Err(fallback);
            }
            cap = Some(role);
        }
        let Some(cap) = cap else {
            return Err(fallback);
        };
        let named = |role: &Role| -> Vec<Origin> {
            let mut names: Vec<_> = members
                .iter()
                .filter(|i| self.definition_role(**i) == role)
                .map(|i| self.faces[*i].origin.clone())
                .collect();
            names.sort();
            names.dedup();
            names
        };
        let neighbour = |ids: &[catallax_recipe::ElementId]| -> Option<Origin> {
            ids.iter()
                .flat_map(|id| {
                    named(&Role::Side {
                        segment: id.clone(),
                        path: None,
                    })
                })
                .min()
        };
        let Some(sweep) = feature.cap_sweep(
            cap == Role::End,
            t,
            matches!(side, ThickenSide::Front),
            &named(&cap),
            &neighbour,
            step,
            dec,
        )?
        else {
            return Err(fallback);
        };
        let first = self.faces.len();
        let node = self.leaf(geometry::Leaf::Feature(sweep.feature), sweep.faces, dec)?;
        for (offset, role) in sweep.roles.into_iter().enumerate() {
            self.definition_roles.insert(first + offset, role);
        }
        Ok(node)
    }
}
