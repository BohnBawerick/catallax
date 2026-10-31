//! Corner membership reuses strict interval states without computing distances.
use super::{BoxState, Classifier};
use crate::{Instruction, Location, Workspace, combine, polynomial_sign};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::BooleanKind;

impl Classifier<'_> {
    pub(super) fn corner_location(
        &self,
        point: &[Rational; 3],
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let plan = self.plan;
        storage.locations.clear();
        for (slot, instruction) in plan.instructions.iter().enumerate() {
            let location = match *instruction {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(i) => match self.states[slot] {
                    // These states certify the same exact point, not a rounded sample.
                    Some(BoxState::Inside) => Location::Inside,
                    Some(BoxState::Outside) => Location::Outside,
                    Some(BoxState::Unknown) | None => plan.leaves[i].location(point, dec)?,
                },
                Instruction::Boolean { kind, a, b } => {
                    combine(kind, storage.locations[a], storage.locations[b])
                }
            };
            storage.locations.push(location);
        }
        let mut location = storage.locations[plan.root];
        if location != Location::Indeterminate {
            return Ok(location);
        }
        if self.on_intersected_boundary(&storage.locations) {
            return Ok(Location::Boundary);
        }
        // Regularised Boolean contacts need the same active faces and exact
        // planar/transverse witnesses as Plan::proven_checked. A leaf boundary
        // alone cannot decide whether a Boolean seam is inside or outside.
        storage.active.clear();
        for (id, face) in plan.faces.iter().enumerate() {
            if face.in_result
                && plan.leaf_registers[face.leaf]
                    .is_some_and(|slot| storage.locations[slot] == Location::Boundary)
                && polynomial_sign(&face.polynomial, point, dec)? == Sign::Zero
                && plan.leaves[face.leaf].on_face(
                    &face.origin,
                    plan.definition_role(id),
                    point,
                    dec,
                )? != Some(false)
            {
                storage.active.push(id);
            }
        }
        if let Some(proven) = plan.planar_location(point, storage, dec)? {
            location = proven;
        } else if plan.boundary_witness(point, storage, dec)? {
            location = Location::Boundary;
        }
        Ok(location)
    }

    /// A point on the boundary of a leaf that an intersection requires is not in the
    /// intersection's interior, which is all a grid corner needs to know. Such a
    /// corner can have no fixed witness ray, as where a gyroid meets a box corner.
    /// Checking this before the contact queries skips the costliest corner work.
    fn on_intersected_boundary(&self, locations: &[Location]) -> bool {
        let mut pending = vec![self.plan.root];
        while let Some(node) = pending.pop() {
            match self.plan.instructions[node] {
                Instruction::Leaf(_) if locations[node] == Location::Boundary => return true,
                Instruction::Boolean {
                    kind: BooleanKind::Intersect,
                    a,
                    b,
                } => pending.extend([a, b]),
                _ => (),
            }
        }
        false
    }
}
