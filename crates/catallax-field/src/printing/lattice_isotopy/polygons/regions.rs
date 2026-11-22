use super::{Box3, Budget, Model, Support, overlaps};
use crate::Instruction;
use catallax_numbers::{Decider, Refusal};
use catallax_recipe::BooleanKind;

fn intersections(
    a: &[Box3],
    b: Option<&[Box3]>,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Vec<Box3>, Refusal> {
    let Some(b) = b else {
        budget.visit_many(a.len(), dec)?;
        return Ok(a.to_vec());
    };
    let mut result = Vec::new();
    for a in a {
        for b in b {
            budget.visit(dec)?;
            if overlaps(a, b) {
                result.push(std::array::from_fn(|i| {
                    [
                        a[i][0].clone().max(b[i][0].clone()),
                        a[i][1].clone().min(b[i][1].clone()),
                    ]
                }));
            }
        }
    }
    result.sort_unstable();
    result.dedup();
    Ok(result)
}

// A bounded positive operand limits where a Boolean can depend on its other
// operand. Preserve separate boxes across unions to keep short lattice limbs
// from introducing cuts on unrelated faces.
pub(super) fn leaf_regions(
    model: &Model<'_>,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Vec<Vec<Box3>>, Refusal> {
    let mut finite: Vec<Option<Vec<Box3>>> = Vec::with_capacity(model.plan.instructions.len());
    for instruction in &model.plan.instructions {
        let bound = match *instruction {
            Instruction::Empty => Some(vec![]),
            Instruction::Leaf(i) => match model.supports[i].as_ref().expect("live support") {
                Support::Block(b) => Some(vec![b.clone()]),
                Support::Cylinder {
                    centre,
                    radius,
                    height,
                } => Some(vec![[
                    [&centre[0] - radius, &centre[0] + radius],
                    [&centre[1] - radius, &centre[1] + radius],
                    [centre[2].clone(), &centre[2] + height],
                ]]),
                Support::Cubic { .. } => None,
            },
            Instruction::Boolean { kind, a, b } => match kind {
                BooleanKind::Subtract => finite[a]
                    .as_deref()
                    .map(|a| intersections(a, None, budget, dec))
                    .transpose()?,
                BooleanKind::Intersect => match (&finite[a], &finite[b]) {
                    (Some(a), b) => Some(intersections(a, b.as_deref(), budget, dec)?),
                    (None, b) => b
                        .as_deref()
                        .map(|b| intersections(b, None, budget, dec))
                        .transpose()?,
                },
                BooleanKind::Union => match (&finite[a], &finite[b]) {
                    (Some(a), Some(b)) => {
                        budget.visit_many(a.len() + b.len(), dec)?;
                        let mut boxes: Vec<_> = a.iter().chain(b).cloned().collect();
                        boxes.sort_unstable();
                        boxes.dedup();
                        Some(boxes)
                    }
                    _ => None,
                },
            },
        };
        finite.push(bound);
    }
    let mut contexts = vec![vec![]; model.plan.instructions.len()];
    let mut leaves = vec![vec![]; model.supports.len()];
    contexts[model.plan.root].push(model.domain.clone());
    for (index, instruction) in model.plan.instructions.iter().enumerate().rev() {
        let mut boxes = std::mem::take(&mut contexts[index]);
        boxes.sort_unstable();
        boxes.dedup();
        match *instruction {
            Instruction::Empty => {}
            Instruction::Leaf(i) => leaves[i].append(&mut boxes),
            Instruction::Boolean { kind, a, b } => {
                let a_mask = if kind == BooleanKind::Intersect {
                    finite[b].as_deref()
                } else {
                    None
                };
                let b_mask = if kind == BooleanKind::Union {
                    None
                } else {
                    finite[a].as_deref()
                };
                contexts[a].extend(intersections(&boxes, a_mask, budget, dec)?);
                contexts[b].extend(intersections(&boxes, b_mask, budget, dec)?);
            }
        }
    }
    Ok(leaves)
}
