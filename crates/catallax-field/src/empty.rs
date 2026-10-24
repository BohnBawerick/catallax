//! Query-only whole-result exclusion using native primitive bounds.

use crate::{Instruction, Plan, geometry::Leaf};
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_recipe::BooleanKind;

type Box3 = [[Rational; 2]; 3];

#[derive(Clone)]
enum Cover {
    Empty,
    Finite(Box3),
    Unbounded,
}
impl Cover {
    fn boolean(kind: BooleanKind, a: &Self, b: &Self) -> Self {
        use Cover::*;
        match (kind, a, b) {
            (BooleanKind::Subtract, _, _) => a.clone(),
            (BooleanKind::Union, Empty, _) => b.clone(),
            (BooleanKind::Union, _, Empty) => a.clone(),
            (BooleanKind::Union, Finite(a), Finite(b)) => Finite(std::array::from_fn(|i| {
                [
                    a[i][0].clone().min(b[i][0].clone()),
                    a[i][1].clone().max(b[i][1].clone()),
                ]
            })),
            (BooleanKind::Union, _, _) => Unbounded,
            (BooleanKind::Intersect, Empty, _) | (BooleanKind::Intersect, _, Empty) => Empty,
            (BooleanKind::Intersect, Unbounded, _) => b.clone(),
            (BooleanKind::Intersect, _, Unbounded) => a.clone(),
            (BooleanKind::Intersect, Finite(a), Finite(b)) => {
                let bounds: Box3 = std::array::from_fn(|i| {
                    [
                        a[i][0].clone().max(b[i][0].clone()),
                        a[i][1].clone().min(b[i][1].clone()),
                    ]
                });
                if bounds.iter().any(|[lo, hi]| lo > hi) {
                    Empty
                } else {
                    Finite(bounds)
                }
            }
        }
    }
}

fn boolean(kind: BooleanKind, a: Option<bool>, b: Option<bool>) -> Option<bool> {
    match kind {
        BooleanKind::Union => match (a, b) {
            (Some(true), _) | (_, Some(true)) => Some(true),
            (Some(false), _) => b,
            (_, Some(false)) => a,
            _ => None,
        },
        BooleanKind::Intersect => match (a, b) {
            (Some(false), _) | (_, Some(false)) => Some(false),
            (Some(true), _) => b,
            (_, Some(true)) => a,
            _ => None,
        },
        BooleanKind::Subtract => boolean(BooleanKind::Intersect, a, b.map(|v| !v)),
    }
}

impl Plan {
    /// Proves an empty regularised result by excluding its entire native box.
    /// Only primitive constraint conjunctions classify boxes. Other definitions
    /// remain unknown, and `false` is not a non-emptiness certificate. This query
    /// neither builds a B-rep nor changes the compiled definition.
    pub fn prove_empty(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        if self.is_proven_empty() {
            return Ok(true);
        }
        let mut covers = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            dec.check_deadline()?;
            covers.push(match *instruction {
                Instruction::Empty => Cover::Empty,
                Instruction::Leaf(i) => self.leaves[i]
                    .display_box(dec)?
                    .map_or(Cover::Unbounded, |b| {
                        Cover::Finite(b.map(|v| [v.lo(), v.hi()]))
                    }),
                Instruction::Boolean { kind, a, b } => Cover::boolean(kind, &covers[a], &covers[b]),
            });
        }
        let bounds = match &covers[self.root] {
            Cover::Empty => return Ok(true),
            Cover::Unbounded => return Ok(false),
            Cover::Finite(bounds) => bounds
                .each_ref()
                .map(|[lo, hi]| Ball::from_rational(lo, 128).hull(&Ball::from_rational(hi, 128))),
        };
        let mut states = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            dec.check_deadline()?;
            let state = match *instruction {
                Instruction::Empty => Some(false),
                Instruction::Leaf(i) => {
                    if !matches!(self.leaves[i], Leaf::Primitive(_)) {
                        states.push(None);
                        continue;
                    }
                    let mut inside = true;
                    let mut outside = false;
                    for constraint in self.leaves[i].constraints() {
                        let value = constraint.polynomial.eval_ball(&bounds, dec)?;
                        if value.lo() > Rational::zero() {
                            outside = true;
                            break;
                        }
                        inside &= value.hi() < Rational::zero();
                    }
                    if outside {
                        Some(false)
                    } else if inside {
                        Some(true)
                    } else {
                        None
                    }
                }
                Instruction::Boolean { kind, a, b } => boolean(kind, states[a], states[b]),
            };
            states.push(state);
        }
        Ok(states[self.root] == Some(false))
    }
}
