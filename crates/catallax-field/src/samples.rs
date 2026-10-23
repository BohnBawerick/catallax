//! Query-only sampling support. These bounds and candidates never build geometry.
use crate::{Instruction, Plan, geometry::Leaf};
use catallax_numbers::{Rational, Refusal};
use catallax_recipe::{BooleanKind, Origin};
use catallax_shapes::{Isometry, RPoint3};

pub(crate) type Bounds = [[Rational; 2]; 3];

pub(crate) fn transform(bounds: &Bounds, pose: &Isometry) -> Result<Bounds, Refusal> {
    let mut result: Option<Bounds> = None;
    for corner in 0..8 {
        let [x, y, z] = std::array::from_fn(|i| bounds[i][(corner >> i) & 1].clone());
        let p = pose
            .apply_rational_point(&RPoint3::new(x, y, z))
            .ok_or_else(placement_error)?;
        let coordinates = [p.x, p.y, p.z];
        if let Some(b) = &mut result {
            for i in 0..3 {
                b[i][0] = b[i][0].clone().min(coordinates[i].clone());
                b[i][1] = b[i][1].clone().max(coordinates[i].clone());
            }
        } else {
            result = Some(coordinates.map(|x| [x.clone(), x]));
        }
    }
    Ok(result.unwrap())
}
pub(crate) fn placement_error() -> Refusal {
    Refusal::new(
        "field.sample-algebraic-placement",
        "the rational candidate API needs a native rational placement",
    )
}
pub(crate) fn fraction(i: u32, n: u32) -> Rational {
    Rational::from_ratio(i64::from(2 * i + 1), i64::from(2 * n)).unwrap()
}
pub(crate) fn circle(n: u32) -> Vec<[Rational; 2]> {
    let mut points = Vec::with_capacity((4 * n) as usize);
    for quarter in 0..4 {
        for i in 0..n {
            let t = fraction(i, n);
            let squared = &t * &t;
            let denominator = &Rational::one() + &squared;
            let x = &(Rational::one() - squared) / &denominator;
            let y = Rational::from_i64(2) * t / denominator;
            points.push(match quarter {
                0 => [x, y],
                1 => [-y, x],
                2 => [-x, -y],
                _ => [y, -x],
            });
        }
    }
    points
}
impl Plan {
    /// Visits rational supporting-face candidates for sketch-based features.
    /// Stock-solid candidates are not included. Candidates are independent of
    /// B-rep construction and display meshes; use `face_query` to establish
    /// surviving boundary membership, names and orientation before accepting one.
    /// `divisions` is in 1..=64. Returning false from `visit` stops cooperatively
    /// and returns `Ok(false)`, never a claim that all candidates were visited.
    /// Through-all walls are sampled over a conservative bound on the field
    /// result, not an arbitrary finite replacement for the unbounded tool.
    /// An algebraic placement without a native rational frame returns
    /// `field.sample-algebraic-placement`; no candidate coordinates are rounded.
    pub fn visit_feature_candidates(
        &self,
        divisions: u32,
        mut visit: impl FnMut(RPoint3, &Origin) -> bool,
    ) -> Result<bool, Refusal> {
        if !(1..=64).contains(&divisions) {
            return Err(Refusal::new(
                "field.sample-options",
                "feature sampling divisions must be in 1..=64",
            ));
        }
        if !self.faces.iter().any(|face| {
            face.in_result
                && self.leaf_registers[face.leaf].is_some()
                && matches!(self.leaves[face.leaf], Leaf::Feature(_))
        }) {
            return Ok(true);
        }
        if self.is_proven_empty() {
            return Ok(true);
        }
        let mut bounds: Vec<Option<Bounds>> = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            let b = match *instruction {
                Instruction::Empty => Some(std::array::from_fn(|_| {
                    [Rational::zero(), Rational::zero()]
                })),
                Instruction::Leaf(i) => match &self.leaves[i] {
                    Leaf::Primitive(p) => Some(p.sample_bounds()?),
                    Leaf::Feature(f) => f.sample_bounds()?,
                },
                Instruction::Boolean { kind, a, b } => match kind {
                    BooleanKind::Subtract => bounds[a].clone(),
                    BooleanKind::Union => match (&bounds[a], &bounds[b]) {
                        (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                            [
                                a[i][0].clone().min(b[i][0].clone()),
                                a[i][1].clone().max(b[i][1].clone()),
                            ]
                        })),
                        _ => None,
                    },
                    BooleanKind::Intersect => match (&bounds[a], &bounds[b]) {
                        (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                            [
                                a[i][0].clone().max(b[i][0].clone()),
                                a[i][1].clone().min(b[i][1].clone()),
                            ]
                        })),
                        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
                        (None, None) => None,
                    },
                },
            };
            bounds.push(b);
        }
        let result = bounds[self.root].as_ref().ok_or_else(|| {
            Refusal::new(
                "field.sample-unbounded",
                "result has no finite feature sampling bound",
            )
        })?;
        for (index, face) in self.faces.iter().enumerate() {
            if !face.in_result || self.leaf_registers[face.leaf].is_none() {
                continue;
            }
            if let Leaf::Feature(feature) = &self.leaves[face.leaf]
                && !feature.sample_candidates(
                    self.definition_role(index),
                    divisions,
                    result,
                    &mut |p| visit(p, &face.origin),
                )?
            {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
