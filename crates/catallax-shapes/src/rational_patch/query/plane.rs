//! A nonzero plane pullback proves different supports, not disjoint supports.

use super::affine::arithmetic;
use super::*;
use std::collections::{BTreeMap, btree_map::Entry};

impl RationalPatch {
    pub(crate) fn proves_off_plane(
        &self,
        plane: &crate::Plane,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        math::validate_point(&plane.origin, dec)?;
        math::validate_vector(&plane.normal, dec)?;
        let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
        let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
        let mut offset = Real::zero();
        for axis in 0..3 {
            let term = arithmetic(normal[axis], origin[axis], true, dec)?;
            offset = arithmetic(&offset, &term, false, dec)?;
        }
        let factors = [
            normal[0].clone(),
            normal[1].clone(),
            normal[2].clone(),
            -offset,
        ];
        for span in &self.data.spans {
            let mut coefficients = BTreeMap::<Vec<u32>, Real>::new();
            for (coordinate, factor) in span.coordinates.iter().zip(&factors) {
                for (monomial, coefficient) in coordinate.terms() {
                    let term = arithmetic(&coefficient, factor, true, dec)?;
                    match coefficients.entry(monomial.exponents().to_vec()) {
                        Entry::Vacant(entry) => {
                            entry.insert(term);
                        }
                        Entry::Occupied(mut entry) => {
                            let sum = arithmetic(entry.get(), &term, false, dec)?;
                            entry.insert(sum);
                        }
                    }
                }
            }
            // W is nonzero, and each admitted span has a rectangle with
            // interior. A nonzero numerator polynomial cannot vanish on it.
            if coefficients
                .values()
                .any(|value| value.rational_ref().is_some_and(|q| !q.is_zero()))
            {
                return Ok(true);
            }
            for coefficient in coefficients.values() {
                if exact_sign(coefficient, dec)? != Sign::Zero {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}
