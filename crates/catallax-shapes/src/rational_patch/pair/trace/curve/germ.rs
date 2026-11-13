//! Outgoing constant-speed acceleration from both retained surface jets.

use super::*;
use crate::rational_patch::differential as d;

fn missing(detail: &str) -> Refusal {
    Refusal::new("intersect.patch-pair-germ", detail)
}

impl PatchPairCurve {
    /// Exact normal acceleration for a nonzero outgoing tangent. At a join,
    /// only the cell side entered by that tangent supplies the surface jets.
    pub fn normal_acceleration(
        &self,
        point: &Point3,
        outgoing: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        if exact_sign(&d::dot(outgoing, outgoing, dec)?, dec)? != Sign::Positive {
            return Err(missing(
                "pair acceleration needs a nonzero outgoing tangent",
            ));
        }
        let mut matched_end = false;
        let mut result = None;
        for piece in self.pieces.iter() {
            for end in 0..2 {
                if !same_point(&endpoint(piece, end), point, dec)? {
                    continue;
                }
                matched_end = true;
                let at = piece.arc().ends[end].parameters();
                let natural = tangent(piece, at, dec)?;
                let along = exact_sign(&d::dot(&natural, outgoing, dec)?, dec)?;
                if matches!(
                    (end, along),
                    (0, Sign::Negative) | (1, Sign::Positive) | (_, Sign::Zero)
                ) {
                    continue;
                }
                let value = acceleration(piece, at, point, outgoing, dec)?;
                if let Some(previous) = &result {
                    if !d::equal(previous, &value, dec)? {
                        return Err(missing(
                            "the outgoing tangent selects unequal cell-side accelerations",
                        ));
                    }
                } else {
                    result = Some(value);
                }
            }
        }
        if matched_end {
            return result
                .ok_or_else(|| missing("the tangent does not enter a retained cell side"));
        }
        let (index, _, at) = self
            .locate(point, dec)?
            .ok_or_else(|| missing("the acceleration point is outside the pair component"))?;
        acceleration(&self.pieces[index], &at, point, outgoing, dec)
    }
}

fn acceleration(
    piece: &PairTracePiece,
    at: &[Real; 4],
    point: &Point3,
    outgoing: &Vector3,
    dec: &mut Decider,
) -> Result<Vector3, Refusal> {
    let arc = piece.arc();
    let first = d::Jet::new(piece.first(), &[at[0].clone(), at[1].clone()], dec)?;
    let second = d::Jet::new(
        piece
            .definition
            .second
            .span(arc.second_chart, arc.second_span),
        &[at[2].clone(), at[3].clone()],
        dec,
    )?;
    let position = Vector3::new(point.x.clone(), point.y.clone(), point.z.clone());
    if !d::equal(&first.point, &position, dec)? || !d::equal(&second.point, &position, dec)? {
        return Err(missing(
            "the retained cell needs exact point correspondence on both sources",
        ));
    }
    let first_velocity = first
        .velocity(outgoing, dec)?
        .ok_or_else(|| missing("the outgoing tangent is outside the first source tangent plane"))?;
    let second_velocity = second.velocity(outgoing, dec)?.ok_or_else(|| {
        missing("the outgoing tangent is outside the second source tangent plane")
    })?;
    let a = first.normal(dec)?;
    let b = second.normal(dec)?;
    let first_second = first.second_along(&first_velocity, dec)?;
    let second_second = second.second_along(&second_velocity, dec)?;
    let ka = d::dot(&a, &first_second, dec)?;
    let kb = d::dot(&b, &second_second, dec)?;
    let aa = d::dot(&a, &a, dec)?;
    let ab = d::dot(&a, &b, dec)?;
    let bb = d::dot(&b, &b, dec)?;
    let determinant = d::sub(&d::mul(&aa, &bb, dec)?, &d::mul(&ab, &ab, dec)?, dec)?;
    if exact_sign(&determinant, dec)? != Sign::Positive {
        return Err(missing(
            "pair acceleration needs independent source normals",
        ));
    }
    // A lies in the normal plane and satisfies n_i.A = n_i.S_i''.
    let alpha = d::div(
        &d::sub(&d::mul(&ka, &bb, dec)?, &d::mul(&kb, &ab, dec)?, dec)?,
        &determinant,
        dec,
    )?;
    let beta = d::div(
        &d::sub(&d::mul(&kb, &aa, dec)?, &d::mul(&ka, &ab, dec)?, dec)?,
        &determinant,
        dec,
    )?;
    let result = d::plus(&d::scale(&a, &alpha, dec)?, &d::scale(&b, &beta, dec)?, dec)?;
    for residual in [
        d::sub(&d::dot(&a, &result, dec)?, &ka, dec)?,
        d::sub(&d::dot(&b, &result, dec)?, &kb, dec)?,
        d::dot(outgoing, &result, dec)?,
    ] {
        if exact_sign(&residual, dec)? != Sign::Zero {
            return Err(missing(
                "the paired-source constant-speed acceleration equations did not prove",
            ));
        }
    }
    Ok(result)
}
