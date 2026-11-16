//! Actual finite pcurve joins, checked in physical space.
use super::{Index, point, uv};
use crate::{
    Sense,
    numerical::{NumericalEdge, NumericalFace, NumericalTrim, bad, distance},
};
use catallax_numbers::{Ball, Decider, NumericalCheck, NumericalSample, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, Surface};

fn endpoint(
    face: &NumericalFace,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    curve: Option<&BSplineCurve>,
    end: bool,
    dec: &mut Decider,
) -> Result<([Ball; 3], Rational), Refusal> {
    dec.charge_work(1)?;
    let slot = usize::from(end);
    let t = &edge.parameters[if trim.sense == Sense::Same {
        slot
    } else {
        1 - slot
    }];
    let at = if let Some(curve) = curve {
        uv(curve, t, dec)?
    } else if matches!(face.surface, Surface::Plane(_)) {
        let p = crate::numerical::curve_at(&edge.curve, t, dec)?
            .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
        let uv = face
            .surface
            .parameters_of(&p, dec)?
            .ok_or_else(|| bad("trim", "planar join endpoint has no parameter projection"))?;
        return Ok((face.surface.evaluate(&uv[0], &uv[1], dec)?, t.clone()));
    } else {
        trim.uv[slot].clone().map(Real::from)
    };
    Ok((point(&face.surface, &at, dec)?, t.clone()))
}

pub(in crate::numerical) fn samples(
    face_index: usize,
    face: &NumericalFace,
    edges: &[NumericalEdge],
    index: &Index<'_>,
    dec: &mut Decider,
) -> Result<Vec<NumericalCheck>, Refusal> {
    if !index.has_face(face_index) {
        return Ok(Vec::new());
    }
    let mut checks = Vec::new();
    for (boundary_index, boundary) in face.loops.iter().enumerate() {
        dec.charge_work(1)?;
        for (use_index, trim) in boundary.iter().enumerate() {
            dec.charge_work(1)?;
            let next_index = (use_index + 1) % boundary.len();
            let next = &boundary[next_index];
            let previous = endpoint(
                face,
                trim,
                &edges[trim.edge],
                index.get(face_index, boundary_index, use_index),
                true,
                dec,
            )?;
            let following = endpoint(
                face,
                next,
                &edges[next.edge],
                index.get(face_index, boundary_index, next_index),
                false,
                dec,
            )?;
            checks.push(NumericalCheck::Samples {
                name: if index.numerical_pair(trim.edge) || index.numerical_pair(next.edge) {
                    "trim.numerical-pair-native-join-position"
                } else {
                    "trim.native-join-position"
                }.into(),
                scope: format!(
                    "face {:?}, loop {boundary_index}, uses {use_index}->{next_index}",
                    face.key
                ),
                samples: vec![NumericalSample {
                    parameters: vec![previous.1, following.1],
                    distance_upper_mm: distance(&previous.0, &following.0),
                }],
            });
        }
    }
    Ok(checks)
}
