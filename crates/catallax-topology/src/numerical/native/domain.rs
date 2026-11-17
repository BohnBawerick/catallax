//! Complete parameter-image hulls for regular native analytic charts.

use super::*;

pub(super) fn check(
    face_index: usize,
    face: &NumericalFace,
    edges: &[NumericalEdge],
    index: &Index<'_>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let first = face
        .loops
        .iter()
        .flatten()
        .next()
        .ok_or_else(|| bad("trim", "a native analytic chart needs a finite boundary"))?;
    let mut bounds: [[Rational; 2]; 2] =
        std::array::from_fn(|axis| [first.uv[0][axis].clone(), first.uv[0][axis].clone()]);
    for (row, boundary) in face.loops.iter().enumerate() {
        for (use_index, trim) in boundary.iter().enumerate() {
            for point in &trim.uv {
                for axis in 0..2 {
                    if point[axis] < bounds[axis][0] {
                        bounds[axis][0] = point[axis].clone();
                    }
                    if point[axis] > bounds[axis][1] {
                        bounds[axis][1] = point[axis].clone();
                    }
                }
            }
            let Some(curve) = index.get(face_index, row, use_index) else {
                continue;
            };
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "native chart boundary names an unknown edge"))?;
            let [from, to] = if edge.parameters[0] < edge.parameters[1] {
                edge.parameters.each_ref()
            } else {
                [&edge.parameters[1], &edge.parameters[0]]
            };
            let mut covered = false;
            for (span, domain) in curve.spans() {
                let lo = domain[0].clone().max(from.clone());
                let hi = domain[1].clone().min(to.clone());
                if lo >= hi {
                    continue;
                }
                covered = true;
                for control in curve.span_bezier_controls_on(span, [&lo, &hi], dec)? {
                    dec.charge_work(1)?;
                    let weight = dec.enclose_at(&control[3], BITS)?;
                    if weight.lo().sign() != Sign::Positive {
                        return Err(Refusal::new(
                            "numerical.native-domain-weight-unverified",
                            "the native parameter hull needs separated positive rational weights",
                        ));
                    }
                    for axis in 0..2 {
                        let value = dec
                            .enclose_at(&control[axis], BITS)?
                            .checked_div(&weight)
                            .ok_or_else(|| {
                                Refusal::new(
                                    "numerical.native-domain-weight-unverified",
                                    "the native parameter hull divisor is not separated from zero",
                                )
                            })?;
                        let (lo, hi) = (value.lo(), value.hi());
                        if lo < bounds[axis][0] {
                            bounds[axis][0] = lo;
                        }
                        if hi > bounds[axis][1] {
                            bounds[axis][1] = hi;
                        }
                    }
                }
            }
            if !covered {
                return Err(bad(
                    "trim",
                    "native trim has no nonempty retained span to bound",
                ));
            }
        }
    }
    analytic::check_native_domain(&face.surface, &bounds)
}
