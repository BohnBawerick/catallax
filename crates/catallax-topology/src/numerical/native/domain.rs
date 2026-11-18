//! Complete parameter-image hulls for native analytic and retained spline charts.

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
    for point in face.loops.iter().flatten().flat_map(|trim| &trim.uv) {
        for axis in 0..2 {
            if point[axis] < bounds[axis][0] {
                bounds[axis][0] = point[axis].clone();
            }
            if point[axis] > bounds[axis][1] {
                bounds[axis][1] = point[axis].clone();
            }
        }
    }
    let cone_turn =
        matches!(face.surface, Surface::Cone(_)).then(|| Ball::pi(BITS).lo().mul_pow2(1));
    let cell = if cone_turn
        .as_ref()
        .is_some_and(|turn| &bounds[0][1] - &bounds[0][0] >= *turn)
    {
        analytic::periodic_cell(face, edges, dec)?
    } else {
        None
    };
    for (row, boundary) in face.loops.iter().enumerate() {
        for (use_index, trim) in boundary.iter().enumerate() {
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
                    if let Some(cell) = &cell
                        && !cell.contains_control(&control, dec)?
                    {
                        return Err(Refusal::new(
                            "numerical.native-cone-domain-unverified",
                            format!(
                                "face {} edge {}: the complete native-curve hull leaves its paired affine fundamental cell",
                                face.key, edge.key
                            ),
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
    if let Some(turn) = cone_turn {
        if bounds[1][0].sign() != Sign::Positive {
            return Err(Refusal::new(
                "numerical.native-cone-domain-unverified",
                "the complete native-curve hull must stay strictly above the cone apex",
            ));
        }
        let width = &bounds[0][1] - &bounds[0][0];
        if width >= turn && cell.is_none() {
            return Err(Refusal::new(
                "numerical.native-cone-domain-unverified",
                "the complete native-curve hull needs a regular angular chart or one paired fundamental cell",
            ));
        }
        return Ok(());
    }
    if let Surface::BSpline(source) = &face.surface {
        let domain = source.carrier().domain();
        if (0..2)
            .any(|axis| bounds[axis][0] < *domain[axis][0] || bounds[axis][1] > *domain[axis][1])
        {
            return Err(Refusal::new(
                "numerical.native-spline-domain-unverified",
                "the complete native-curve hull leaves its retained spline carrier domain",
            ));
        }
        Ok(())
    } else {
        analytic::check_native_domain(&face.surface, &bounds)
    }
}
