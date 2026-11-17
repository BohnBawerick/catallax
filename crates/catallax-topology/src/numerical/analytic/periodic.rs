use super::*;

pub(super) fn check_loops(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let (axes, family, code) = match &face.surface {
        Surface::Torus(_) => (2, "torus", "numerical.torus-periodic-domain-unverified"),
        Surface::Sphere(_) => (1, "sphere", "numerical.sphere-periodic-domain-unverified"),
        _ => {
            return Err(bad(
                "unsupported",
                "periodic chart checks require a sphere or ring torus",
            ));
        }
    };
    for axis in 0..axes {
        let mut values = face
            .loops
            .iter()
            .flatten()
            .flat_map(|trim| &trim.uv)
            .map(|uv| &uv[axis]);
        let first = values
            .next()
            .ok_or_else(|| bad("trim", "a periodic face has no finite parameter boundary"))?;
        let (mut lo, mut hi) = (first, first);
        for value in values {
            dec.charge_work(1)?;
            lo = lo.min(value);
            hi = hi.max(value);
        }
        let width = hi - lo;
        if width < Ball::pi(BITS).lo() * Rational::from_i64(2) {
            continue;
        }
        if width != period() || !paired_seam(face, edges, axis, dec)? {
            return Err(Refusal::new(
                code,
                format!(
                    "a {family} chart must stay below one turn in each periodic parameter or retain a rectangular paired single-period cell with strictly interior holes"
                ),
            ));
        }
    }
    for (index, first) in face.loops.iter().flatten().enumerate() {
        for second in face.loops.iter().flatten().skip(index + 1) {
            dec.charge_work(1)?;
            if first.edge == second.edge
                && !(0..axes).any(|axis| translated_seam(first, second, axis))
            {
                return Err(Refusal::new(
                    code,
                    format!(
                        "a repeated {family} edge must have opposite uses translated by one complete parameter period"
                    ),
                ));
            }
        }
    }
    Ok(())
}

fn translated_seam(first: &NumericalTrim, second: &NumericalTrim, axis: usize) -> bool {
    if first.edge != second.edge || first.sense == second.sense {
        return false;
    }
    let other = 1 - axis;
    let shift = &second.uv[1][axis] - &first.uv[0][axis];
    shift.abs() == period()
        && &second.uv[0][axis] - &first.uv[1][axis] == shift
        && second.uv[1][other] == first.uv[0][other]
        && second.uv[0][other] == first.uv[1][other]
}

fn paired_seam(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    axis: usize,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let other = 1 - axis;
    let Some(boundary) = face.loops.first() else {
        return Ok(false);
    };
    // Only a rectangular paired cell is covered by this certificate. A seam
    // on one interval cannot authorize other identified contacts after wrapping.
    if boundary.len() != 4 {
        return Ok(false);
    }
    for (i, first) in boundary.iter().enumerate() {
        for (j, second) in boundary.iter().enumerate().skip(i + 1) {
            dec.charge_work(1)?;
            if first.edge != second.edge || first.sense == second.sense {
                continue;
            }
            if edges.get(first.edge).is_none() {
                return Err(bad("trim", "periodic seam names an unknown edge"));
            }
            if translated_seam(first, second, axis)
                && first.uv[0][axis] == first.uv[1][axis]
                && second.uv[0][axis] == second.uv[1][axis]
                && first.uv[0][other] != first.uv[1][other]
                && boundary.iter().enumerate().all(|(index, rim)| {
                    index == i
                        || index == j
                        || (rim.uv[0][other] == rim.uv[1][other]
                            && (rim.uv[0][other] == first.uv[0][other]
                                || rim.uv[0][other] == first.uv[1][other])
                            && (&rim.uv[1][axis] - &rim.uv[0][axis]).abs() == period())
                })
            {
                let low = (&first.uv[0][axis]).min(&second.uv[0][axis]);
                let high = (&first.uv[0][axis]).max(&second.uv[0][axis]);
                if face
                    .loops
                    .iter()
                    .skip(1)
                    .flatten()
                    .flat_map(|trim| &trim.uv)
                    .any(|uv| &uv[axis] <= low || &uv[axis] >= high)
                {
                    return Ok(false);
                }
                return Ok(true);
            }
        }
    }
    Ok(false)
}
