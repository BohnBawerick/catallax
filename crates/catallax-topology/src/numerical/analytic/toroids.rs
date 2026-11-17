use super::*;

pub(super) fn check_loops(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for axis in 0..2 {
        let mut values = face
            .loops
            .iter()
            .flatten()
            .flat_map(|trim| &trim.uv)
            .map(|uv| &uv[axis]);
        let first = values
            .next()
            .ok_or_else(|| bad("trim", "torus face has no finite parameter boundary"))?;
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
                "numerical.torus-periodic-domain-unverified",
                "a torus chart must stay below one turn in each parameter or retain an explicitly paired single-period seam",
            ));
        }
    }
    Ok(())
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
    for (i, first) in boundary.iter().enumerate() {
        for second in &boundary[i + 1..] {
            dec.charge_work(1)?;
            if first.edge != second.edge || first.sense == second.sense {
                continue;
            }
            if edges.get(first.edge).is_none() {
                return Err(bad("trim", "torus seam names an unknown edge"));
            }
            let shift = &second.uv[1][axis] - &first.uv[0][axis];
            if shift.abs() == period()
                && &second.uv[0][axis] - &first.uv[1][axis] == shift
                && second.uv[1][other] == first.uv[0][other]
                && second.uv[0][other] == first.uv[1][other]
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
