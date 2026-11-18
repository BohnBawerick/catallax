use super::*;
use catallax_numbers::Real;

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
        if width != period() || paired_seam(face, edges, axis, dec)?.is_none() {
            return Err(Refusal::new(
                code,
                format!(
                    "a {family} chart must stay below one turn in each periodic parameter or retain a paired single-period fundamental cell with strictly interior holes"
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

pub(in crate::numerical) struct Cell {
    axis: usize,
    origin: [Rational; 2],
    height: Rational,
    drift: Rational,
    bounds: [Rational; 2],
}

impl Cell {
    /// Positive-weight homogeneous controls stay in the translated affine
    /// cell exactly, without projecting a sloping cell onto a wider u box.
    pub(in crate::numerical) fn contains_control(
        &self,
        control: &[Real; 4],
        dec: &Decider,
    ) -> Result<bool, Refusal> {
        let offset =
            |axis: usize| &control[axis] - &(&control[3] * &Real::from(self.origin[axis].clone()));
        let across = &(&offset(self.axis) * &Real::from(self.height.clone()))
            - &(&offset(1 - self.axis) * &Real::from(self.drift.clone()));
        let lower = &across - &(&control[3] * &Real::from(self.bounds[0].clone()));
        let upper = &across - &(&control[3] * &Real::from(self.bounds[1].clone()));
        Ok(dec.enclose_at(&lower, BITS)?.lo().sign() != Sign::Negative
            && dec.enclose_at(&upper, BITS)?.hi().sign() != Sign::Positive)
    }
}

pub(super) fn paired_seam(
    face: &NumericalFace,
    edges: &[NumericalEdge],
    axis: usize,
    dec: &mut Decider,
) -> Result<Option<Cell>, Refusal> {
    let other = 1 - axis;
    let Some(boundary) = face.loops.first() else {
        return Ok(None);
    };
    if boundary.len() < 4 {
        return Ok(None);
    }
    for (i, first) in boundary.iter().enumerate() {
        for (j, second) in boundary.iter().enumerate().skip(i + 1) {
            dec.charge_work(1)?;
            if first.edge != second.edge || first.sense == second.sense {
                continue;
            }
            let Some(edge) = edges.get(first.edge) else {
                return Err(bad("trim", "periodic seam names an unknown edge"));
            };
            if matches!(face.surface, Surface::Cylinder(_) | Surface::Cone(_))
                && !matches!(edge.curve, Curve::Line(_))
            {
                continue;
            }
            if !translated_seam(first, second, axis) || first.uv[0][other] == first.uv[1][other] {
                continue;
            }
            // Subtract the affine seam direction before measuring across the
            // cell. Decimal frame rounding may give both seam branches a slope.
            let height = &first.uv[1][other] - &first.uv[0][other];
            let drift = &first.uv[1][axis] - &first.uv[0][axis];
            let across = |point: &[Rational; 2]| {
                &(&point[axis] - &first.uv[0][axis]) * &height
                    - &(&point[other] - &first.uv[0][other]) * &drift
            };
            let shift = across(&second.uv[0]);
            let low = Rational::zero().min(shift.clone());
            let high = Rational::zero().max(shift);
            if !boundary.iter().enumerate().all(|(index, rim)| {
                if index == i || index == j {
                    return true;
                }
                let a = across(&rim.uv[0]);
                let b = across(&rim.uv[1]);
                [(&a, &rim.uv[0]), (&b, &rim.uv[1])]
                    .iter()
                    .all(|(value, point)| {
                        **value >= low
                            && **value <= high
                            && (**value != low && **value != high
                                || first.uv.contains(point)
                                || second.uv.contains(point))
                    })
                    && !(a == b && (a == low || a == high))
            }) {
                continue;
            }
            if face
                .loops
                .iter()
                .skip(1)
                .flatten()
                .flat_map(|trim| &trim.uv)
                .any(|point| {
                    let value = across(point);
                    value <= low || value >= high
                })
            {
                return Ok(None);
            }
            return Ok(Some(Cell {
                axis,
                origin: first.uv[0].clone(),
                height,
                drift,
                bounds: [low, high],
            }));
        }
    }
    Ok(None)
}
