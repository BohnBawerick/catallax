//! A univariate cut on a proved monotone source parameter.

use super::*;

pub(super) fn points(
    piece: &PairTracePiece,
    span: &Span,
    cut: &RealPoly,
    side: usize,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let terms = cut.terms();
    let Some(axis) = (0..2).find(|&axis| terms.iter().all(|(m, _)| m.exponents()[1 - axis] == 0))
    else {
        return Ok(None);
    };
    let arc = piece.arc();
    let other = if side == 0 {
        piece
            .definition
            .second
            .span(arc.second_chart, arc.second_span)
    } else {
        piece.first()
    };
    let tangent: [RealPoly; 3] = std::array::from_fn(|i| {
        &span.first[1 - axis][i] * &span.coordinates[3]
            - &span.coordinates[i] * &span.first[1 - axis][3]
    });
    let targets = if side == 0 { [0, 1] } else { [2, 3] };
    let others = if side == 0 { [2, 3] } else { [0, 1] };
    let velocity = (0..3).fold(RealPoly::zero(4), |sum, i| {
        sum + lift(&tangent[i], &targets, 4) * lift(&other.normal_polynomials()[i], &others, 4)
    });
    let domain = [
        arc.first_cell[0].clone(),
        arc.first_cell[1].clone(),
        arc.second_cell[0].clone(),
        arc.second_cell[1].clone(),
    ];
    let bits = 128.min(dec.policy().max_precision_bits);
    let box_: Vec<_> = domain
        .iter()
        .map(|[lo, hi]| Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits)))
        .collect();
    if bounds::eval(&velocity, &box_, dec)?.contains_zero() {
        return Ok(None);
    }
    let column = 2 * side + axis;
    let ends = [
        &piece.endpoint_parameters(0)[column],
        &piece.endpoint_parameters(1)[column],
    ];
    let order = exact_sign(&(ends[1] - ends[0]), dec)?;
    if order == Sign::Zero {
        return Err(missing(
            "strict source monotonicity conflicts with equal endpoints",
        ));
    }
    let equation = RealPoly::from_terms(
        1,
        terms
            .into_iter()
            .map(|(m, c)| (vec![m.exponents()[axis]], c))
            .collect(),
    );
    let mut points = Vec::new();
    for root in equation.real_roots(dec)? {
        let low = exact_sign(&(&root - ends[0]), dec)?;
        let high = exact_sign(&(ends[1] - &root), dec)?;
        if low == Sign::Zero {
            points.push(endpoint(piece, 0));
        } else if high == Sign::Zero {
            points.push(endpoint(piece, 1));
        } else if low == order && high == order {
            let mut equations: Vec<_> = piece
                .system
                .polynomials
                .iter()
                .map(|p| {
                    RealPoly::from_terms(
                        4,
                        p.terms()
                            .into_iter()
                            .map(|(m, c)| (m.exponents()[..4].to_vec(), c))
                            .collect(),
                    )
                })
                .collect();
            equations.push(RealPoly::var(4, column) - RealPoly::constant(4, root));
            // Strict source monotonicity and the endpoint bracket prove one
            // root on this owned arc. Reconstruction selects that joint root.
            let at = RealPoly::simple_algebraic_root_in_reconstructed(
                &equations,
                &domain,
                &[0, 1],
                bits,
                dec,
            )?;
            points.push(piece.first().point(&[at[0].clone(), at[1].clone()]));
        }
    }
    Ok(Some(points))
}
