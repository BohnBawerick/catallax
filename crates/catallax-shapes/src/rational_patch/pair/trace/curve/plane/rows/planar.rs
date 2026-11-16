//! Solve an owned source row directly against a proved planar support.

use super::*;

pub(super) fn point(
    piece: &PairTracePiece,
    span: &Span,
    other: &Span,
    side: usize,
    axis: usize,
    row: &Real,
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    let middle = other
        .domain
        .each_ref()
        .map(|[lo, hi]| Real::from((lo + hi).mul_pow2(-1)));
    let [u, v] = other.derivatives(&middle);
    let plane = Plane {
        origin: other.point(&middle),
        normal: u.cross(&v),
        ref_dir: u,
    };
    if !equation(other, &plane).is_zero(dec)? {
        return Ok(None);
    }
    let cut = equation(span, &plane).substitute(axis, &RealPoly::constant(2, row.clone()));
    if cut.is_zero(dec)? {
        return Ok(None);
    }
    let cut = RealPoly::from_terms(
        1,
        cut.terms()
            .into_iter()
            .map(|(m, c)| (vec![m.exponents()[1 - axis]], c))
            .collect(),
    );
    let arc = piece.arc();
    let (source, unit, domain, source_domain) = if side == 0 {
        (
            &piece.definition.second,
            arc.second_chart,
            &arc.second_cell,
            &arc.first_cell,
        )
    } else {
        (
            &piece.definition.first,
            arc.first_chart,
            &arc.first_cell,
            &arc.second_cell,
        )
    };
    let mut found = None;
    for root in cut.real_roots(dec)? {
        let at = if axis == 0 {
            [row.clone(), root]
        } else {
            [root, row.clone()]
        };
        if !contains(source_domain, &at, dec)? {
            continue;
        }
        let point = span.point(&at);
        let mut incident = false;
        for (_, patch) in source.query_charts(unit) {
            if let Some(other_at) = patch.parameters_at(&point, dec)? {
                incident |= contains(domain, &other_at, dec)?;
            }
        }
        if incident {
            if found.is_some() {
                return Err(missing(
                    "a monotone source row acquired two owned plane roots",
                ));
            }
            found = Some(point);
        }
    }
    found
        .map(Some)
        .ok_or_else(|| missing("an interior monotone source row lost its proved plane root"))
}
