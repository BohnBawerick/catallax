//! Original tensor spans with their retained, certified atlas chart cover.

use super::*;
use crate::{RationalAtlas, SplineBoundary};
use std::collections::BTreeMap;

const SIDES: [SplineBoundary; 4] = [
    SplineBoundary::UMin,
    SplineBoundary::UMax,
    SplineBoundary::VMin,
    SplineBoundary::VMax,
];

#[derive(Debug)]
pub(super) struct OriginalUnit {
    pub(super) span: Span,
    pub(super) physical: [bool; 4],
    pub(super) charts: Vec<usize>,
}

#[derive(Debug)]
pub(super) struct OriginalAtlas {
    pub(super) atlas: RationalAtlas,
    pub(super) units: Vec<OriginalUnit>,
}

fn failure(detail: &'static str) -> Refusal {
    Refusal::new("intersect.original-span-ownership", detail)
}

impl OriginalAtlas {
    fn new(atlas: &RationalAtlas, dec: &mut Decider) -> Result<Arc<Self>, Refusal> {
        let source = atlas.source();
        let knots = source.knots();
        let mut covers: BTreeMap<[usize; 2], Vec<usize>> = BTreeMap::new();
        for (index, chart) in atlas.charts().iter().enumerate() {
            dec.check_deadline()?;
            let bounds = chart.patch().domain();
            let mut indices = [0; 2];
            for axis in 0..2 {
                let Some(span) = knots[axis]
                    .partition_point(|k| k <= &bounds[axis][0])
                    .checked_sub(1)
                else {
                    return Err(failure(
                        "a certified chart starts outside its original knot domain",
                    ));
                };
                let Some(end) = knots[axis].get(span + 1) else {
                    return Err(failure(
                        "a certified chart has no containing original knot interval",
                    ));
                };
                if knots[axis][span] >= *end || &bounds[axis][1] > end {
                    return Err(failure(
                        "a certified chart crosses an original tensor-span boundary",
                    ));
                }
                indices[axis] = span;
            }
            covers.entry(indices).or_default().push(index);
        }
        let mut units = Vec::with_capacity(covers.len());
        for (indices, charts) in covers {
            dec.check_deadline()?;
            let domain: Domain = std::array::from_fn(|axis| {
                [
                    knots[axis][indices[axis]].clone(),
                    knots[axis][indices[axis] + 1].clone(),
                ]
            });
            let mut area = Rational::from_i64(0);
            let mut physical = [None; 4];
            for &index in &charts {
                let chart = &atlas.charts()[index];
                let bounds = chart.patch().domain();
                let size = &(&bounds[0][1] - &bounds[0][0]) * &(&bounds[1][1] - &bounds[1][0]);
                area = &area + &size;
                for (side, flag) in physical.iter_mut().enumerate() {
                    if bounds[side / 2][side % 2] != domain[side / 2][side % 2] {
                        continue;
                    }
                    let value = chart.neighbor(SIDES[side]).is_none();
                    if flag.is_some_and(|previous| previous != value) {
                        return Err(failure(
                            "an original span edge has nonuniform physical or periodic ownership",
                        ));
                    }
                    *flag = Some(value);
                }
            }
            // The atlas constructor supplies disjoint chart interiors. Exact
            // area equality proves this retained component owns the full span;
            // a partial logical face is never silently extended to its source.
            let full = &(&domain[0][1] - &domain[0][0]) * &(&domain[1][1] - &domain[1][0]);
            if area != full {
                return Err(failure(
                    "the retained chart component does not completely cover an original span",
                ));
            }
            let [Some(a), Some(b), Some(c), Some(d)] = physical else {
                return Err(failure(
                    "the complete original-span cover has no proved edge ownership",
                ));
            };
            let coordinates = source.span_polynomials(indices, dec)?;
            units.push(OriginalUnit {
                span: Span::new(indices, domain, coordinates),
                physical: [a, b, c, d],
                charts,
            });
        }
        if units.is_empty() {
            return Err(failure(
                "an original-span meeting needs a nonempty certified source cover",
            ));
        }
        Ok(Arc::new(Self {
            atlas: atlas.clone(),
            units,
        }))
    }
}

impl RationalAtlas {
    /// Complete meeting on original tensor spans, retaining every certified
    /// atlas chart as cover, regularity and endpoint-ownership authority.
    /// The existing field, normal-cone and root search is unchanged. Partial
    /// span covers or nonuniform edge ownership refuse by name.
    ///
    /// Source-unit indices in returned arcs refer to original spans; inverse
    /// point queries still visit their actual certified embedding charts.
    pub fn patch_meeting_with_original_spans(
        &self,
        other: &Self,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        let budget = catallax_numbers::operation_budget_active();
        let start = if budget {
            catallax_numbers::current_operation_work().used
        } else {
            dec.stats().work
        };
        let first = OriginalAtlas::new(self, dec)?;
        let second = OriginalAtlas::new(other, dec)?;
        let mut answer =
            prepared::meeting_cached([Source::Original(&first), Source::Original(&second)], dec)?;
        answer.search.work = if budget {
            catallax_numbers::current_operation_work().used
        } else {
            dec.stats().work
        } - start;
        Ok(answer)
    }
}
