//! A sheet is the regularized difference of two existing exact networks.
use crate::{Plan, geometry, gyroid::Gyroid};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{BooleanKind, Grade, GyroidKind, Number, Pose};

/// Compile both sides, retaining the network evaluator's coverage restrictions.
pub(crate) fn sides(
    period: &Rational,
    level: &Grade,
    frame: &Pose,
) -> Result<[Gyroid; 2], Refusal> {
    let negate = |n: &Number| Number::new(-n.value());
    let positive = |n: &Number| n.value() > &Rational::zero();
    let negative = match level {
        Grade::Constant(c) if positive(c) => Grade::Constant(negate(c)),
        Grade::Linear {
            from,
            to,
            start,
            end,
        } if positive(start) && positive(end) => Grade::Linear {
            from: from.clone(),
            to: to.clone(),
            start: negate(start),
            end: negate(end),
        },
        Grade::Distance { .. } => {
            return Err(Refusal::new(
                "gyroid.unsupported",
                "distance-graded sheets need a certified grade evaluator",
            ));
        }
        _ => {
            return Err(Refusal::new(
                "gyroid.impossible",
                "a sheet needs a strictly positive level everywhere",
            ));
        }
    };
    Ok([
        Gyroid::certified(period, level, GyroidKind::Network, frame)?,
        Gyroid::certified(period, &negative, GyroidKind::Network, frame)?,
    ])
}
impl Plan {
    /// Only live leaves introduced by an explicit sheet operation need the new
    /// print guard. Existing network Boolean recipes retain their prior route.
    pub(crate) fn has_sheet_band(&self) -> bool {
        self.leaves.iter().enumerate().any(|(i, leaf)| {
            self.leaf_registers[i].is_some()
                && matches!(leaf,geometry::Leaf::Gyroid(g) if g.is_sheet_side())
        })
    }
    pub(crate) fn sheet(
        &mut self,
        region: usize,
        period: &Rational,
        level: &Grade,
        frame: &Pose,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let [mut upper, mut lower] = sides(period, level, frame)?;
        if upper.saturated() == Some(true) && lower.saturated() == Some(false) {
            return Ok(region);
        }
        upper.mark_sheet_side();
        lower.mark_sheet_side();
        let a = self.leaf(geometry::Leaf::Gyroid(Box::new(upper)), vec![], dec)?;
        let b = self.leaf(geometry::Leaf::Gyroid(Box::new(lower)), vec![], dec)?;
        let band = self.boolean(BooleanKind::Subtract, a, b, dec)?;
        self.boolean(BooleanKind::Intersect, region, band, dec)
    }
}
