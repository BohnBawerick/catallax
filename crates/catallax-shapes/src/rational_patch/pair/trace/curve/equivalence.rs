//! Component identity uses a shared parameter lift, never nearby world boxes.

use super::*;

impl PatchPairCurve {
    /// Exact point-set equivalence when retained component ownership proves it.
    /// A shared parameter lift joins maximal components of the same complete
    /// meeting. A missing world witness proves different sets. Otherwise the
    /// available certificates do not decide the query.
    pub fn same_component(&self, other: &Self, dec: &mut Decider) -> Result<Option<bool>, Refusal> {
        if self.same_definition(other) {
            return Ok(Some(true));
        }
        let a = &self.pieces[0].definition;
        let b = &other.pieces[0].definition;
        let swapped = if a.first.same_definition(&b.first) && a.second.same_definition(&b.second) {
            Some(false)
        } else if a.first.same_definition(&b.second) && a.second.same_definition(&b.first) {
            Some(true)
        } else {
            None
        };
        if let Some(swapped) = swapped {
            let seed = self.pieces[0].arc();
            let parameters = seed.ends[0].parameters();
            let (charts, first, second) = if swapped {
                (
                    [seed.second_chart, seed.first_chart],
                    [parameters[2].clone(), parameters[3].clone()],
                    [parameters[0].clone(), parameters[1].clone()],
                )
            } else {
                (
                    [seed.first_chart, seed.second_chart],
                    [parameters[0].clone(), parameters[1].clone()],
                    [parameters[2].clone(), parameters[3].clone()],
                )
            };
            for piece in other.pieces.iter() {
                let arc = piece.arc();
                if [arc.first_chart, arc.second_chart] == charts
                    && contains(&arc.first_cell, &first, dec)?
                    && contains(&arc.second_cell, &second, dec)?
                {
                    return Ok(Some(true));
                }
            }
        }
        // Atlas charts are individually embedded; unrelated charts need not
        // have a unique global lift. A failed lift match alone proves nothing.
        for (a, b) in [(self, other), (other, self)] {
            for point in a.witnesses() {
                if !b.contains(&point, dec)? {
                    return Ok(Some(false));
                }
            }
        }
        Ok(None)
    }
}
