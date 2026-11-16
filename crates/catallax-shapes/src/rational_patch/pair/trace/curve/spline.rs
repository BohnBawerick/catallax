//! Sufficient path clearance against either complete retained source cover.

use super::*;
use crate::BSplineCurve;
use crate::rational_patch::{SplineParameterHits, spline_parameters};

impl PatchPairCurve {
    /// Proves that the finite path avoids a complete source cover of this
    /// component. Only its initial native parameter may be permitted to meet.
    /// False is inconclusive: contained/unsupported queries and meetings with
    /// unrelated parts of a source cover never become an empty point set.
    pub fn clears_spline_sources(
        &self,
        spline: &BSplineCurve,
        allow_initial_contact: bool,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut partial = None;
        for side in 0..2 {
            let mut patches = Vec::new();
            for piece in self.pieces.iter() {
                dec.charge_work(1)?;
                let arc = piece.arc();
                let (source, chart) = if side == 0 {
                    (&piece.definition.first, arc.first_chart)
                } else {
                    (&piece.definition.second, arc.second_chart)
                };
                for (_, patch) in source.query_charts(chart) {
                    dec.charge_work(1)?;
                    patches.push(patch);
                }
            }
            let answer = match spline_parameters(&patches, spline, dec) {
                Ok(answer) => answer,
                Err(error)
                    if error.code == "shapes.rational-patch-unproven"
                        && error.detail
                            == "partial spline containment needs finite interval ownership" =>
                {
                    partial = Some(error);
                    continue;
                }
                Err(error) => return Err(error),
            };
            if let Some(SplineParameterHits::Points(parameters)) = answer {
                let mut clear = true;
                let initial = Real::from(spline.domain()[0].clone());
                for parameter in parameters {
                    dec.charge_work(1)?;
                    if !allow_initial_contact
                        || exact_sign(&(parameter - &initial), dec)? != Sign::Zero
                    {
                        clear = false;
                        break;
                    }
                }
                if clear {
                    return Ok(true);
                }
            }
        }
        if let Some(error) = partial {
            return Err(error);
        }
        Ok(false)
    }
}
