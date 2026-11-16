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
        let trace_start = catallax_numbers::current_operation_work().used;
        eprintln!("PAIR_CLEAR enter work={trace_start}");
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
            eprintln!("PAIR_CLEAR side={side} patches={} before={} delta={}", patches.len(), catallax_numbers::current_operation_work().used, catallax_numbers::current_operation_work().used-trace_start);
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
            eprintln!("PAIR_CLEAR side={side} queried work={} delta={} kind={}", catallax_numbers::current_operation_work().used, catallax_numbers::current_operation_work().used-trace_start, match &answer {Some(SplineParameterHits::Points(p))=>p.len() as i64,Some(SplineParameterHits::Contained)=>-1,None=>-2});
            if let Some(SplineParameterHits::Points(parameters)) = answer {
                let mut clear = true;
                let initial = Real::from(spline.domain()[0].clone());
                for parameter in parameters {
                    eprintln!("PAIR_CLEAR compare work={}",catallax_numbers::current_operation_work().used);
                    dec.charge_work(1)?;
                    if !allow_initial_contact
                        || exact_sign(&(&parameter - &initial), dec)? != Sign::Zero
                    {
                        clear = false;
                        break;
                    }
                }
                if clear {
                    eprintln!("PAIR_CLEAR clear work={}",catallax_numbers::current_operation_work().used);
                    return Ok(true);
                }
            }
        }
        eprintln!("PAIR_CLEAR inconclusive work={}",catallax_numbers::current_operation_work().used);
        if let Some(error) = partial {
            return Err(error);
        }
        Ok(false)
    }
}
