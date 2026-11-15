//! Complete isolated path meetings from either retained source cover.

use super::*;
use crate::{BSplineCurve, SectionHits};

impl PatchPairCurve {
    /// Complete isolated meetings of a finite native spline and this component.
    /// `None` means no source cover supplies a point-only query. Containment is
    /// never replaced by selected endpoints or an empty point set.
    pub fn spline_points(
        &self,
        spline: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Point3>>, Refusal> {
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
            let answer = match crate::rational_patch::spline_points(&patches, spline, dec) {
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
            if let Some(SectionHits::Points(candidates)) = answer {
                let mut points = Vec::new();
                for point in candidates {
                    dec.charge_work(1)?;
                    if self.contains(&point, dec)? {
                        points.push(point);
                    }
                }
                return Ok(Some(points));
            }
        }
        if let Some(error) = partial {
            return Err(error);
        }
        Ok(None)
    }
}
