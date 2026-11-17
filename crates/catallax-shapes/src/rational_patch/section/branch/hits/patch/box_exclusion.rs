//! Exact empty-set certificates over a bounded Bezier subdivision tree.

use super::*;

const MAX_BOX_DEPTH: u32 = 8;
const MAX_BOX_CELLS: usize = 128;

struct CutBoxSearch<'a> {
    cut: &'a SectionSurface,
    patch: &'a RationalPatch,
    bounds: &'a [Interval; 3],
    visited: usize,
}

impl CutBoxSearch<'_> {
    fn cell(
        &mut self,
        span: [usize; 2],
        domain: &Domain,
        depth: u32,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        dec.charge_work(1)?;
        if self.visited == MAX_BOX_CELLS {
            return Ok(false);
        }
        self.visited += 1;
        {
            let controls = self.patch.carrier().span_bezier_controls_on(
                span,
                domain.each_ref().map(|[lo, hi]| [lo, hi]),
                dec,
            )?;
            if let Some(hull) = homogeneous_control_hull(&controls, dec)?
                && (0..3).any(|axis| {
                    hull[axis].hi() < self.bounds[axis].lo()
                        || self.bounds[axis].hi() < hull[axis].lo()
                })
            {
                return Ok(true);
            }
            if !self.cut.retains_span(&controls, dec)? {
                return Ok(true);
            }
            let (values, _) = self
                .cut
                .net(&controls, self.patch.carrier().degrees(), dec)?;
            if uniform(&values, dec)?.is_some() {
                return Ok(true);
            }
        }
        if depth == MAX_BOX_DEPTH {
            return Ok(false);
        }
        let [low, high] = split(domain, &self.patch.data.widths);
        if !self.cell(span, &low, depth + 1, dec)? {
            return Ok(false);
        }
        self.cell(span, &high, depth + 1, dec)
    }
}

impl SectionSurface {
    /// Proves that no point of `patch` on this cut lies in the closed box.
    ///
    /// Each original knot rectangle is covered by exact positive-weight Bezier
    /// subrectangles. A cell disappears only through strict outward Cartesian
    /// hull separation, retained-sheet exclusion, or an exact single-signed cut
    /// net. Every cell must disappear. The finite depth/cell limits return
    /// `false` when inconclusive, never a numerical empty-set assertion.
    pub fn excludes_patch_in_box(
        &self,
        patch: &RationalPatch,
        bounds: &[Interval; 3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        self.validate(dec)?;
        let mut search = CutBoxSearch {
            cut: self,
            patch,
            bounds,
            visited: 0,
        };
        for span in &patch.data.spans {
            if !search.cell(span.indices, &span.domain, 0, dec)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
