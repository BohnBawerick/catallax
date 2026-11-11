//! The carrier's display grid, cached per tolerance.

use super::*;
use crate::rational_sheet::sheet_grid;

impl RollingEnvelope {
    /// Grid with `chord` at most `tolerance`. Doubling keeps every carrier knot
    /// on a grid line. A tolerance out of reach after 2^10 cells per span
    /// direction refuses instead of widening the report.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        let local = self.local_grid(tolerance, max_points, dec)?;
        let frame = self.data.definition.frame(self.data.net.bits, dec)?;
        Ok(EnvelopeGrid {
            columns: local.columns,
            rows: local.rows,
            points: local.points.iter().map(|p| frame.world(p)).collect(),
            chord: local.chord.clone(),
        })
    }

    fn local_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<Arc<EnvelopeGrid>, Refusal> {
        // A worker leaves the shared cache unread and unfilled for the writer.
        let mut cache = catallax_numbers::shared_caches()
            .then(|| self.data.net.grid.lock().expect("no panic holds this lock"));
        #[cfg(feature = "work-span")]
        let span_owner = std::ptr::from_ref(&self.data.net.grid) as usize as u64;
        if let Some((key, grid)) = cache.as_deref().and_then(Option::as_ref)
            && key == tolerance
        {
            // A hit skips the grid's charge, so it waits for the fill.
            #[cfg(feature = "work-span")]
            catallax_numbers::work_span::read_at("shared cache: envelope grid", span_owner, 0);
            return Ok(grid.clone());
        }
        let grid = Arc::new(sheet_grid(
            Patch::all(&self.data.carrier, &self.data.net.controls, DEGREE),
            tolerance,
            max_points,
            dec,
        )?);
        if let Some(cache) = cache.as_deref_mut() {
            *cache = Some((tolerance.clone(), grid.clone()));
            #[cfg(feature = "work-span")]
            catallax_numbers::work_span::write_at("shared cache: envelope grid", span_owner, 0);
        }
        Ok(grid)
    }
}
