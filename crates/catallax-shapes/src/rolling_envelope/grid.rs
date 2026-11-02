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
        let mut cache = self.data.net.grid.lock().expect("no panic holds this lock");
        if let Some((key, grid)) = cache.as_ref()
            && key == tolerance
        {
            return Ok(grid.clone());
        }
        let grid = Arc::new(sheet_grid(
            Patch::all(&self.data.carrier, &self.data.net.controls, DEGREE),
            tolerance,
            max_points,
            dec,
        )?);
        *cache = Some((tolerance.clone(), grid.clone()));
        Ok(grid)
    }
}
