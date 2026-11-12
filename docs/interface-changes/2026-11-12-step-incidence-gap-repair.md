## 2026-11-12: measured STEP incidence and bounded planar derivatives

Pieces: 8, 11.
Change: `measure_incidence` returns an `IncidenceGap` enclosure for a saved source-incidence refusal, with its finite-feature or enclosing-carrier scope.
`ImportRefusal::IncidenceGap` carries the measurement, the missing repair certificate and any later strict solid refusal.
Spline-bearing planar solids with line or spline boundaries may build an exact derivative within the existing 0.01 mm import cap and one tenth of the certified neighbouring feature size.
Attained vertex and retained-curve displacement bounds are recorded in the imported body's distance marks.
Why: exact decimal STEP sources can disagree at their boundaries even when their intended joins are much closer than the imported-geometry cap.
Migration: handle the new serialized `incidence_gap` refusal and preserve the body's marks when consuming imported results.
Only exact repaired geometry can create exact incidence facts.
The measurement is a source-incidence distance, not a whole trimmed-face Hausdorff certificate or an output-error estimate.
Unsupported support fitting, boundary coverage and solid checks remain typed refusals.
