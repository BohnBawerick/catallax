## 2026-10-22: consume temporary dyadic magnitudes

Piece: 1
Change: ball products and midpoint-rounding residuals consume fresh dyadic values when taking absolute values. Midpoint/radius normalization also retains owned, already-precision-sized results instead of copying their integer mantissas. Precision reduction of borrowed balls keeps a borrowed path.
Why: radius formulas and interval normalization copied temporary results before immediately discarding them. An odd normalized mantissa changes whenever significant bits are discarded, so the representable-center case needs no rounded copy or error subtraction.
Migration: public signatures, arithmetic, rounding and work charges are unchanged. A 57,456-coefficient public jet workload retains identical midpoints, radii and precisions; independent GMP/MPFR tests check containment.

