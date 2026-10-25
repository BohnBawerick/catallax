## 2026-10-22: avoid rational differences in diamond-chart comparisons

Piece: 5
Change: cylindrical projection compares each diamond phase against shifted anchor thresholds instead of normalizing a rational difference. Corner bounds are computed once, and four phases use a fixed array.
Why: chart subdivision repeatedly normalized fractions only to compare their order.
Migration: phase quotients, corner order, strict wrap cutoffs and outward hull rounding are unchanged. A before/after smoke retains identical results for 2,916 scaled boxes at 8, 64, 128 and 256 bits. Permanent checks independently bound all four quadrants, a periodic seam and an axis, and reject a box containing the projection pole.

