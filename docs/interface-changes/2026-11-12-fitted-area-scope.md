## 2026-11-12: a fitted blend's area says what it encloses

Pieces: 6 (measurement)
Change: `volume_area_centroid` and `area` on a body with rolling-envelope, rolling-saddle or sphere-band faces add a `measure.area-of-definition` mark (`catallax_measure::AREA_OF_DEFINITION`, `BoundKind::Distance`). Its bound is the worst fitted carrier's measured distance to its definition. The mark states that the area enclosure contains the exact definition's area only, and that the carriers' own area is not bounded. Volume and centroid are unchanged: their widths already enclose the carriers. `blend_definition_properties` asks for the definition alone and carries no such mark.
Why: PLAN M12 asks for each quantity of a body with fitted faces to carry its own stated bound. The area was the definition's area under the represented scope, without saying so, so its enclosure could be read as covering the fitted carriers.
Migration: none. Readers of marks see one more mark on fitted-blend area answers; values are unchanged.
