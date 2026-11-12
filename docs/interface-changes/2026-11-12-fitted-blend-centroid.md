## 2026-11-12: the centroid of a fitted blend answers with its stated bound

Pieces: 6 (measurement)
Change: `catallax_measure::centroid` on a body with rolling-envelope, rolling-saddle or sphere-band faces answers through the fitted-blend junction, with the same enclosure `volume_area_centroid` gives: it contains the exact definition's centre and the fitted carriers' moment slack. Before, it refused `measure.unsupported-surface` ("surface kind has no measurement definition"). `catallax_measure::inertia` on such a body now refuses before any work with the detail `mass_properties` already uses ("second moments of a rolling-envelope blend are not covered").
Why: PLAN M12 asks for the centre of gravity of bodies with fitted faces with its own stated bound. One entry point stated it and the other refused.
Migration: none. A refusal becomes an answer; the inertia refusal keeps its code.
