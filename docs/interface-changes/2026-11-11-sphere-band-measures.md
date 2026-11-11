## 2026-11-11: bodies with sphere bands measure volume, area and centroid

Pieces: 6
Change: `volume`, `area`, `volume_area_centroid` and `blend_definition_properties` answer for a fillet whose ring bands lie between a sphere and a cylinder that passes through it off its centre. They refused `measure.unsupported-surface` ("a band between a sphere and a cylinder has no certified change integral"). The result is the recorded unfilleted base plus each band's pocket and area change. Volume and centroid also hold the fitted carrier, through the same tube bound as rolling-saddle bands; area is the definition's. `mass_properties`, face and edge measures on the band and its rails, `bounding_box` and meshing still refuse.
Why: M12 asks for area, volume and centre of gravity for bodies with fitted faces, each with its own stated bound. Sphere bands were the remaining ring band without a change integral.
Migration: none. A refusal becomes an answer; no proven answer changes.
