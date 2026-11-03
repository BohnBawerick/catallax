## 2026-11-03: volume, area and centroid of a body with a loop fillet

Pieces: measure
Change: `volume`, `area`, `centroid` and `volume_area_centroid` now answer for a body whose fillet rounds the loop of two crossing cylinders, where they refused with `measure.unsupported-surface`. As on a through ring, the volume and centroid enclose both the exact definition and the fitted carrier, and the area is the exact definition's.
Why: this is the second solid piece of the loop slice. The through ring's integrals stand columns along the narrow axis on the wide cylinder, which a loop crosses at its outline, so a loop's change takes the divergence theorem over the band and two wall strips instead.
Migration: code that treated a loop body's volume as a refusal sees an enclosure. `mass_properties`, face measures, STEP and the certificate still refuse on a loop.
