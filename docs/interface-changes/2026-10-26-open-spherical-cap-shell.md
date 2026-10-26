## 2026-10-26: Planar opening below one spherical cap

Pieces: recipe features, independent field, application, verification consumers.
Change: A convex coaxial cylinder/sphere cell with a complete planar base selected open now produces an exact distance-defined shell when the base is on the sphere-centre plane or towards the cap before the join. This bounded family upgrades `distance.cylinder-domain` to a shell answer in both kernels; other open circular families retain coverage refusals.
Why: A whole-domain nearest-projection proof supplies witnesses on the finite kept cylinder and spherical cap, including their join and the open rim, before constructing inward supports.
Migration: No signature, recipe, format, defaults or consumer dispatch changes. Callers consume the existing `CircularSkin::cavity` result. See [bounded capability, proof and controls](../OPEN-SPHERICAL-CAP-SHELL.md); this does not authorise unlisted release-replay upgrades or exporter changes.
