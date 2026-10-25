## 2026-10-23: saved first-hit feature instances

Pieces: application, features
Change: `expand_feature` transports explicit directions for `extrude_to` and `extrude_next`, supplies the mapped authored normal for reflected default directions, and reuses explicit input/origin binding for target and limit references.
Why: a saved reflected capped extrusion must query its reflected target along the original physical sweep direction, not the repaired sketch frame's opposite normal.
Migration: no schema change; save the ordinary version 4 steps returned by expansion.

