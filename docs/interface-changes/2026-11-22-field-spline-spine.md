## 2026-11-22: field spline spine certificate

Pieces: catallax-field
Change: `spline_spine::Spine::prove` certifies the centre path of a constant-radius ball against two whole, single-span rational patches. The result carries span models, error bounds and both border roots. Invalid definitions, unproved paths and unattainable tube widths have distinct refusal codes.
Why: the field needs an independent path for spline fillet queries and two-kernel checks.
Migration: existing callers need no changes. New clients supply the radius, contact sides, the continuation parameter, a seed station and part tolerance through `SpineDefinition`. They must establish that the certified component belongs to the selected recipe edge before using it to change a solid.
