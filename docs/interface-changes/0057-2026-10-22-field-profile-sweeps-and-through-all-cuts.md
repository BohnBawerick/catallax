## 2026-10-22: field profile sweeps and through-all cuts

Pieces: field, features, recipe
Change: `catallax_field::Plan::compile` accepts sketch profiles, signed `extrude`, two-sided `extrude_through` and rational-frame full/partial `revolve`. Sketch validation uses `catallax_features::Profile::new`; feature supports and trims are evaluated without calling a B-rep builder. Source face candidates retain segment aliases, split domains, reversed cut orientations and 180-degree merged-cap names.
Why: both recipe evaluators must build the same feature definitions and report their creating origins.
Migration: no Rust signature, recipe format or command-line protocol changes. Sweep values use `Value::Residual` and `Approximate::true_distance_formula = false`; unsigned distance upper bounds remain unbounded off the certified boundary. Callers must not treat a sweep residual as true Euclidean distance. The specific unfinished `features.unsupported-*` codes from RECIPE.md section 7.14 now apply to field compilation as well.

