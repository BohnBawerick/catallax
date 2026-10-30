## 2026-10-30: terminal network gyroid

Pieces: 1, 5 (catallax-field `lattice_stage`, `printing`; catallax-recipe)
Change: `Operation::Gyroid` now compiles in the field kernel for a network with a constant level and an identity frame. Its leaf is intersected with the region. `lattice_stage` adds `GyroidSpec` and the `TerminalField` trait, and `Stage` becomes `Stage<S = LatticeSpec>`, so `Stage::new(body, GyroidSpec { .. }, settings, dec)` builds a gyroid stage. `Stage::return_before_lattice` returns the stage's own spec. `Recipe::check_lattice_terminal` and `check_no_cad_after_lattice` now also treat a gyroid as terminal, with the same error and refusal codes. `RecipeError::LatticeNotTerminal` keeps its name, code and fields; its message now names a terminal field.
Why: V7 needs one gyroid through the terminal-field stage: unproved display, certified export at the shipped `1/10` mm, and labelled STL. A second operation or spelling would duplicate the recipe definition.
Migration: none for lattice callers. A signature that names `Stage` keeps meaning a lattice stage. An application that adds a gyroid passes a `GyroidSpec`, must treat refusal `gyroid.unsupported` as a deferred function and `gyroid.impossible` as invalid input, and must label the level as dimensionless.

### Supported domain

- `kind` is `network`, `level` is `{"constant": c}`, `frame` is the identity. Any finite rational `c` is valid, including zero and negative values. `period` must be positive.
- `g` is a dimensionless residual. `Query::value` is `Value::Residual`, and `Approximate::true_distance_formula` is false. The leaf's scalar is `(g - c)/L` with `L = 22/a`, so its magnitude bounds distance to the level set from below.
- Point queries are exact. Mass, B-rep evaluation, STEP, drawings and copies by pattern or mirror are not supported. `mass::properties` refuses with `field.print-coverage`, and copies refuse with `gyroid.unsupported`.
- The certified display refuses with the existing code `field.display-lattice`, so the existing adapter fallback to `lattice_stage::preview_plan` applies. The preview stays `Proof::Unproved`.

### Proof and budgets

`Stage::export` certifies a two-sided `1/10` mm bound by the method in `docs/FIELD.md`. `printing::mesh_plan` uses the same certificate for a gyroid at an accuracy of one hundredth of its period or finer, and every other plan keeps its uniform-grid route. Caps and work limits are unchanged: 2,000,000 cells, 1,000,000 vertices and 100M work units. An export of a truly disconnected clipped solid still refuses with `field.print-disconnected`, including a small corner piece. At period 10 and level `1/2`, the box `[-10,10] × [-10,10] × [-5,15]` is such a solid.
