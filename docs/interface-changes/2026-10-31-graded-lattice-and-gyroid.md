## 2026-10-31: graded lattice radius and gyroid level

Pieces: 1, 5 (catallax-field `grade_range`, `lattice_stage`, `printing`, `mass`; catallax-features `evaluate`, `lattice`; catallax-cli `catallax-stand-in/1`)
Change: One clamped linear grade along a part axis now drives both the cubic lattice radius and the network gyroid level through the terminal stage, certified print, field moments and the explicit mass stand-in. `Gyroid` compiles that linear level; `LatticeSpec` and `GyroidSpec` validation names the graded refusals; the B-rep kernel carries a gyroid by its region envelope, as it already carries a lattice.
Why: V7d needs graded density on both field families: upstream edits, native proof, labelled export and required mass evidence, without a second position function or field family.
Migration: An application may pass a linear grade in `LatticeSpec::radius` or `GyroidSpec::level`. It must label the radius in millimetres and the level as dimensionless, and treat the new grade codes below as coverage refusals of valid input. Callers matching `features.unsupported-operation` for a B-rep gyroid now receive the field-envelope codes listed below.

### Supported domain

- A grade is `{"constant": c}` or `{"linear": {"from": f, "to": g, "start": c0, "end": c1}}`, read in part coordinates. Equal endpoints compile to the constant they equal and take every constant route, byte and limit.
- `Stage::new` certifies linear grades whose anchors differ along exactly one part axis by a nonzero distance. A lattice radius needs endpoints in `(0, size/4]` and slope `K = |c1 - c0|/|g - f| <= 1/2`. A gyroid level needs endpoints in `[-1/4, 1/4]` and `K <= 1/10` per millimetre.
- New stage refusals: `lattice.grade-function-unsupported` and `gyroid.grade-function-unsupported` for anchors off one axis; `lattice.grade-range-unsupported` and `gyroid.grade-range-unsupported` for an endpoint or slope outside the range. Existing codes stay: `lattice.impossible` and `gyroid.impossible` (including coincident anchors), `lattice.unsupported-grade` and `gyroid.unsupported` (distance grades, sheets, frames), `lattice.settings`.
- Lattice point queries keep every valid linear grade, including oblique anchors. A gyroid leaf compiles only the certified subset and refuses others with `gyroid.unsupported`.
- Leaf scalars are labelled residuals: `(β - ρ)/(1 + K)` and `(g - c)/(22/a + K)`. Neither is a signed distance.

### Public additions

- `catallax_field::grade_range(grade, lo, hi, dec)` returns the exact clamped range of a constant or linear grade over a closed rational box. Refusals: `field.grade-box`, `field.grade-unsupported`, `field.grade-impossible`, and the decider's work and deadline outcomes.
- `catallax_field::printing::classify_box(plan, lo, hi, dec)` returns the print certificate's interval `BoxState`; `printing::BoxState` is now public.
- `catallax_features::lattice::validate_gyroid(period, level)` checks the gyroid parameters both kernels read.

### Behaviour

- `printing::mesh_plan` and `Stage::export` use the adaptive certified candidate mesh for every axis-graded lattice or gyroid, at every minimum detail. A graded gyroid keeps a coarse candidate pitch of at most one thirty-second of its period. Constant fields keep their routes and bytes.
- `mass::properties` measures an intersection of axis blocks or axis prisms with one gyroid, constant or graded, or one axis-graded cubic lattice, by monotone-graph integration. A constant gyroid previously refused with `field.print-coverage`. Constant lattices keep analytic columns and admitted-box integration.
- `catallax_features::evaluate` records a gyroid in `Built::stand_ins`, so a plain recipe STEP or IGES refuses with `field.stand-in-export` and a B-rep print with `print.brep-field-definition`, as for a lattice. The `field.stand-in-export` detail now names a field envelope.
- `catallax-stand-in/1` compares the exact surrogate and the STEP-rounded exchange with the field enclosures before it writes `answer.step`. Codes are unchanged: `stand-in.inertia-not-matched`, `stand-in.properties-not-matched`, `stand-in.exchange-not-enclosed`.
- The sampled preview of an axis cubic lattice uses its smallest radius as its pitch, graded or constant.
