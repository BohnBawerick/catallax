## 2026-11-18: numerical embedding moves to the features crate

Pieces: 7, features, application, verification CLI
Change: `catallax_step::numerical::build_from_exact_model`, `definition_from_exact_model` and `embeds` are gone. `catallax_features::numerical_embedding::embeds` answers as before, and `definition_from_exact_model` there gives the same embedding of an exact model as a `NumericalDefinition` with its native UV trims. The caller admits it with `catallax_topology::numerical::build_with_native_trims` and writes it with `catallax_step::numerical::build`. This corrects the entry point named in [the numerical STEP embedding fragment](2026-11-16-numerical-step-embedding.md); the embedded family, tolerance and manifest described there stand.
Why: the rolling-ball fillet embeds the body it rounds, and features cannot depend on the STEP writer. Embedding now has no production dependency on it.
Migration: replace a `build_from_exact_model` call with the three calls above, in that order, inside one work request.
