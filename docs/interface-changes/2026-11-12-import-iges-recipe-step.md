## 2026-11-12: the `import_iges` recipe step

Pieces: recipe, 8, features, field, compare, verification exchange

Change: recipes gain `import_iges`, with the parameters of `import_step`
(`file`, `pose`, `decisions`). `catallax_recipe::Operation::ImportIges` carries
them; its faces take the existing `imported` role, keyed by IGES directory-entry
sequence numbers such as `D17`. `catallax_import` gains `SourceFormat`,
`apply_source` and `propose_source`; `apply` and `propose` are their STEP
forms and are unchanged. The features evaluator, the field compiler, the
comparison's imported-part checks and `import_edit::replace_step` accept the
new step the way they accept `import_step`. A repaired IGES import's derivative
is a STEP file. The verification exchange's `catallax-import-recipe/1` wrapper
carries IGES text the way it carries STEP text. See RECIPE.md section 4.18.2.

Why: PLAN M24 asks for IGES import usable by the application. Builds this step
admits were refused or impossible before, so no existing answer changes.

Coverage, raw evidence retained privately: a Catallax-written and an
OpenCASCADE-written 2 by 3 by 5 block build placed, with volume exactly 30 and
`D` keys; a shell that opens the top face by its imported key builds; the field
compiles the import; replacing the Catallax file by the OpenCASCADE file with an
explicit face map rewrites the shell's reference and rebuilds the same volume,
and an unmapped key fails with `import.replace-unmapped`.

Migration: exhaustive matches on `Operation` add the `ImportIges` arm. Callers
that need IGES proposals or repairs use `propose_source` and `apply_source`.
No recipe format version, saved answer, tolerance or work allowance changes.
