# M24: multiple direct STEP solids and explicit replacement

Baseline: `371774eaf3ffb8321ad54311c7596a25c8a0cce3`.

This is a bounded first cut of PLAN M24 and RECIPE 4.18.1. No recipe format,
face-origin grammar, tolerance, proof allowance or saved-answer inventory changes.

## Multiple solids

The recipe B-rep evaluator admits every closed direct `MANIFOLD_SOLID_BREP`
when certified enclosing boxes are pairwise strictly separated. It collects the
existing lumps, in source order, without a boolean, face renaming or material
duplication. The recipe pose applies to the complete body. `Built::import_solids`
retains each solid's source entity key and its placed lump membership, including
across an unchanged incremental prefix. As before, `ADVANCED_FACE` entity keys
are the imported face identities, fixed by the file hash.

No solid can be omitted because its repair was declined. No partially built
recipe is published. Contact, overlap and disjoint shapes whose enclosing boxes
overlap refuse `import.multisolid-separation`; this is a coverage refusal, not a
clash verdict. Multiple-solid mapped occurrences and transformed representation
relationships refuse before construction because their occurrence placements
are not supported. Assemblies and void-bearing STEP entities retain the reader's
existing refusals. This does not extend STEP product-label or assembly semantics.

The field shares source admission. Its existing planar component route can admit
multiple planar solids. The public two-block fixture builds a field display, but
its requested volume comparison refuses `field.volume-definition`; this is not
a passing two-kernel volume comparison. Other multiple-solid analytic domains
may still refuse.
This change adds no field, mass or FEA algorithms.

## Replacement transaction

`catallax_features::import_edit::replace_step` accepts an immutable recipe,
`StepReplacement`, hash-indexed immutable source bytes and a caller decider. The
proposal names the import step, expected old hash, new file name/hash, new saved
decisions and an explicit old-face-to-new-face relation. It returns a fully
rebuilt `Built`; `Built::recipe()` is the candidate revision to save atomically.
The caller retains the old recipe and old bytes for undo. Failure never mutates
either revision or source bytes, and no file is opened using a recipe label.

Both source hashes and all supplied map keys are checked. Every referenced old
key needs an explicit singleton target, even if the new file has the same
numeric entity key. Nested `of`, edge and vertex origins, `touching` and ordinary
face selections are visited through the typed recipe's serialized schema.
Re-parsing restores canonical names and rejects collapsed vertex names; collapsed
edge selectors are rejected explicitly. A complete rebuild checks downstream
reference resolution and supported feature geometry before returning a candidate.

Splits, deleted targets and step-wide origin wildcards return
`import.replace-selection`; this entry does not expand lists or choose one member
of an ambiguous relation. Many-to-one maps are allowed only when references
remain expressible and the rebuilt recipe succeeds. Geometry is never used to
guess correspondence. The proposal is editorial intent, not evidence that old
and new faces are geometrically equivalent. Existing unfiltered selections keep
their saved meaning; they contain no old entity key to rename.

Old decisions are never copied. The supplied new decisions pass the existing
new-source question-id and repair validation. Other failures distinguish stale
edits, missing sources, hash mismatches, invalid keys, unmapped references,
invalid structured references and rebuild failures with `import.replace-*`
codes. A rebuild failure retains the underlying error in its detail.

The API creates no UI undo stack and writes no files itself. Applications must
commit the returned recipe/build pair together. Unsupported downstream geometry
remains a typed refusal, including unsupported imported spline certification.

## Validation scope

Rust tests use public synthetic STEP boxes and the existing public application
fixture. They check separate membership, exact volume, placement and face keys;
overlap/contact/duplicate-material refusal; mapped occurrence refusal; explicit
renumbering and rebuild; stale/missing/hash/unmapped/split/decision failures;
incompatible edge references; nested names and merged-edge refusal; rollback;
and identical recipe/STEP bytes at 1, 2 and 4 workers.

These checks do not qualify all of M24. No C++/OCCT oracle or private saved-answer
corpus is run locally. Hosted exact-head integration and preservation gates remain
required; local Rust evidence is not an overall PASS.
