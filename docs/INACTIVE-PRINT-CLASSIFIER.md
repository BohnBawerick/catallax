# Inactive print and mass classifier definitions

The classifier prepares only leaves whose existing `leaf_registers` entry is
present. An inactive entry occupies its original definition slot without
coefficient enclosure or coverage dispatch. Every executed leaf instruction
still indexes the same original leaf and its prepared definition. Compaction,
the root Boolean tape and finite named-face domains are unchanged. This changes
preparation, not the represented point set or the supported active leaf kinds.

## Reproductions and independent controls

The frozen public base is `866c9a30c4730ebe906cba352f61bf350b251652`.
Before the production edit, the dedicated mass and print regressions both
failed with `field.print-coverage` for retained but inactive triangular extrusion
definitions. The active triangular extrusion control retained that refusal.
The frozen verifier SHA-256 is
`660956d7076d6229f367398f6271d31762df55fa5fe587bb56cddf59566c4742`.
The real print adapter answered the compact box and refused the otherwise
equivalent box with an unused Feature definition. Source inspection alone is
not counted as reproduction.

Dedicated tests put distinct inactive Feature definitions before, between and
after live leaves. Two overlapping boxes of side `2s`, displaced by `s` along x,
have independently derived union volume `12s³` and intersection/subtraction
volumes `4s³`, with `s=1/128`. Ordinary `mass::properties` is called at default
Accuracy and caps. Compact and retained plans must give equal exact volume,
centroid, inertia, classified cells and work. Exact interior, exterior and tie
points check the root Boolean semantics; finite named-face queries include
positive membership and an off-domain negative.

The printable single-box control uses the default `0.01` mm tolerance and print
caps. Its retained and compact versions must produce actual nonempty meshes,
identical STL/OBJ bytes, cell counts, work and certified deviation. Boolean print
controls preserve the complete compact decision, including any remaining
refinement refusal. They do not count such refusals as printable output.

## Public thicken fixture and remaining gate

The explicit public fixture handoff is
[`collinear-neighbours.json` at PR119's frozen head](https://github.com/BohnBawerick/catallax/blob/1c187dec82a078a83d8dcc8947fe6612deae8804/crates/catallax-features/tests/fixtures/boolean-face-thicken/collinear-neighbours.json).
Its finite thicken result is `[0,4] × [0,4] × [3,4]`, volume 16. The original
ordinary library mass probe has no operation-budget wrapper. Its default
volume/centroid/inertia Accuracy is relative `1e-12`, absolute `1e-30`, and its
cell cap is 2,000,000. Do not describe that original probe as a 100M request.

Mandatory printing uses tolerance `1/100` mm, 2,000,000 cells and 1,000,000
vertices. The equivalent compact box already returns `field.print-budget`
(`print interval refinement cell limit reached`) on the frozen main verifier.
It consumes 2,000,000 work under a separate 100,000,000 verification allowance,
without exhausting that work allowance. No STL or OBJ is published. This is a
real consumer blocker beyond preparation; the required print smoke remains
open. A smaller synthetic success does not satisfy this required smoke.

PR119 landed as `68f80436a4a9804c5d6a15435c4dfecc6119d758` before integration
and migration of its single obsolete mass-refusal assertion. The original
unwrapped library call now requires exact volume 16. A dedicated independent
oracle also requires centroid `(2,2,7/2)` and centroidal inertia
`diag(68/3,68/3,128/3)`. Its publication implementation and all other assertions
are preserved. Frozen-head integration probes remain separate from final-main
evidence. The draft PR records exact final source/runtime identities, commands
and terminal checks.

The required handed-off print success test was introduced by this change at
`172636abda9b108587726c09db2d78ab1504f448`; it did not exist on main. A scope
split moves that added test to a separate follow-up without changing its enabled
`#[test]`, default controls or success assertions. The [original enabled test](https://github.com/BohnBawerick/catallax/blob/172636abda9b108587726c09db2d78ab1504f448/crates/catallax-field/tests/inactive_print_classifier.rs#L281)
has SHA-256 `92a8509bc4aa71bd52ac222908abfd3ce1308b5c6c8105771e7dacb11455f1a8`.
The preserved restoration patch has SHA-256
`76bc91febd45f6a728dc63b4f0d8faee98ca71f4c20d0603fad5cd25204574ed`.
No existing main test is removed or ignored, and the new exact mass oracle and
active unsupported refusal control remain in this classifier change. The print
limit remains an open consumer defect; moving its test does not establish print
acceptance. Neither an expected-refusal assertion nor a smaller print replaces
the unchanged follow-up success gate.

## Defect and attempt ledger

| Unit or attempt | Observed result | Disposition |
| --- | --- | --- |
| Inactive Feature preparation | Original dedicated mass and print regressions fail with `field.print-coverage`; active Feature refusal passes | Skip preparation using existing liveness, with indexed empty slots |
| First synthetic Boolean print expectation | A default-cap print reaches `field.print-budget` after preparation is repaired | Preserve compact decision and keep a separate enabled actual-output test; no cap increase |
| Exact-size compact thicken equivalent | Real default-cap CLI print reaches `field.print-budget` | Required thicken print gate remains open; report the independent limit |
| PR119 mass-refusal assertion | Dependency landing verified at `68f80436`; original call now answers exact volume 16 | Migrate only that assertion; add independent exact moments |
| Required retained thicken print | Both ordinary library and actual CLI reach `field.print-budget` after preparation repair | Preserve the unchanged enabled success test and restoration patch for separate follow-up qualification; no print acceptance claim |

Verification uses the checkout's `verify-catallax` skill, public synthetic
inputs, preserved baseline executable and separate evidence directories outside
the repository. No private proof inputs are read. The 1M scoreboard, 20M
application and 100M verification allowances remain distinct. Default mass and
printing caps and all tolerances are unchanged. Shared-cloud timings are
diagnostic; no wall-time speedup or complete target acceptance is claimed.
