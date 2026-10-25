# Convex polygon stock-quadric first-hit caps

## Scope and provenance

Base/main inspected before edits: `357be43e3ea5ba00df9fd449ea2afca09260ae24`.
Branch: `cloud/convex-polygon-quadric-first-hit`. The delivery head and exact-head
qualification are recorded in the draft PR and external proof identity. No
private corpus or held-out input is used. The new public fixtures derive from
the existing `curved-limits/{sphere,cylinder,cone}-entry-next.json` definitions,
replacing only the authored profile with rational triangle, pentagon and rotated
square loops. Directions, both limit modes, entry/exit starts and radius edits
are explicit test variants.

The initial checkout was clean. No `AGENTS.md` was present in the checkout or its
ancestor paths. Repository verification-skill, feature-map, protocol and fragment
instructions were read. PRs 93, 101, 103, 107, 109, 112 and P3/108 were inspected
before edits. The pre-edit reservation and commands are retained under
`/workspace/convex-quadric-cap-proof` in the saved environment.

Owned production changes are `AuthoredCap` compilation/side ids and
`CurvedCapped` domain queries in `end_condition/curved.rs`; construction/build
dispatch and shared prism radius calculation in `curved/build.rs`; the new
`curved/polygon.rs` validator/supports/column topology; and `CurvedField` side
membership/origins in field `end_condition.rs`. Existing Plan hooks, projected
quadratic certification, stock-scene compilation, annular torus path, Boolean
engine, display, revolution and distance-feature production files are untouched.

## Defect ledger

| Defect | Before repair | Repair and evidence |
| --- | --- | --- |
| Four-line rectangle gate | Fresh B-rep and field CLI runs reject triangle, pentagon and rotated diamond with `features.unsupported-curved-limit`; no STEP. | Retain the whole rectangular enclosure proof; construct and query the intersection of exact polygon side halfplanes. Both kernels consume the same finite analytic column definition. |
| Adjacent collinear authored sides | The added control produced `boolean.merge-ambiguous` because adjacent coplanar prism faces did not close after Boolean merging. | Merge straight topology sides using exact rational collinearity and union their authored origins. Analytic cells retain equivalent halfplanes. The control now passes without Boolean-engine edits. |
| Enclosure too broad | A diamond of circumradius 3/2 has sphere hits everywhere, but its enclosing square crosses the radius-2 silhouette. | Retain `features.unsupported-curved-limit`; the test independently establishes the polygon lies strictly within the silhouette. This is unfinished coverage, not an impossible missing cap. |
| Native spherical STEP volume consistency | OpenCASCADE imports valid solids, but the imported entry/exit volume sum differs from the exact polygon integral. | Still open. Gauss-Kronrod integration requested at 1e-12 reports discrepancies of about 6.869e-8 (pentagon), 1.334e-7 (rotated square) and 1.274e-6 mm³ (triangle). The new 1e-7 identity assertion remains enabled and fails for the latter two profiles. The cause is not isolated; exact Catallax volume checks do not clear this independent gate. No export/import production repair is attempted in this lane. |

Nonconvex, self-crossing, degenerate/repeated-vertex, trim-crossing and exhausted
requests retain typed refusals. Tangency, changing branches and general target
graphs remain outside this lane. Construction never drops higher-degree ray
terms, increases a cell cap or changes a tolerance or request allowance.

## Focused checks

`cargo test --release -p catallax-features --test polygon_end_conditions -j 2 -- --test-threads=2`
passes ten tests on the implementation before final qualification:

- 36 profile/stock/mode/branch combinations build in B-rep and field, retain
  every authored side origin and distinguish exact polygons from their bounding
  rectangles. Rational and irrational off-axis roots are checked independently.
- An independently authored rational inequality oracle checks 1,764 field
  points, including polygon boundaries and points beyond the far target root.
- Clockwise profiles, proper world placement and positive/negative raw direction
  scaling preserve entry/exit roots and physical side origins.
- All nine radius edits resolve original face origins and physical-edge names.
- All 36 consumer combinations produce B-rep display and STEP with unchanged
  options, no marks and one solid.
- Requested default-accuracy B-rep volumes enclose an independent exact integral
  for the cylinder triangle. All nine profile/stock entry-plus-exit integrals
  enclose exactly five times the authored polygon area.
- The public native STEP-import regression covers both modes and branches;
  independent import validity/volume checks remain distinct from certified
  Catallax property enclosures.

Initial local Rust compile/borrow errors and the red collinear control are
retained in external logs; they are not passing qualification. Full-suite,
real-CLI, exact-head CI and consumer limits are reported separately in the draft
PR after execution. No private final477 run is performed here; the parent owns
that gate. No matched timing experiment or speedup claim is made.

Native field display was attempted on all nine new profiles with default
options and the unchanged 100M verifier allowance. Each independently compiled
definition reaches `field.display-unsupported`, with detail `display chart needs
a rational unit frame`. This separate held display path is not a prerequisite
for the construction/query slice and is not repaired here. Full default field
mass and final-head checks retain their actual outcomes in external evidence and
the PR; no success is inferred from B-rep properties.

## Reproduction commands and evidence

Use the installed workspace toolchain with
`PATH=/workspace/.cargo/bin:$PATH`, `CARGO_HOME=/workspace/.cargo`, and
`RUSTUP_HOME=/workspace/.rustup` in this environment. The initial release build
command was `nice -n 10 cargo build --release -p catallax-cli --bins -j 2`.
The frozen base verifier SHA-256 is
`25f91d05a6f591042bc47e759d4ed062b10e7b3ced45be2b017c3904c4e40208`;
doctor verified `catallax-request/1` and the embedded base commit.

```sh
.claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary /workspace/convex-quadric-cap-proof/baseline-verify --json
.claude/skills/verify-catallax/scripts/control-catallax.py run /workspace/convex-quadric-cap-proof/reproduction-inputs --source recipe --binary /workspace/convex-quadric-cap-proof/baseline-verify --out /workspace/convex-quadric-cap-proof/baseline-reproduction --json
.claude/skills/verify-catallax/scripts/control-catallax.py run /workspace/convex-quadric-cap-proof/reproduction-inputs --source field --binary /workspace/convex-quadric-cap-proof/baseline-verify --out /workspace/convex-quadric-cap-proof/baseline-field-reproduction --json
```

External command records preserve argument arrays, UTC start/end, exit code,
base/head, branch, source provenance and diagnostic wall time. Each CLI run
retains config, ordered inputs, request/response, operation logs, manifests and
actual STEP/OBJ side effects. Fresh retained-domain recipe, field and comparison
runs use the same nine public rectangle/torus inputs and unchanged 100M verifier
allowance. Comparisons use the existing driver exclusions and strict artifacts.
Agreement is sampled evidence, work is not latency, and native field display or
mass refusals are separate from successful analytic construction/query checks.
