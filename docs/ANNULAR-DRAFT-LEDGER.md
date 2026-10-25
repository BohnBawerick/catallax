# Bounded axial annular draft

This lane constructs one coaxial annular cell with two perpendicular finite
planar caps. Either cylindrical wall, or both walls, may move. The inner wall
keeps reversed material sense. The hole is excluded throughout construction,
queries, mass integration, B-rep meshing and STEP export.

## Identity and source provenance

- Inspected and frozen base/main: `357be43e3ea5ba00df9fd449ea2afca09260ae24`.
- Branch: `cloud/axial-annular-draft-20261025`. The final PR description and
  external proof identity record name the exact delivery head and CI runs.
- Original worktree was clean. No `AGENTS.md` existed in the checkout, its
  parent directories or the repository file inventory. Read the repository
  verification skill, feature map, answer workflow and interface fragment rules.
- Baseline release verifier SHA-256:
  `25f91d05a6f591042bc47e759d4ed062b10e7b3ced45be2b017c3904c4e40208`.
- All new inputs are public synthetic recipes under
  `crates/catallax-features/tests/fixtures/draft-annular/`. The source is a
  two-circle sketch, radii 8 and 2 mm, extruded 2 mm. There is no private corpus,
  held-out fixture or private proof-input access in this lane.
- External commands, identities, logs and CLI artifacts are retained in
  `/workspace/annular-draft-proof/`. Ownership was recorded before production
  edits in `ownership-before-edits.txt`; polling status is in `CHECKPOINT.md`.

## Ownership and construction proof

The new feature `draft_annular.rs` owns `AnnularWall` and `AnnularCell`.
`draft.rs` registers its module. `draft_build.rs` owns the bounded recognizer,
annular assembler and `Cell` dispatch for membership, anchoring, bounds and
origins. The recognizer requires four faces, each with two complete circular
ring boundaries, opposite wall material senses and two opposing cap planes.
No untrimmed support, convex hull or largest-component selection is substituted.

The field's new `draft_annular.rs` recognizes the four carriers from an
independently compiled finite feature, checks its annular and hole membership,
and evaluates distances to the four finite boundary segments in the meridian
plane. Its `draft.rs` owns direct annular cells, constraints, bounds, transforms,
fast caches and publication. Small unheld hooks are:

- `geometry.rs::Constraint::reversed` negates both exact and cached polynomial
  terms for the inner wall; `Leaf::draft_witness` uses a strict annular point.
- `samples.rs::Region::sample_candidates` samples both walls and cap annuli.
- `draft_source.rs::Source::atom` retains the curved-source coverage refusal
  in its exhaustive match.

The source and replacement radii and their difference are affine in axial
height. Strict positivity at both finite ends certifies the entire interval.
For the inner wall, the authored positive angle increases radius in the pull
direction. This is the reversed source material normal, not a sign inferred
from the cylinder's axis direction. The caps retain their source origins and
each has an outer and inner ring. Default retention uses positive-area neutral
intersection; explicit picks require strict annular material interior.

Active production ownership for PR93, PR101, PR103, PR107, PR109, PR112 and P3
was inspected before edits. Their files, native display, Boolean, numeric,
shapes storage, shared `features.rs` and `lib.rs` are unchanged. No second draft
arrangement/dispatch worker was used. The aggregate interface log is unchanged;
the additive fragment is
[`2026-10-25-axial-annular-draft.md`](interface-changes/2026-10-25-axial-annular-draft.md).

## Defect ledger

1. Current-main recipe and field CLI runs independently returned
   `draft.unsupported-arrangement` on the public `both.json` geometry. The B-rep
   guard accepted only one cylinder; the field's single-loop witness refused
   the annular source. Neither baseline produced an answer STEP or display.
2. An uncommitted initial candidate represented the field annulus as subtraction
   of two existing round cells. Five of six feature tests and all three then
   written field tests passed. An expanding outer wall at its top cap corner
   remained `Indeterminate`: the general Boolean witness did not settle this
   boundary. This is a failed predecessor, not a passing head. Replaced it with
   the direct annular field cell. No Boolean witness or shared hook changed.
3. Two local test-compilation attempts exposed ordinary new-test type errors
   (`Real` addition ownership and converting `Rational` test coordinates to
   binary64). They were corrected; they are not geometric refusals or passing
   runs. The later failed-test compilation log is retained externally.
4. Workspace Clippy initially stopped while building its GMP test dependency,
   before linting, because `m4` was missing. Restored Debian GNU m4 1.4.19 under
   the external proof directory, without a system or repository change. The
   failed setup log is retained separately from the restored run.
5. The direct-cell focused suite passes 12 tests: eight feature tests and four
   field tests. Subsequent final-head checks and any failures are recorded in
   the PR and external logs, without relabelling predecessor results.

## Independent oracles and retained behavior

For `both.json`, with `0 <= z <= 2`, the exact radii are `R=8-z` and `r=2+z`.
Integrating annular disks independently gives volume `80*pi`, area
`(80+40*sqrt(2))*pi`, centroid `(0,0,5/6)`, centroidal inertia
`Izz=2400*pi`, `Ixx=Iyy=(11020/9)*pi`, and zero off-diagonals. The requested
B-rep mass enclosures contain these values at unchanged default accuracy.
Single-wall volumes are `272*pi/3` and `328*pi/3`. Signed 30-degree checks use
`(60 -/+ 10/sqrt(3))*pi` on the one-mm cell. Rational exact section controls
cover both caps and walls; posed, scaled, reflected and scalar-edited cases
preserve material and names.

Field true-distance controls include the hole axis, an axial point below the
annulus (`sqrt(5)`), a point above whose closest point is on the finite inner
wall (`sqrt(25/2)`), and a strict material point whose closest cap is 1 mm away.
Inner-wall normals point into the hole. Cap queries reject hole and outer
exterior points. B-rep mesh checks require all four origin-bearing faces,
paired opposite edge uses, outward winding and cap triangles that preserve
the hole within the unchanged 0.005-mm requested chord bound.

Retention tests include both cap anchors, an off-body neutral plane, explicit
empty picks, duplicate interior picks, hole/boundary pick errors, zero-angle
identity and a neutral plane tangent along a generator (no area). Edits rebuild
radius laws and fail rather than retarget a saved pick. Existing draft tests
provide retained-domain coverage for the earlier planar and single-wall paths.

## Exact commands and remaining gates

Shell prerequisites were restored by setting `CARGO_HOME=/workspace/.cargo`,
`RUSTUP_HOME=/workspace/.rustup` and prepending `/workspace/.cargo/bin` to `PATH`.
No account, credential or protected setting changed. The initial missing-Cargo
command did not execute a build. Commands run from the repository root include:

```sh
nice -n 10 cargo build --release -p catallax-cli --bins -j 2
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py doctor \
  --binary /workspace/annular-draft-proof/base-bin/catallax-verify --json
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run \
  /workspace/annular-draft-proof/annulus.json --source recipe \
  --binary /workspace/annular-draft-proof/base-bin/catallax-verify \
  --out /workspace/annular-draft-proof/base-recipe --json
# Same frozen input and binary, --source field, --out .../base-field.
nice -n 10 cargo test --release -p catallax-features -p catallax-field \
  --test draft_annular -j 2 -- --test-threads=2
nice -n 10 cargo test --release -p catallax-field --test drafts -j 2 -- --test-threads=2
nice -n 10 cargo clippy --workspace --all-targets -j 2 -- -D warnings
```

Applicable completion checks are workspace release tests and doc tests, format,
Clippy, the edge-hunt/release-policy/verification-driver suites, generated page
check, public CLI comparison and preserved-domain replay, native STEP imports
with CI-pinned `cadquery-ocp==8.0.1.1.0`, and automatic final-head CI. Exact
commands, exits, counts and head identity belong to the retained external
command log and factual PR description; listing a gate here does not assert it
has passed.

Crossings, touching walls, axis folds, noncoaxial sources, oblique nonzero-angle
neutral curves, more general curved arrangements and ambiguous annular nearest
references remain honest coverage refusals. The independent field recognizer
currently requires one finite source feature; general Boolean curved-source
extraction and redrafting conical walls retain their existing refusals. Native
field display and field mass integration are separate consumers and are not
implied by the B-rep display/property checks. No allowances, tolerances, cell
caps, comparison coverage or typed fail-closed policy changed. No speedup or
latency qualification is claimed. Public correctness and green CI do not
replace full acceptance: final477, application interaction, readiness and merge
remain parent-owned gates. This lane does not run final477.
