# N02 P2 endpoint-sharing repair

Repair base: `fec7e989c4ae0c2a9b2f96ed1835e036eaa0a7bc`, verified as the local,
remote and open PR #605 head before editing. The branch remains
`dot/n02-numerical-trim-sewing`; base target remains main. The independent review
is [review 5907293492](https://github.com/BohnBawerick/catallax/pull/605#pullrequestreview-5907293492).
This author implements only that P2, not an independent re-review.

The executable source witness matches the review prediction. The square tube is
scaled by 0.001 in its vertex positions, line carriers, plane origins and UV
trims. Two disjoint pairs of corresponding outer/inner top corners are joined.
Each squared displacement is exactly 0.000002 mm², giving a gap about
0.001414 mm under the unchanged 0.01 mm local tolerance. The regression confirms
that no edge collapses. The clean scaled shell has Euler characteristic 0.

Before the production repair, the new refusal regression was executed against the
old implementation and failed as intended (test exit 101). Transcribed output:

```text
N02 P2 observed=accepted vertices=14 edges=24 faces=10 euler=-2 work=168608
outer and hole loops must not share walked endpoint IDs after sewing
```

Thus the source prediction was not merely assumed: the sewn body was accepted,
its Euler characteristic changed to -2, and the input-preservation assertion
passed before the expected-refusal assertion failed. The new production check
reads both endpoint IDs walked by each numerical coedge loop after the explicit
sewing remap. Outer/hole ID-set intersection refuses before model publication.
It neither changes the global topology checker nor relies on its Euler test.
Successful two-loop admission records `trim.boundary-endpoints-disjoint`
separately from the retained UV disjointness check.

After repair, all four `numerical_trim` tests pass, including the new regression,
original material/orientation and clean tube, zero and 0.005 mm gap success,
excessive gap and malformed-loop controls, invalid pairs and allowance-1 refusal.
Actual final transcribed summaries:

```text
N02 P2 observed=refused code=numerical.invalid-trim work=344 input=unchanged
N02 gap=0 result=numerical work=168644 exhausted=false repairs=0
N02 gap=0.005 result=numerical work=168644 exhausted=false repairs=1
```

One active 20M operation allowance surrounds sewing and rebuilding. The new
checks charge that same counter; there is no nested/reset allowance. The old
success count of 168612 remains historical evidence at the repair base, not the
new measured count. No tolerance, cap or expected answer was changed.

Portable commands, from the repository root with the existing authorized tool
roots/cache selected as environment variables:

```sh
cargo test --offline --release -p catallax-topology --test numerical_trim sewing_outer_and_hole_endpoints -- --nocapture
CATALLAX_N02_EVIDENCE="$N02_P2_EVIDENCE_DIR" cargo test --offline --release -p catallax-topology --test numerical_trim -- --nocapture
ci/run clippy -p catallax-topology
rustfmt --edition 2024 --check crates/catallax-topology/src/numerical.rs crates/catallax-topology/tests/numerical_trim.rs
git diff --check
```

The first command is the before-fail regression run prior to the production edit;
the second is the successful four-test final capture. Set `N02_P2_EVIDENCE_DIR`
to the absolute local destination for this repository-relative evidence directory
because Cargo tests run from the crate directory. Actual commands used the
existing Cargo and Rustup roots and reusable target cache; no installation or
permission change was needed. The original dual-bin manifest warning is
unchanged. Scoped official Clippy passes all targets with warnings denied;
changed-file rustfmt and diff hygiene pass. No broad campaign, application tests,
hosted gates or worker comparison was rerun for this narrowly scoped correction.

`near-gap.txt` and the refusal text files here are actual repaired test captures;
`shared-endpoints-refusal.txt` is the new P2 refusal. They contain no machine paths.
The original sibling N02 capture files remain byte-for-byte unchanged at their
original locations and in the old commit. Original report/validation/hash metadata
were also copied to authorized private scratch before revision. Root N02 artifact
digests identify the final public source, documentation and both capture sets.

Limits remain N01's numerical-only contract and pointwise sampling. This repairs
one explicit finite-loop topology rejection; it does not establish general 3D
separation, exact/bounded answers, STEP authority or global topology completeness.
No Answer/BoundedEvidence or generic topology code changed. No C++, private input,
companion repository, N03, new PR, comment or review posting, merge or deployment was
performed by this correction. Independent re-review remains with the supervisor.
