# R11 first step: public refusal census

**11 November 2026: 51 fresh requests completed; all four observed shell spans
answered.** Offset and thicken have no coverage in this input set. This replaces
the earlier missing-current-log finding with a bounded, newly generated census.
It does not establish the surviving-face refusal share needed to authorize
[PLAN.md R11 step two](PLAN.md#126-robustness-and-many-core-steps).

## Source, interval and coverage

Run baseline: current fetched main `439aa03061fef1078bc82688463fe7f627fde4a0`.
The clean detached source checkout produced `catallax-verify` with SHA-256
`13927727f7edc400ac1ae4b3769645720d6d9872fb1b316707140cac1458a289`.
Build started **2026-11-11 04:29:18 UTC**, completed successfully in 3m38s, using
rustc 1.99.1 (58cf8f9c3, 2026-10-13) and cargo 1.99.1.
The collection interval was **2026-11-11 04:33:23.747865–04:37:33.439355 UTC**,
not the whole day. All four fresh worker groups exited successfully; no cache
hits, observation cutoffs, timeouts, missing responses or internal errors occurred.
This is not a latency or independent 1/2/4-worker qualification.

[Preserved evidence](evidence/r11-2026-11-11.json) contains the complete emitted
operation-log objects, full refusal diagnostics, comparison outcomes, input and
raw response/log digests, per-group UTC intervals, build identity and eligibility
inventory. Evidence SHA-256:
`db54adfbf0df28103c725ea157b85b05b3482cd9fb69d3586a400cf76d87becb`.
Raw driver artifacts, doctor output and build log are retained outside the checkout.
The public repository's tracked inputs at the frozen revision establish provenance;
no external file imports or oracle adapters were admitted.

| Input group | Unique files | Executed requests | Answered responses | Refused responses | Work cap per request |
| --- | ---: | ---: | ---: | ---: | ---: |
| Complete P0–P5 recipes, B-rep (`recipe`) | 14 | 14 | 5 | 9 | 1,000,000 |
| Same recipes, field | 14 | 14 | 6 | 8 | 1,000,000 |
| Same recipes, sampled compare | 14 | 14 | 6 | 8 | 1,000,000 |
| Public saved `catallax-question/1` replay | 9 | 9 | 9 | 0 | 100,000,000, preserved from each input |
| Total, deduplicated files | **23** | **51** | **26** | **25** | No cap changes |

The 14 part files are eight `parts/P0/*.json` controls, P1 bracket and drafted
bracket, and P2–P5. The nine exported replays are the app's
[bracket-width-140 question](../crates/catallax-app/tests/fixtures/bracket-width-140.question.json)
and eight CLI `soliddraft-*` questions (flush/gui minimal and bracket, narrow,
narrow just-under/just-over, and slot-cut). All names and hashes are in evidence.
The existing [driver skill](../.claude/skills/verify-catallax/SKILL.md) defines the
1M parts profile and preservation of exported replay allowances. No eligible
saved request uses 20M; no additional profile or retry at a larger cap was added.

Discovery covered tracked JSON/JSONL in **`parts`, `crates`, and `.claude` only**:
255 files = 23 selected + **232 excluded**, individually listed with reasons.
These exclusions include bare recipe fixtures outside complete P0–P5, metadata,
P6, prefix cases and collections that are not exported question documents.
Other public roots were not searched for replays: this is explicitly not the whole
public test suite. The protected 477-case inputs were neither opened nor
enumerated; their private provenance is documented in
[RELEASING.md](RELEASING.md#merge-checks). C++/external oracle checks were
excluded entirely. No saved-answer inventory was read or modified for this census.

## Shell, offset and thicken by actual cause

Count only each execution's immediate `operation_log.operations` array.
Deduplicate by (group, driver artifact key, sequence); conflicting copies are an
audit error. The artifact key binds input, allowance and binary; sequence
identifies the kernel/operation span. Manifest/root copies add no executions.
Distinct adapters and the two kernels within compare remain separate executions.
There are **567 distinct executed spans**, of which **4** match the selected
operation names, all in **one unique case**, [P4](../parts/P4-shelled-tray.json).
An unexecuted span is neither a success nor a refusal.

| Operation | Executed | Answered | Unsupported scope | Invalid selection | Proof/work limit | Downstream measurement | Internal failure | Unknown refusal |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Shell | 4 | 4 | 0 | 0 | 0 | 0 | 0 | 0 |
| Offset | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Thicken | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Stable shell identifiers are `(parts-recipe, P4, brep, tray, sequence 1)`,
`(parts-field, P4, field, tray, 1)`, and compare's `(field, tray, 1)` and
`(brep, tray, 5)`. B-rep shell used 63,289 work units per execution; field shell
used 2,894. No shell refusal was observed (0/4 spans, 0/1 cases); the offset and
thicken rates are **undefined (0/0)**. There is no measured surviving-face failure
share. Conservative refusals elsewhere are not correctness defects.

## Refusals elsewhere in the same logs

**26 refused spans across 10 unique cases** occur outside the selected operations.
Diagnostics, not just code names, support these classifications:

| Actual cause/stage | Refused spans | Representative public identifiers and diagnostic |
| --- | ---: | --- |
| Boolean proof/work limit | 12 | Six cases in both recipe and compare, including P4: `numbers.operation-budget`, “operation work allowance exhausted” |
| Fillet proof/work limit | 1 | P5 recipe: same code and detail |
| Downstream field display work limit | 7 | P0 crossing-pipes/general-cone/general-torus/near-miss, both P1 variants, P4; display diagnostic preserves unresolved regional-tape operands, including empty lists |
| Unsupported certified intersection route | 4 | P0 general-cone and general-torus in recipe and compare: `intersect.traced-pair-required`, non-coaxial intersections require certified branch tracing |
| Unsupported fillet continuation | 2 | P5 field and compare: `fillet.unsupported-overflow`, contact region reaches another support and needs trimmed continuation |
| Invalid selection, downstream measurement, internal failure, unknown **span** cause | 0 | None logged in these categories |

These span counts intentionally differ from the **25 refused responses across
11 unique cases**: 21 root work-limit refusals and 4 root unsupported-scope
refusals. Two compare responses are `answered` with **`comparison.verdict =
unresolved`**, preserving the cone/torus B-rep refusals; only four compare
responses report sampled agreement. Agreement is not a global certificate.
Conversely, P2 compare exhausts its allowance after six answered construction
spans, with no refused span or stage attached. Its work-limit cause is known,
but its exact failing stage is **unclassified**, not inferred to be measurement.

P4's shell answers in all adapters, then recipe/compare refuse a later boolean
and field refuses display. Attributing those three final request refusals to
shell would misstate the evidence. No downstream measurement refusal is identified
by this run; zero observed does not establish measurement coverage.

## Reproduce and validate

Run from the census branch with Rust, Python 3 and normal build prerequisites.
Set `RUN`, `BUILD`, and `EVIDENCE` to fresh external directories; their values
below are placeholders, not repository paths. The collector invokes the existing
per-operation-log driver with four workers, preserving saved allowances.

```sh
MAIN=439aa03061fef1078bc82688463fe7f627fde4a0
RUN='<clean-source-directory>'
BUILD='<build-directory>'
EVIDENCE='<fresh-evidence-directory>'
git worktree add --detach "$RUN" "$MAIN"
(cd "$RUN" && CARGO_TARGET_DIR="$BUILD" nice -n 10 cargo build --release -p catallax-cli --bins -j 2)
sha256sum "$BUILD/release/catallax-verify"
python3 -B tools/research/r11_census.py --source-root "$RUN" --binary "$BUILD/release/catallax-verify" --out "$EVIDENCE" --workers 4
python3 -B tools/research/r11_census.py --count "$EVIDENCE/evidence.json"
# Count the committed observation without running geometry:
python3 -B tools/research/r11_census.py --count docs/evidence/r11-2026-11-11.json
```

Collector SHA-256: `ad00c8e31ceaf49b62c789de851c8bc5f4a6b9a056500520595e3f555cc84bba`;
frozen driver SHA-256: `7a3ed6daab44ec1062f1bed765c9345149348e5bb77e0a68588f7d0d08cf5cff`.
All 51 input, response and operation-log digests were independently checked against
retained artifacts; coverage, caps, ordered membership, arithmetic, clean source
revision and worker exits were checked. Log wall times and binary hashes can vary
on another build/machine. A collector doctor-parsing error was corrected before
any geometry request; its diagnostic is retained separately, outside the 51 runs.
No Rust test suite, full CI, saved-answer regeneration or C++ execution accompanied
this research run. No production source, tolerance or allowance changed.

## Historical control and next questions

The [22 October snapshot](checkpoints/2026-10-22-parts.json), documented in its
[checkpoint note](checkpoints/2026-10-22-verification.md), remains a separate
historical control: 39 logs, 240 spans, two shell refusals in one P4 case, both
unsupported-operation before work. Offset/thicken had no spans. Its revision was
`d98ac6a00552760399e2cd5f58b55a4a90179c0b`, snapshot SHA-256
`b3d1eef01215a9a679be1790badf3e3cb4ede49c6911fe63959cb5935fce7eee`;
no precise UTC interval exists. It is not included in today's denominators, and
its embedded historical oracle results are not used.

1. Preserve the construction/display distinction: P4 shell is not the observed
   bottleneck; its later boolean and display work limits warrant separate study
   within existing caps, without attributing them to face survival.
2. Broader offset/thicken and growing-thickness coverage is needed before deciding
   R11 step two. This run supplies no such cases and authorizes no prototype.
3. Keep unsupported tracing and fillet continuation separate from ownership and
   numeric refusal. The P2 compare root-only refusal also needs stage attribution
   before assigning a more specific cause. None justifies raising an allowance.
