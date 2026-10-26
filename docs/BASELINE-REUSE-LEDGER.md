# Baseline evidence reuse qualification ledger

Assignment: `9a2fead6b25ac7de`. Scope is the verification driver, its strict
comparator integration, public synthetic tests and documentation. Base is
`e6a93ac32e08519b990057cdb9fcb0d852f23987`. PR138 reporting at
`1089120e14f777a082c9f09308c10dd81d90a74f` remains an ancestor.

The checkout and its ancestors have no available `AGENTS.md`. The
`.agents/skills/verify-catallax` instructions, verification protocol, replay
contract and current check workflow were read. Open PR filenames were reconciled
before editing. No owned file overlaps the 12 then-open PRs. In particular,
the four read-only qualifications for 130/142/48/140 and the separately owned
finite-cone implementation were not altered or interrupted. Private proof
inputs and held-out data were not read. No real 477 replay was run here.

## Defects and dispositions

| ID | Observation / evidence | Disposition |
| --- | --- | --- |
| Q1 | Ordinary cache/resume entries do not establish an earlier fresh six-run qualification or complete build/environment provenance. | They cannot seed this opt-in archive. Unknown identity always selects six; only a fresh qualified six can publish an archive. |
| Q2 | Old run records omit worker exit codes. Generic manifest comparison does not independently audit complete work/log evidence. | New run records retain worker policy and exits. Qualification audits all exchange files and same-build non-time logs/work, plus the current strict comparator. Old evidence is ineligible. |
| Q3 | Normal JSON float parsing can collapse distinct numeric endpoints below binary64 spacing. | Qualification reads numeric JSON decimals exactly and retains full interval strings. A public regression demonstrates endpoints that ordinary float parsing merges. |
| Q4 | A copied historical baseline's original `cache_hit: false` alone does not describe its current use. | Qualification and copied run records explicitly label reuse and archive provenance. `report` labels historical timing separately. Original responses/work remain preserved. |
| Q5 | A hash-intact archive can still fail a changed validator or reveal a baseline defect. | Stop before new jobs and require diagnosis plus explicit requalification. Never retry automatically or trust an old `equal` flag. |
| Q6 | Independent review suspected incomplete operation records could pass if coherently changed in all three baseline logs. Executed public reproduction against `69bfe0c` confirmed an archive omitting `kernel` in all baseline logs/response copies was accepted after resealing. | Require every operation field, correct types, complete unsigned decision counts, status-dependent decisions/refusal, complete marks and bound/certificate records. All-worker omission/type mutations now fail on schema completeness, independently of determinism. |
| Q7 | Independent review suspected an adapter returning `[]` could qualify without workers. Executed public reproduction against `69bfe0c` confirmed `qualified`, a published archive and six runs containing zero cases. | Reject empty interpreted inputs before creating output or launching jobs, for fresh and reuse entry points; reject empty archived qualification evidence too. Ordinary adapter loading is unchanged. |
| Q8 | The new real face-measurement positive control exposed an owned over-rejection: face-check display observations have no exported OBJ, unlike field/print output. | Preserve that explicit existing schema by request format while retaining delivered-file presence/name checks for field and print. A real public torus face-measurement triplet passes alongside all export controls. |
| Q9 | Bounded recheck suspected coherently deleted whole entries could evade field validation. Executed public mutations against `60e0755` confirmed acceptance of empty histories, a deleted final entry and a deleted middle entry with renumbering across all baseline runs. The corresponding same-executable counter mutation also exposed a gap across the two triplets. | Compare ordered operation identity/status histories, including their counts, against the independent counterpart retained by the fresh six-run archive. Identical executables also require full non-time log/work equality across both triplets. All-baseline deletions and same-executable counter changes now fail. A real public comparison preflight returns legitimate empty logs and still qualifies/reuses with three fresh candidate jobs. Different-build work/counter changes remain allowed; changed execution histories require separate audit. |

Build provenance is a trusted build-time operator attestation bound to executable
bytes and the binary's source commit. Its schema requires complete source,
toolchain executable hashes and versions, resolved/native dependencies, build options/environment/config and
feature information. It does not authenticate a malicious attestation or infer
missing compiler inputs. Tests use explicitly synthetic attestations and never
claim those are measured build provenance.

The strict gate blocks changed certificates, newly answered refusals and changed
exports for the parent's separate preservation/geometry/export audits. It does
not implement accepted-difference exceptions. SINGLEPAIR and existing running
qualifications are unchanged. Any key mismatch, unknown provenance, missing or
corrupt evidence, or explicit fresh audit selects six. Candidate 1/2/4 always
executes fresh. An intact semantic baseline/validator defect instead stops for
diagnosis. No global cache, private path dependency or geometry source changes
were introduced.

## Measurements and validation

All measurements below are observations on a shared cloud environment. They
are separate scopes, not terms of an aggregate makespan or speedup claim.

| Completed observation | Result |
| --- | --- |
| Isolated locked release CLI build, two build jobs | Exit 0; 136.74 s elapsed. |
| Real CLI doctor | Exit 0; request protocol `catallax-request/1`; executable SHA-256 `14239b7e825157a79f3c7784903340bc4cf5bc43c54bbe6a367fbee76dc2948e`. |
| Public one-block explicit fresh-six smoke | All six fresh; strict gate passed; command span 447.42 ms, final gate audit 6.39 ms. |
| Final complete driver suite after review corrections, including PR138 reporting and public mutation tests | 51 passed in 31.376 s. |
| Focused qualification mutation coverage | 21 tests included in the final passing driver suite; includes coherent all-baseline field/type/entry deletions, same-executable counter histories, explicit nullable schemas, empty fresh/reuse adapter rejection, legitimate empty preflights and real face-measurement certificates. |
| STEP regressions with CI's pinned OCP 8.0.1.1.0 | 15 passed in 69.708 s. |
| Edge hunt / release evidence policy | 15 / 21 passed. |
| Formatting / generated comparison page / diff whitespace | Passed. |
| Workspace Clippy, all targets, warnings denied | Passed after prerequisite repair; 138.35 s elapsed. |
| Existing local workspace release suite, including doctests | Exit 0; 1,237 passed, 0 failed, 11 existing ignored across 165 targets. Command span including its serialized build: 2,129.75 s. Rust sources/binaries were from the pinned original base; current Python checks are recorded separately above. |

Initial Rust test/Clippy attempts stopped in GMP configuration because `m4`
was missing. The Debian trixie `m4` 1.4.19-8 package was extracted into isolated
scratch tooling after checking SHA-256 against Debian's package index:
`221b3e224708c1da2a1a88b6824380593579b2254ae7d914afa14cbd5d343334`.
Builds were serialized in an isolated target; active qualifications' binaries
were not replaced. No source/budget/assertion change was used to repair a tool
prerequisite. The first STEP attempt used the saved OCP 7.8.1.1.post1 and
completed with 22 subtest errors in 15 tests, including the unavailable
`TopoDS.Face` API. An isolated installation of CI's pinned OCP version resolved
that prerequisite mismatch; no geometry code changed. The workspace release
suite's terminal result is recorded with the PR.

The review reproductions retain the published `69bfe0c` implementation and
public synthetic archives outside the repository. Both source-only suspicions
were executed and confirmed before repair; their red test logs and the later
passing tests are preserved. The earlier head's hosted check and target-parts
workflows both completed successfully, but that result does not clear the
review defects or qualify the corrected head. No existing run was cancelled
or replaced to obtain these results.

The same saved environment briefly disconnected while starting the Q9
reproduction. The test-file edit persisted; the command had not started and no
log existed. After connection recovery, process/worktree reconciliation found
the earlier workspace suite terminal with exit 0 and no active local validation.
Only the previously unstarted focused reproduction resumed. No replacement
environment or duplicate workspace suite/build was created.

There is no saved-time claim: this implementation has not reused a real
matching future baseline. Synthetic archive reuse tests prove the control path,
not real-corpus performance. Parent-owned **FULL SIX after green CI** remains
mandatory for this driver/comparator change. Independent geometry and
changed-export audits remain outside this strict synthetic qualification.
