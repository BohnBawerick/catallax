# P02 frozen M10/M8 exchange evidence

Partial Phase 1 evidence; not tangent or release qualification.

Baseline: current remote `main` observed as
`f4db345f2e7c51835d1ebeba9f6446d117144117` on 2026-11-14,
matching the original clean checkout. Isolated worktree branch:
`phase1/P02-freeze-f4db345`. All changes are confined to P02 paths.

## Provenance and method

The public baseline generator `tools/step-regressions/test_m10_spline_boolean.py`
provides M10-001..014. `cases()`/`write_questions()` were imported unchanged to
retain the exact 14 question/1 wrappers. It reads the two public loft fixtures
and embeds the existing `m27-loft-occt.step` fixture bytes with their SHA-256.
That is an existing public import input, not an oracle execution.
All 13 existing `crates/catallax-features/tests/fixtures/p5-later-cuts/*.json`
recipes were sorted by filename and wrapped in `catallax-step/1` then
`catallax-question/1`, source `recipe`, max_work 20,000,000.
No expected recipe was missing or synthesized.

`evidence/phase1/P02/` retains both exact JSONL inputs, inventory, source hashes,
commands, doctor, driver exit codes, and complete fresh run directories with
config, inputs, requests, responses, manifests, operation logs, portable
questions and artifacts. The stock driver used workers 1 then 4, with
`CATALLAX_THREADS` respectively 1 and 4, separate outputs and `--fresh`.
Config's unused default `max_work` is 100M; each portable question and actual
response retains its explicit 20M request allowance. No cache evidence is used.

### Portable publication copies

Files listed in `evidence/phase1/P02/publication.json` are sanitized portable
copies, not original raw logs. Original captures are retained privately outside
the repository. `${P02_CHECKOUT}` replaces the capture's absolute evidence
checkout location. `${P02_BINARY_ROOT}` replaces the separate checkout location
containing the observed binary. Expand these tokens to local paths before using
the recorded command arrays. Request `out` paths retain their historical
temporary directory names under the portable checkout token.

`publish.py` applies these replacements from a private capture, records every
sanitized file, runs `analyze.py` and regenerates `SHA256SUMS`. Supply the private
capture directory and its original `--checkout-root` and `--binary-root` paths.
It does not execute the driver or kernel. Comparison digests and `SHA256SUMS`
describe the published bytes. Case keys retain their original driver identities;
source, binary, input and STEP hashes remain unchanged. No host names were
present in the captured records. Responses, work and geometry bytes are unchanged.

A fresh build was attempted and failed with exit 127 because Cargo is absent.
The pre-existing release binary identifies baseline f4db345 and passes doctor
with request contract `catallax-request/1`. Its SHA-256 is
`1dfe54e616aee1aa353b3f6a07933ab4451f44a8791cf53ec4c748ffd6f48128`.
This establishes observed binary identity, not independent build reproducibility.

`analyze.py` compares complete responses after removing only response and
operation diagnostic `wall_ms` and manifest `stats`; work and operation decisions
remain compared. It checks actual artifact files and hashes STEP bytes.
`comparison.json` records every result, full refusal, charged work, comparison
verdict, STEP presence and per-case file hashes. Absent STEP on both sides is
recorded as absence, not successful export.

## Limits and remaining work

See [controls.md](controls.md) for tiny-loop, border-touch and tangent controls.
Their sources are present, but no compiled control executable or Rust toolchain
is available; none was freshly executed. No higher allowance, algorithm change,
Result adapter change, tolerance change or saved-answer change was made.

The tangent prototype's historical results remain source inference. Its source
still refuses border and curve touches and has no Boolean consumer. This wave
cannot qualify tangent handling. No OCCT import of newly emitted STEP was run;
STEP presence/strict equality alone does not establish external validity.
No latency, broad regression, hosted gate, sanitizer or memory-fault claim is
made. Refused budget prefixes are not completed proof costs.

## Executed results

All 27 requests completed freshly at each worker/thread setting (54 executions).
Full comparable responses, charged work and strict STEP bytes agree at 1/1 and
4/4. All four driver runs exited zero; zero crashes were reported.

| Case | Status / refusal | Work | STEP present (both) |
| --- | --- | ---: | --- |
| M10-001 | step.unsupported-geometry | 8,219,519 | False |
| M10-002 | step.unsupported-geometry | 8,219,568 | False |
| M10-003 | step.unsupported-geometry | 8,214,023 | False |
| M10-004 | numbers.operation-budget | 19,999,916 | False |
| M10-005 | numbers.operation-budget | 19,999,916 | False |
| M10-006 | numbers.operation-budget | 19,999,916 | False |
| M10-007 | step.unsupported-geometry | 2,527,978 | False |
| M10-008 | numbers.geometric-gap-unproven | 1,802,159 | False |
| M10-009 | intersect.parametric-coincidence | 115,984 | False |
| M10-010 | intersect.parametric-coincidence | 115,984 | False |
| M10-011 | intersect.parametric-coincidence | 115,984 | False |
| M10-012 | intersect.unsupported-pair | 74,693 | False |
| M10-013 | intersect.unsupported-pair | 74,693 | False |
| M10-014 | intersect.unsupported-pair | 74,693 | False |
| k01-wall-hole | numbers.operation-budget | 19,999,980 | False |
| k02-top-slot | intersect.parametric-surface | 2,180,807 | False |
| k03-rim-notch | numbers.operation-budget | 20,000,000 | False |
| k04-torus-pocket | numbers.operation-budget | 19,999,976 | False |
| k05-plane-notch | numbers.operation-budget | 19,990,133 | False |
| k06-cylinder-drill | numbers.operation-budget | 19,999,951 | False |
| k07-cone-drill | numbers.operation-budget | 19,926,846 | False |
| k08-sphere-dimple | numbers.operation-budget | 19,995,283 | False |
| k09-ring-groove | numbers.operation-budget | 19,988,811 | False |
| k10-axis-section | numbers.operation-budget | 19,999,276 | False |
| k11-offset-section | numbers.operation-budget | 19,998,787 | False |
| k12-far-disjoint | answered | 8,267,193 | True |
| k13-fillet-first | intersect.unsupported-pair | 3,471,214 | False |

All 14 M10 requests refused and left no STEP. Later cuts k01..k11 and k13
refused; only k12 far-disjoint answered and retained a seamed STEP file.
The later-cut total is 213,818,257 units per round, including refused prefixes.
