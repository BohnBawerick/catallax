# Threads inside one field request

Baseline: `eee8fd123ce8a97e1115013d436f9db7e51d3043`.

## Question and candidates

R4 asks whether one immutable field plan can use the available processors while
keeping input order, output bytes and deterministic work at every worker count.
The earlier research probe reported identical bits at 1, 2, 4 and 8 threads and
5.7 to 6.0 times speedup. Raw evidence retained privately.

Candidate A is a shared atomic allowance. Reject it because scheduling would
choose which input first runs out of work. Candidate B is a shared worker pool
with ordered admission. Independent arithmetic runs in parallel; the request
coordinator charges reservations and appends geometry in the original order.
Checked numeric batches use fixed waves of at most 256 input items. Split the
remaining request allowance equally between inputs, with the first remainder
items receiving one extra unit. Shares are single-use, have cold proof caches,
and record actual reservations. Commit every share's accepted work in input
order before admitting results. Retry incomplete items in input order against
the remaining request cap. Failed attempts stay charged; unused units are not
charged. Both computed numeric work and admitted work stay inside the one cap.
The partition and retry order do not depend on workers.
Existing scalar mesh decisions retain their request caches and charging order.

## Frozen checks and decision rule

One implementation and measured round, with at most one repeat for an invalid
measurement. Compare P1 to P4 display and print at 1, 2, 4, 8 and default workers, using the canonical scoreboard's unchanged 1,000,000-work allowance and the default 0.01 mm accuracy and mesh limits.
Compare artifact bytes, logical work, cell counts and typed refusals. Include
empty input, non-finite inputs and allowances exhausted inside a batch item.
Keep all 477 saved replay answers and their allowances unchanged. Record one
1-thread versus all-core wall-clock observation as information, without a speed
gate. GPU evaluation and the exact B-rep kernel are outside this experiment.

Promote only when the correctness checks pass. The direction is to use every
available processor by default. Serial proof admission and ordered mesh assembly
remain necessary where later work depends on an earlier decision.

## Verdict: promote

R4 is adopted. The shared pool defaults to all available processors. One round
compared 50 case-runs at 1, 2, 4, 8 and 12 CPUs against the baseline: identical
artifact bytes, work, cell counts and typed refusals. No tolerance or allowance
was raised. The small-block controls produced matching display and print files,
with 442 and 4,494 work units. Six allowance tests, two worker tests, 47 selected
optimized regressions and the FEA query-count check passed. The numeric-cap
probe admits seven units and computes six arithmetic units under a ten-unit
request, then refuses the next item identically at every worker setting.

The table records whole-request seconds from the one measured round. This was a
shared machine with concurrent validation, not a timing pass. The 12-CPU total
was slower in this observation; speed is not an adoption gate. `budget` means
`numbers.operation-budget`; `coverage` means `field.print-coverage`. These are
baseline capability gaps, with no partial output accepted.

| Request | Outcome | Work | 1 CPU (s) | 12 CPUs (s) |
| --- | --- | ---: | ---: | ---: |
| P1 display | budget | 1,000,000 | 5.506 | 7.946 |
| P2 display | answered | 98,447 | 4.097 | 6.193 |
| P3 display | answered | 888,736 | 56.446 | 115.286 |
| P4 display | budget | 1,000,000 | 12.051 | 30.999 |
| P1 print | budget | 1,000,000 | 33.439 | 67.476 |
| P2 print | coverage | 3,299 | 0.034 | 0.007 |
| P3 print | coverage | 3,622 | 0.098 | 0.083 |
| P4 print | budget | 1,000,000 | 17.041 | 12.996 |

The original 477-case saved-answer replay remains a delivery gate at its
unchanged 100,000,000-work allowance. Its inputs and accepted answers are not
modified. Raw evidence retained privately. An initial baseline attempt at the
CLI's 100M default was interrupted without an answer and replaced by the
canonical 1M scoreboard run before this candidate round.
