## 2026-11-11: threaded unproved processor preview

Pieces: catallax-field

Change: `Plan::unproved_preview_batch(points, workers, storage)` returns
`Result<&[UnprovedValue], Refusal>` using reusable `UnprovedBatch` storage.
It shares the immutable plan across scoped processor threads and gives each
worker private registers. Results retain input order and the scalar F4 bits,
including canonical zero, on the same machine at any worker count. There is
no cross-point arithmetic or schedule-dependent reduction. This implements
the threaded batch part of PLAN M14 F5. It adds no regional pruning.

Failure: zero workers returns `field.preview-workers`, including on empty
input. Non-finite points and arithmetic overflow keep F4's refusals. All
workers finish before return; an evaluation failure from the earliest input
chunk wins over later failures. Thread creation failure or a worker panic
returns `field.preview-thread`. No failure retains partial results. Empty
input with positive workers returns an empty slice without starting threads.
At most the requested number of workers run, capped by the number of points;
one worker runs directly on the caller.

Allowance: F4's processor point evaluation performs no counted proof
arithmetic. F5 calls the same evaluator and
neither creates nor enlarges a request allowance. No shared numeric allowance
is needed for this path. Tests compare a subsequent proven answer and its
work count against a request with no preview, at 1, 2, 4 and 8 workers.

Fence: results remain `UnprovedValue`, with `PreviewHardware::Processor`.
There is no conversion to a proven result. The source fence allows the new
batch module and checks `UnprovedBatch` as well as the existing preview names.
Exports, proven measurements and proof routes remain forbidden callers.

Migration: none. Existing scalar, deterministic and proven entry points keep
their signatures and behavior. The F4 benchmark now also compares every batch
value against scalar F4 on the same seeded 2,000,000 uniform points per part
and near-surface points, at 1, 2, 4 and 8 workers. Its shared-machine times are
informational. The starred 60 ns P3 / 150 ns P1 timing pass still needs a quiet
machine under the timing rules in ARCHITECTURE.md. Raw evidence retained privately.
