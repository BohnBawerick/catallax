# How the exact kernel uses every core

The exact kernel's STEP export uses every available core by default, and so do the two checks that take most of the time in a Boolean of flat and round faces. The rest of the Boolean can use them behind a setting and is slower that way today. An answer, a refusal, a work count and every output byte are the same at any thread count, one included, with the setting on or off. This page says how that is done, what a thread may not touch, where the threads work and where they do not pay yet. It is the first step of decision L1 of [PLAN.md](PLAN.md#124-named-for-later). The measurement that preceded it is R5.

## 1. The settings

The kernel uses as many threads as the process may run on. To force a count:

- `catallax-verify --threads N stream` and `catallax --threads N answer ...`, before the command;
- the `CATALLAX_THREADS=N` environment variable, for any program that links the kernel;
- `catallax_numbers::set_threads` from code.

With one thread no worker is started and no record is made. The kernel then runs the code it ran before this change.

`CATALLAX_RUN_AHEAD=loops`, or `catallax_numbers::set_run_ahead_loops(true)`, also hands the Boolean's loops to the other threads. [Section 5](#5-where-the-threads-work) says why that is off by default.

## 2. Why the loops cannot simply be split

One request has one work counter and one set of proof caches. A proof found in a cache costs less than a cold one, and it is charged less. The first decision that needs a root enclosure pays for it; every later one reuses it for one unit.

Split a loop across threads with a cache each, and every thread pays for the proofs it needs. The total grows, and it depends on which thread got which item. Share the caches instead, and the total depends on who finished first. Either way the work count stops being a property of the request. The count decides the refusal at the work allowance, so the answer would depend on the schedule too.

Model handles have the same problem. A handle is the position of an entity in its arena, so two threads that add entities at once would number them by arrival.

## 3. Run-ahead

So one thread still does everything, in order. It is the request's own thread, called the writer here. The other threads work ahead of it on copies and hand back records. A record is a result with a tape of the work that produced it. The writer adopts a record at the place where it would have done the work, and replays the tape against its own caches. That charges it what the work would have cost it cold. If there is no record, the writer does the work as it always has.

There are two kinds of record.

### A whole computation

`catallax_numbers::record` runs a computation on a helper thread and returns its value with its `Charges`. The writer adopts them in its own order and skips the computation, the work between its decisions included. This is for a computation whose result depends on nothing but values both threads share. The STEP export uses it for the spans of a pcurve fit.

A loop of such computations uses `values_ahead!`. The writer visits the items in order, and other threads compute later ones whole. Each item is computed once: the writer adopts a finished value with its charges and statistics, waits for one a worker has begun, and computes the rest itself. A search that stops at the first success works the same way; what the workers computed past the stop is dropped.

### One decision

A loop of independent proofs can hand its later items to worker threads with `run_ahead!`. A worker runs an item on its own copy of the model, with its own work counter and its own proof caches, and throws the copy away. It keeps a record of each decision that took it a quarter of a millisecond or more: the question, the outcome and the tape. When the writer reaches a decision and a record of the same question is ready, it adopts the record instead of doing the arithmetic.

A quick decision gets no record on purpose. Finding a record and proving that it answers the same question costs about as much as a quick decision does. A worker that has been at one decision for a millisecond enters it as pending. A thread that reaches a pending decision waits for the record, for as long again as the worker has already had, and then decides for itself. A worker the machine has set aside cannot hold the request up.

### The same question

A decision's work depends on which expression nodes are shared: one evaluation charges each distinct node once. Two questions are the same only if their expression graphs match node for node, under the same policy. `structure.rs` proves the match by walking both graphs. A fingerprint of the top few nodes finds candidates and proves nothing.

Nodes and roots made before the loop was handed over are shared by the writer and the workers, so they may only match themselves. Each node carries the count of hand-overs at the time it was made, which tells the two kinds apart without a walk through the shared part.

### The same charge

A decision's outcome depends only on the question and the policy. Caches change what a proof costs, never what it proves. That rule is older than this change, and the replay relies on it. The charge is the part that depends on the caches, so the tape is replayed against the writer's own caches:

- Work outside the proof caches is replayed as the amount recorded.
- Each algebraic root enclosure goes through the writer's real cache protocol, `enclose_root`. The protocol chooses between a reuse, a certificate the root holds and a cold proof exactly as it does without workers. Only the cold bisection is taken from the record.
- Each checked pi enclosure goes through `pi_enclosure` the same way.

A worker keeps the cold arithmetic of every site, including the ones it reused itself, so the writer has it whichever way its own caches fall.

A worker starts from the root certificates the writer held when it handed the work over, so it does not prove them again. A record that reused one of those says so. The writer adopts it only while it still holds that certificate, and the replay is then a reuse, as it would be cold.

### The same refusal

Each record knows the most any replay of it can charge. If that might not fit the remaining allowance, the writer ignores the record and does the work cold. The allowance therefore runs out at the same charge, and the request reports the same work, as with one thread.

A deadline is kept the same way. A helper works to the writer's deadline. The writer waits for a helper's result no longer than its deadline, and adopts nothing once the deadline has passed. It then does the work itself and refuses at its own next check.

A reserved allowance run inside recorded code starts a counter from zero and sets the proof caches aside, which no tape describes. The thread stays a worker inside it and leaves the shared caches alone, and whatever it was recording is dropped.

## 4. What a worker may not touch

Nothing a worker does reaches the writer except through a record. A cache that lives on a shared object would break that: a worker would fill it, the writer would find it filled, skip the charge and report less work than one thread does. So a worker neither reads nor fills these, and a recorded computation that meets one is dropped and done again by the writer.

| Cache | Where | A worker |
| --- | --- | --- |
| Certificates on an algebraic root | `AlgebraicRoot::enclosures` | uses its own record of cold proofs |
| Request proof caches: root certificates, root lists, enclosures, pi | `operation.rs` | has its own, started for each item from the writer's root proofs |
| A decider's root lists | `Decider::root_cache` | has its own decider |
| Taylor panels of a traced curve | `intersect` `lib.rs`, `jets` | computes without the cache |
| Constants of a trace | `intersect` `chart.rs`, `balls` and `interval` | computes without the cache |
| Controls of a spline span | `intersect` `span_hull.rs` | computes without the cache |
| Grid and moments of a rolling envelope, grid of a rolling saddle | `shapes` | computes without the cache |
| Bézier spans of a rational patch | `shapes` `rational_patch.rs` | fills a second cell the writer never reads |
| Intersection plans | `intersect` `pair.rs`, per thread | has its own; a recorded computation that reads them is dropped |

Cells whose fill is not charged, such as a node's coarse enclosure or a curve's equations, are shared freely. A worker that fills one saves the writer time and changes no count.

The rule for new code: a cache on a shared object, or one a thread keeps for itself, whose miss spends decision work must ask `catallax_numbers::shared_caches()` before it is read or filled, or call `private_cache()` when it is read. `crates/catallax-boolean/tests/shared_state.rs` lists every such place and fails when one is added without a line in the list.

## 5. Where the threads work

### By default, the STEP export

The fit of a pcurve cuts the curve into spans until each meets the tolerance. Helper threads compute the spans the cut is going to reach, as whole computations, and the fit adopts them in its own order. A fit that finishes within 2 milliseconds starts no helper.

This is where the time was. Seven of the P0 to P5 recipes take more than a second and export. For them the export was between 62% and 100% of the request at one thread, and its pcurve fits were most of that. [Section 7](#7-speed-on-p0-to-p5) has the times.

### By default, two checks in the Boolean

Two checks read the model and change nothing, so their loops hand over whole items with `values_ahead!`:

- whether a probe's path on a face is clear of each boundary edge, one item for each edge;
- whether a new intersection point is one of the points already known, one item for each known point.

In P1's Boolean they were about two fifths of the time at one thread. A loop hands over once its items take a fifth of a millisecond each and what is left adds up to two milliseconds.

Measured on the same laptop, in alternating runs of one binary at one thread and at 12. Work counts are the same in every run.

| Part | Machine load | Boolean and features, 1 thread | 12 threads |
| --- | ---: | ---: | ---: |
| P1 bracket | about 7 | 3.8 s, 4.8 s, 5.8 s | 3.5 s, 4.4 s, 4.5 s |
| P1 bracket | about 32 | 20.1 s, 22.8 s | 12.2 s, 12.9 s |
| P4 shelled tray | about 7 | 0.23 s to 0.29 s in four runs | 0.21 s to 0.25 s |
| P0 crossing pipes | about 6 | 0.54 s | 0.45 s |

On a quiet machine the gain on P1 is about a seventh, and one pair of the three went the other way. It is larger when the machine is busy. The other three fifths of P1's Boolean still run on one thread.

### Behind `CATALLAX_RUN_AHEAD=loops`, the rest of the Boolean

With the setting on, each of these loops hands its later items to workers once the loop has run for 20 milliseconds, or at once if its recent runs were that slow. A request under a work allowance wakes no worker in its first 20,000 units.

| Stage | Loops |
| --- | --- |
| Intersect | face pairs, carrier hits, carrier pairs, point faces, carrier pieces, isolated points |
| Imprint | boundary edge splits, section edge splits, face pieces, face cells |
| Classify | pieces; and for each point query, its faces and its ray crossings |
| Sew | shell orientation, void nesting |
| Searches that stop at the first success | probe arcs, probe depths, barrier equations, Bernstein coefficient signs; comparisons against known carriers and edges |

In a search the writer still tries the candidates in order and stops at the first that succeeds. The workers try later candidates in case the earlier ones fail. Their work is never charged unless the writer comes to need it.

The answers are right and the Boolean is slower than without the setting, so it is off by default. Measured on a 12-thread laptop that other jobs were using, one run each:

| Part | Boolean and features, 1 thread | With the loops handed over, 12 threads |
| --- | --- | --- |
| P1 bracket | 6.0 s | 13.5 s |
| P0 crossing pipes | 1.4 s | 1.8 s |
| P4 shelled tray | 0.36 s | 0.52 s to 0.93 s; 0.83 s with one worker |

Two things explain it. A record of one decision saves the writer only that decision's arithmetic. In P1's Boolean, decisions are about half the time, and those long enough to be worth a record are about a quarter of it, so workers that cost nothing could make it a third faster at most. And on this machine a busy thread costs. Other jobs already kept its six cores full, and the tray's Boolean at one thread went from 0.4 s to 0.6 s beside one busy process that shares nothing with it, and to 1.1 s beside eleven. A worker that runs a whole loop item to produce a few records is such a thread, and most of what it computes the writer computes again.

## 6. What they do not reach yet

- Most of the Boolean. The items that carry the rest of its time, the cells of a face in imprint, the classification of a piece, a carrier's hits and sections, add entities and facts to the model as they go, and later stages depend on what they leave. Run on scratch copies, with only their values kept, they make sewing refuse on P1, the seam, the pipes and the tray. Handing them over needs a record of what an item appends and which facts it reads and writes, replayed in order with the new handles renumbered. That is the large step of L1.
- A fillet request turns on enclosure reuse, a cache keyed by expression node. `enclose_at`, `enclose_all_at` and an unbounded `enclose_budgeted_at` are then decided cold, and a recorded computation that calls one is dropped.
- The STEP export outside the pcurve fit, the summary, and the certification of a fitted result run on the writer alone. The last is most of P5.
- Recipe steps run one after another. No P0 to P5 recipe has two independent expensive steps.

## 7. Speed on P0 to P5

One run of each P0 to P5 recipe through `catallax-verify`, at one thread and at the default, which was 12 threads on this laptop. Other jobs were loading the machine, with a load average near 20, so each time is one sample. The times are information and gate nothing. Every response, operation log, work count and file byte was the same in both runs.

| Part | Work units | Request, 1 thread | Request, all cores | Export, 1 thread | Export, all cores |
| --- | ---: | ---: | ---: | ---: | ---: |
| P0 block and sphere | 2,406,937 | 33.5 s | 6.64 s | 33.2 s | 6.45 s |
| P0 crossing pipes | 6,783,002 | 62.5 s | 22.2 s | 60.5 s | 21.0 s |
| P0 equal cylinder and sphere | 4,513 | 0.01 s | 0.01 s | 0.00 s | 0.00 s |
| P0 general cone (refused) | 2,785 | 0.01 s | 0.00 s | 0.00 s | 0.00 s |
| P0 general torus (refused) | 4 | 0.00 s | 0.00 s | 0.00 s | 0.00 s |
| P0 line tangent | 133,128 | 0.04 s | 0.03 s | 0.00 s | 0.00 s |
| P0 near-miss | 7,076 | 0.02 s | 0.01 s | 0.00 s | 0.00 s |
| P0 seam | 6,360,050 | 8.26 s | 1.87 s | 7.95 s | 1.63 s |
| P1 bracket | 44,824,339 | 44.7 s | 14.7 s | 36.6 s | 6.75 s |
| P1 drafted bracket | 45,051,525 | 65.1 s | 16.7 s | 51.8 s | 7.30 s |
| P2 turned collar | 7,199 | 0.02 s | 0.01 s | 0.01 s | 0.01 s |
| P3 filleted boss | 275,859 | 21.6 s | 6.01 s | 21.5 s | 5.98 s |
| P4 shelled tray | 1,733,492 | 1.53 s | 0.83 s | 0.95 s | 0.39 s |
| P5 flanged housing | 1,737,135 | 13.1 s | 13.9 s | 0.00 s | 0.00 s |
| All fourteen | | 250.5 s | 82.9 s | 212.6 s | 49.5 s |

What is left at all cores is mostly the Boolean and feature work on the one thread. For P1 that part took between 5 s and 14 s in runs of the same binary at the same setting, so a difference of that size between the columns is noise. P5 does not gain because most of it is the certification of a fitted result.

## 8. How it is checked

- `crates/catallax-numbers/tests/run_ahead.rs`: decisions over shared and fresh expressions, root isolation, recorded computations and loops of whole items give the same outcomes, statistics and work at 1, 2, 4 and 8 threads, with every record adopted and with none, at the full allowance and at allowances that run out part-way.
- `crates/catallax-cli/tests/threads.rs`: `catallax-verify` answers recipes at 1, 2, 4 and 8 threads and at the default, with the loops handed over and without. Responses, operation logs, work counts and file bytes must match, for full answers and for runs cut short by the allowance.
- `CATALLAX_RUN_AHEAD=wait` hands the loops over, makes every loop wait for its workers before the writer goes on and keeps a record of every decision, so every record is adopted. The tests use it to compare complete replays without depending on timing.
- `crates/catallax-boolean/tests/shared_state.rs` fails when a cell, lock or static a worker can reach has no line saying why it cannot change the writer's work.
- The same numbers test file and the unit tests in `crates/catallax-step/src/uv.rs` cover the edges: a nested allowance inside recorded code, a deadline that passes while the writer waits or before it adopts, and a writer that panics while helpers wait for it.
- The 477 saved answers replay in CI with two kernel threads, so their exports run with a helper thread. They do not replay with the loops handed over.
