# How the exact kernel uses every core

The exact kernel's STEP export uses every available core by default, and so do four checks that take about half the time in a Boolean of flat and round faces. The rest of the Boolean can use them behind a setting and is slower that way today. An answer, a refusal, a work count and every output byte are the same at any thread count, one included, with the setting on or off. This page says how that is done, what a thread may not touch, where the threads work and where they do not pay yet. It is the first step of decision L1 of [PLAN.md](archive/PLAN-2026-11.md#124-named-for-later). The measurement that preceded it is R5.

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

A loop of such computations uses `values_ahead!`. The writer visits the items in order, and other threads compute later ones whole. Each item's value and charges are taken once: from a worker when its record is usable and fits the allowance, otherwise from the writer's own run. The writer waits a while for an item a worker has begun, no longer than its own deadline, and then runs the item itself. A search that stops at the first success works the same way; what the workers computed past the stop is dropped.

### Work that adds to the model

Most loop items in a Boolean are not values. The cell of a face in imprint, or the classification of a piece, adds points, coedges, faces and facts to the model while it runs, and later stages read what it left. Such an item can still run on a copy, with a journal.

A model only grows: entities are appended, facts recorded, names attached. `Model::journal_begin` starts a journal on a copy, and from then on every adding method writes what it was asked to add. `Model::replay` does the same additions to the model the copy was taken from, in the same order and through the same methods, so the tables those methods keep come out as if the work had run there: which coedges an edge has, which face owns a coedge, which facts follow from a new one. If the model has grown since the copy was taken, the handles the work made are moved past what it gained, and `Relocate` moves the handles in the work's result the same way.

The request's own thread runs imprint cells and piece classification this way, each on a copy taken at that moment and replayed at once. That is the same work in the same order as before, and it gives the same model, work counts and bytes.

Two changes are not journalled, because they belong to building and moving whole bodies and no loop item makes them: rewriting a body's face origins, and recording chart expansions. Either one on a journalled copy is a programming error and panics.

### What the work read

A worker's copy is older than the model its journal lands on. By the time the writer reaches an item, the items before it have added to the model. What a handle names never changes, but some answers do:

- whether a fact is in the fact book, which class a subject is in, what a class lies on, and what is tangent to a subject;
- the coedges of an edge and the owner of a coedge;
- the names kept beside faces, edges and bodies;
- how many surfaces, curves and points there are, and which bodies.

`Model::watch` starts a list of every such question the work asks its copy, the store's and the fact book's included. Copies of a watched copy write to the same list. The fact book lists what its own bookkeeping reads too: the classes a query walks through, and the sizes and ties a union compares.

Before the writer takes an item's result, `Model::unchanged_where_read` checks the list against the model as it is now: each entry the work read must be what it was in the model the copy was taken from. A question about something the work made itself always passes, since nobody else can have touched it. When every entry holds, the work would have found the same answers on the writer's model and so done the same there. The writer then adopts the item's charges and replays its journal. When one entry does not hold, the result is dropped with its charges and the writer runs the item itself. That happens when an earlier item recorded a fact this one would have used, for example.

The check is by entry, not by answer. A changed entry drops the result even where the answer would have come out the same.

### A cache the loop keeps

Some loops keep a cache of their own beside the model. The cells of one face share the traces of the face's boundary edges, and the graph's point and section loops share traces and point enclosures the same way. A hit charges nothing and a miss is computed and charged, so what an item costs depends on the items before it.

A worker runs an item with a copy of the cache as it was when the loop was handed over, and the copy counts what it did not hold. The writer takes the item only when that count is zero. The writer's cache holds everything the copy began with, so an item that only found values finds them here too, charges the same and leaves the cache as it was. An item that missed a value may find it in the writer's cache, filled by an earlier item, and would then have charged less than the worker did. The writer runs such an item itself.

The rule gives up the items whose misses the writer would have missed too. Today that loses nothing: a worker that evaluates a traced curve meets the guarded caches of section 4 and its record is dropped, so no item that misses a trace could be taken under any rule.

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

A deadline is kept the same way. A helper or a pool worker works to the writer's deadline. The writer waits for a helper's result no longer than its deadline, and adopts nothing once the deadline has passed. It then does the work itself and refuses at its own next check.

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
| Equations, nets and shared children of the spline pair search | `shapes` `rational_patch/pair` | never fills them; the writer fills them before each level, and a worker that finds one empty drops its record |
| Intersection plans | `intersect` `pair.rs`, per thread | has its own; a recorded computation that reads them is dropped |

Cells whose fill is not charged, such as a node's coarse enclosure or a curve's equations, are shared freely. A worker that fills one saves the writer time and changes no count.

The rule for new code: a cache on a shared object, or one a thread keeps for itself, whose miss spends decision work must ask `catallax_numbers::shared_caches()` before it is read or filled, or call `private_cache()` when it is read. `crates/catallax-boolean/tests/shared_state.rs` lists every such place and fails when one is added without a line in the list.

## 5. Where the threads work

### By default, the STEP export

The fit of a pcurve cuts the curve into spans until each meets the tolerance. Helper threads compute the spans the cut is going to reach, as whole computations, and the fit adopts them in its own order. A fit that finishes within 2 milliseconds starts no helper.

This is where the time was. Seven of the P0 to P5 recipes take more than a second and export. For them the export was between 62% and 100% of the request at one thread, and its pcurve fits were most of that. [Section 7](#7-speed-on-p0-to-p5) has the times.

### By default, four checks in the Boolean

Four checks read the model and change nothing, so their loops hand over whole items with `values_ahead!`:

- whether a probe's path on a face is clear of each boundary edge, one item for each edge;
- whether a new intersection point is one of the points already known, one item for each known point;
- whether each graph point lies on the arc of a boundary edge being split in imprint, one item for each point;
- whether each graph point lies on a section piece, one item for each point.

In P1's Boolean the first two were about two fifths of the time at one thread, and the last two about a twentieth. A loop hands over once its items take a fifth of a millisecond each and what is left adds up to two milliseconds. Each loop also remembers how its items went with the workers' help against the writer alone. While that has not been quicker, as on a machine whose cores are all busy, the loop stays with the writer and tries again one time in sixteen.

Measured on the same laptop in six alternating pairs of runs of one binary, at one thread and at 12, while other jobs held the load between 10 and 27. Work counts are the same in every run.

| Part | Boolean and features, 1 thread | 12 threads |
| --- | ---: | ---: |
| P1 bracket | 4.3 s to 19.2 s, median 6.0 s | 5.4 s to 6.6 s, median 5.5 s |
| P4 shelled tray | 0.27 s to 0.38 s, mean 0.31 s | 0.24 s to 0.36 s, mean 0.30 s |

At 12 threads P1's Boolean took about the same time whatever else the machine was doing. At one thread it followed the load. In the two quietest pairs one thread was the quicker, 4.3 s and 5.0 s against 6.4 s and 5.5 s, so on an idle machine expect no gain from this step alone. The times above are with the first two checks. The last two made no difference that alternating runs could show: 4.5 s against 4.8 s mean at 12 threads over five pairs. A little over half of P1's Boolean still runs on one thread.

### By default, five loops whose items add to the model

Five loops whose items add to the model hand over whole items, with `values_ahead!` and the journal of section 3:

- the cut of each face into pieces in imprint, one item for each face of the two bodies;
- the cell of each boundary walk of a face, one item for each walk;
- the class of each piece against the other body, one item for each piece;
- the faces that hold each point of the intersection graph, one item for each point;
- the section pieces of each carrier curve, one item for each carrier.

A worker runs an item on a copy of the model as it was when the loop was handed over. It hands back the item's value, its journal, the list of what it read and its charges. The writer checks the list when it reaches the item, adopts the charges, replays the journal and moves the handles in the value. An item that refused is never taken: the writer runs it, so a refusal comes from its own run.

The cells, the points and the carriers each share a region cache across their loop, and a worker's item is taken only when its copy of that cache held all the item asked for (section 3). A face that another thread cuts finds its cells in order on that thread. The cells of a face the writer cuts itself are handed over in turn.

A carrier whose first cut depends on a carrier the graph skipped has to ask where that carrier crosses, and it adds the answer to lists the carriers after it read. Another thread's result for such a carrier is not taken, and the writer cuts it itself.

At one thread on a 12-thread laptop, P1's sixteen Booleans took 6.5 s in one instrumented run: 2.5 s in intersect, of which 1.2 s in the carriers' hits and 1.1 s in their sections, 2.7 s in imprint, 1.2 s in classification and under 0.1 s in sewing. About seven tenths of that is in these five loops.

In an instrumented run on that laptop at 12 threads, with faces and pieces handed over and other jobs keeping every core busy, P1 took 30 face results and 128 piece results from workers, and the tray 9 and 37. None failed the read check. On each part three piece results came with a charge tape the writer could not use, and the writer ran those pieces itself, as it did the five on P1 and two on the tray that it stopped waiting for.

With all five loops handed over, an instrumented build ran P1's build three times on a 24-thread desktop. Each run took 10 to 117 face results, 95 to 174 cell results, 117 to 192 piece results and 12 to 15 carrier results from workers. None failed the read check or the cache condition. The writer ran itself the 38 to 47 items whose charge tape it could not use, which include every item that evaluated a traced curve, and the 8 or 9 it stopped waiting for. The work count was the same in every run.

That is not a faster Boolean. On the same desktop, at a load near 7, three binaries ran P1's build eight times each in alternating order: main before any loop that adds to the model was handed over, main with faces and pieces handed over, and all five loops. Work counts are the same in every run.

| P1 bracket, build only, 24-thread desktop | Before | Faces and pieces | All five loops |
| --- | ---: | ---: | ---: |
| 1 thread | 2.3 s to 2.9 s, median 2.6 s | 2.4 s to 2.9 s, median 2.6 s | 2.4 s to 3.0 s, median 2.6 s |
| 24 threads | 2.2 s to 2.6 s, median 2.4 s | 2.2 s to 2.5 s, median 2.3 s | 2.2 s to 3.0 s, median 2.4 s |
| CPU time over wall time at 24 threads | 1.0 | 1.1 | 1.3 |

P1's Boolean is about a tenth shorter at 24 threads than at one, and it was so before these loops were handed over. The workers do deliver results, and with all five loops they add about a quarter to the CPU time, but the request is no shorter for it. P1 is sixteen Booleans of about an eighth of a second each, made of many loops that last a few milliseconds. Why that leaves nothing to gain is not measured. Two things would explain it and neither is confirmed: a hand-over has a fixed cost in waking a worker and copying the model for it, which is large beside items of a millisecond; and a worker starts each item with the proof caches as they were at the hand-over, so it repeats proofs the writer's warmer caches would reuse.

The export is where every core pays: on the same desktop P1's export took 8.7 s at one thread and 1.1 s at 24.

### Behind `CATALLAX_RUN_AHEAD=loops`, the rest of the Boolean

With the setting on, each of these loops hands its later items to workers once the loop has run for 20 milliseconds, or at once if its recent runs were that slow. A request under a work allowance wakes no worker in its first 20,000 units.

| Stage | Loops |
| --- | --- |
| Intersect | face pairs, carrier hits, carrier pairs |
| Imprint | boundary edge splits, section edge splits |
| Classify | for each point query, its faces and its ray crossings |
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

The loops that hand over whole items hold about seven tenths of P1's Boolean at one thread: the faces of imprint with their cells about a third, the pieces of classification about a fifth, the sections of the carriers about a sixth and the faces of the graph points a fiftieth. What the threads do not reach:

- Three kinds of loop in the Boolean that merge into lists in order. They are the end of what this design can hand over without changing an answer's work count.
  - A carrier's hits, about a fifth of P1's Boolean at one thread (1.2 s of 6.5 s). Each hit is merged into the list of known points, and a new point is compared with every point found before it. An item run against an older list makes fewer comparisons, so it charges less than the same item does in order. Handing these over means separating the queries from the merge and running the merge afterwards in order. That changes the order of charged work at one thread too, and with it the work counts of saved answers, so it is a decision and not a step of this design.
  - The face pairs and carrier pairs of the intersect stage, about a hundredth of P1's Boolean together. They merge carriers and points into the same lists.
  - The splits of boundary and section edges in imprint, under a tenth of P1's Boolean. A split merges its cuts into the builder's lists, which the next split reads. The point checks inside them are already handed over as whole items.
- Sewing, under a hundredth of P1's Boolean.
- An item of the five loops above that evaluates a traced curve on a worker. The trace's panels and constants are caches on a shared object whose fill is charged (section 4), so the worker's record is dropped and the writer runs the item. Faces of crossing cylinders are such items. Taking them would need those caches to be recorded like root proofs are.
- A shorter Boolean. Section 5 has the measurement: the loops above hand over hundreds of items on P1 and its Boolean takes as long as before. The loops that could not be handed over are a quarter of its time, so they are not the reason.
- The next measurement is why P1's Boolean does not gain beyond one thread: on a 24-thread desktop its build takes 2.6 s at one thread and 2.4 s at 24, and the five loops add about 25% to its CPU time there for no time saved.
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
- The journal has unit tests in `crates/catallax-shapes/src/store.rs` and `crates/catallax-topology/src/model/journal.rs`: a journal replayed onto a store or model that has grown gives exactly what the same work gives when run there directly, moved handles included, and a journal kept inside another becomes part of it.
- The list of what work read has unit tests in the same two files and in `crates/catallax-numbers/src/facts.rs`. Additions the work never asked about pass the check. A fact the work could have used, a new coedge on an edge it listed, a new owner of a coedge it asked about, a count it took and a new body each fail it. Work whose reads hold replays onto the grown model exactly as it runs there directly.
- `a_value_that_no_longer_holds_is_dropped_with_its_charges` in the numbers test file: a worker's value the writer rejects is computed again by the writer and charged once, at full and at short allowances.
- `an_item_that_missed_the_loops_cache_is_the_writers_to_run` in `crates/catallax-boolean/src/region.rs`: three loops over one region cache through the real hand-over. Answers, statistics, work and cache contents are one thread's at 2, 3 and 4 threads and at allowances that run out part-way, and exactly the items whose copy of the cache held everything are taken. Without that condition the test fails on the work count.
- The same numbers test file and the unit tests in `crates/catallax-step/src/uv.rs` cover the edges: a nested allowance inside recorded code, a deadline that passes while the writer waits or before it adopts, and a writer that panics while helpers wait for it.
- The 477 saved answers replay in CI with two kernel threads, so their exports run with a helper thread. They do not replay with the loops handed over.
