# How parallel the exact kernel could be

Research item R5, adopted from the Csanády interview study: measure how parallel the exact B-rep could be, and find out whether "a B-rep cannot use many cores" holds for Catallax. The decision of 11 November 2026 in [ARCHITECTURE.md](ARCHITECTURE.md#write-once-storage-handles-coordinates-as-a-cache) adopts parallel execution in both kernels now, including L1 (plan in parallel, append in order) for the B-rep; this report is the map of where that parallel work goes. This page is a source observation at `7765406` plus the `work-span` observer added with it. The observer runs the kernel on one thread; run-ahead workers do not change the writer's order of work or its work count. It changes no kernel decision, work charge, allowance, tolerance or cap. Raw evidence is retained privately.

## 1. Answer

At the grain of the kernel's existing loops, a rebuild of P0 to P5 exposes a parallelism of 1.0 to 3.3: the speed-up from running the loops' iterations side by side as the code is written, on the dependency graph section 2 describes. Allowing speculative candidate searches, removing the waits on shared bookkeeping and merging results in order (L1) raises the bound to between 1.0 and 8.2, with a median of 2.3 over the 14 cases. Only the P1 brackets (8.2) pass 4. As the code stands, Csanády's claim holds: the cost is proof arithmetic whose result does not depend on order, but the work sits in a few long proof chains, one Boolean after another. The sections below show which chains to break first.

Four things bound the span:

- **One task dominates its loop.** In line tangent, the search for a certified probe point in the imprint's face cells is 94% of the critical path; in P4 it is 45% as written and 49% under L1. Running the search's candidates, barriers and Bernstein coefficients side by side is what lifts P1 from 3.30 to 5.61 and seam from 2.39 to 3.19.
- **Single proofs are indivisible at this grain.** Crossing pipes spends 475,533 units in one root isolation (`solve.rs:266`) inside one ray crossing during void nesting. Its build cannot go below that call: 3.90M / 0.48M = 8.2 is a hard ceiling at decision grain. Finer parallelism needs the isolation's own subdivision run in parallel.
- **Recipes are chains.** No P0 to P5 recipe has two expensive steps that are independent: every Boolean consumes the previous result, and tool bodies cost nothing. Feature-level parallelism is exactly 1.00 in all 14 cases.
- **Shared bookkeeping and ordered merges serialise what remains.** Fact-book unions, the boundary-trace caches and the request proof caches account for the gap between the speculative and unshared columns (crossing pipes 1.60 to 2.32, P1 5.61 to 6.02). Point and carrier merges, each an exact comparison against every earlier entity, are why L1 helps P1 (6.02 to 8.17) and no other case.

STEP export is outside the rebuild and was not divided into tasks. It is a serial 43% of crossing pipes' request work, 57% of block and sphere's, 33% of seam's and 92% of P3's. The pcurve fit inside it runs ahead on every core since #414, which the observer does not split into tasks; its work count is the same at any thread count.

The 477 saved replay inputs give the same picture (section 4.4): by summed work 1.31 as written and 1.64 at best, median 1.31 and 1.43 per input.

## 2. Method

### The observer

The `work-span` feature of `catallax-numbers` replays one request on an ideal machine with unlimited workers and free scheduling, in request work units. Code runs in its normal order and the observer only keeps clocks, so the measurement needs no quiet machine. Work charged on the current strand advances its clock. A loop marked as parallel starts each task at the time the loop body reached it rather than after the previous task. A read of a tracked resource moves the reader's clock to the time its last writer wrote it. Work, reads and writes in a parallel loop but outside its tasks join every task so far. The span is the finishing time of the whole request; work divided by span bounds the speed-up of any schedule of the observer's task graph. Every run also checks that the observed work equals the request counter and that the critical path's segments sum to the span.

The observer reads the request counter from `operation::charge`. It never charges work, asks a numeric question or changes an outcome, and ordinary builds compile none of it.

```sh
cargo build --release -p catallax-cli --features work-span --example work_span
target/release/examples/work_span parts/P0/seam.json --out "$PROOF/seam" --max-work 100000000
```

The output gives the request work, the build's work (the recipe step loop), each model's span, its parallelism and critical path by region, every marked loop with its tasks and per-model span, the tracked dependencies that delayed a strand, and the largest top-level decider calls.

### Five models

One run keeps five clocks. They differ only in which loops run their tasks side by side and which reads wait.

| Model | Side by side | Reads that wait | Extra work it ignores |
| --- | --- | --- | --- |
| Steps | Recipe steps only | All | None |
| Independent | Recipe steps and every loop whose iterations all run in the sequential build | All | None |
| Speculative | Also searches that stop at the first accepted candidate | All | Candidates after the accepted one |
| Unshared | As speculative | All except fact books, indexes and caches | Recomputing what a shared cache or fact would have saved |
| Ordered | As unshared; each task's merge into the loop's shared list leaves the task for one commit chain per loop, in task order | As unshared | As unshared |

Independent is what the code allows as written. Ordered is decision L1 applied to the loops below. Speculative, Unshared and Ordered are optimistic: running a search side by side does work the sequential build skips, and not sharing a cache or fact means repeating its proof. Neither cost is counted, so those columns are bounds, not predictions.

### What is a task

| Loop | Where | Kind |
| --- | --- | --- |
| Recipe steps | `catallax-features` `evaluate_impl`; question steps in `catallax-cli` `build_into` | Steps |
| Face pairs, carrier hits, carrier pairs, point faces, carrier pieces, isolated points | `catallax-boolean` `graph.rs` | Parallel |
| Boundary splits, section splits, face pieces, face cells | `imprinting.rs` | Parallel |
| Classify pieces | `stages.rs` | Parallel |
| Point sides, ray crossings | `containment.rs` | Parallel |
| Orient shells, nest voids | `sewing.rs` | Parallel |
| Envelope nodes, bound spans, clearance spans | `catallax-shapes` `rolling_envelope.rs` | Parallel |
| Probe candidates, probe depths | `points.rs` | Search |
| Probe barriers, barrier equations, Bernstein signs | `probe_path.rs` | Search |
| Point, carrier and edge merges | `canonical_point`, `add_curve`, `Builder::edge` | Search inside an ordered merge |

Each parallel loop was read for values carried from one iteration to the next. Those that change results are tracked as resources; collections that take one entry per task, such as the classification map, are allowed. The Boolean stages run in sequence.

### What is tracked

| Resource | Read | Written |
| --- | --- | --- |
| Arena entities: every surface, curve, point, vertex, edge, coedge, loop, face, shell, lump and body, keyed by the copy of the arena that appended them | `Arena::get` | `Arena::push` |
| Arena length, one per copy of an arena | `len`, `is_empty`, `handles` | Every push |
| Topology side tables: edge coedges, coedge owner, face origins, edge constructions, blend base | Their accessors | Their mutators |
| Fact book: component parent, size, incidence sets; tangents and the fact list per subject | `holds`, `facts_about`, `iter` | `record`, unions |
| Conic support index | `has_conic_incidence` | `record_fact` |
| Recipe step outputs and the fillet, draft and chamfer registries | Step inputs | Step end |
| Request proof caches: root certificates, root lists, enclosures, pi, intersection plans | Hit | Insert |
| Boolean caches: face bounds, trim bounds, edge traces, trace projections, point coordinates, surface pairs | Lookup | Insert |
| Envelope fit node rows reused by the next span | Next span | Span |

Every resource names one entry of one store, cache or copy. Each arena and cache has its own number, and every clone of an arena gets a fresh one: a copy shares the entries it inherited, with their writers, and appends under its own number. Entries of tables that can be rewritten (the fact book, the topology side tables and the conic index) take a fresh generation on every write, and each copy keeps the generations it inherited, so a later write to either copy, including a rewrite of an inherited entry, is invisible to the other. Probe copies of the model therefore no longer alias the model they came from. Cache keys are dense entry numbers, never hashes. A read waits for the last write in program order. Handle numbers are allocation order, not a dependency: an ordered append reproduces them.

The graph still over-states dependencies in three named places, each a conservative choice: a fact-book read waits for every change to a component it walks, not only the changes its answer depends on; a read of an arena copy's length waits for that copy's latest append; and work between the tasks of a loop joins every task before it. So the measured parallelism is a bound for the observer's graph and a lower estimate for the kernel at this grain, and finer grain can expose more. A dependency carried by a value outside every tracked resource would make it an over-estimate; the loop audit found none.

Not tracked, with the reason: `Decider::root_cache` hits charge the full certificate cost again, so they cannot change work; lazily filled `Real` metadata, curve-equation cells, jet programs and Gröbner bases are not charged; chart and face expansion maps are written only while their body is built. A value passed between tasks through a local variable that is none of the above would be missed; the loop audit found none.

## 3. Ordinary builds are unchanged

Every hook is behind `#[cfg(feature = "work-span")]`. With default features, this change charges the same work as `7765406` on all 14 P0 to P5 cases, from 4 units (general torus) to 44,017,730 (drafted P1). With 1, 2 and 4 workers, and in the observed run, each case gave the same answer digest as `7765406`. Main's own changes since the first version of this report (#370, #371 and #414) lowered the work of block and sphere, crossing pipes, seam and both P1 brackets; every figure below is re-measured on them. Observed runs of all 477 saved replay inputs gave the accepted answer digest and observed exactly the request counter's work. The `rounded STEP speed` check compares native work between base and head on every pull request that touches the numbers, Boolean, intersection or STEP crates.

## 4. Results

Work in units. Parallelism is work divided by span, so 1.00 is fully serial.

### 4.1 Build (the recipe steps)

| Case | Outcome | Request work | Build work | Steps | Independent | Speculative | Unshared | Ordered | Largest decision |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| crossing pipes | answered | 6,783,002 | 3,897,430 | 1.00 | 1.60 | 1.60 | 2.32 | 2.32 | 475,533 |
| line tangent | answered | 133,128 | 131,052 | 1.00 | 2.02 | 2.02 | 2.03 | 2.02 | 73 |
| equal cylinder and sphere | answered | 4,513 | 3,457 | 1.00 | 2.29 | 2.43 | 2.43 | 2.43 | 73 |
| block and sphere | answered | 2,164,259 | 938,819 | 1.00 | 2.29 | 2.40 | 2.48 | 2.34 | 44,372 |
| seam | answered | 1,268,420 | 846,990 | 1.00 | 2.39 | 3.19 | 3.50 | 3.01 | 14,740 |
| near miss | answered | 7,076 | 4,202 | 1.00 | 2.22 | 2.23 | 2.25 | 2.25 | 73 |
| general cone | `intersect.traced-pair-required` | 2,785 | 2,785 | 1.00 | 1.78 | 1.78 | 1.78 | 1.78 | 73 |
| general torus | `intersect.traced-pair-required` | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1 |
| P1 bracket | answered | 43,790,544 | 42,152,303 | 1.00 | 3.30 | 5.61 | 6.02 | 8.17 | 70,654 |
| P1 drafted bracket | answered | 44,017,730 | 42,377,532 | 1.00 | 3.30 | 5.62 | 6.03 | 8.18 | 70,654 |
| P2 collar | answered | 7,199 | 2,751 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 284 |
| P3 filleted boss | answered | 275,859 | 22,542 | 1.00 | 1.65 | 1.79 | 1.80 | 1.75 | 6,429 |
| P4 shelled tray | answered | 1,733,492 | 1,554,975 | 1.00 | 2.62 | 2.94 | 2.94 | 2.93 | 15,708 |
| P5 flanged housing | answered | 1,737,135 | 1,077,799 | 1.00 | 1.75 | 1.75 | 1.78 | 1.83 | 8,420 |

The largest decision is the largest single top-level decider call. Without parallelism inside it, no schedule finishes the build in less.

Ordered can fall below Unshared. It sends every merge through one chain per loop, even merges that would not have waited for each other; block and sphere's, seam's and P3's did not.

### 4.2 Whole request, including summary and STEP export

| Case | Steps | Independent | Speculative | Unshared | Ordered |
| --- | ---: | ---: | ---: | ---: | ---: |
| crossing pipes | 1.00 | 1.27 | 1.27 | 1.49 | 1.49 |
| line tangent | 1.00 | 1.99 | 1.99 | 1.99 | 1.99 |
| equal cylinder and sphere | 1.00 | 1.76 | 1.82 | 1.82 | 1.82 |
| block and sphere | 1.00 | 1.32 | 1.34 | 1.35 | 1.33 |
| seam | 1.00 | 1.63 | 1.85 | 1.91 | 1.81 |
| near miss | 1.00 | 1.48 | 1.49 | 1.49 | 1.49 |
| general cone | 1.00 | 1.78 | 1.78 | 1.78 | 1.78 |
| general torus | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| P1 bracket | 1.00 | 3.04 | 4.78 | 5.07 | 6.44 |
| P1 drafted bracket | 1.00 | 3.04 | 4.80 | 5.08 | 6.45 |
| P2 collar | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| P3 filleted boss | 1.00 | 1.03 | 1.04 | 1.04 | 1.04 |
| P4 shelled tray | 1.00 | 2.25 | 2.45 | 2.45 | 2.44 |
| P5 flanged housing | 1.00 | 1.36 | 1.36 | 1.37 | 1.39 |

### 4.3 By Boolean stage, measurement and export

The table [PLAN.md](PLAN.md#126-robustness-and-many-core-steps) asks for: total work, the longest chain that must run in order, and their ratio, for the three slow P0 cases and P1 to P5. Work by stage is summed over every Boolean in the build.

| Case | Intersect | Imprint | Classify | Sew | Other build | Measurement and certification | STEP export |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| block and sphere | 234,821 | 477,201 | 226,353 | 444 | 0 | 0 | 1,225,440 |
| crossing pipes | 16,589 | 1,208,906 | 62,820 | 2,609,115 | 0 | 0 | 2,885,572 |
| seam | 101,910 | 549,153 | 182,002 | 13,925 | 0 | 8 | 421,422 |
| P1 bracket | 17,432,253 | 18,042,019 | 6,667,736 | 9,585 | 710 | 0 | 1,638,241 |
| P1 drafted bracket | 17,628,244 | 18,064,462 | 6,671,875 | 9,817 | 3,134 | 0 | 1,640,198 |
| P2 collar | 0 | 0 | 0 | 0 | 2,751 | 0 | 4,448 |
| P3 filleted boss | 3,722 | 5,435 | 3,156 | 6,724 | 3,505 | 0 | 253,317 |
| P4 shelled tray | 462,478 | 879,723 | 208,631 | 309 | 3,834 | 32 | 178,485 |
| P5 flanged housing | 145,216 | 28,281 | 75,613 | 8,609 | 820,080 | 659,336 | 0 |

Other build is work outside the Boolean stages: P5's fillet and pattern, P2's revolve and chamfer. P5's certification of the fitted result writes its STEP file inside the certification column. Measurement and export run once, after the build, and were not divided into tasks: their longest chain is their work.

Longest chain by stage, with work over it in brackets:

| Case | Intersect, Independent | Intersect, Ordered | Imprint, Independent | Imprint, Ordered | Classify, Independent | Classify, Ordered | Sew, Independent | Sew, Ordered |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| block and sphere | 67,862 (3.46) | 78,335 (3.00) | 270,379 (1.76) | 269,272 (1.77) | 70,519 (3.21) | 53,370 (4.24) | 391 (1.14) | 248 (1.79) |
| crossing pipes | 8,779 (1.89) | 9,642 (1.72) | 200,043 (6.04) | 198,930 (6.08) | 16,862 (3.73) | 11,062 (5.68) | 2,219,260 (1.18) | 1,462,223 (1.78) |
| seam | 41,398 (2.46) | 71,501 (1.43) | 209,925 (2.62) | 184,355 (2.98) | 89,808 (2.03) | 20,551 (8.86) | 13,916 (1.00) | 4,863 (2.86) |
| P1 bracket | 8,675,516 (2.01) | 3,065,270 (5.69) | 3,016,062 (5.98) | 1,747,012 (10.33) | 1,061,749 (6.28) | 339,690 (19.63) | 8,837 (1.08) | 4,292 (2.23) |
| P1 drafted bracket | 8,767,190 (2.01) | 3,085,796 (5.71) | 3,018,370 (5.98) | 1,747,316 (10.34) | 1,061,625 (6.28) | 339,530 (19.65) | 8,893 (1.10) | 4,348 (2.26) |
| P3 filleted boss | 1,596 (2.33) | 1,700 (2.19) | 1,008 (5.39) | 557 (9.76) | 827 (3.82) | 426 (7.41) | 6,704 (1.00) | 6,673 (1.01) |
| P4 shelled tray | 160,347 (2.88) | 148,556 (3.11) | 372,309 (2.36) | 365,902 (2.40) | 55,608 (3.75) | 12,595 (16.56) | 278 (1.11) | 252 (1.23) |
| P5 flanged housing | 105,782 (1.37) | 87,743 (1.66) | 1,072 (26.38) | 455 (62.16) | 9,756 (7.75) | 9,136 (8.28) | 8,440 (1.02) | 8,369 (1.03) |

A stage's longest chain is the time it is open on its strand, summed over the build's Booleans. P1's 35 steps run one after another, so its stage chains add up to its build's. Ordered can lengthen a stage, as in seam's and block and sphere's Intersect: its merges wait for one commit chain even when they would not have waited for each other.

### 4.4 Saved replay inputs

All 477 saved replay inputs were observed at the replay allowance of 100M units: 458 answered and 19 refused, 5.27G units in total. Every observed answer matched its accepted digest. Parallelism here is for the whole request.

| Model | Work over summed span | Median case | 90th percentile | Largest case |
| --- | ---: | ---: | ---: | ---: |
| Steps | 1.00 | 1.00 | 1.00 | 1.00 |
| Independent | 1.31 | 1.31 | 2.35 | 7.61 |
| Speculative | 1.38 | 1.38 | 2.81 | 14.35 |
| Unshared | 1.64 | 1.43 | 3.41 | 14.36 |
| Ordered | 1.63 | 1.42 | 3.46 | 9.38 |

No saved input has two independent steps: every step after the first consumes the one before, which is why Steps is 1.00 throughout. The single-Boolean families (unite 85, subtract 89, intersect 86 inputs) stay at or below 1.05 in every model by summed work. The 217 mixed chains carry 73% of the work and reach 1.46, 1.58, 2.08 and 2.08 in the last four models. Only 37 inputs pass 4 in Unshared and 32 in Ordered. The largest values all come from one 48,881-unit intersection (replay index 109); among the 320 inputs above 1M units the largest is 8.20 Unshared.

The heaviest answered inputs, all mixed chains of three or four steps:

| Replay index | Work | Independent | Speculative | Unshared | Ordered |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 356 | 93,099,262 | 1.39 | 1.40 | 1.63 | 1.61 |
| 191 | 88,793,433 | 1.57 | 2.18 | 3.49 | 3.46 |
| 337 | 83,638,495 | 1.44 | 1.45 | 2.14 | 2.15 |
| 294 | 83,118,235 | 1.57 | 1.61 | 2.21 | 2.21 |
| 359 | 82,729,459 | 1.46 | 1.47 | 2.63 | 2.60 |

### 4.5 Where the span concentrates

Shares of the request span in the Independent and Ordered models, by innermost region; the Boolean step and stage name the feature.

| Case | Independent: largest regions on the critical path | Ordered: largest regions on the critical path |
| --- | --- | --- |
| crossing pipes | STEP export 54.2%; Sew, void nesting, ray crossings 23.3%; void nesting, point sides 12.5% | STEP export 63.2%; void nesting, ray crossings 12.7% (one root isolation of 475,533); point sides 12.5% |
| block and sphere | STEP export 75.0%; imprint face cells, probe depths 14.2% | STEP export 75.3%; imprint face cells, probe depths 14.3% |
| seam | STEP export 54.3%; imprint face cells, probe depths 16.9%; classification, Bernstein signs 9.3% | STEP export 60.0%; imprint face cells, probe depths 24.8%; point merge 7.7% |
| line tangent | Imprint face cells, probe depths 93.7% | The same, 93.9% |
| P1 bracket | Final `bracket` cut's carrier hits, point merge 12.0%; STEP export 11.4%; slot cuts' point merges 4.5% to 5.0% each | STEP export 24.1%; final cut's point merge 6.8%; slot cuts' imprint probe depths 3.0% each |
| P3 filleted boss | STEP export 94.9% | STEP export 95.2% |
| P4 shelled tray | Imprint face cells, probe depths 44.9%; STEP export 23.2%; carrier hits 11.7% | Imprint face cells, probe depths 48.7%; STEP export 25.2%; point merge 9.4% |
| P5 flanged housing | Certification, display and STEP export of the fitted result after the build 51.6%; the fillet step outside its marked loops 34.2% | The same, 52.8% and 35.0% |

In P1, 18 steps waited for an earlier step's output: the cuts form one chain.

Dependencies that delayed a strand most, in the Speculative model (\* shared bookkeeping, ignored by Unshared and Ordered):

| Case | Largest delays, in units summed over waits |
| --- | --- |
| crossing pipes | trace projection cache\* 6.7M; fact links\* 5.5M; fact class\* 1.3M; root list cache\* 0.6M |
| seam | trace projection cache\* 315k; root certificate cache\* 42k; fact list\* 30k |
| P1 bracket | fact class\* 919M; recipe step output 22.9M; edge entities 21.7M; trace projection cache\* 8.9M |
| P4 shelled tray | fact class\* 4.4M; trace projection cache\* 3.0M; edge entities 1.4M |
| P5 flanged housing | fact class\* 797k; enclosure cache\* 213k; root list cache\* 170k; envelope node rows 162k |

Delays add up over many waits on overlapping intervals, so they rank causes; they are not shares of the span.

## 5. What the three single-threaded places would have to change

The rule that holds everything else in place: an answer, a refusal and the work count must be the same at any worker count ([ARCHITECTURE.md](ARCHITECTURE.md#what-stays-deterministic-and-proved)). Work decides refusals at the allowance boundary, so a task's work must not depend on what ran beside it.

### 5.1 The single writer

Every operation takes `&mut Model`. The geometry and topology arenas are append-only and persistent (`im::Vector`), handles are issued in order, and the fact book is a union-find whose unions rewire incidence sets in place.

- **Plan against a snapshot, commit in order.** A clone of the model is cheap, so each task can plan against the snapshot at its loop's start and return candidate entities under local names. A commit in task order appends them, which reproduces the sequential handle numbers, and runs the merges that are dedup comparisons today: `canonical_point`, `add_curve` and `Builder::edge`. That is the Ordered model. Its commit chain is where P1's remaining span sits, so the merge itself must stop comparing every new point against every earlier one by exact proof. Separating candidates first by certified boxes, then running the remaining comparisons side by side, is what the point-merge search measured.
- **Fact-book reads must see exactly the facts the sequential order would.** A fact proved by an earlier task can save a later task a proof. If the later task ran on the snapshot, it would do that proof itself, its work would change, and so could a refusal at the boundary. Either the reader waits for the earlier task's commit, which is the measured delay, or it runs speculatively and is replayed when a fact it could have used lands first. Unions should be recorded as commit-time events rather than in-place rewiring, so a snapshot's components never change under a reader.

### 5.2 The per-thread allowance

`with_operation_budget` keeps the request's limit, used work and exhaustion latch in a thread-local `ACTIVE`, and `charge` latches exhaustion at the first reservation that does not fit. `OperationAllowance::split_off` divides a limit ahead of time, which would move where refusals happen.

- **One counter per task, decided in order after the join.** Each task charges its own counter with the request's remaining limit as a cap. After the join, prefix sums in task order give the sequential total, the first task that would have overflowed, and the exact `used` figure the sequential run reports. Tasks after that point are discarded. This needs every task's work to be independent of scheduling, which section 5.3 provides.
- **Speculation charges only what the sequential run would have tried.** A search run side by side keeps a counter per candidate and charges candidates up to and including the first accepted one; the rest is discarded and never reaches the request counter.
- **Wall-clock deadlines stay as they are.** The legacy 50-second deadline was never deterministic and decides nothing in the verification exchange.

### 5.3 The per-thread proof caches

`ROOT_PROOFS`, `ROOT_LIST_PROOFS`, `ENCLOSURE_PROOFS` and `PI_PROOFS` (`operation.rs`) and the intersection `PLANS` (`pair.rs`) are thread-local. A hit on any of the first four skips the arithmetic it saved and its charge; a plan hit skips the plan's sign decisions. A task on another thread with its own cache would miss, charge more, and could change a refusal. The Boolean's stage caches (face bounds, trim bounds, edge traces, trace projections, point coordinates and surface pairs) have the same effect inside one Boolean.

- **Request-scoped, order-aware caches.** Each entry records the task that produced it. A task later in the order counts it as a hit and charges the reuse cost, waiting for the value if it is not ready or computing it without charge. An earlier task never sees it. Work then depends only on task order, not on the schedule.
- **Deterministic eviction.** `PLANS` keeps 32 entries first in, first out; the enclosure cache clears at 8,192 entries. Eviction by arrival order makes hits depend on the schedule. Keep entries for the whole request, or evict only at stage boundaries.
- **`Decider::root_cache` can stay per decider.** A hit charges the full certificate cost again, so it saves time but never changes work.

## 6. What this does not measure

- Parallelism inside a decider call. Root isolation, enclosure refinement and resultant elimination are single tasks here. Crossing pipes needs that to pass 8.2.
- Rewriting the recipe. A chain of subtractions such as P1's slot cuts could subtract one union of tools, which would turn the chain into a tree. That changes the construction and its proofs, so it is a recipe question, not a scheduling one.
- STEP export, summary and fitted-recipe certification after the build. They run as one serial region each.
- Scheduling overhead, memory traffic and wall time. Work units only; speed is reported, not gated.

## 7. Where the parallel work goes first

The 11 November decision adopts L1, speculative probe searches and unshared caches and facts now. The measurement orders them by effect. Today the longest chain is still 12% (P1) to 55% (P5) of the build's work across the five builds above 1M units even in the Ordered model, and 70% of the median saved input's request, so each step below shortens the chain only where it names:

1. Speculative search over probe candidates and their certification (P1 3.30 to 5.61, seam 2.39 to 3.19), with the per-candidate charging of section 5.2.
2. Order-aware shared caches and commit-ordered facts (sections 5.3 and 5.1), worth up to a further 0.7 on crossing pipes and 0.4 on P1.
3. Parallel root isolation and enclosure refinement for single proofs above about 10^5 units; crossing pipes needs it to pass 8.2.
4. Plan in parallel, append in order, with merges that do not compare every new point against every earlier one by exact proof. On its own it raises P1's bound from 6.02 to 8.17, lowers seam's from 3.50 to 3.01, and moves no other case by more than 0.2.

STEP export, summary and fitted-recipe certification after the build are the next serial regions to divide. Every step keeps the same answer, refusal and work count at every thread count; wall time is reported, not gated.
