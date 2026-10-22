# Catallax roadmap

Written 22 October 2026. [PRD.md](PRD.md) defines the product. This roadmap replaces the previous implementation direction and orders the work. It was written after reading Hanno Braun's retrospective ["Shutting Down Fornjot"](https://archive.hannobraun.com/fornjot/blog/shutting-down-fornjot/) and measuring the project against it, lesson by lesson. Where the old plan conflicted with a lesson, this plan pivots. Where a recommendation would change the PRD, it is listed under [Proposed PRD changes](#proposed-prd-changes) and waits for the project owner.

The companion documents are unchanged: [ARCHITECTURE.md](ARCHITECTURE.md) and [DESIGN.md](DESIGN.md) govern the B-rep implementation, [RECIPE.md](RECIPE.md) defines the shared language, [FIELD.md](FIELD.md) the field evaluator, and [COMPARE.md](COMPARE.md) the two-way display comparison and its measurements, which used to live in this file.

## 1. The premise

Catallax is the kernel for SolidDraft, and SolidDraft is the application that runs on Catallax. Kernel work is chosen by what SolidDraft's target parts and features need, in the order its users need them. Generalising the kernel for other callers comes later, once it carries one application well. The two stay separate repositories with separate releases.

This is a change of stance, not of mathematics. The old plan aimed at "reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels" and put SolidDraft integration last, after every operation was covered. Fornjot's author cut his application to focus on the kernel and watched the vision blur; he now says a kernel can still be a tool for specific use cases rather than a building block for everyone. We take that at face value. The kernel's target is five parts rebuilding in one application at click speed with no wrong answer, not parity with another kernel.

SolidDraft, as its own README describes it, is a SolidWorks-style desktop CAD application in Rust with a feature tree, an AI assistant, a ribbon and a 3D viewport. It runs today on NVIDIA's SMLib kernel through an interface the application owns, sketches on planes and faces, extrudes, cuts, fillets and chamfers, makes drawings, and has a headless command line and a browser build. Catallax is not plugged in. The project owner's direction is that SolidDraft becomes the application that runs only on Catallax, with Catallax designed for the features its users drive. What that means for SolidDraft's own code is decided in its repository, not here.

Everything kept from the PRD and the design rules still holds: the same question gives the same answer on every machine and at every thread count; no unsafe code; derived geometry stored by definition; a proven answer is never wrong; only impossible geometry is an acceptable final refusal.

## 2. Fornjot's lessons, and what each one changes here

Hanno Braun worked on Fornjot from July 2020 to June 2026 and did not reach a usable kernel. His retrospective names six mistakes. The table says where Catallax stood on each one on 22 October 2026 and what this roadmap does about it.

| Lesson | Where Catallax stood | What changes |
|---|---|---|
| Early progress does not extrapolate. Primitives and transversal intersections come fast; tangent, coincident and tolerance-driven cases are a cliff. | Exact arithmetic and certified decisions avoid wrong answers. The cliff reappears as speed and refusals: cylinder pairs with seams or equal radii, and general cone and torus poses, mostly refuse at the 50-second deadline. Extrude and revolve came in a day; fillets, shell and offset are the cliff ahead. | A cliff corpus becomes a standing gauge, measured and published at every checkpoint. Progress is counted as target parts that build at click speed, never as operations that work on easy poses. The arithmetic design itself is put to a measured test (the floating-point filter prototype) instead of being assumed. |
| Funding before a product sells a dream. | The kernel site named topology optimisation and instant FEA beside live sponsor links. | A works-today list stands beside every sponsorship ask, built from `main`, and the pull request that changes a capability updates it. Aspirations stay on the site only with an "in design" label. Each checkpoint publishes a dated note of what changed. |
| Incremental change alone cannot fix an architecture; both extremes are wrong. Prototype constantly, integrate only what proved itself. | Designs were written on paper by the architect and implemented straight into the main code. Speed was being ground down refusal by refusal. The shared algebraic foundation is a large redesign that must land whole. | A prototype lane with a fixed shape: question, measurement, time box and a decision rule written before the spike starts. Spike code is never merged; a winning idea is re-implemented under the design rules. Two spikes start now. |
| Half-measures in a crisis. Each time, the options were full commitment or quitting, and he took neither. | AI workers keep the code moving without the owner's hours. The analogous risk is direction churn: three kernel efforts in about a week before Catallax. | Dated commit-or-park checkpoints every two weeks. At each one, every item in flight is either committed, with a landing date before the next checkpoint, or parked with a written state. A change of kernel direction happens only at a checkpoint and only with a record. |
| A generic kernel is a vaguer target than an application. | The target was parity with other kernels; SolidDraft integration was the last item. | Section 1. SolidDraft's target parts and features choose the kernel work. Integration moves from last to first. |
| Prototyping came too late and the redesign never landed; the codebase sat half-migrated. | One representation-wide migration is in flight (the algebraic foundation) with four lanes waiting to rebase onto it. | One migration at a time. The foundation lands whole by the first checkpoint or is parked. No lane introduces a competing representation. A lane that cannot rebase by its checkpoint is parked, not carried. |

## 3. What "useful" means: target users and target parts

The first users are SolidWorks users making machined and printed parts, which is who SolidDraft is built for. For them a kernel is useful when the part they are drawing rebuilds after every click, exports a STEP file their shop can open, and never shows a wrong body.

Five target parts define that. They are chosen so that each adds one capability the previous one did not need, and together they cover what a SolidWorks user does in a first session. Each is written as a `catallax-recipe/3` file in a `parts/` directory the next items create, with its dimensions fixed so timings compare across revisions.

| Part | What it exercises | Needs beyond today |
|---|---|---|
| P1 Bracket. A plate with a flange along one edge, a boss with a through bore, four bolt holes on a rectangular pattern and a slot through the flange. | Extrusions and booleans only, including a union whose operands share a face, and a chain of about a dozen subtractions. | A sound coplanar-face boolean in both kernels. Rebuild under the click budget. |
| P2 Turned collar. A rectangle revolved a full turn into a ring, chamfered on both outer edges. | Revolution, torus and cone surfaces, chamfer on circular edges. | Chamfer. |
| P3 Filleted boss. A block with a cylindrical boss on top, the junction edge filleted, the block's top edges filleted. | Plane-to-plane and plane-to-cylinder fillets, concave and convex, with ends that meet transversally. | Fillet, per RECIPE.md section 4.5, with the end and corner rules of section 7.2 settled for these cases. |
| P4 Shelled tray. A block shelled to a thin wall with its top face open, then intersected with a cylinder to round its corners. | Shell with an open face, then a boolean on the shelled body. | Shell, per section 4.7, with section 7.3 settled. |
| P5 Flanged housing. A revolved housing with a bolt-hole pattern of countersunk holes and fillets where the flange meets the body. | Pattern, countersink cones in booleans, fillet against a revolved face. | Pattern, hole and countersink, fillet on a torus or cone neighbour. |

Beside them stands P0, the cliff corpus: the kernel's own hard cases, kept small and fixed. Equal-diameter pipes crossing at right angles, two cylinders tangent along a line, a cylinder and sphere with equal radii, a block and sphere union, poses along a seam, and the general cone and torus poses that refuse today. P0 is a gauge, not a target. It measures the cliff; nobody is asked to flatten it on a schedule.

A target part counts as built when, on the reference machine (one core, release build, `nice 10`, no other load):

- the B-rep builds it from the recipe with no refusal and no mark above the geometric limit;
- the field evaluator builds the same recipe and the two-way comparison agrees at display accuracy, or the recipe uses an operation the field does not yet build and that is stated;
- a full rebuild finishes inside the click budget below;
- its STEP file opens in OpenCASCADE with the same face count and a volume within the exported deviation;
- a parameter edit (plate thickness, boss radius, fillet radius) rebuilds with every recipe name still resolving to the same face.

The click budget is the one the verification exchange records: 100 ms for a stock solid, 500 ms for one boolean, 1,000 ms for a complete recipe, 2,000 ms for the complete sampled comparison. These are the first numbers, not a proof of what the design can reach. P1 is the first part to meet them. Whether it does by the second checkpoint is the first real test of exact-first arithmetic as the interactive path, and the floating-point filter prototype exists to give that test a second candidate.

## 4. Process

### 4.1 Prototype lanes

A prototype answers one question that a paper design cannot. Every spike is set up in a short file under `docs/spikes/<name>.md` before any code is written, with these fields:

- the question, in one sentence;
- the measurement: which target parts and which P0 cases, which numbers come out, and which tool measures them. The verification exchange is the default instrument; a spike that needs another one says so;
- the baseline: the same measurement on `main` at a named revision;
- the decision rule: the numbers that promote the idea, written before the spike starts;
- the time box, in checkpoints. One checkpoint by default, two at most.

Spike code lives on a `spike/<name>` branch and is never merged. It may break every design rule except two: it may not change the measurement tool, and it may not report a wrong answer as a right one. When the time box ends the spike file gets its result section, the branch is left where it stopped, and the decision is one of three:

- promote: a build item re-implements the idea on `main` under the design rules, with the spike's measurements as its exit check;
- park: the idea stays recorded with its numbers, no further work;
- repeat once: only when the measurement itself was wrong, with a new decision rule.

Two spikes start now.

**Float filter.** Question: does a floating-point filter in front of the exact decisions, with a certified error bound that falls back to exact arithmetic whenever the filter cannot decide, bring P1 inside the click budget and shorten the P0 refusals, without any decision ever being made by the filter alone? Measured on P1 and P0: wall time, work units, refusal count, and zero changes in any proven answer. Baseline: `main` after the algebraic foundation lands. Promotion rule: P1 under 1,000 ms and no P0 case slower than before. If the filter does not reach that, the result still answers the question nobody has asked yet, whether the arithmetic design is the bottleneck.

**Fillet, B-rep-first or field-first.** Question: for P3's fillets, which kernel should build the fillet the user sees? The B-rep candidate constructs the rolling-ball surface by definition with analytic contact curves. The field candidate evaluates the definition as a distance expression and the B-rep imports only at export. Measured on P3: time to first display, exactness against the section 4.5 definition at sampled points, names surviving an edit of the fillet radius, and STEP export of the result. Promotion rule: the candidate that meets the click budget and the definition on P3 owns fillet display; if both do, B-rep, because STEP needs it anyway; if neither does, the fillet definition is wrong for interactive use and section 4.5 is revisited.

### 4.2 Checkpoints

Checkpoints fall every two weeks. At each one a short dated note is appended to this file: what landed, what parked, the P0 and target-part numbers, and what the next checkpoint must show. The same note goes to sponsors. A checkpoint can be missed by a day for a build to finish; it is not moved for a design argument.

Every item in flight gets one of two verdicts at a checkpoint. Committed means merged, or a named landing date before the next checkpoint and a rebase already done. Parked means the branch stays, a state file under `docs/parked/<name>.md` says what works, what does not and what it would take to resume, and no further work happens until a later checkpoint un-parks it. There is no third verdict. An item that is "almost there" twice is parked.

A change of kernel direction, meaning a new representation, a new arithmetic path as the default, or a new kernel effort, is made only at a checkpoint, only with a spike result or a measured failure as its reason, and only by a written decision in ARCHITECTURE.md. This is the rule against direction churn.

| Checkpoint | Date | Must show |
|---|---|---|
| C1 | 5 November 2026 | The algebraic foundation merged with the rational fixtures byte-identical, or parked. The verification exchange merged. The `parts/` corpus with P0 to P5 recipes and a scoreboard report that builds, compares and times them. The two spike files written and both spikes started. Site aspirations relabelled. |
| C2 | 19 November 2026 | P1 builds inside the click budget with the coplanar-face boolean sound, or the float-filter spike result says why not and what the design change is. Catallax answers SolidDraft's kernel interface for P1's operations, and P1 displays in the application through it. Both spikes decided. |
| C3 | 3 December 2026 | P2 and P3 built: chamfer and fillet on the B-rep as the fillet spike decided. STEP of P1 to P3 opens in OpenCASCADE. The direction question answered in writing: exact-first with a filter, or a changed default path. |
| C4 | 17 December 2026 | P4 built: shell with section 7.3 settled. Names survive the parameter edit on P1 to P4. The error-policy certificate implemented for the operations P1 to P4 use. |
| C5 | 14 January 2027 | P5 built: pattern, hole, countersink. All five parts inside the click budget on the reference machine. The P0 gauge republished. The next five parts chosen with SolidDraft's users. |

### 4.3 Landing work in flight without a half-migrated codebase

ARCHITECTURE.md already sets the merge order for the algebraic foundation: land the decision, land the conic preparation, land one coordinated foundation with every consumer adapted and the rational subset byte-identical, then rebase each lane once. This roadmap adds the rules that keep that from becoming Fornjot's transitional state:

1. One representation-wide migration at a time. The foundation is that migration until it merges or parks. No second one starts before C3.
2. The foundation merges only whole. Its gate is the existing fixture set: every successful rational exchange fixture keeps byte-identical STEP and manifest content, and every enabled test passes on GitHub's runners. If that cannot be shown by C1 the branch parks and `main` stays rational; the lanes then land on rational `main` and the foundation returns as a spike question.
3. Each lane rebases exactly once, onto the merged base, and lands inside one checkpoint of that rebase or parks.
4. A lane that needs a type the base does not have asks for a recorded decision; it does not add a parallel type.
5. Independent fixes that do not widen an interface may land before the base and are folded into it.
6. Every merge updates the in-flight ledger below. The ledger is the single list of what is open; a branch not in it is parked by definition.

### 4.4 Sponsorship tied to shipped capability

- A works-today list stands beside every sponsorship ask, in README.md and on the site. It names only what `main` does, with a "not yet" line for what it does not.
- The pull request that adds, removes or materially speeds up a supported operation updates both lists in the same change. A reviewer checks for it.
- The site may name an aspiration only with an "in design" label, and never as a capability. The field kernel card and the manifesto's instant-FEA and topology-optimisation lines are relabelled as a next item.
- The checkpoint note is the sponsor update: dated, with the numbers, including the ones that went the wrong way.
- No new funding ask goes out between checkpoints. The question "would you pay for what works today" is answered by the list, not by the roadmap.

### 4.5 How work is produced

Most Catallax code is written by AI coding agents, directed and checked by a human. The process rules above are written so that an orchestrator can dispatch them: a spike is a lane with a file, a time box and a decision rule; a checkpoint is a dated review of a ledger; a target part is a recipe and a scoreboard row. The agents' throughput makes two-week checkpoints realistic where a single developer needed months. It also makes direction churn cheap, which is why the checkpoint rule is strict.

## 5. Verdicts on the current plan

Each item of the previous plan, each piece of work in flight, and each queued item gets one of three verdicts: stays, changes, or dropped. Dropped means parked under section 4.2, with its state recorded; nothing is deleted.

### 5.1 The previous plan's items

| Previous item | Verdict | Reason |
|---|---|---|
| 1. Continue independent testing of the B-rep half and repair refusals of possible geometry | Changes | The independent tester and the exchange contract stay. Refusal repair is no longer a sweep over a corpus; a refusal is fixed when a target part or a spike metric hits it. The P0 gauge measures the rest and is published, not scheduled. |
| 2. Settle RECIPE.md section 7 and widen exact geometry beyond rational frames | Changes | Only the open items the target parts need are settled now: 7.2 fillet ends and corners for P3's cases, 7.3 shell for P4, 7.9 feature patterns and 7.8 one-sided through-all for P5. The rest park. Widening exact geometry happens through the algebraic foundation and stops there until a target part needs a new kind. |
| 3. Define and enforce the composable geometric error policy | Changes | The policy is defined and stays the contract. Enforcement is staged: piece 1's certificate types land with the foundation; whole-result certificates are required for the operations P1 to P4 use by C4; the full producer and consumer audit follows P5. It does not gate click speed work. |
| 4. Extend field classification of degenerate zeros and finite closest-point bounds before claiming complete queries | Changes | The coplanar-face part of it is P1's blocker and moves to the front. The rest waits for the fillet spike, which decides how much the field does interactively. The field does not need to match the B-rep operation for operation before P1 reaches the application. |
| 5. Finish the section 7.14 geometry bugs in both evaluators, then holes, fillets, chamfers, shell, offset, draft, sweep, loft, pattern, mirror | Changes | Reordered by target part: chamfer (P2), fillet (P3), shell (P4), hole, countersink and pattern (P5). Offset, draft, sweep, loft and mirror park until a chosen part needs them. The 7.14 algebraic items land with the foundation; elliptic cylinders and spindle tori park unless a part produces them. |
| 6. Integrate display meshes, then STEP import, conversions, hand-off checks, export accuracy levels and stand-ins | Changes | Display meshes with shared edges and face origins are what the application consumes and stay at the front. STEP import, on-demand conversions, export accuracy levels and the large-part stand-in park; the application exports nothing but drawings today and STEP export is the kernel's first output. |
| 7. Integrate into SolidDraft once the required operations are covered | Pivot | Integration moves to C2, with P1's operations only. The kernel is designed against what the application calls. Field-only features, printing, FEA queries and interference checks park. |
| Open policy conflict: the 50-second deadline against the same-answer rule | Changes | Resolved by direction: the verification exchange's deterministic work allowances replace the deadline as the outcome rule; wall-clock becomes a cancellation signal that can only stop work, never choose an answer. The legacy command keeps its contract until a versioned change. |
| Measure useful latency | Stays, promoted | It is now the central metric. The scoreboard reports it per target part at every checkpoint. |
| Background comparison prerequisites and whole-check measurements | Stays, moved | Documentation of `catallax-compare`, not plan. Moved verbatim to [COMPARE.md](COMPARE.md). |
| Evidence for completion | Stays | Mathematical cases, independent oracles, scale and placement variations, repeatability and measured timings remain the evidence. Target parts are added to the list. |

### 5.2 Work in flight

| Item | State on 22 October 2026 | Verdict | Reason and condition |
|---|---|---|---|
| Shared algebraic foundation | Branch with three commits, about 144 files and ten thousand added lines, one commit behind `main` | Stays, lands first | It is the common base every other lane is told to rebase onto, and it already carries the right gate: the rational subset byte-identical. It is also the project's one Fornjot-shaped risk, a redesign that must land whole. Condition: merged by C1 with the fixture gate green on GitHub's runners, or parked with `main` staying rational. |
| Refusal fix batch 5 | Branch with six performance and contact-proof commits, two commits behind `main` | Changes | Land the commits that rebase cleanly onto the base and improve a P0 or P1 number; drop the corpus-driven ordering. Later batches are chosen by target parts and by the float-filter result, not by sweeping refusals. Condition: rebased and measured by C1. |
| New verification exchange and command line (pull request 12) | Open, checks green, conflicts with `main`, waiting on its own full-corpus evidence | Stays, becomes the instrument | Deterministic work allowances replace the deadline, replay caching and timing reports are what spikes and checkpoints measure with, and its interactive targets are the click budget. Condition: rebase and merge by C1 with its own tests; the full-corpus evidence runs as C1 checkpoint data rather than blocking the merge. |
| Recipe open-item definitions | In progress on paper | Changes | Settle 7.2 for P3's cases, 7.3 for P4, 7.8 and 7.9 for P5, each as a definition both kernels can implement, each before its build item starts. Park 7.1, 7.4 to 7.7 and 7.10 to 7.13 with their options recorded. |

### 5.3 Queued items

| Item | Verdict | Reason and condition |
|---|---|---|
| Display meshes for new kinds (elliptic cylinders, open conics) | Changes | After the foundation lands, and only for kinds a target part produces. The application needs edge polylines and face roles from the mesh more than new kinds; those come first. |
| Primitive and feature name equivalence | Stays | The application's persistent references depend on it, and the target-part exit check (an edit with every name still resolving) tests it directly. Lands by C4. |
| Comparison-found gaps | Changes | Split. The coplanar subtraction indeterminate query is P1's blocker and moves to the front. The same-body boolean assertion is small and lands with it. The torus meshing budget parks. |
| Error-policy enforcement | Changes | Staged as in 5.1 item 3: types with the foundation, certificates for P1 to P4 operations by C4, the full audit after P5. |
| Codex cloud offload | Stays, outside the product roadmap | Capacity, not product. Spikes are a good fit for offloaded workers because their code is throwaway and their decision rule is written down. It changes nothing in this plan's order. |

### 5.4 Parked now

STEP import. Export accuracy levels and the large-part stand-in. Field-only features: lattice, gyroid, graded density, smooth blend, mesh and scan import. FEA and interference queries. Sweep, loft, offset, draft, thicken, mirror. Elliptic cylinders and spindle torus sheets unless a target part needs them. Torus meshing allowance. RECIPE.md items 7.1, 7.4 to 7.7 and 7.10 to 7.13. Jev. Each keeps its definition and its refusal code; none is declared impossible.

## 6. Next items

Numbered in priority order. Each names its deliverable, its exit check and whether it is a prototype or a build. The project owner picks from this list at each reconvene; items 1 to 3 are the C1 set.

1. **Land the shared algebraic foundation.** Build. Deliverable: the foundation branch merged onto `main` with every consumer adapted. Exit: every successful rational exchange fixture keeps byte-identical STEP and manifest content excluding elapsed time; every enabled test passes on GitHub's runners; each of the four waiting lanes has rebased or parked. By C1.
2. **Land the verification exchange.** Build. Deliverable: pull request 12 rebased and merged, with its stream workers, deterministic request allowances, replay cache and timing report. Exit: a block question that refuses through the legacy deadline answers through the new exchange; one, two and four workers give identical outcomes and work counts on the existing matrix; the legacy commands and manifest are unchanged. By C1.
3. **Target-part corpus and scoreboard.** Build. Deliverable: `parts/` with P0 to P5 as `catallax-recipe/3` files and a `parts` command in the verification exchange that builds each part on both kernels, compares, times and writes one report. Exit: the report runs on `main` and lists, per part, build time, refusal or marks, comparison outcome and STEP deviation; the first report is committed under `docs/checkpoints/`. By C1.
4. **Float filter.** Prototype. Deliverable: `docs/spikes/float-filter.md` with the question, measurement, baseline and decision rule, and the `spike/float-filter` branch. Exit: the result section filled in with P1 and P0 numbers against the post-foundation baseline and one of promote, park or repeat. Decided by C2.
5. **Fillet, B-rep-first or field-first.** Prototype. Deliverable: `docs/spikes/fillet-owner.md` and the `spike/fillet-owner` branch with both candidates on P3. Exit: time to display, agreement with the section 4.5 definition, names after a radius edit and STEP export for each candidate, and the ownership decision. Decided by C2.
6. **Coplanar-face booleans.** Build. Deliverable: a sound union, subtraction and intersection when operands share a face or a face lies in another's plane, in both kernels, with the field's degenerate zero classified rather than indeterminate. Exit: P1's flange-on-plate union and the two-block reproducers in COMPARE.md agree with no unresolved observation; the same-body boolean reaches the set rule, not the assertion. By C2.
7. **Catallax answers SolidDraft's kernel interface for P1.** Build. Deliverable: on the Catallax side, a converter from the application's recipe to `catallax-recipe/3` following RECIPE.md section 6.2, and a library entry that builds a recipe and returns the body, its display mesh with edge polylines, mass properties and per-step face origins in the shapes the application's interface asks for. Exit: P1 built through that entry displays in the application through a backend for its interface, with faces picked and named; the backend itself lives in the application's repository. By C2.
8. **Chamfer.** Build. Deliverable: RECIPE.md section 4.6 implemented on the B-rep for plane-plane and circular edges, with names. Exit: P2 builds inside the click budget; STEP opens in OpenCASCADE; the field builds the same recipe or states the gap. By C3.
9. **Fillet.** Build. Deliverable: section 4.5 implemented where the fillet spike put it, for plane-plane and plane-cylinder edges with transversal ends, and section 7.2 settled for those cases. Exit: P3 builds inside the click budget; a radius edit keeps every name; STEP opens in OpenCASCADE. By C3.
10. **Deadline out of outcomes.** Build. Deliverable: library paths take a cancellation signal that can only stop work; no wall-clock value selects an answer, a mark or a refusal in any library call; the legacy command keeps `numbers.deadline` under its contract, with the versioned exchange decision recorded in INTERFACE-CHANGES.md. Exit: the same question under load and idle gives the same outcome through the library; the independent tester's contract is unchanged. By C3.
11. **Shell.** Build. Deliverable: section 4.7 implemented with section 7.3 settled, for planar and cylindrical walls with one open face. Exit: P4 builds inside the click budget and compares. By C4.
12. **Name equivalence and edit survival.** Build. Deliverable: primitive and feature names equivalent across both kernels, and the parameter-edit check in the scoreboard. Exit: the edit check passes on P1 to P4 with every reference resolving to the same face. By C4.
13. **Error-policy certificates for P1 to P4 operations.** Build. Deliverable: piece 1's certificate types, the two reserved refusal codes, and whole-result certificates for extrude, revolve, boolean, chamfer, fillet and shell. Exit: P1 to P4 publish a certificate below 1 nm or refuse with the specific code; a forged mark is rejected. By C4.
14. **Hole, countersink and pattern.** Build. Deliverable: sections 4.3 and 4.13 implemented, with 7.8 and 7.9 settled for P5. Exit: P5 builds inside the click budget; the pattern's names resolve per copy. By C5.
15. **Site and README relabel.** Build. Deliverable: the field kernel card and the manifesto name their aspirations under an "in design" label; the SolidDraft card says what its README says; the works-today lists are the ones a reviewer checks on every capability change. Exit: no capability is named on the site that `main` does not have. By C1.
16. **Checkpoint notes.** Build, recurring. Deliverable: a dated note under `docs/checkpoints/` at each checkpoint with the ledger verdicts and the scoreboard. Exit: the note exists on the date and the same text went to sponsors.

## 7. In-flight ledger

Updated at every merge and every checkpoint. A branch not listed here is parked.

| Lane | Branch | Verdict at C0 (22 October 2026) | Lands by |
|---|---|---|---|
| Algebraic foundation | `fm/catallax-algebraic-foundation` | Committed | C1 |
| Verification exchange | `fm/catallax-verify-cli`, pull request 12 | Committed | C1 |
| Refusal fix batch 5 | `fm/catallax-refusal-batch5` | Committed, measured subset | C1 |
| Recipe open items 7.2, 7.3, 7.8, 7.9 | none yet | Committed | each before its build item |
| Statement and works-today list | `fm/catallax-roadmap-fornjot-statement`, pull request 15 | Committed | C1 |
| This roadmap | `fm/catallax-roadmap-fornjot` | Waits for the project owner | C1 |

## 8. Proposed PRD changes

These recommendations change the PRD's direction. They are proposals; docs/PRD.md is unchanged until the project owner decides.

1. **Section 1, what Catallax is.** "It is the kernel for SolidDraft and for anyone else who wants it" becomes "It is the kernel for SolidDraft, designed for the features SolidDraft's users drive. Anyone else may use it; generalising it for them comes after it carries one application well."
2. **Section 3.8, target.** Replace "reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels" with "the target parts in section 3 of the roadmap build in SolidDraft at click speed with no wrong answer. Comparisons with other kernels are measurements we publish, not the target."
3. **Section 10, direction in order.** Move "Integrate Catallax into SolidDraft" from step 6 to step 2, limited to the operations the first target part needs, and make "Build the field kernel for everything the B-rep half already builds" follow the target parts rather than precede widening. The field keeps its role as the second opinion in the background check and as the owner of field features; it stops being a gate on each B-rep operation.
4. **Section 3.2, refusals.** Keep the rule. Add one sentence: a refusal that arrives after the click budget is a bug twice, once for refusing and once for taking that long, and refusal latency is reported with every refusal.
5. **Section 3.5, checks.** Keep all three. State that the in-app background check runs only for operations both kernels build, and that a part whose recipe uses a B-rep-only operation is still a valid part with that fact recorded, not a rebuild error.
6. **Section 4, features table.** Reorder by target part and mark each row with the checkpoint that first needs it. The table's claims about the field kernel for fillets, shells and sweeps become "to be decided by the fillet spike" until it has run.
7. **Section 9.4, SolidDraft.** "The kernel remains usable independently of that application" stays true as a licence fact. Add the direction that SolidDraft runs on Catallax only, and that SolidDraft's interface and recipe format converge on RECIPE.md through the mapping in its section 6 rather than staying a swappable slot.
8. **Section 11, open items.** Close the deadline item by the direction in next item 10, and add: the click budget numbers; which five parts follow P1 to P5.

## 9. Open questions for the project owner

Nothing below blocks items 1 to 3. Each is a decision to take at the reconvene.

1. **The five target parts.** Section 3 proposes P1 to P5. Confirm, swap or add; the recipes are written after the choice. A sixth candidate, a sheet-metal style bent plate, was left out because it needs sweep.
2. **Click budget.** Adopt the verification exchange's 100, 500, 1,000 and 2,000 ms targets as the budget, or set different ones. They are stated as first numbers, not as proven reachable.
3. **Field kernel parity.** Proposed PRD change 3 defers field parity until the target parts build. The alternative keeps the PRD order and accepts that P1 in the application waits for the field to cover every operation. The roadmap assumes the deferral.
4. **Sponsorship during the pivot.** Keep the sponsor links live with the works-today list beside them, as this roadmap assumes, or pause asks until P1 builds in the application. Fornjot's author would not sell a dream again "without an actual product to deliver alongside"; the list is our product-alongside.
5. **Checkpoint cadence.** Two weeks is proposed because agent throughput makes it realistic. Four weeks is the alternative if the owner's review time is the constraint.
6. **Spike ownership.** Whether a spike may be run by a worker with no prior context of the repository, which tests the design against a fresh reading, or only by workers who know the design rules.

## Evidence for completion

Use mathematical cases, independent numeric oracles, scale and placement variations, repeatability checks, measured timings and the target-part scoreboard. Keep refusals, marked answers and wrong proofs distinct. An independent tester built without seeing the kernel's code uses the documented exchange contract. No unpublished tester design or acceptance system forms part of this plan.
