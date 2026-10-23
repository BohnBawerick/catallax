# Fillet construction and ownership

Protocol recorded before prototype code. Work bound: one measured round of both candidates against the baseline, with at most one repeat, only if the measurement itself was wrong. No calendar deadline. Decide this prototype before the general-fillet build it informs starts; promotion is a separate implementation on `main`.

## Question and candidates

Which kernel builds the general fillet the user sees? Can a B-rep face built from the offset-face intersection meet the 100 nm target within the 1,000 ms recipe-feature click budget?

- **B-rep first:** construct the rolling-ball surface from its retained definition, trim the input faces and display the resulting B-rep. The field measures deviation against that definition.
- **Field first:** display the same rolling-ball definition through the field; defer B-rep construction to requested STEP export and drawings. Deferred construction must still produce the named fillet, not a different blend.

Both use RECIPE.md section 4.5, including the selected edge, material side, contact arc and end treatment. A smooth minimum or whole-body morphological opening is not a substitute. Neither may report an uncertain or sampled answer as proven. The verification CLI and its driver stay unchanged. Prototype code lives only on local `prototype/fillet-ownership`, never pushed, proposed for merge or merged.

## Fixed cases and baseline

The B-rep phase baseline is `main` at `b05584bf5b967d74c887e43bb0d804b465df2f89`, before candidate code; its executable SHA-256 is `a44170db307fe8df5eb36f5c789112204bd6d92721a7e870cf5b27ec36e69403`. The field phase uses the subsequently landed native-display CLI at `84d453037a49e0aa3fb50697ba485256697e596a`, as recorded below. These phase baselines are not a same-build speed comparison. Earlier observations at `4feea2361fba8d75766df60f181d196aaa98c00e` are preliminary only. Recipe bytes, ordered inputs, hashes, work allowance and selections were frozen before runs. Dimensions below are millimetres. General cases use radius 1, then 3/2 for the edit; primitive poses and edge selections stay fixed.

| Case | Construction and selected edge | Role |
| --- | --- | --- |
| G1 crossed cylinders | Union of radius-8, height-32 cylinders centred at the origin along z and x; round their side/side intersection, using both named wall origins. | General offset-cylinder intersection; decides ownership. |
| G2 oblique cone/plane | Subtract a 60×60×30 block occupying z=-30..0 from a centred frustum with radii 10 and 4, height 24 and axis direction proportional to (3,0,4). Round the exposed cone wall/cutter-end edge. | Oblique conic edge, not the coaxial circular reduction; decides ownership. |
| G3 variable-curvature edge | Union of a radius-12 sphere at the origin and a radius-5, height-40 z-cylinder centred at (4,0,0). Round the sphere-side/cylinder-wall intersection. | Noncircular edge with varying curvature; constant fillet radius, not a variable-radius feature; decides ownership. |
| P3 exact reference | Unchanged `parts/P3-filleted-boss.json`; edit junction radius 2 to 5/2 and top-edge radius 1 to 3/2, retaining IDs and selections. | Exact plane/plane and perpendicular plane/cylinder reductions, built directly in stage 3; cannot decide general ownership. |

Validate that each general selection reaches the intended nonempty input edge before interpreting fillet numbers. A refusing input boolean is an upstream coverage limit, not a fillet failure or an impossible fillet. Keep the failed request and its operation log.

## Measurement and decision

Run the unchanged verification CLI through `.claude/skills/verify-catallax/scripts/control-catallax.py`, with independent `recipe`, `field` and `compare` requests for every original and edited case. Use the same 100,000,000-unit allowance for baseline and candidates. Record refusal codes and latency separately from successful geometry. Keep commands, responses, operation logs and STEP files outside the repository. A missing response or panic is incomplete, never a geometric refusal.

| Required number | Evidence and acceptance |
| --- | --- |
| Time to first display, ms | Fillet construction and tessellation on an already-built input body, starting at the first changed feature. Both candidates use the same immutable input bodies from the frozen prefixes. Input construction is recorded separately and is not charged to the fillet click. Record requested comparison/export costs separately. Feature budget 1,000 ms; complete background comparison budget 2,000 ms. |
| Definition deviation, nm | Field-measured, two-sided bound over every whole fillet face against its own rolling-ball definition. Target 100 nm, hard cap 1 µm, also at most one tenth of the shortest neighbouring edge or narrowest face. Samples, empty marks and STEP curve deviation are not this certificate. |
| Names after radius edit | Resolved selected input origins and resulting blend/corner origins before and after the fixed edit; report retained, lost and ambiguous counts. Matching strings without resolved geometry do not prove identity. |
| STEP import | Requested export opens in OpenCASCADE without repair, has valid topology and matching face count. Record actual counts and export latency. Missing export or bindings is not a pass. |

Use doctor-selected pinned cores at nice 10 with at most two build jobs and two local test threads. Retain ambient machine load. Any changed outcome or claimed speed difference gets alternating baseline/candidate pairs with one warmup and three measured pairs. Run earlier-answer corpus and failure-pack replay, compare complete manifests and strict STEP bytes, and check independent one-, two- and four-worker runs before promotion. Use the application loop when available; explicitly mark missing application observations unmeasured. An unavailable application loop does not waive the required whole-face deviation certificate.

Prerequisite: build the field's two-sided per-face distance to each face's own definition, the measuring stick in PLAN.md section 3.7, on `main` as a separate task and expose it in the verification CLI. Wait for that implementation to land before retaking the baseline or writing candidate code. Then measure the baseline and both candidates with the same unchanged CLI and driver; no supplemental observer is used. Point samples, whole-result boundary distance and display chord error cannot substitute for the required whole-face certificate. The missing measuring stick is a prerequisite, not a measured failure of either candidate and not grounds for a park verdict.

The candidate meeting the click budget and 100 nm target on the general cases owns fillet display. If both do, choose B-rep first because STEP and drawings need it anyway. If neither does, revisit RECIPE.md section 4.5 for interactive use. Missing certificates, valid topology, edit identity or earlier-answer preservation prevent promotion. Verdict is **promote**, **park**, or **repeat once** only when the measurement itself was wrong. P3 numbers remain a reference, never evidence that general fillets work.

Measure B-rep first. If it meets both the click budget and the 100 nm target, field-first display timing cannot change the tie-breaking outcome and is explicitly not measured. If B-rep first fails either criterion, pause for a separate `main` task to expose field-native display timing in the CLI before measuring field first. B-rep mesh time and field compilation time are not substitutes for field-first display time.

Cached-input timing is the experiment boundary. The measured G1 and G3 prefixes already consume more than the click budget, but SolidDraft rebuilds from the first changed feature and keeps earlier bodies valid; input construction is not the operation being compared. The unchanged CLI's fillet operation span starts after its input body exists. Report that span directly. For a completed B-rep display, subtract the preceding input-construction operation spans from `first_mesh_ms` to obtain a conservative cached-input first-display upper bound. Parsing and other unassigned overhead remain in that bound. A bound below budget proves the timing criterion; a bound above budget alone does not prove failure. The fillet span itself supplies a lower bound. Keep full-recipe prefix costs as a separate finding.

The first B-rep construction probe uses the equal-radius perpendicular-cylinder run in G1. Its two offset-cylinder intersections are ellipses. Cover those exact centre paths by rational stations, form an envelope of exact rolling-ball spheres with ordinary B-rep booleans, and trim the added region against the offset supports and input body. The centre-path cover targets 100 nm; it is not a face-deviation certificate. Any completed approximation remains explicitly non-exact with a conservative whole-body distance bound until the unchanged CLI measures its actual face deviation. Stop on the first failed B-rep criterion before expanding the probe to the remaining cases.

Field-display measurement uses the separately landed CLI at named `main` revision `84d453037a49e0aa3fb50697ba485256697e596a`, preserved before field-candidate code with executable SHA-256 `aeae0d1db4afb27f5c78debf901e604fc5ae3407e6e6cb0183644326c31a130b`. The earlier B-rep baseline and observations remain historical, not a same-tool speed comparison. The field probe constructs G1's continuous rolling-ball contact strips directly from the two exact ellipse centre paths, with the remaining cylinder walls and caps trimmed analytically. Its field-native mesh uses parameter interpolation bounds and outward coordinate rounding, not B-rep faces or a sampled maximum. Both radii and all named cases are requested; out-of-domain cases remain explicit coverage refusals. Display certificates do not become whole-face fillet certificates, input-source picking does not become resolved blend naming, and absent deferred STEP does not pass export.

## Results

### Historical input and B-rep observations

The pre-fillet input controls used the preserved B-rep baseline in both slots, core 2, one warmup and three alternating pairs per case. All 12 measured pairs preserved complete manifests and strict STEP bytes. These input costs are outside the cached-input fillet click.

| Input case | Samples | Construction minimum, ms | Median, ms | Maximum, ms |
| --- | ---: | ---: | ---: | ---: |
| G1 | 6 | 1,220.766 | 1,384.829 | 2,229.054 |
| G2 | 6 | 294.528 | 335.370 | 615.676 |
| G3 | 6 | 1,075.773 | 1,191.638 | 1,371.575 |
| P3 reference | 6 | 39.065 | 54.853 | 60.841 |

Priority qualification: these historical launchers requested nice 10 without recording the inferior's actual niceness. The B-rep round nested two `nice -n10` increments, which would reach nice 19 from a normal-priority launcher. Its clocks are not canonical nice-10 click measurements. They remain historical observations, not a speed comparison or proof that general fillets inherently cost minutes.

The pre-code cached-input amendment is `08df463ee0f2b3da488c194221763062d6d1d4a7`; measured B-rep code is `f7bea7fb694383b2ee4731b679dd1cf199c75fa8`, executable SHA-256 `02069343cacc396f27bd38703dea912e56e53171ca532086a9a4395b27b8ccd4`. The CLI, driver and face measuring stick were unchanged from that phase's baseline. All three completed G1 radius-3/2 samples exhausted the 100,000,000-unit allowance: 3,515,134 input units plus 96,484,842 fillet units. No body, display or STEP resulted. The historical fillet-span median was 243,995.007 ms, **until refusal**, with the priority qualification above.

The original-radius warmup was censored by the batch guard, not geometrically refused. The initial B-rep construction stopped at the recorded gate. A diagnostic-symbol replay at actual nice 10 found 12 envelope Boolean unions accounting for 97.7% of its marked wall interval; fitting and independent face certification were not reached. These are limitations of this finite-ball construction, not of every B-rep method.

### Field-first continuous G1 display

The pre-code field protocol amendment is `eeb7855`; measured implementation is `736770ef5110d7909fcbb54e9d09188bdc10409d`, executable SHA-256 `cdb2e3ba7434dce42c706c18059e2ba16062cf45a5d1a07164d8292775e51dc9`. The CLI, driver and independent face measuring stick stayed unchanged from named field baseline `84d453037a49e0aa3fb50697ba485256697e596a`.

One warmup and three alternating baseline/candidate pairs completed for **all eight original/edited inputs**, on core 2 with observed launcher/driver nice 10. Children inherit that priority without nested increments. The shared host's one-minute load was 5.44 at round start and 7.16 at completion; it was not idle.

G1 uses continuous rolling-ball contact strips, analytically trimmed cylinder walls and unchanged caps. Whole-chart derivative bounds and outward-rounded coordinate enclosures certify its native display mesh. No B-rep builder or mesher is used. The table is the conservative cached-input first-mesh upper bound, obtained by subtracting preceding input-construction spans from the CLI clock.

| G1 radius | Minimum, ms | Median, ms | Maximum, ms | Vertices | Triangles | Complete CLI response median, ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 297.877 | 322.948 | 454.773 | 527,432 | 786,432 | 41,817.607 |
| 3/2 | 495.785 | 520.295 | 588.300 | 790,600 | 1,310,720 | 61,338.194 |

**G1 passes the kernel first-mesh timing criterion at both radii.** Complete-response time includes writing the large OBJ; the frozen first-mesh clock excludes that I/O, transport and UI presentation. No application-frame observation was made. Counts, finite coordinates and triangle indices were checked in the actual OBJ files; all three measured repeats had identical strict OBJ bytes for each radius.

The returned bound is `3/400` mm, or **7,500 nm, for display triangles against the defined boundary**. It is not the required 100 nm whole-fillet-face construction certificate. `face_measurements.faces` remains empty and `edit` null. Neither `exact: true` for the retained definition nor an empty face list measures zero construction error.

G2, G3 and P3 remain outside this native constructor's domain. The baseline returns `field.unsupported-operation`; the candidate returns `prototype.field-fillet-domain`. The following are diagnostic **whole-request refusal medians**, not display times or speedups. No mesh is published for these cases; they do not prove impossible geometry.

| Case | Radius / radii | Baseline request, ms | Candidate request, ms |
| --- | --- | ---: | ---: |
| G2 | 1 | 0.693 | 0.663 |
| G2 | 3/2 | 0.938 | 0.686 |
| G3 | 1 | 0.606 | 0.573 |
| G3 | 3/2 | 0.667 | 0.582 |
| P3 reference | 2, 1 | 0.913 | 0.734 |
| P3 reference edit | 5/2, 3/2 | 0.911 | 0.590 |

### Names, deferred STEP and comparison

Native OBJ groups carry six input-source handles, not resolved blend/corner origins. Authoritative point/face queries for this display-only probe explicitly refuse rather than answer with the unfilleted input union. Retained/lost/ambiguous **blend identity counts are unmeasured**, not zero; matching source strings is not an edit-identity pass.

Deferred STEP was actually requested through `recipe` for every original/edited input at the field phase's named baseline and native candidate. G1 exhausted the request allowance at both radii; G2/G3/P3 returned `prototype.fillet-domain`. The single G1 export-contract observations below are not a repeated or paired speed benchmark.

| G1 radius | Fillet span, ms | Whole export request, ms | Accepted request work | Outcome |
| --- | ---: | ---: | ---: | --- |
| 1 | 15,219.234 | 15,921.593 | 99,995,512 | `numbers.operation-budget` |
| 3/2 | 16,626.479 | 17,661.443 | 99,999,976 | `numbers.operation-budget` |

No candidate STEP file exists. OpenCASCADE import, valid imported topology and matching face count therefore do **not** pass and were not attempted on nonexistent files. Both baseline and candidate `compare` requests returned `unresolved` for all eight inputs; that is not agreement or a global certificate.

Four independent pre-fillet control replays against the field phase's baseline preserved complete manifests and strict STEP bytes. Full earlier-answer corpus/failure-pack replay and independent one-/two-/four-worker promotion proofs were not run; no promotion is proposed. Commands, doctor records, binary/input hashes, complete responses, operations, OBJ files, refusal cleanup checks and censored evidence are retained outside the repository.

### Verdict: park

**Park both tested constructions; no fillet-display owner is promoted.** B-rep first produces no fillet within the work allowance. Field first demonstrates sub-budget G1 meshing, but does not cover the other general cases, provide the required whole-face 100 nm measurement, resolve blend identity after the edit, or produce deferred STEP. The decision does not rest on a missing display CLI or on the historical priority error.

The native G1 result is a useful latency reference, not evidence that field first satisfies all ownership gates. Revisit the interactive construction/representation choice in RECIPE.md section 4.5 before the general-fillet build, without weakening the recipe definition or proof requirements. A later trial can test compact contact-curve/patch construction with certified local bounds. Stage-3 exact plane/cylinder fillets remain direct builds and P3 remains a reference only.

No measurement repeat is requested: deterministic construction/coverage/export failures and missing edit identity already prevent promotion; the qualified historical clocks are not used as a canonical timing decision. The local `prototype/fillet-ownership` branch is never pushed or merged. Only this protocol/results file is delivered; promotion would require a separate proper implementation on `main`.
