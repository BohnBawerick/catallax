# Fillet construction and ownership

Protocol recorded before prototype code. Work bound: one measured round of both candidates against the baseline, with at most one repeat, only if the measurement itself was wrong. No calendar deadline. Decide this prototype before the general-fillet build it informs starts; promotion is a separate implementation on `main`.

## Question and candidates

Which kernel builds the general fillet the user sees? Can a B-rep face built from the offset-face intersection meet the 100 nm target within the 1,000 ms recipe-feature click budget?

- **B-rep first:** construct the rolling-ball surface from its retained definition, trim the input faces and display the resulting B-rep. The field measures deviation against that definition.
- **Field first:** display the same rolling-ball definition through the field; defer B-rep construction to requested STEP export and drawings. Deferred construction must still produce the named fillet, not a different blend.

Both use RECIPE.md section 4.5, including the selected edge, material side, contact arc and end treatment. A smooth minimum or whole-body morphological opening is not a substitute. Neither may report an uncertain or sampled answer as proven. The verification CLI and its driver stay unchanged. Prototype code lives only on local `prototype/fillet-ownership`, never pushed, proposed for merge or merged.

## Fixed cases and baseline

Preliminary coverage observations use `main` at `4feea2361fba8d75766df60f181d196aaa98c00e`; retain that executable and its evidence, but do not use them as the candidate round's baseline. Retake the baseline on `main` after the measuring-stick prerequisite below lands, record its full revision here, and preserve its executable, SHA-256 and doctor output outside the repository before candidate code changes. Freeze recipe bytes, ordered inputs, hashes, work allowance and selected origin names before candidate runs. Dimensions below are millimetres. General cases use radius 1, then 3/2 for the edit; primitive poses and edge selections stay fixed.

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
| Time to first display, ms | First completed visible fillet, including required construction and tessellation. Record recipe-build and requested comparison/export costs separately. CLI response time without a rendered fillet is not first display. Feature budget 1,000 ms; complete background comparison budget 2,000 ms. |
| Definition deviation, nm | Field-measured, two-sided bound over every whole fillet face against its own rolling-ball definition. Target 100 nm, hard cap 1 µm, also at most one tenth of the shortest neighbouring edge or narrowest face. Samples, empty marks and STEP curve deviation are not this certificate. |
| Names after radius edit | Resolved selected input origins and resulting blend/corner origins before and after the fixed edit; report retained, lost and ambiguous counts. Matching strings without resolved geometry do not prove identity. |
| STEP import | Requested export opens in OpenCASCADE without repair, has valid topology and matching face count. Record actual counts and export latency. Missing export or bindings is not a pass. |

Use doctor-selected pinned cores at nice 10 with at most two build jobs and two local test threads. Retain ambient machine load. Any changed outcome or claimed speed difference gets alternating baseline/candidate pairs with one warmup and three measured pairs. Run earlier-answer corpus and failure-pack replay, compare complete manifests and strict STEP bytes, and check independent one-, two- and four-worker runs before promotion. Use the application loop when available; explicitly mark missing application observations unmeasured. An unavailable application loop does not waive the required whole-face deviation certificate.

Prerequisite: build the field's two-sided per-face distance to each face's own definition, the measuring stick in PLAN.md section 3.7, on `main` as a separate task and expose it in the verification CLI. Wait for that implementation to land before retaking the baseline or writing candidate code. Then measure the baseline and both candidates with the same unchanged CLI and driver; no supplemental observer is used. Point samples, whole-result boundary distance and display chord error cannot substitute for the required whole-face certificate. The missing measuring stick is a prerequisite, not a measured failure of either candidate and not grounds for a park verdict.

The candidate meeting the click budget and 100 nm target on the general cases owns fillet display. If both do, choose B-rep first because STEP and drawings need it anyway. If neither does, revisit RECIPE.md section 4.5 for interactive use. Missing certificates, valid topology, edit identity or earlier-answer preservation prevent promotion. Verdict is **promote**, **park**, or **repeat once** only when the measurement itself was wrong. P3 numbers remain a reference, never evidence that general fillets work.

## Results

Not measured yet. No ownership decision or promotion is claimed by this protocol.
