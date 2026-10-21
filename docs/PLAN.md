# Catallax implementation direction

[PRD.md](PRD.md) defines the product. This plan separates the inherited B-rep implementation from the work still required. The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels. The former 1-in-10,000 programme and grid-cell acceptance goals are retired; keep their useful mathematical regressions.

## Keep the B-rep starting code

Carry forward exact arithmetic and certified decisions with deterministic work budgets, definition-based analytic geometry, append-only curved topology, the five-stage boolean engine, bounded mass properties and STEP AP242 export. Preserve the public exchange contract unchanged, including its legacy names. `catallax-features` adds named B-rep recipe stock solids and rationally represented sketch features. The independent `catallax-field` evaluator covers stock solids and booleans; see [FIELD.md](FIELD.md) for its guarantees.

The existing tests include ignored cases. Their source records the reasons; a green default suite does not mean those cases were repaired. Existing work still in the source repository must arrive as a separately reviewed export delta, never by silently replacing this tree.

## Next work, in order

1. Continue independent testing of the B-rep half and repair refusals of possible geometry. General cone and torus pairs, ambiguous imprints and classification, precision limits and missed measurement accuracy remain unfinished work. Only impossible geometry is an acceptable final refusal.
2. [RECIPE.md](RECIPE.md) defines each operation and its origin-based naming. Named B-rep recipe builds now preserve origins through splits and maximal merges. Settle the remaining section 7 items and widen exact geometry beyond rational analytic frames. Arena indices are not persistent names.
3. Define and enforce a composable geometric error policy. A bounded decision may stand in for exact coincidence only after proving the gap smaller than 1 nm. Existing scalar, angle and distance marks do not yet establish that invariant over a whole result.
4. The field evaluator now covers recipe blocks, cylinders, spheres, cone frusta, tori and booleans, with named face candidates, exact predicates, certified primitive distances and conservative Boolean distance bounds. Extend certified classification of degenerate regularised zeros and finite closest-point bounds on surviving trimmed faces before claiming complete distance/proximity queries. Add development agreement checks and the background comparison only when the prerequisites below exist.
5. Finish the B-rep sketch, extrusion and revolution geometry bugs listed in RECIPE.md section 7.14, then bring those features to the field evaluator and compare their sets and names. Continue through holes, fillets and chamfers, shell and offset, draft, sweep and loft, pattern and mirror. Every operation needs its shared definition before two implementations can agree.
6. Integrate the shared-edge display meshes from `catallax-mesh`, then add STEP import, on-demand conversions, output hand-off checks, export accuracy levels and large-part stand-ins. Seam curves require a versioned exchange-contract change. Measure export interoperability rather than assuming a syntactically valid file works in every CAD tool.
7. Integrate the kernel into SolidDraft once the required operations are covered. Add field-only features, printing, FEA queries and interference checks with the hand-off guarantees in the PRD. FEA solver design remains separate.

## Background comparison prerequisites

The B-rep half has no display meshes yet. Build a deterministic display tessellator with explicit mesh-to-exact-surface deviation bounds and face origin sets. Add a B-rep recipe adapter for the same stock definitions and booleans, retaining operation names through trimming and merging. Neither is part of the field crate.

The application then needs revision-tagged edit regions, stable point-to-face associations and cancellation that cannot publish stale or partial results. Sample B-rep meshes against certified field queries, use bounds rather than raw residual magnitudes, and compare origin sets as well as locations and directions. The reverse check needs field boundary samples projected onto the appropriate exact B-rep faces, not only one-way vertex checks. Unresolved queries or unbounded upper distances must remain unresolved, never count as a 1 nm agreement certificate. Measure the whole two-way check on representative edits before claiming an interactive budget.


## Open policy conflict

The CLI currently uses a 50-second wall-clock deadline and returns `numbers.deadline` if it passes. Identical questions can therefore answer on an idle machine and refuse under load, contrary to the same-answer rule across machines and thread counts. Decide how to separate cancellation from deterministic geometric outcomes, then review any required exchange-contract change. No deadline behaviour changes in this publication.

## Measure useful latency

Record rebuild and query time per representative edit, including refusal latency and agreement-check cost. Speed per click matters. Optimize measured costs without weakening proofs. Safe background threads and cancellation must preserve deterministic answers and must never expose a half-built result to the application.

## Evidence for completion

Use mathematical cases, independent numeric oracles, scale and placement variations, repeatability checks and measured timings. Keep refusals, marked answers and wrong proofs distinct. An independent tester built without seeing the kernel's code uses the documented exchange contract. No unpublished tester design or acceptance system forms part of this plan.
