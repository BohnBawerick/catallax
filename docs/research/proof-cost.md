# Proof-cost research

Outcome: recorded. Its findings stand. The machinery it led to was rejected.
When: 13 to 14 November 2026.

## What we were trying to achieve

An honest answer to four questions:

1. What does it cost to prove geometry, and is 20M a real boundary between parts that can and cannot be built?
2. How do SMLib and OpenCASCADE handle the same cases?
3. Is "sample, fit, then zoom in" a sound method?
4. Can a GPU rescue a case that runs out of work?

## Why it was proposed

On 13 November the owner ruled that kernel work must be general, and that when a general method cannot meet a strict requirement the right move is an exception, not more tuning: drop or relax the proof, accept a measured threshold, or change strategy. He named SMLib and OpenCASCADE as the yardstick: if they solve a case and Catallax cannot, Catallax must catch up or its requirement is wrong. He asked for research on the cost of proof and on the exception options so he could decide.

## What was tried

Reading, not benchmarking: seven earlier reports, a week of merged changes, the source of SMLib and OpenCASCADE at pinned versions, and the primary papers on certified geometry. The only new computation was an exact counterexample to sampling. No kernel was run.

## Result

- No fixed allowance covers every input. Spline inputs have no limit on size, degree, coefficient precision or how badly conditioned they are.
- Nothing found says ordinary CAD is too dear to certify. Several papers show large savings from better algorithms. None shows every part fits one constant.
- Sampling alone is unsound. "Sample, fit, bound the unsampled region, then refine" is an old and sound family. "Refine only where samples look bad" is a guess: a tiny loop or a narrow spike can leave every sample looking good. The report builds an exact example.
- A GPU does not rescue a work refusal. Work units count arithmetic, so faster hardware does not lower the count. A GPU helps only after a bottleneck that can be batched has been found.
- The reference kernels do not prove these things. They work to tolerances and checks, and a valid file from them is not a proof of the original design.

## What was decided from it

The report recommended, among other things, separately requested proof jobs that could be resumed. The owner asked for a specification of that and of a bounded modelling contract. On 14 November the owner rejected the specification as over-engineering: the point of a kernel is to create the geometry.

The owner then chose the minimum option in each of four areas:

- Extra proof work: better algorithms. No proof-job system.
- Numerical modelling: editable numerical B-reps under stated tolerances, with meaningful checks.
- Certification: report what the checks establish. No general certificate system.
- Release and use: a numerical STEP file that opens and edits elsewhere, marked approximate.

Those four choices are in `VISION.md` and are milestone 1.1 on the roadmap. The research findings above are why the roadmap treats refusals as work for better methods, and why the speed phase ranks cheaper proofs first.

## Sources

- `VISION.md`, sections "Numerical where needed, honest everywhere" and "Efficiency, not machinery"
- The research report "What proof should cost in Catallax", 13 November 2026, the rejected specification, and the options report "Twelve product choices for Catallax", 14 November 2026. All three are kept outside this repository.
- The bounded answer kind that already existed: https://github.com/BohnBawerick/catallax/pull/576
- Related: [cheaper-proof-p1.md](cheaper-proof-p1.md), [body-checker-audit.md](body-checker-audit.md)
