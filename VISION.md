# Vision

`Catallax` exists so that anyone can design, edit, analyse and deliver a mechanical part without a kernel vendor's permission.
It serves SolidDraft first, and anyone who copies the code to build their own product.
It turns a recipe, the feature tree, into a part the user can edit again, measure, analyse and send to manufacturing as STEP.
It owns exactly one thing: the part's geometry and the answers computed from it.

## Finish the part

The test of a product change is whether a real part gets finished, edited again, analysed or delivered.
A refusal of possible geometry is unfinished work, never a product boundary.
Only impossible geometry, such as a non-positive radius or a self-touching shell, is a legitimate final refusal.
Faster refusals, more refusal codes and better labels do not count as coverage.
Coverage is measured on whole edit chains on real target parts, against a reference kernel doing the same job on the same part.
A result built by another route, such as the field, is welcome when the part it yields can be edited again and delivered, with its deviation reported.

## Research is how the kernel is built

Catallax is a research kernel: experiments and prototypes are how it finds the best way to build each capability.
A prototype does not need a blocked part to justify it.
Experimental code may land in the shipping crates, switched off by default, so SolidDraft can try it on real parts.
Each experiment records what it measured, and is adopted, parked or dropped on that evidence.
A prototype stays labelled as one until it is adopted, and never appears as landed capability.
Compatibility with earlier versions is never a requirement: old recipes rebuild with the better algorithm, and no old algorithm is kept alive for them.

## Numerical where needed, honest everywhere

Catallax builds editable numerical B-reps under stated tolerances where exact construction is unavailable or not worth its cost.
Tolerance values, work allowances and caps are engineering choices that change on evidence; the commitment is that each one is stated.
Every numerical result carries its tolerance, the checks it passed and any repair it received.
Each answer claims what its checks establish and nothing more.
Proven, bounded and refused answers stay separate, and nothing turns a bounded answer into a proven one.
A proven answer is never wrong, and a wrong proof stops work on that method until it is repaired.
Existing exact proofs keep their named scope, and a numerical edit does not wait for a whole-part theorem.
When an exact STEP file cannot be written, Catallax writes an approximate one that opens and edits in another CAD tool, and marks it approximate in the export report.
Catallax never silently returns the input, fills a gap, moves geometry or substitutes a mesh.

## Efficiency, not machinery

Hard cases get better algorithms and shared methods, not more infrastructure around slow ones.
A larger work allowance is an explicit choice by the caller, and the operation restarts honestly under it.
Catallax does not build durable proof jobs, proof services, independent certificate systems or research into building them.
A new subsystem earns its place only when it helps finish, edit, analyse or deliver a part, or establishes a claim the user asked for.
Work is counted in deterministic units, so the machine never changes the answer: one version gives the same outcome, topology and bounds on every machine and at every thread count.
Speed per click comes first from algorithms that need less work, then from cores and the GPU running that work in fewer seconds.

## One recipe, geometry by definition

The recipe is the master, and every operation has one written mathematical definition.
Surfaces, edges and vertices are stored by definition; coordinates and fitted curves are caches.
Faces and edges are named by the operation that made them, never by their order in a result.
The signed distance field evaluates the same recipe, checks the B-rep on demand, owns field-only features such as lattices, and carries analysis such as stress without a volume mesh.
An unchecked field answer is labelled unproved and never stands in for an export or a measurement.
Catallax borrows methods from other kernels, such as OpenCASCADE and SMLib, but never links or calls another kernel as a fallback.
Ported code proposes candidates, and Catallax's own checks decide what is published.

## Free to copy

The code is MIT licensed safe Rust, with no licence fees, seats or royalties.
Copying, changing, selling and competing with Catallax are welcome.
Adapted third-party code keeps its notice and attribution, and published ideas are credited without copying their source.

## Scope

Catallax is not a CAD application; SolidDraft owns the interface, documents and user workflow.
It is not a proof service, a job scheduler or a certification authority.
Public text describes landed work and its limits; plans and experiments stay in the roadmap and research notes.

A change aligns when it finishes, edits, analyses or delivers a real part that failed before, makes an existing answer cheaper without changing it, or runs a measured experiment toward one of those, and its claims match its checks.
A change should be resisted when it adds machinery, labels or refusals instead of geometry, claims more than its checks establish, makes an answer depend on the machine or the clock, keeps old behaviour alive for compatibility, calls another kernel as a fallback, or presents a prototype as landed capability.
