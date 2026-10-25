## 2026-10-24: stack-safe region sampling and bounded exact display samples

Pieces: 5, display mesh, application

Change: certified face-region sampling has a maximum subdivision depth of 32
and visits at most 4096 cells per chart build across all boundary pieces.
Exhaustion returns `boolean.sampling-budget`, which is not retried with greater
numeric precision. The inherited numeric precision limit still applies.
Display meshing retains at most 4096 exact traced planar boundary samples per
body and refuses with `mesh.resource-limit` before constructing another graph.
The existing vertex allowance and every requested display tolerance still apply.

Why: using numeric precision as the recursion allowance overflowed thread stacks
and made tiny nonzero coordinates much slower than exact zero. A vertex-only
display allowance did not constrain the expensive exact boundary expressions.
The new bounds govern work and retained data, not permissible coordinate sizes.
Linux application regressions isolate every case under a 512 MiB address-space
limit, an 8 MiB main-stack limit, 12 CPU seconds and 15 wall seconds.
These test limits are not a library-wide allocator or wall-clock policy.

Sub-tolerance inputs stay exact. DESIGN R9 forbids snapping merely because a
coordinate is below 1 nm. No whole-effect distance certificate or scoped
decide-once record authorizes replacing these coordinates with zero, so an
unsettled classification refuses instead. A resource refusal is unfinished
coverage, not a claim that the input geometry is impossible.

Migration: callers must preserve the two new typed resource refusal codes.
No public type, request format, numeric budget, modelling tolerance, STEP
representation or successfully sampled geometry changes. The guards do not
introduce timing, measurement or proof reporting into an application build.


