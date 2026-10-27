# Retained-body print API

`catallax_field::printing::brep::mesh(&Built, &Options, &mut Decider)`
qualifies an already-built exact body. `mesh_certified` accepts a `CertifiedBuilt`
and verifies its geometry-bound certificate. Neither call evaluates a recipe,
uses a preview, or requests a measurement. The CLI retains recipe construction
and route selection, then calls this same adapter.

The implementation moves the existing binary32 motion, topology, shell sense,
nesting and finite-trim bound checks from the CLI; it does not add another mesh
or proof implementation. Existing options/defaults remain the single controls.
The native-field wrapper also retains its selected options and never retries
a refused B-rep as a field request.

The returned `Output` has private fields and read-only geometry/bound accessors.
It retains a copy of the proof's options. `Output::stl()` creates the existing
metadata-bearing STL entirely in memory and returns the conservatively rounded
certificate alongside it. It takes no replacement options. The CLI's existing
streamed publication uses the same header encoder. Failed proof/encoding returns
no bytes for publication. Application publication is a separate responsibility:
revalidate current body identity, revision, proof status and selected settings
before committing any destination.

## Public evidence and open defects

Base: Catallax `5c6d6b25d55c9fdb52295a832616a7b592063787`.
Run from the repository root:

```sh
cargo test -p catallax-cli --test printing --locked
cargo test -p catallax-field --lib printing::brep --locked
cargo clippy -p catallax-field -p catallax-cli --all-targets --locked -- -D warnings
cargo fmt --all --check
```

The retained-body regression compares API and CLI STL bytes for a planar box,
decimal-coordinate box and cavity. It also mutates the caller's settings after
proof, verifies the result retains the original settings, and checks an
insufficient vertex cap refuses. Existing print regressions retain topology,
cavity orientation, rounding, work/vertex caps and no-file refusal assertions.
These are public functional checks, not application click timings or corpus
qualification.

The listed validation completed with eight CLI print tests and four moved
encoding/topology unit tests passing, formatting passing, and all-target Clippy
for both affected crates passing with warnings denied. Cargo retains its existing
manifest warning about two binary names sharing `src/main.rs`; no lint waiver was
added.

| Item | Status |
| --- | --- |
| Reusable current-body print proof and STL encoding | Implemented; validation recorded with the candidate head |
| SolidDraft panel, revision/settings publication gate and certified 3MF | Pending application integration |
| U1/U3 rendered export acceptance and U6 over-cap request | Not exercised by this kernel slice |
| Periodic STEP finite-trim failures, including P3 | Unchanged; separate writer ownership |
| PR35 CLI over-tolerance guard | Separate SolidDraft branch; not superseded here |

No repaired STEP-output policy, dependency revision, budget, tolerance or test
expectation changes are included. Qualification routing remains with Firstmate.
