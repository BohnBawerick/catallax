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

The listed validation completed with ten CLI print tests and four moved
encoding/topology unit tests passing, formatting passing, and all-target Clippy
for both affected crates passing with warnings denied. Cargo retains its existing
manifest warning about two binary names sharing `src/main.rs`; no lint waiver was
added.

### Field-envelope admission correction

Independent review found that the initial shared API at
`ec89725e870372f43eb11f742bac0c1365cb7653` omitted the CLI's field-definition
guard. A public 4 x 4 x 4 box followed by a cubic lattice (cell size 2, strut
radius 0.25) reproduces the defect: `evaluate` retains the region as a stand-in,
and the original shared API wrongly returned a 684-byte solid STL. The new
`retained_body_refuses_result_reachable_lattice_before_stl_encoding` regression
failed with that result before the correction. Earlier passing API/CLI checks
did not cover this admission gap and did not establish readiness.

The existing result-reachable `native_field` predicate now resides in the shared
printing module. Both CLI route selection and the retained-body guard use it.
The guard checks `Built::recipe()` and retains `print.brep-field-definition`;
it never substitutes the solid envelope or retries the field route. A separate
regression retains an unused lattice stand-in while selecting the ordinary box
result and verifies successful API/CLI byte agreement. This deliberately does
not reject every nonempty `stand_ins` list. No claim is made that certified
construction supports lattice definitions.

| Item | Status |
| --- | --- |
| Reusable current-body print proof and STL encoding | Implemented; validation recorded with the candidate head |
| Result-reachable field envelope admitted as exact B-rep | Reproduced on initial API head; shared reachability guard and regressions added |
| SolidDraft panel, revision/settings publication gate and certified 3MF | Pending application integration |
| U1/U3 rendered export acceptance and U6 over-cap request | Not exercised by this kernel slice |
| Periodic STEP finite-trim failures, including P3 | Unchanged; separate writer ownership |
| PR35 CLI over-tolerance guard | Separate SolidDraft branch; not superseded here |

No repaired STEP-output policy, dependency revision, budget, tolerance or test
expectation changes are included. Qualification routing remains with Firstmate.
