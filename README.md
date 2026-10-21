# Catallax

Geometry should not require a vendor's permission. Catallax is a geometry kernel in safe Rust, free for everyone, for ever.

## Why this is free

Nearly every commercial CAD program runs on one of a handful of closed geometry kernels. Their owners licence them, meter them and can withdraw them. The maths inside is a century old and the algorithms are published. What is closed is the code, and the code is kept closed by law.

That law is not property. Property exists because things are scarce. If you copy this kernel, we still have it, and now you have one too. Nobody is poorer. Programs are not scarce, so the "intellectual property" over them has to be manufactured by the state: a grant that forbids you from arranging your own bits in a pattern someone else arranged first. Murray Rothbard wrote, "Patents, therefore, are grants of exclusive monopoly privilege by the State" in [Man, Economy, and State, chapter 10, section 7](https://mises.org/online-book/man-economy-and-state-power-and-market/10-monopoly-and-competition/7-patents-and-copyrights). Stephan Kinsella extended the point to copyright: "Ideas are not naturally scarce" in Against Intellectual Property, Journal of Libertarian Studies 15:2, 2001. We stand with that reading.

Geometry is mathematics and belongs to everyone. A public kernel lets a one-person CAD product enter the market without buying permission from a kernel vendor. Copy this code. Change it. Sell a product built on it. Compete with us.

Catallax puts AI coding agents to work on exact geometry, directed and checked by a human. Most of the inherited implementation was written that way. Whether this approach works at the scale a useful CAD kernel demands is an open question. We will find out in public.

## About the licence

Copyright attaches to code by default. We use the legal instrument available to let you use the code without asking us. MIT grants permission to use, copy, modify, publish and sell it, with the notice kept on copies and no warranty.

Licensed under MIT. See [LICENSE](LICENSE).

Catallax was started by Pedro De Paiva.

## One recipe, two kernels

One recipe, the feature tree, is the master. An exact B-rep kernel and a signed distance field kernel build each part from it, using one mathematical definition for each operation. Faces and edges will be named by the operation that created them, not their position in a result list.

The kernels must agree exactly where mathematics allows, and within 1 nm only where necessity forces a tolerance. Impossible geometry is the only acceptable final refusal. Every other refusal is unfinished work. The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels, not a claim to beat them.

[The PRD](docs/PRD.md) sets that direction. The B-rep starting code remains intact. [RECIPE.md](docs/RECIPE.md) defines the shared language, and `catallax-field` now evaluates its stock solids and booleans with exact predicates, bounded distances and named face queries. [FIELD.md](docs/FIELD.md) states the query guarantees and unresolved cases. Cross-kernel agreement checks are not built. The inherited B-rep mark APIs do not yet enforce the new 1 nm geometric limit, and its 50-second deadline conflicts with the same-answer rule. These gaps remain explicit in [the plan](docs/PLAN.md).

## What works today

The code is early, and supported operations can still refuse. Measure it by its results, not the length of its roadmap.

- `catallax-numbers` supplies exact arithmetic, certified enclosures and decisions, polynomial roots, construction facts and marks.
- `catallax-shapes` stores exact rational poses, analytic geometry and derived points by definition, with proven parameter enclosures.
- `catallax-topology` builds blocks, cylinders, spheres, cones and tori in append-only curved topology and checks structural faults.
- `catallax-intersect` handles plane, cylinder and sphere pairs, curve-surface meetings, branch-aware queries and certified traces. General cone and torus boolean coverage is missing.
- `catallax-boolean` implements unite, subtract and intersect through intersect, imprint, classify, keep and sew. Ambiguous and unsupported cases remain unfinished work.
- `catallax-measure` computes volume, area and centroid with proven bounds for the body actually built.
- `catallax-step` exports deterministic STEP AP242 with analytic geometry, certified spline fits for intersection branches and a deviation report. The legacy contract has interoperability limits, including no seam curves.
- `catallax-mesh` builds shared-edge display meshes of analytic faces and certified intersection traces, with picking IDs, proven chord bounds and OBJ export.
- `catallax-recipe` reads, checks and writes recipes in the canonical `catallax-recipe/2` JSON form. It evaluates no geometry.
- `catallax-field` compiles recipe stock solids and booleans into a fixed-order tape, keeps names outside that tape, and supplies tentative and proven scalar/batch queries. Boolean residuals carry conservative distance bounds, not a true-distance claim.
- `catallax-cli` builds `catallax` and the legacy `exact-brep` program from the same main. They read a question, build and operate, measure, and write a manifest and optional STEP file.

Sketch features, fillets, chamfers, shell, offset, draft, sweep, loft and STEP import are still missing. The kernel does not yet implement the whole PRD. Existing ignored regressions remain visible; publishing the code does not resolve them.

## Building

```sh
cargo build --release -p catallax-cli --bins
target/release/catallax about
target/release/catallax check question.json
target/release/catallax answer question.json --out answer
```

`check` validates without evaluating geometry. `answer` writes `answer/manifest.json` and, for a nonempty solid when requested, `answer/answer.step`. Refused and invalid questions also write a manifest and exit 0. [JUDGE-INTERFACE.md](docs/JUDGE-INTERFACE.md) defines the exact command and file contract, including its unchanged legacy names.

Local checks:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
python3 -B -m unittest discover -s tools/edge-hunt -v
```

GitHub-hosted CI runs formatting, clippy, all enabled release tests and doctests. Crate `TESTS.md` files describe extra opt-in workloads. Timing benchmarks report measurements, never correctness gates.

## Display meshes

`catallax_mesh::mesh_body(&model, body, &options, &mut decider)` returns an
`Answered<Mesh>`. `MeshOptions::chord_tolerance` is a positive exact length in
millimetres. `max_vertices` limits deterministic refinement, including temporary
periodic chart copies. The default is 0.01 mm and one million vertices.

`Mesh::positions` is shared across faces. Each `Triangle` carries its input
`FaceId`; each `FaceMesh` identifies a contiguous triangle range and its proven
maximum chord error. `EdgeMesh` records the one polyline used by both incident
faces, including its proven curve-to-polyline error. Face IDs belong to the
input model, not to rebuild-stable recipe names.

Planes, cylinders, cones, spheres and tori use their exact definitions.
Intersection edges use piece 4's certified traces. Surface bounds use Taylor's
theorem and analytic second-derivative bounds over parameter triangles.
Bounds include replacement of boundary vertices and binary64 coordinate
rounding. They describe chord error against the analytic face and edge
definitions, not a proof that polygonal trims equal the exact curved region.
All covered definitions have proven bounds, so no sampled estimate is presented
as a proof. Input and trace marks are preserved.

The implementation is serial with fixed traversal, subdivision and output
orders. Pure-Rust angle projection and dyadic periodic coordinates avoid
platform trigonometric libraries and inconsistent seam copies. Concurrent
callers do not share state. Tests compare OBJ bytes with 1, 2 and 4 workers,
count every mesh edge twice with opposite directions, and sample triangle
interiors against independent analytic distance formulas.

Write an OBJ for Blender without changing the judge CLI:

```sh
cargo run --release -p catallax-mesh --example obj -- torus torus.obj 0.05
cargo run --release -p catallax-mesh --example obj -- cross intersection.obj 0.05
```

The example also accepts `block`, `cylinder`, `cone`, `frustum`, `sphere`,
`drill` and `lens`. Library callers write any built body with `Mesh::write_obj`.
OBJ groups are named `face_<id>`, and their triangles share vertex indices.
Resource or coordinate-precision failures return `mesh.precision-or-budget`;
invalid controls return `mesh.invalid-options`. No over-tolerance or
non-watertight mesh is returned as a successful answer. A caller-supplied
wall-clock deadline can still cancel meshing through `Decider`; the default
library decider has no deadline.


## Where to start

- [PRD](docs/PRD.md), [plan](docs/PLAN.md) and [research](docs/research/) explain the direction.
- [Recipe language](docs/RECIPE.md) defines operations, names and the recipe file.
- [Field queries](docs/FIELD.md) defines distance semantics, proof limits, selection and batch storage.
- [Design](docs/DESIGN.md) sets the rules; [architecture](docs/ARCHITECTURE.md) connects the B-rep crates.
- [Piece contracts](docs/pieces/) describe the current interfaces and mathematical tests.
- [Interface changes](docs/INTERFACE-CHANGES.md) records migrations.
- [Provenance](docs/PROVENANCE.md) records the publication audit and dependency licence boundary.

## Contributing

Read the design rules first. A guessed topological decision, an unproven construction fact, an `unsafe` block or a silent fallback is a defect. Keep API and format changes explicit, update callers and document the migration. Bring mathematical evidence for geometry changes and report refusals honestly.

Agent-written code is welcome. Say how it was produced and stay accountable for what it does. This repository's independent tester was built without seeing the kernel's code.

## Other products

[check out my other products](https://zigrawbohnbawerick.org/)

## Sponsorship

Sponsors can help pay for work on Catallax. If you want to sponsor the project, open a GitHub issue to get in touch. There is no sponsor page yet.
