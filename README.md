# Catallax

Geometry should not require a vendor's permission. Catallax is a free geometry kernel in safe Rust, built for [SolidDraft CAD](https://github.com/BohnBawerick/soliddraft).

We support free markets and oppose intellectual property laws. Copying code does not take it away from anyone. A public kernel lets a solo engineer or a new CAD product enter the market without buying permission from a kernel vendor. Copy this code. Change it. Sell a product built on it. Compete with us.

Catallax is developed by Zigrawbohnbawerick Pty Ltd. No licence fees, seats or royalties. [Why this is free](https://kernel.zigrawbohnbawerick.org/why.html).

[Website](https://kernel.zigrawbohnbawerick.org/) · [Code and downloads](https://github.com/BohnBawerick/catallax/releases) · [Other products](https://zigrawbohnbawerick.org/)

## One recipe, two representations

An exact boundary representation (B-rep) and a signed distance field evaluate a shared subset of one feature tree. Geometry is stored by definition. Faces and edges keep names tied to the operations that created them. The field can check the B-rep through a separate, sampled comparison. This is not a whole-boundary proof.

In design: general tolerance-policy composition, broader agreement checks between the two representations, and complete industrial CAD coverage. The [roadmap](ROADMAP.md) separates landed work from those goals.

## What works today

Catallax 0.5.0 includes:

- Exact blocks, cylinders, spheres, cones and tori; unite, subtract and intersect in supported geometric configurations.
- Sketch extrusions and revolutions, holes, patterns, mirrors, sweeps, lofts, fillets, chamfers, offsets, shells, thicken and draft within their documented domains.
- Proven bounds for supported measurements, display meshes and STEP AP242 exports. Analytic STEP import, planar repair and spline source previews have separate admission and proof limits.
- A field evaluator, sampled two-kernel comparisons, cubic lattices and network gyroids with labelled print output.
- An application entry used by SolidDraft, revision-bound construction display and deferred proof, plus a deterministic verification CLI with replay and worker-count checks.
- Tolerant geometry beside exact geometry: an edge gap is admitted at a measured local tolerance or refused, and approximate results are labelled.
- STEP assemblies with shared parts, placements, names, colours and supported PMI; clash, clearance and mass checks; and drawing curves for sections, silhouettes and hidden lines.
- First builds of bend, twist, curvature-continuous fillets, mesh Booleans and mesh decimation.

This is alpha software. Valid cases still refuse or take seconds. Checkpoints 3 and 4 measured the 0.5.0 release against OpenCASCADE and SMLib:

- Edit chains: Catallax finished 11 of 13 test parts and 9 of 10 census parts. OpenCASCADE finished 11 and 10, SMLib 7 and 8.
- STEP census: 767 of 3,000 public files (25.57%) export as valid solids. OpenCASCADE's healed reference is 2,701 (90.03%).
- Coverage: 31 of 61 reference-kernel lines (50.8%), with working evidence on every line and the full bar met on one.

Catallax still refuses the general cone and torus intersections, all three phase-2 targets (the threaded boss, the lofted duct and the curved housing) and one census-part drill. Some answers cost more than the 20M work target.

See [Checkpoint 3](docs/checkpoints/checkpoint-3.md), [Checkpoint 4](docs/checkpoints/checkpoint-4.md), the [coverage recount](docs/checkpoints/checkpoint-3-coverage.md), the [0.5.0 release scope](docs/milestones/M5.md#050-phase-4-alpha-release) and [implementation limits](docs/CAPABILITIES.md).

## Download or build

Download an alpha CLI archive for Linux, macOS or Windows from [GitHub Releases](https://github.com/BohnBawerick/catallax/releases). Archives include `catallax`, `catallax-verify` and the `exact-brep` compatibility command. See [platform requirements and verification](docs/RELEASING.md#downloads). The binaries are not platform code-signed or notarized. Crates are not published to crates.io.

To build from a checkout, use a Rust toolchain supporting edition 2024 and run from the repository root:

```sh
cargo build --release -p catallax-cli --bins
target/release/catallax about
```

On Windows, run `target/release/catallax.exe about`. See [building and verification](docs/BUILDING.md) for question inputs, local checks and replay commands. For the CAD application, use [SolidDraft's releases](https://github.com/BohnBawerick/soliddraft/releases).

## Documentation

- [Recipe language](docs/RECIPE.md).
- [Field queries](docs/FIELD.md) and [sampled comparison](docs/COMPARE.md).
- [Verification CLI](docs/VERIFY-INTERFACE.md) and [legacy question interface](docs/JUDGE-INTERFACE.md).
- [Design rules](docs/DESIGN.md), [architecture](docs/ARCHITECTURE.md) and [interface migrations](docs/INTERFACE-CHANGES.md).
- [Contributing](docs/BUILDING.md#contributing), [provenance](docs/PROVENANCE.md) and [third-party notices](docs/THIRD-PARTY.md).

## Licence and support

[MIT](LICENSE). Use, copy, modify, publish and sell the code with the notice kept on copies and no warranty.

Support development through [GitHub Sponsors](https://github.com/sponsors/BohnBawerick) or [Ko-fi](https://ko-fi.com/pedropaiva56920). The scope and limits above describe the current code.

## Thanks

Catallax draws on openly published research, standards and engineering accounts, including Fornjot, Fidget, libfive and OpenCASCADE. Their published ideas informed the work; their source code was not copied into the kernel. No affiliation or endorsement is implied. See [special thanks](https://kernel.zigrawbohnbawerick.org/thanks.html) and the [research notes](docs/research/).
