# Catallax

Catallax is an MIT-licensed geometry kernel in safe Rust for [SolidDraft CAD](https://github.com/BohnBawerick/soliddraft). One recipe describes a part. The boundary representation builds its faces and edges; the signed distance field evaluates the same definitions and builds field-only features such as lattices and gyroids.

This is alpha software. See [capabilities and release evidence](docs/CAPABILITIES.md) for what builds, imports and exports today. Valid cases can still refuse.

## Download or build

Get the command-line tools for Linux, macOS or Windows from [GitHub Releases](https://github.com/BohnBawerick/catallax/releases). See [platform requirements and download verification](docs/RELEASING.md#downloads). Catallax does not publish crates to crates.io.

From the repository root, with a Rust toolchain that supports edition 2024:

```sh
cargo build --release -p catallax-cli --bins
target/release/catallax about
```

On Windows, use `target/release/catallax.exe about`. For the application, see [SolidDraft releases](https://github.com/BohnBawerick/soliddraft/releases).

## Documentation

- [Roadmap](ROADMAP.md) and [glossary](docs/GLOSSARY.md) for the joint release plan and terms.
- [Architecture](docs/ARCHITECTURE.md) and [product scope](docs/PRD.md) for where changes belong and what a result must prove.
- [Building and contributing](docs/BUILDING.md), [recipe format](docs/RECIPE.md) and [verification exchange](docs/VERIFY-INTERFACE.md).
- [Code provenance](docs/PROVENANCE.md) and [third-party notices](docs/THIRD-PARTY.md).

[MIT licence](LICENSE). You can copy, change, sell or reuse the code with the licence notice intact.
