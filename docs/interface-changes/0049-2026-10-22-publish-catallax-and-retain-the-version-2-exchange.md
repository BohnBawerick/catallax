## 2026-10-22: publish Catallax and retain the version 2 exchange

Pieces: all, cli, exchange
Change: workspace packages, Rust crate paths and directories use `catallax-` and `catallax_`. `catallax-cli` builds `catallax` and `exact-brep` from the same main. Benchmark and dev-only deadline environment identifiers use `CATALLAX_`. Protocol strings `exact-brep-question/2`, `exact-brep-manifest/2` and `kernel.name = exact-brep` remain unchanged.
Why: the B-rep implementation becomes the starting code for a new public MIT kernel with one recipe and two evaluators. Rename packages without breaking the existing exchange.
Migration: Rust callers use the new crate names. Existing exchange clients may keep invoking `exact-brep`. A future protocol rename requires a versioned contract change. The new refusal rule, 1 nm necessity bound and operation-derived persistent names are documented requirements with explicit implementation gaps, not geometry changes in this publication.

