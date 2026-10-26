# Code provenance and licences

Catallax's own code is licensed under MIT, copyright 2026 Zigrawbohnbawerick Pty Ltd. This public tree starts from the existing B-rep implementation without importing its Git history. The publication changes package names and documentation, adds a binary name, and replaces the CI deployment; it does not import a reference implementation or change geometry algorithms.

## Reference implementation audit

The B-rep pieces were built from written contracts and mathematical tests, then compared with existing implementations. The source comparison records distinguish scratch adapters and reference checkouts from production code. Those records state that reference source did not enter the production workspace. Public reference projects considered in that work include Nacre, ogeom-rs and monstertruck. Algorithmic ideas and mathematical formulas are not a claim that their source code was copied.

No copied third-party file or function was found in this export. The publication audit checked the tracked source inventory, all source attribution and SPDX/copyright markers, package and lockfile dependencies, relevant source commit messages and each retained internal comparison record's source-use statement. No vendored reference source, submodule or reference-project production dependency was found. This is an audit of the available repository evidence, not proof of the absence of unattributed copying.

Because no imported implementation was identified, this publication adds no third-party source exception, NOTICE or per-file third-party SPDX header. If copied code is identified later, preserve its original licence and notices in NOTICE and its files, and correct this record before publishing it. The project's MIT licence does not override a dependency's terms.

## External dependencies

Cargo resolves registry dependencies under their own licences, as declared in their package metadata. Production uses the num arithmetic crates and serde JSON crates. The independent mathematical test oracles also use `rug` and `gmp-mpfr-sys`, with LGPL terms for those libraries and their native components. They are development dependencies, not a relicensing of the kernel.

Display meshing uses `spade` for constrained triangulation and its `robust`
predicate dependency, both under MIT OR Apache-2.0, and `libm` under MIT for
platform-independent angle projection. These are registry dependencies, not
vendored sources. No dependency type crosses the mesh library's public API.

Review the resolved dependency licences when distributing binaries or native test tools. This repository does not vendor those dependency sources or replace their licence files with the project's MIT text.

## Publication boundary

The public export omits operational infrastructure and internal comparison/replay records. It keeps code, mathematical tests and their useful technical descriptions. Research notes are published separately under [research/](research/), retaining source links and uncertainty labels. A cited paper is linked rather than republished as extracted full text.
