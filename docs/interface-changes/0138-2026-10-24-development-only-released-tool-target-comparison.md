## 2026-10-24: development-only released-tool target comparison

Pieces: verification tooling and target scoreboard

Change: the driver adds `head-to-head-reference` and `head-to-head`. Format `catallax-head-to-head/1` records pinned alternating P0–P5 construction/edit pairs, separate display and requested measurement, work, operations rebuilt, refusals, independent property checks, name survival and OpenCASCADE STEP import checks. SolidCraft v0.21.0 is pinned by archive and executable SHA-256; scripts are authored from the target definitions.

Why: compare complete edit workloads and incorrect results without adding timing, measurement or reporting to application builds.

Migration: existing kernel, recipe, application and verification exchange interfaces are unchanged. Build the opt-in `catallax-cli` example `head_to_head` and follow the benchmark rerun guide. Missing measurements remain unmeasured; sampled checks and numeric independent-oracle comparisons are not certificates.

