## 2026-11-02: one-solid USD BRep export

Pieces: features, exchange verification

Change: The opt-in `catallax-features` example `usd_brep` writes `analytic-square-through-hole/1` as a USD `BrepArray` with a companion display mesh, recipe and `catallax-usd-golden/1` metadata. It supports the 2 mm baseline, a 3 mm identity-preserving height edit and one rigid placement. `tools/checks/usd_brep.py` decodes the actual USD through OpenUSD and checks the analytic fixture and eight negative controls.

Why: Establish one-way analytic exchange, incidence and identity evidence before adding more geometry or an application-facing command. NVIDIA consumer interoperability remains unproved.

Migration: No existing caller changes. The example requires a new output directory, height, placement and a caller-supplied full source commit. Retain `part.usda`, `metadata.json` and `recipe.json` together. Do not treat the display mesh, interchange tolerance or bare USD file as a lossless recipe or general certificate transport. Run the commands under [Optional USD BRep and Omniverse exchange](../archive/PLAN-2026-11.md#optional-usd-brep-and-omniverse-exchange).
