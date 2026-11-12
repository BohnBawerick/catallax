## 2026-11-12: STEP export refuses an undecided recognition instead of writing another layout

Pieces: 7
Change: band recognition, the band chart frame search and the cone axial-ring test refuse with the `numbers.*` code when a test they need is undecided. Before, an undecided test declined the special route, and the file was written in the legacy charted layout, on another chart frame or with legacy cone topology. The equator-split retry and the exporter's charted retry no longer answer a refusal that keeps an undecided decision or carries a work or precision code.
Why: the written layout depended on the work and precision allowance. With less work, a fillet that writes two band faces wrote the legacy layout instead (cone-top at 1,000 units a decision). See [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 2.4.
Migration: none for callers that already handle `write_step` refusals. A request that wrote the legacy layout only because a band test was undecided now refuses with the `numbers.*` code.
