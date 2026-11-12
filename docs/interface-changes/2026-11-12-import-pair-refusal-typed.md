## 2026-11-12: an undecided imported face pair keeps its decision

Pieces: import, 5
Change: `ImportRefusal::SelfIntersectionNotVerified` gains `refusal: Option<Box<Refusal>>`, the engine's refusal for the pair: a `numbers.*` refusal that keeps the undecided decision and its bound, or `import.pair-undecided` (outside what the curved pair check covers) or `import.pair-setup` (the faces could not be set up). It is skipped on the wire. `reason` is now that refusal's `code: detail` text; a refusal from the engine used to be written in Rust debug form. `import.self-intersection-not-verified` refusals from recipe evaluation now keep the undecided decision.
Why: the curved pair check kept only a string, so a work or precision limit could not be told apart from unsupported geometry, and a later tolerance rule had no bound to read. See [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 2.4.
Migration: code that builds or fully destructures the variant adds `refusal` (or `..`). The JSON form is unchanged apart from the `reason` text.
