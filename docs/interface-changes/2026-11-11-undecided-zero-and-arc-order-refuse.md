## 2026-11-11: equal-crossing saddle tests and imported arc order refuse when undecided

Pieces: 2, import
Change: `RollingSaddleDefinition::crossing`, `RollingSaddle::is_edge` and `RollingSaddle::centres_on` refuse when an exact-zero test they need is undecided (a `NotExact` sign, or a zero proven only by a separation bound). Before, they read it as "not equal". A part proven nonzero still answers no. In import, an undecided or refused order along an open edge makes a curved face pair not verified instead of "not on the arc".
Why: three of the ten unsafe paths in [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 2.2. "Not equal" picked another saddle family, told measurement a band's centres were elsewhere, and could drop a real cut from an imported face pair or invent a boundary contact.
Migration: callers already handle these functions' refusals. A refusal now carries the `numbers.*` code, or `numbers.geometric-gap-unproven` for a separation-bound zero, instead of an answer derived from it.
