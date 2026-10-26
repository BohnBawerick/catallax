## 2026-10-19: a step may name one body twice

Pieces: judge interface, cli
Change: a step's `args` may name the same operand or step twice. `X ∪ X` and `X ∩ X` are `X`, and `X − X` is empty. `catallax-cli` answers such a step by that set rule, without calling piece 5's `boolean`, whenever both arguments resolve to the same `BodyId`. The result carries `X`'s marks once, and the step makes no decision. The question format is unchanged, so its version stays at 2.
Why: the question format already accepted repeated names, but piece 5 treats one `BodyId` passed as both operands as a programming error, so such a valid question made the program panic.
Migration: the CLI builder applies the set rule. Pieces 1 to 7 are unchanged.

