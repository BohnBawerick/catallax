## 2026-11-11: a refusal's Debug text ignores its decision

Pieces: all
Change: `Refusal`'s `Debug` prints only the code and detail, whether or not an undecided decision is attached. This corrects [2026-11-11-refusal-keeps-undecided-decision.md](2026-11-11-refusal-keeps-undecided-decision.md), which printed the decision when one was attached.
Why: answers embed refusal `Debug` text. The STEP import self-intersection check writes `format!("{r:?}")` into the `reason` of `import.self-intersection-not-verified`, and that refusal reaches recipe and field details through `Debug` again. With the decision printed, an undecided intersection inside that check would change an answer's text.
Migration: read the decision with `Refusal::undecided()`, not from `Debug` output.
