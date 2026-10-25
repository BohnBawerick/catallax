## 2026-10-18: ND8 questions carry no labels

Pieces: judge interface, cli
Change: the question `id` and all names are opaque; questions carry no pose, family or cell label and the kernel reads none. The examples now use opaque ids. This confirms the judge's provisional behaviour.
Why: an answer must come from the geometry alone.
Migration: none.

