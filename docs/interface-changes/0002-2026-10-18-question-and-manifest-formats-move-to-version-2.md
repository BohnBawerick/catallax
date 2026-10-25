## 2026-10-18: question and manifest formats move to version 2

Pieces: judge interface, cli
Change: the question format is now `exact-brep-question/2` and the manifest format `exact-brep-manifest/2`, carrying the ND1 to ND11 decisions below. `catallax_cli::QUESTION_FORMAT` and `catallax_cli::MANIFEST_FORMAT` change to match. Version 1 is withdrawn.
Why: the decisions change the default deadline, when a STEP file is written and how `result` counts touching lumps, and `JUDGE-INTERFACE.md` requires a version bump for any format change.
Migration: the `catallax-cli` creator writes and reads only version 2; `about` lists only `exact-brep-question/2`.

