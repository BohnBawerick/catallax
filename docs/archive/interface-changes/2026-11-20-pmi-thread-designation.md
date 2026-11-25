Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-20: a thread designation on a diameter

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::Size` gains `thread: Option<String>`, the designation of a threaded feature such as `M6x1 - 6H`. `write_step_with_pmi` writes it as a user defined attribute `thread` (`pmi::THREAD`) on the feature the diameter is on, and `read_step_pmi` reads it back. A thread on a radius, or a designation that is not plain text, refuses with `step.pmi-invalid`.
Why: the Checkpoint 2 drawings and PMI item (step D2) puts threads on the part as intent while modelling, and a STEP file has to carry them. AP242 edition 1 has no thread entity, so the form is the one already used for a through hole.
Migration: a caller that builds `Size` adds `thread: None`. Files written before read with `thread: None`.
