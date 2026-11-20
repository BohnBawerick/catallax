## 2026-11-20: feature recognition on a body

Pieces: 6; SolidDraft adapter
Change: `catallax_measure::recognition::recognise(model, body, dec)` returns `Recognised`: the body's round holes (through or blind, with their faces, diameter and depth), groups of equal holes, round bosses and round-ended slots. It reads planes and circular cylinders and decides every condition exactly. It adds no refusal code; a decision that cannot be made returns the decider's refusal.
Why: the Checkpoint 2 drawings and PMI item ([decisions, section 4.8](../checkpoints/checkpoint-2-decisions.md#48-drawings-and-pmi-the-p1-bracket-test), step D3) needs the features of an imported part with no manufacturing information, so an interview can ask only for what is missing.
Migration: none. The function is new and reads the body as it stands.
