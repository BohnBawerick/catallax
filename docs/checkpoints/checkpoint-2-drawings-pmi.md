# Checkpoint 2 drawings and PMI

Adopted by [Checkpoint 2 decisions, section 4.8](checkpoint-2-decisions.md#48-drawings-and-pmi-the-p1-bracket-test), Q3, Q20 and Q26. Copied from the review report's bracket test and D1 to D6 plan. The bracket geometry is in [parts/P1-bracket.json](../../parts/P1-bracket.json).

## P1 bracket acceptance test

Part: P1 bracket (`parts/P1-bracket.json`: plate 60x40x5, flange 60x5x25, boss Ø16 h10 at (30,25) with Ø8 bore, four Ø4 holes at x10/50 y12/32, flange slot 20x6 with R3 ends at z18).
Manufacturing information added in SolidDraft (14 items): datums A (plate underside), B (flange back face), C (left end face); bore Ø8 H7, Ra 1.6; bore position Ø0.05 to A|B|C; boss Ø16 ±0.1, top 15 ±0.1 from A; 4x Ø4 +0.1/0 through; bolt pattern position Ø0.2 to A|B|C; slot 20x6 +0.1/0, R3; slot centre 18 ±0.1 from A; flatness 0.05 on A; perpendicularity 0.1 of B to A; ISO 2768-m; 6061-T6, break sharp edges.
Checks, all in the browser app:
- C1 round trip: STEP carries all 14 as semantic PMI on the right faces; NIST STEP File Analyzer no errors; OCCT reads all 14; re-import, save, reopen keep all 14 on the same faces. Runs first (the dot's advice).
- C2 edit: bore Ø8 to Ø10; dimension reads Ø10 H7, position tolerance still attached, C1 passes again.
- C3 one button: STEP plus drawing (PDF, DXF), front/top/side and a section through the bore; every item exactly once; holes, boss and slot located from datums; no overlaps; title block with name, material, general tolerance, revision; no hand moves.
- C4 nothing guessed: nothing on the drawing absent from the file; unplaceable items listed by name.
- C5 bare import: PMI stripped, recognition finds bore, boss, 4 holes, slot; the interview asks only the missing items with proposals; result passes C3 and C4.

## D1 to D6 implementation steps

The steps below preserve the review's technical scope and dependencies. Q25 and Q26 supersede its earlier split staffing and estimates: one Fable owner across kernel and app, moving to Claude Opus only if Fable quota runs short. The dot reviews. Work runs in parallel only where it does not touch 3.0A or 3.0B code; roadmap order is just behind 3.0B. No implementation starts before owner approval of the roadmap PR, Q29.

| # | Item | After |
|---|---|---|
| D1 | AP242 round trip: a hole diameter and tolerance on the part survives export, re-import, save and reopen (dot's advice) | now |
| D2 | Intent while modelling: threads, fits, general tolerance; drawing tolerances move to the part | SolidDraft 0.7.1 |
| D3 | Feature recognition on imports: half-day quiddity spike, then build or adopt | now |
| D4 | Dimension layout | D1, D3 started |
| D5 | The interview (rules, no AI) | D1, D2 |
| D6 | One button: STEP with PMI plus drawing | D4, D5 |
