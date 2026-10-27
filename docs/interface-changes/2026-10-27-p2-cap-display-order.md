# Preserve certified cap triangulation across the P2 chamfer edit

The public `parts/P2-collar.json` edit from distance 1 to 0.75 changes the
rounded lexicographic ordering of points on the upper annular cap. Previously,
its exact hull check rejected that floating proposal and sent the entire cap
through bridge/ear triangulation. At the existing 0.01 mm chord tolerance,
certified display then exhausted its 1,000,000-work allowance on face 8.

The planar mesher now repairs inversions in the proposed order with metered,
exact lexicographic comparisons. Coincident exact points still reject the
candidate. The exact hull, positive triangle orientation, oriented boundary
incidence, shared rails, returned-coordinate controls and fallback remain.
No point is merged, moved or omitted; no numeric/core cache changes, tolerance
changes or work-allowance increases are introduced.

Public release-build observations on main `7a347113d34131a54dc59af2c6f073da2af9099f`
and this repair, with a separate fresh allowance for each certified mesh:

| Controls | Before original / edit work | After original / edit work |
| --- | --- | --- |
| 0.01 mm chord, 1M allowance | 810916 / exhausted 1000000 | 812186 / 815854 |
| 0.00036 mm chord, 0.03490658503988659 rad angle, 20M allowance | 4179760 / 17260891 | 4186734 / 4203393 |
| 0.01 mm chord, 0.2617993877991494 rad angle, 1.5 mm edge, 20M allowance | 18191481 / exhausted 19999992 | 18192751 / 18196419 |

The dedicated P2 regression checks both complete recipes, closed oriented
incidence, reported chord/angular bounds, exact rational lengths of returned
binary64 triangle edges, shared boundary incidence and revision-local preview
names. A unit regression covers distinct exact coordinates that round to the
same binary64 value, and retains rejection of coincident exact points.

These are deterministic kernel-work observations, not SolidDraft click times.
Real application revolution, edit, picking, drawing and proof-gated output
acceptance remain a separate application slice. Isolated tessellation
qualification tier: SINGLEPAIR, to be routed through the coordinator; not
executed by this implementation lane.
