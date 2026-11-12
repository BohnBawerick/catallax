## 2026-11-13: direct solids with overlapping enclosing boxes

Pieces: 8, recipe import and explicit replacement

Previously, two disjoint solids whose enclosing boxes overlapped refused
`import.multisolid-separation`. The recipe importer now asks the existing exact
boundary-intersection stage for a complete, unmarked graph. Zero section edges
and zero isolated meeting points establish that the boundaries do not meet.
An exact point on every connected shell must then classify strictly outside
the other solid, in both directions. This second check excludes containment,
including a small solid wholly inside another without a boundary meeting.

Only a completed proof admits the pair. Contact, overlap, bounded decisions,
missing intersection coverage and exhausted work still refuse
`import.multisolid-separation`. The detail names the two source solids and gives
the section-edge, isolated-point and mark counts when a graph was completed,
or retains the underlying proof refusal. Zero boundary meetings with a shell
witness inside the other solid is reported separately. The existing strict
box proof still runs first and spends the same work as before.

No Boolean union or healing occurs. The importer appends the original lumps
in source order, preserving face names, membership and the recipe placement.
STEP and IGES recipe steps share the admission rule. Explicit STEP replacement
uses the same complete rebuild and keeps its existing reference and rollback
rules. Product assemblies, mapped occurrences and void-bearing STEP entities
retain their reader refusals. Unsupported boundary meetings remain coverage
refusals rather than clash verdicts.

Public synthetic checks cover a sphere in a torus hole, three parallel oblique
blocks with overlapping boxes, contact, overlap, duplicate material and nested
spheres in both source orders. The positive checks retain source membership,
checked topology, exact block volume and matching exported STEP, repair records
and deterministic work at 1, 2 and 4 workers. Raw evidence retained privately.
No saved replay answer, modelling tolerance or work allowance is changed.
