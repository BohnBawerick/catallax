## 2026-10-24: retain region ownership through lattice clipping

Pieces: field

Change: lattice clipping inherits the region snapshot's named-face ownership,
including predicates from earlier Boolean results.
Why: retaining the region's source-face registry without its ownership predicates
could restore a buried name on a surviving coincident cap.
Migration: no public signature, geometry, tolerance, budget or B-rep change.
The field-only lattice caller follows the same ownership rule as other trims.

