## 2026-10-25: share retained conic contact certificates

Pieces: 2, 5
Change: `GeomStore::has_conic_incidence(curve, point) -> bool` exposes the
read-only retained plane/quadric certificate for a circle or ellipse. Probe
paths use it to preserve source contact before expanding barrier equations.
Why: the same exact certificate must survive from carrier membership into
contact coefficients and straight-path restrictions.
Migration: `false` means no certificate, never geometric nonmembership.
Consumers using the proof as a decision call `Decider::note_fact_decision`;
ordinary membership does this itself. Whole-path exclusion and intersection
component checks remain required. Allowances, tolerances and codes are unchanged.
See [the supporting conic theorem](2026-10-25-retained-planar-conic-incidence.md).
