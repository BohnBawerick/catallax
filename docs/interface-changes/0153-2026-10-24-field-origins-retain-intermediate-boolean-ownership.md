## 2026-10-24: field origins retain intermediate Boolean ownership

Pieces: field, comparison, field display

Change: named field queries retain the boundary predicates of intermediate
Boolean results. Closed-contact sectors and cached display queries apply the
same predicates. Copies preserve the seed snapshot's ownership independently
of scalar geometry sharing.

Why: testing only the source leaf and final boundary could restore a buried
or clipped-away name when a later operand supplied a coincident carrier.
The B-rep correctly omitted these names.

Migration: no signature or format change. Consumers receive fewer spurious
coactive origins. Legitimate surviving paths and coincident aliases remain
selectable. B-rep geometry, STEP export, tolerances and work limits are unchanged.

