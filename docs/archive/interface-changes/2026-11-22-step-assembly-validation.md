Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-22: STEP assembly exchange validation

Pieces: 7, 8

Change: Assembly export includes a file-level `complete_body_bound=unestablished` warning when a part report lacks a complete deviation bound. It refuses void shells with `step.assembly-unsupported-part`. Assembly import refuses unsupported representation items, including mapped leaves, and requires solid colours to resolve equally on both sides.

Why: A mapped leaf can carry its own motion, and a one-sided style cannot be recovered as a uniform part colour. The file must retain the part writer's bound warning and the assembly writer must respect the reader's void-shell limitation.

Migration: Use direct manifold B-rep items for assembly leaves. Provide a single two-sided colour or matching colours on both sides of each solid. The ordinary part writers and readers keep their contracts.
