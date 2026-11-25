Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-22: IGES sphere and torus export

Pieces: 7; verification exchange

Change: `write_iges` and the `catallax-iges/1` request accept whole spheres and ring tori that meet the exact writer's existing representation and tolerance rules. The files use analytic IGES 196 and 198 supports with shared circular seam edges.

Why: roadmap item 3.2 calls for exporting both solid types in IGES.

Migration: sphere and ring torus requests can now return an IGES artifact and its deviation report. Count the added vertices and seam edges from the report. Trimmed spherical and toroidal faces still return a typed refusal.
