## 2026-10-26: require nonempty complete baseline evidence

Pieces: verification driver and comparator.
Change: correct the [explicit baseline reuse interface](2026-10-26-explicit-baseline-evidence-reuse.md)
to reject empty interpreted qualification corpora and validate every required
operation-log field, type and conditional null/certificate form before reuse.
Why: public reproductions confirmed that coherently resealed incomplete logs
and zero-case adapter runs could pass the earlier validation.
Migration: no producer protocol changes. Current complete
`catallax-operation-log/1` records, including interrupted and refused forms,
remain valid; incomplete older layouts cannot seed reuse. Preserve the existing
face-check internal-display schema separately from delivered field/print OBJ
records. Ordinary adapter loading and SINGLEPAIR remain unchanged. Requalify
after diagnosing rejected evidence; do not retry until green or trust stored
equality flags.
