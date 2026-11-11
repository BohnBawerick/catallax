## 2026-11-11: exact spline section imprint queries

Pieces: 5.

Change: Boolean curve meetings consume certified `PatchSection` definitions through complete cutting-surface queries and exact component membership. Finite rational patch and atlas body queries check incidence before ray parity, without requesting a global implicit side. The new refusal `boolean.patch-section-meetings` names a missing finite curve query or overlap proof.

Why: A box against the smooth three-section loft reached certified section branches, then stopped at conic-only curve elimination or an off-carrier atlas side query. These queries must keep the finite carrier and selected component as part of the definition.

Migration: No signatures, formats, tolerances or work allowances change. Callers handling Boolean refusals must accept the new code. Spline/cutter polynomial meetings require rational homogeneous coefficients. Whole-span containment without a finite second constraint and partial-span containment refuse by name. These query additions do not claim complete spline face imprinting, a sewn loft Boolean or STEP output for section branches.
