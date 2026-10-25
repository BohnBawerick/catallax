## 2026-10-23: closed trim contacts and chart-local display proofs

Pieces: field, comparison
Change: `Plan::closed_face_query_checked` refines conservative named-face queries using exact transverse support sectors and planar fans. Reverse comparison uses these closed-patch proofs. Tangent, singular, ambiguous-partition and larger curved contacts remain conservative. Native field display restricts rational constraints to its own rational charts, proves open-cell signs with exact Bernstein coefficients, and verifies rational trim rails before splitting. A positive constraint excludes an analytic inequality `Cell`; profile `Feature` winding semantics remain unchanged.
Why: valid closed circular rails were unresolved during reverse comparison, while world-coordinate interval dependency forced excessive native display refinement.
Migration: no request or manifest format, tolerance or budget changes. Native OBJ bytes, cells, vertices and measured work may change; retain the reported two-sided bound. The ordinary named-face query remains available with its conservative behavior. These proofs do not construct a B-rep or promote sampled agreement into a whole-result certificate.

