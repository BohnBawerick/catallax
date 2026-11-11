## 2026-11-12: Keep spline incidence facts within their proved scope

Pieces: 2, 3, 4, 8

Change: spline incidence distinguishes an off-surface edge, a certified finite edge and a certified whole retained curve. The STEP builder records `Incident(Curve, Surface)` only for the whole-curve result.

Why: the [spline admission route](2026-11-12-step-spline-certification.md) checked only spans crossed by the physical edge. Unused retained spans could leave the surface, so that proof could not authorize a global incidence fact. A full-curve intersection query previously consumed the fact and returned an incorrect exact `Contained`.

Migration: public signatures and serialized refusal variants remain. A finite-edge-only certificate does not authorize a full-curve containment answer. Queries or body checks that need unsupported full-curve evidence may refuse. Whole-curve certificates, exact line reductions and checked complete patch boundaries retain their existing global facts.

The regression builds a planar face with a degree-one spline spanning parameters zero to two, while its physical edge uses only zero to one. A point at parameter three-halves is proved off-plane. The retained curve stays untrimmed; its full-curve query refuses and no global incidence fact is recorded. A separate control verifies the whole-curve certificate when every retained span is proved on a surface.
