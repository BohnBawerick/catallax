## 2026-10-24: analytic distance features and query-only trim candidates

Pieces: shapes, features, field, measure and compare
Change: `Surface::parallel_support` constructs exact plane, circular-cylinder, cone, sphere and ring-torus parallels. A consumed radial support returns `None`; an unsupported derived or spindle carrier retains a typed refusal. Features expose bounded convex plane arrangements, complete cylinder and radial band definitions, and face/edge/vertex distance pieces shared with the independent field compiler. Ordinary recipe evaluation now handles the analytic offset, shell and thicken families listed in RECIPE.md section 7.15, retaining origins, input marks, empty erosions and consumed cavities.
Why: distance features need actual trim witnesses and explicit edge/corner geometry, not sharp support movement or residual thresholds.
Migration: recipe and exchange formats do not change. Unsupported legal families remain explicit coverage refusals. The existing certified elliptic-offset prism entry and its local tolerance policy are unchanged.

Change: `Plan::visit_distance_candidates` adds independent, query-only trim-aware candidates. Comparison invokes it only for recipes containing offset, shell or thicken. `Plan::select` prefers a live replacement support with the requested origin before an older registry candidate.
Why: narrow open rims, inner walls and clipped circular patches can have no sample on a fixed primitive grid.
Migration: callers must certify candidate trim membership, origin sets and directions. Sampling does not prove whole-boundary coverage and adds no work to part construction.

Change: property queries on finite curved distance-feature trims return `measure.budget-exhausted` before entering uncovered quadrature. Exact planes, whole analytic supports and full-ring reductions retain their existing integrals. The measurement crate uses typed recipe provenance to identify this property-proof boundary, including copied distance faces.
Why: an offset's distance and incidence proofs do not prove its volume or moment integral. A typed refusal preserves the valid geometric result without a stalled property request.
Migration: treat this refusal as unfinished property coverage, not invalid geometry or permission to accept a sampled value. No work limit, accuracy, modelling tolerance or approximation cap is relaxed.

