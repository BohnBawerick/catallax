## 2026-10-25: output-bounded planar boundary samples

Pieces: mesh and application

Change: planar boundary samples retain a shared trace and an exact parameter,
not an independent exact-coordinate expression graph. Interval enclosures feed
the existing planar filters; ambiguous predicates construct the same exact
coordinates on demand. Boundary storage uses the existing `MeshOptions.max_vertices`
limit, including topological and trace-piece endpoints. The fixed aggregate
4096-sample refusal is removed.

Why: the sample count did not measure retained memory or arithmetic cost and
rejected ordinary fine-chord drawings. Shared trace definitions keep boundary
storage proportional to the generated vertices without retaining every exact
construction graph.

Migration: no API, mesh tolerance, vertex allowance or arithmetic allowance
changes. Callers must continue handling resource refusals. A large-circle
application build now exhausts the existing aggregate arithmetic allowance and
returns `application.work-limit`, rather than the removed
4096-sample `mesh.resource-limit`. Tiny-coordinate region refusals and the
isolated regression process caps remain unchanged.

