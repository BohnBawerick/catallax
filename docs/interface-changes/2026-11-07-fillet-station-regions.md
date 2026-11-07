## 2026-11-07: exact fillet station regions

Pieces: features, recipe consumers, topology, Boolean, STEP.

Change: the B-rep evaluates the `keep_neighbours` syntax introduced by [the schema prerequisite](2026-11-07-fillet-overflow-schema.md) for one straight plane/plane run with explicit limits at both endpoints. It builds the exact affected-region prism, subtracts it for a convex edge or adds it for a concave edge, and retains the original rolling-ball definition. The cylinder and two station-plane closures share their exact boundary edges. Closures use the original endpoint `corner.vertex` names; the cylinder keeps `blend.edge`.

Why: accepting the syntax alone did not construct the requested stationary closures or preserve unselected neighbours. The existing exact tangent-triangle hull encloses the whole affected region. Strict separation from every other input support proves that those neighbours keep their complete trimmed point sets. The station planes remain exactly normal to the input edge, at the authored distances. No fit, radius change, moved station, work-allowance change or roll-over fallback is introduced.

Refusals: `fillet.incompatible-setbacks` reports a limit on the wrong edge or endpoint, an ambiguous endpoint name, conflicting distances at one endpoint, or stations leaving no positive interval. Missing limits, curved/multiple runs and unproved neighbour clearance retain `fillet.unsupported-overflow`. `fillet.invalid-station-region` reports a failed topology check on the unpublished affected region. These are not a general impossible-geometry verdict.

Migration: no new syntax or caller argument is required. Existing default and explicit `roll_over` paths remain unchanged. The independent field compiler still explicitly refuses `keep_neighbours`; local whole-face field measurements are separate from independent whole-recipe compilation. The corpus changes only the demonstrated setback expectation and adds material witnesses. Its dimensions, radius, station distances, edge and OCCT window are unchanged.

Evidence: the selected ordinary recipe writes a raw-valid nine-face solid at 1M work and matches the independent station-window volume of 3996.566370614358 mm^3. Its cylindrical blend and both closure planes have measured whole-face deviation zero, and all three names resolve after a radius edit to 1.5 mm. Material regressions cover convex removal, concave addition, asymmetric/reversed station records, original end walls and incompatible constraints. Full-corpus results belong to the clean published head's report.
