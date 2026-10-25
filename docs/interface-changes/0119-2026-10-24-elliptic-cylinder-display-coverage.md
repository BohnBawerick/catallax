## 2026-10-24: elliptic-cylinder display coverage

Pieces: display mesh

Change: elliptic cylinders now use their exact periodic/linear chart through the existing shared-edge display path. Inverse projection scales the transverse coordinates by the semi-axes. The major semi-axis bounds chart acceleration; the major/minor ratio bounds normal rotation. Facet-normal and edge-length controls retain their existing certificate semantics. Two angularly monotone boundary loops can form an axial triangle strip without interior samples; every candidate triangle must pass the existing chord and quality bounds, or the periodic triangulator supplies refinement.

Why: SolidDraft must display all six analytic support kinds without replacing exact geometry. Finite parabola and hyperbola edges already use certified traces; new closed-solid regressions exercise both forms alongside elliptic caps, oblique circular caps and rigid placements.

Migration: no API, face or edge picking IDs, recipe role, tolerance, budget, manifest or STEP changes. Edge polylines still contain the single shared vertex sequence used by adjacent faces. Display bounds remain bounds against the analytic definition, not a claim that polygonal trims equal curved trims.

