## 2026-11-08: fillets over concave blends at mixed vertices

Pieces: recipe definition, features.

Change: RECIPE section 4.5.1 now defines a vertex where selected convex and concave edges meet. The concave blends form first and the convex sheets roll over them. Each concave blend leaves a boundary arc against the vertex's other faces, and the convex sheets continue along the tangent chain this arc forms. The patch that rolls over the concave blend is the vertex's junction sector and is named `corner.vertex`. When every edge at each mixed vertex is selected, the B-rep builds the operation in a concave pass and then a convex pass. Both passes reuse the ordinary strip construction and its proofs. The convex pass takes the selected convex edges and the arcs at mixed vertices.

The strip guard also proves a face clear of a strip when the face is a plane bounded only by straight edges and a hull face separates its corners from the strip's hull. A plane face lies in the convex hull of its corners, so this is exact. Before, the guard tested only the unbounded support, which refused `l-mixed-all` because the wall-side plane crosses the hull of a base edge's strip far from the wall-side face itself.

Why: the definition left a mixed vertex open. Under the plain set reading, the concave addition and the convex removals meet at a knife edge, which no CAD kernel produces. OCCT's simultaneous answer for `l-mixed-all` equals forming the concave edge first and then rounding every remaining sharp edge. The raw volumes agree to 8e-10 mm^3. The rule needs no new syntax and keeps one constant radius per operation.

Refusals: a mixed vertex that leaves an incident edge unselected keeps the strip construction's refusal. An edge the construction cannot classify as convex or concave also keeps its earlier route.

Migration: no syntax or caller changes. Selections without a mixed vertex take the unchanged single pass. The torus at each mixed vertex of `l-mixed-all` carries `corner.vertex`, not a blend name that refers to a same-step blend.

Evidence: the unchanged `l-mixed-all` recipe constructs 38 exact faces (18 cylinders, 10 spheres, 2 tori, 8 planes) at 1M work, then STEP export exhausts 1M, as for `box-all`. At an explicit 20M allowance the complete request uses 2,997,397 work and writes one raw-valid 38-face solid whose volume is 7.8e-10 mm^3 from OCCT's. Ten material witnesses pass in the native classifier and in OCCT's independent construction. Full-corpus results belong to the clean published head's report.
