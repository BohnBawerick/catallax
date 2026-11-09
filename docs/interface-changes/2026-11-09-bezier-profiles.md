## 2026-11-09: polynomial spline profiles

Pieces: shapes, recipe, features, topology, intersections, Boolean operations, mesh, measurement, application

`catallax_recipe::Segment::Bezier` adds a quadratic or cubic polynomial span. Its start is the preceding segment endpoint. `controls` contains one or two interior control points and `to` is the endpoint. The SolidDraft converter accepts `type: "bezier"` with explicit `from`, `controls` and `to` coordinates. Consecutive endpoints retain the existing application join tolerance.

Closed profiles may combine these spans with lines, including nested holes. Exact control-hull subdivision proves that interiors do not cross or touch. Exact polynomial winding determines nesting and Green integrals determine signed area. Collinear monotone spans become lines. Other curve families in the same profile return `sketch.bezier-curve-mix`; an exhausted separation proof returns `sketch.bezier-proof`.

Normal extrusion retains the authored polynomial as a ruled B-spline surface. Each span must have distinct endpoints, nonzero endpoint tangents and a proved monotone chord projection. A span without that embedding proof returns `shapes.bezier-extrusion-embedding`. Subdivide such a span into regular monotone spans before submitting it. Oblique extrusion returns `extrude.bezier-oblique`.

Finite patch charts support exact membership, normal queries, rigid transport and reflection. Planar caps normal to the extrusion produce exact finite section curves. Section-to-line roots include tangencies and endpoints. Section-to-section queries prove complete separation except at shared original endpoints; interior intersections and unresolved subdivision return `intersect.bezier-curve-pair`. Parallel section-to-wall queries use the same proof and check finite axial bounds. Other wall directions return `intersect.bezier-curve-surface`. Bounded de Casteljau subdivision can prove an oblique plane disjoint even when the original control hull crosses it. Unproved oblique plane sections return `intersect.bezier-plane-oblique`.

Boolean subtraction supports a spline prism contained laterally within a planar block, including through cuts, pockets and reversed extrusion directions. Axial wall partitions require complete span boundaries, exact interval coverage and opposite uses of every interior section. More general Boolean intersections remain subject to the named finite-query refusals above.

Display meshes use certified chord bounds and shared boundary samples. Body volume, centroid and inertia use exact polynomial integrals; surface area retains enclosed integration. Bounds use polynomial extrema. Complete finite patch equality is proved from matching retained boundaries, or disproved by a corner outside the other patch. Matching corners alone return `shapes.bezier-extrusion-equality`.

Revolution, path sweep and first-hit extrusion limits return `features.bezier-profile`. Field compilation returns `field.bezier-profile`. STEP export of these new finite patches remains unsupported and returns `step.unsupported-geometry`. No operation replaces an unsupported exact profile with a polygon.

Exhaustive recipe segment and spline support matches must handle the new variants. Existing line, circle, arc and ellipse formats remain valid.
