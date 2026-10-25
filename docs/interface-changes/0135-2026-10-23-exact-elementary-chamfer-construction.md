## 2026-10-23: exact elementary chamfer construction

Pieces: features, field, comparison

Change: recipe chamfers construct exact planar strips on plane-plane lines and conical bands on perpendicular plane-cylinder rings. They reuse the edge-blend boundary replacement machinery from elementary fillets. `Built::chamfers` retains each `ChamferDefinition`, including ordered supports, Euclidean setbacks, original edge, material sense, contacts and analytic carrier. Body patterns and mirrors transform these definitions with their persistent origins. The field independently resolves sharp material contacts and constructs triangular extruded or revolved affected regions. Comparison now attempts chamfer recipes.

Unequal setbacks follow the two input origin patterns, even when their order differs from topology storage. Ambiguous ordering and conflicting overlapping references refuse explicitly. Separation checks cover the finite affected region, including the original sharp edge, rather than the unbounded bevel carrier. Changed rails and inherited edges must retain their direction; uncertain overflow refuses before publishing a body.

Full revolutions of convex line polygons use certified section half-planes and signs of radical expressions for field classification and face membership. Convexity and the side of the axis are established from the exact profile; other profiles retain the general winding evaluator. This avoids isolating the same radial square root at every display query without changing geometric or comparison tolerances. Comparison uses the upstream certified trimmed-face distance queries; residuals retain their existing meaning.

Migration: no recipe schema, command, work allowance, modelling tolerance, legacy answer or STEP-option changes. New failures are `chamfer.non-positive-distance`, `chamfer.ambiguous-order`, `chamfer.conflicting-distances`, `chamfer.lost-edge`, `chamfer.unsupported-contact`, `chamfer.unsupported-topology`, `chamfer.unsupported-overflow` and `field.chamfer-contact`. General roll-over, three-edge vertex patches, trimmed circular arcs, filtered field references and nonrational field sections remain implementation gaps. SolidDraft's positional-selector converter is unchanged. Sampled comparison is separate from a whole-result geometric certificate.

