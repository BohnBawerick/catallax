# M10: finite spline imprint queries

The box against a smooth three-section loft now passes the native curve comparisons and outgoing germ queries. It still refuses while finding a strict interior probe for a face piece. The loft solid, its classification and sewing remain unproven.

All capability requests below use the unchanged 20,000,000-unit allowance. No tolerance or proof-search allowance changed.

| Loft cut | Operations | Work | Result |
| --- | --- | ---: | --- |
| Top | subtract, union, intersect | 5,743,555 | `boolean.spline-face-probe` |
| Side | subtract, union, intersect | 8,276,808 | `boolean.spline-face-probe` |

The refusal names the missing exact parameter-space path and boundary-clearance certificate. Retained branch definitions and exact radial order do not by themselves supply that face-interior proof.

## Proved queries

Native spline-line meetings cover every knot span and verify both line constraints. Shared-knot meetings are deduplicated exactly; an overlap interval names its missing ownership proof. Native cap comparisons reuse the patch-pair engine's strict closed-box exclusion certificate.

Outgoing acceleration uses homogeneous derivatives through second order, including both weight derivatives. Native curves require a unique finite preimage or an exposed-endpoint certificate. Plane-section branches check every retained chart preimage and select the outgoing cell side. Ambiguous selectors and unsupported cutters refuse by name.

The fitted Boolean's existing apart-body proof remains available before a new native germ query spends its allowance. A private certificate binds the original bodies and carries the main call's paid proof to the fallback. Speculative workers own separate contexts. The disjoint fitted cut therefore retains its original carriers and blend lineage.

## M8 cuts after a fitted fillet

These are the requests in [the M8 list](M8-LATER-CUTS.md), followed by canonical P5. Each built request exports STEP. Work is the charged prefix; a budget refusal can leave units when the next reservation exceeds the remainder.

| Request | Work | Result or missing step |
| --- | ---: | --- |
| k01 wall hole | 20,000,000 | `numbers.operation-budget`: operation work allowance exhausted |
| k02 top slot | 19,995,626 | `numbers.operation-budget`: operation work allowance exhausted |
| k03 rim notch | 2,850,911 | `intersect.parametric-surface`: circle with fitted rolling-envelope spline: this intersection needs a finite parametric surface solver |
| k04 torus pocket | 3,282,910 | `intersect.traced-pair-required`: skew cylinder/torus intersection requires certified branch tracing |
| k05 plane notch | 3,424,555 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k06 cylinder drill | 3,425,416 | `intersect.unsupported-pair`: fitted rolling-envelope spline with cylinder: surface pair is outside piece 4 coverage |
| k07 cone drill | 3,927,165 | `intersect.traced-pair-required`: non-coaxial cone/cylinder intersection requires certified branch tracing |
| k08 sphere dimple | 2,828,404 | `intersect.traced-pair-required`: non-coaxial cone/sphere intersection requires certified branch tracing |
| k09 ring groove | 2,830,766 | `intersect.unsupported-pair`: fitted rolling-envelope spline with torus: surface pair is outside piece 4 coverage |
| k10 axis section | 2,825,497 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k11 offset section | 3,419,097 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k12 disjoint tool | 12,916,559 | solid, 22 faces |
| k13 fillet before the countersunk pattern | 19,997,484 | `numbers.operation-budget`: operation work allowance exhausted |
| Canonical P5 | 1,737,135 | solid, 22 faces |

## Validation

The Shapes, Intersect and Boolean library suites pass all 199 tests. Both later-cut integration tests and all three box/loft integration tests pass. Formatting and Clippy pass for the changed crates and their test targets.

The six loft requests and fourteen M8/P5 requests have identical comparable manifests, exported bytes and work with worker/kernel counts 1/1, 2/2, 4/4 and 1/12. The last configuration uses all available kernel threads. Both previously proven M8 answers and STEP files match their saved versions exactly.

Every built request exports STEP and opens in OpenCASCADE with the expected 22 faces. The raw, unhealed shapes are also valid and retain 22 faces. The loft cases export no solid.

Raw evidence retained privately.
