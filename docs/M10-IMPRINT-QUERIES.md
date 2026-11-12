# M10: finite spline imprint queries

## Retained implicit boundary traces

At `5183cb22d0c1e445917a2ed582d9b8419027a28c`, retained section curves have complete traces in their own certified cells. Clipped endpoints retain exact component ownership. Monotone root bounds enclose every parameter interval; implicit-function recurrences supply Taylor jets, and exact point queries use the caller's decider. Analytic and native spline trace definitions keep their existing evaluation paths.

Spline probes may follow a certified internal chart join without treating that join as a physical face boundary. Planar paths use complete finite meetings for retained sections and native spline boundaries. Algebraic third-surface constraints use exact elimination and check a strict transverse Jacobian. Linear constraints eliminate one parameter directly. None of these changes raises the allowance or changes a tolerance.

| Loft cut | Operation | Work | B-rep result |
| --- | --- | ---: | --- |
| Side | subtract | 10,942,803 | solid |
| Side | union | 10,942,852 | solid |
| Side | intersect | 10,934,705 | solid |
| Top | subtract, union, intersect | 19,997,386 | `numbers.operation-budget` |

The top cut remains open. Complete algebraic meetings advance it beyond the rational-coefficient refusal, but the exact proof still exceeds the unchanged 20M allowance. The all-solid target remains separately ignored, not recorded as passing.

The side solids are feature-builder results, not exported STEP evidence. Ordinary recipe exchange refuses genuine spline geometry under manifest/3. The documented opt-in certified exchange also refuses these cut lofts with `features.unsupported-offset-family`, because they have no supported complete prism certificate. No STEP files were produced and no OpenCASCADE validity is claimed.

A side cut through the public M7 quadratic spline prism independently refuses `intersect.parametric-coincidence` at 16,018 work units. Its missing exact branch-component equivalence is separate from the top cut's proof cost. Saved-corpus preservation is not established by these capability runs.

## Chart-local interior probe prototype

Commit `226f3e1af6a65a00f6d1ad71042c2251025f9cbb` adds an exact path inside one retained rational chart and one knot span. The path starts strictly left of a boundary arc. Substituted cutting equations and native-curve control hulls certify clearance from every boundary edge for the complete open-start path. Positive spline weights keep its homogeneous denominator positive. A separation-bound coincidence cannot certify the path.

The first probe now succeeds in all six box/loft requests. Classification next asks for a certified trace of a retained `PatchSection` boundary. Native spline traces already retain finite knot-span ranges, but `arc::trace` has no corresponding retained-section case and reaches `query::charts`, which refuses:

| Loft cut | Operations | Work | Result |
| --- | --- | ---: | --- |
| Top | subtract, union, intersect | 6,063,036 | `intersect.spline-chart` |
| Side | subtract, union, intersect | 8,596,235 | `intersect.spline-chart` |

Fresh CLI runs at worker/kernel counts 1/1, 2/2, 4/4 and 1/24 have identical refusal manifests, complete operation decision records and work after excluding diagnostic wall times. All 24 requests retain the unchanged 20M allowance, use no cache and export no STEP. Two focused regression tests pass: a witness remains inside its finite curved patch, and an intermediate native boundary hit rejects the path even when its endpoint is beyond that boundary. The separately ignored solid target test was exercised and still fails at the named trace gap. This step has not established classification, sewing, STEP export, an M7 cut or preservation of the saved corpus.

`TracePiece` currently represents explicit algebraic expressions and their jets. A retained section is an implicit graph selected by a complete cell cover. Its checked approximation supplies a distance bound, not an exact boundary trace. The next step extends the shared trace contract with certified implicit-graph pieces, clipped endpoint ownership and certified root bounds, rather than adding a second Boolean boundary-subdivision implementation. A floating-point trace alone supplies none of those proofs.

## Finite imprint query baseline

The box against a smooth three-section loft now passes the native curve comparisons and outgoing germ queries. It still refuses while finding a strict interior probe for a face piece. The loft solid, its classification and sewing remain unproven.

All capability requests below use the unchanged 20,000,000-unit allowance. No tolerance or proof-search allowance changed.

The recipes and complete section checks are in [`m10_box_loft.rs`](../crates/catallax-features/tests/m10_box_loft.rs). The initial reproduction at `53da1a599cf3b50ded26bdb1cf3eca4d49491603` stopped at `boolean.imprint-ambiguous` for the side box (1,074,776 units) and `shapes.parametric-query` for the top box (2,449,104 units). The tables below record query source `e56b8d2e1f942cb9852a234d2225087eb8c89675`, after merging main at `166dde09daa5b87178dd97de95689873c17a1d92`.

| Loft cut | Operations | Work | Result |
| --- | --- | ---: | --- |
| Top | subtract, union, intersect | 5,743,555 | `boolean.spline-face-probe` |
| Side | subtract, union, intersect | 8,276,808 | `boolean.spline-face-probe` |

The refusal names the missing exact parameter-space path and boundary-clearance certificate. Retained branch definitions and exact radial order do not by themselves supply that face-interior proof.

## Proved queries

Retained section meetings use the other curve's complete intersections with the cutting surface, then prove membership in the selected finite component. Native spline cuts substitute the implicit equation into each span's homogeneous polynomials, clear only positive denominators and retain every exact root in the closed span. Shared-knot roots are deduplicated. Containment needs a second complete analytic constraint or returns `boolean.patch-section-meetings`; a finite sample never replaces an overlap. Coefficient products and accumulations reserve the caller's exact-arithmetic cost before execution.

Native spline-line meetings cover every knot span and verify both line constraints. Shared-knot meetings are deduplicated exactly; an overlap interval names its missing ownership proof. Native cap comparisons reuse the patch-pair engine's strict closed-box exclusion certificate.

Outgoing acceleration uses homogeneous derivatives through second order, including both weight derivatives. Native curves require a unique finite preimage or an exposed-endpoint certificate. Plane-section branches check every retained chart preimage and select the outgoing cell side. Ambiguous selectors and unsupported cutters refuse by name.

The fitted Boolean's existing apart-body proof remains available before a new native germ query spends its allowance. A private certificate binds the original bodies and carries the main call's paid proof to the fallback. Speculative workers own separate contexts. Successful journals are reused only after checking their model reads. A refused plan runs on the main call, which produces and carries its own separation certificate. The disjoint fitted cut therefore retains its original carriers and blend lineage.

## M8 cuts after a fitted fillet

These are the requests in [the M8 list](M8-LATER-CUTS.md), followed by canonical P5. Each built request exports STEP. Work is the charged prefix; a budget refusal can leave units when the next reservation exceeds the remainder.

| Request | Work | Result or missing step |
| --- | ---: | --- |
| k01 wall hole | 19,999,984 | `numbers.operation-budget`: operation work allowance exhausted |
| k02 top slot | 19,994,736 | `numbers.operation-budget`: operation work allowance exhausted |
| k03 rim notch | 19,999,997 | `numbers.operation-budget`: operation work allowance exhausted |
| k04 torus pocket | 3,282,020 | `intersect.traced-pair-required`: skew cylinder/torus intersection requires certified branch tracing |
| k05 plane notch | 3,423,665 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k06 cylinder drill | 3,424,526 | `intersect.unsupported-pair`: fitted rolling-envelope spline with cylinder: surface pair is outside piece 4 coverage |
| k07 cone drill | 3,926,275 | `intersect.traced-pair-required`: non-coaxial cone/cylinder intersection requires certified branch tracing |
| k08 sphere dimple | 2,827,514 | `intersect.traced-pair-required`: non-coaxial cone/sphere intersection requires certified branch tracing |
| k09 ring groove | 2,829,876 | `intersect.unsupported-pair`: fitted rolling-envelope spline with torus: surface pair is outside piece 4 coverage |
| k10 axis section | 2,824,607 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k11 offset section | 3,418,207 | `intersect.unsupported-pair`: fitted rolling-envelope spline with plane: surface pair is outside piece 4 coverage |
| k12 disjoint tool | 12,915,677 | solid, 22 faces |
| k13 fillet before the countersunk pattern | 19,997,226 | `numbers.operation-budget`: operation work allowance exhausted |
| Canonical P5 | 1,736,245 | solid, 22 faces |

## Validation

The Shapes, Intersect, Boolean and Topology library suites pass all 224 tests. All three later-cut integration tests and all three box/loft integration tests pass. Six application carrier contracts, five native curve contracts and the shared-state inventory also pass. Formatting and Clippy pass for the changed crates and their test targets.

The six loft requests and fourteen M8/P5 requests have identical comparable manifests, exported bytes and work with worker/kernel counts 1/1, 2/2, 4/4 and 1/24. The last configuration uses all available kernel threads in this run. Both previously proven M8 answers and STEP files match their saved versions exactly.

Every built request exports STEP and opens in OpenCASCADE with the expected 22 faces. The raw, unhealed shapes are also valid and retain 22 faces. The loft cases export no solid.

## Saved replay

Before integrating main at `01c0c95d4181248b7115cb41e82e8d2b9f747719`, query commit `7bc0a96e69cfc58820cabd4ba952be763a2791e2` and its base `8aa58c98d6504b175c923b0acbd5e3c334e1fb52` agree on all 477 saved cases: comparable manifests, exported bytes and work records. Both runs have 460 answers, 8 `boolean.self-touching-shell` refusals, 9 budget refusals and no crashes.

The three older budget refusals already answer on that base: case 194 uses 46,980,708 units and has 21 faces, case 330 uses 47,901,314 and has 6 faces, and case 378 uses 72,438,464 and is empty. Both solids pass normal and raw OpenCASCADE checks with the same face counts. All 457 earlier proven answers and exported files remain unchanged.

The saved corpus retains its existing 100M request allowance. The loft and M8 capability checks use 20M.

Raw evidence retained privately.
