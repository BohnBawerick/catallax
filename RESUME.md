# P1 resume state

## Closeout decision

Further P1 bottleneck work is paused pending a new kernel plan. This archival branch, `fm/catallax-p1-prototype`, contains the unqualified canonicalization prototype and this record. It has no pull request and is not production code. It does not claim P1 or milestone M2 complete.

## Reviewable code

- https://github.com/BohnBawerick/catallax/pull/144 merged the prepared-region cache reuse and compatibility-preserving circle/cylinder exclusions.
- https://github.com/BohnBawerick/catallax/pull/148 contains the next exact circle-plane and support-cylinder/plane exclusions. Its geometry head is `5f80f2906ddb243cfba15d2153f5a4118fd66e20`, based on `c8d301b7de1b84680c70f332c180a112614cc238`. Full-six qualification was pending at closeout; no result is claimed here.
- The latter code passed 183 release intersection/Boolean tests, with one existing long-run qualification ignored, and formatting checks.
- Its 41-case public replay preserved all 35 prior answered manifests, including exact property endpoints, after excluding kernel identity/statistics. All 35 STEP files were byte-identical and passed calibrated standard-reader validity with expected face and solid counts. This is not a raw no-healing or volume-containment claim.
- All six unsuccessful public requests remained refused. Case `eh-20261020-0000354` reached `s2` sewing at 73,232,266 accepted work and reported `boolean.self-touching-shell` instead of exhausting 100M; it published no body or output. The other five refusals were unchanged.

The requested complete P1 recipe and the boss-radius 8-to-9 edit remain unchanged. Construction still refuses before the remaining slot cuts and bridge complete.

| Measured code | Original accepted work | Radius-9 accepted work | Failure |
| --- | ---: | ---: | --- |
| Reviewable exclusions, real CLI | 999,987 | 999,828 | `slot-cut-1` |
| Unqualified canonical-pool prototype | 999,670 | 1,000,000 | `slot-cut-1`, `Intersect` |

Neither run delivered a complete P1 body. Edit naming checks, selectable unproved preview delivery and background proof remain unreached. Preview wait and proof time are unmeasured, not zero. Development evaluation waits varied substantially with logging and shared-host load; they are not latency qualification.

## Measurements and accounting

The [proof-cost report](docs/milestones/2026-10-26-p1-proof-costs.md), [baseline stages](docs/milestones/2026-10-26-p1-baseline-stages.csv) and [plane-exclusion stages](docs/milestones/2026-10-26-p1-plane-exclusions-stages.csv) retain the shipped measurements. Refused stages were captured with temporary `OperationSpan` probes before error propagation; successful-only `StageReport` data is insufficient.

Accepted operation work and cumulative attempted numeric reservations are different counters. An unaffordable arithmetic reservation can increase numeric statistics without becoming accepted work. No allowance, tolerance, geometric definition or output accuracy was changed.

The exact plane exclusions reduced matched `slot-cut-0` intersection work from 279,884 to 131,766 in the original, and from 275,831 to 127,713 in the edit. They retained support-pair enumeration and nonempty RUR witnesses. The initially suggested duplicate planar-incidence proof accounted for only 1,035 work across the reached chain and was not the first optimization selected.

## Unmerged prototype: observed, not accepted

Two related private changes were investigated:

1. Reuse certified coordinate enclosures in graph, curve-meeting and imprint canonical pools. Keep first-existing witness order, meter cached comparisons, retain literal-fact and scalar fallbacks, and restrict each cache to one append-only model.
2. For points with six proven current-store incidences to the same three rational cylinder supports, use the existing budgeted `Decider::certify_root` uniqueness certificate on a hull containing both points. The hull padding is a certificate interior, not a coincidence tolerance. No new point expressions are created.

An isolated fixture reproduced the original branch/cylinder and line/cylinder routes for four points at each boss radius. Scalar equality cost 26,042–26,109 accepted work per equivalence. Existing uniqueness certification cost 11,027–11,094. These were independent profiling requests, not a split or reset of P1's allowance; no wall-time win was qualified.

In the full P1 prototype, four such identities completed. Reached canonicalization work was 183,093 original and 175,995 edited, despite reaching later queries. Refined decisions in the refused intersection stage fell to 174 and 176. Guarded enclosure setup increased the prefix before `slot-cut-1` by 934 work. This is later progress within the same refusal, not completed construction.

The focused regression run was **2 passed, 1 failed**:

- Passed: independent RUR routes retained the four existing witnesses under a bounded request.
- Passed: two incidence facts did not discharge a third cylinder for a nearby non-root.
- Failed: `two_subnanometre_roots_in_the_same_certificate_hull_stay_distinct` exhausted the unchanged 1M allowance. The failure was a resource refusal, not an observed wrong merge. The prototype therefore has not passed the full suite or prior-manifest/STEP audit and must not be treated as landable.

Latest refusal attribution:

- Original: point `#2639`, source `CurveSurface(#383, #14)`, failed during canonicalization after 1,727 accepted work in that intern. The complete request used 999,670.
- Edit: curve pair `#380/#390` refused after 4,335 accepted work in that pair; the request used 1,000,000.
- These handles identify one retained run, not stable geometry names. Reconstruct their definitions before choosing another optimization.

## Recovery material

The complete seven-file unqualified prototype, including development probes and the failing regression, is archived on `fm/catallax-p1-prototype`, based on geometry head `5f80f2906ddb243cfba15d2153f5a4118fd66e20`. The qualified review branch is unchanged. The original local stash snapshot is `aa5763dc7d792a1ab05ccd1037f102709d933821`. The task evidence bundle also retains the portable patch `closeout/canonical-prototype.patch`, recovery metadata `closeout/prototype.json` and a copy of this record. A fresh clone can fetch the archive branch; the local stash object is not assumed to exist remotely.

Relevant evidence names within that bundle:

- `qualification/isolated-equality-baseline.json`
- `qualification/isolated-equality-uniqueness.json`
- `qualification/identity-route-p1.json` and `qualification/identity-route-profile/`
- `qualification/canonical-pool-p1.json` and `qualification/canonical-pool-profile/`
- `qualification/canonical-pool-sources-p1.json`
- `qualification/rebased-exclusions-versus-qualified-audit/`
- `qualification/rebased-exclusions-versus-qualified-export-gate.json`

The portable development fixture is retained with the evidence. Earlier full point-solving and finite-support-pruning experiments changed prior enclosures, including wider intervals; do not restore them as the starting point.

## Next rungs, only after a new plan

1. Resolve the near-singular certificate resource regression before considering the prototype for production review. Preserve strict non-merging of distinct nearby roots, missing-incidence rejection, original witness choice, typed refusals and caller accounting. Do not repin the failure as success, reset an exhausted allowance or raise a cap.
2. Profile the newly reached canonicalization and curve-pair sites from their actual definitions. Separate repeated elimination, coordinate enclosure and equality work. Reusing an exact negative or complete system result is only a research candidate; completeness, component coverage and unchanged witness expressions must be proved before skipping queries. Carrier incidence alone proves neither trim membership nor orientation.
3. Any candidate must pass the focused boundary regressions, full affected suites, identical prior successful manifests/property endpoints and calibrated STEP validity/count/byte checks before another code PR. A cheap local fixture or a later refusal is not full-P1 acceptance.
4. Then resume the entire original recipe, configured edit and naming checks, complete revision-bound selectable preview, and requested background proof. Report preview wait and proof phases separately under the [26 October contract](docs/ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26).

The application now has `WorkAllocation`, revision-bound previews and a retained-build `ProofJob` in `catallax-app/src/lifecycle.rs`. It reserves the explicit total before preview and does not rebuild or grant a new allowance for proof. Its documented proof scope is topology validation plus certified display, not a whole-model or second-kernel certificate. Keep requested output obligations explicit; comparison remains on demand. No lifecycle path was exercised for P1 because its construction never completed.
