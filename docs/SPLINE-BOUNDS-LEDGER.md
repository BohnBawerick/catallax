# Certified finite spline-prism bounds

The reproduction base is main
`866c9a30c4730ebe906cba352f61bf350b251652`. The initial checkout was clean. No
`AGENTS.md` existed in the checkout or workspace ancestors. The checkout's
verification skill, feature map, PLAN, ARCHITECTURE, relevant RECIPE contracts
and interface-fragment instructions govern this change. Only public synthetic
fixtures are used. Final source/runtime identities and completed qualification
are recorded in the draft PR and external evidence.

Production changes are confined to measurement bounds and its new private
spline module. Dedicated measurement tests consume the existing
[certified recipe](../.claude/skills/verify-catallax/references/examples/certified-recipe.json)
through `catallax_app::build_certified`. Application, shapes, topology,
classification, distance, display and exchange implementations are unchanged.
The measurement development dependencies add the real application fixture;
they do not add a production dependency on the application.

The intended integration base is
`b4617972468bc69cccdcae9e59df60e784704350`, after the disjoint planar-face
thickening and application-loop changes landed. The original reproduction and
its binary remain separately identified above; qualification is repeated on
the integrated source rather than attributed to that earlier executable.

| Unit defect or control | Executed before repair | Repair / qualification |
| --- | --- | --- |
| Closed finite spline rims rejected as open edges | A current-base application build and `SplinePrism::check` succeed; `bounding_box` and all three `face_bounding_box` calls return `measure.unbounded-face`, “open edge lacks finite endpoints”. The regression fails with exit 101. | Admit only the complete checked prism domains, retaining endpoint-free closed rings. The same body/face regression must answer without accuracy marks. |
| Missing represented spline extrema | The base refusal prevents any extrema answer. Source inspection alone is not counted as reproduction. | Split nonzero knot spans, use homogeneous interval de Casteljau subdivision and attained endpoint witnesses, and require all six extremum widths. A loose control hull cannot terminate the query. |
| Excess query work | Intermediate implementations exhaust the unchanged 20M request allowance, including a fixed 256-bit control-net pass. These failed controls are retained in external evidence. | Choose working precision from the requested accuracy, use a decider-aware seam jet and carry interval control nets/witnesses through subdivision. Width checks remain authoritative. The same 20M regression answers at 1,765,165 work units for absolute accuracy 1e-8. This is not a 1M result or a wall-time speed claim. |
| Finite cap trim | Each base cap refuses independently. | Bound each cap by its own complete rim, rather than both ends of the carrier; the lower cap excludes the upper carrier's coordinate extrema. |
| Unproved face trims | No new coverage is inferred from a carrier's finite active knot interval. | Bare side, partial side and detached duplicate cap controls refuse. A failed or exhausted query preserves the immutable built body. |
| Reflected raw topology copy | An exact reflection passes structural `Model::check`, but the unchanged `SplinePrism::check` rejects its cap-ring travel with `numbers.geometric-gap-unproven` before bounds. A proper rotation/translation passes that checker and the independent extrema oracle. | Preserve the reflected-copy refusal as an explicit negative control. Extending the prism orientation proof belongs to its producer/checker; no shapes or topology repair is included here. |

The frozen base verifier SHA-256 is
`660956d7076d6229f367398f6271d31762df55fa5fe587bb56cddf59566c4742`.
Doctor reports the exact base commit and `catallax-request/1`. All three public
verification examples answer on that base, with their original requests.
The separately preserved base library regression executable SHA-256 is
`50c13a26decfb2f67a475c6f5f796998f54573174cc0fcb36ceddae1afee8a2f`.
Its command records identify the release libraries used before production edits.

The first Cargo test attempt hit a missing `m4` build prerequisite for the
existing GMP test dependency. A Debian build tool was extracted outside the
checkout; no system configuration or dependency version was changed. The
infrastructure failure is distinct from the executed geometric refusal.

## Proof and independent controls

`SplinePrism::check` proves one convex periodic represented sweep, its complete
shared physical rings, two perpendicular cap disks, finite domains and material
orientation. A face-only query recovers the three incident faces and checks a
scratch shell containing those original faces. It creates no incidence facts,
seams, vertices or new geometry. Domains are not guessed from endpoint samples.

The admitted carrier is composite degree-(3,1) Bezier geometry. The measurement
consumer verifies knot multiplicities before treating active controls as Bezier
controls, skips zero-width knot spans and preserves each finite endpoint. With
positive weights, a rational Bezier segment lies in the Cartesian convex hull
of its homogeneous control net. Interval de Casteljau subdivision encloses the
exact restricted net. The hull supplies outer bounds; attained endpoints and a
seam jet supply the opposite bounds on each extremum. Only panels that can
affect an unresolved global minimum or maximum are refined. Every accepted
minimum/maximum enclosure satisfies the existing absolute/relative `Accuracy`
rule; insufficient work refuses rather than publishing an unconverged hull.
Working precision respects the caller's maximum and is selected from its
requested tolerance with guard bits. The query never weakens `Accuracy` to
compensate for rounding or exhausted work.

Every side coordinate is affine in the checked v sweep, so its extrema occur
on the two rims. A linear coordinate on either compact planar disk attains its
extrema on its own boundary, or is constant. Thus cap, side and body enclosure
and extremum-width proofs follow from the rim proofs without display padding.
The generic traced-branch extrema algorithm is unchanged.

The independent test oracle converts each cubic Bernstein coordinate to its
power basis, solves its quadratic derivative algebraically and evaluates every
admitted root plus both endpoints at 512-bit precision. It never calls spline
evaluation, jets, tracing or measurement bounds. It checks enclosure and
requested extremum error separately. Controls include a global extremum strictly
inside a knot span, repeated knots, positive nonconstant rational weights,
negative and zero coordinates, exact placements, scale changes, an offset edit,+a shorter axial trim and a finite interval excluding the carrier maximum.
Original face origins and physical rings remain present.

The fitted extrema are also compared with the independently written support
function of the retained exact inward parallel ellipse. A generic placement
distinguishes the represented and exact-target extrema. Their discrepancy is
bounded by the existing local correspondence certificate, not by silently
substituting that target law for the spline. The additive
[interface fragment](interface-changes/2026-10-26-certified-spline-prism-bounds.md)
states the explicit composition rule for a consumer requesting target bounds.

## Commands and limits

```sh
cargo build --release -p catallax-cli --bins -j 2
cargo test --release -p catallax-measure --test spline_bounds -j 2 -- --nocapture --test-threads=2
cargo test --release -p catallax-measure --lib bounds::spline -j 2
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release -j 2 -- --test-threads=2
```

CLI preservation uses the checkout verification driver on the three public
examples, retains base/candidate doctor identities, manifests and actual
STEP/OBJ bytes, and compares fresh one/two/four-worker and warm-cache runs.
The existing consumer smoke remains separate from the new public library
bounds regression; no new verification request adapter is claimed.

The 1M scoreboard, 20M application and 100M verification allowances remain
distinct. Modelling tolerance, local fitting targets/caps, native field
accuracy/cell limits, printing/mass defaults and independent vertex limits are
unchanged. Shared-cloud timings are diagnostic, not a quiet pinned AB/BA speed
claim or complete-click acceptance. The 100/500/1000/2000 ms engineering targets
are unchanged. Private release preservation and readiness remain external
qualification gates; this public fixture work does not replace them.

## Cross-body sharing correction

Frozen source `eda40cf1e94c6b06409b252cfeeeac0a5cdb5bb8` incorrectly used the
model-wide coedge owner graph even when the caller supplied a body. The earlier
qualification is retained as historical evidence and does not qualify this
correction. The failure is an acceptance gap, not a false enclosure.

The public regression starts with checked prism A, adds a fresh coedge, loop
and lower-cap face on the same ring and plane with the same senses, then forms
prism B from that cap and A's existing upper cap and side. Both bodies pass the
structural and complete geometric prism checks. On the frozen source, both
body-bound queries return `measure.unsupported-surface`; the regression exits
101 before repair. Its original source, executable, command, output and hash
remain separate from the corrected-head evidence.

Body queries now read actual face/loop/coedge uses within that body's inventory
to build their adjacency map. They do not filter model-wide first-owner links:
even a complete loop and coedge may be shared across bodies, with a different
first owner in the model's lookup table. A query's model, optional body and
Accuracy are immutable throughout its cached face/edge extrema. The same
complete-domain checker, closed rims and unchanged interval extremum proof are
used after the scoped walk. No geometry, incidence facts or storage change.

Four additional application regressions check the valid A/B bodies against the
independent rim-extremum oracle, shared loops/coedges in all six face orders,
unowned foreign coedges, missing-cap/side refusals, model immutability and
alternating loose/tight Accuracy requests. A face-only query still has no body
selector and explicitly refuses the ambiguous four-face model-wide cover.
Incomplete body inventories cannot borrow a globally available cap or side.
No budget, tolerance, target-law guarantee or existing refusal is weakened.
See the [additive correction contract](interface-changes/2026-10-26-body-scoped-spline-prism-bounds.md).
