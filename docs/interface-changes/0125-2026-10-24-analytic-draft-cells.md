## 2026-10-24: analytic draft cells

Pieces: recipe, features, field, application, comparison

Change: `Draft` accepts optional `keep_regions`. Omission keeps neutral-plane-connected material; an explicit empty list selects no material; nonempty lists select components by strict interior membership. Existing version-3 recipes serialize without a new property. Feature-pattern placement transforms saved region points along with the pull and neutral plane.

`DraftDefinition` retains exact pull, neutral plane and algebraic angle values. Plane replacements follow the neutral-line generators. Pull-parallel cylinders reduce to a single cone nappe only when the neutral normal is parallel to pull. Oblique cylindrical neutral planes return `draft.unsupported-oblique-neutral`, an unfinished-coverage refusal, not an impossible-input refusal. The neutral-curve definition is unchanged.

The cell constructors rebuild convex planar arrangements from exact support intersections instead of preserving old adjacency. Redundant and zero-area patches disappear. Axial round cells trim their caps or replace a swallowed cap with the cone apex, without adding a face. Surviving faces retain their source origins. `Built::consumed` and `Plan::consumed_faces` expose removed-face provenance separately from selectable names. Field name lookup prefers the result's replacement candidate when the origin also occurs in construction history.

The field compiles the same analytic definitions without constructing B-rep topology. Its cell distance uses projection to trimmed planar facets and their edges, or projection to the round cell's meridian and cap disks. Comparison now attempts draft builds instead of classifying every draft operation as unsupported.

Coverage: convex planar cells and axial cylindrical cells with complete perpendicular caps. Nonconvex, mixed curved and Boolean-input arrangements and filtered field face references still have explicit unfinished-coverage refusals. These constructors do not by themselves complete roadmap item 19 or establish the drafted P1 latency and comparison exits.

Prior art: [OpenCASCADE's cylindrical draft construction](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKOffset/Draft/Draft_Modification_1.cxx) explicitly checks neutral-normal/pull parallelism before reducing to a cone. Its topology-preserving restriction is not adopted for swallowed planar neighbours or round caps. No source code is copied.

Migration: callers matching `Operation::Draft` exhaustively must bind or ignore `keep_regions`. No version, budget, tolerance, legacy question or STEP-option change.

