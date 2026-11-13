# M19 adaptive printing verdict

Baseline: `07d545b1b63aa39a09cc422c5540e56e1dd07d32`, local main on
13 November 2026. This is a bounded release verdict, not M19 completion.
General material-topology proof, sheet printing, full P6 printing and calibrated
cost estimates remain unestablished. The [adaptive failure report](../interface-changes/2026-11-13-adaptive-print-progress.md)
now retains actual candidate and certificate progress. No acceptance rule,
tolerance, capacity, work allowance or saved-answer approval changes.

## What main already supplies

[The adaptive increment](../interface-changes/2026-11-11-adaptive-print-and-sheets.md)
landed in [the original change](https://github.com/BohnBawerick/catallax/pull/405).
It supplies exact sheet point queries, an opt-in adaptive request with retained
cell bounds, staged failures and one admitted nonempty geometry family, a
single convex block. Proven-empty plans are also admitted. The
[shared worker change](https://github.com/BohnBawerick/catallax/pull/413) changed
worker reuse, not this material-topology restriction. Those capabilities are
retained, not rebuilt here.

The adaptive candidate has conforming coordinate-plane refinements and two
independent distance checks. It permits four candidate rounds, the initial
round and three refinements. Its final material check still refuses every
nonempty root other than the supported block family with
`field.print-detail-topology` at `MaterialTopology`. Earlier unresolved cells
keep their own refusal and stage.

## Material-topology proof

Parked as research. A closed candidate with complete two-sided distance
coverage does not establish that its material has the same components,
cavities or handles as the definition. The retained guard must not be removed
on that evidence alone.

An exact counterexample uses the control block
`[-0.5, 3.5] x [-0.5, 0.5] x [-0.5, 0.5]` and a spherical void centred at
`(0.13, 0.11, 0.44)` mm, radius `0.03125` mm. The void is entirely inside the
block. Its gap to the top face is `0.02875` mm. Every point on the void boundary
is within `0.5 - 0.44 + 0.03125 = 0.09125` mm of that face. Filling the void
keeps the complete outer boundary and has two-sided boundary distance at most
`0.09125 < 0.125` mm. Yet it removes one cavity and changes the number of
boundary components from two to one. Even an exact distance certificate at
0.125 mm cannot authorize that change.

The exact point-query control checks that the near-surface centre is material
in the filled block and void in the hollow one, and that both z-axis sphere
endpoints are real cavity boundary points. It does not export that example.
The adaptive small-void control below keeps the earlier centre
`(0.13, 0.11, 0.09)` mm. It is a separate refinement refusal. The disjoint-block
control exercises the final guard after both distance checks, rather than
calling every early refusal a material-topology result.

Accepted cells retain the conservative global distance bound, not a tighter
bound inferred from their widths. For boundary points in one cell it bounds
distance to the complete mesh, not necessarily that cell's own facets.
Meshing and interval certification remain serial. The adaptive worker argument
controls final facet representation checks, which a material-topology refusal
does not reach. Worker equality here does not establish meshing speed-up.

Resuming requires a certificate that establishes source-to-candidate material
correspondence inside each mixed cell and on shared faces. It must preserve
interior components, cavities and handles; account for exact contacts; and
compose over conforming unequal cells. Independent closed-mesh validation
remains necessary. Neither a genus count alone nor a component-size filter
provides this certificate. Keep existing outputs unchanged and require complete
failure and geometry equality at 1, 2 and 4 workers.

## Sheet printing

Parked pending the material proof. The supported sheet definition is the band
`|g| <= c`, not a millimetre wall thickness. The control has period 4 mm and
constant level `c = 0.25`. Its two boundary sheets require separate complete
coverage and material correspondence through clipping and joins.

The tests exercise both 0.125 mm and 0.01 mm public print requests. These take
different existing meshing paths, so a refusal on one does not qualify the
other. The adaptive control uses 0.125 mm, 100,000 cells and 100,000 vertices.
Exact membership tests still cover both sides and the linear grades. No sheet
STL is admitted by this change.

Resuming requires printable native-field sheet output with both sheets intact,
complete two-sided error bounds, material orientation and preserved small voids.
Sheet queries, a network export or a B-rep envelope do not close that exit.

## P6 at recorded controls

The full public recipe is the 14-step `parts/P6-lattice-bracket.json`. Its
12 x 8 x 8 mm bracket has 3 mm legs, a 1 mm shell, 2 mm cubic infill cells with
0.3 mm beam radius, a 0.25 mm overlap at the cavity walls and a radius-0.5 mm
through hole. No feature or step is removed for these controls.

The bounded diagnostic uses the default print settings, 0.01 mm accuracy,
2,000,000 cells and 1,000,000 vertices, with the existing 1,000,000-unit
scoreboard allowance. It is not qualification at the larger verification
allowance. The explicit 0.02 mm request is a separate diagnostic choice with
the same limits.
It is not an automatic retry, not a request estimated to fit and not a 0.01 mm
export. The earlier 256-cell experiment remains a small-cap control only.

### Why the current uniform route cannot fit

This is an exact lower bound from the recipe and the current grid rule, not a
measured traversal total or a prediction of the adaptive route. The bottom
face contains the full rectangle `0 <= x <= 5`, `0 <= y <= 8`, `z = 0` mm.
The cavity starts at `z = 1`, and the bore has `x >= 5.5`, so neither removes
this 40 mm² patch. Strict interior/exterior classification cannot discard a
cell touching its boundary.

The uniform route chooses a dyadic grid edge `h` with `4h <= tolerance`.
At 0.01 mm, `h = 1/512` mm; at 0.02 mm, `h = 1/256` mm. The grid aligns to
the integer patch boundaries, so even one layer along this single patch needs:

| Explicit accuracy | Grid edge | Boundary-cell lower bound | Cell cap |
| --- | --- | ---: | ---: |
| 0.01 mm | 1/512 mm | `5 * 512 * 8 * 512 = 10,485,760` | 2,000,000 |
| 0.02 mm | 1/256 mm | `5 * 256 * 8 * 256 = 2,621,440` | 2,000,000 |

Parents, other faces, infill and later refinement add work. These counts are
therefore lower bounds only, already above the unchanged cap at both settings.
Doubling the accuracy value does not establish a fitting choice. This identifies
the uniform boundary subdivision as a capacity cause before any claim about
the lattice's material topology. Increasing a cap is not the implementation
of adaptive cells or their missing proof.

Resuming requires the complete native-field recipe to export within the chosen
recorded controls, with preserved material topology and a complete file bound.
P6 also needs its separate mass-equivalent STEP stand-in, including volume,
centre of gravity and inertia. That acceptance is not demonstrated by the print
controls here. No stand-in or two-kernel property agreement is claimed.

## Calibration and coarser choices

Parked pending qualified output and M26's measurement window. There are zero
throughput calibration samples and zero fitted throughput models in this
increment. The adaptive failure reports actual Decider work and charged request
work, plus actual cumulative visited cells, the current candidate's vertices
and facets, its zero-based round and unresolved cells observed so far. These
partial counts do not predict total demand or peak memory. Its explicit
alternative doubles the requested tolerance while retaining
both capacity limits. The label remains `UNVALIDATED`; its estimate states that
fit, cells, triangles, peak memory and calibrated time are unavailable.

`unavailable` is the estimate status, not a fabricated runtime refusal code.
The geometric or capacity refusal remains the original typed cause reported
below. A calibrated throughput cannot resolve missing material topology.

Resuming requires bounded coarse traversal to estimate unresolved-cell demand,
identified refinement overhead, live-storage accounting for peak memory and
throughput samples tied to an implementation and measurement environment.
Separate actual counts, estimates, lower bounds and certificates. Measure
prediction error on held-out requests before claiming a coarser request is
estimated to fit. Retain the original allowance while estimating. Timing work
belongs to M26 and no elapsed-time target is evaluated here.

## Fresh execution evidence

The six original adaptive controls passed freshly on the baseline above.
The disjoint-block guard passed on
`b8c18ce9993e40ac7b6202886f419254c3a9df22`. Field work below measures meshing
after plan compilation, not the complete CLI request. Raw evidence retained
privately.

| Field control | Observed result | Charged request work |
| --- | --- | ---: |
| 4 x 1 x 1 mm block, 0.125 mm, 1/2/4/default workers | Same 4,368 facets and 13,710 visited cells; independent closed-mesh check passed | 63,799 |
| Period-4 sheet, level 0.25, 0.125 mm | `Refinement`, `field.print-unresolved`, 436 unresolved homogeneous cells | 693,843 |
| Radius-0.03125 mm cavity at the earlier centre | `Refinement`, `field.print-unresolved`, 5 unresolved homogeneous cells | 97,015 |
| Two disjoint blocks, 0.125 mm, 1/2/4 workers | `MaterialTopology`, `field.print-detail-topology`; both distance checks completed | 74,587 each |
| P6 and V7c at the deliberately small 256-cell cap | `Preparation`, `field.print-budget`; not default-cap qualification | 1 each |

The real CLI also completed the following requests at one worker. Its print
work includes compilation. These are completed observations within longer
matrices, not a completed default-cap comparison at 1, 2 and 4 workers.

| Public CLI request | Observed refusal | Request work |
| --- | --- | ---: |
| Full P6, 0.01 mm, default print caps, 100M allowance | `field.print-budget`, uniform interval-refinement cell limit | 2,000,001 |
| Full P6, explicit 0.02 mm, same caps and allowance | `field.print-budget`, same limiting stage | 2,000,001 |
| Both P6 accuracies, same print caps, 1M allowance | `numbers.operation-budget` | 1,000,000 each |
| Sheet control, 0.125 mm and 0.01 mm | `field.print-detail-topology`, before meshing | 0 each |

The 100M public-CLI worker matrix did not finish its observation window. No
unfinished request is counted as a geometric refusal or a worker comparison.
The native adaptive full-budget P6 exit remains unestablished. The uniform
lower bounds above stand independently of these execution windows.

The field controls compare complete failure values, including attached
undecided evidence, request settings and work, at 1, 2 and 4 workers. The CLI
control starts the real `catallax-verify` command at each worker count, compares
complete responses except clocks, checks the on-disk manifest and verifies
that stale `print.stl`, `display.obj` and `answer.step` files are removed.

The routine suite includes bounded 1,000-unit P6 and sheet controls and the
ordinary 0.125 mm sheet refusal. The larger 1M diagnostics remain explicit
ignored tests and do not qualify P6. Reproduction from the repository root:

```sh
cargo test --release -p catallax-field --test adaptive_print -- --nocapture --test-threads=1
cargo test --release -p catallax-cli --test m19_print_verdict -- --nocapture
cargo fmt --all --check
git diff --check
```

To attempt the separate 1M diagnostics, add `--ignored` after `--` to the
respective test command. An incomplete attempt is not a passing test result.

The 477 saved replay inputs are not part of this focused run. The new counts
do not feed a geometric decision, change a work allowance or approve a saved
answer difference. An exact release-candidate replay remains the release
procedure's gate.
