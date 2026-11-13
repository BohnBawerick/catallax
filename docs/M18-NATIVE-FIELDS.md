# M18 native field construction

M18 remains open. This increment adds exact field construction for reduced sweep
and loft families, circular smooth-blend point queries and an on-demand table of
independent recipe-prefix builds. It does not replace the authored P5 verdict.
The radius-2 overflow and its geometry are already recorded on main in
[the P5 note](M18-P5-VERDICT.md), delivered by
https://github.com/BohnBawerick/catallax/pull/542.

## Construction coverage

| Definition | Field construction | Remaining limit |
| --- | --- | --- |
| Empty path | Empty solid | No sketch relocation |
| Rational straight path | Exact endpoint displacement, either orientation | A displacement tangent to the sketch has zero thickness |
| Fixed line chain | Regularized union of translated profile prisms | Curved spans remain separate work |
| Rotation-minimising straight chain | Collinear forward spans | A tangent change needs the stationary rotation region |
| One rotation-minimising cardinal arc | Existing field revolution definition, with path-qualified side names | Meridional profile and the revolution domain; angles 90, 180 or 270 degrees |
| Translated rational-chart loft | Exact prism; ruled stations are ordered on a line | Smooth stations must be uniformly spaced translations |
| Two coaxial full-circle sections | Exact cone or cylinder, including different radii | Positive axial separation and corresponding circle charts |
| Circular smooth blend | Exact signed stock distances, exact membership and terminal boundary origins | No certified implicit display/print patch or downstream consumer |

The field reads recipe definitions and constructs its own membership predicates.
These paths do not call the B-rep evaluator or import its generated boundary.
Ordinary extrusions keep their existing path and arithmetic order. No tolerance
or work allowance changes.

Sweep `M(0)` is identity. A path starting at `(100, 200, 300)` does not move a
sketch to that point. The first line's endpoint displacement is applied directly,
without normalizing it into a length. Intermediate exposed span caps use the
`envelope` path role; only the first and last caps use `start` and `end`.

## Circular blend values and names

The field implements RECIPE section 4.17 for `union`, `intersect` and `subtract`.
It computes exact signed Euclidean distances for placed blocks, cylinders,
cones, spheres and ring tori. Intersect negates both signed distances before
the circular union formula and negates the result. Subtract first negates the
tool distance. Nonpositive radii refuse with `smooth.non-positive-radius`.
Composite operands refuse with `smooth.distance-unsupported`; a Boolean
minimum or maximum residual never substitutes for true distance.

The returned `Value::Residual` and fast value are `f/2`, where `f` is the
recipe formula. In the active circular region the derivatives with respect
to the two signed distances have Euclidean norm one and sum at most
`sqrt(2)`. Each signed distance is 1-Lipschitz in model space. Elsewhere the
formula reduces to a minimum or maximum of those distances. Thus `f/2` is
at most 1-Lipschitz globally, and its absolute value is a certified lower
bound on distance to its zero set. It is not signed distance. The returned
distance upper bound remains infinity away from a certified boundary point.

`Plan::smooth_origins_checked` returns `smooth` exactly where both circular
terms are active on the boundary. At the unchanged boundary it returns only
the contributing stock origins. It returns no names away from the boundary.
The existing analytic `faces`/`select` registry has no carrier for the new
implicit patch; zero registered smooth faces is an explicit limit. Normals
remain ambiguous. Display and print preparation refuse with
`smooth.boundary-unsupported`, reporting zero admitted patches. The smooth
step must be the recipe result; other consumers refuse `smooth.not-terminal`.

The right-angle test uses `k=5 mm`, signed distances `a=2 mm`, `b=1 mm` and
the exact 3-4-5 circle, so the new union boundary passes through `(2,1,0)`.
Sign-reflected points test intersect and subtract independently.

## Independent step table

```sh
catallax-verify --threads 1 steps parts/P1-bracket.json --max-work 20000000
```

The command writes `catallax-independent-steps/1` JSON. It attempts every
body-valued recipe prefix independently in both kernels, in fresh threads and
with the same supplied allowance per kernel and prefix. Each row records work,
the limit, result face candidates and empty status, or the original typed
refusal and its first blocked step. A later row blocked by an earlier operation
does not claim its own operation was reached. Sketches, paths and scans are
reference rows, validated by the recipe and subsequent independent body builds.

An envelope-only B-rep stand-in reports `compare.envelope-only`. The table
proves construction coverage, not display, measurement, exchange or geometric
agreement. Existing sampled `source: compare` remains the request for a
two-way boundary check. Its sweep and loft coverage now comes from the actual
compilers instead of an unconditional unsupported-operation list.

Run the public evidence collector with an external proof directory:

```sh
python3 tools/checks/m18_field.py --binary target/release/catallax-verify --out "$PROOF"
```

It drives native field display fixtures at 1, 2 and 4 workers, compares complete
manifests and display bytes, and captures every P1-P5 prefix at 1, 2 and 4
workers and kernel threads. Raw evidence retained privately.

## Work still required

- General spline-path translations, fixed circular translations and joined
  rotation-minimising arc/corner motions. Unsupported curved spans report
  `field.sweep-motion-unsupported` with their index and total span count;
  unsupported corners report `field.sweep-corner-unsupported` with the join.
- Non-translated rational lofts, including the three-circle natural lofts.
  `field.loft-correspondence-unsupported` reports the first section outside
  the checked correspondence and the total section count. No general loft
  patch is admitted by this increment. Spline-profile validation retains its
  original field refusal.
- Smooth blends of composite bodies need certified true distance to their
  trimmed result, followed by implicit boundary coverage for display, printing,
  normals, selection handles and later features. There are zero admitted
  composite-distance definitions and zero certified smooth boundary patches.
- P5 still needs the four mouth-contact envelope regions described in its
  existing verdict. This increment does not revisit or widen that admission.

These limits keep M18 open. They are construction/proof work, separate from
the field speed work assigned to M26.
