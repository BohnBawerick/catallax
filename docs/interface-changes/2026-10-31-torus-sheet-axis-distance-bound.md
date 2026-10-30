## 2026-10-31: torus-sheet distance bounds reach the axis points

Pieces: 2

Change: a `NotExact` from `side`, `side_of` or `point_on_curve` on a spindle
`TorusSheet` now converts the residual bound `b` to `b/h`, with
`h^2 = minor^2 - major^2`, not `b/minor`. A horn sheet keeps `b/minor`. A sheet
with irrational radii whose `h` cannot be separated from zero refuses with
`numbers.divisor-not-separated`.

Why: `b/minor` holds only where the nearest tube point is retained. In the
triangle between the tube centre and the axis points `(0, +-h)` the nearest
retained point is an axis point. The centre of the `(1, 4)` sheet is `sqrt(15)`
from it, not `15/4`.

Migration: none. Bounds are larger, never smaller; proven sides, coincidences
and other surfaces are unchanged.
