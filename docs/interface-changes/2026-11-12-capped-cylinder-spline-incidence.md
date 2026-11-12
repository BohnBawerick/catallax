## 2026-11-12: circular spline incidence on capped cylinders

Pieces: production STEP body import and local source-deviation marks.

A measured spline incidence refusal can now produce an exact circular edge
when its complete support neighbourhood is one unchanged cylinder and two
perpendicular planar caps. The source interval and target finite arc need an
exact one-to-one parameter correspondence, complete nonoverlapping ring
coverage, unchanged shared vertices and a proved source feature margin. The
existing strict solid certificate still runs on the derivative.

Each changed physical edge adds `import.incidence-curved-repair`, a distance
mark. Its bound is an outward upper bound on the maximum chosen pointwise
parameter-correspondence displacement. That upper bound also bounds both
finite-curve Hausdorff directions. The detail JSON contains:

- `solid`, `entity`, `entity_kind` and `source_geometry_sha256`.
- `definition`, explicitly naming the parameter-correspondence metric.
- `deviation_lower_mm` and `deviation_upper_mm`, outward decimal bounds on
  that maximum displacement. The lower endpoint is not a Hausdorff bound.
- `approved_tolerance_mm`, the unchanged 0.01 mm import cap.
- `feature_size_squared_mm2`, an exact rational lower bound for the source
  neighbouring edge lengths and face widths, squared.

The attained displacement also respects a tenth of the source feature lower
bound. Zero movement adds no tolerance mark. Every new arithmetic operation
uses the existing shared work allowance. Larger shells, moved vertices,
physical generator edges and spline supports still refuse with a named
certificate step. Unproved bounds retain `incidence_gap` with its measured
source gap and repair reason. Exhaustion remains `resource_limit`.

No enum, proposal or repair-record version changes. Readers that display
local import deviations should retain this additional body mark and its
metric definition. It does not authorize later edits or a looser tolerance.
Existing exact admissions and the planar repair path keep their answers,
marks, work and bytes. This slice adds no curved user-override protocol.

The proof, limitations and frozen census comparison are in
[the R8 report](../milestones/R8-curved-incidence.md).
