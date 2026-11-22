# Independent sweep and loft shells

Existing recipes and comparison options keep their format. Catallax constructs
more sweep and loft bodies in the field from the saved path, sketch and loft
correspondence. The B-rep uses its existing construction route.

Fixed spline translations use the same tensor-product construction as
SMLib's `N_CreateTransSweepSrf`. Cardinal line and arc chains use successive
translations and meridional sectors with a carried frame. OpenCASCADE's
`BRepOffsetAPI_MakePipeShell` distinguishes a fixed trihedron from transported
frames and supplies transition modes at spine breaks. Catallax uses the
stationary cardinal turns defined by its recipe.
[SMLib NLib](https://docs.nvidia.com/smlib/manual/nlib/functions/index.html) ·
[OpenCASCADE pipe shell](https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___make_pipe_shell.html)

Lofts follow the saved section interpolation. The smooth case uses the
natural cubic equations from Carl de Boor's *A Practical Guide to Splines*
with the recipe's unit station intervals. Compatible
profile intervals come from the saved segment charts. OpenCASCADE's
`BRepOffsetAPI_ThruSections` also offers ruled and smoothed constructions.
Catallax retains its authored start segments and common rational breaks.
Parallel Bézier roofs use the existing SMLib-derived Coons construction for
their planar caps. Nonparallel affine polygon pairs use a symmetric Jacobian
test over a convex source box to prove the loft is embedded.
[OpenCASCADE through sections](https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___thru_sections.html) ·
[Catallax NURBS toolkit](../NURBS-TOOLKIT.md)

`Plan::result_patch_shell` exposes a sole closed patch shell to a caller
that needs its unsigned boundary range. `catallax-compare` uses the range
for B-rep display vertices of terminal sweep and loft shells. It adds each
vertex's incident display bound before applying the existing comparison
accuracy. Reverse named-face and normal checks for these shells remain
unresolved. There is no change to the 20M work target or 100M allowance.

## Exact recipe inventory

Fixed spline caps orient each ring according to nesting and use an outer ring
for the cap normal. The coaxial loft builds compatible circular sections at
their saved angular starts, including loops authored as quarter arcs.

Both kernels built 19 of the 27 saved or reconstructed sweep and loft
prefixes within 20M work. This count covers construction, including an
empty coplanar loft. The results and work matched at 1, 2 and 4 threads.

At the unchanged 100M allowance, the comparison completed 18 of those
19 inputs. It found no forward discrepancies across 441,885 B-rep
display vertices. Five cases reported agreement, including the empty
loft. The other 13 returned unresolved because the patch shells lack
reverse named-face samples and normals. The largest sampled upper bound
was below 0.009968 mm against the existing 0.01 mm comparison accuracy.

**Curved housing loft.** Display comparison carried to 3.5. The field
build uses 110,867 units and the B-rep build uses 1,371,690. The existing
0.005 mm B-rep mesher exhausts 100M total units before the comparator
can obtain a display mesh or run a field boundary query.

| Prefix | First stopped B-rep step | Exact runner result |
|---|---|---|
| P8 thread | `thread_path` | `sweep.requires-numerical-evaluation` |
| P8 spring | `thread_path` | `sweep.requires-numerical-evaluation` |
| P9 outer loft | `outer` | `loft.admissibility-unproven` |
| P9 inner loft | `outer` | `loft.admissibility-unproven` |
| Schema fixture tube | `holes` | `intersect.trace-native-window` |
| Schema fixture column | `holes` | `intersect.trace-native-window` |
| Crossed linear spline | `body` | `sweep.spline-admissibility-unproven` |
| Half turn with no saved axis | `body` | `sweep.missing-corner-axis` |

The P8 rows indicate routing to a fitted evaluator. The fitted P8 part has
its own recorded refusal at the threaded union. The standalone fitted helix,
guide and twist comparisons are the next part of 4.2.
[Fitted sweep evidence](../RECIPE.md#work-recorded-for-phase-5)
