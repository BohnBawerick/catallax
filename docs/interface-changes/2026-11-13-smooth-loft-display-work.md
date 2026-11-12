# Smooth loft display work

The certified application display of three centered square sections, widths
2, 4 and 2 mm at z = 0, 3 and 6 mm in Smooth mode, now fits the unchanged
20,000,000-unit application allowance. Its exact volume is 4596/70 mm^3.

## Work attribution

The development `app_work_profile` observer reproduced the refusal on main
`33d8cf0729d774db39700a64227da37fbe13f178` at the default 0.01 mm chord tolerance.
Its `--exact-only` run completed in 34,234 units. The complete application
entry exhausted 20,000,000 with `application.work-limit` after the exact loft
had already built.

Required shared stations were inserted into each rational display-grid axis
with a linear scan. The observer charged 4,196,350 exact sign-filter calls,
16,785,400 units, to that scan before refusal. This was 83.93% of the entire
application allowance.

Display-grid insertion now uses bisection on the already sorted axis and the
same exact sign predicate. It retains exact ordering and duplicate removal,
including algebraic stations. Section-curve station insertion keeps its
existing implementation.

| Measurement | Main | With bisection |
| --- | ---: | ---: |
| Exact build | 34,234 | 34,234 |
| Complete application build | 20,000,000, refused | 7,435,243, built |
| Display work after exact build | 19,965,766, unfinished | 7,401,009, complete |
| Display-grid insertion | 16,785,400, unfinished | 467,472, complete |
| Display positions | unavailable | 16,388 |
| Display triangles | unavailable | 32,772 |
| Remaining application allowance | 0 | 12,564,757 |

The requested chord tolerance remains 0.01 mm, with the same complete-cell
chord bounds, embedding and rounding proofs, shared boundary checks, vertex
cap and optional controls. The exact loft, mass and STEP definitions are
unchanged. These are deterministic work measurements, not a timing pass.

## Regression checks

`natural_three_square_loft_builds_inside_the_default_allowance` exercises the
application entry with its default guard. It checks every face's chord bound,
physical picking keys, 6 faces, 12 edges and 8 vertices, volume 4596/70 and
centroid [0, 0, 3]. A separate public grid test covers unordered required
stations, an algebraic duplicate and both domain endpoints.

The application test also retains an exact checkpoint and requests its display
with a fresh decider under a separate 20M allowance. Its mesh must equal the
eager application's mesh, including all face and edge certificates.

The development observer also reports `display_obj_sha256`, binding the
returned display coordinates, triangles, face groups and chord bounds. Fresh
runs at 1, 2 and 4 kernel threads have identical complete observer records:

| Kernel threads | Work | Display OBJ SHA-256 |
| ---: | ---: | --- |
| 1 | 7,435,243 | `8846fee09dd0f323bb8799c7359e59a4eb13b84c2d57e3f87aa2854a233f0fe7` |
| 2 | 7,435,243 | `8846fee09dd0f323bb8799c7359e59a4eb13b84c2d57e3f87aa2854a233f0fe7` |
| 4 | 7,435,243 | `8846fee09dd0f323bb8799c7359e59a4eb13b84c2d57e3f87aa2854a233f0fe7` |

Raw evidence retained privately.
