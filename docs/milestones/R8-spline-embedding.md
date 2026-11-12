# Exact spline embedding after the interval filter

The source cubic quarter-circle approximation can have a strictly positive endpoint-chord derivative while its power-basis interval enclosure includes zero. Import used to stop at `import.spline-embedding-unproven` in that case, before measuring source incidence.

The existing interval check remains the first path. Its successful answers and work counts stay unchanged. On that specific unproved result, an exact rational Bernstein certificate now tries the complete active source interval:

- Extract homogeneous Bezier spans by charged rational blossoming, including nonclamped source knots.
- Require positive weights and exact C0 agreement between adjacent spans.
- Use one common direction, from the actual first endpoint to the actual last endpoint.
- Form the Bernstein coefficients of the projected rational derivative numerator `S' W - S W'`. Strictly positive coefficients place its hodograph in an open half-space. Exact midpoint subdivision can establish positivity when the initial hull crosses zero.

One common projection increases strictly across every span and join. It therefore proves a regular, one-to-one curve, and also proves embedding of any retained finite subinterval. Different directions on individual spans would not give that global conclusion and are not used.

Every new rational operation, comparison and subdivision is charged before execution. The certificate stops at 4096 visited cells, depth 24 or the caller's existing work allowance. No tolerance or work allowance changes. A zero endpoint derivative, discontinuous join, closed projection or curve without a proved common monotone projection remains refused. Failure keeps the existing typed embedding refusal; operation exhaustion remains `resource_limit`.

Unclamped source controls can start and end at the same stored point while their actual active endpoints differ. The fallback uses those active endpoints. The line-collapse shortcut also requires a nonzero source-control chord, so a certified curved arc cannot become a zero-direction line. Its regression checks the returned spline carrier and exact endpoint parameters.

The production STEP test now reaches the next measured refusal for a cubic quarter-circle approximation: its exact polynomial carrier is not exactly on its authored cylinder. The measured maximum gap is between 0.000272529946630 and 0.000272530425428 mm. Curved repair is a separate step. A genuine cubic self-intersection, with two distinct interior parameters mapping exactly to the origin, still refuses at embedding.

Focused controls cover the old fast path's identical work count, rational weighted span extraction against the independent kernel, discontinuous spans, exact subdivision, zero and partial work budgets, and larger operand charging. Raw evidence retained privately.

## Frozen census comparison

The same 728 model IDs were imported with the base production implementation and the candidate, using the existing 20-million work allowance. The control compiles the base import entry and spline checks from `ac63201c`; all other compiled production files equal that base. The candidate producer is `5fa525d7e2cbc07ecc835efe6d62c8445162a5b1`. Binary and source hashes are in [the manifest](R8-spline-embedding.json).

| First result | Base | Candidate |
| --- | ---: | ---: |
| `ADVANCED_FACE:boundary_certification` | 3 | 3 |
| `ADVANCED_FACE:support_admission` | 108 | 108 |
| `CONICAL_SURFACE` | 93 | 93 |
| `EDGE_CURVE` | 62 | 62 |
| `EDGE_CURVE:curve_embedding` | 12 | 6 |
| `ELLIPSE` | 12 | 12 |
| `admitted` | 1 | 1 |
| `incidence_gap:edge_surface_incidence` | 133 | 133 |
| `incidence_gap:spline_curve_surface_incidence` | 33 | 39 |
| `incidence_gap:spline_endpoint_incidence` | 125 | 125 |
| `incidence_gap:vertex_curve_incidence` | 119 | 119 |
| `resource_limit` | 27 | 27 |

Newly admitted files: 0. The previously admitted file keeps identical answers, local marks, work and raw result bytes. 6 files reach a different first refusal. The remaining causes are recorded per model in [the CSV](R8-spline-embedding.csv). No observation is incomplete.

The 14-model subset covers all previous embedding refusals plus the admitted and work-limited controls. Its results, marks, work and output bytes match the full census at 1, 2 and 4 process and kernel workers. The two synthetic STEP controls also match at 1, 2 and 4 workers. All 213 import tests and strict import Clippy pass. Saved-answer replay and STEP regression results are recorded in the PR checks.
