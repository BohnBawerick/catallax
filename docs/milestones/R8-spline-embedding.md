# Exact spline embedding after the interval filter

The source cubic quarter-circle approximation can have a strictly positive endpoint-chord derivative while its power-basis interval enclosure includes zero. Import used to stop at `import.spline-embedding-unproven` in that case, before measuring source incidence.

The existing interval check remains the first path. Its successful answers and work counts stay unchanged. On that specific unproved result, an exact rational Bernstein certificate now tries the complete active source interval:

- Extract homogeneous Bezier spans by charged rational blossoming, including nonclamped source knots.
- Require positive weights and exact C0 agreement between adjacent spans.
- Use one common direction, from the actual first endpoint to the actual last endpoint.
- Form the Bernstein coefficients of the projected rational derivative numerator `S' W - S W'`. Strictly positive coefficients place its hodograph in an open half-space. Exact midpoint subdivision can establish positivity when the initial hull crosses zero.

One common projection increases strictly across every span and join. It therefore proves a regular, one-to-one curve, and also proves embedding of any retained finite subinterval. Different directions on individual spans would not give that global conclusion and are not used.

Every new rational operation, comparison and subdivision is charged before execution. The certificate stops at 4096 visited cells, depth 24 or the caller's existing work allowance. No tolerance or work allowance changes. A zero endpoint derivative, discontinuous join, closed projection or curve without a proved common monotone projection remains refused. Failure keeps the existing typed embedding refusal; operation exhaustion remains `resource_limit`.

The production STEP test now reaches the next measured refusal for a cubic quarter-circle approximation: its exact polynomial carrier is not exactly on its authored cylinder. The measured maximum gap is between 0.000272529946630 and 0.000272530425428 mm. Curved repair is a separate step. A genuine cubic self-intersection, with two distinct interior parameters mapping exactly to the origin, still refuses at embedding.

Focused controls cover the old fast path's identical work count, rational weighted span extraction against the independent kernel, discontinuous spans, exact subdivision, zero and partial work budgets, and larger operand charging. Raw evidence retained privately.
