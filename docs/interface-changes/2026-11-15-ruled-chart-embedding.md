# Rational ruled chart embedding

`RationalPatch::new` can prove a single quadratic-u, linear-v rational Bezier chart after the centre-tangent projection fails its root interval check.
The controls must be rational, with positive weights repeated across v.
The carrier, knots, native domain and controls stay unchanged.
The fallback proves the existing complete symmetric-Jacobian lower bound; samples do not supply this proof.
Other layouts retain the interval subdivision path.
Atlas cell admission keeps its existing interval criterion.
This preserves the chart sizes used by mass-property integration.

For either Cartesian projection, the u derivative numerator has quadratic Bernstein controls `2*w0*w1*(p1-p0)`, `w0*w2*(p2-p0)`, and `2*w1*w2*(p2-p1)`.
Their minimum and maximum, divided by the appropriate squared weight bounds, enclose the derivative over both v rows and the complete unit rectangle.
The v derivative is a positive weighted average of the three projected row differences.
These formulas also bound nonparallel rows caused by Cartesian rounding.

Keep the first centre-dual projection f and construct a second projection g from the lower row's control plane, normalized by its first v displacement.
The bounds give `a <= f_u <= A`, `|f_v| <= M`, `|g_u| <= N`, and `b <= g_v <= B`.
Require positive a and b, choose `k = 1 + M²/(a*b)`, and require `D = a*k*b - (M+k*N)²/4 > 0`.
The symmetric Jacobian of `[f,k*g]` has determinant at least D and trace at most `A+k*B`.
Therefore `D/(A+k*B)` is a positive global lower eigenvalue bound.
Strong monotonicity proves regularity and injectivity on the convex rectangle, with the existing diagonal-congruence rule for restrictions.

This is an exact support certificate only for the retained rational definition.
A chart built from numerical native controls retains the source's numerical category and tolerance; this proof adds no exact relation to the original analytic surface or an intersection curve.
The default exact cone/cylinder window gate is unchanged.

OCCT [splits finite native B-spline spans into Bezier patches](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/GeomConvert/GeomConvert_BSplineSurfaceToBezierSurface.cxx#L89).
The P04 SMLib port uses the same finite parametric supports through `SmCrvOnSurf`, recorded in [the port scope](../M25-INTERSECTION-PORT-SCOPE.md).
The numerical consumer keeps that representation.
This small rational certificate satisfies Catallax's additional global-chart requirement without asking the numerical source for an exact analytic-incidence proof.

The `numerical_hole_carriers` example accepts `--check-supports` after the worker count.
It retains each finite span domain, complete chart lower bound, projection and measured work alongside its numerical intersection report.
Run it with the same actual recipes and fresh output directories at 1, 2 and 4 workers.
The unit regressions cover wide cones, nonunit domains, affine transport, rounded nonparallel rows, unique finite inverse queries and collapsed-apex refusal.
