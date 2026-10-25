# Piece 1: numbers and decisions (`catallax-numbers`)

Status: complete for blind test writing and blind building.

This is the only input the test writer and the creator get. Read `docs/DESIGN.md` and `docs/ARCHITECTURE.md` for the rules around it; nothing else is needed.

## Approved algebraic-definition extension

The shared foundation implements the 2026-10-22 architecture decision below. The rest of this document describes the retained rational baseline unless this section overrides it.

`Real` stays the existing lazy algebraic expression type. Do not add transcendental nodes or a second scalar type. A stored geometry definition must have proven-valid domains for division and square roots. Constructing an expression remains infallible; validating a geometry definition is a separate decision.

Add `RealPoly`, a polynomial in a fixed number of variables with `Real` coefficients. Use `Monomial`'s existing order and a deterministic expression order when combining terms. Arithmetic, partial derivatives, substitution and exact evaluation build expressions without making decisions. Expose `nvars`, ordered `terms`, `eval_real`, `partial`, arithmetic and conversion from `MPoly`. `as_rational() -> Option<MPoly>` succeeds only when every coefficient's existing `Real::as_rational` succeeds. It must not round or invoke elimination. The rational conversion feeds existing solvers unchanged.

`rational() -> Option<&MPoly>` borrows the rational backend without copying coefficients; `into_rational() -> Option<MPoly>` consumes it. `total_degree`, `is_zero`, `equals`, `resultant` and `real_roots` accept the caller's mutable decider; `eval_ball` accepts its enclosure context. `Decider::rational_value` can certify a rational value for a correlated algebraic expression. Its `None` result is inconclusive and must not be used as an irrationality or domain certificate.

`RealPoly::gcd` takes two univariate polynomials and the caller's mutable
decider. It returns their monic greatest common divisor in the existing exact
coefficient field, with `gcd(0, 0) = 0`. It validates all coefficient domains
before a constant shortcut and certifies algebraically zero leading terms.
The sparse Euclidean remainder path does not eliminate algebraic coefficients
into a larger rational polynomial; numeric refusals and work limits propagate.

A syntactic term bound is not the mathematical degree. Algebraically zero coefficients may remain until a decider proves zero. Degree selection, leading-coefficient division, polynomial equality and zero-polynomial tests must therefore be certified operations; no caller may use list length or expression identity as those proofs. `RealPoly` has no mathematical `Eq`, `Ord` or `Hash`. Enclosure evaluation includes coefficient error as well as input error, using the caller's precision and decision policy.

For non-rational coefficients, solvers must either operate in a certified common algebraic extension or lower the expression graph to rational polynomial equations with auxiliary variables and certified root/domain selectors. Lowering retains shared leaves and their correlations, positive square-root choices and nonzero divisors. Elimination gives candidate roots; it must not admit other conjugates or denominator zeros. Certify and order candidates against the original equations and component selectors. Rational `UPoly`, `MPoly`, `PolySystem` and `IsolatedRoot` keep their current contracts. Lifted systems need their own certificates when singular; the existing nonsingular `IsolatedRoot` constructor cannot certify them by assertion.

Add a shared budgeted `sin_cos_degrees(angle: &Rational, dec: &mut Decider) -> Result<[Real; 2], Refusal>`, returning sine then cosine. Reduce the exact angle modulo 360, select the intended roots of unity using exact sector/root certificates, and retain those certified algebraic values. Standard rational/quadratic angles use exact fast paths. No floating-point trigonometry or rationalised answer is allowed. Arbitrary rational degree angles may have high algebraic degree; failure to construct their certificate is a temporary coverage/resource failure, never impossible geometry. General `sin(Real)` and `cos(Real)` are not introduced. Continuous trigonometric chart evaluation remains in `Ball`.

Use the existing deterministic work policy for coefficient validation, root selection and elimination. A failed proof keeps its `numbers.*` reason. Cache reuse must preserve outcomes and deterministic work charges across cold/warm and serial/parallel runs. The alternative of eagerly replacing every expression by one isolated number simplifies arithmetic boundaries but pays elimination before it is needed and loses shared construction identities.

## Geometric proof records

`ModellingTolerance::default()` records `1/1_000_000` mm before proof work starts. Its value has no mutation API. `check` accepts only nonnegative millimetre bounds strictly below that value. Equality refuses with `numbers.geometric-gap-exceeded`; negative bounds and other units refuse with `numbers.geometric-gap-unproven`. This setting is not an export accuracy or an entity-local approximation cap.

`EvidenceBound` holds an exact rational upper bound and an `EvidenceUnit`. `enclose_magnitude` uses the caller's decider to enclose a shared `Real`, retaining the caller's declared units without converting a residual into geometry. `add` adds nonnegative bounds with matching units. The geometric caller must prove that addition is the appropriate transfer rule.

`CertificateIdentity<K>` names the immutable recipe revision, representation, exact target and complete `CertificateScope<K>`. Scopes distinguish a whole-result cover, complete named entities, finite parameter boxes and individual queries. `GeometricEvidence<K>` retains that binding, an optional geometric bound, `GeometricProof<K>`, ordered dependencies and `ProofContext<K>` with producing piece, step and stage. A missing metric bound is not zero displacement. A conversion additionally retains `DecisionEvidence<K>` with its original code, residual enclosure, units, affected entities and selected branch. An `Approximation` proof names the retained exact definition and its correspondence theorem.

`GapCertificate<K>` is a scoped metric claim under the strict 1 nm decision or legacy-output rule. `GeometricCertificate<K>` instead lists every `RegionCertificate<K>` in a complete exact/exception cover, converted bounded choices and a whole-result evidence root. An exact region requires zero-displacement evidence. Every changed region, including an exact support with a moved trim, needs a `LocalException<K>`. An optional uniform metric summary needs its own complete oriented proof; it is not another sub-1-nm gate or a maximum inferred from face bounds.

`LocalToleranceSettings::default()` records a 100 nm kernel-fit target, inclusive 1 micrometre kernel cap and inclusive 10 micrometre import cap. `LocalException` records the attained enclosure, source, applicable limit and neighbouring edge/face size lower bounds. Both tenth-of-feature limits also apply. `NoNeighbour` explicitly claims absence; the geometric verifier must check adjacency rather than use it for an unknown size. A fit above target but within every hard limit passes arithmetic checking and reports `kernel_target_met=false`. An import's declared source uncertainty is separate from its file-geometry deviation; a looser limit requires a file-specific `ImportToleranceAcceptance` whose authorisation the consumer verifies. Decision-only exceptions retain the strict 1 nm rule and cannot be relabelled as fits.

`EvidenceBook<K>` is append-only and ordered by stable definition keys. Equal shared records deduplicate; conflicting keys refuse. `check_gap_structure` checks a gap's binding, units, strict cap, domain and reachable dependency graph. `check_structure` checks a whole-result certificate against the caller's actual region inventory, including every exception's source and feature-size caps, converted decisions, identities and reachable proofs. Missing dependencies, cycles, stale domains, unlisted changed regions and unproved metric summaries refuse. Neither method verifies geometric theorems, adjacency, topology, coverage or material orientation. Public records and structural success do not authorise publication; the geometric consumer verifies those obligations against immutable definitions. No mark-log heuristic supplies that proof.

`DecisionBook<K>` records a `RecordedChoice<K>` containing the chosen branch and supporting evidence key. `decide_once` calls its closure only for an unsettled key and checks cancellation before lookup and after a new decision. Failed or interrupted decisions record no choice. `record` accepts an identical repeated choice and refuses a conflicting branch or evidence key. Consumers must validate imported evidence before trusting it. Keys include immutable definition revisions, operands, scope and branch selectors, never pointer addresses or worker completion order. Iteration is in canonical key order.

Record and graph operations use the existing decider cancellation checkpoints and request work accounting. There is no second allowance, and no numeric decision path or legacy `Answer::is_exact` behavior changes merely because these types exist.

## What this piece is for

Every other piece of the kernel asks this piece two kinds of question.

1. "Give me exact numbers." Rationals of any size, exact real algebraic numbers, and enclosures (intervals that are proven to contain a true value) at any precision, including enclosures of pi, square roots, sines and logarithms for measurement.
2. "Decide this." What is the sign of this exact quantity? Which of these two exact numbers is larger? Where are the real roots of this polynomial? Is there exactly one solution of this polynomial system in this box?

Every decision ends in one of four outcomes: proven, proven coincident (exactly zero), not exact with a proven bound, or a named refusal. A decision never returns a sign it has not proven.

The piece also holds the shared vocabulary for time budgets, decision statistics, recorded construction facts and marks (the record of every "not exact" answer).

It has no geometry. Points, surfaces and curves belong to piece 2.

## Outcomes, and the one rule that must never break

`Outcome<T>` has four variants:

- `Proven(T)`: the strict answer `T` is certified. For a sign, the quantity is proven non-zero and `T` says which side of zero.
- `ProvenCoincident(Coincidence)`: the quantity is proven exactly zero. `Coincidence::proof` names the proof: `ExactArithmetic` (exact rational evaluation gave zero), `Algebraic` (an exact algebraic method showed it), `SeparationBound` (an enclosure narrower than a proven root-separation bound), or `ConstructionFact` (a recorded fact; this piece never produces that one itself, piece 2 does).
- `NotExact(NotExact)`: the budget ran out. `NotExact::bound` is a non-negative rational `b` with a proof that the true quantity lies in `[-b, b]`. The sign is unknown.
- `Refused(Refusal)`: the decision has no answer. `Refusal::code` is a stable code from the list below.

The rule: **a `Proven` or `ProvenCoincident` outcome is never wrong.** A true zero must never come back as `Proven`. A true non-zero must never come back as `ProvenCoincident`. This rule is checked by every test in this piece, and a single violation stops all work on the piece (DESIGN.md, R2).

The second rule: a `NotExact` bound is always a true bound. If the outcome says `[-b, b]`, the true value is in it.

## Numbers

### Integer

Arbitrary-precision signed integer. Exact. Values of thousands of digits must work.

### Rational

Arbitrary-precision exact rational in normal form: the denominator is positive and shares no factor with the numerator. Zero is `0/1`. Two equal rationals are `==`, `Ord`-equal and hash the same.

Parsing (`Rational::parse`, also `FromStr`) is exact and accepts exactly this grammar, with no surrounding whitespace:

```
rational  = sign? digits ( "." digits )? ( ("e" | "E") sign? digits )?
          | sign? digits "/" digits
sign      = "+" | "-"
digits    = one or more of 0-9
```

- `"1.1"` is 11/10. `"0.1"` is 1/10. `"1e-7"` is 1/10000000. `"2.5E3"` is 2500. `"-6/8"` is -3/4. `"-0"` is 0. `"007"` is 7.
- Rejected with `ParseNumberError`: `""`, `"."`, `"1."`, `".5"`, `"1.1.1"`, `" 1"`, `"1 "`, `"1_000"`, `"0x10"`, `"inf"`, `"nan"`, `"1e"`, `"1/-2"`, `"1/0"` (kind `ZeroDenominator`), and any exponent whose absolute value exceeds 100,000 (kind `ExponentTooLarge`).

Error kinds: `Empty` only for the empty string; `ZeroDenominator` only for a string that matches the fraction form and whose denominator digits are all zero (`"1/0"`, `"0/000"`); `ExponentTooLarge` only for a string that matches the decimal form and whose exponent's magnitude exceeds 100,000; `InvalidSyntax` for everything else (including `"   "` and `"1/0/2"`).

`Integer::parse` accepts `sign? digits` only, and returns only `Empty` or `InvalidSyntax`.

Formatting:

- `Display` for `Rational` prints `n` for integers and `n/d` otherwise, in lowest terms, with a leading `-` for negatives.
- `0^0` is 1 for both `Integer::pow` and `Rational::pow`.
- `to_decimal_exact` returns the exact decimal string when the denominator has no prime factors other than 2 and 5, and `None` otherwise.
- `to_decimal_bounds(significant_digits)` returns `(lo, hi)`, two decimal strings with at most that many significant digits, such that `lo <= value <= hi` and both are as tight as that digit count allows (each is the nearest such decimal on its side). When the value is exactly representable, `lo == hi`. `significant_digits` of 0 is treated as 1. Each string is written from its digits `d1 d2 ... dk` (trailing zero digits dropped, `k >= 1`) and its decimal exponent `E`, the value being `d1.d2...dk * 10^E`:
  - zero is `"0"`;
  - when `-6 <= E <= 20`, plain notation: the digits with the decimal point placed and zeros padded as needed, no decimal point when the value is an integer (`1.2 * 10^4` is `"12000"`, `1.2 * 10^-6` is `"0.0000012"`, 100 is `"100"`);
  - otherwise scientific: `d1`, then `.` and the remaining digits if `k > 1`, then `e`, a sign `+` or `-`, and the exponent without leading zeros (`"1e+30"`, `"2.5e-7"`);
  - a negative value has a leading `-`.

Conversions:

- `from_f64_exact` gives the exact value of a finite double (`0.1_f64` becomes 3602879701896397/36028797018963968), `None` for NaN and infinities.
- `to_f64_nearest` gives the nearest double, rounding ties to even, with overflow to infinity. It is for display and caches only.

Division of a `Rational` by zero through the `/` operator panics; `checked_div` and `recip` return `None` instead. `pow` with a negative exponent of zero returns `None`.

### Sign

`Sign` (`Negative`, `Zero`, `Positive`) is the exact sign of an exact type (`Integer`, `Rational`). Exact types never need a decision to know their sign.

### Interval

A closed interval of doubles, `[lo, hi]`, that is proven to contain the true value. All arithmetic rounds outward, so the result of every operation contains every value the exact operation could give on the inputs. Endpoints may be infinite; NaN never appears in an interval (operations that would produce one give the whole line instead).

- `Interval::from_rational(r)` contains `r`. When `r` is a double, the interval is that single point. Otherwise its endpoints are the two neighbouring doubles of `r`.
- `sign()` returns `Some(StrictSign)` only when zero is not in the interval.
- Division by an interval containing zero gives the whole line.
- `sqrt` of an interval entirely below zero returns `None`; otherwise the negative part is cut off first.

Edge rules:

- endpoints `-0.0` are stored as `+0.0`, so `point(-0.0) == point(0.0)`;
- `sign()` of an interval with an endpoint at zero is `None` (zero is in the closed interval);
- `mid()` of the whole line is 0; of `[lo, +inf]` with finite `lo` it is `f64::MAX`; of `[-inf, hi]` it is `-f64::MAX`;
- `[0, 0]` times any interval, bounded or not, is `[0, 0]`; `[0, 0] / [0, 0]` is the whole line.

Intervals are the fast filter. They must be cheap: no allocation.

### Ball

An arbitrary-precision enclosure: a midpoint and a radius, both exact binary fractions, representing `[mid - rad, mid + rad]`. Every ball carries a working precision in bits; arithmetic on two balls uses the larger of the two precisions. The radius of a result accounts for every rounding, so the true result is always inside.

Balls have proven enclosures of:

- `sqrt`, `square`, `recip`, `pow_u32`;
- `pi(precision)`;
- `exp`, `ln`, `sin`, `cos`, `atan`, `asin`, `acos`, `atan2`.

Each function returns `None` outside its domain, as the domain rules below define. At precision `p`, for an input of radius zero, the result radius should be within a small factor of `2^-p` times the magnitude of the result; that is a quality goal that is measured, not a correctness rule. Containment is the correctness rule.

Precision rules: a precision of 0 is treated as 1. `from_rational(r, p)` rounds `r` to a midpoint of `p` significant bits and puts the rounding error in the radius (radius 0 when `r` fits). `with_precision(p)` with `p` below the current precision rounds the midpoint to `p` bits and widens the radius to keep containment; with a larger `p` it only changes the working precision. `hull` uses the larger of the two precisions.

`from_ratio(numerator, denominator, p)` encloses an integer ratio directly, without first reducing it to a canonical `Rational`. It accepts either denominator sign and returns `None` only for a zero denominator. Common factors do not change the midpoint or radius. Rounding and precision rules are the same as `from_rational`.

`is_exact_zero()` is true only for the singleton ball `[0, 0]`; it is false for a nonzero-width ball that contains zero. It reads the binary fields without allocating or copying them. `mid`, `rad`, `lo` and `hi` still return canonical exact `Rational` values.

Domain rules: `sqrt` returns `None` only for a ball entirely below zero; otherwise it encloses the square root of the ball's non-negative part, like `Interval::sqrt`. `ln` needs `lo > 0`. `asin` and `acos` need `-1 <= lo` and `hi <= 1` (the ends may touch `±1`). `atan2` returns `None` when both balls contain zero. `recip` and `checked_div` return `None` when the divisor contains zero.

These functions are the kernel's own. They must not call the platform's floating-point library, so results are the same on every machine (DESIGN.md, R10).

Sine and cosine range reduction use a stored 4096-bit Machin enclosure of pi,
outward rounded to the required guard precision. Higher precision recomputes the
series. Loading and truncating the immutable dyadics charge no numeric work,
independently of first use or worker count. `Ball::pi` and `pi_checked` retain
their existing series construction and enclosure representation.

### Polynomials

- `UPoly`: a univariate polynomial with `Rational` coefficients, stored lowest degree first, with trailing zero coefficients removed (the zero polynomial has no coefficients and no degree). Operations: arithmetic, exact evaluation at a rational, evaluation on an `Interval`, a `Ball` and a `Real`, derivative, division with remainder, monic greatest common divisor, square-free factorisation (monic factors, each with its multiplicity, in increasing multiplicity; the zero polynomial and constants give an empty list), composition, and conversion to primitive integer coefficients with a positive leading coefficient (lowest degree first, empty for the zero polynomial). `UPoly::constant(0)` is the zero polynomial and `from_roots(&[])` is 1.
- `MPoly`: a polynomial in a fixed number of variables with `Rational` coefficients. `Monomial` holds one exponent per variable. Terms are kept in a fixed canonical order, which is also the `Ord` of `Monomial`: lower total degree first, and within one total degree the exponent vectors compared lexicographically, smaller first. Iteration order is therefore deterministic. Operations: arithmetic, exact evaluation at rationals, evaluation on intervals, balls and reals, partial derivatives, substitution of a polynomial for a variable, the resultant with respect to a variable (the Sylvester resultant; when one polynomial has degree 0 in that variable, the result is it raised to the other's degree; when either is zero, the result is zero), and conversion to `UPoly` when only one variable appears.

The resultant keeps the same number of variables; `var` simply no longer appears. Calling an `MPoly` operation with mismatched variable counts, a variable index of `nvars` or more, or an evaluation point of the wrong length is a programming error and may panic; so is `IsolatedRoot::coordinate` with an axis out of range.

### Algebraic

A real algebraic number. Either an exact rational, or the unique root of a square-free polynomial with integer coefficients inside a rational isolating interval `(lo, hi)` that contains no other root.

- `as_rational()` returns the exact value when the number is rational. Every rational root found by `Decider::real_roots` is returned in this exact form.
- `defining_polynomial()` is `q x - p` for a rational `p / q` in lowest terms. For an irrational number it is a square-free polynomial with primitive integer coefficients and a positive leading coefficient that has the number as a root; it need not be the minimal polynomial.
- `isolating_interval()` returns `(lo, hi)`; for a rational number both are the number itself. For an irrational number `lo < hi`, neither end is a root of the defining polynomial, and the number is the only root strictly between them.
- `refined(max_width)` returns the same number with an isolating interval no wider than `max_width`, keeping those properties. `max_width` must be positive; zero or negative is a programming error and may panic.
- `to_ball(precision)` encloses it.

`Algebraic` has no `PartialEq`: equality of two algebraic numbers is a decision. Turn them into `Real` values and use the `Decider`.

### Real

The exact real number type that the rest of the kernel computes with. It is an immutable, cheaply cloned expression over exact leaves:

- leaves: a `Rational`, an `Algebraic`, or one coordinate of an `IsolatedRoot` (a certified solution of a polynomial system, see below);
- operations: `+`, `-`, `*`, `/`, unary `-`, `sqrt`, `square`, `pow_u32`.

Building a `Real` never decides anything and never fails. Questions about its sign are asked through the `Decider`.
`Real` has no `PartialEq` or `PartialOrd`.

- `as_rational()` returns `Some` when the expression contains only rational leaves and the operations `+ - *`, `square`, `pow_u32` and `/` with divisors that are non-zero rationals; in that case it returns the exact value. It does no refinement.
- `approx()` returns an `Interval` that contains the value, computed once and cached. When the value is undefined (division by an interval containing zero, square root of an interval below zero) it may be the whole line.
- `to_f64()` returns a double near the value, for display only.

There are no transcendental operations on `Real`. `pi` exists only as a `Ball`. This keeps every zero question decidable in principle: every `Real` is algebraic.

Undefined values:

- A divisor proven exactly zero makes any decision that needs the quotient return `Refused("numbers.division-by-zero")`.
- A divisor that cannot be separated from zero within the budget makes the decision return `Refused("numbers.divisor-not-separated")`.
- A square root of a value proven negative gives `Refused("numbers.sqrt-of-negative")`.
- A square root of a value that cannot be proven non-negative within the budget gives `Refused("numbers.sqrt-not-separated")`. A square root of a value proven exactly zero is zero.
- Definedness is required of every subexpression, whatever its surroundings: `0 * (1 / 0)` refuses with `numbers.division-by-zero`, it does not simplify to zero. A decision is never `Proven`, `ProvenCoincident` or `NotExact` unless every square-root argument in the expression is proven non-negative and every divisor proven non-zero.
- When several subexpressions are undefined, the refusal reported is the one for the first undefined node met in a depth-first, left-operand-first, post-order walk of the expression.

### Polynomial systems and certified roots

- `PolySystem`: `n` equations in `n` variables (`MPoly` values that all have `n` variables, `n >= 1`).
- `RootBox`: a closed box given by rational lower and upper corners.
- `IsolatedRoot`: a system together with a box that is proven to contain exactly one solution of the system, at which the Jacobian is invertible. `coordinate(i)` returns that solution's `i`th coordinate as a `Real` leaf, so it can be combined with other numbers and its signs decided. `enclose(precision)` returns a ball for each coordinate.

Only the `Decider` creates an `IsolatedRoot`, by certifying a candidate box.
Finding candidate boxes is the caller's job (piece 4); this piece only certifies them.

## Decisions

### Decider

A `Decider` holds a `BudgetPolicy`, an optional wall-clock deadline and the running `DecisionStats`. Every decision creates a fresh per-decision budget from the policy.

Nonlinear root isolation may reuse complete exact certificates held by that decider in a bounded memo. A lookup compares the full polynomial, not only a hash. Reuse preserves the certificate's original deterministic work charge and still checks the operation deadline; it does not cache failed or partial isolation.

Clones of an irrational algebraic number share an immutable polynomial/interval snapshot. Its bounded enclosure memo keys the exact working precision and replays the original work cost; it never substitutes a higher-precision enclosure. A changed interval or defining polynomial starts a fresh memo. An insufficient remaining allowance uses ordinary computation so budget-limited results and partial-work statistics are preserved.

- `sign(x)`: the sign of a `Real`. `Outcome<StrictSign>`.
- `compare(a, b)`: the sign of `a - b`, as `Outcome<StrictOrder>`; `ProvenCoincident` means equal. It must give exactly the outcome, statistics and work of `sign(&(a - b))`.
- `enclose(x, max_width)`: a proven ball around `x`. A ball's width is `hi - lo`, twice its radius. It refines until the width is at most `max_width` or the budget runs out, and returns the tightest ball it reached either way. `Err(Refusal)` only when `x` is undefined (for the same reasons as above) or the deadline has passed. The caller checks the width. It counts as one decision: in `refined` when the width was reached, in `not_exact` when it was not, in `refused` on an error. `max_width` must be positive; zero or negative is a programming error and may panic.
- `enclose_at(x, precision)`: a ball around `x` evaluated at working precision `precision`, whatever width that precision reaches. It takes `&self`, is not a decision and does not change decision statistics. A square-root argument whose ball reaches below zero gives the square root of the ball's non-negative part, as `Ball::sqrt` does, so `x` must be known to be defined (every constant of a returned trace is). `Err(Refusal)` when that precision cannot separate a divisor of `x` from zero (the caller asks again at a higher precision) or the deadline has passed. Piece 4's trace pieces use it to enclose their constants at the precision each call asks for.
- `real_roots(p)`: all real roots of a `UPoly`, in increasing order, each with its multiplicity in `p`, each isolated from every other root. Rational roots come back exact. A root that is not rational comes back with a defining polynomial irreducible over Q whenever the square-free factor holding it has degree at most 4: a quartic that splits into two rational quadratics is split, so a quadratic irrational root is known as one and decisions on it can use exact arithmetic in its quadratic field. The zero polynomial gives `Refused("numbers.zero-polynomial")`. A non-zero constant has no roots. The outcome is `Proven` (the root list) or `Refused("numbers.budget-exhausted")`; it is never `NotExact` or `ProvenCoincident`, because a partial root list has no meaning. Narrowing an isolating interval, whether to separate roots or to decide whether a root is rational (down to the width at which at most one rational whose denominator divides the leading coefficient fits), uses steps that converge quadratically and account for both the root's magnitude and nearby root separation (safeguarded or interval Newton steps, for example), each kept only when exact signs confirm the root is still bracketed, and falls back to bisection only when a step is rejected. The work then grows with the logarithm of the bits needed, so large coefficients and scaled variables stay within the default budget.
- `certify_root(system, candidate)`: either `Proven(BoxRoots::One(root))`, a proof that the closed box holds exactly one solution and that it is non-singular; or `Proven(BoxRoots::None)`, a proof that the box holds no solution; or `Refused` with `numbers.root-not-certified` when neither can be shown within the budget. That includes a box holding a singular solution (a tangent contact or a double root) and a box holding more than one solution. This method never returns `NotExact` or `ProvenCoincident`. A root on the boundary of the closed box counts as inside it; such a box may be refused, since a refusal is never wrong. The certified `root_box()` may be a sub-box of the candidate. When the work or precision limit is reached the code is `numbers.root-not-certified`; when the deadline has passed it is `numbers.deadline`. `system` and `candidate` of different dimensions are a programming error and may panic.
- `note_fact_decision()`: counts a decision that piece 2 settled from a recorded fact, adding one to `decisions` and one to `coincident_fact`, so the statistics cover every decision.
- `set_deadline(Some(instant))`: decisions started after that instant return `Refused("numbers.deadline")` without doing any work. Decisions in progress check during work and refuse when the deadline passes; they never retain a partial result as `NotExact`.
- `check_deadline()`: returns `Err(Refusal)` with `numbers.deadline` after that instant, otherwise `Ok(())`. It does not change decision statistics.
- `stats()` and `take_stats()`.

How a sign is decided is the creator's choice, within these rules:

1. **Filter first.** Most decisions on ordinary inputs must finish from a cached `Interval` enclosure with no allocation.
2. **Then more precision.** When the filter cannot decide, refine with `Ball` arithmetic at growing precision, up to `max_precision_bits`.
3. **Then an exact method.** To prove zero, use one of:
   - exact rational evaluation, when the expression is rational;
   - exact arithmetic in a tower of real quadratic extensions, when every leaf is rational, a square root or a root of degree 2: the value is computed exactly in the field `Q(sqrt(r1))(sqrt(r2))...` that its square roots generate, each radicand first tested for being a square in the field built so far. This decides zero and sign exactly, and a decision tries it before building any annihilating polynomial or separation bound;
   - an exact algebraic method;
   - a proven root-separation bound: a number `B > 0` computed from the expression's structure (degrees, heights, leaf polynomials, system sizes) such that a non-zero value has magnitude at least `B`. Published constructive bounds of this kind exist for expressions with square roots and algebraic leaves (for example the bounds of Burnikel, Fleischer, Mehlhorn and Schirra, and the Li-Yap degree-measure bound). The creator picks the bound, cites its source in the code, and must make sure it is valid for every leaf kind, including system-root coordinates.
4. **Then stop.** When the budget runs out, return `NotExact` with the tightest proven bound reached. The code names the limit that stopped the decision first: `numbers.precision-limit` when the next step would need more than `max_precision_bits`, `numbers.budget-exhausted` when the next step would need more than `max_work`, whichever comes first in the decision's sequence of steps. An exact method of rule 3 is an attempt beside refinement, never a stop: when building its proof (an annihilating polynomial, a separation bound, an algebraic identity) would pass `max_work`, the decision skips that attempt and keeps refining. So a decision ends `NotExact` only when the next refinement step would pass a limit, and a non-zero value whose sign refinement within both limits can show is always `Proven`.

A fixed cutoff such as "smaller than 2^-180 means zero" is not a proof and is forbidden.

### Budgets and work units

`BudgetPolicy` has two limits per decision: `max_precision_bits` and `max_work`.
With `max_work = 0`, no decision work runs, not even the filter, and a sign returns `NotExact` with a bound from the cached enclosure.

A work unit is counted the same way on every machine, from the sizes of the numbers involved and not from the clock:

- an operation on doubles (one interval operation) costs 1 unit;
- adding, subtracting or comparing two exact numbers of `a` and `b` 64-bit limbs costs `a + b` units;
- multiplying or dividing them costs `a * b` units, whatever multiplication algorithm is actually used;
- ball operations at precision `p` count as multi-limb operations on numbers of `ceil(p / 64)` limbs.

The creator may refine this count (for example, square roots), but it must be deterministic, documented in the code, and never read the clock.

`BudgetPolicy::default()` is `max_precision_bits = 16384` and `max_work = 1_000_000_000`. These values may be retuned later by the architect; changes are logged in new `docs/interface-changes/` fragments.

`with_operation_budget` retains root-enclosure certificates already paid for in
that request. The key includes the exact defining polynomial, original rational
isolating interval and requested working precision. Reuse returns the same
published ball without reserving its arithmetic again against the whole-request
allowance. Cancellation checks and the per-decision cold-proof charge still
apply. A cached certificate from an earlier request must be paid for again.
Receipts are cleared on normal return and panic unwinding; they never transfer
between requests or threads. `DecisionStats::work` retains per-decision
accounting; `OperationWork::used` counts request-wide reservations and checkpoints.

The request also retains complete nonlinear root lists across deciders and
bounded-memo eviction. Their keys include the full polynomial and both
`BudgetPolicy` limits: policy-bounded exact shortcuts can change the isolation
certificate and its cold fee. Reuse preserves root values, ordering,
multiplicities, per-decision cold work and deadline checks. The request retains
the certificate data, not only a marker waiving another arithmetic reservation.
Failed or partial isolation is not recorded.

On complete rational-search branches, a small integer quadratic may obtain its
Descartes variations directly from the three exact transformed coefficients.
Integer admission proves all intermediates fit the native signed arithmetic;
it is not a geometric tolerance or a sign estimate. Other inputs retain the
filtered Taylor-shift predicate and its exact fallback. The zero/one-variation
proof rules and coarse interval selectors are unchanged.


### Statistics

`DecisionStats` counts every decision a `Decider` has made, by how it ended:

- `filter`: proven by the `f64` filter alone;
- `refined`: proven after a precision increase;
- `exact`: proven non-zero by an exact method;
- `coincident_exact`, `coincident_separation`, `coincident_fact`: proven zero, by proof kind (`coincident_exact` covers `ExactArithmetic` and `Algebraic`);
- `not_exact` and `refused`;
- `max_bits`: the highest precision any decision used;
- `work`: total work units.

`decisions` is the sum of all outcome counts. `real_roots` and `certify_root` count as one decision each: in `exact` when `Proven`, in `refused` when `Refused`.

## Construction facts

`FactBook<K>` stores facts about objects named by keys of type `K` (piece 2 uses its geometry ids). It is append-only and deterministic.

- `FactKind`: `Coincident` (same point set), `Coaxial`, `Concentric`, `Parallel`, `Tangent`, `Incident` (the first object lies on the second).
- `Fact<K>` holds the kind, the two subjects `a` and `b`, an optional third key `at` (for example the curve or point along which two surfaces are tangent), and `FactOrigin`: `Construction(name)` with a static name of the constructing routine, or `Proof(ZeroProof)`.
- `record(fact)` returns `true` when the fact is new and `false` when an equal fact is already present (it keeps the first). Two facts are equal when they have the same kind, the same `at`, and the same subjects, where the subjects of the symmetric kinds (all but `Incident`) are compared as an unordered pair; the origin is ignored. A fact with `a == b` is not stored and `record` returns `false`.
- `holds(kind, a, b)` answers:
  - for `Coincident`, `Coaxial`, `Concentric`, `Parallel`: true when `a` and `b` are linked by a chain of facts of that kind, in either direction (these kinds are equivalence relations; every key is related to itself);
  - for `Tangent`: true when a tangent fact names `a` and `b` in either order (so false when `a == b`);
  - for `Incident`: true when there is a chain `a` incident on `x1`, `x1` incident on `x2`, and so on up to `b`, where any step may also cross a `Coincident` link (false when `a == b` with no such chain);
  - no other kind is closed through `Coincident` links: `Coincident(a, a2)` and `Parallel(a2, b)` do not make `holds(Parallel, a, b)` true.
- `facts_about(k)` lists the recorded facts that name `k` as `a`, `b` or `at`, in recording order.
- `iter()` lists all facts in recording order.

Equivalence queries use union-by-size components instead of walking every
historical relation. Incidence queries traverse the directed graph of
`Coincident` components. Merging two components rewires incoming and outgoing
edges and retains any incidence self-edge. Coincidence alone never proves
incidence. The original facts, their recording order and provenance stay intact;
the index does not remove geometry, merge handles or change decision accounting.


## Marks

A mark is the record of one "not exact" answer, carried on results.

- `BoundKind`: `Distance` (the bound is a length in model units), `Angle` (radians), `Value` (the bound applies to the decided quantity itself and has no geometric meaning).
- `Mark` holds a code, a free-text detail, a non-negative rational bound and its kind. `Mark::from_not_exact(ne, kind, bound, detail)` keeps the `NotExact` code; the caller supplies the bound converted to the stated kind (when the kind is `Value`, pass `ne.bound`).
- `MarkLog` is an ordered list of marks: `push`, `extend` (appends another log), `iter`, `len`, `is_empty`, and `worst(kind)`, the largest bound of that kind if any.

`Answer<T>` holds a value and a `MarkLog`. `Answered<T>` is `Result<Answer<T>, Refusal>`. `Answer::exact(value)` builds one with an empty log, and `is_exact()` says whether the log is empty.

## Refusal and not-exact codes

| Code | When |
| --- | --- |
| `numbers.budget-exhausted` | Work limit reached (in `NotExact`, or in `Refused` from `real_roots`) |
| `numbers.precision-limit` | Precision limit reached (in `NotExact`) |
| `numbers.division-by-zero` | A divisor is proven zero |
| `numbers.divisor-not-separated` | A divisor cannot be separated from zero within the budget |
| `numbers.sqrt-of-negative` | A square root of a proven negative value |
| `numbers.sqrt-not-separated` | A square-root argument cannot be proven non-negative within the budget |
| `numbers.zero-polynomial` | `real_roots` of the zero polynomial |
| `numbers.root-not-certified` | `certify_root` could prove neither one root nor no root |
| `numbers.deadline` | The operation's wall-clock deadline passed |
| `numbers.geometric-gap-unproven` | Missing or invalid geometric evidence, units, binding, scope or dependencies; a geometric consumer also uses it for an unproved conversion or theorem |
| `numbers.geometric-gap-exceeded` | The available certified geometric upper bound is greater than or equal to the strict recorded modelling cap; this does not prove the true gap exceeds the cap |

## Exact public interface

The crate re-exports the certificate APIs above and the numeric baseline below at its root. Arithmetic operators are implemented for owned values and for references (`&a + &b`) unless noted. Fields not listed are private.

```rust
// Exact integers and rationals
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Integer;
impl Integer {
    pub fn zero() -> Self;
    pub fn one() -> Self;
    pub fn from_i64(v: i64) -> Self;
    pub fn from_u64(v: u64) -> Self;
    pub fn parse(s: &str) -> Result<Self, ParseNumberError>;
    pub fn to_i64(&self) -> Option<i64>;
    pub fn sign(&self) -> Sign;
    pub fn is_zero(&self) -> bool;
    pub fn abs(&self) -> Self;
    pub fn bit_length(&self) -> u64;              // bits of |self|; 0 for zero
    pub fn pow(&self, exp: u32) -> Self;
    pub fn gcd(&self, other: &Self) -> Self;      // non-negative; gcd(0, 0) = 0
    pub fn div_rem_floor(&self, d: &Self) -> Option<(Self, Self)>; // None when d = 0
    pub fn mul_pow2(&self, k: u64) -> Self;
}
// Integer: Add, Sub, Mul, Neg; Display (decimal); FromStr; From<i64>, From<u64>, From<i32>

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ParseErrorKind { Empty, InvalidSyntax, ZeroDenominator, ExponentTooLarge }
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ParseNumberError { pub kind: ParseErrorKind }
// ParseNumberError: Display, std::error::Error

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Rational;
impl Rational {
    pub fn zero() -> Self;
    pub fn one() -> Self;
    pub fn from_integer(n: Integer) -> Self;
    pub fn from_i64(v: i64) -> Self;
    pub fn new(numerator: Integer, denominator: Integer) -> Option<Self>;
    pub fn from_ratio(numerator: i64, denominator: i64) -> Option<Self>;
    pub fn parse(s: &str) -> Result<Self, ParseNumberError>;
    pub fn numerator(&self) -> &Integer;
    pub fn denominator(&self) -> &Integer;       // always positive
    pub fn sign(&self) -> Sign;
    pub fn is_zero(&self) -> bool;
    pub fn is_integer(&self) -> bool;
    pub fn abs(&self) -> Self;
    pub fn recip(&self) -> Option<Self>;
    pub fn checked_div(&self, other: &Self) -> Option<Self>;
    pub fn pow(&self, exp: i32) -> Option<Self>;
    pub fn floor(&self) -> Integer;
    pub fn ceil(&self) -> Integer;
    pub fn mul_pow2(&self, k: i64) -> Self;
    pub fn from_f64_exact(v: f64) -> Option<Self>;
    pub fn to_f64_nearest(&self) -> f64;
    pub fn to_decimal_exact(&self) -> Option<String>;
    pub fn to_decimal_bounds(&self, significant_digits: u32) -> (String, String);
}
// Rational: Add, Sub, Mul, Div (panics on zero), Neg; Display; FromStr;
//           From<i64>, From<i32>, From<Integer>

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Sign { Negative, Zero, Positive }

// Enclosures
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Interval;
impl Interval {
    pub fn new(lo: f64, hi: f64) -> Option<Self>;  // None if NaN or lo > hi
    pub fn point(v: f64) -> Option<Self>;          // None if NaN
    pub fn entire() -> Self;
    pub fn from_rational(r: &Rational) -> Self;
    pub fn from_integer(n: &Integer) -> Self;
    pub fn lo(&self) -> f64;
    pub fn hi(&self) -> f64;
    pub fn contains(&self, v: f64) -> bool;
    pub fn contains_zero(&self) -> bool;
    pub fn is_finite(&self) -> bool;
    pub fn width(&self) -> f64;                    // rounded up
    pub fn mid(&self) -> f64;
    pub fn sign(&self) -> Option<StrictSign>;
    pub fn square(&self) -> Self;
    pub fn sqrt(&self) -> Option<Self>;
    pub fn hull(&self, other: &Self) -> Self;
}
// Interval: Add, Sub, Mul, Div, Neg (by value; Interval is Copy)

#[derive(Clone, Debug)]
pub struct Ball;
impl Ball {
    pub fn from_rational(r: &Rational, precision: u32) -> Self;
    pub fn from_ratio(numerator: &Integer, denominator: &Integer, precision: u32) -> Option<Self>;
    pub fn from_integer(n: &Integer, precision: u32) -> Self;
    pub fn precision(&self) -> u32;
    pub fn with_precision(&self, precision: u32) -> Self;
    pub fn mid(&self) -> Rational;
    pub fn rad(&self) -> Rational;
    pub fn lo(&self) -> Rational;
    pub fn hi(&self) -> Rational;
    pub fn contains(&self, r: &Rational) -> bool;
    pub fn contains_zero(&self) -> bool;
    pub fn is_exact_zero(&self) -> bool;
    pub fn sign(&self) -> Option<StrictSign>;
    pub fn to_interval(&self) -> Interval;
    pub fn hull(&self, other: &Self) -> Self;
    pub fn square(&self) -> Self;
    pub fn sqrt(&self) -> Option<Self>;
    pub fn recip(&self) -> Option<Self>;
    pub fn checked_div(&self, other: &Self) -> Option<Self>;
    pub fn pow_u32(&self, exp: u32) -> Self;
    pub fn pi(precision: u32) -> Self;
    pub fn pi_checked(precision: u32, dec: &Decider) -> Result<Self, Refusal>;
    pub fn exp(&self) -> Self;
    pub fn ln(&self) -> Option<Self>;
    pub fn sin(&self) -> Self;
    pub fn cos(&self) -> Self;
    pub fn atan(&self) -> Self;
    pub fn atan_checked(&self, dec: &Decider) -> Result<Self, Refusal>;
    pub fn asin(&self) -> Option<Self>;
    pub fn acos(&self) -> Option<Self>;
    pub fn atan2(y: &Self, x: &Self) -> Option<Self>;
}
// Ball: Add, Sub, Mul, Neg (no Div; use checked_div)

// Polynomials
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct UPoly;
impl UPoly {
    pub fn new(coefficients: Vec<Rational>) -> Self; // lowest degree first
    pub fn zero() -> Self;
    pub fn constant(c: Rational) -> Self;
    pub fn x() -> Self;
    pub fn from_roots(roots: &[Rational]) -> Self;   // monic product of (x - r)
    pub fn coefficients(&self) -> &[Rational];
    pub fn degree(&self) -> Option<usize>;          // None for the zero polynomial
    pub fn is_zero(&self) -> bool;
    pub fn leading(&self) -> Option<&Rational>;
    pub fn eval(&self, x: &Rational) -> Rational;
    pub fn eval_interval(&self, x: Interval) -> Interval;
    pub fn eval_ball(&self, x: &Ball) -> Ball;
    pub fn eval_real(&self, x: &Real) -> Real;
    pub fn derivative(&self) -> Self;
    pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>; // None for a zero divisor
    pub fn gcd(&self, other: &Self) -> Self;         // monic; gcd(0, 0) = 0
    pub fn square_free_factors(&self) -> Vec<(UPoly, u32)>;
    pub fn compose(&self, inner: &Self) -> Self;
    pub fn primitive_integer(&self) -> Vec<Integer>;
}
// UPoly: Add, Sub, Mul, Neg (references)

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Monomial;
impl Monomial {
    pub fn new(exponents: Vec<u32>) -> Self;
    pub fn exponents(&self) -> &[u32];
    pub fn degree(&self) -> u32;
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct MPoly;
impl MPoly {
    pub fn zero(nvars: usize) -> Self;
    pub fn constant(nvars: usize, c: Rational) -> Self;
    pub fn var(nvars: usize, index: usize) -> Self;
    pub fn from_terms(nvars: usize, terms: Vec<(Vec<u32>, Rational)>) -> Self; // like terms are summed
    pub fn nvars(&self) -> usize;
    pub fn is_zero(&self) -> bool;
    pub fn total_degree(&self) -> Option<u32>;
    pub fn terms(&self) -> Vec<(Monomial, Rational)>; // canonical order, no zero coefficients
    pub fn scale(&self, c: &Rational) -> Self;
    pub fn eval(&self, x: &[Rational]) -> Rational;
    pub fn eval_interval(&self, x: &[Interval]) -> Interval;
    pub fn eval_ball(&self, x: &[Ball]) -> Ball;
    pub fn eval_real(&self, x: &[Real]) -> Real;
    pub fn partial(&self, var: usize) -> Self;
    pub fn substitute(&self, var: usize, value: &MPoly) -> Self;
    pub fn resultant(&self, other: &Self, var: usize) -> Self;
    pub fn to_upoly(&self, var: usize) -> Option<UPoly>; // None if another variable appears
}
// MPoly: Add, Sub, Mul, Neg (references)

// Exact real numbers
#[derive(Clone, Debug)]
pub struct Algebraic;
impl Algebraic {
    pub fn from_rational(r: Rational) -> Self;
    pub fn as_rational(&self) -> Option<&Rational>;
    pub fn defining_polynomial(&self) -> UPoly;     // square-free, primitive integer coefficients
    pub fn isolating_interval(&self) -> (Rational, Rational);
    pub fn refined(&self, max_width: &Rational) -> Self;
    pub fn to_ball(&self, precision: u32) -> Ball;
}

#[derive(Clone, Debug)]
pub struct RealRoot { pub value: Algebraic, pub multiplicity: u32 }

#[derive(Clone, Debug)]
pub struct Real;
impl Real {
    pub fn zero() -> Self;
    pub fn one() -> Self;
    pub fn from_rational(r: Rational) -> Self;
    pub fn from_i64(v: i64) -> Self;
    pub fn from_algebraic(a: Algebraic) -> Self;
    pub fn from_root_coordinate(root: &IsolatedRoot, axis: usize) -> Self;
    pub fn sqrt(&self) -> Self;
    pub fn square(&self) -> Self;
    pub fn pow_u32(&self, exp: u32) -> Self;
    pub fn as_rational(&self) -> Option<Rational>;
    pub fn approx(&self) -> Interval;
    pub fn to_f64(&self) -> f64;
}
// Real: Add, Sub, Mul, Div, Neg; From<Rational>, From<i64>, From<Algebraic>

#[derive(Clone, Debug)]
pub struct PolySystem;
impl PolySystem {
    pub fn new(equations: Vec<MPoly>) -> Option<Self>; // None unless n equations in n variables, n >= 1
    pub fn equations(&self) -> &[MPoly];
    pub fn dimension(&self) -> usize;
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RootBox;
impl RootBox {
    pub fn new(lo: Vec<Rational>, hi: Vec<Rational>) -> Option<Self>; // None unless same length >= 1 and lo <= hi
    pub fn lo(&self) -> &[Rational];
    pub fn hi(&self) -> &[Rational];
    pub fn dimension(&self) -> usize;
}

#[derive(Clone, Debug)]
pub struct IsolatedRoot;
impl IsolatedRoot {
    pub fn system(&self) -> &std::sync::Arc<PolySystem>;
    pub fn root_box(&self) -> &RootBox;
    pub fn coordinate(&self, axis: usize) -> Real;
    pub fn enclose(&self, precision: u32) -> Vec<Ball>;
}

#[derive(Clone, Debug)]
pub enum BoxRoots { None, One(IsolatedRoot) }

// Outcomes
#[must_use]
#[derive(Clone, Debug, PartialEq)]
pub enum Outcome<T> {
    Proven(T),
    ProvenCoincident(Coincidence),
    NotExact(NotExact),
    Refused(Refusal),
}
impl<T> Outcome<T> {
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Outcome<U>;
    pub fn proven(&self) -> Option<&T>;
    pub fn is_proven(&self) -> bool;
    pub fn is_coincident(&self) -> bool;   // true for ProvenCoincident only
    pub fn is_not_exact(&self) -> bool;
    pub fn is_refused(&self) -> bool;
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum StrictSign { Negative, Positive }
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum StrictOrder { Less, Greater }

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum ZeroProof { ExactArithmetic, Algebraic, SeparationBound, ConstructionFact }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Coincidence { pub proof: ZeroProof }

#[derive(Clone, PartialEq, Debug)]
pub struct NotExact { pub bound: Rational, pub code: &'static str, pub detail: String }

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Refusal { pub code: &'static str, pub detail: String }
impl Refusal { pub fn new(code: &'static str, detail: impl Into<String>) -> Self; }
// Refusal: Display, std::error::Error

#[derive(Clone, Debug)]
pub struct Answer<T> { pub value: T, pub marks: MarkLog }
impl<T> Answer<T> {
    pub fn exact(value: T) -> Self;
    pub fn with_marks(value: T, marks: MarkLog) -> Self;
    pub fn is_exact(&self) -> bool;
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Answer<U>;
}
pub type Answered<T> = Result<Answer<T>, Refusal>;

// Budgets and the decider
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct BudgetPolicy { pub max_precision_bits: u32, pub max_work: u64 }
// BudgetPolicy: Default (16384 bits, 1_000_000_000 units)

#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct DecisionStats {
    pub decisions: u64,
    pub filter: u64,
    pub refined: u64,
    pub exact: u64,
    pub coincident_exact: u64,
    pub coincident_separation: u64,
    pub coincident_fact: u64,
    pub not_exact: u64,
    pub refused: u64,
    pub max_bits: u32,
    pub work: u64,
}

#[derive(Debug)]
pub struct Decider;
impl Decider {
    pub fn new(policy: BudgetPolicy) -> Self;
    pub fn policy(&self) -> &BudgetPolicy;
    pub fn set_deadline(&mut self, deadline: Option<std::time::Instant>);
    pub fn check_deadline(&self) -> Result<(), Refusal>;
    pub fn sign(&mut self, x: &Real) -> Outcome<StrictSign>;
    pub fn compare(&mut self, a: &Real, b: &Real) -> Outcome<StrictOrder>;
    pub fn enclose(&mut self, x: &Real, max_width: &Rational) -> Result<Ball, Refusal>;
    pub fn enclose_at(&self, x: &Real, precision: u32) -> Result<Ball, Refusal>;
    pub fn real_roots(&mut self, p: &UPoly) -> Outcome<Vec<RealRoot>>;
    pub fn certify_root(&mut self, system: &std::sync::Arc<PolySystem>, candidate: &RootBox) -> Outcome<BoxRoots>;
    pub fn note_fact_decision(&mut self);
    pub fn stats(&self) -> &DecisionStats;
    pub fn take_stats(&mut self) -> DecisionStats;
}
// Decider: Default (BudgetPolicy::default(), no deadline)

// Facts
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum FactKind { Coincident, Coaxial, Concentric, Parallel, Tangent, Incident }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum FactOrigin { Construction(&'static str), Proof(ZeroProof) }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Fact<K> { pub kind: FactKind, pub a: K, pub b: K, pub at: Option<K>, pub origin: FactOrigin }

#[derive(Clone, Debug)]
pub struct FactBook<K>;
impl<K: Copy + Ord> FactBook<K> {
    pub fn new() -> Self;
    pub fn record(&mut self, fact: Fact<K>) -> bool;
    pub fn holds(&self, kind: FactKind, a: K, b: K) -> bool;
    pub fn facts_about(&self, key: K) -> Vec<&Fact<K>>;
    pub fn iter(&self) -> std::slice::Iter<'_, Fact<K>>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
}
// FactBook: Default

// Marks
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BoundKind { Distance, Angle, Value }
#[derive(Clone, PartialEq, Debug)]
pub struct Mark { pub code: &'static str, pub detail: String, pub bound: Rational, pub kind: BoundKind }
impl Mark {
    pub fn new(code: &'static str, kind: BoundKind, bound: Rational, detail: impl Into<String>) -> Self;
    pub fn from_not_exact(not_exact: &NotExact, kind: BoundKind, bound: Rational, detail: impl Into<String>) -> Self;
}
#[derive(Clone, PartialEq, Debug, Default)]
pub struct MarkLog;
impl MarkLog {
    pub fn new() -> Self;
    pub fn push(&mut self, mark: Mark);
    pub fn extend(&mut self, other: &MarkLog);
    pub fn iter(&self) -> std::slice::Iter<'_, Mark>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn worst(&self, kind: BoundKind) -> Option<&Rational>;
}
```

The crate in `crates/catallax-numbers` implements this interface. The compiled crate is authoritative if a signature listing becomes stale; review and correct the listing with any API change.

## What defines correct behaviour

The test writer builds tests from the cases below. Every case has an answer that maths gives, so tests compare against that answer, never against a run of the code.
Independent oracles are allowed in tests: a third-party big-number crate used only as a dev-dependency, or tables of values generated by a separate program, as long as nothing of this piece's own code computes the expected answer.

### Parsing and exact arithmetic

1. Every accepted and rejected string listed under "Rational" above.
2. `parse("0.1") + parse("0.2") == parse("0.3")` exactly, while the same sum in doubles is not 0.3.
3. Normal form: `6/4` equals `3/2`; `-6/-4` is rejected by the grammar, but `Rational::new(-6, -4)` gives `3/2`; `0/5 == 0/1` and both hash the same; the denominator is always positive.
4. Field laws on random rationals with numerators and denominators up to 2^256: associativity, distributivity, `a / b * b == a` for non-zero `b`, `a - a == 0`.
5. `from_f64_exact(0.1)` is `3602879701896397 / 36028797018963968`. `from_f64_exact` then `to_f64_nearest` gives back every finite double, including subnormals and `-0.0` (which becomes 0).
6. `to_decimal_bounds`: `1/3` with 5 digits gives `("0.33333", "0.33334")`; `2/3` gives `("0.66666", "0.66667")`; `5/4` with 2 digits gives `("1.2", "1.3")`; `5/4` with 3 digits gives `("1.25", "1.25")`; `-1/3` with 3 digits gives `("-0.334", "-0.333")`; `10^30` with 3 digits gives `("1e+30", "1e+30")`; `1/10^7` with 2 digits gives `("1e-7", "1e-7")`. For random rationals, parsing both strings back gives `lo <= value <= hi`.
7. Integers of 10,000 digits: `(a * b) / b == a`, `gcd` matches the Euclidean definition, `div_rem_floor` returns `q = floor(a / d)` and `r = a - q * d`, with `0 <= |r| < |d|` and non-zero `r` having the sign of `d`. For `a = 7`, `d = -3`, `(q, r) = (-3, -2)`.

### Intervals

8. `Interval::from_rational(1/10)` has `lo < 1/10 < hi` with `lo` and `hi` adjacent doubles. `from_rational(1/2)` is the single point 0.5.
9. The sum of ten copies of `from_rational(1/10)` contains 1.
10. Containment under random operations: for random rationals `x`, `y` and each operation, the interval result contains the exact rational result.
11. `[1e300, 1e300] * [1e300, 1e300]` has sign `Positive` and an infinite upper end; it never turns into NaN.
12. `from_rational(1e-320)` (a subnormal range) contains the value and does not claim `Positive` unless its lower end is above zero.
13. Division by `[-1, 1]` gives the whole line. `sqrt([-2, -1])` is `None`. `sqrt([-1, 4])` is `[0, 2]`.

### Balls and elementary functions

14. `Ball::pi(p)` for `p` = 64, 1000 and 10,000 contains the first 3,000 decimal digits of pi as published (both the truncation and the truncation plus one unit in the last digit place lie on the correct sides).
15. `sin(pi)` contains 0; `cos(pi)` contains -1; `atan(1) * 4` overlaps `pi`; `exp(ln(2))` contains 2; `asin(1)` contains pi / 2; `atan2(1, -1)` contains 3 pi / 4.
16. `exp(1)` contains e to 1,000 digits; `ln(2)` contains ln 2 to 1,000 digits (published constants).
17. Containment under random inputs for every function, checked against an independent high-precision oracle.
18. Domain edges: `ln` of a ball containing 0 is `None`; `asin` of a ball reaching past 1 is `None`; `atan2` of two balls containing 0 is `None`; `sqrt` of a ball entirely below zero is `None`.
19. The same ball computed twice, and on two machines, has the same midpoint and radius.

### Signs that must be proven, with known answers

20. The classic near-degenerate orientation test: for points `p = (0.5 + i * 2^-53, 0.5 + j * 2^-53)` with `i, j` in `0..256`, `q = (12, 12)`, `r = (24, 24)`, the sign of `(q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)`, built as a `Real` from the exact rationals, must match exact rational evaluation for all 65,536 points. A plain double evaluation gets many of these wrong; this piece must get all of them right.
21. The same with the 3D orientation determinant of four points on random exact grids near a common plane.
22. `sign(parse("0.1") + parse("0.2") - parse("0.3"))` is `ProvenCoincident` with proof `ExactArithmetic`.
23. `sign(sqrt(1000001) + sqrt(999999) - 2 * sqrt(1000000))` is `Proven(Negative)` (the value is about -2.5e-10). With 10^20 in place of 10^6 the value is about -2.5e-31; still `Proven(Negative)` under the default budget.
24. `sign(sqrt(10^40 + 1) - 10^20)` is `Proven(Positive)` (about 5e-21).
25. Large exact cancellations: with `a = 2^1100` as a rational, `sign(a * a - a * a + 1)` is `Proven(Positive)`, and `sign(a * a - a * a)` is `ProvenCoincident`.

### Hidden zeros that must be proven coincident

26. `sqrt(2) * sqrt(2) - 2`.
27. `sqrt(2) + sqrt(3) - sqrt(5 + 2 * sqrt(6))`.
28. `sqrt(3 + 2 * sqrt(2)) - 1 - sqrt(2)`.
29. With `g = (1 + sqrt(5)) / 2`: `g * g - g - 1`.
30. For random non-negative rationals `x`, `y`: `sqrt(x) * sqrt(y) - sqrt(x * y)`.
31. For the root `a` of `x^2 - 2` in `(1, 2)` from `real_roots`: `compare(Real::from_algebraic(a), sqrt(2))` is `ProvenCoincident`; `compare` with `parse("1.41421356")` is `Proven(Greater)`.

Each must end `ProvenCoincident` under the default budget, never `Proven`.

### Budget behaviour

32. With `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }`, case 24 must end `NotExact` (or `Proven(Positive)` if an exact method fits in that budget), never `Proven(Negative)` and never `ProvenCoincident`. If it is `NotExact`, its bound must be at least the true magnitude, about 5e-21.
33. For every hidden zero in cases 26 to 30 under a tiny budget: the outcome is `ProvenCoincident` or `NotExact`, never `Proven`.
34. For random expressions of depth up to 6 over `+ - * / sqrt` with small rational leaves, under random budgets: every `Proven` agrees with an independent high-precision oracle that also proves the sign, every `NotExact` bound contains the oracle's value, and an expression the oracle finds undefined (a negative square-root argument or a zero divisor) never ends `Proven`, `ProvenCoincident` or `NotExact`. A refusal is always allowed.
35. A `Decider` whose deadline is in the past returns `Refused("numbers.deadline")` for every decision and does no work.
36. Repeating any decision gives the identical outcome and adds identical work to the statistics. Running the same list of decisions on eight threads (one `Decider` each) gives identical outcomes and statistics.

### Undefined values

37. `sign(1 / (sqrt(2) * sqrt(2) - 2))` is `Refused("numbers.division-by-zero")`.
38. `sign(sqrt(-1))` is `Refused("numbers.sqrt-of-negative")`.
39. `sign(sqrt(sqrt(2) * sqrt(2) - 2))` is `ProvenCoincident` (the square root of a proven zero is zero). Under the default budget, `sign((sqrt(2) + sqrt(3) + sqrt(5) + sqrt(7))^2 - (17 + 2 * (sqrt(6) + sqrt(10) + sqrt(14) + sqrt(15) + sqrt(21) + sqrt(35))))` and `sign(sqrt(3 + 2 * sqrt(2)) - 1 - sqrt(2))` are `ProvenCoincident`, and `sign(sqrt(2) + sqrt(3) + sqrt(5) - 5.3824)` is `Proven(Negative)`.

### Roots of univariate polynomials

40. `real_roots(x^2 - 2)` gives two simple roots, near -1.41421 and 1.41421, neither rational. Under the default budget, `real_roots((N^2 x^2 - 2) * (N^2 x^2 - 3))` with `N = 10^30 * (2^607 - 1)` gives four simple roots, `-sqrt(3) / N`, `-sqrt(2) / N`, `sqrt(2) / N` and `sqrt(3) / N`, none rational. With `e = 10^-30`, `a = 4e - e^2` and `b = 8e - e^2 - 12`, `real_roots((x^2 - a) * (x^2 - b))` (expanded) gives the two roots `-sqrt(a)` and `sqrt(a)`, each with a defining polynomial of degree 2, and for each root `r`, `sign(1 - (sqrt(4 - r^2) - 3 + e)^2)` is `ProvenCoincident` under the default budget (`4 - a = (2 - e)^2`).
41. `real_roots` of the product of `(x - i)` for `i` in 1..=20 (Wilkinson's polynomial) gives exactly 20 simple roots, each exactly the rational `i`.
42. `real_roots((x - 1)^3 * (x + 2)^2 * (x^2 + 1))` gives `-2` with multiplicity 2 and `1` with multiplicity 3, both exact.
43. Mignotte polynomials `x^n - 2 * (a * x - 1)^2`: for `a = 10` and `n` = 20 and 21, the two roots near `1 / a` are closer together than `2 * a^(-(n + 2) / 2)`; they must both be found and isolated from each other. The full root counts are given in "Values checked for this description" below.
44. Chebyshev polynomials `T_n` for `n` up to 64: exactly `n` simple roots, each inside `(-1, 1)`. The list is in increasing order, so its `k`th entry (counting from 1) overlaps `cos((2(n - k) + 1) pi / (2n))` computed with `Ball`.
45. Random polynomials with known integer and rational roots multiplied by irreducible quadratics with no real roots: every listed root found exactly, nothing else.
46. `real_roots(0)` is `Refused("numbers.zero-polynomial")`; `real_roots(5)` is `Proven` with no roots.

### Certified roots of systems

47. System `x^2 + y^2 - 1 = 0`, `x - y = 0`; box `[1/2, 1] x [1/2, 1]`: `Proven(One)`, and `compare(coordinate(0), sqrt(1/2))` is `ProvenCoincident`.
48. Same system, box `[-1/10, 1/10] x [-1/10, 1/10]`: `Proven(None)`.
49. Same system, box `[-1, 1] x [-1, 1]` (two roots inside): `Refused("numbers.root-not-certified")`.
50. System `x^2 + y^2 - 1 = 0`, `y - 1 = 0` (a tangent, double root at `(0, 1)`); a box around `(0, 1)`: `Refused("numbers.root-not-certified")`, never `Proven(One)`.
51. Three spheres of radius 5 centred at `(0, 0, 0)`, `(6, 0, 0)` and `(0, 6, 0)`: their two common points are `(3, 3, ±sqrt(7))`. Boxes around each give `Proven(One)`, and the `z` coordinate compared with `sqrt(7)` (or `-sqrt(7)`) is `ProvenCoincident`.

### Facts

52. `record(Coaxial(a, b))` then `holds(Coaxial, b, a)` is true. Adding `Coaxial(b, c)` makes `holds(Coaxial, a, c)` true. `holds(Coaxial, a, a)` is true for any key.
53. `Tangent(a, b)` and `Tangent(b, c)` do not make `holds(Tangent, a, c)` true.
54. `Incident(p, c)` and `Incident(c, s)` make `holds(Incident, p, s)` true but not `holds(Incident, s, p)`. `Incident(p, c)` and `Coincident(c, d)` make `holds(Incident, p, d)` true.
55. Recording a fact twice returns `true` then `false`, and `len` grows by one. `iter` and `facts_about` follow recording order.
56. A book with 100,000 facts answers `holds` for chains of length 1,000.

### Marks and answers

57. `MarkLog::worst(Distance)` returns the largest distance bound and ignores other kinds; an empty log returns `None`. `extend` keeps order. `Answer::exact(v).is_exact()` is true.

### Scale

58. Scale: for cases 20 and 21, multiplying every point coordinate by `10^k` for `k` in `-30..=30` leaves every sign unchanged (the determinants are homogeneous). For cases 22 to 31, multiplying the whole expression by the rational `10^k` (building `Real::from_rational(10^k) * E`) leaves every outcome kind unchanged.

### Volume of tests

The foundation is tested hardest (DESIGN.md, proof level). The generated families (cases 4, 10, 17, 20, 21, 30, 34, 45 and 58) together must run at least 30,000 maths-answered sign or containment questions per test run, from fixed seeds so runs are repeatable.

## Values checked for this description

The architect checked these numbers with an independent high-precision tool when writing this file:

- case 23: `sqrt(1000001) + sqrt(999999) - 2000` is about `-2.50000000000078e-10`; with `10^20` it is about `-2.5e-31`;
- case 24: `sqrt(10^40 + 1) - 10^20` is about `5e-21`;
- cases 27 and 28 evaluate to zero at 80 digits;
- case 43: the root counts and close-root gaps for `a = 10`, `n = 20` and `n = 21` are listed below.

| Polynomial | Real roots | Roots near `1 / a` | Gap | `2 * a^(-(n + 2) / 2)` |
| --- | --- | --- | --- | --- |
| `x^20 - 2 * (10x - 1)^2` | 4: about -1.3529322, 0.0999999999929289, 0.1000000000070711, 1.3306538 | 2 | about 1.414e-11 | 2e-11 |
| `x^21 - 2 * (10x - 1)^2` | 3: about 0.0999999999977639, 0.1000000000022361, 1.3106253 | 2 | about 4.472e-12 | about 6.325e-12 |

For `n >= 3` and integer `a >= 3`, Descartes' rule of signs together with the sign changes at `0`, `1/a` and `1` gives exactly four real roots for even `n` (one negative) and exactly three for odd `n` (none negative). The test writer may extend the family with an independent root counter.

## Speed: measured, never gated

No test fails for speed. A benchmark target (`cargo bench` or a `benches/` binary) measures and prints:

1. the share of decisions that finish in the `f64` filter on random general-position inputs of cases 20 and 21 (the expectation is above 99%, but it is reported, not asserted);
2. time per decision for those inputs, split by proof path;
3. time and precision reached for the near-degenerate families (cases 20, 23, 24) and for hidden zeros (cases 26 to 30);
4. time of `real_roots` against degree for Wilkinson, Mignotte and Chebyshev polynomials up to degree 100;
5. time of `certify_root` for case 51;
6. time to produce `pi` to 1,000 and 10,000 digits;
7. `FactBook::holds` time at 100,000 facts.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Points, vectors, surfaces, curves and poses (piece 2).
- Finding candidate boxes for system roots (piece 4).
- Geometric residual conversions, topology and whole-solid theorem verification. This piece records the modelling setting and proof claims and checks their arithmetic and dependency structure.
