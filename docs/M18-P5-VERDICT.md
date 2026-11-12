# Authored P5 at radius 2: exact field overflow verdict

Baseline: `4d51b5c36bdad29632a448585290cfaa627098ba` (13 November 2026).
The input is the public `parts/P5-flanged-housing.json` with its authored fillet
radius **2**. This note does not admit the recipe and changes no compiler
acceptance, tolerance, work allowance, face registry or answer. It records where
the field overflows, with the exact enclosure and figures, and what a build must
supply.

## Verdict

The field's `fillet.unsupported-overflow` at radius 2 is correct for the current
representation, and the geometry is not impossible. The plain torus fillet's
plane rail, radius `R + r = 14 mm`, runs through every countersink mouth, whose
near rim lies at `o - m = 13 mm`. The rail crosses each mouth circle transversally
at two exact points and reaches 1 mm inside the 4 mm throat. Emitting the plain
torus would place the plane contact, and with it part of the fillet surface, over
the hole. The correct solid needs the four rim-contact envelope faces that the
B-rep builds; the field has neither the compiled incidence nor the material and
named-boundary integration for them.

The last radius the field's certified clearance proof admits is `r < 1`. At
`r = 1` the rail is exactly tangent to every mouth, and at `r = 2` it crosses.

## Exact figures

**What overflows:** the plain torus's plane rail and the material it carries,
not the rolling ball, the numeric precision or the work allowance. **At which
bound:** the strict clearance bound `d > R + r + m = 18 mm`, missed by exactly
1 mm, or equivalently the untrimmed-torus radius bound
`r <= o - m - R = 1 mm`, exceeded by 1 mm at the authored radius.

The housing is a revolution about `y`; the flange top is the plane `y = 5 mm`.
Each of the four countersink mouths is a circle of radius `m = 4 mm` centred
`o = 17 mm` from the housing axis in that plane. The body radius where the
flange-top meets the body-wall is `R = 12 mm`.

| Quantity | Exact value |
| --- | --- |
| Body radius `R` | 12 mm |
| Authored ball radius `r` | 2 mm |
| Mouth radius `m` | 4 mm |
| Mouth-axis offset `o` | 17 mm |
| Plain torus plane rail `R + r` | 14 mm |
| Mouth near rim `o - m` | 13 mm |
| Rail penetration at the nearest point | `14 - 13 = 1 mm` |
| Rail/mouth crossings per mouth | 2: `(469/34, 5, ±21*sqrt(15)/34)` mm |
| Crossing cosine on the rail | `cos(theta0) = 469/476 = 67/68` |
| Crossing half-angle | `theta0 ~= 9.8382 deg` |
| Untrimmed rail arc inside one mouth | `2*theta0 ~= 19.6765 deg` |
| Four-mouth rail arc | `8*theta0 ~= 78.7058 deg`, 21.86% of the ring |
| Torus pieces left between mouths | four arcs of `90 - 2*theta0 ~= 70.3235 deg` |
| Ball-centre dip at the deepest station | `2 - sqrt(3) ~= 0.267949 mm` below the rail |

Derivation. The rail circle is `x^2 + z^2 = 196`; a mouth circle is
`(x - 17)^2 + z^2 = 16`. Subtracting gives `34x - 289 = 180`, so
`x = 469/34` and `z^2 = 196 - (469/34)^2 = 6615/1156`. The same rational
appears in the rolling region's transition formula
`(C*C + o*o - m*m) / (2*C*o)` with `C = R + r = 14`: both give `67/68`. Because
`|14 - 4| = 10 < 17 < 18 = 14 + 4`, each pair of circles has exactly two
crossings; the strict inequality is the transversal case, not a tangency.

The B-rep's `RollingEnvelopeDefinition` transitions are exactly these crossings.
At `u = ±1` its mouth foot is `(469/34, 5, ±21*sqrt(15)/34)`; at `u = 0` the
ball centre is `(14, 5 + sqrt(3), 0)` and it touches the mouth's nearest point
`(13, 5, 0)`. The plain torus holds the ball centre at height `5 + r = 7` over
the whole ring, so over each mouth arc it is up to `2 - sqrt(3) mm` too high;
that is where its material overflows.

## Why the clearance proof cannot cover it

The finite-round clearance certificate (introduced by PR 409, implemented as
`finite_round_clear` in `crates/catallax-field/src/fillet.rs`) separates the
whole finite round primitive from the closed radial band `[R, R + r]` only when
either

- `d > R + r + m`, or
- `R > m` and `d < R - m`,

where `d` is the distance between the two parallel axes. Both are strict; an
equality is reported as unresolved, not clear.

Enclosure. The certificate encloses the material-change region in the closed
radial band `[R, R + r] = [12, 14] mm`, or `[144, 196]` squared. The countersink
cone's finite primitive is contained in the tube of radius `m = 4 mm` about its
axis at offset `o = 17 mm`; its nearest material point is the mouth rim at radial
`o - m = 13 mm`, or `169` squared. Since `144 < 169`, the two closed enclosures
intersect, and the strict disjunction above is the only exit. This is a
containment argument, not a sampled measurement.

For the authored P5, `d = o = 17` and `R + r + m = 18`: the first bound fails by
exactly 1 mm (`289` against `324` in squared form, a deficit of `35`). The second
bound would need `d < R - m = 8`, and 17 is 9 mm beyond it. Therefore the
certificate returns false at the countersink cone, and `fillet_domain` reports
`fillet.unsupported-overflow` with the "needs trimmed continuation" detail.

The same numbers give the radius threshold. For `r = 1`, `R + r + m = 17`, an
exact equality and external tangency at `(13, 5, 0)`; for `r = 1/2`,
`R + r + m = 16.5 < 17`, and the proof holds. That is why the radius-1/2 edit
builds and every `r > 1` that reaches the mouth refuses.

## What a field build requires

The missing work is the same three obligations the rolling-region note carries,
now with the specific figures:

1. **Recipe incidence.** Discover the four mouths from the recipe's Boolean
   material and construct the four envelope definitions (body 12, mouth 4,
   offset 17, radius 2, edge `[1, 1]`). Those definitions are valid at the
   authored radius; the existing region and this verdict's test construct them.
2. **Material coverage.** Prove that the four envelope charts cover exactly the
   arcs between the rail crossings (`2*theta0` each) and that the four torus
   pieces fill the remaining `70.3235`-degree arcs, including the transition
   stations, mouth arcs and both rails. Angular coverage alone is not material
   coverage.
3. **Named boundary integration.** Register the envelope charts with exact
   origins, rays and distance coverage. The field's `Face` registry holds an
   analytic support polynomial and trim predicates; a fitted carrier or the
   plain torus cannot stand in, and no unregistered envelope leaf may answer a
   named-subset query.

The B-rep already builds this solid: `crates/catallax-features/tests/rim_contact_fillets.rs`
checks the authored radius 2, the `5/2` edit and translated, scaled and rotated
definitions. The built result has 22 faces: 4 rolling envelopes, 4 torus pieces,
3 planes, 4 cones and 7 cylinders, one blend name, a closed shell and Euler
genus five. Reusing that B-rep as the field answer would not be an independent
field build and is not claimed here. At `7/2` even the B-rep refuses, because the
ball would leave the mouth edge for the cone face.

## Executed evidence

```sh
cargo test -p catallax-field --test m18_p5_verdict --test m18_clearance --test rolling_region --test rolling_sectors
cargo test -p catallax-field --test fillets
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
git diff --check
```

The new `m18_p5_verdict` test pins, from the public recipe and the compiler:

- the authored dimensions `R = 12`, `r = 2`, `m = 4`, `o = 17`, plane `y = 5`
  and pattern count 4;
- the admitted radius-1/2 blend torus: `major - minor = 12`, `minor = 1/2`;
- the 1 mm rail/rim overlap and the two exact crossings `x = 469/34`,
  `z^2 = 6615/1156`, both satisfying the mouth circle exactly;
- `cos(theta0) = 67/68`, equal to the rolling region's transition formula;
- the exact envelope stations at `u = -1, 0, 1`, including the
  `2 - sqrt(3)` dip and the mouth foot `(13, 5, 0)`;
- compilation outcomes: `1/2` and `3/4` build; `1` and `2` refuse with
  `fillet.unsupported-overflow`.

Executed counts: `m18_p5_verdict` 3, `m18_clearance` 1, `rolling_region` 6,
`rolling_sectors` 3 and the field `fillets` target 6 tests, all passed. Strict
field Clippy, formatting and whitespace checks pass. The B-rep
`rim_contact_fillets` suite is the existing control cited above; this change
does not modify it and it was not rerun here. The existing duplicate CLI
binary-target manifest warning remains. No production code changed, so the 477
replay answers, proven answers, work counts, tolerances and worker-count behavior
are untouched; no new parallel route is added.

## Limits

No authored-P5 admission, no M18 exit, no section 4.17 smooth blend, no wider
sweeps or lofts, no STEP export, no two-kernel comparison and no mass authority
are claimed. This note is a verdict and a pin, not coverage.
