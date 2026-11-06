# Circular latitude jets for fillet STEP charts

## Scope

STEP chart export can now use affine UV jets for a native circular latitude
of a toroidal fillet. An actual written torus or torus-sheet face must have a
recorded `blend` origin to select this route. Unused fillet history does not
select it. The existing chart, pole, seam and material proofs still run.

The shortcut proves its conditions against the rounded geometry actually
written to STEP, not the exact source placement. For written surface centre
`o`, axis `a`, major radius `M`, and written conic centre `c`, normal `n`,
semi-axes `r` and `s`, it requires:

- `r = s > 0` exactly as rational written values;
- `n` parallel to `a`, proved by an exact rational cross product;
- `c - o` orthogonal to both written radial basis directions;
- the tube vector `(r - M, height)` is nonzero.

No proximity test or tolerance can replace those equalities. An ellipse,
offset circle, tilted circle, zero radius or tube-centre circle keeps the
existing generic path. Meridian seams also keep that path.

In the written orthonormal frames the parameter inverse is
`u(t) = sign(dot(n,a)) t + phase` and
`v(t) = atan2(height, r - M)`. Thus `u'` is exactly `+1` or `-1`, `v'` is zero,
and all higher Taylor coefficients are zero. This avoids interval dependency
between sine and cosine in repeated generic fits. Constant values still use
the existing angular lifts and interval inverse. Existing Taylor controls,
endpoint coupling, rounding bounds, fitting reservations, arithmetic metering
and final file-bound checks remain active. No 3D carrier or fillet changes.

## Observed coverage

At the unchanged 1M CLI request allowance, `sphere-dome` now exports a
raw-valid one-solid/five-face file at 945,971 work. The same fresh-main recipe
requires 2,214,090 work at an explicit 20M request. Both use 17 digits and
0.0000001 mm file tolerance. The candidate reports
5.656937652983913246864836148956529112911e-9 mm deviation and zero fitted
3D curves. Its existing export-only pole and torus-equator cuts remain.

Ten material witnesses check the cap and removed fillet material. An
independent washer integral gives 47.3569106407326 mm3; raw STEP gives
47.35691064032189 mm3. The independent OCCT construction gives
47.356912004414255 mm3, within the unchanged comparison allowance. These
checks do not prove global face distance or attribute another kernel's error.

The corpus report becomes 20 checked exports, 22 refusal records and the
unchanged oblique-cylinder failed volume comparison. Three refusal records
remain schema gaps. No recipe, radius, material side, overflow rule, caller
accuracy, work allowance or accepted earlier-answer policy changes.

`unequal-ordered` and `mixed-ordered` remain checked but intentionally use
fewer UV fit spans and different STEP bytes. Their 3D curves, topology,
material and reported maximum bounds remain unchanged. Other 37 executable
corpus outputs retain strict manifests and STEP bytes against fresh main.
This is not a latency claim.

## Checks and limits

Rational-frame tests reject the ineligible conics and cover both normal
directions. Jet tests compare the generic inverse's derivatives against the
affine certificate. The standard 1M API export now includes `sphere-dome`.
The stricter local-fitting-budget test still refuses without touching output,
unchanged. Existing precision and no-partial-output tests remain.

Four extra controls pass raw/reference validity, ten witnesses and their
independent washer integrals at 1M: a larger ball, a lower cap, an intermediate
radius and a translated quarter-turn. Fresh main refuses those same inputs
at 1M. The other ten tested parameter/pose controls remain explicit refusals.

Parameter coverage is not universal. Smaller radii, larger spheres and the
tested general placements can still exhaust 1M. The tested spherical
horn/spindle scenes still refuse STEP. There is no new spherical constructor,
arbitrary-trim route, general-pose promise or complete-band bypass.
