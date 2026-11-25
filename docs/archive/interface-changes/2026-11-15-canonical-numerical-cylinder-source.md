Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Canonical numerical cylinder sources

`catallax_boolean::numerical::NumericalCylinderSource::new` constructs a finite
source independently of a cone cut. It accepts `&Cylinder`, an increasing
`[Rational;2]` interval, `NumericalOptions` and `&mut Decider`. The original
cylinder is retained. Coordinates measure signed millimetres along its unit
axis from its original origin, including algebraic compiled axes.

The descriptor contains the full numerical B-spline source, all four embedded
rational patches, source/control evidence and measured source/chart work. Quarter
IDs are 0 through 3. Keys are descriptor-local `cylinder/quarter-i`; callers
supply their own assembly identities. Each patch retains native u=[i,i+1] and
v=[0,1]. No cone face key enters the source or quarter identities.

`axis_parameters_mm()` records signed millimetres at native v=0 and v=1.
`signed_mm_at_v` applies this affine mapping and refuses v outside [0,1].
The original origin, axis and seam direction remain available through `original()`.
Descriptor and quarter fields are private. Read-only `source()`, `quarters()`,
`evidence()`, parameter and work accessors keep the constructed definitions,
interval and actual checks together. Quarter accessors expose `id()`, `key()`,
`surface()` and `axis_parameters_mm()` for assembly without mutating the source. The same carrier
and interval produce the same full source as the finite cone/cylinder query.

`catallax_intersect::numerical_native_cylinder_source` supplies the lower-level
retained source. Both factories reuse the merged native `Frame`, circle control
net, Cartesian rounding and analytic carrier checks. No second parameterization,
UV fitting or curve reparameterization is introduced. Quarter admission uses the
existing retained rational chart embedding check. The prescribed tolerance and
all source checks stay numerical. Exact embedding of a retained patch does not
prove exact incidence on the original analytic cylinder.

Every call uses the caller's active work allowance and deadline. Standalone calls
keep the existing 20M allowance. Invalid radii, frame directions, collapsed or
reversed intervals and nonpositive tolerances refuse. Source controls and
quarter patches do not by themselves form Boolean topology, an admitted body,
a STEP result or a complete body deviation bound. Hole assembly retains the
Entry/Wall provenance and performs its own complete shell admission/export.
