# Circular first-hit caps on stock spheres

The public fixtures contain one circular disk of radius 3, a stock sphere of
radius 5 centred 10 mm above the original sketch, and each of the two limit
modes. Entry starts outside; exit starts at the sphere centre. These are
synthetic definitions. The initial clean source base was
`e9d362f7eb8fa41597ca5f16bf5a40b208fa0a09`, containing merged #117. Its preserved
verifier SHA-256 is
`322a1801958de0b28b489d639af1ee8f08e818971409c68fda741a9ffc4928bf`.
The delivery integrates `49ac2af8a6ae99b6a1c476ad7981918048e5c689`;
the PR records exact final source/runtime identities and qualification results.
Earlier proof bases remain historical observations.

## Defect ledger

| Unit | Reproduction | Result and proof |
| --- | --- | --- |
| Circular profile coverage | Four fresh `recipe` requests and the same four `field` requests on the initial base return `features.unsupported-curved-limit`. No STEP or OBJ is published. Five new enabled Rust tests also fail before implementation. | A distinct whole-disk radial proof replaces the polygon-only guard for coaxial stock spheres. Exact discriminant, start-state and positive finite root bounds cover the axis, every radial coordinate and the closed rim. The interface fragment states the inequalities. |
| Disk domain | A square enclosure contains points outside the circle, and a radius-4 disk on a radius-5 sphere has a valid first hit even though its square enclosure crosses the silhouette. | The field and B-rep retain the disk itself: one cylindrical side, a plane and a spherical cap. Rational field-oracle checks include square corners, start extensions, poles, rims and the remote spherical root. No circle is polygonised. |
| Invalid circle radius | A new zero/negative-radius regression reproduces a different refusal from direct cap compilation than from the sketch consumer. | Reuse `sketch.non-positive-radius` in all three entry paths; no new refusal code is introduced. |
| Spherical ring chart | The first implementation's transverse spherical chart imports valid solids, but all 24 enabled independent native volume cases fail the unchanged `1e-8 mm³` assertion, by roughly `5e-5 mm³`. The failing executable, exports and logs are retained separately. | The new circular builder aligns its sphere chart with the sweep axis, so the circular rim is a constant-latitude trim. All 24 original read-back cases then pass. The exact sphere, circle, sense, finite point set and numerical assertion are unchanged. No exporter or repair pass is changed. |

The first test-development compile errors and two setup corrections are also
retained in external logs: the volume-width assertion now explicitly requests
the tighter `1e-11 mm³` absolute accuracy it checks, and the missing-reference
fixture uses a structurally valid shared-edge query. No numerical assertion was
relaxed and no failing test was disabled.

## Geometry and independent controls

The [interface fragment](interface-changes/2026-10-26-circular-sphere-first-hit.md)
specifies the full admitted domain, uniform root bounds, entry/exit sets and
remaining gaps. Both evaluators compile operation definitions independently.
The field never reads B-rep topology to obtain geometry, and its direct query
uses implicit support signs and the ray derivative rather than the height
evaluator used by the exact builder.

Focused geometry tests cover both limit modes and branches, 8,100 rational
field-oracle queries, exact membership and origin unions at the pole and rim,
positive and negative raw direction scaling, proper noncardinal world frames,
uniform length scaling, and radius/start-height edits with physical-edge names.
Requested B-rep volume enclosures contain the independent disk integral
`2 pi (R³-(R²-r²)^(3/2))/3`, combined with the start-height cylinder term.
The width assertion is `1e-11 mm³`. Cases include rational and irrational rim
heights and a disk whose enclosing square crosses the silhouette.

The real application tests request default display, picking maps and an
incremental target-radius rebuild. Both ring edges and all three extrusion-owned
faces retain their names. Caller budget exhaustion stays a typed refusal.
Independent STEP tests use raw OpenCASCADE reads with no topology healing,
and check validity, three analytic faces, material samples, axis/rim controls
and the analytic volume in both modes, both branches, radius/sphere/start
edits, reversal and a proper rigid placement. Import validity and numerical
properties remain distinct from certified Catallax property bounds.

Native field display, requested field mass, sampled comparison and complete
application timings have separate outcomes in the external run evidence and
PR. A successful independent field construction is not a mass or whole-boundary
distance certificate. Shared-environment timings are diagnostic; this change
makes no speedup claim.

## Commands and preservation

Reproduce the public circle fixtures through the checkout verification driver:

```sh
cargo build --release -p catallax-cli --bins -j 2
.claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary "$BINARY" --json
.claude/skills/verify-catallax/scripts/control-catallax.py run crates/catallax-features/tests/fixtures/circular-sphere-limits --source recipe --binary "$BINARY" --out "$EVIDENCE/recipe" --fresh --json
.claude/skills/verify-catallax/scripts/control-catallax.py run crates/catallax-features/tests/fixtures/circular-sphere-limits --source field --binary "$BINARY" --out "$EVIDENCE/field" --fresh --json
.claude/skills/verify-catallax/scripts/control-catallax.py run crates/catallax-features/tests/fixtures/circular-sphere-limits --source compare --binary "$BINARY" --out "$EVIDENCE/compare" --fresh --json
cargo test --release -p catallax-features --test circular_sphere_end_conditions -j 2 -- --test-threads=2
cargo test --release -p catallax-app --test circular_sphere_limits -j 2 -- --test-threads=2
python3 -B -m unittest discover -s tools/step-regressions -p 'test_circular_sphere_limits.py' -v
```

The final qualification uses the full release workspace suite, doctests,
format/clippy checks, verification-driver tests, all independent STEP
regressions and CI on the recorded final head. Public polygon, rectangle and
annular-torus controls retain strict manifests and STEP/OBJ bytes under the
existing comparison exclusions. Fresh worker-count checks retain independent
output directories and the original request allowances.

No aggregate capability file or historical interface fragment is edited.
Production changes stay in the curved-cap dispatcher, new circle helper and
builder, the two circular construction dispatches/private cell constructor,
and field cap membership/origins. Polygon algorithms, field root hooks, native
display, print routes and exporters are unchanged. The 1M, 20M and 100M
allowances, 0.01 mm field accuracy, 2M native display-cell limit, independent
vertex limits, mass/print caps and modelling/local tolerances are unchanged.
