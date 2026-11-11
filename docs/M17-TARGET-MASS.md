# Target mass and centre: bounded increment

This is not completion of M17. The complete authored P5 still refuses; P6 has
certified but wide field enclosures, not a matching B-rep solid. No existing
measurement, allowance, tolerance, recipe or saved answer is changed.

Baseline: `cdc0a3884c59b2bea538c230e37639bd74cdb9a2`. The new opt-in
`catallax_field::target_mass::measure` compiles the complete native recipe and
uses the merged F7 classifier/grid. It requests only mass and first moments,
not inertia. This avoids the older all-properties route's unsupported classifier
families and inertia requirement without weakening that route's contract.

Density is an explicit positive rational in g/mm³, uniform over all material.
The evidence uses `157/20000` g/mm³ as a control, not an inferred material.
Recipe placements apply before integration; centres are in world millimetres.
Multi-density assemblies, material inference, imported file maps and inertia are
outside this increment. Callers retain rational intervals, not their midpoint.

## Exact-source acceptance matrix

Every row uses the complete unchanged file in `parts/`; source digests are in
[provenance](evidence/m17-target-mass/provenance.txt). All recorded runs use one
100M operation allowance including compile, prepare, jobs and reduction.

| Target / exact source | Field at 64³: mass interval (g) | B-rep agreement at explicit 16³ | Remaining limit |
|---|---|---|---|
| [P1](../parts/P1-bracket.json) | [89709957/655360, 216340191/1310720] | Common certified mass and centre; 44,258,952 work | Field bounds alone remain wide |
| [P2](../parts/P2-collar.json) | [17384139/640000, 25344981/640000] | Common certified mass and centre; 1,329,565 work | Field bounds alone remain wide |
| [P3](../parts/P3-filleted-boss.json) | [43820279577/524288000, 38919115121/262144000] | Common certified mass and centre; 1,351,147 work | Fillet field bounds especially wide |
| [P4](../parts/P4-shelled-tray.json) | [1651797/51200, 266832333/6553600] | Common certified mass and centre; 3,540,192 work | Field bounds alone remain wide |
| [P5](../parts/P5-flanged-housing.json) | `fillet.unsupported-overflow`, 3,367 work | Same compile refusal | M18 authored radius-2 finite transition/material coverage; no other lane incorporated |
| [P6](../parts/P6-lattice-bracket.json) | [106355097/40960000, 266175759/81920000] | `target-mass.brep-definition` | B-rep is an envelope stand-in, not lattice material; field centre uncertainty remains |

Exact centre endpoints, first moments, geometry/material revision hashes and
complete work records are retained in [field64](evidence/m17-target-mass/field64.txt),
[field16](evidence/m17-target-mass/field16.txt) and
[agreed16](evidence/m17-target-mass/agreed16.txt). The records are newly generated
public recipe evidence, not a saved-answer inventory. No finite width is silently
called an exact scalar. No target accuracy threshold is introduced or relaxed.

## Certificate and agreement boundary

The F7 grid partitions the entire compact definition bounds into equal rational
boxes with disjoint interiors. Each occupied cell volume lies in its fraction
enclosure times exact cell volume. Full-measure cells have exact box first moments.
For other cells, multiplying the occupied-volume interval by the signed coordinate
interval encloses its first moment: this remains sound across zero and with
negative coordinates. Sum in canonical cell order, then divide by a strictly
positive total-volume lower bound. Intersect the centre with the global bounds
(a positive-mass barycentre lies in their convex hull). A nonempty result without
positive lower volume refuses; certified zero volume has no centre. Boundaries
of measure zero neither add nor remove volume. No mesh or point sampling enters.

`measure_agreed` builds the **same recipe** independently through B-rep features,
rejects all stand-ins and marked result bodies, and requires unmarked volume and
centroid measurements. Those two exact-definition routes enclose the same quantity;
their rational interval intersection is therefore one justified result. Arbitrary
interval overlap is not used to establish definition equivalence. Disjoint
intervals or conflicting emptiness refuse. Original first moments remain field
enclosures. Existing B-rep and field measurements are not modified. A B-rep refusal
is returned, never hidden behind a successful field result. Field-only measurement
remains available separately. The 16³ agreement configuration is explicit, not
automatic coarsening after a refusal; 64³ plus B-rep is not claimed within 100M.

## Work, workers and validation

The existing F7 fixed 128-cell jobs reserve their allowances once, independent of
workers. The new serial ordered reduction charges 64 units before each cell's
rational moment arithmetic and 32 before division. It does not give each worker
a fresh request allowance. Defaults use all available cores bounded by nonempty
jobs. All target records compare full typed answers/refusals (including undecided
evidence), actual classifier work, charged operation work and exhaustion at
1/2/4/default workers. Four cores were available in this run. Fixed job assignment,
local state and canonical reduction explain any-count numeric invariance; finite
tests alone are not a universal proof. Spawn failures and wall-clock cancellation
are environmental refusals, not deterministic numeric disagreements. Deadlines
are checked by the existing query and before every reduction cell.

Five focused regressions cover signed exact box moments/density, analytic translated
sphere containment, full records at 1/2/3/4/default/99 workers and low budgets,
expired deadlines, and independent exact box agreement. The 9³ control creates
six nonempty canonical jobs, exercising uneven and overprovisioned scheduling.
Tests inspect full records rather than Debug strings. Evidence printing uses Debug
only for persistence, after typed equality assertions.

Runs completed on 2026-11-11 UTC. Reproduce from this branch's committed source:

```sh
cargo test -p catallax-field --release --test target_mass
cargo test -p catallax-field --release --test mass --test assembly
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
cargo build -p catallax-field --release --example target_mass
target/release/examples/target_mass 16 > field16.txt
target/release/examples/target_mass 64 > field64.txt
target/release/examples/target_mass 16 agreed > agreed16.txt
sha256sum target/release/examples/target_mass
```

Focused checks and retained target evidence are local Rust evidence only. Hosted
gates, C++/OCCT, protected replay and full-suite qualification were not run.
