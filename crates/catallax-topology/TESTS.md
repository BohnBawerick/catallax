# Piece 3 blind tests

The suite follows `docs/pieces/03-topology.md`. Expected geometry comes from its tables and exact rational formulas, not from a topology implementation. Production files and public signatures are unchanged.

## Run

All 39 blind contract tests were enabled after comparison. These historical results do not replace running the current suite. This publication changes no mathematical assertions or enabled test bodies.

```sh
cargo test -p catallax-topology --release -- --test-threads=1
```

The evaluation command includes ignored tests explicitly, so later exclusions cannot disappear from a comparison:

```sh
cargo test -p catallax-topology --release -- --include-ignored --test-threads=1
```

Compile tests and benchmarks against the public interface:

```sh
cargo test -p catallax-topology --release --all-targets --no-run
```

The eight-thread determinism test launches its own workers even with `--test-threads=1`. Generated families run 11,487 builder-and-check questions with fixed seeds: 7,700 posed/scaled questions, 427 powers-of-ten questions, and 3,360 determinism questions. Each posed/scaled question checks exact fields, counts, orientation, Euler characteristic, construction facts and a fault-free result. Determinism compares structured records from two sequential runs and eight threaded runs, including handles, geometry, entity fields and fact order. No assertion uses elapsed time.

## Coverage

| Contract cases | Tests | Expected answers |
| --- | --- | --- |
| 1–5 | `builders.rs` | All primitive faces in specified order, exact surfaces and curves, point definitions, directed cycles up to rotation, coedge senses, counts and empty marks |
| 6 | `builders.rs` | Empty body, empty queries/maps, no arena growth except one body |
| 7 | `builders.rs` | Exact growth of all eleven entity/geometry arenas, repeated construction without reuse |
| 8–9 | `builders.rs` | Named quarter-turn coordinates, rational poses, independently evaluated implicit signs and block left-hand rule |
| 10 | `common/assertions.rs`, called for every checked builder | Required incidence and parallel facts with construction provenance; the two forbidden block relations remain false |
| 11 | `builders.rs` | All named refusal codes, validation precedence, zero/negative dimensions, unchanged arenas, facts and existing entities on refusal |
| 12 | `builders.rs` | Every primitive at all scales `10^k`, `k = -30..=30`, with exact scaled fields and invariant directions/slopes |
| 13–14 | `bodies.rs` and builder assertions | Creation-order back-links, orphan ownership, first-added-face ownership, first-sight query order, reversed coedge endpoints, deduplication and apex vertex |
| 15–16 | `bodies.rs` | Two hand-built blocks as separate lumps, line contact with distinct topology, both equal and shared line definitions |
| 17 | `bodies.rs` | Void shell with reversed face senses, reversed coedge order and senses; Euler 2 per shell |
| 18–19 | `bodies.rs` | Equator-split sphere and meridian-cut torus, including one edge used twice by one face |
| 20 | `bodies.rs` | Block with cylindrical through-hole, opposite ring uses, exact kind counts and Euler 0 |
| 21–32 | `faults.rs` | Every planted fault checks membership of its exact specified code; repeated checks also compare deterministic fault order |
| 33 | `builders.rs` | Equal handle/entity/geometry/fact sequences across separate models and eight threads |

Additional boundary cases cover empty shells and coedge loops, both one-ended edge variants, a closed curve with one vertex, invalid one-vertex line edges, unbounded cylinders/cones, ring coedges mixed in a loop, repeated coedges, vertices or edges shared across shells, missing endpoint/vertex-loop incidence, and distinct handles at equal coordinates. A deliberate false-fact fixture verifies that `check` trusts recorded incidence rather than making geometric decisions. Fault tests require specified codes but allow additional faults, as the contract requires.

`geometry.rs` also covers intersection curves whose closedness is unknown, meeting-point vertices without explicit coordinates, incidence resolved through coincident facts, a third distinct coedge use, and marks preserved across shared bodies and model clones.

`common/mod.rs` contains fixture-local geometry and boundary descriptors. Hand-built sound bodies do not call topology builders. Arena probes append one sentinel per arena to a clone and inspect the returned indices; they neither change the original model nor require private accessors. Point comparisons match `PointDef::Explicit` and compare its rational point because `PointDef` itself has no equality implementation.

## Speed measurements

The benchmark is opt-in, following the earlier pieces' convention. Without the environment variable it prints the activation command and exits before calling `Model::new`.

```sh
CATALLAX_TOPOLOGY_RUN_BENCH=1 cargo bench -p catallax-topology --bench speed
```

It prints nanoseconds per call and sample counts for:

- Every primitive, including both frustum directions and the apex cone, at identity and a fixed-seed rational pose. Each sample starts with a fresh model; model creation and destruction are outside the measured interval.
- `check`, `counts`, `faces_of_body`, `edges_of_body` and `vertices_of_body` on hand-built height-one prisms with `n = 3, 10, 100, 1,000, 10,000` polygon vertices. Output includes face count.
- `coedges_of_edge` and `owner_of_coedge` across the largest prism's edges.

The prism has `n + 2` faces, `3n` edges and `2n` vertices, so `V - E + F = 2`. Increasing samples of `y = x²`, closed by the last-to-first chord, form a counterclockwise convex polygon. The bottom cycle is reversed, the top cycle is forward, and a side's outward normal is its bottom edge direction crossed with positive z. Fixture construction and initial correctness checks are excluded from query timings. Timing never gates a test or benchmark.

## Verification boundary

The release evaluation passed all 39 blind contract tests and all 10 implementation-private unit tests. The 11,487 generated builder-and-check questions completed without an unexpected refusal, and the planted-fault cases found their required codes. Every blind test now runs in default CI.

The opt-in benchmark's enabled path also completed, including the 10,002-face prism. The comparison report records all timings. Those observations establish this run's results, not a general proof of correctness or a speed threshold.
