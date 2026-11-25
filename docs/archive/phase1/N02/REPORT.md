Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# N02: one numerical trim and endpoint-sewing family

P2 correction to PR #605 at `fec7e989c4ae0c2a9b2f96ed1835e036eaa0a7bc`:
numerical admission now rejects outer/hole loops sharing walked endpoint IDs
after sewing. The review's scaled two-join witness was executed and confirmed:
before repair it accepted a shell with Euler characteristic -2; after repair it
refuses without mutating input. See [before/after evidence](../../../../evidence/phase1/N02/P2/VALIDATION.md).

The public square-tube fixture builds with two-loop top and bottom faces. A split
endpoint 0.005 mm away sews under the explicitly local 0.01 mm tolerance; a
0.02 mm gap refuses. Crossing/touching/outside holes, reversed orientation and
self-crossing trims refuse without changing the input. Results use N01's separate
numerical type, including zero-gap/zero-residual controls.

## Baseline, scope and prerequisites

Repository: `BohnBawerick/catallax`. Remote main was read, fetched and frozen to
`f26c898c94582dd640c754b210cca796071416d4` on 2026-11-14 UTC. N01
[PR #603](https://github.com/BohnBawerick/catallax/pull/603) was verified merged at
that SHA. P02 [PR #599](https://github.com/BohnBawerick/catallax/pull/599) was
verified merged at `416f4b46780a7cba576a44395da343f0e1be137d`, an ancestor of the
baseline. The isolated branch is `dot/n02-numerical-trim-sewing`.

Read N01's merged public producer/interface records and numerical implementation,
the [P02 baseline and limitations](../P02/README.md), repository CI commands and
available checkout guidance. No applicable `AGENTS.md` or `.agents/skills/SKILL.md`
was found. The verify skill was read for context; no CLI/export qualification
workflow was invoked. No standalone N02 card file was found in the tracked plan;
The previous roadmap supplies the numerical
modelling requirement supplies context, not authority to widen the result API.

The open-PR check found #604's separate face-mounted extrusion work and no N02
branch/PR. This change touches only numerical topology admission/sewing, one
focused test and N02 documentation/evidence. P04, B01, B02, P02 and their branches
and artifacts are preserved. Neither #604 nor its review is part of this work.

## Implementation and contract

[The interface record](../../interface-changes/2026-11-14-numerical-trim-sewing.md)
describes the bounded API. [Numerical admission](../../../../crates/catallax-topology/src/numerical.rs)
now stores both loops, checks their opposite orientations, rejects inclusive
boundary intersections and checks strict hole containment using rational
parameter-space predicates. A separate check also rejects shared walked endpoint
IDs across the two loops, including IDs unified by sewing. It retains N01's outer-loop checks and patch-domain
containment. A duplicate outer loop still refuses as `numerical.invalid-trim`.
Only one hole is admitted; no general region/Boolean proof machinery is added.

`numerical::sew(&definition, &[[keeper, removed]], options, &mut decider)` measures
and records endpoint movement, unifies only the explicit disjoint vertex pairs,
then runs N01 admission on a private clone. The input is borrowed and remains
unchanged. Edge carriers must already be shared and explicitly identified;
independent edge or face matching, mesh sewing and automatic gap search are not
implemented. Every edge and every coedge of both loops retains N01's actual
endpoint and nine-location distance checks. These are sampled checks, not bounds
between samples or over the body.

The successful type is `NumericalBody` with `NumericalEvidence`, never a new exact
or bounded result. The exact checker still reports missing incidence. No
`Answer`, `BoundedEvidence`, N01 adapter, PRD, budget, global tolerance, exporter or
saved answer changes. N01's existing public summary can consume the unchanged
type; this card does not add a wire adapter or serialize an exact result.

## Public accept/refuse examples

[The focused Rust fixture](../../../../crates/catallax-topology/tests/numerical_trim.rs)
is independently constructed from integer-coordinate rectangles. It has outer
cross-section [0,4]², hole [1,3]² and height 1 mm. The retained shell has 16 vertices,
24 edges, 10 faces and 12 loops. Top material normal is +z; inner-wall normals
point toward the hole. The two top loops have signed areas +16 and -4 mm²,
respectively, so the independently known material footprint is 12 mm². These
orientation and area witnesses do not claim a general numerical mass or material
query API.

A source endpoint is split into a separate vertex used by one edge. Unsewn
admission refuses structural closure. Explicit sewing restores the shared graph.
The tests compare the borrowed input's complete debug representation before and
after success and refusal. Nonzero movement retains the original point, keeper
position, removed key, upper movement distance and mapping check. Zero-gap sewing
records the mapping but no fictitious movement repair.

| Example | Observed outcome |
| --- | --- |
| Two-loop square tube | Numerical body, expected counts and material orientation |
| Split endpoint, gap 0 mm | Numerical body; no vertex movement repair |
| Split endpoint, gap 0.005 mm; local tolerance 0.01 mm | Numerical body; one measured vertex movement |
| Gap 0.02 mm with the same local tolerance | `numerical.repair-tolerance-exceeded`; no returned body; input unchanged |
| Crossing, touching, outside, reversed or self-crossing hole | `numerical.invalid-trim`; no returned body; input unchanged |
| Reused, identical or invalid vertex indices | `numerical.invalid-definition` |
| Request allowance 1 | Refusal and exhausted accounting; no nested/reset allowance |

The near-gap [retained definition and evidence](../../../../evidence/phase1/N02/P2/near-gap.txt)
are actual repaired Rust debug records captured by the test's opt-in evidence sink, not a
new stable wire schema. Refusal records in the same directory are actual debug
outputs. [Validation](../../../../evidence/phase1/N02/VALIDATION.md) separates executed
checks from prior evidence and remaining gates.

## Evidence boundary

P02's frozen M10/M8 refusals motivate keeping approximation explicit. Those
requests were not rerun and none is newly claimed successful by this planar
family. N01's existing sampling/result rules are reused; there is no C++ oracle,
private corpus, companion repository, broad duplicate campaign or new global
proof layer. The implementation is serial; no worker comparison or production
parallelism is claimed. The caller's allowance spans endpoint sewing and the
whole rebuild. Numerical STEP export, whole-body deviation/separation,
self-intersection freedom in 3D, general patch sewing and application edit chains
remain outside this evidence.
