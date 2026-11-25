Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# N03: checked approximate STEP for the N02 planar family

The numerical two-loop tube and its N02 endpoint-repaired control export as
explicitly approximate STEP. Both actual disk files were reread and checked,
then independently imported raw in OCCT 7.8 without healing and cut successfully.
Source admission, emitted decimal-value admission, byte identity and reference
geometry checks have separate records. None promotes numerical evidence to an
exact or bounded answer.

## Baseline, prerequisites and ownership

Actual origin/main was observed, fetched and frozen on 2026-11-15 UTC to
`ecad1812b339bb759d83ba8f057c927d69d393ec`. Git ancestry checks confirm N02 #605
merge `3e2a929e709c7916bfa8519e6cf5f7179cbc7cd5` and the supplied E01 #606 merge
`ecad1812b339bb759d83ba8f057c927d69d393ec` are ancestors. Branch:
`phase1/N03-approximate-step`; no competing local/remote branch or open N03 PR was
found before work. Previous branches, checkouts and evidence are preserved.

Read repository CI/verification guidance, N01/N02's numerical contract and
public N02 fixture, the STEP writer/seam contracts and supported raw reference
checker. No applicable AGENTS.md or .agents/skills/SKILL.md was found. The
verification skill was read for context; no broad CLI campaign was invoked.
No standalone tracked N03 card file was found; the explicit N03 assignment and
ROADMAP 1.1 provide this bounded scope.

Provenance discrepancy: the supplied #606 merge is titled “Round mixed fillet
vertices while retaining unselected edges”; no separate E01 pcurve-planning
record was found at this baseline. This planar family does not need curved
pcurves. It reuses only existing graph/decimal/AP242 framing code and the
compatible planar seam check, with exact-model uncertainty claims removed from
the numerical file. It inherits no exact/pcurve proof authority.

Only the new STEP numerical module, its declaration, directly scoped tests,
reference script and N03 documentation/evidence change. No bored_shell,
planar_shell or other E03 geometry code, topology implementation, exact writer,
E01 feature, N01 adapter, PRD, tolerance/cap or saved answer is changed. Writing
`CLOSED_SHELL` serializes the already retained shell; it is not a new shell
construction algorithm or inferred topology.

## Result and file contract

[Public contract](../../interface-changes/2026-11-15-numerical-step-export.md).
[Implementation](../../../../crates/catallax-step/src/numerical.rs).
`numerical::build` takes N01/N02's immutable numerical body and returns a separate
`NumericalStep`. No `Answer`, `BoundedEvidence` or exact `StepReport` is produced.
The generating definition, source numerical checks/tolerance/repair history,
emitted decimal-value checks and scoped scalar-rounding bound are retained.
N01's generic-constructor trust caveat is not solved by this separate producer.

The first family is rational axis-aligned planar faces and finite line carriers,
with explicit shared vertex/edge identities and outer/inner loop orientation.
Merely tolerated source gaps refuse; the exporter performs no healing. Unsupported
carriers refuse rather than becoming a mesh or inferred plane. The decimal model
is separately admitted through N01/N02 with the same IDs; collapsed edges and
invalid loop relationships after rounding refuse. The actual file must match the
checked bytes; this is distinct from source-body validation and from raw OCCT
shape validation.

The file header states `catallax numerical approximate`, source name, source
local tolerance, check/repair counts, scoped scalar error and
`complete_body_bound=unestablished`. It omits the exact writer's global uncertainty
entity/context, so the 0.01 mm sampled source tolerance is not a declared complete
file/body bound. Zero residual and zero scalar rounding still remain numerical.

## Public inputs and executed checks

[Fixture](../../../../crates/catallax-step/tests/common/n03_tube.rs) copies N02's
public tube construction without geometry changes; N02 source is untouched.
The tube has outer [0,4]², through-hole [1,3]² and height 1 mm. A second control
splits an endpoint by 0.005 mm and uses N02's explicit sewing at the unchanged
local 0.01 mm tolerance. That result retains one actual source vertex movement.
“Raw import” below refers only to import of its exported file, not to a claim
that this repaired source had no prior repair.

[Focused Rust tests](../../../../crates/catallax-step/tests/numerical_step.rs) pass:

- Both bodies export, are physically written, reread and byte/seam checked; each
  has 10 faces, 24 edges and 16 vertices. Source evidence and emitted-value
  evidence are retained separately. Emitted repair histories are empty.
- Tolerated source gap refuses `step.numerical-gap` without mutation/healing.
- A tiny translated body that collapses at the requested decimal precision
  refuses during emitted-value admission, despite passing source admission.
- A flat spline carrier refuses `step.numerical-unsupported`; no planar shortcut.
- Invalid options and allowance-1 controls refuse. One combined fixture,
  source-sew and export request succeeds at 349,106 work under unchanged 20M.
- Exact legacy unit block retains all 6,360 bytes and 2,400 work units; both legacy
  and required seam options preserve the same bytes. The existing independent
  exact block contract test also passes unchanged.
- Compile-fail test rejects implicit conversion from `NumericalStep` to an exact
  `Answer<StepReport>`. This does not prevent explicit misuse of generic wrappers.

[Raw reference script](../../../../tools/step-regressions/check_numerical_step.py)
reads the exact files below with OCCT/OCP 7.8.1.1, Python 3.12.14. It uses the
supported direct `StepToTopoDS_Builder`, `read.stdsameparameter.mode=0`,
`read.maxprecision.mode=0`, precision/max tolerance 1e-7 mm and explicit identity
millimetre/radian factors. No TransferRoots, ShapeFix, sewing, SameParameter or
file rewrite is used. The older shared helper's two-argument Init call initially
failed at the binding interface; N03's script supplies 7.8's required
`StepData_Factors` to the same builder, without changing the shared checker,
reference kernel or healing policy.

Both imports are raw-valid single solids: 10 faces, 24 edges, 16 vertices, two
2-loop faces, face-hole-adjusted Euler 0 and volume 11.999999999999998 mm³.
Subtracting the reference box with corner (-1,-1,-1) and size (1.5,6,3) removes
the x≤0.5 wall portion. The result is valid, one solid, volume 10 mm³ with the
same topology counts and through-hole. Material probes show (0.25,2,0.5) changes
from inside to outside; (0.75,2,0.5) remains inside; (2,2,0.5) remains outside in
the hole; the new x=0.5 boundary is identified. Input bytes remain unchanged.

| File | SHA-256 |
| --- | --- |
| two-loop.step | `7f3778332cb35917fbe5c4f840a6cd0d97f8e98cb069f207c1c28e528a456738` |
| repaired.step | `715713ed356758392c51639f87ad3852618e7e50e033d7dc46e91f9474e413b7` |
| legacy-block.step | `7148692d5d9b552c6c93484699740718d056ea6654889792f2b1f7d7ee3d3c2f` |

Files, full source/emitted evidence records, numerical summaries and
[raw reference outcomes](../../../../evidence/phase1/N03/raw-reference.json) are
retained under `evidence/phase1/N03/`. Commands and source/artifact hashes are in
[validation](../../../../evidence/phase1/N03/VALIDATION.md).

## Acceptance and limits

| N03 condition | Evidence status |
| --- | --- |
| Public N02 two-loop/control export, approximate category | Established for the two named planar controls |
| Source/tolerance/repair/bound provenance retained separately | Established; complete body bound explicitly unestablished |
| Actual emitted files reread and checked | Established by Rust disk reread plus hash-bound reference checks |
| Raw OCCT import with healing off | Established for the two retained hashes |
| Subsequent cut with relevant topology/material checks | Established for that reference box cut on both files |
| Exact legacy bytes and unchanged policies | Established for retained unit-block bytes/work and independent block contract control; no global replay claim |
| General numerical/curved STEP, nonzero tolerated geometric gaps | Not established; typed refusal in this first writer |
| Whole-body deviation, global self-intersection/separation certificate | Not established by numerical samples or scalar rounding bounds |
| Native application edit chain, universal/multiworker qualification | Not run; no production-parallelism claim |

This is a bounded N03 delivery, not general STEP qualification. No held-out/private
input, personal machine, different reference oracle, donor code, manual hosted CI,
merge, mark-ready or deployment is involved. Remaining broader numerical export
families need their own retained contracts and file/reference checks.
