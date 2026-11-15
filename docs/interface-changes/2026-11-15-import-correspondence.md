# Verified source replacement correspondence

`catallax_import::correspondence::compare(old, new, dec)` compares two completed
`Imported` results and returns `catallax-import-correspondence/1`. Both source
digests come from producer-owned snapshots of the actual imported bytes. The
report contains verified old/new face and physical edge addresses, unresolved
old entries with reasons and ambiguity candidates, and unmatched new entries.
Each address includes kind, source solid and source entity key.

A verified pair has unique exact equality of its complete accepted definition
in both imports. Face equality includes the surface, orientation, every finite
loop, boundary carrier, endpoint and sense. Physical edges retain their carrier,
endpoints and natural sense; coedge direction is excluded and repeated uses are
deduplicated. Entity numbers, labels, source order and geometric proximity do
not establish correspondence. A reused number on changed geometry remains
unresolved. Duplicate geometry in either file prevents an automatic match.

The importer retains the actual planar, seam and cylinder derivative definitions
after its recorded repairs. These accepted definitions are compared rather than
the inconsistent raw source. Solids left open by declined decisions remain
unresolved. The strict builder's internal spline repair currently returns marks
without its resulting finite definition; such marked solids also remain
unresolved. No raw source definition is substituted for that missing result.

Equality is deliberately conservative. Equivalent reparameterizations, loop
rotations, support-frame changes, splits, merges and changed geometry can remain
unresolved. `changed_or_removed` means that no equal verified new definition was
found, including when a new definition is unavailable. `unmatched_new` likewise
does not claim an addition. All old face and physical-edge addresses appear
exactly once in either `verified` or `unresolved`. Refusal returns no partial
report. The existing active allowance applies, with 20M as the default when no
request is active. Complete canonical representations are charged by byte size,
including rational and spline data, and compared directly rather than by hash.

`Imported` now contains private producer-owned comparison data. Callers continue
obtaining it through `apply`, `apply_source` or `rebuild`; struct literals must
use these constructors instead. Editing its public repair record cannot alter
the retained definitions or actual digests. Serialized correspondence is a report
that callers can edit, not authority to bypass import or document-rebuild checks.
The existing editorial `StepReplacement` API retains its contract. Consumers
must bind a correspondence to both actual sources and rebuild the complete
replacement before atomically updating a document.

## Approach

Before implementing the comparator, inspected OCCT's
[`BRepTools_ReShape`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingData/TKBRep/BRepTools/BRepTools_ReShape.cxx)
and the [P04 SMLib port/consumer map](../phase1/P04/PORT-CONSUMERS.md).
OCCT records an explicit old/new replacement relation and reapplies it to shape
structure. That relation does not itself discover equivalent entities in two
unrelated files. This step follows explicit correspondence at the consumer
boundary and uses complete-definition equality to generate a conservative
relation. It needs no tracer, projection or spline intersection port.

This is the correspondence step of milestone 1.2. Numerical import families,
the full census rerun and SolidDraft's replacement UI retain their separate
acceptance gates.
