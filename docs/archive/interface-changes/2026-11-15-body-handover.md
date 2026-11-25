Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-15: transfer built bodies between kernel instances

`catallax_app::handover` exports four functions:

```rust
pub fn write_body(body: &catallax_features::Built, decider: &mut Decider)
    -> Result<Vec<u8>, Refusal>;
pub fn read_body(bytes: &[u8], decider: &mut Decider)
    -> Result<catallax_features::Built, Refusal>;
pub fn write_numerical(body: &NumericalBody, decider: &mut Decider)
    -> Result<Vec<u8>, Refusal>;
pub fn read_numerical(bytes: &[u8], decider: &mut Decider)
    -> Result<NumericalBody, Refusal>;
```

A worker can build and mesh a recipe, send the body bytes and its display mesh,
then let the page read the body into independent stores. Reading does not build
the recipe, mesh the result, or import STEP. The existing desktop proof completion
publishes an `Arc<IncrementalApplicationBody>` through revision authority in
`catallax_app::lifecycle`. This format supplies the owned state needed across a
worker boundary where that pointer cannot cross. The caller still rejects stale
edit completions.

The exact entry retains `Built`, including intermediate bodies, recipe, fillet
and chamfer definitions, removed-face history, field stand-ins, source solid
keys and import repair records. The model retains every topology and geometry
arena, persistent names, construction/proof facts, chart covers and marks.
Numerical transfer retains its separate type, original/current definitions,
tolerance, checks, repair evidence and prior states. It grants no extra proof or
export authority. Picking handles keep their indices and belong to the new model.
Incremental evaluation caches and display meshes travel separately.

The design follows OCCT's [BinTools shape-set serialization](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingData/TKBRep/BinTools/BinTools_ShapeSet.cxx):
number shared records, write definitions once, and restore their links. Catallax
also needs exact expression/root sharing and retained evidence. A geometry-only
exchange would lose those. [BRepTools](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingData/TKBRep/BRepTools/BRepTools.cxx)
uses the same shape-set pattern with textual storage.

Format version 1 has an eight-byte `CXBODY\0\0` magic, a little-endian `u32`
version, the root kind, a type-name table, numbered records in dependency order,
the root record, and a SHA-256 digest of all preceding bytes. Lengths use
little-endian `u64`; handles and shared record numbers use `u32`. Import checks
platform widths. Arbitrary-precision integers use signed little-endian bytes.
Rationals retain their canonical numerator and positive denominator, with no
float conversion. Floating bounds retain their IEEE bits. Both graph walks are
iterative, so expression depth does not become decoder stack depth.

`handover.version-mismatch` refuses another format version.
`handover.invalid-data` refuses damaged envelopes, truncation, wrong root kinds,
unknown records/variants, invalid lengths/references, noncanonical rationals,
duplicate keys and handles outside their arenas. The digest detects damage; it
does not authenticate an untrusted producer or prove its geometry. Consumers
should accept transfers from their own kernel build. Changes to record fields,
variant order, or type names require a format-version bump.

Transfer charges the active request allowance through `Decider::charge_work`.
The encoder charges before variable-length walks, integer conversion, buffer
growth, record copies and hashing. Metadata JSON streams through the same
meter, with checks before string escaping and collection visits. It uses the
existing lengths and integer bit lengths, without a size-estimation pass.
Without an active allowance it uses the application default. Input scanning,
decoding and rational restoration consume work; import does not obtain a fresh
allowance inside an existing request. Journals, read observers and thread fork
bookkeeping reset because they belong to the old execution. Retained geometric
definitions, evidence and numeric enclosures survive. Existing public fact/mark
codes use static strings; import interns their vocabulary once with a 64 KiB
process limit to bound storage for custom codes.

The change adds codec derives to retained types and a small arena/handle hook.
It does not change the STEP curved numerical exporter or analytic import dispatch.
