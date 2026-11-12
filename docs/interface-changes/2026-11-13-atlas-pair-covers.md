# Atlas covers in the exact patch-pair engine

Change: `RationalAtlas::partition` takes `[bool; 2]` periodic flags instead of a u-only `bool`. Pass `[true, false]` for existing u-periodic sources and `[false, false]` for nonperiodic sources. `periodic_v` exposes the second flag. Both directions require the same exact position and G1 row proof; the flags never close a gap.

`RationalAtlas::patch_meeting(&RationalAtlas, &mut Decider)` uses the existing Yang and Jia cell-pair search over all retained charts. Per-span physical-border ownership and exact periodic identifications reach its exclusion and branch-join stages. The result retains the first atlas for evaluating certified endpoints. `PairArc` adds `first_chart` and `second_chart`; its span indices are local to those charts. Single-patch answers keep chart index zero and their existing span indices.

The join stage identifies fixed periodic row coordinates through the atlas's proof, then uses the existing exact root uniqueness certificates. It does not tolerance-merge endpoints or replace the cell certificates with traced candidates. Original root parameters and chart identities remain available in the result.

Migration: all repository constructor calls use the two flags. Existing `RationalPatch::patch_meeting`, offset meetings, allowances and tolerances are unchanged. No new refusal code, Boolean curve form, exchange schema or STEP route is introduced. See [the measured atlas proof](../M10-SPLINE-SPLINE.md#joined-atlas-inputs).
