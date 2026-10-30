## 2026-10-30: incremental exact builds without eager display

Pieces: catallax-app, SolidDraft Catallax adapter

Change: Added `build_incremental_without_display` and `EvaluatedApplicationBody`. The entry returns the checked exact checkpoint and picking names without a certified mesh. `EvaluatedApplicationBody::display` requests that mesh later with the original marks and refusal codes. Existing eager entries and their return types are unchanged.

Why: Tool bodies and intermediate rebuild results do not need certified display meshes before a caller requests them.

Migration: Existing callers need no changes. The SolidDraft companion can opt into the new entry, retain its checkpoint and names, and cache the result of `display` only when certified display data is needed. Both incremental entries accept the same checkpoints. Each new build or display call spends an active caller allowance unchanged, or uses its own 20M default application guard when no allowance is active. Default exhaustion reports `application.work-limit`. A later display refusal leaves the exact checkpoint usable. The companion owns the kernel pin, display cache, scenario work guards and hosted click-time acceptance; this kernel change does not remove any budget exception.
