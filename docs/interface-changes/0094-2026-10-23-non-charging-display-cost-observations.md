## 2026-10-23: non-charging display cost observations

Pieces: numbers, field, verification CLI
Change: `current_operation_work` reads accepted request work without charging it; outside a budget it returns zero work and no exhaustion. `FaceId::index` exposes the owning field registry index. `mesh_plan_reported` retains `Diagnostics`, `Cost`, `FaceWork` and `SupportWork` for an executed display prefix, including refusals. Field verification manifests add `catallax-display-diagnostics/1` with per-support cells, maximum depth, decisions and work, and per-face boundary-query costs.
Why: an exhausted shared numeric allowance must retain which face/support work consumed it, without changing the geometry's allowance.
Migration: no existing signature, tolerance, allowance or refusal code changes. Diagnostic reads perform no geometric query and consume no work. Shared support costs include face query costs and must not be added to them. Numeric-budget refusal manifests preserve only these observations, never partial geometry or display certificates. Work fields are zero outside a request budget.

