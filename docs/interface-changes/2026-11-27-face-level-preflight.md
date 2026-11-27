## 2026-11-27: exact per-face display preflight

Pieces: catallax-mesh, catallax-app, SolidDraft display
Change: `plan_face_levels` prepares the shared trim samples and face topology once. `FaceLevelPlan::point_count` reports the distinct mesh vertices a logical face references at each level, including its share of a sewn edge. `emit_selected` consumes the plan with one level per logical face; `emit`, `emit_bounded` and `emit_preview` select a uniform level.
Why: A viewport can choose levels before allocating its render buffers. Coarser rounds use points retained by the finer rounds. Preflight runs the geometric proof and is not a cheap estimate.
Migration: Use the planner when a device budget needs counts before emission. Handle a typed refusal if the checked preflight cannot establish a schedule. Existing preview and certified mesh entry points remain available; definition-native first paint is unchanged. The `display --levels` benchmark records count equality, direct certified mesh equality, and preflight work and time by carrier.
