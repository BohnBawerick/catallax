## 2026-10-27: connected print output for lattices

Pieces: 5 (catallax-field printing)
Change: adds `printing::Mesh::component_count` and `printing::mesh_plan_connected`, which runs `mesh_plan` and refuses a closed mesh with more than one component as `field.print-disconnected`. `mesh_plan` is unchanged.
Why: lattice export must be one connected closed solid; separate struts must be refused, not published.
Migration: none for existing callers. A lattice export entry calls `mesh_plan_connected` instead of `mesh_plan`.
