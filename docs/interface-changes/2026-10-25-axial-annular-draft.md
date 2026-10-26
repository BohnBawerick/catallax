## 2026-10-25: bounded axial annular draft cells

Pieces: recipe features, independent field evaluator
Change: `draft::annular::{AnnularWall, AnnularCell}` describes one coaxial annulus with two perpendicular caps. Draft now accepts either or both complete cylindrical walls when their affine radii remain strictly positive and separated over the finite height; source names and the hole survive.
Why: the previous single-cylinder cell guard returned `draft.unsupported-arrangement` for this bounded hollow section.
Migration: no recipe syntax or caller changes. Consumers retain inner-wall reversed material sense, existing retention rules and finite annular cap domains. Crossings, folds and general curved arrangements retain coverage refusals; no tolerance, allowance, cell cap or display interface changes.
