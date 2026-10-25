## 2026-10-22: recipe format `catallax-recipe/1` and crate `catallax-recipe`

Pieces: recipe (new); the judge interface is unchanged
Change: `docs/RECIPE.md` defines the recipe language, version 1: operation definitions, origin-based face and edge names, references, rebuild failures and the canonical JSON file `catallax-recipe/1`. The new crate `catallax-recipe` holds its data model: `Recipe`, `Step`, `Operation` and their parameter types, `Origin`, `Role`, `EdgeName`, `VertexName`, `OriginPattern`, `FaceRef`, `EdgeRef`, `EdgeSelection`, `Number`, `RecipeError` with codes `recipe.<reason>`, `Recipe::parse`, `Recipe::check` and `Recipe::to_text`. It depends only on piece 1, and nothing depends on it yet.
Why: PRD section 10 step 3. Both kernels and SolidDraft need the same definitions and names before either kernel implements an operation.
Migration: none. No piece, the CLI or the exchange contract changes. The question format stays the judge's input.

