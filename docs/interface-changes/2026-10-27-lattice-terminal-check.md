## 2026-10-27: opt-in terminal lattice check

Pieces: 1, 5 (catallax-recipe)
Change: `Recipe::check_lattice_terminal` and `RecipeError::LatticeNotTerminal` (code `recipe.lattice-not-terminal`) report the first step that reads a lattice step. `Recipe::check` and `Recipe::parse` are unchanged and still accept a lattice as an input.
Why: an application that makes the lattice its last feature must block CAD after it. The check is opt-in because the public every-operation fixture unites lattices into later steps.
Migration: callers that match `RecipeError` exhaustively add the new variant. An application that treats the lattice as terminal calls the check on the history it is about to build.
