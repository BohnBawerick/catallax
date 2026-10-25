## 2026-10-22: composable geometric error decision

Pieces: 1 to 7, recipe evaluation, field, mesh, compare
Change: [ARCHITECTURE.md](../ARCHITECTURE.md#composable-geometric-error-decision-of-2026-10-22) defines a whole-result geometric certificate, typed local conversions, operation transfer proofs and a strict bound below 1 nm. It reserves `numbers.geometric-gap-unproven` and `numbers.geometric-gap-exceeded` for implementation.
Why: inherited mark logs record local uncertainty without proving its effect over the final solid, a recipe chain or an exported representation.
Migration: implement the decision in the documented piece order on the shared algebraic foundation. Audit every producer and consumer before certifying marked results; keep measurement and display uncertainty separate. This is a documentation contract, with no runtime API or code change. The command-line contract, legacy manifest meanings and work-budget policy remain unchanged. A later implementation must log its concrete API migration; certificate exposure through the CLI needs a separate versioned exchange decision.

