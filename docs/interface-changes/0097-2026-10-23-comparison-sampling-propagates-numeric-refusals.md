## 2026-10-23: comparison sampling propagates numeric refusals

Pieces: comparison, verification CLI
Change: reverse comparison sampling uses the checked field-query API and records numeric refusals as unresolved observations instead of unwrapping them. Feature sampling stops on the first such failure. Request-wide exhaustion still returns `numbers.operation-budget` from the verification CLI.
Why: a valid field face handle does not guarantee that a trim or orientation proof fits the caller's arithmetic allowance.
Migration: no interface signature, tolerance, work allowance or successful answer changes. Callers receive the existing unresolved or numeric-budget outcome instead of `cli.internal-error` when a comparison exhausts its proof allowance.

