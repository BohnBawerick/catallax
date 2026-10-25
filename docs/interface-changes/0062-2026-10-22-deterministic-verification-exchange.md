## 2026-10-22: deterministic verification exchange

Pieces: numeric accounting and CLI verification
Change: add `catallax-verify stream`, `catallax-request/1`, `catallax-response/1` and `catallax-manifest/3`. Question, recipe/3 and two-kernel comparison sources share a deterministic whole-operation allowance. Add `with_operation_budget` for synchronous request scopes. Manifest/3 defines distinct elliptic-cylinder, parabola and hyperbola tags.
Why: the legacy 50-second deadline lets machine load choose between an answer and a refusal. Verification needs deterministic outcomes and separate latency evidence.
Migration: existing `catallax` and `exact-brep` commands, question/2 and manifest/2 remain unchanged. Adopt the new binary explicitly. Question adapters parse but ignore legacy `deadline_ms`; `max_work` is the new allowance. Do not dispatch numeric work across threads inside a request scope. Parallelize independent requests. See [VERIFY-INTERFACE.md](../VERIFY-INTERFACE.md) for streaming, cache identity and cancellation behavior.

