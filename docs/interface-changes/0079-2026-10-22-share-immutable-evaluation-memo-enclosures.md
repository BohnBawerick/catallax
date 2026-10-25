## 2026-10-22: share immutable evaluation-memo enclosures

Piece: 1
Change: one real-expression evaluation shares immutable ball values between memo entries and their callers. The memo is dropped before returning the owned root ball. Square-root clipping replaces only its local value.
Why: memo insertions and hits copied arbitrary-precision midpoint and radius mantissas.
Migration: memo scope, expression traversal, domain checks, rounding and work charges are unchanged. No public ownership or thread-safety contract changes. The complete numbers suite passes, including cold/warm budget boundaries and known-defined square-root clipping.

