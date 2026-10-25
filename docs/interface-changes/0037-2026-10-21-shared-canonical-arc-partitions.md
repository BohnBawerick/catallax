## 2026-10-21: shared canonical arc partitions

Piece: 5
Change: imprint splits section arcs and operand boundaries with the same canonical cut pool before identifying shared edges. A graph section can therefore map to several common topology edges; their union is exactly the original section.
Why: splitting an operand edge without splitting a coincident section left overlapping arcs with different endpoints. Their duplicate outgoing germs had no strict radial order. Applying one partition removes the duplicate arcs instead of adding a tie-breaker. This resolves the two overlap-driven imprint-order refusals; the remaining fourth-order contact requires the derivative extension below.
Migration: no public Rust signature changes. Consumers of `Imprinted::common_edges` must use the complete common-edge set, not assume a one-to-one correspondence with graph sections. Exact tangent and curvature predicates, genuine ambiguity refusals and default budgets are unchanged.

