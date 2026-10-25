## 2026-10-22: skip zero and repeated Taylor products

Piece: 4
Change: interval Taylor multiplication and division skip coefficient ranges proved identically zero. Squared expressions and square-root recurrences compute symmetric products once, while retaining the original order of rounded additions and subtractions.
Why: high-order traces used by face membership and measurement repeatedly multiplied the same coefficient pairs or scanned zero tails.
Migration: no public signatures, trace definitions, orders, precision or rounding rules change. A before/after public-API probe retained identical midpoint, radius and precision fingerprints for 57,456 coefficients across circles and quartic branches, Taylor orders zero through twenty and four working precisions.

