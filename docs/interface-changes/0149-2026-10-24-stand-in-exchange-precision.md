## 2026-10-24: stand-in exchange precision

Pieces: CLI stand-in exchange
Change: The explicit `catallax-stand-in/1` adapter writes 32 significant digits instead of inheriting ordinary STEP's 17-digit default.
Why: Once field inertia reaches its unchanged 1e-9 width, 17-digit STEP rounding can exceed the retained field enclosures. Tightening the exchange removes that error without widening a field or surrogate bound.
Migration: No input or format change. Ordinary recipe STEP options and bytes are unchanged. The adapter still refuses unless every exported property enclosure fits the original field bounds.


