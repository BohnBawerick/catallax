## 2026-10-18: ND2 numbers are decimal or p/q strings

Pieces: judge interface, cli
Change: question numbers use piece 1's grammar, a decimal with an optional exponent or a fraction `p/q`, and nothing else. Every number in the manifest is a decimal in the notation of `to_decimal_bounds`, never a fraction.
Why: both forms parse exactly; writing only decimals keeps the manifest readable by any decimal parser.
Migration: `render_manifest` writes bounds, mark bounds and `max_deviation` with `to_decimal_bounds`.

