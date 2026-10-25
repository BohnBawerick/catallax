## 2026-10-22: omit Taylor coefficients unused by panel certificates

Piece: 6
Change: integration programs lazily compile center and remainder dependency plans into fixed-size boxed slices, keeping both plans independently lazy without spare-capacity metadata. Symmetric center integration needs even output coefficients; its whole-range certificate needs only the order-20 output coefficient. Products, variable quotients and roots retain every lower coefficient needed by their recurrences. Constant divisors retain their zeroth coefficient and requested numerator coefficients.
Compiled node, output-index and exact-constant definitions are immutable boxed slices too. The builder counts retained nodes while marking dependencies and allocates their final buffer once; compaction adds no allocation or coefficient arithmetic. Prepared planar area state retains only its numerator and squared normal. Optional open-torus program bundles are indirect, so their payload is not embedded in every curved face.
Why: evaluating every output coefficient forced arithmetic that never entered the integral or its remainder bound.
Migration: requested coefficient arithmetic keeps its order and outward rounding; definedness checks and nonsymmetric angular integrals remain intact. No public signature, requested accuracy, default allowance or deadline changes. Regressions use exact antiderivatives and keep undefined zeroth-order operations visible.

