## 2026-10-22: seed polynomial content from a small coefficient

Piece: 1
Change: primitive integer normalization starts its coefficient GCD from a nonzero coefficient of minimum bit length, stops when the content is one, and avoids dividing an already primitive positive-leading polynomial by one.
Why: annihilator construction spent time computing content from very large constant terms even when a unit coefficient already certified primitive content.
Migration: no public signature, policy or normalized coefficient changes. Existing polynomial and root suites retain their rational-content, sign, multiplicity and independent-oracle checks. A 65,536-bit monic-polynomial smoke probe retains every coefficient.

