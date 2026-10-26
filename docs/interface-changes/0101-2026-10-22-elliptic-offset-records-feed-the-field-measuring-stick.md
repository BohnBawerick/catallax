## 2026-10-22: elliptic-offset records feed the field measuring stick

Pieces: shapes, field and offset construction consumers
Change: `catallax_field::measurement::FaceDefinition` implements `From<&EllipticOffsetDefinition>`, retaining the exact original supporting surface and nonnegative offset magnitude required by the landed checker.
Why: construction records and the independent measuring stick must use the same definitions without repeating unsigned-distance conversion at each caller.
Migration: construct the checker's input from a borrowed validated record and keep the original record for its inward/outward side. The field metric is unsigned; it does not certify that branch, a finite construction domain or a fitted B-rep carrier. No existing measurement, exchange or geometry path changes.

