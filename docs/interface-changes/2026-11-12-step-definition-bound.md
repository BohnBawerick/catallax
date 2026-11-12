## 2026-11-12: STEP states its bound against the exact definition

Pieces: 7
Change: `StepReport` gains `definition_deviation: Option<Rational>`. `write_step`, and the writers built on it, state the file's distance to the exact definition: `max_deviation` plus the worst measured fit bound of the fitted rolling-envelope, rolling-saddle and rolling-sphere faces written. For a file without fitted faces it equals `max_deviation`. Other writers leave it `None`.
Why: PLAN M12 asks for STEP of bodies with fitted faces with its own stated bound. `max_deviation` bounds the file only against the represented carriers, and only the certified recipe path composed it with the fit.
Migration: callers that build `StepReport` literals add the field or use `..Default::default()`.
