Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Optional browser timings

`OperationEntry.wall_ms` is now `Option<f64>`, and display `Cost.elapsed` is `Option<Duration>`. Both are `None` on `wasm32-unknown-unknown`, where the standard library has no clock. Native measurements and the native operation debug dump keep their existing values and format. Native Rust callers that require a measurement must unwrap the option.

Browser run-ahead uses one thread, including when a caller requests more. STEP fitting skips the timed helper-thread threshold on this target. Geometry, decisions and deterministic work counts do not depend on these observations.

Run `ci/wasm-clocks` with the wasm target and Node installed. It builds SolidDraft's browser dependencies, then executes the wasm logging, run-ahead and display checks without clock imports.
