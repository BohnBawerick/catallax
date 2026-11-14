# Publication sanitation

The `raw` directory contains sanitized copies of the original P03 execution records. Original bytes are privately preserved in existing authorized scratch, not republished here.

Only machine-specific locations were replaced: the executable uses `target/release/catallax-verify`, scratch paths use `<PROOF>`, and the report uses a portable checkout description. `<PROOF>` is a placeholder in JSON metadata, not a literal execution destination. In shell examples `$PROOF` must be set to a new authorized scratch directory. Restricted toolchain-location prose was generalized.

The six refusal outcomes, five successful control artifacts, source identities, exact inputs, 20M allowance, work counts, timings and original binary identity are unchanged. No new geometry execution, toolchain probe or test occurred.

Each `artifact.json` now hashes its published case files. `SHA256SUMS.json` hashes the published evidence files excluding itself. Original case keys are retained as execution identifiers; they were not recomputed as sanitized-content identities. Changed metadata has new digests; original digests do not authenticate sanitized bytes. Unchanged STEP, response, manifest and portable-question bytes retain their existing digests.

This follow-up removes paths from the current published branch tree. It does not erase the earlier commit from Git history; no history rewrite or force push was performed.
