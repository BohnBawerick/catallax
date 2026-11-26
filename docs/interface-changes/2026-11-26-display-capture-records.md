## 2026-11-26: display capture records for display comparisons

Pieces: catallax-mesh, catallax-field, tools/display-diff
Change: `catallax-mesh` and `catallax-field` gain a `display-capture` feature, off by default. With it, `mesh_body`, `mesh_body_reusing`, `display_mesh_body`, `preview_mesh_body` and field display write one `displaycapture` line per call to standard error: record format 1, thread name, function, spline carriers, options, work counters and a SHA-256 digest of an explicit encoding of the result. The encoding is in `catallax_mesh::capture`; it does not use `Debug` text. `tools/display-diff/check.py capture` uses this interface by default; `--collector probe` keeps the old source rewrite for revisions without it, and `migrate` matches the two collectors on one revision.
Why: the probe rewrote four function bodies by regular expression and hashed `Debug` text, so a renamed local or a formatting change could break capture or move a digest without moving a triangle.
Migration: none for library callers; normal builds do not compile the feature. Compare captures from one collector only. An exception pinned against probe rows does not match interface rows, so re-pin it against an interface capture of both sides.

### Display answer check

The `display answers` workflow uses the interface collector when a pull request changes display code or its tests or capture tools. It compares the exact PR head at 1 and 4 workers and uses the PR's main base commit as its baseline. The job summary reports elapsed time, and the `display-captures` artifact contains the JSONL records and work counters.

When a mesh change is intentional, list the part and reason in the PR and put one JSON object per changed call in `tools/display-diff/exceptions.jsonl`. Copy `key`, `main.sha256` and `candidate.sha256` from the comparison into `key`, `main_sha256` and `candidate_sha256`; add a nonempty `reason`. An exception already on main cannot authorize another change. New or edited exceptions must match a change in this PR exactly. No exception can waive a lost call, a new refusal or a difference between workers.

Format 2 adds a final `ok` or `refused` field without changing the result digest. The gate accepts format 1 as a baseline, including main `6182a6141`. If a changed result refuses and that baseline does not expose its outcome, the gate fails: it cannot prove the prior result refused too. To run the CI comparison locally, capture main at 1 worker and the candidate at 1 and 4 with `capture --threads 1` and `capture --threads 1 --threads 4`, respectively, then run `compare <main-dir> <candidate-dir> --gate`.
