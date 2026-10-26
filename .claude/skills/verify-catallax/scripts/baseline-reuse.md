# Explicit baseline evidence reuse

`control-catallax.py qualify INPUT --baseline BASELINE_BINARY --candidate
CANDIDATE_BINARY --out NEW_EVIDENCE --json` runs fresh baseline and candidate
jobs at 1, 2 and 4 workers. It never resumes a qualification. Use the ordinary
`run` and `compare` commands for separate diagnostic experiments. `pair`
(including the operational SINGLEPAIR workflow) is unchanged.

To permit reuse, supply `--baseline-identity BUILD_RECORD.json` and explicitly
select `--reuse-baseline EARLIER_EVIDENCE`. Only an earlier fresh six-run
qualification produced by this command can supply baseline evidence. No archive
search, global cache or automatic migration of old proof bundles exists.
Candidate 1/2/4 always executes fresh with no cache or resume, even when both
executables have identical bytes. Reused qualifications cannot seed another
archive. Already running commands are unaffected.

Qualification requires a nonempty interpreted corpus, including when a trusted
adapter is selected. An empty adapter result is rejected before any output
directory, worker job or archive is created; selecting reuse cannot bypass this.

`--fresh-audit` forces six. A missing archive, unknown build provenance,
changed complete key, missing evidence or failed integrity check also selects
six. An **intact** archive that fails current semantic validation stops before
starting jobs: diagnose the baseline/validator defect and explicitly
requalify. Fresh-job defects stop qualification without retries. Preserve the
failed evidence; do not rerun until green or relabel it as a geometric refusal.

## Build provenance

The build record is a trusted operator attestation captured at build time,
bound to the executable hash. This command does not infer flags, native
dependencies or compiler provenance from a binary. Do not reconstruct missing
information with placeholders. Missing or unknown information means six runs
and no reusable archive. The source commit must also match the executable's
`about` record. Exact required fields are:

```json
{
  "format": "catallax-build-evidence/1",
  "source": {
    "commit": "<40 lowercase hex digits>",
    "tree": "<40 lowercase hex digits>",
    "dirty": false,
    "submodules": []
  },
  "executable_sha256": "<64 lowercase hex digits>",
  "build": {
    "toolchain": {
      "rustc": {"version": "<full version>", "sha256": "<64 lowercase hex digits>"},
      "cargo": {"version": "<full version>", "sha256": "<64 lowercase hex digits>"},
      "linker": {"version": "<full version>", "sha256": "<64 lowercase hex digits>"}
    },
    "dependencies": {
      "cargo_lock_sha256": "<64 lowercase hex digits>",
      "resolved_sha256": "<64 lowercase hex digits>",
      "native_sha256": "<64 lowercase hex digits>"
    },
    "options": {
      "command": ["<complete build command and arguments>"],
      "env_sha256": "<64 lowercase hex digits>",
      "config_sha256": "<64 lowercase hex digits>",
      "target": "<target triple>",
      "profile": "<profile>",
      "features": []
    }
  }
}
```

Toolchain hashes identify the actual compiler, Cargo and linker executable bytes;
version strings alone are insufficient. Submodules, if any, contain `path` and exact `commit`. Dependency hashes cover
the complete resolved dependency sources and native build inputs, not merely
package names. Build environment and config hashes cover all effective
compiler/linker flags, Cargo configuration, generated build inputs and options.
Empty arrays explicitly mean none; use a hash of an explicitly empty native
dependency inventory when appropriate. These are integrity/provenance records,
not signatures or protection against a malicious attesting operator.

The reuse key contains this complete record plus the executable identity,
ordered effective requests including duplicate occurrences, interpreted source,
adapter hash, default and per-request work allowances, embedded accuracy/options
and caps, driver/comparator hashes, worker counts/dispatch/core policy and
determinism environment. The latter includes a digest of the entire environment
(values are not published), Python identity, OS/machine identity, CPU identity,
affinity, resource limits and current dynamic-library bytes. Changing even an
irrelevant environment variable conservatively prevents reuse. Live CPU MHz
and bogomips are excluded as observations, not execution policy.

## Validation and limits

The directory archive retains both earlier triplets and a complete file
inventory. Every selection rechecks inventory, file hashes, absence of
symlinks, fresh provenance, worker exit codes, exact case order/completeness,
saved config, requests, original portable allowances, responses, manifests,
operation logs and STEP/IGES/OBJ/STL presence and bytes. JSON decimal numbers
are parsed exactly; interval strings are never shortened. Only established
time/build exclusions apply to geometric manifests. Same-build checks also
compare complete work and non-time logs. Current comparisons are recomputed
for both historical triplets; no stored `equal` verdict is authoritative.

Every `catallax-operation-log/1` entry must retain all protocol fields, including
nullable fields and observed `wall_ms`. Answered/refused entries require the
complete unsigned decision-counter record; interrupted entries require null
decisions. Refused entries require code/detail; answered/interrupted entries
require null refusal. Marks retain their typed complete schema and exact bound.
Successful `face_measure` local certificates must match the manifest's complete
face record and exact interval; other operations have null local tolerances.
Failed whole requests may retain null face certificates. These explicit nullable
forms preserve the current schema; omitted fields in older layouts are never
silently upgraded. Additive fields remain present in equality checks.

The existing face-check schema records internal B-rep display measurements
without an OBJ file. That specific request format is kept distinct from field
and print display records, which explicitly declare a delivered `display.obj`.

The gate deliberately requires strict manifest and artifact preservation
between builds. It cannot automatically approve a refusal becoming an answer,
new certificate or changed export. Those remain blocked for the parent's
existing preservation, independent geometry and changed-export audits. This
command is an opt-in strict operational subset, not a replacement for that
qualification policy. It adds no accepted-difference list or relaxed gate.

`qualification.json` labels every triplet member `fresh` or `reused` and retains
the selected archive descriptor hash. Historical baseline elapsed time is
explicitly historical. `elapsed_ms` measures the actual command span;
`archive_audit_ms`, `gate_audit_ms` and per-run observations describe their
own scopes. Do not sum overlapping jobs or historical timings into a makespan.
The existing preview `not-measured`, opt-in comparison/null and timing reports
remain unchanged. Reuse evidence does not establish fresh baseline latency.

Run the public synthetic mutation suite with the normal driver test discovery.
Its build attestations are deliberately synthetic fixtures. The implementation
alone establishes no saved time; a real future matching baseline is required.
Any driver/comparator change, including this one, still requires parent-owned
FULL SIX qualification after green CI. No private corpus is needed by these
tests, and none is included here.
