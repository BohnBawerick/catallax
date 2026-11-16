# 2026-11-16: cheaper body adoption

`read_body` now caches validated rational encodings for the duration of one
transfer. The first occurrence still checks the positive denominator and runs
`Rational::new`, then compares the normalized result with the encoded integers.
Only an identical pair of length-prefixed integer encodings can reuse that
result. New encodings, including unreduced fractions and invalid denominators,
pass through the existing checks. One decode owns the cache and lends it
exclusively to its reader views. The cache does not survive the request.

Integer decoding also borrows its byte slice instead of allocating and metering
a temporary vector one byte at a time. It charges the entire slice before
materializing the integer. Cache hashing, copying and growth charge the active
request allowance before doing that work. Reading never obtains a new allowance
inside an existing request.

Shared records now use a dense indexed table instead of hashing their consecutive
IDs. The decoder checks the minimum record width and charges table allocation
before reserving it. Duplicate, out-of-range, unresolved and mistyped references
still refuse. Fixed arrays decode into stack storage instead of temporary vectors.
Small byte and scalar readers permit cross-crate inlining; their work and deadline
checks remain in place.

The wire layout and version 1 bytes are unchanged. Checksums, graph references,
handle bounds and typed refusals keep their existing behavior. This adds no
geometric algorithm. It retains the existing numbered shared-record design
described in [the original interface note](2026-11-15-body-handover.md).

## Measurements

The fixture is `crates/catallax-cli/tests/fixtures/soliddraft-gui-bracket.json`.
Its finished body occupies 3,779,427 bytes before and after this change.
Construction, output destruction and display meshing stay outside read timing.

The temporary exclusive-time profiler found the following baseline costs.
Native profiling ran on the work PC; wasm profiling used Node 22 on the Blade.
Instrumentation adds overhead, so these rows identify the hotspot rather than
replace the uninstrumented end-to-end measurements.

| Stage | Native ms | wasm32 ms in Node |
| --- | ---: | ---: |
| Envelope and container decoding | 15-19 | 27.9 |
| SHA-256 | 1.8 | 16.9 |
| Shared graph decoding | 7.7-8.3 | 16.9 |
| Retained JSON decoding | 1.0 | 2.8 |
| Rational decoding outside validation | 3.5-3.7 | 8.5 |
| Rational canonicalization | 86-92 | 254.9 |
| Integer binary decoding | 26-28 | 58.4 |
| Allocation and freeing | 23-25 | 44.1 |

One read visits 23,694 shared records and 20,922 rationals. The profiler counted
483,200 allocation/free calls natively and about 525,175 on wasm32. Canonicalizing
repeated rationals was the largest stage; retained JSON was not the bottleneck.

Uninstrumented release measurements on the work PC:

| Runtime | Read before ms | Read after ms | Write before ms | Write after ms |
| --- | ---: | ---: | ---: | ---: |
| Native, median of seven reads | 138.0 | 39.1 | 46.8 | 44.1 |
| wasm32, Node 24.11.0 median of seven warm reads | 146.7 | 49.9 | 34.9 | 38.1 |
| wasm32, Chrome 143 median of seven warm reads | 431.4 | 169.1 | 137.8 | 152.4 |

The encoder is unchanged; its native timing variation is not an encoder speedup.
The final standalone Node probe reports a 123.0 ms cold read. Chrome cold reads
fall from 432.8 ms to 187.6 ms. Warm measurements follow three unrecorded reads.
Import work changes from 12,705,134 to 13,830,917 units because
cache operations also consume work. Export work remains 23,872,737 units.

The standalone Chrome kernel result does not meet the requested margin below
150 ms. That target remains open. SolidDraft was not rebuilt or retested here;
its performance comparison will happen after the app can pin this kernel revision.

The Blade was overloaded during its later Chrome measurements, with load
averages above 37 on a 12-thread host. Those runs are not a controlled timing
comparison.

## Reproduce

Run the builds on the designated build host, from the repository root. The native
example accepts an optional second argument to save the finished body:

```sh
cargo run --release -p catallax-app --example handover_profile -- \
  crates/catallax-cli/tests/fixtures/soliddraft-gui-bracket.json /tmp/bracket.handover
cargo build --release --target wasm32-unknown-unknown -p catallax-app --lib
```

Compile the standalone probe against those wasm libraries:

```sh
target_dir=${CARGO_TARGET_DIR:-target}
deps=$target_dir/wasm32-unknown-unknown/release/deps
externs=()
for package in catallax_app catallax_features catallax_numbers; do
  lib=$(find "$deps" -name "lib${package}-*.rlib" -print -quit)
  externs+=(--extern "$package=$lib")
done
rustc --edition 2024 -C opt-level=3 --target wasm32-unknown-unknown --crate-type cdylib \
  crates/catallax-app/examples/handover_browser/probe.rs "${externs[@]}" \
  -L "dependency=$deps" -L "dependency=$target_dir/release/deps" -o /tmp/handover-probe.wasm
node crates/catallax-app/examples/handover_browser/measure.mjs \
  /tmp/handover-probe.wasm /tmp/bracket.handover
```

For Chrome, serve the wasm file, body and `measure.mjs` from a temporary localhost
directory. Import `measure` in the page and pass the two fetched ArrayBuffers.
Use an unthrottled browser. The function reports the first read separately from
seven warm read/write samples. It never builds geometry in the page.

All 13 number-codec tests, 15 application handover tests and 9 numerical-trim
tests pass. The existing
codec and early-abort tests are unchanged. New cases cover repeated rational
round trips, allowance exhaustion on cached reads, and noncanonical or invalid
fractions after a cache hit, invalid fixed-array elements, and duplicate or
out-of-range shared IDs. The CI-selected all-target clippy check and
`ci/run fmt` also pass.

The run-ahead source-state inventory guard passes under nextest 0.9.148's CI
profile with the exclusively borrowed cache. A native-generated bracket payload
also decodes in wasm32 and re-encodes to the same 3,779,427 bytes, with read work
13,830,917 and write work 23,872,737. The pre-fix archived shard reproduced the
inventory failure; an additional graded-field test reached the local 600-second
limit despite passing hosted CI in 152.8 seconds.
