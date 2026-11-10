# Three-kernel comparison data

Run from the checkout to collect real data without touching `site/`:

```sh
python3 -B tools/kernel-compare/collect.py \
  --out /tmp/catallax-comparison-release \
  --scratch /tmp
```

The output directory must be new and outside the checkout. Requirements are Linux x86-64, Rust, `uv`, and enough disk space for a release build and the pinned Python wheels. The collector installs the existing five-tool benchmark's hash-locked Python 3.11 dependencies, builds with `-j 6`, pins execution to the first allowed CPU, and runs kernels sequentially at nice level 10. Rust, OpenMP and BLAS worker counts are limited to one during measurement. It removes its own temporary build, Python installation and download cache on exit. Evidence stays in `--out`.

For a tagged release, add `--source-tag v0.1.5`. This requires a clean, committed branch. The collector temporarily checks out that exact tag in this worktree, builds and preserves the producer binary, then restores the branch before building the observer and running cases. It creates no additional checkout. The driver checks the binary's embedded commit against the tag, and records the source tag, kernel revision and separate harness revision. A rebased tooling commit is not relabelled as the release binary.

The collector writes `results.json` with full observations, `data.json` with a smaller reviewable projection, `summary.md`, and per-case commands, stdout, stderr, producer STEP files and observer reports. To retain data on a comparison branch without editing site pages:

```sh
python3 -B tools/kernel-compare/report.py /tmp/catallax-comparison-release/results.json \
  --data data/kernel-compare-runs/release.json \
  --summary data/kernel-compare-runs/release.md
```

The projection accepts only a completed inventory with all three engines present for every case. Raw per-sample intervals remain in the evidence directory. Their digest and the raw result digest bind the projection to that run. Nothing deploys, regenerates site pages, or imports into `site/comparison-results.json`.

After an adapter correction, pass `--supplement /path/to/narrow-rerun/results.json` to the projection. It accepts unchanged case definitions, kernel binaries and completion policy only. Latest observations enter the matrix; superseded observations and each run's source, script hashes and raw digest remain in the data. The 9 November dataset retains the original 375 producer attempts plus a five-case, fifteen-attempt adapter rerun. That rerun fixes split circular contours and trimmed-face witness selection, and distinguishes released-package refusals from adapter errors. It does not change any requested geometry or tolerance.

## Inputs

The current inventory has 125 case variants:

- All 45 recipes in `parts/fillets/cases.jsonl`, including station setbacks and linear radius laws.
- All 13 five-tool targets, with their authored construction and scalar-edit recipes, for 26 variants. Their original recipe hashes must match the existing benchmark pins. This runs the requested three kernels, not SolidCraft or Pali Ijo.
- The 43 distinct recipe files referenced by `include_str!` in `crates/catallax-compare`.
- All 11 representative recipes from `catallax-compare`'s `measure` example, exported by `measure --export` rather than copied into another case list.

The fillet adapter reuses the independently authored OCCT geometry and edge witnesses. Its existing ruled-section linear-radius reference and explicit station-window construction remain labelled in the raw native report. SMLib executes its public native radius-law overload. No constant-radius substitution, smaller-radius retry, healed STEP, or successful prefix stands in for the requested part.

## One completion policy

`policy.py` owns every kernel's geometric thresholds and solve deadline. The sweep is 0.1, 0.01, 0.001, 0.0001 and 0.00001 mm. The solve deadline is 120 seconds for the whole producer process, including startup, complete construction and STEP serialization. Observer costs are separate. A timeout or an answered process over the deadline is `unsolved-time-limit`, never a geometric refusal.

Each unchanged producer result is checked at every threshold. This is a threshold slider over recorded geometric error, not five hidden construction retries or changes to a kernel's modelling tolerance.

`done-sampled` requires all of the following:

1. The producer answered within the solve deadline and exported the complete result.
2. Direct STEP import, without repair, has valid topology, at least one solid and positive finite volume.
3. All authored material witnesses pass.
4. Both boundary-sample checks finish, with an upper distance within the chosen threshold.

The reference is the definition-native Catallax field compiler, run independently of B-rep construction. Four divisions produce rational reference candidates; trim queries retain only surviving boundary points. Each producer's imported STEP boundary is also sampled on a five-by-five trimmed face grid, five points per edge and its vertices. The reverse observer computes point-to-definition distance intervals under a fixed work allowance. The forward observer measures reference points against the complete imported face boundary, with a recorded 0.0000001 mm numerical allowance.

This is finite sampled evidence, not a global Hausdorff, tangency, self-intersection or correctness certificate. The field compiler shares recipe semantics and authorship with Catallax B-rep. OCCT both competes and imports every STEP file. These dependencies remain in the data; agreement cannot establish their independence.

An unavailable reference definition, missing upper bound, unresolved trim candidate or observer timeout is `unresolved`. It cannot become a pass. A proven sample outside the threshold, invalid raw topology or failed witness is `not-done`. Maximum topological tolerance and analytic volume differences are separate diagnostics, never substituted for geometric distance in millimetres.

Native refusals, crashes, protocol errors and adapter gaps remain separate outcomes. An adapter gap is not attributed to the kernel. SMLib entries come only from real invocations of the installed released package. This collector does not label another tool's output as SMLib.

## Provenance and interpretation

Results retain source revision, the Catallax binary's embedded identity, binary and observer hashes, script and input hashes, package versions, installed native-library hashes, release pins, host CPU and operating system, CPU affinity, load averages and UTC dates. Shared-host solve times are observations, not quiet-host timing claims or cross-kernel speed rankings.

For a narrow smoke run, repeat `--case` with exact IDs, for example `fillets/box-single`, `five-tool/P2-collar/construction` and `comparison-measure/sphere`. `run.py` also accepts prepared binaries and Python paths to avoid rebuilding during development. Use a fresh output directory for every run.

Check the shared acceptance boundaries with:

```sh
python3 -B -m unittest discover -s tools/kernel-compare -p test_policy.py -v
```

The 9 November smoke run exercised all three producers on those three cases. All nine results passed at 0.001 mm. At 0.00001 mm, SMLib's P2 collar failed with an observed distance interval of approximately 0.0000250544 to 0.0000252544 mm; Catallax and OCCT passed. This smoke is not the full corpus result.
