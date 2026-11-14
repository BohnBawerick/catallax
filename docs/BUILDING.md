# Building and verification

Run these commands from the repository root with a Rust toolchain supporting edition 2024. Save the [complete legacy question example](JUDGE-INTERFACE.md#a-complete-example) as `question.json` before running `check` or `answer`. Recipe questions use the separate verification CLI.

```sh
cargo build --release -p catallax-cli --bins
target/release/catallax about
target/release/catallax check question.json
target/release/catallax answer question.json --out answer
```

`check` validates without evaluating geometry. `answer` writes `answer/manifest.json` and, for a nonempty solid when requested, `answer/answer.step`. Refused and invalid questions also write a manifest and exit 0. [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md) defines the exact command and file contract, including its unchanged legacy names.

Local checks:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
python3 -B -m unittest discover -s tools/edge-hunt -v
```

CI selects checks by changed paths. Full-suite runs include formatting, clippy, enabled release tests and doctests. Crate `TESTS.md` files describe extra opt-in workloads. Timing benchmarks report measurements, never correctness gates.

## Verification

`catallax-verify` adds a persistent, deterministic exchange for questions, recipes and sampled two-kernel comparisons. Whole-request work budgets choose outcomes; wall time is reported separately. The legacy `catallax` and `exact-brep` commands retain their existing deadlines and file formats.

```sh
cargo build --release -p catallax-cli --bins
.claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary target/release/catallax-verify --json
.claude/skills/verify-catallax/scripts/control-catallax.py --help
```

The [verification skill](../.claude/skills/verify-catallax/SKILL.md) covers corpus replay, caching, resume, strict STEP comparison and matched-pair timing. The [versioned contract](VERIFY-INTERFACE.md) defines work accounting and the expanded geometry kind tags. Pass private inputs as runtime paths and keep proof artifacts outside the repository.


## Contributing

Read the design rules first. A guessed topological decision, an unproven construction fact, an `unsafe` block or a silent fallback is a defect. Keep API and format changes explicit, update callers and document the migration in a new [interface-change fragment](interface-changes/README.md). Bring mathematical evidence for geometry changes and report refusals honestly.

Catallax uses AI coding agents, directed and checked by a human. Most of the inherited implementation was written that way. Whether that approach scales to a useful CAD kernel remains an open question.

Agent-written code is welcome. Say how it was produced and stay accountable for what it does. This repository's independent tester was built without seeing the kernel's code.
