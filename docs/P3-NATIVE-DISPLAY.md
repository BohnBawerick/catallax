# P3 native field display

The unchanged `parts/P3-filleted-boss.json` completes native display with its existing plane, cylinder and torus supports at the fixed 1,000,000-work scoreboard allowance and default 0.01 mm display tolerance. The public display controls, vertex/cell limits, typed resource refusals and elementary fillet construction remain unchanged.

The evaluator represents certified quarter-arc and ring-fillet contact regions with exact polynomial predicates. Sharing those predicates across Boolean operands preserves the dependence between the boss and adjoining fillet at buried seams. Exact chart restrictions, regional signs, material germs and the existing deferred-cell cover establish the display's whole-boundary distance certificate. General algebraic charts and new fillet carriers are outside this slice.

The public P3 display stage uses 885,114 work, with 75,811 vertices, 108,962 triangles and a reported 0.01 mm bound. Independent tests check all thirteen surviving origins against the B-rep, analytic surviving domains, outward triangle winding, planar area controls and valid B-rep topology. Strict P3 STEP bytes are unchanged. Zero-work and exhausted cell limits return refusals without a mesh.

Development observations remain raw timing misses against the 1,000 ms complete-click target. The latest recorded complete CLI wait is about 17.49 seconds; earlier 35- and 38-second observations remain misses. These are unqualified development measurements, with no speed claim or replacement baseline.

The paired comparison retains its original accuracy (0.01 mm), B-rep chord tolerance (0.005 mm), divisions (4) and normal cosine (999/1000). At the separate verification default of 100,000,000 work, it completes in 4,503,230 work with 10,098 forward vertices, 332 reverse samples, no findings and no refusals. This is sampled two-way agreement, not a whole-result certificate.

At the unchanged 1,000,000-work scoreboard contract, comparison still exhausts its allowance with unresolved observations. It remains a typed refusal and does not establish scoreboard agreement. The larger verification allowance is not substituted into the scoreboard. The 2,000 ms background-comparison target remains unchanged; development comparison timings also exceed it. Complete target-part acceptance is still pending.

Reproduce field display through the verification skill:

```sh
cargo build --release -p catallax-cli --bins -j 2
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run \
  parts/P3-filleted-boss.json --source field --max-work 1000000 \
  --binary target/release/catallax-verify --out "$PROOF/p3-field" --json
```

Run the same public recipe separately with `--source compare` and inspect its refusal, forward/reverse coverage and all findings. The 1M scoreboard, 20M application default and 100M verification default remain separate contracts.
