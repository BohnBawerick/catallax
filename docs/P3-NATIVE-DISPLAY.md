# P3 native field display

The unchanged `parts/P3-filleted-boss.json` completes native display with its existing plane, cylinder and torus supports at the fixed 1,000,000-work scoreboard allowance and default 0.01 mm display tolerance. The public display controls, vertex/cell limits, typed resource refusals and elementary fillet construction remain unchanged.

The evaluator represents certified quarter-arc and ring-fillet contact regions with exact polynomial predicates. Sharing those predicates across Boolean operands preserves the dependence between the boss and adjoining fillet at buried seams. Exact chart restrictions, regional signs, material germs and the existing deferred-cell cover establish the display's whole-boundary distance certificate. General algebraic charts and new fillet carriers are outside this slice.

The public P3 display stage uses 885,114 work, with 75,811 vertices, 108,962 triangles and a reported 0.01 mm bound. Independent tests check all thirteen surviving origins against the B-rep, analytic surviving domains, outward triangle winding, planar area controls and valid B-rep topology. Strict P3 STEP bytes are unchanged. Zero-work and exhausted cell limits return refusals without a mesh.

Development observations remain raw timing misses against the 1,000 ms complete-click target. A recorded complete CLI wait is about 17.49 seconds; earlier 35- and 38-second observations remain misses. These are unqualified development measurements, with no speed claim or replacement baseline.

The paired comparison retains its original accuracy (0.01 mm), B-rep chord tolerance (0.005 mm), divisions (4) and normal cosine (999/1000). After integrating the generated-band row reuse, the unchanged 1,000,000-work scoreboard request completes in 993,304 work, with all 9,842 forward vertices, 332 reverse samples, no findings and no refusals. Its exact sampled upper bound is unchanged. The separate 100,000,000-work verification control checks the same complete result. This is sampled two-way agreement, not a whole-result certificate.

Historical evidence remains distinct: the frozen pre-integration runtime needed 1,008,614 work and refused at 1M; the earlier projector candidate used 993,013 work with a different 10,098-vertex mesh. Neither is a replacement baseline for the current result. Generated-band row reuse preserves mesh positions, triangles, rails, face bounds and all geometric checks while avoiding the separate intermediate-strip traversal.

The larger verification allowance is not substituted into the scoreboard. Exhausted work still returns a typed refusal with no incomplete agreement. Rational-input plane, cylinder and ring-torus proposals use the same parametrisations and binary rounding as the previous projector, sharing the radial-root enclosure across torus phases. Existing unsupported-frame, algebraic-input and low-precision fallbacks remain intact; the original trim and actual-distance checks still certify each candidate. The 2,000 ms background-comparison target remains unchanged; a recorded complete CLI wait of 21.87 seconds remains a raw miss. Complete target-part acceptance and local timing qualification are still pending.

Reproduce field display through the verification skill:

```sh
cargo build --release -p catallax-cli --bins -j 2
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run \
  parts/P3-filleted-boss.json --source field --max-work 1000000 \
  --binary target/release/catallax-verify --out "$PROOF/p3-field" --json
```

Run the same public recipe separately with `--source compare` and inspect its forward/reverse coverage, complete findings and resource status. The 1M scoreboard, 20M application default and 100M verification default remain separate contracts.
