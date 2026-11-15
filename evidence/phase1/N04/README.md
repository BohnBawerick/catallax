# N04 numerical edit evidence

Source tested: `822847575290fd42676f7eb95a5f079efb6d28b2`.

The baseline regression at `7a90e83e` failed after a no-op `replace_surface`:
`a later edit discarded the source sewing checks`. The final regression checks the
retained source evidence through the typed history API and confirms that old
sewing checks do not appear as checks on the edited body.

`recipe.json` is the generated save of the N02 square tube with a 0.005 mm split
endpoint, explicit sewing under 0.01 mm, height changed from 1 to 2 mm under
0.002 mm, then a cut removing material above z=1.5 under 0.002 mm. The result
retains all ten source face origins and adds `edit:trim/cap` to the cut cap. One
source repair remains historical. Two earlier checked states keep their own
controls and actual checks.

`edited.step` is the actual approximate STEP from that replay. Fresh example
runs at `CATALLAX_THREADS=1`, `2` and `4` produced identical STEP bytes and identical
`replay-summary.json` records, including 698,253 shared work units for parse,
source admission, edits, export and written-byte readback. The independent
feature test measured 506,206 work units for construction and edits alone.

`raw-reference.json` records OCCT 7.8.1.1 reading these exact bytes with the raw
builder and healing disabled. It checks one valid solid, ten faces, 24 edges,
16 vertices, two faces with holes and volume 18 mm3. Material probes check the
kept wall, removed cap region, hole, outside and new cap. A subsequent OCCT cut
returns a valid 15 mm3 solid.

Reproduction:

```sh
cargo run --release --locked -p catallax-features --example numerical_replay -- evidence/phase1/N04/recipe.json /tmp/n04.step
python3 -B tools/step-regressions/check_numerical_edit.py /tmp/n04.step
cargo test --release --locked -p catallax-topology -p catallax-features -p catallax-step -p catallax-app --test numerical_trim --test numerical_edits --test numerical_step --test numerical -- --nocapture
```

The scoped release tests passed, including the N03 exact legacy byte/work control.
Adverse controls cover an over-tolerance current residual, invalid dimensions and
cuts, duplicate edit IDs, modified source gaps, a low precision cap and a shared
work limit of one. A refused CLI replay removed a pre-existing stale STEP file.
No tolerance, precision cap or deterministic allowance was increased.

Coverage remains rational planar coordinate-frame edits and cuts across a
straight two-layer prism. Curved or general numerical Booleans and app checkpoint
integration remain outside this API. Numerical samples do not establish a
complete body bound.
