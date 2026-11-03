## 2026-11-03: STEP and a certificate for a body with a loop fillet

Pieces: step, features
Change: `write_step` and `write_rolling_junction_step` now answer for a body whose fillet rounds the loop of two crossing cylinders, where they refused with `step.unsupported-geometry`. `certify_junction` answers too, where it refused with `numbers.geometric-gap-unproven`, so `catallax-verify` answers the loop boss. Each wall's chart cut goes where its rail never does. The chart code accepts a rail span that lies wholly behind the cylinder's axis as clear of its cut.
Why: a part with a loop fillet needs a file, and `catallax-verify` needs that file and a certificate to answer.
Migration: a SolidDraft flow that showed the STEP or certificate refusal for these bodies now gets an answer. `mass_properties` and face measures still refuse with `measure.unsupported-surface`.
