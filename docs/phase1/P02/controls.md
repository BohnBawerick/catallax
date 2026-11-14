# P02 controls: source inventory, not executed evidence

Frozen source: `f4db345f2e7c51835d1ebeba9f6446d117144117`.

These controls are present in the public repository. They were not executed in
this wave: Cargo/rustc and compiled test/example executables were unavailable.
No substitute fixtures, increased allowance, C++ oracle or private evidence
were used. Source hashes are retained in `source-provenance.json`.

| Control | Existing source / intended bounded command | Source inference only |
| --- | --- | --- |
| Tiny closed loops and tangent plane | `cargo test --release -p catallax-intersect --test patch_plane complete_search_finds_every_open_arc_loop_and_tangent_contact -- --nocapture` | Existing test uses 20M per case, checks tiny closed loop and named tangent refusal, and prints charged work. |
| Cylinder/sphere tiny loops and tangencies | `cargo test --release -p catallax-intersect --test patch_quadric complete_search_finds_open_arcs_missed_loops_and_tangent_contacts -- --nocapture` | Existing test covers tiny loops and named cylinder/sphere tangent refusals under 20M. |
| Border touch, tiny crossing and transformed controls | `cargo run --release -p catallax-intersect --example tangent_touch` | Existing example has 17 cases, five placements, workers 1/2/4 and 20M per request. It asserts answer/work equality across worker counts. |

The tangent example's `boundary` case reaches the solver's named
`prototype.tangent-boundary-contact` refusal path. A tiny crossing in that
example tests stationary-contact enumeration; it does not replace the tiny
closed-loop coverage tests. The prototype has no application/Boolean consumer.
The current research note says the prototype is parked: interior touches
classify, border and curve touches remain unresolved. Its historical 255-request
round and historical work table are source reports, not fresh P02 measurements.
Cylinder is unmeasured; cone, sphere and torus are outside the prototype's
implemented scope. A refusal/report string in that example is not a measured
surface request.

Firstmate needs an authorized Rust toolchain (or baseline-identified compiled
control executables) to complete these controls. No sanitizer or memory-fault
probe is necessary for this evidence gap, and none was attempted.
