# Retry a certified interior sphere pole on a perpendicular STEP frame

The STEP writer retries its perpendicular spherical frame when a trimmed face
has a certified interior pole but none of the attempted latitude caps clears its
trims. Previously that preparation result ended the chart attempt and kept the
legacy encoding. Earlier answer 431 followed that path and was raw-invalid.

The retry is export-only. The built body, exactness, marks, measured properties
and native topology do not change. An uncertified pole still keeps the legacy
encoding. Unsupported charts, non-budget preparation errors and the existing
short-edge guard retain their fallback behaviour. Request budget and deadline
failures still propagate. A charted result must stay within the requested STEP
deviation before the writer uses it.

The public regression intersects a radius-five ball with a block whose lower x
face is at -0.001. Its volume is `261.87792761444194`. Both direct no-healing
translation and the standard STEP reader must yield one valid solid with that
volume within the existing area-times-tolerance comparison.

Merge correction metadata now admits the approved `invalid-to-valid-step` class
only for an index on the base revision's active known-raw-invalid list, at its
trusted approved STEP hash. The declaration pins before and after digests.
For new corrections, CI builds the accepted base revision and replays only the
changed indices. The base manifest must match its accepted full comparable hash
before the checker derives the protected-field digest from it. The candidate's
commitment and protected fields must match that authenticated base value; a
self-supplied commitment cannot authorize a new exactness, mark or property.
The direct STEP audit also checks the corrected solid's volume and maximum
deviation. The merge-check contract is documented in
[RELEASING](../RELEASING.md#merge-checks).

This merge corrects 431 only. Firstmate authorized splitting the original V8
merge A under its E-4 rule. Answer 262 keeps its approved legacy bytes and remains
on the known-raw-invalid list. Its touching-trim chart correction follows as
merge A2. No application pin or release-correction file changes here.
