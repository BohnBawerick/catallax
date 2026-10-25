## 2026-10-24: correction to sampling depth and aggregate application work

Pieces: 1, 5, application

This corrects the preceding "stack-safe region sampling and bounded exact display
samples" entry. A fixed depth of 32 rejected previously answered near-contact
and scaled geometry. Region sampling now uses an explicit heap worklist in the
same left-before-right order as the old recursion. The inherited numeric depth
guard remains, but thread-stack use no longer grows with subdivision depth.
The unchanged 4096-cell allowance is shared across every piece, projection
choice and precision retry of one signature query, rather than renewed for each
chart build. Its exhaustion still returns `boolean.sampling-budget`.

An independent phase probe also found a valid displaced-hole build using
88 million arithmetic units before meshing. A sampling-only bound did not cap
that complete application edit. `catallax_app::DEFAULT_BUILD_WORK` is now
20,000,000 units for one construction/display request. The normal bracket uses
9.8 million units. Both exact and certified application entries share this
allowance and discard any result after exhaustion with the new
`Error::ResourceLimit`, code `application.work-limit`. Numeric per-decision
policies remain in force.

`catallax_numbers::operation_budget_active` reports whether the caller already
owns an aggregate work scope. Application builds preserve that explicit
allowance rather than nesting, resetting or enlarging it. In particular, the
verification exchange retains its existing per-request allowance and outcomes.
The default application limit is not an allocator cap or a wall-clock deadline.
The 4096 exact display sample bound and unchanged isolated process caps from the
preceding entry remain in force. No coordinate-size cutoff, geometric snapping
or tolerance increase is introduced.

Migration: consumers of the non-exhaustive application error type must retain
`application.work-limit` as a resource refusal. Callers with an explicit
operation budget keep control of their request allowance. There is no recipe,
manifest or STEP format migration.

