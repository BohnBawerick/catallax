# Spline-face offset refusal reason

Nonzero ordinary B-rep offsets of a body containing a B-spline face now return
`distance.spline-offset-unsupported`. The detail names M10's missing certified
offset and trim arrangement. This replaces the generic envelope refusal for
those requests. The operation still returns no offset body.

Zero-offset and empty-body identity rules run first. Other analytic carrier
paths and cap-only spline shell operations retain their existing behavior.
