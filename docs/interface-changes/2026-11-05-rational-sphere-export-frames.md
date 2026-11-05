# Rational spherical STEP frames and chart budget refusals

For triangular sphere trims, STEP reframing now tries normalized 3/5 and
4/5 turns before the legacy diagonal frames. Other trims keep their original
candidate list. The stored-axis permutations remain first, and every old
diagonal remains available as a fallback. On rational unit bases, these
turns avoid the extra square root introduced by a 45-degree diagonal.

This changes only the exported parameter frame. Sphere centres, radii, senses,
trim curves, recipe origins and model geometry stay unchanged. The default
17 significant digits and 1/10000000 mm file tolerance also stay unchanged.
`StepReport.max_deviation` still bounds the written file; this is not a new
fitted fillet or a relaxed construction tolerance.

The radius-1 top loop followed by radius-1 upright rounds on the corpus's
20 by 20 by 10 block now exports within the 1M request allowance. OCCT's own
fillet operation refuses this case, so its regression checks raw import,
five material witnesses and the independent rounded-box section integral:

`V = 9 * (396 + pi) + 324 + 56*pi/3 = 3888 + 83*pi/3 mm^3`.

Failed local work, precision and decimal certificates during chart export now
remain refusals. They no longer fall back to a possibly unseamed legacy file.
The regression uses the sphere-dome recipe with separate 1M per-decision and
request allowances. Before the fix, `write_step` reported success but raw OCCT
rejected the file. The fixed call returns `step.budget-exhausted` without
changing its output buffer. This is unfinished proof work, not impossible
geometry. Unsupported encodings still retain their existing fallback rules.
The closed-form conic boolean test previously assumed that legacy STEP
fallbacks were valid. A vertical half-torus union disproves that assumption:
its old raw file fails OCCT validity and the independent volume control.
Construction and exact mass assertions remain; uncertified exports must
refuse rather than restore that file.

No public API signature, recipe format, radius law, overflow policy, display
mesh, SolidDraft plumbing or dependency pin changes. The wider construction,
whole-face proof and commercial-kernel comparison gaps remain open.
