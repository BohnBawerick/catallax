# Adaptive print failure progress

`printing::adaptive::Failure` adds a `progress: Progress` field. It retains
actual counts when an opt-in adaptive request stops, including inside candidate
construction or either interval certificate. The original typed refusal,
request, stage, arithmetic work and charged request work remain intact.

`candidate_round` is zero-based and absent before construction. `cells` counts
visited candidate and certificate cells across attempts, using the existing
counter without clamping it to the cap. `vertices` and `triangles` describe
the current rejected candidate. `unresolved_cells` counts homogeneous cells
observed so far; it is absent before traversal. A partial count of zero is not
a proof that no unresolved cells remain.

These are observations, not predicted total demand, peak memory, calibrated
time or permission to use a partial mesh. A preparation refusal has zero
visited cells and no candidate round. A failure during a later candidate keeps
cumulative cells and that candidate's partial geometry counts. The report
captures certificate work before propagating its original failure.

No certificate, material-topology guard, tolerance, work allowance or refinement
limit changes. Legacy print requests do not acquire these fields or a new
meshing route. Successful output and its cell certificates keep their existing
meaning. The [M19 verdict](../parked/m19-printing.md) records the measured
controls and remaining proof and calibration work. Raw evidence retained
privately.
