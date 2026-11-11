# M20 aligned field cantilever: run policy

UNPROVED numerical FEA; not suitable for safety sign-off. This bounded first
result does not complete M20's plate, sliver or imported-part exits.

Frozen baseline: cdc0a3884c59b2bea538c230e37639bd74cdb9a2.
The public recipe cyclically rotates a 5 x 5 x 100 mm native block so the
beam occupies x=[0,100], y=z=[-2.5,2.5] mm. Start/End remain the support/load
names. E=200000 N/mm², nu=0.3, traction=(0,0,-0.04) N/mm², force=-1 N.

Policy fixed before runs: Q1 displacement, eight-point Gauss volume integration,
consistent face loads, direct fixed-DOF elimination, diagonal-preconditioned CG,
recomputed residual <=1e-9, at most 10000 iterations. Retain the original solver
caps of 80 cells per axis and 5120 total. Refinement grids in internal
(y,z,x) order: [2,2,20], [4,4,40], [8,8,80]. Each run has the unchanged
100M operation allowance shared by geometry and solver. Charge one unit per
assembled cell and one per sparse matrix row evaluation; also report actual
matrix-entry products. Postprocessing reserves one unit per cell centre and
eight units per path/station (the maximum incident traces). No clock selects a numerical answer. Caller cancellation
through numeric work checks is distinct from convergence failure. Parallel row
arithmetic keeps each row's order and all reductions canonical; all available
cores by default. Geometry preparation and authoritative queries remain serial.

F7 exact cell integration must return fraction=[1,1], a connected native box
and one certain homogeneous material piece matching each cell. The conservative
closed-box state is retained, including Unknown on boundaries. F7's retained
axis-box certificate supplies coverage, never samples. Named box patches must
survive, have area 25 mm² and the prescribed centroids. Point membership and
material queries are additional checks, not coverage proofs.

Tip is the consistent area mean. Centreline and both outer-fibre traces retain
every incident cell's tensor and von Mises value; station summaries give ranges.
Cell-centre tensors are retained separately. Pictures use unsmoothed computed
cell traces, never beam-formula stresses. Stress errors are comparisons only;
there is no separate stress percentage pass. The displacement acceptance is
abs((-tip)/0.032-1)<=0.02. Clamp/load end regions and corner concentrations are
not beam-theory stress gates. SVG colours and finite trace samples are
numerical display evidence, not global stress bounds.
