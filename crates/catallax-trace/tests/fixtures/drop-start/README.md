# Singular drop-start gate

`cases.txt` records `SmSurfaceDropCurve::ComputePointValues` from SMLib
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`. The test probe initializes the
working interval and domain explicitly before calling the method. It records
the surface point, derivatives, normal, projected curve tangent, UV point,
UV tangent and two drop distances. Cargo never builds or links this C++.

The cases follow parameter curves from collapsed sides in both parameter
orientations, at corners and interior side points, close to the singularity
and at regular points. A collapsed point curve has no usable tangent and
returns a typed failure. At a singular start, a point 0.0001 along the
projected tangent is minimized onto the surface at the source's 0.001
search tolerance. That drop chooses the parameter of the inward first step.
These are candidate search constants, with no authority to prove geometry.

The Rust unit gate can access the private tracer without adding a public
API. It compares all fields within the existing 16-ULP candidate allowance
and checks identical candidate bits with 1, 2 and 4 workers.

Build `oracle.cpp` with the flags in `../foundation/README.md`, then run
`python3 generate.py "$ORACLE" "$RAW_DIR"`. Raw evidence retained privately.
