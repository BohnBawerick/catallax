## 2026-10-26: direct finite conical face distance

Pieces: 2 (shapes), 6 (measure)

Change: `distance::Finite::recognize_oriented` additionally recognizes complete
positive-nappe circular cone bands between perpendicular coaxial rims.
`Finite::recognize_conical` accepts the same oriented rim inventory and an optional
single exact apex vertex representing the collapsed lower boundary.
`PreparedFace` adapts complete circle loops and a cone's apex loop to this route.

Why: the nearest point on the entire retained band is the aligned-azimuth
projection onto its closed radial/axial meridian segment, including both finite
rims and the apex. Repeated parameter subdivision is unnecessary for this domain.

Migration: existing callers and options are unchanged. New conical callers must
supply the complete boundary inventory, not a selected subset: two ring loops,
or one ring and one apex vertex loop. Each ring boolean must combine edge,
coedge and face senses as in the existing spherical contract. Natural travel
relative to the outward cone normal retains increasing height when its circle
normal agrees with the cone axis. An apex is exactly the carrier apex; it cannot
stand for an omitted finite rim. Unoriented recognition still declines cones.

Recognition proves positive height, circular incidence, coaxiality, perpendicular
planes, valid reference directions, and compatible lower/upper retained sides.
The opposite nappe and both carrier continuations are excluded. Arcs, angular
sectors, oblique/elliptic rims, extra loops, incompatible orientations and all
other unrecognized trims keep the existing complete search. This adds no new
unsupported-input refusal. Spherical metadata and effective sense semantics are
unchanged.

Direct answers have `cells == 0`, exact admitted witnesses and two-sided bounds
under the unchanged width, cell and work controls. If fixed-precision bounds miss
the requested width, the complete fallback still runs with the original options.
The existing `Branch` variants remain: interior projection includes endpoint
ties, an outside projection reports `Rim` even when that endpoint is the apex,
and an axial query reports `Axis`. No new branch variant is needed.

This library query does not enable shell/offset/thicken construction, native
field distance/display, or an application latency claim. The default 1M
scoreboard, 20M application and 100M verification allowances, output accuracy,
cell and vertex caps are unchanged.
