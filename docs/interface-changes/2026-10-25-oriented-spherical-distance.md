## 2026-10-25: oriented finite spherical distance

Pieces: 2, 6

Change: `catallax_shapes::distance::Finite::recognize_oriented(surface, rims, dec)`
accepts complete boundary circles paired with their effective forward traversal.
`true` means natural circle travel with the carrier's outward normal; combine
edge, coedge and face senses before calling. It recognizes a spherical cap or a
complete coaxial band only after proving exact circle incidence and compatible
retained sides. `PreparedFace` supplies this metadata for single-ring loops.
`Branch::Axis` also identifies spherical axis/centre cases; stored rim references
resolve circle ties. Recognized point queries have `cells == 0`.

Why: one unoriented rim on a sphere bounds two complementary regions. Carrier
recognition alone cannot determine the finite retained face. For a recognized
height interval, maximizing the point's dot product over the sphere admits the
radial sphere point only inside that interval; otherwise the maximum lies on a
finite rim. The sphere centre ties everywhere and axial queries include admitted
poles or rim-circle ties. Exact witnesses and both distance bounds remain required.

Migration: the measurement adapter combines edge/coedge travel with face sense
and uses the new entry point. Existing `Finite::recognize` callers remain valid;
that unoriented method still declines spheres. Non-spherical recognized families,
widths, cell/work limits, refusal codes and all unrecognized-trim searches are
unchanged. Degenerate or incompatible rings are not promoted to spherical faces.
The per-piece interface review and final acceptance gates remain PR review work.
