# Surface and wire bodies

Milestone 2.3 adds library-level surface edits, free spline wires and constructed-face sewing. These are geometry APIs, not new recipe operators or SolidDraft commands. Numerical bodies remain separate from exact models.

## Body roots and STEP

`BodyKind::Surface` contains face shells without asserting material. `BodyKind::Wire` contains nonempty ordered `Wire::coedges` paths and no lumps. Solid and surface bodies have empty `Body::wires`; empty bodies have neither wires nor lumps. Every existing body constructor supplies the new field.

Wire paths may be open or closed. Consecutive directed endpoints must share a vertex; a ring edge is a complete one-edge path. Coedges and edges cannot be repeated within the body's wire paths. Crossings do not imply a material or face-regularity claim. Body traversal, counts, journal replay and isometries preserve wires. Reflection preserves a wire's directed traversal rather than applying the outward-face loop reversal rule.

The `CXBODY` envelope now writes version 2 for certified and numerical handovers. Readers also accept the supported pre-wire version-1 layout and supply empty wire lists. Existing solid, surface and empty discriminants retain their values; wire is appended. Earlier numerical record layouts with obsolete positional native trims still refuse through their existing per-record tags.

STEP preserves the actual kind:

- Surface bodies use `SHELL_BASED_SURFACE_MODEL`; actual boundary uses decide whether each shell is open or closed.
- Wire bodies use `EDGE_BASED_WIREFRAME_MODEL` with one `CONNECTED_EDGE_SET` per path.
- Solids retain their closed material representation. Mixed exports use the generic shape representation.

Numerical export reconstructs and checks the emitted body with the same kind instead of forcing an open sheet through solid admission.

## Public operations

| Operation | API | Contract |
|---|---|---|
| Free spline wire | `catallax_topology::spline_wire` | Exact active-domain clamping and rational endpoint joins; no snapping or input mutation. Returns a fresh model and wire body. |
| Finite spline sheet | `numerical::surface::SplineSurfaceBody::new` | Retains the rational generating patch and builds a separate numerically admitted surface body. |
| Trim | `SplineSurfaceBody::trim` | Exact rational extraction of a nonempty parameter rectangle inside the current domain; rebuilds boundary topology and retains history. This is parameter trimming, not surface/surface intersection trimming. |
| Extend | `SplineSurfaceBody::extend` | Extends one chart side to an external parameter value. The numerical candidate retains the original parameter frame; an independent 9-by-9 sample grid checks the old patch within the caller's tolerance. |
| Thicken | `SplineSurfaceBody::thicken` | Retains the original rational cap unchanged. The signed unit-normal offset uses cubic interpolation with simple knots on every source-span subdivision and Greville nodes. Independent 7-by-7 checks per subdivision include both one-sided source-knot limits. Two caps and exact common-basis ruled rims share their vertices and edges explicitly. Both signs retain the source-face key and history. |
| Sew constructed faces | `catallax_import::numerical::sew_constructed_surfaces` | Takes native analytic face definitions, not STEP bytes. Reuses the import sewer, source checks and finite-parameter mapper; returns a separately checked numerical solid and the complete sewing/source-pair evidence. |
| Ruled patch | `RationalSplineSurface::ruled` | Exact Euclidean-linear ruling between continuous clamped rational profiles with the same active domain. A common homogeneous product basis retains continuity multiplicities rather than introducing redundant C0 joins. |

Extension and normal-offset checks are sampled. They establish neither a whole-domain approximation bound nor global injectivity or absence of self-intersection. Numerical STEP stays marked numerical. The exact generating definition is not an exact-body or exact-offset certificate.

Thickening does not relax the existing strong-monotone bounded-chart admission. A geometrically valid ruled rim can still refuse if that chart certificate is unavailable. The localized regression keeps a 66.7 nm bump within parameter support `[0.270, 0.274]` and checks both caps at 1 nm tolerance on a 1000 mm panel with 0.1 mm thickness. This aspect ratio admits the rims; it does not claim coverage of thinner rims outside that chart certificate.

Sewing uses the existing finite line-segment and identical finite circle/ellipse correspondence rules. Its displacement cap and the numerical reconstruction tolerance are independent explicit controls. Open, ambiguous or disconnected candidates never gain solid authority. The borrowed faces remain unchanged on failure. The report's constructed-definition fingerprint is a diagnostic identity for versioned in-memory text, not an imported STEP digest or a persistence codec.

## Algorithm choices

The donor pin remains SMLib `c9c5979f04cea8c38213004d4f04df0dcfb40d05`; existing port attribution stays in [THIRD-PARTY.md](THIRD-PARTY.md).

| Operation | Compared approaches | Choice |
|---|---|---|
| Wire topology | [SMLib wire shells and edge uses](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/.agents/skills/smlib-brep-model/references/brep-format.md); [OCCT MakeWire](https://occt3d.com/dev/doc/refman/html/class_b_rep_builder_a_p_i___make_wire.html) | Separate free directed uses, with exact ordered endpoint joins. No tolerance welding or partial-result publication. |
| Rectangular trim | SMLib homogeneous knot splitting; [OCCT rectangular trimmed surface](https://occt3d.com/dev/doc/refman/html/class_geom___rectangular_trimmed_surface.html) | Existing exact tensor extraction, not a fitted replacement or an intersection solver. |
| Extension | SMLib `N_SrfExtendByParamDist`, `NL_G1R`, in [NL_SrfAdv.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_SrfAdv.cpp); [OCCT ExtendSurfByLength](https://occt3d.com/dev/doc/refman/html/class_geom_lib.html) | Existing `catallax-trace::surface_extend` candidate, parameter remapping and independent retained-source samples. No new donor code or OCCT code is copied. |
| Normal offset | SMLib `SmOffsetSurface` unit-normal law and tensor interpolation; [OCCT OffsetSurface](https://occt3d.com/dev/doc/refman/html/class_geom___offset_surface.html) and [BSplSLib prescribed-knot interpolation](https://occt3d.com/dev/doc/refman/html/class_b_spl_s_lib.html) | Keep the original cap. Evaluate the signed normal law and use OCCT's prescribed-knot approach with the existing exact rational tensor solver. Source-span knots prevent a narrow source feature from falling between global nodes; simple cubic knots retain C2 joins. Independent per-span checks remain sampled. |
| Ruled rims | SMLib `N_CreateRuledSrf` blends in [NL_SrfAdv.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_SrfAdv.cpp); [OCCT two-curve BSpline filling](https://occt3d.com/dev/doc/refman/html/class_geom_fill___b_spline_curves.html) | Use a common spline basis, not a four-boundary Coons sum. Exact homogeneous products give Euclidean-linear rulings even when profile weights differ. Knot multiplicities retain input continuity so STEP decimal rounding does not depend on redundant C0-join constraints. No numerical admission rule changes. |
| Constructed-face sewing | SMLib stitch-then-solid validation; [OCCT Sewing](https://occt3d.com/dev/doc/refman/html/class_b_rep_builder_a_p_i___sewing.html) | The existing import lane's mutually unique boundary matcher and subsequent checks. Its private algorithms and byte-import entry paths are unchanged. |

## Reproduction and observed proof

Run on a build host, with scratch directories outside the checkout:

```sh
cargo test -p catallax-topology --test wires --test surface_thickening
cargo test -p catallax-app --test handover version_one_bodies_keep_geometry_history_and_empty_kind
cargo run -p catallax-step --example surface_bodies -- /tmp/cx-surface-bodies
cargo run -p catallax-step --example surface_bodies -- /tmp/cx-back-thicken -2/5
cargo run -p catallax-step --example wire_bodies -- /tmp/cx-wire-bodies
cargo run -p catallax-import --example constructed_sew -- /tmp/cx-constructed-sew
```

| Scenario | Construction and export work units | Observed result |
|---|---:|---|
| Curved sheet, exact trim, tangent extension, positive thickening | 15,269,091 | Three surface STEP files and a six-face solid; positive volume 12.808125482551887 mm³. |
| Same chain with negative thickening | 15,269,554 | Correctly oriented six-face solid; positive volume 12.792125695277566 mm³. |
| Open, multi-edge closed and single-edge closed spline wires | 9,455 | Three STEP wireframes, correct unique edge/vertex counts, no faces or solids. |
| Six independently constructed planar faces | 408,897 | Twelve recorded seam matches, a six-face solid of volume 24 mm³, and refusal of the open five-face sibling. |

The wire regression also verifies actual directed geometry after replay onto a grown model and reflection. Artifact imports used `cadquery-ocp==8.0.1.1.0`, direct `StepToTopoDS_Builder` surface/solid/wireframe entry points, `read.stdsameparameter.mode=0`, and `read.maxprecision.mode=0`. All written artifacts passed `BRepCheck_Analyzer` without `TransferRoots`, `ShapeFix` or an explicit healing pass. These checks do not establish global self-intersection absence.
