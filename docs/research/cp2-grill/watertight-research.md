# Watertight-by-design B-rep: precedent, cost, verdict

Question from the Catallax owner: why can't watertight-by-construction methods (Sederberg 2008, Urick 2019/2020, the hybrid Bezier-triangle paper) just be adopted, what would a kernel that is watertight all the way take (including STEP in and out), and has anyone shipped it?

Source labels:
- [read] means I fetched and read the primary text. Downloaded copies are in `scratchpad/dl/`.
- [snippet] means search-result text or a summary of a fetched web page.
- [unread/memory] means my own background knowledge, not checked this session.

## 0. The core reason it is not easy (the math)

- An exact surface-surface intersection cannot in general be a shared low-degree spline edge.
  - The image of a cubic trim curve on a bicubic patch has degree at most 18 and genus zero.
  - The intersection of two bicubic patches has degree 324 and is generally not genus zero.
  - So conventional trimmed NURBS must leave gaps. [read: Sederberg et al. 2008, sec. 2]
- It follows that every watertight method perturbs at least one surface near every intersection edge.
  - Sederberg 2008: "involves a perturbation of the original surfaces." [read]
  - Marussig & Hughes 2018: "approximations are introduced at least in the vicinity of trimming curves. The number of control variables increases..." [read: ICES 17-03, sec. 3.4 summary]
- Analytic-on-analytic cases are not exempt.
  - Two cylinders with skew or general axes meet in a quartic of genus 1, which has no rational parametrisation. [unread/memory]
  - So even a cylinder-cylinder edge cannot be shared exactly unless at least one face is rebuilt as a spline.
  - Inference: "watertight" and "keep my analytic faces exact" conflict on most Boolean edges.
  - Exception: some low-degree cases share boundaries exactly today. Plane with cylinder gives a conic. A planar profile extruded gives side faces that share the profile curve exactly.
- Mainstream kernels therefore chose tolerance.
  - Li et al. 2024 (Cocasoft) say T-splines are "not currently adopted by mainstream CAD systems", are not an ISO standard, and raise development and maintenance cost. [read: arXiv 2402.10216]
  - The same paper says Urick's reparameterisation cannot accurately control its effect on the original surface shape. [read: arXiv 2402.10216]

## 1. Precedent: products and kernels

| System | Watertight by design? | STEP in / out | What it loses |
|---|---|---|---|
| Autodesk T-Splines, then Fusion Form | Yes, inside a Form (T-spline) body only. Solids are ordinary tolerant B-rep. | Export: T-spline bodies must be converted to B-rep, otherwise they are dropped from STEP/IGES [snippet: Autodesk support article]. Import: no conversion from a B-rep body to T-spline; only "BRep Face to T-Spline", one face at a time, often messy [snippet: Autodesk forum accepted answer]. | Editability after conversion; the B-rep is tolerant again. |
| Rhino SubD | Yes for the SubD object (limit surface, no trims). | Export through ToNURBS as packed or unpacked NURBS patches. | Extraordinary vertices can only be approximated: G0, G1, G1x or G2 options, each trading deviation for smoothness [snippet: Rhino docs, Grasshopper docs]. One forum benchmark: packed conversion of a huge SubD took about 66 min, unpacked about 9 s [snippet]. |
| NX Realize Shape / Solid Edge subdivision | Subdivision cage converted to B-rep [unread/memory]. | Same pattern as Rhino. | Same pattern. |
| Coreform | U-splines are watertight. They are built from a hex mesh in Cubit, not converted from trimmed CAD [read: Coreform IGA page]. | Exports BEXT (Bezier extraction) for LS-DYNA [snippet]. Flex immerses trimmed CAD in a background spline mesh and does not convert it [read: Coreform IGA page]. | The analysis company co-founded by the IGA and T-spline people chose to tolerate trimmed CAD for import. |
| nVariate WatertightCAD (Urick) | Implements WBO. Planned as a plugin inside existing CAD, not a kernel [snippet: NIST SBIR award text]. | QIF bracket STEP converted to an untrimmed multi-patch NURBS STEP [read: OSTI 1635598]. | See section 3. No shipping product found. Patents granted through 2026 [snippet: Justia]. |
| nTop (implicit) | Yes. No trims; field-defined. | STEP/Parasolid out by meshing then Quadrangulate then "CAD Body from Quad Mesh", or by "CAD Body from Implicit Body" (beta) [snippet: nTop support]. | Complex lattices fail or explode in size: a 600 MB STEP from a model under 30 kB in .ntop [snippet]. Analytic faces only survive through a "Preserved CAD" sewing input. |
| Hyperganic / voxel tools | Watertight as fields or voxels [unread/memory]. | Mesh export [unread/memory]. | No exact B-rep. |
| Parasolid convergent modelling | No. Mixes facet B-rep with classic B-rep under tolerant modelling [snippet: Siemens v32/v33 posts]. | Classic B-rep. | Not applicable. |
| SpaceClaim / Discovery | No. Direct modelling on a tolerant ACIS-family B-rep [unread/memory]. | Not applicable. | Not applicable. |
| Mesh kernels (e.g. Manifold) | Yes, as guaranteed-manifold meshes [unread/memory]. | No exact surface export. | Exactness. |

Conclusion for 1 (inference from absence):
- I found no shipping general-purpose solid kernel that is watertight by design and round-trips STEP.
- Every shipping case keeps the watertight representation inside one body type (T-spline, SubD, implicit, mesh).
- Each converts to tolerant trimmed NURBS B-rep at the boundary with solids and with STEP.

## 2. What the research was for, and where adoption stopped

- Sederberg 2008 was written for design.
  - Goals: watertight fillets, feature curves that need not follow isoparameter lines, and gap-free NURBS export. [read]
  - Limits it reports:
    - Bicubic only.
    - Extraordinary points are patched to G1 (Peters) or G2 (Loop).
    - The least-squares fit converges linearly, can take tens of seconds, and can hit local minima.
    - Merging two T-splines with mismatched parametrisations kinks the surface, which "most pairs" of converted surfaces have.
    - Example perturbation error: 2.5e-4 of patch width. [read]
- Urick, Marussig, Hiemstra and Hughes worked for isogeometric analysis (IGA) and analysis-suitable CAD. [read: IGA2019 abstract; Hiemstra et al. 2020]
  - nVariate later added additive manufacturing (patent US20230030783) and digital-thread use. [snippet]
- Song & Cohen 2019 (IDETC, "Making Trimmed B-Spline B-Reps Watertight With a Hybrid Representation"):
  - Keeps the original surfaces except in narrow bands around trims.
  - Seals the bands with Bezier-triangle regions that match the model-space trim curve exactly. [snippet]
  - Xia & Qian (rTBS) is a similar triangular-Bezier conversion. [read: summary in Marussig & Hughes]
  - I found no adoption of either.
- Autodesk T-Splines timeline:
  - Autodesk bought T-Splines in Dec 2011.
  - The SolidWorks plugin (tsElements) was dropped about three years later.
  - The Rhino plugin was dropped on 7 Jan 2017.
  - The technology now lives in Fusion and Alias Form tools. [snippet: SolidSmack, Autodesk press release]
  - The T-spline patent US 7,274,364 expired in 2024. [snippet: Wikipedia only, not verified at USPTO]
- Marussig & Hughes 2018 conclude that "there is no canonical way to deal with trimmed models, neither in analysis nor in design." [read]
  - Most IGA publications use local, immersed methods that leave the geometry untouched. [read/snippet]
  - Coreform Flex follows the same pattern commercially.

## 3. What a kernel that is watertight all the way would take

### Surface types

- Needed: an unstructured spline type with T-junctions, mixed degree and extraordinary points (T-splines, U-splines or LR-splines). Plain multi-patch NURBS is also possible, but then continuity is hard to maintain:
  - Hiemstra: "cumbersome to maintain geometric continuity between patches unless continuity is built into the basis." [read]
- Analytic faces are converted.
  - In the QIF bracket test, 152 cylinders, 98 planes, 52 tori, 8 conics and 2 spheres all became splines. [read: OSTI 1635598]
  - Hiemstra's fitting criteria list "conic shapes should be represented exactly" and "flat faces should remain flat" as goals that can conflict with the rebuild. [read]

### Booleans and fillets (Urick WBO)

- Three stages [read]:
  1. Analyse the parameter space and split it into u-type and v-type regions.
  2. Reparameterise.
  3. Update the model-space surface so it interpolates the approximated intersection curve.
- WBO runs after a conventional Boolean, reusing its SSI output. So it is a post-pass on a tolerant kernel, not a replacement for one. [read: CAD 2019 abstract]
- Degree grows with the reparameterisation function.
  - Marussig & Hughes: "the degree of the resulting patches may become large." [read]
  - Li et al.: an m x n patch becomes m x (m*p + n). [read]
- Control-point growth, from the only industrial number found (QIF bracket: 318 faces, 834 curves, tolerance 0.000515 in) [read: OSTI 1635598]:
  - Original, cast to pure NURBS: 6,468 surface control points + 5,016 model-space trim control points + about 5,016 parameter-space trim control points, about 16,500 in total.
  - Watertight result: 20,429 surface control points, a 23.8% increase.
  - About half of those are duplicated shared-boundary points that a proper data structure would store once.
  - Raising everything to one degree would add "significant" further growth.
  - The run was done by "non-automated application of scripts".
  - Extraordinary points of valence 3 and 5 were added to avoid degenerate slivers.
- Fillet continuity was lost.
  - WBO uses only tangent data from the intersection curve, so faces from different source surfaces end up only G0.
  - Example: a fillet between two planes and a cylinder loses its tangency. Normal and binormal handling was listed as future work. [read]
- Runtime: no published wall-clock numbers found for WBO. Sederberg 2008 reports tens of seconds per fit (2008 hardware). [read]

### STEP import

- Every trimmed face has to be rebuilt.
- Hiemstra's global method needs a valid triangulation of the B-rep first, and "current CAD systems do not have all the necessary tools" for complex or defective models. [read]
- Fusion refuses whole-body B-rep to T-spline conversion. [snippet]
- Inference: real-world STEP needs healing before conversion, and the conversion success rate on dirty supplier files is unknown.

### STEP export

- Option A: untrimmed multi-patch NURBS (B_SPLINE_SURFACE faces with coincident boundaries). nVariate did this for the QIF bracket. [read]
  - Geometry is effectively gap-free.
  - T-junction and extraordinary-point structure, editability and analytic types are all lost.
  - The face count rises.
- Option B: AP242 edition 2 (2020) has `locally_refined_spline_surface`, covering T-splines, LR B-splines and hierarchical splines. [snippet: steptools AP242e2 notes, ISO 10303-42:2021 schema]
  - This corrects the earlier note: STEP has had such an entity since the 2018 Part 42 update.
  - I found no translator that reads or writes it, and SINTEF says support is "limited". [snippet: arXiv 2010.01980]
  - nVariate says no public format handles their unstructured splines; .tsm is bicubic and uses T-NURCC extraordinary points. [read: OSTI 1635598]
- Option C: AP242 tessellated. Watertight but faceted.
- Option D: convert back to trimmed NURBS. This reintroduces the receiver's tolerances.

### Round trip

- Lossy in practice (inference from the evidence above).
  - Analytic types do not come back.
  - Reimporting your own multi-patch export means rebuilding again or accepting a tolerant B-rep.
  - Geometry stays within tolerance, but semantics (cylinder axes for holes, planes for mates and CAM) are gone unless re-detected.

## 4. Failure modes and open problems

- Extraordinary points:
  - Sederberg's extraordinary points are G1 or G2 patched. [read]
  - Rhino approximates them on conversion to NURBS. [snippet]
  - Converted T-splines distort near extraordinary points (GA boundary-optimisation follow-up paper). [snippet]
- Continuity: G0 only between rebuilt faces from different sources (WBO). [read]
- T-spline merges kink when the parametrisations do not match. [read: Sederberg]
- Linear dependence of T-spline bases; mixed degrees across a model. [read: Marussig & Hughes]
- Editing after conversion:
  - Sederberg and WBO both claim editability. [read]
  - In practice the model is a control net. A user can no longer "change the hole diameter" without feature history. [inference]
- Analytic shapes are lost. [read: OSTI]
- File size:
  - About +24% control points with a good data structure. [read]
  - Implicit to STEP can blow up by orders of magnitude. [snippet: nTop]
- Robustness on dirty input is unproven at scale. [read: Hiemstra; Marussig & Hughes: "the step from academic examples to practical multipatch models is perhaps the most challenging task"]
- No mainstream data structure or format exists for the result. [read: OSTI]

## 5. Verdict for Catallax

Evidence:
- No shipping precedent for option (a).
- The best industrial data point is one bracket, done by hand-run scripts.
- It cost about 24% more control points, all analytic faces, and fillet tangency.
- No STEP reader consumes the native entity.
- The IGA founders' own company immerses trimmed CAD rather than converting it.

Inference:
- **(a) Watertight all the way: no.**
  - It requires inventing the data structure, the Boolean, the fillet with tangency and the import healer, all of which are open research.
  - Every STEP exchange still degrades to either multi-patch NURBS or tolerant trimmed NURBS.
  - It discards the analytic faces that CAM, PMI and mating need.
- **(b) Tolerant B-rep with watertight construction where the kernel owns both sides: yes, as a targeted feature.**
  - Cheap exact cases:
    - Extrude and revolve of planar profiles.
    - Lofts and sweeps whose rails or profiles are the shared edges.
    - Fillets on analytic supports. Composing a pcurve with a plane or cylinder gives a low-degree exact curve; a bicubic support gives degree up to 18. [read: Sederberg bound]
  - Record an "exact-shared" flag per edge so downstream code can skip tolerance logic.
  - Optionally offer WBO-style export of a multi-patch NURBS solid for IGA or additive-manufacturing users, as an export mode, not the core representation.
- **(c) Tolerant only:** the industry baseline. Combined with Catallax's SDF evaluator, robust inside/outside queries need not depend on watertightness. Generalised winding numbers over trimmed NURBS are an existing method. [snippet: arXiv 2504.11435 title only]

Recommendation: (b), built on (c).

## Could not verify

- Full text of Urick et al. CAD 2019 and JCISE 2020 (paywalled). Read only the abstract and the OSTI report.
- Song & Cohen 2019 full text.
- Exact USPTO expiry date of the T-spline patent.
- Whether any CAD translator implements `locally_refined_spline_surface`.
- Current state of Fusion Form (still shipping as of the 2017-2024 sources).
- WBO runtime figures.
