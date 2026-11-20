# Tolerance numbers research (Catallax, units mm)

Labels: [read] = I read the primary text (source file, PDF, vendor doc page) in this session.
[read-excerpt] = seen only as a search-engine summary of the page, not the page itself.
[unread/memory] = from memory or inference, not checked in this session.

Local copies of downloaded sources: `scratchpad/dl/` (vda41.txt, vda3.txt, jackson.txt, ps_prec.html.txt, acis_tol.htm.txt, arx.txt), `scratchpad/usdbrep/repo` (SMLib), `scratchpad/occt/occt` (OCCT sparse checkout).

## 1. Published PDQ guidance

| Source | Finding | Label |
|---|---|---|
| VDA 4955 v4.1 (Dec 2006, 93 pp, German) | Defines G-ED-TI tiny edge, G-LO-LG large edge gap, G-FA-EG edge-face gap, G-FA-VG vertex gap, G-FA-NA narrow face, G-FA-TI tiny face, G-SH-FR free edge. **No numeric thresholds anywhere in the text.** 3.2.1.3 O-CM-AP says the accuracy parameters (gap, tangency angle, tiny element) must match the "Unternehmensstandard" (company standard). The tiny-edge rationale: edges under a minimum size "can lead to invalid elements on exchange into a system of lower accuracy" and loss of trims. Only number in the document: a knot-tolerance example of 0.0001. | [read] |
| VDA 4955 v3 (2002) | Same criteria, same lack of numbers. | [read] |
| SASIG PDQ = ISO/PAS 26183:2006 (withdrawn) | Same criterion codes. Excerpt describes "tiny" as a user-set accuracy above a "zero length" accuracy. No numbers found. | [read-excerpt] |
| JAMA/JAPIA PDQ | Nothing found online. | not found |
| ISO 10303-59 (PDQ-S), ISO/TS 8000-311 | Schema for criteria, measurement and inspection results; thresholds are parameters carried in the data, not fixed by the standard. | [read-excerpt] |
| CAx-IF, LOTAR, NIST | No numeric model-tolerance guidance found. CAx-IF rec-pracs cover PMI and validation properties, not gap thresholds. | not found |
| CADInterop (vendor) | SASIG/VDA/JAMA share core criteria with "numerical thresholds proper to each standard"; micro-edge/micro-face "business threshold (often 0.1 mm)" because they break meshers. | [read-excerpt] |
| Siemens community | Tolerant-edge tubes "typically 0.001 mm to 0.01 mm". | [read-excerpt] |
| Memory | Commonly cited OEM/Q-Checker defaults around 0.01 to 0.02 mm for gaps and tiny edges. Not verified. | [unread/memory] |

Conclusion: the standards name the checks and leave the numbers to each OEM. The de facto numbers come from receiving systems (section 3), which all model at 0.01 mm.

## 2. Kernel defaults and caps

### OCCT (source read, master)
- Precision::Confusion = 1e-7 [read]. Header comment: tuned for mm; "1/10000 micron, rather difficult to measure" [read].
- Precision::Approximation = Confusion x 10 = 1e-6. "You may use a smaller tolerance in an approximation algorithm, but this option might be costly." [read]
- Precision::Intersection = Confusion x 0.01 = 1e-9; PConfusion = 1e-9; Angular = 1e-12 [read].
- BOPAlgo_Options fuzzy value default = Confusion (1e-7); SetFuzzyValue floors at Confusion [read].
- ShapeFix_Root: myPrecision = myMinTol = myMaxTol = Confusion by default; if MaxTol < Precision it is raised to Precision [read].
- ShapeFix_ShapeTolerance::LimitTolerance(shape, tmin, tmax, type): any vertex/edge/face tol above tmax is set to tmax, below tmin set to tmin; no re-measurement, so it can make a shape invalid [read].
- STEP/IGES resource files: FixShape.MinTolerance3d = 1e-7; Tolerance3d = Runtime.Tolerance; MaxTolerance3d = Runtime.MaxTolerance [read].
- XSTEP statics: read.precision.mode = 0 (File: use the file's uncertainty), read.precision.val = 1e-3, read.maxprecision.mode = 0 (Preferred), read.maxprecision.val = 1.0 (mm), write.precision.mode = 0 (Average), write.precision.val = 1e-3 [read for values; meaning of mode names from memory].
- "Preferred" means the max may be exceeded by healing; "Forced" applies LimitTolerance [unread/memory].
- BRepLib_ValidateEdge: tol = distance x 1.00001 [read]. BRepCheck_Edge returns sqrt(tolCal) x 1.05; BRepCheck_Vertex sqrt(Tol x 1.05) [read]. Approx_SameParameter 1.05 x deviation [given in brief, not re-read].

### Parasolid (V35 docs, chapter "Session and local precision", read)
- Session linear precision 1e-8 model units, angular 1e-11 rad; units strongly recommended as metres, size box 1000 units (1 km) [read]. In mm: 1e-5 mm, box +-500 000 mm.
- Exact (non-tolerant) edges and vertices are treated as having half session precision: 5e-9 m = 5e-6 mm [read-excerpt of XT format reference; consistent with Jackson 1995 "default 5.0e-9 metres", read].
- Separation guideline: vertices not meant to coincide more than 100 x precision apart (1e-6 m = 1e-3 mm); edges not meant parallel more than 100 x angular precision [read].
- Performance: "Precision in the range 1.0e-6 to 1.0e-8 modeller units is likely to result in the generation of SP-curves with large numbers of segments" [read]. In mm that band is 1e-3 to 1e-5 mm.
- PK_EDGE_optimise lowers an existing tolerance towards the max deviation between attached curves; short edges (all points inside a vertex bubble) are skipped by default [read].
- Setting a smaller precision without replacing SP-curves can make the model fail checking; "we recommend that applications do not decrease precision in this manner" [read]. So lowering must re-fit, not just rewrite the number.
- PK_VERTEX_set_precision: vertex precision not less than the max of its edges [read].
- Loop gap closing default 1e-6 (m) = 1e-3 mm [read].
- Scaling imports by up to factor 2 is safe [read].
- Solid Edge (Parasolid) import: stitch at 1e-5 m = 0.01 mm; models from looser systems arrive with many tolerant edges, which a simplify step rebuilds to session precision [read-excerpt].

### Jackson 1995, "B-rep modelling with local tolerances" (read)
- Old global Parasolid tolerance 5e-9 m, "chosen to be below the smallest realistic tolerance for most mechanical parts".
- Single small tolerance forces surfaces with "a large amount of defining data" (blends); single large tolerance blocks small features.
- "It is not possible to split a short edge, of length less than 4 times the tolerance, without creating edges which are too short."
- Vertex tol >= edge tol >= face tol; entities intersect if within the sum of their tolerances.
- Section 7.3 "Preventing tolerance growth": compression (merging) grows local tolerances; "models defined to a certain tolerance should unite to create models of the same or similar tolerance ... tolerance growth should be contained where possible." No numbers.

### ACIS (R17 docs, read)
- SPAresabs = 1e-6 (coincidence), chosen with "at least an order of magnitude guard band".
- SPAresnor = 1e-10 (ratio smallest/largest; largest representable = 1e4 units). Do not change.
- SPAresfit = 1e-3; should be a "reasonable multiple (e.g. 1000)" of resabs.
- SPAresmch = 1e-11.
- SAT files store resabs/resnor; mismatch on restore gives a warning [read-excerpt, RobotStudio forum].

### SMLib (usd-brep source, read)
- `SmBrep::ValidateAndUpdateTolerances` (source/SMLib/src/SmBrep.cpp ~33405-33736): `dSize = boundingBox.GetSize().Length()` (diagonal). Edge tol bumped to `2.0 * dMaxEdgeFaceGap` only `if (dMaxEdgeFaceGap < dSize/800.0)`. Header comment: "Edge gaps larger than BrepBoundingBoxSize/800 are not fixed. It's assumed such big gaps are symptoms of more serious Brep problems." No derivation of 800 given. Rebuilds the UV trim when the Edge/Surface gap is good but Edge/UVTrimCurve gap is bad.
- SmTol.h ZoneTol3d by model size estimate: 0.00005-0.005 -> 1e-7; 0.005-0.5 -> 1e-6; 0.5-500 -> 1e-5 (default); 500-50000 -> 1e-4; 50000-5e6 -> 1e-3; outside: 10^floor(log10(size)/2) x 1e-5. I.e. tol grows as sqrt(size).
- XSectTol3d = ZoneTol1 + ZoneTol2 (= 2 x Zone); ApproxTol3d = Zone/2 (fitted geometry sits inside half the zone).
- Planned gap model (comment): EdgeZoneTol = max(Def, connected faces, GapGain x MaxGap_VertexEdge/2, GapGain x MaxGap_VertexFace - Def). GapGain value not found.
- "Pinhole in a battleship": user flags small topology to get smaller tolerances.

### CGM
- Nothing public found.

### Other
- arXiv 2310.10351 (Boolean result repair): per-intersection-edge tolerance t = k x d(f1) x d(f2) x d(e0) x t0, where d = polynomial degree and k = max(fit error / t0) of the inputs; fit errors often exceed t0 by over 1000 x [read]. Empirical, but it ties tolerance to the input's measured error, not to model size.
- Autodesk patent US 11016470: "keeping the gap size < 10 eps renders tolerant modeling operations as reliable as precise model operations" [read-excerpt].
- Fang, Bruderlin, Zhu 1993: tolerances updated dynamically during Booleans to keep relations consistent (planar/quadric) [read-excerpt].

## 3. Receiving systems

| System | Default | Notes | Label |
|---|---|---|---|
| Rhino | 0.01 mm absolute | "works best ... absolute tolerance around 0.01 to 0.001", small features >= 10 x tolerance, model size <= 100000; "smaller than 0.0001 will noticeably slow some intersection and fitting processes"; import does not change doc tolerance | [read] |
| Rhino STEP import | uses tolerance written in the STEP file, with a floor (formerly 0.0001); CATIA mm STEP files usually declare 0.005; McNeel: exporters often declare a tighter tolerance than the model was built to, so join fails; fix is override to e.g. 0.001 | [read-excerpt] |
| Siemens NX | Modeling distance tolerance 0.01 mm since NX 9 (was 0.0254 mm), chosen as automotive de facto standard; angle 0.5 deg; B-spline approx tol default 0.0254, tighter makes export slower and larger; tolerant vertex can reach sqrt(2) x 0.01 = 0.01414 | [read-excerpt] |
| Creo | Absolute accuracy default 0.01 mm (0.00039 in) since Creo 7.0 | [read] |
| Creo relative | 0.0012 default, range 0.01-0.0001 (lower bound configurable 1e-6..1e-4); STEP/IGES import at relative accuracy, Parasolid/ACIS tolerant and CATIA V5 at absolute | [read-excerpt] |
| CATIA V5 | Join/healing merging distance 0.001 mm default, max 0.1 mm; heal tangency 0.5 deg | [read-excerpt] |
| SolidWorks | Parasolid inside (session 1e-5 mm); no official knit-tolerance default found; third-party says 0.001 mm tight / 0.01 mm general | [read-excerpt] |
| Onshape | Parasolid inside; no import tolerance setting documented | [read-excerpt] |
| Cross-kernel | CADInterop: Parasolid data at 1e-6 mm can show micro-gaps in ACIS at 1e-4 mm | [read-excerpt] |

Fail behaviour: no source shows a receiver rejecting a file because it is tighter than its own tolerance. The documented failures are the other way: real gaps larger than the declared or receiver tolerance give open shells or "faulty topology". Tighter source data mostly costs size and speed (Parasolid SP-curves, NX B-spline export, Rhino fitting).

## 4. Choosing tolerances: what the evidence says

- Absolute vs relative. Kernels use an absolute resolution with a fixed size box (Parasolid 1e11 ratio, ACIS 1e10, OCCT 1e-7 in mm). SMLib and Creo-relative scale with model size. Every receiving CAD system ships an absolute 0.01 mm modelling tolerance today, and Creo moved from relative to absolute in 7.0 and recommends absolute for multibody. Relative tolerances make the same feature valid in one part and invalid in a larger assembly.
- Local feature size. Jackson: an edge shorter than 4 x tol cannot be split cleanly. Rhino: features >= 10 x tolerance. Parasolid: short edge = entirely inside vertex bubbles (an error). These all tie tolerance to the smallest local feature, not to the bounding box.
- Margins. OCCT 1.00001 (exact distance computation) and 1.05 (sampled deviation); SMLib 2.0 (edge) and 2 x max (vertex); NX users report edge tolerances at 2 x to 4-5 x modelling tolerance [read-excerpt]. None is derived from an error bound in the sources.
- Growth through chained ops. Coincidence is "within sum of tolerances" (Jackson, SMLib XSect = Zone1 + Zone2), so a new vertex or edge made from two tolerant edges can need t1 + t2. With a 2 x margin at each re-validation, n steps can grow tolerances like 2^n unless re-measured. Parasolid's answer is PK_EDGE_optimise (re-measure and lower, with re-fit). Jackson asks for containment but gives no rule.
- Cost of tight fitting. For cubic spline approximation, error ~ h^4, so segment count ~ tol^(-1/4): 1e-6 -> 1e-7 is about 1.8 x more segments, 1e-3 -> 1e-7 about 10 x [unread/memory, standard approximation theory]. Parasolid and OCCT docs confirm "costly" without numbers [read].

## 5. Can we improve on SMLib?

1. Replace bbox/800 with an absolute cap plus a local cap. For a 1000 mm part, /800 allows 1.25 mm; for a 5 mm part, 0.006 mm. That is scale-dependent and much looser than any receiver's 0.01 mm on big parts. Proposal: cap = min(C_abs, L_local / 10), with L_local = shortest edge in the entity's neighbourhood or narrowest adjacent face width; hard error if tol > L_local / 4 (Jackson).
2. Replace fixed 2 x with measured gap + measurement error bound, times 1.05 to 1.1. Measure the max gap by sampling plus local maximisation (Newton on the gap function at the best samples). The bound delta is the change in max between the last two refinements. tol = 1.1 x (gap + delta), floored at the class floor. This keeps OCCT's 1.05 style tightness while being honest about sampling error, and stops 2 x per-step growth.
3. Separate classes: own exact, own fitted, good import, poor import, each with its own floor and cap (below).
4. After every operation, run an optimise pass: re-measure, lower to measured + margin, re-fit the pcurve first (Parasolid warns that lowering without re-fitting breaks checking).
5. Vertex: max(incident edge tols, 1.1 x measured vertex gaps), not 2 x.
6. Write STEP uncertainty as the max entity tolerance (not the average, which OCCT writes by default), so receivers never see a declared tolerance smaller than real gaps.

On 1e-7 mm for own fitted geometry: it is 50 x below Parasolid's exact-edge threshold (5e-6 mm) and 10 x below ACIS resabs (1e-6 mm with mm units), so no receiver can see the difference from 1e-6. It costs roughly 1.8 x more spline segments than 1e-6 (theory), and it equals OCCT Confusion, so measured gaps sit at the noise floor of iterative routines (OCCT Intersection 1e-9). No interop failure from being tight was found; the risk is only size and speed.

## 6. Proposed values (mm)

| Class | Value | Evidence | Confidence |
|---|---|---|---|
| Session resolution eps (coincidence; implicit tol of exact entities) | 1e-7, size box about 1e4 mm | OCCT Confusion 1e-7; Parasolid ratio 1e11 and ACIS 1e10 give 1e-7 to 1e-6 at a 1e4 box | medium |
| Own exact builds | no stored tol; behave as eps (Parasolid gives exact entities half session precision) | Parasolid, OCCT | high |
| Own fitted builds | fit target 1e-6; stored tol = 1.1 x (gap + delta), floor 1e-6, cap 1e-5 (Parasolid session precision) | OCCT Approximation 1e-6; below Parasolid exact threshold 5e-6 and ACIS resabs 1e-6, so edges arrive exact in Parasolid/ACIS receivers; 1e-7 costs ~1.8 x segments for no visible gain | medium-high |
| Good imports | file uncertainty, floor eps, per-edge measured; auto cap 1e-3 | Parasolid separation 100 x precision = 1e-3 mm; CATIA 0.001 merging; Parasolid SP-curve cost band 1e-5..1e-3 | medium |
| Poor imports | auto cap 0.01; up to 0.1 only by explicit user opt-in; above that, refit or report | NX, Creo, Rhino all default 0.01; Solid Edge stitches at 0.01; CATIA join/heal max 0.1; OCCT default max 1.0 is too loose | high for 0.01, medium for 0.1 |
| Margin | 1.1 x (measured max gap + measurement error bound) | OCCT 1.05 / 1.00001; SMLib 2 x has no stated basis and compounds | medium |
| Local cap | warn if tol > L_local/10, error if > L_local/4 | Rhino features >= 10 x tol; Jackson 4 x tol split limit | medium |
| Growth cap | per op: re-measure and lower afterwards; never raise above class cap; flag any op that raises a tol by > 10 x its inputs' max | Jackson 7.3 (qualitative); Parasolid PK_EDGE_optimise; patent "< 10 eps as reliable as exact" | low |
| Short edge | error if whole edge inside vertex bubbles; error if length < 4 x max(vertex tol); PDQ warning < 0.1 mm, strong < 0.01 mm | Parasolid short-edge definition; Jackson; VDA tiny-edge rationale; receivers collapse below 0.01; CADInterop 0.1 mesher threshold | medium (validity), low (PDQ numbers) |
