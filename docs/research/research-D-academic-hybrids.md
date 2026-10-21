# Research D: academic hybrids of B-rep and F-rep / implicit / SDF

Date of search: 2026-10-21. About 52 tool calls used (over the 40 budget). Gaps marked [GAP].
Tags: [SRC] = stated by the cited source (or its abstract as relayed by search). [UNVERIFIED] = from search-tool summary or general knowledge, not read in the primary text.

---

## 1. Hybrid models

### 1.1 Pasko, Adzhiev, Fryazinov and colleagues

- **HFRep, Hybrid Function Representation** - Tereshin, Pasko, Fryazinov, Adzhiev, "Hybrid function representation for heterogeneous objects", Graphical Models 114 (2021) 101098, DOI 10.1016/j.gmod.2021.101098; preprint https://arxiv.org/abs/2012.15176 (Dec 2020). [SRC] A "unifying functionally-based hybrid representation called HFRep" that gives "a continuous smooth distance field in Euclidean space". Covers paired combinations of FRep with SDF, ADF (adaptive distance fields) and IDF (interpolated distance fields), generating geometry and attributes. It is FRep + sampled distance field, not B-rep + FRep.
  - Earlier short paper: "Hybrid Function Representation with Distance Properties", Eurographics 2019 short papers, 17-20. http://eprints.bournemouth.ac.uk/33493/
  - PhD thesis: Tereshin, "Hybrid Modelling of Heterogeneous Volumetric Objects", Bournemouth 2020. https://eprints.bournemouth.ac.uk/35542/
  - Code: https://github.com/teshaTe/HybridFrep
- **BSP-fields (exact polygonal B-rep to FRep)** - Fryazinov, Pasko, Adzhiev, "BSP-fields: An exact representation of polygonal objects by differentiable scalar fields based on binary space partitioning", Computer-Aided Design 43(3) (2011) 265-277, DOI 10.1016/j.cad.2010.11.011. https://eprints.bournemouth.ac.uk/18561/ [SRC] Requirements: function is zero exactly on the boundary, "no extra zero-value isosurfaces", C1 everywhere. Uses BSP partitioning by face planes; independent of genus and number of components. Keywords include "Exact conversion". Compares min/max, R-functions and SARDF for combining the BSP halfspaces. Earlier tech report TR-NCCA-2008-03 (2008). This is the closest Pasko-group result to "B-rep as Boolean tree of per-face implicits", restricted to polygonal faces.
- **Constructive hypervolumes** - Pasko, Adzhiev, Schmitt, Schlick, "Constructive Hypervolume Modeling", Graphical Models 63(6) (2001) 413-442, DOI 10.1006/gmod.2001.0560. [SRC] FRep for geometry and attributes "represented independently using real-valued scalar functions"; each function "evaluated at the given point by a procedure traversing a constructive tree structure". HyperFun is the modelling language.
- **Hybrid cellular-functional (heterogeneous objects)** - Adzhiev, Kartasheva, Kunii, Pasko, Schmitt, "Hybrid cellular-functional modeling of heterogeneous objects", ASME JCISE 2 (2002, issue printed 2003) 312-321, DOI 10.1115/1.1559580; conference version SMA '02. [SRC] "hybrid model of geometry and attributes combining a cellular representation and a functionally based constructive representation of dimensionally non-homogeneous entities". Example: "conversion of a 3D implicit complex to the cellular representation".
- **SARDF** - Fayolle, Pasko, Schmitt, Mirenkov, "Constructive heterogeneous object modeling using signed approximate real distance functions", ASME JCISE 6(3) (2006) 221-229, DOI 10.1115/1.2218366. [SRC] "smooth approximation of the min/max operations ... for maintaining an approximate signed distance function in constructive shape modeling". LNCS 4889 chapter (2008), DOI 10.1007/978-3-540-68443-5_5, proves C1 continuity. Code: https://github.com/fayolle/sardf
- **Bounded blending** - G. Pasko, A. Pasko, T. Kunii, "Bounded blending for function-based shape modeling", IEEE CG&A 25(2) (2005) 36-45, DOI 10.1109/MCG.2005.37. [SRC per search summary] Blends built from R-functions and displacement functions with a localised area of influence set by a bounding solid; supports blends on single edges, partial edge blends, recursive blends.
- **HyperFun** - https://en.wikipedia.org/wiki/HyperFun ; Pasko et al., "Function-based shape modeling: mathematical framework and specialized language" https://fab.cba.mit.edu/classes/S62.12/docs/Pasko_frep.pdf [SRC] tools include a polygonizer and a POV-Ray plug-in; HyperFun compiles to C/C++/Java. GraphiCon 2000 invited talk said the team "planned to develop an advanced computer-aided design system based on several geometric representations including F-rep and voxels" https://www.graphicon.ru/html/2000/Reports_Invited_by_IPC/Pasko_Adzhiev_Savchenko.pdf
- **Commercial CAD plug-in: Symvol for Rhino (Uformia, 2011)** - [SRC] Skoltech lab page lists "Symvol for Rhino, a plug-in with FRep to a well-known software package" and MeshUp https://sci.skoltech.ru/page16523582.html . An FRep slicer paper says Symvol "uses F-rep as a flexible backend for computational solid geometries" https://www.sciencedirect.com/science/article/abs/pii/S0167839618300268 . Vilbrandt (Uformia CTO): "We are not a voxel modeler or a blob-modeler, or a CSG modeler, or surface modeler (however our underlying platform can convert meshs and voxels)" https://www.fabbaloo.com/2011/11/uformia-saving-the-world-from-polygons-html (2011). [GAP] No source found describing how Symvol consumed Rhino B-rep geometry (whether exact face-by-face conversion or mesh conversion). No Maya or SolidWorks FRep plug-in confirmed.

### 1.2 Vadim Shapiro and colleagues

- **R-functions survey** - Shapiro, "Semi-analytic geometry with R-functions", Acta Numerica 16 (2007) 239-303, DOI 10.1017/S096249290631001X. [SRC] Rvachev called R-functions "logically charged functions" because they "encode complete logical information within the standard setting of real analysis".
- **Real functions for solids** - Shapiro, "Real functions for representation of rigid solids", CAGD 11(2) (1994) 153-175. https://www.sciencedirect.com/science/article/abs/pii/0167839694900302 [UNVERIFIED content].
- **Implicit functions from B-rep and CSG, SAGE** - Shapiro & Tsukanov, "Implicit Functions With Guaranteed Differential Properties", SMA '99, DOI 10.1145/304012.304038. [SRC] "R-functions provide the methodology for constructing exact implicit functions for any semi-analytic set"; implemented "in the SAGE (Semi-Analytic Geometry Engine)" at UW Madison, and "tested composite functions built automatically from boundary and CSG representations". https://www.researchgate.net/publication/2514675
- **SAGE architecture** - Shapiro & Tsukanov, "The Architecture of SAGE - A Meshfree System Based on RFM", Engineering with Computers 18 (2002) 295-311, DOI 10.1007/s003660200027. [SRC] meshfree system where geometry "neither conforms to, nor is restricted by a spatial discretization"; based on Rvachev's Function Method.
- **Approximate distance fields** - Biswas & Shapiro, "Approximate distance fields with non-vanishing gradients", Graphical Models 66(3) (2004) 133-159, DOI 10.1016/j.gmod.2004.01.003. [SRC] exact Euclidean distance fields "are expensive to compute and may not be sufficiently smooth"; quality of an approximate field is characterised locally near the boundary by "its order of normalization". Companion: Biswas, Shapiro, Tsukanov, "Heterogeneous material modeling with distance fields", CAGD 21(3) (2004) 215-242. Also Freytag, Shapiro, Tsukanov, "Field modeling with sampled distances", CAD (year not confirmed) [GAP].
- **Representation-agnostic / query-centric** - Hoffmann, Shapiro, Srinivasan, "Geometric interoperability via queries", CAD 46 (2014) 148-159, DOI 10.1016/j.cad.2013.08.027, free PDF https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=914283 . [SRC] "we advocate a query-centric approach"; queries are "formally specified computable functions that are widely supported by different systems"; criticises data-centric exchange for "lack of consistent supporting mathematical models".
- **Intact Solutions** (Shapiro founder/CTO) - [SRC] Immersed Method of Moments described as "a geometry-agnostic FEA methodology" working on "B-reps, STL/PLY, voxels, implicit, g-code" https://intact-solutions.com/immersed-method-of-moments-imm/ . No published "Intact geometry kernel" found. Earlier product Scan&Solve (Rhino plug-in), 2013 interview https://lifeupfrontblog.wordpress.com/2013/06/18/scan-and-solve-intact-solutions-vadim-shapiro/

### 1.3 Convergent modelling

- [SRC] Siemens Parasolid fact sheet: B-rep "extended to include facet mesh as a new surface type"; facet meshes "fully integrated into the core model and not a parallel representation" https://www.plm.automation.siemens.com/media/global/de/siemens_parasolid_with_convergent_modeling_fact_sheet_tcm53-70396.pdf (NX 11, ~2016). Convergent modelling is mesh + B-rep, not implicit + B-rep.
- [GAP] No peer-reviewed academic paper titled or framed as "convergent modelling" of B-rep with implicit/F-rep was found.

### 1.4 Systems where a construction tree is evaluated into both B-rep and field

- [SRC] Hoffmann & Shapiro, "Solid Modeling", ch. 57 in Handbook of Discrete and Computational Geometry, 3rd ed., CRC 2017 (preprint https://www.csun.edu/~ctoth/Handbook/chap57.pdf): PADL "evaluated into an internal format, usually CSG expressions plus a Brep evaluation for efficient rendering"; "CSG expressions may be used as implicit representations [Bloomenthal & Bajaj 1997]"; modern constructions "are usually evaluated into Breps, mesh representations, or spatial discretizations".
- SAGE built implicit functions automatically from B-rep and CSG (see 1.2).
- NH-Rep builds a Boolean tree from a B-rep and evaluates it as a field (see 2.3).

---

## 2. B-rep to SDF

### 2.1 Analytic faces

- [SRC] Point-to-circle distance in 3D has a coordinate-free closed form (project P onto circle plane, closest point on ray from centre) - Eberly, "Distance to Circles in 3D", Geometric Tools https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf ; code https://www.geometrictools.com/GTE/Mathematics/DistPoint3Circle3.h . Degenerate case: point on axis, "entire circle is equidistant".
- [UNVERIFIED, follows from the above] Distance to a full torus = |distance to core circle - tube radius|; closed form. Plane, sphere, infinite cylinder are textbook closed forms. Cone distance reduces to a 2D problem in the meridian plane. [GAP] Primary source for cone/frustum not retrieved.
- [GAP] Quilez's primitive list (iquilezles.org/articles/distfunctions) returned 404 / redirect; not read.
- [UNVERIFIED] For a trimmed face, distance to the face is not the distance to its carrier surface; edges and vertices must also be tested. No primary source retrieved stating this.

### 2.2 NURBS faces

- [SRC] Ma & Hewitt, "Point inversion and projection for NURBS curve and surface: Control polygon approach", CAGD 20(2) (2003) 79-99, DOI 10.1016/S0167-8396(03)00021-9: subdivide into Bezier patches, choose candidates from control net, then iterate. A later paper notes control-polygon elimination "may exclude the correct solution"; improved by Selimovic, CAGD 23(5) (2006) 439-445.
- [SRC per search summary] Upreti/Subbarayan algebraic level sets on low-degree NURBS (CAD 2016 https://www.sciencedirect.com/science/article/abs/pii/S0010448516301129 ; CAD&A 10(3) 2013 https://www.cad-journal.net/files/vol_10/CAD_10(3)_2013_427-443.pdf) avoid iterative projection but give approximate signed distance measures, not Euclidean distance.
- [SRC per search summary] Krishnamurthy et al. GPU minimum-distance / closest-point on trimmed NURBS with bounding-box hierarchies.
- Closed form for general NURBS: none found. Projection is numeric (subdivision plus Newton).

### 2.3 Sign (inside/outside)

- [SRC] Spainhour & Weiss, "Robust Containment Queries over Collections of Trimmed NURBS Surfaces via Generalized Winding Numbers", ACM TOG 45(3) art. 26 (2026), DOI 10.1145/3797957, arXiv 2504.11435: described as the first algorithm to compute GWNs "over unstructured collections of trimmed NURBS patches"; robust to non-watertightness; uses Stokes' theorem to a boundary integral plus line-surface intersection when near; trims handled by 2D GWN in parameter space; "accurate to arbitrary precision at arbitrary distances". Curve version: Spainhour, Gunderman, Weiss, TOG 43(4) (2024), DOI 10.1145/3658228.
- Follow-ups (2026 arXiv): "Spatially Accelerated Winding Numbers for Curved Geometry" arXiv 2605.19200; "The Antipodal Method" arXiv 2605.01536. Closed-form GWN for rational curves: TOG 44(4) (2025), DOI 10.1145/3730886.

### 2.4 Composition: min/max is not an exact SDF

- [SRC] Quilez, "Interior Distance" https://iquilezles.org/articles/interiordistance/ (announced 2020): min/max of SDFs give "just a lower bound to the real distance"; union's interior is "totally broken", not a small deviation.
- [SRC] Marschner, Sellan, Liu, Jacobson, "Constructive Solid Geometry on Neural Signed Distance Fields", SIGGRAPH Asia 2023, DOI 10.1145/3610548.3618170: booleans by min/max produce "Pseudo-SDFs", "eikonal yet not true distance functions"; proposes a closest-point loss to train toward exact SDF.
- [SRC] Curv docs: offset of an exact field gives rounded offset; smooth blends and chamfer are "approximate" https://github.com/curv3d/curv/blob/master/docs/shapes/Distance_Field_Operations.rst

### 2.5 B-rep as Boolean tree of per-face implicits

- **NH-Rep** - Guo, Liu, Pan, Guo, "Implicit Conversion of Manifold B-Rep Solids by Neural Halfspace Representation", ACM TOG 41(6) art. 276 (SIGGRAPH Asia 2022), DOI 10.1145/3550454.3555502, arXiv 2209.10191, code https://github.com/guohaoxiang/NH-Rep . [SRC] "a Boolean tree built on a set of implicit functions represented by the neural network"; "efficient algorithm to extract the Boolean tree from a manifold B-Rep solid"; each leaf is a B-rep (sub)patch with a neural implicit; root "approximates the signed distance field". Tested on ~10,000 models with NURBS. Limitations: halfspace functions "hard to edit", slow learning, no non-manifold B-reps. [GAP] Abstract does not say how the tree is chosen (convex/concave edge rule not confirmed from text read).
- **BSP-fields** (1.1) is the exact, non-neural analogue for polygonal B-reps (2011).
- **Boundary-sampled halfspaces** - Du, Zhou, Carr, Ju, ACM TOG 40(4) art. 53 (SIGGRAPH 2021), DOI 10.1145/3450626.3459870. [SRC] shape defined "by sparsely placed samples on the boundary of each halfspace" instead of Boolean expressions; halfspaces are primitives or free-form implicits; released code extracts boundary with Marching Cubes on a grid.
- **BR-DF** - arXiv 2511.14870 (2025) https://arxiv.org/pdf/2511.14870 [SRC per summary] SDF plus per-face UDFs, "equivalent to and mutually convertible with a faceted B-Rep model" (segmented triangle mesh, not exact surfaces).
- **DualBrep** - arXiv 2606.31579 (2026) [SRC per summary] SDF for geometry plus UDF with Voronoi partition for topology; neural rebuilder outputs watertight B-rep; evaluated on ABC.
- Jones et al., "Self-Supervised Representation Learning for CAD", CVPR 2023: per-face parametric surface with face boundary as an SDF in (u,v) https://openaccess.thecvf.com/content/CVPR2023/papers/Jones_Self-Supervised_Representation_Learning_for_CAD_CVPR_2023_paper.pdf

---

## 3. SDF / implicit to B-rep

### 3.1 Implicitisation and parametrisation theory

- [SRC] Sederberg, Anderson, Goldman, "Implicit representation of parametric curves and surfaces", CVGIP 28(1) (1984) 72-84, DOI 10.1016/0734-189X(84)90140-3: implicitisation of rational parametric curves/surfaces has closed-form solutions using only +,-,*,/; a rational "bi-m-ic" patch has implicit degree 2m^2 (bicubic: 18). Direction parametric -> implicit always works for rational parametric input.
- [SRC] Schicho, "Rational Parametrization of Surfaces", J. Symbolic Computation 26(1) (1998) 1-29: algorithm "decides if such a parametric representation exists, based on Castelnuovo's rationality criterion" (arithmetic genus and second plurigenus both zero); computes it via adjoints if yes. So implicit -> rational parametric is possible only for rational surfaces, and decidable. Real case: Schicho, ISSAC '98, 302-308.
- [SRC] Pratt, "Cyclides in computer aided geometric design", CAGD 7 (1990) 221-242: Dupin cyclides "generalise ... the plane, cylinder, cone, sphere and torus", with "a rational biquadratic representation for cyclide surface patches". So these surfaces admit exact rational (NURBS) patches.
- [SRC] Offsets: Pottmann & Peternell (Laguerre geometry, CAGD 15 (1998) 223-249, https://www.geometrie.tuwien.ac.at/peternell/laggeo2.pdf): polynomial/rational curves "in general not possess rational offsets", so offsets "have to be approximated in rational B-spline form". Farouki & Sakkalis, "Pythagorean hodographs", IBM J. Res. Dev. 34(5) (1990) 736-752: PH curves are the class with rational offsets.
- [UNVERIFIED] Generic smooth quartic surfaces are not rational (K3); arbitrary non-algebraic fields have no exact NURBS form. Not checked against a primary source in this run.

### 3.2 Isosurface extraction with sharp features

- [SRC] Ju, Losasso, Schaefer, Warren, "Dual contouring of Hermite data", SIGGRAPH 2002, DOI 10.1145/566570.566586: uses "exact intersection points and normals" on grid edges; QEF minimiser per cell; avoids explicit feature identification; octree simplification.
- Kobbelt, Botsch, Schwanecke, Seidel, Extended Marching Cubes, SIGGRAPH 2001 [UNVERIFIED citation details]. Later work notes EMC, DC and dual-grid MC "preserve sharp features" but "can potentially generate intersecting faces".
- [SRC] Plantinga & Vegter, "Isotopic Approximation of Implicit Curves and Surfaces", SGP 2004, DOI 10.2312/SGP/SGP04/251-260, https://www.math.rug.nl/~gert/documents/2004/pv-sgp-2004.pdf : interval arithmetic plus octree gives an output "regularly isotopic" to the true surface; for non-singular surfaces. Topology guarantee for a mesh, not an exact B-rep.
- [SRC] Keeter, "Massively Parallel Rendering of Complex Closed-Form Implicit Surfaces", TOG 39(4) 141 (2020), DOI 10.1145/3386569.3392429: interval arithmetic to skip regions and prune expressions; only needs C0. libfive: https://libfive.com

### 3.3 B-rep / CSG recovery

- [SRC] Shapiro & Vossler, "Separation for boundary to CSG conversion", ACM TOG 12(1) (1993) 35-55, DOI 10.1145/169728.169723: needs "separating half-spaces" that do not contribute to the boundary; linear separators suffice for planar-edged solids and for solids bounded by general quadrics; general case "remains poorly understood". Exact, not learned. Also Shapiro & Vossler, "Construction and optimization of CSG representations", CAD 23(1) (1991) 4-20.
- [SRC] Hoffmann & Shapiro 2017 handbook: B-rep -> CSG for curved faces "is mostly solved for second-degree surfaces, but for higher degree surfaces it is open". Example: a solid bounded by two spherical faces and one cylindrical face cannot be written using only those three halfspaces; two extra separators D, E are needed; optimised result A + CDE + B.
- [SRC] Buchele & Crawford, "Three-dimensional halfspace CSG tree construction from implicit boundary representations", CAD 36(11) (2004) 1063-1073, DOI 10.1016/j.cad.2004.01.006 (SM '03 DOI 10.1145/781606.781629): input "a valid boundary representation ... of piecewise implicit surfaces", output halfspace CSG; disjoint decomposition makes minimisation "generally not necessary".
- [SRC] Du et al., "InverseCSG", ACM TOG 37(6) 213 (2018), DOI 10.1145/3272127.3275006: inverse CSG as program synthesis; finds "simple programs which are close to the ground truth" (approximate, from meshes).
- CSGNet: Sharma et al., CVPR 2018 [venue per Semantic Scholar listing; details not read].
- [SRC] CAPRI-Net: Yu et al., CVPR 2022, 11768-11778: quadric primitives assembled by CSG in fixed order; excludes "cones and tori". Follow-up D2CSG, NeurIPS 2023.
- [SRC] Friedrich & Fayolle, "On the Complexity of the CSG Tree Extraction Problem", arXiv 2103.06139 (2021): about the search-space complexity given a point cloud and fitted primitives [specific result not read].
- [SRC] Point2CAD, Liu et al., CVPR 2024, arXiv 2312.04962: segment, fit primitives and neural freeform surfaces, intersect neighbours for edges and corners; claims "topologically correct reconstructions"; says ComplexGen outputs "are often invalid, in the sense that the predicted geometry and topology are inconsistent". ComplexGen: Guo et al. 2022, arXiv 2205.14573.
- **Tolerance guarantee**: [GAP / negative finding] No method found that takes an arbitrary field and outputs an exact B-rep with a proven distance tolerance. Guarantees found are: exact B-rep -> CSG for quadric/planar B-reps (Shapiro-Vossler 1993, Buchele-Crawford 2004), topological isotopy for meshes (Plantinga-Vegter 2004), and tolerance-controlled approximation of rolling-ball blends (see 4). Learned methods (CSGNet, CAPRI-Net, InverseCSG, Point2CAD, DualBrep) report empirical metrics only.

---

## 4. Fillets and blends with true metric radius

### 4.1 Definition via offsets

- [SRC] Rossignac & Requicha, "Constant-Radius Blending in Solid Modeling", Computers in Mechanical Engineering 3(1) (1984) 65-73: rolling-sphere blends "defined mathematically in terms of offset solids, and represented computationally in an extended form of constructive solid geometry that includes offset operators"; blend surfaces approximated with piecewise circular curves. Also "Offsetting operations in solid modelling", CAGD 3 (1986) 129-148.
- [SRC] Hoffmann & Shapiro 2017 handbook: ball centre "lies on the intersection of the offsets, by r, of both f and g"; round = locally convex blend, fillet = locally concave.
- [SRC, GitHub issue, low authority] rolling-ball fillet = morphological closing (concave) / opening (convex) by a ball and needs true distance; 1-Lipschitz bounds do not suffice https://github.com/bmander/geomsolver/issues/66

### 4.2 Implicit blends and their radius behaviour

- [SRC] Hoffmann & Hopcroft, "The Potential Method for Blending Surfaces and Corners", in Farin (ed.), Geometric Modeling: Algorithms and New Trends, SIAM 1987, 347-365. Rockwood & Owen, "Blending Surfaces in Solid Modeling" (1987). Rockwood, "The displacement method for implicit blending surfaces in solid models", ACM TOG (1989), DOI 10.1145/77269.77271.
- [SRC] Curv docs: Hoffmann-Hopcroft with a=b=r, lambda=0 is "The Elliptic Blend"; "At 90 degrees the fillet is a quarter-circle of radius r", at other angles "the fillet deforms to an ellipse"; "might be bad for engineering, if you need a constant radius fillet". Smooth blends listed as "approximate" distance.
- Bounded blending (Pasko 2005) and SARDF (2006) are R-function blends; [GAP] neither source read claims a metric rolling-ball radius.
- [GAP] No academic paper found giving an angle-independent exact rolling-ball fillet operator directly on SDFs.

### 4.3 Exact B-rep counterparts

- [SRC] Peternell & Pottmann, "Computing Rational Parametrizations of Canal Surfaces", J. Symbolic Computation 23 (1997) 255-266, DOI 10.1006/jsco.1996.0087, https://www.dmg.tuwien.ac.at/geom/ig/peternell/canalsurf.pdf : "any canal surface to a rational spine curve m(t) and a rational radius function r(t) possesses rational parametrizations".
- [SRC] Dahl & Krasauskas, "Rational fixed radius rolling ball blends between natural quadrics", CAGD 29(9) (2012) 691-706, DOI 10.1016/j.cagd.2012.06.008: classifies all configurations with exact rational fixed-radius rolling-ball blends; highlights list four pairs: "cone/plane, two cones touching at one or two points, and cone/sphere touching at one point", with minimal bi-degree parametrisations. [UNCERTAIN] Whether trivial rotational cases (e.g. plane perpendicular to cylinder, giving a torus patch) are excluded from that list was not confirmed; the abstract page returned 403.
- [SRC] Karciauskas & Krasauskas, "Rational rolling ball blending of natural quadrics", Math. Modelling and Analysis 5(1) (2000) 97-107: "apart from simple cases" fixed-radius blends "usually turn out to be irrational" and "cannot be exactly represented as Bezier or B-spline surfaces"; uses variable radius to get rationality.
- [SRC] Farouki-style approximation, "Approximation of rolling-ball blends for free-form parametric surfaces", CAD (1996) https://www.sciencedirect.com/science/article/abs/pii/0010448596000085 : G1 sequences of rational patches to "a prescribed precision", tracing the offset-intersection spine; isoparametric curves "only approximately circular".
- [SRC] Rational parametrisation condition (Lu & Pottmann 1996, as cited): needs rational alpha, beta with |C'(t)|^2 = alpha^2 + beta^2.
- [SRC] Cyclides: Pratt 1990 gives cylinder-plane cyclide blend ("only the pillar type of blend is possible" per Garnier et al.); Garnier et al., "Blending Planes and Canal Surfaces Using Dupin Cyclides" https://www.researchgate.net/publication/221014325 ; Allen & Dutta, "Cyclides in pure blending I", CAGD 1997; Shene 1998 cone-cone via offset trick "does not always work properly".

---

## 5. Recipe-as-master kernels with several evaluators

- [SRC] Requicha, "Representations for Rigid Solids: Theory, Methods, and Systems", ACM Computing Surveys 12(4) (1980) 437-464, DOI 10.1145/356827.356833: framework for representation schemes; Section 1.4 formal properties (validity, unambiguity etc.). [UNVERIFIED] exact wording not read.
- [SRC] Shapiro, "Maintenance of Geometric Representations through Space Decompositions", IJCGA 7 (1997) 21-56 (errata/reprint 7(4) 383-418, DOI 10.1142/S0218195997000247): space decompositions as a graph ("Finite Discrete Lattice"), nodes are parts of space, links are inclusion.
- [SRC] Raghothama & Shapiro, "Consistent updates in dual representation systems", SMA '99 (DOI 10.1145/304012.304019), CAD 32(8-9) (2000) 463-477, DOI 10.1016/S0010-4485(00)00036-1: parametric and B-rep "remain consistent if they continue to model the same solid"; gives "necessary and sufficient conditions for consistency" of B-rep and CSG; two algorithms using "sign-invariant space decompositions"; also applicable to "feature attachments". This is the main consistency-check result found.
- [SRC] Related: Raghothama & Shapiro, "Boundary representation deformation in parametric solid modeling", ACM TOG 17(4) (1998); "Topological framework for part families", JCISE 2(4).
- [SRC] Hoffmann & Shapiro 2017 handbook: feature/constraint interfaces need infrastructure "capable of supporting and maintaining consistency of both types of representations within a single system [RS00]"; conversions follow a "generate and test" paradigm (generate candidates, test against source, assemble); conversions split into "deterministic evaluations" and "nondeterministic comprehensions" [HSS11]; open problem: maintaining multiple views "is a generalization of a representation conversion problem where multiple representations of the same shape must be maintained simultaneously".
- [SRC] Hoffmann, Shapiro, Srinivasan 2014: query-centric interoperability (representations interoperate through shared queries, not translation).
- Related CSG-equivalence check: "Solving CSG equations for checking equivalency between two different geometric models", CAD (2003/04) https://www.sciencedirect.com/science/article/abs/pii/S0010448503002185 [content not read].
- [GAP] No recent (2020-2026) academic proposal found for a kernel whose master is a feature/CSG tree with simultaneous exact-B-rep, field and mesh evaluators and explicit cross-evaluator consistency checks. Closest: SAGE (B-rep/CSG -> R-function field, 1999-2002), Raghothama-Shapiro dual-representation consistency (1999/2000), NH-Rep and BR-DF/DualBrep (learned conversions).

---

## Not found or uncertain

1. Symvol for Rhino internals: whether it converted Rhino B-reps face-by-face into FRep or via meshes.
2. Academic "convergent modelling" paper combining B-rep with implicits (only Siemens mesh+B-rep material found).
3. NH-Rep's exact tree-extraction rule (convex/concave edge classification) not read from paper text.
4. Dahl & Krasauskas 2012: whether rotational trivial cases (torus fillets for perpendicular plane-cylinder) are outside the four-pair list.
5. Exact angle-independent rolling-ball blend operator on SDFs: no paper found.
6. Any method giving an exact B-rep from an arbitrary field with a proven distance tolerance: none found.
7. Primary sources for point-to-cone closed form and trimmed-face edge/vertex distance handling not retrieved.
8. Friedrich & Fayolle 2021 complexity result not read.
