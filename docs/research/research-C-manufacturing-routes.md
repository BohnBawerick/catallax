# Research C: routes from implicit / SDF models to manufacturing

Compiled 2026-10-21. Every claim carries its source. "Vendor claim" marks statements from a vendor, reseller or vendor-authored article. "Unverified" marks items the sources did not confirm. Search-result summaries were cross-checked against fetched pages where marked "(fetched)".

Key question per route: does a circle come out as a true circle (exact cylinder face), what tolerance is held, which tools accept the output.

---

## 1. 3D printing

### 1.1 STL and 3MF meshes (the default route)

- A mesh never holds a true circle. Chord height (chordal deviation) is "the maximum distance that your software will allow between the surface of the original 3D model and the surface of the STL file." Markforged recommends 0.1 mm chord, 1 deg angular. https://markforged.com/resources/blog/how-to-create-high-quality-stl-files-for-3d-prints
- Rule of thumb cited by BigRep: chord deviation about 1/20 of layer thickness, not below 1 micron; example 0.02 mm at 0.4 mm layers. https://bigrep.com/posts/stl-file-best-practices/
- Other published recommendations: 0.013-0.025 mm and about 5 deg for FDM/SLS (projects-house) http://www.projects-house.com/articles/3d-printing/stl-file-preparation-mistakes/ ; 0.01-0.001 mm (CMAC) https://www.cmac.com.au/blog/5-vital-things-about-stl-file-format-3d-printing ; 0.001 mm and 15 deg (FacFox) https://facfox.com/docs/kb/3d-printing-stl-files-a-step-by-step-guide
- CAD export defaults vs recommended print values (Inn3D Studio): SOLIDWORKS default about 0.15 mm / 30 deg, recommended 0.02 mm / 5 deg; Creo 0.03 mm / 1 deg; NX 0.02 mm (default 0.08 mm). https://inn3dstudio.com/en/blogs/%E6%9C%80%E6%96%B0%E6%B6%88%E6%81%AF/%E5%A6%82%E4%BD%95%E5%9C%A8-cad-%E4%B8%AD%E6%AD%A3%E7%A2%BA%E5%BB%BA%E7%AB%8B-3d-%E5%88%97%E5%8D%B0%E7%94%A8-stl-%E6%AA%94
- Halving the tolerance roughly quadruples face count on a curved surface; a coarse export can turn a 6.35 mm (0.250 in) hole into a twelve-sided pocket that a dowel will not enter; an SLS case at 0.5 mm default left 0.6 mm facets, 0.02 mm made a 210 MB file, 0.05 mm chosen at 42 MB. (Summarised across the guides above; the hole and SLS examples come from projects-house.)
- Arithmetic (mine, not a source): for a circle of radius r tessellated with chord error e, the largest segment angle is 2*acos(1 - e/r). For r = 5 mm and e = 0.01 mm that is about 7.2 deg, so about 50 segments.

Typical printer accuracy (service-bureau design guides, all vendor claims):
- Xometry: MJF +/-0.3% (+/-0.2 mm), SLS +/-0.3% (+/-0.3 mm), SLA +/-0.5% (+/-0.2 mm), FDM +/-0.5% (+/-0.5 mm), Carbon DLS +/-0.1% (+/-0.1 mm), DMLS +/-0.2% (+/-0.1-0.2 mm), PolyJet +/-0.1 mm. Its MJF section contradicts its own table (+/-0.3 mm lower limit). https://xometry.pro/en/articles/3d-printing-tolerances/
- MakerVerse: SLA/DLP +/-0.1 mm or 0.2%; SLS/MJF +/-0.3 mm or 0.3%; DMLS +/-0.1-0.2 mm or 0.2%. https://www.makerverse.com/resources/3d-printing/mastering-3d-printing-tolerances-for-precision/
- RapidDirect: SLA X/Y +/-0.05-0.1 mm; SLS +/-0.3 mm or 0.3%; MJF +/-0.3 mm first 25 mm then 0.1%. https://www.rapiddirect.com/blog/3d-printing-tolerances-how-accurate-is-3d-printing/
- Sovol blog (Protolabs-style "first inch" figures, origin unconfirmed): DMLS +/-0.076 mm X/Y, SLA +/-0.051 mm X/Y. https://www.sovol3d.com/blogs/news/optimize-3d-printing-tolerance
- Read across: printer error (about 0.05-0.3 mm) usually exceeds the recommended chord error (0.001-0.05 mm). Markforged's 0.1 mm chord is the loosest value and would consume a +/-0.1 mm budget (this comparison is from the guide summary, not from one source).

### 1.2 3MF Volumetric and Implicit Extension

Status and dates:
- Volumetric Design extension first announced 2021-11-17 as v0.5c for public feedback. https://www.prnewswire.com/news-releases/3mf-consortium-announces-volumetric-design-support-grows-adoption-of-3d-printing-specification-301426769.html
- Implicit namespace presented by Jan Orend (EOS) at CDFAM NYC, 2024-10-02/03. https://3mf.io/news/2024/09/interview-with-jan-orend-on-the-volumetric-and-implicit-extension-release/ , https://cdfam.com/3mf-consortium-volumetric-and-implicit-extension-release/
- Spec document in GitHub: header "Version 1.0.0", "Status Published"; Change History lists 1.0 "First published version" dated 2025-11-13. (fetched) https://raw.githubusercontent.com/3MFConsortium/spec_volumetric/master/3MF%20Volumetric%20Extension.md
- Conflict: the 3mf.io spec table (page "Last updated 2026-09-17") still lists Volumetric Extension v0.8.0 dated 02/27/2025. (fetched) https://3mf.io/spec/ . The 3mf.io homepage lists "Volumetric v1.0", "Voxel-based and implicit (Signed Distance Function) representations", and Core "v2.0". (fetched) https://3mf.io/
- 3MF core is an ISO standard: ISO/IEC 25422:2025. https://3mf.io/

What the implicit part defines (fetched spec):
- Closed-form functions as acyclic node graphs (`<i:implicitfunction>`): arithmetic, trig, vector, matrix, function-call nodes; `mesh` node gives signed distance to a mesh, `unsignedmesh` unsigned, `beamlattice` signed distance to a beam lattice.
- `<levelset>`: the zero level of a scalar function output (input vector "pos") is the shape boundary.
- A mesh is still required: "Since fields can be evaluated in an unbounded way, a closed mesh is required to enclose any levelset element"; evaluation is limited to that mesh's extent (`meshid`), or its bounding box with `meshbboxonly`.
- Precision: functions "define a field with arbitrary precision"; single-precision IEEE 754 or better RECOMMENDED; a consumer that cannot resolve a feature below `minfeaturesize` SHOULD warn.
- So a circle defined analytically in the function stays exact in the file; the printed accuracy depends on how the consumer samples it. No tolerance is mandated by the spec.

Implementations found:
- lib3mf 2.5.0 writes levelsets; open bug (issue #472, 2026-10) says it omits the volumetric namespace from requiredextensions, making every lib3mf levelset file non-conforming to 1.0.0; fix PR #476 opened 2026-10-18. https://github.com/3MFConsortium/lib3mf/issues/472 , https://github.com/3MFConsortium/lib3mf/pull/476
- lib3mf issue #473: doubles written with 6 significant digits, so implicit constants and levelset settings cannot round-trip; names KForge as an affected producer. https://github.com/3MFConsortium/lib3mf/issues/473
- Gladius (3MF Consortium): experimental viewer/editor; imports volumetric 3MF (implicit graph or Image3D); exports 3MF volumetric, OpenVDB, STL, SVG contours, CLI contours; OpenCL compute. https://github.com/3MFConsortium/gladius/blob/main/README.md
- Go package lestrrat-go/3mf/volumetric implements a 2022/01 schema. https://pkg.go.dev/github.com/lestrrat-go/3mf/volumetric
- Spherene joined the 3MF Consortium (2026-05) and frames the Volumetric Extension as the way to carry internal geometry to slicers and printers; this is not a statement of shipping support in nTop or Fusion. https://3mf.io/news/2026/05/spherene-joins-the-3mf-consortium-to-bring-intelligent-internal-structures-to-additive-manufacturing/
- NOT FOUND: any production slicer, printer, or build processor (EOSPRINT, Magics, PrusaSlicer, Bambu, Cura, HP, Stratasys) documented as reading 3MF levelset/implicit data. The 3mf.io "Compatibility Matrix" contents were not retrieved.

### 1.3 Printing directly from an implicit without a mesh (proprietary formats)

nTop + EOS:
- nTop Implicit File (.implicit) loaded by the nTop Plugin for EOSPRINT; EOSPRINT slices it. Requires EOSPRINT 2.14+; export via "Export Implicit Body (Beta)" block introduced in nTop 4.3; announced June 2023 as public beta. Current beta status not confirmed. https://www.ntop.com/resources/product-updates/ntop-plugin-eosprint/ , https://support.ntop.com/hc/en-us/articles/17319219880467-nTop-Plugin-for-EOSPRINT-Overview , https://support.ntop.com/hc/en-us/articles/17495967403283-nTop-4-3-What-s-New
- Vendor claims: files a few MB vs GB meshes, up to 60% faster load, "no loss of accuracy"; Siemens Energy heat exchanger implicit file about 1 MB, 99% smaller, 500x faster, printed on EOS M 290. Same sources.
- SME article (2023-06-13, written by an nTop marketing manager): "no loss of geometric precision"; files "up to 10,000 times" smaller. No dimensional accuracy results reported. (fetched) https://www.sme.org/technologies/articles/2023/june/am-without-making-a-mesh/
- TCT Magazine on nTop/EOS "Implicit Interop", and stated intent to standardise via 3MF. https://www.tctmagazine.com/ntopology-and-eos-bring-implicit-interop-to-3d-printing/

nTop + Materialise:
- Partnership announced around RAPID 2024: nTop Core integrated into Magics and NxG Build Processor. https://www.voxelmatters.com/materialise-partners-with-ntop-to-accelerate-part-design-for-am/ , https://www.design-engineering.com/materialise-ntop-partner-on-additive-manufacturing-1004042943/
- Magics 2025 (RAPID+TCT 2025, apparently Magics 29): "seamless processing of nTop implicit geometries without the need for mesh conversion", paired with next-generation Build Processors. Case studies: Wartsila cylinder head, DMG MORI early access (vendor claims). https://develop3d.com/3d-printing/magics-2025-new-ntop-integration-and-more-announced-by-materialise/ , https://www.ntop.com/resources/product-updates/ntop-5/
- No source describes how the Build Processor samples the field or which machine vendors' build processors support it.

Other consumers of .implicit (vendor claims):
- Inkbit embedded the nTop Core library to read and manufacture from .implicit with no meshing. https://www.ntop.com/software/capabilities/integrations/
- Autodesk Fusion has an add-in to import implicit files (stated in nTop / trade coverage). https://www.ntop.com/software/capabilities/integrations/
- nTop also exports slice formats directly: EOS CLI, CLI/CLF, SLC and others, bypassing STL. https://www.digitalengineering247.com/article/skip-the-mesh-print-from-cad/pre-processing-and-meshing

---

## 2. CNC machining

### 2.1 CAM from meshes

- Autodesk Fusion: STL/OBJ are treated as mesh, not parametric bodies; 3D toolpaths (e.g. 3D Adaptive, Parallel) run on meshes; 2D operations do not (user report, about 9 years old). https://support.bantamtools.com/hc/en-us/articles/115001671574-Fusion-360 , https://community.carbide3d.com/t/stl-to-fusion360-cam/6127 . No current official Autodesk list of mesh-capable operations found.
- Mastercam: STL can be opened as mesh and toolpaths applied like surfaces since X3/X4 (forum); X5 added eight triangular-mesh multiaxis toolpaths that accept surface, solid and STL input (CTE Magazine). Mastercam 2022 added mesh creation from solids/faces, normals on mesh for tool-axis control, Project on meshes. Later mesh smoothing, trim, split, offset. https://www.emastercam.com/forums/topic/36722-machining-stl-files/ , https://www.ctemag.com/products/mastercam-multiaxis , https://www.mastercam.com/community/blog/5-new-mastercam-2022-features-you-need-to-know/ , https://ctemag.com/articles/meshing-cad-cam/
- Autodesk PowerMill: imports IGES, VDA and STL. Official help: "the accuracy of the toolpath produced by PowerMill is dependent on the accuracy of the model read into the program. Therefore, it is essential that the original model is produced to a suitable tolerance." Forum: Mesh.MaxTriangleSize option limits triangle length where faceting shows; PowerMill builds an internal triangle model for toolpath calculation. https://help.autodesk.com/cloudhelp/2023/ENU/PWRM-ReferenceHelp/files/GUID-D434293D-CFD5-4EEF-B146-D9881C9BB5DA.htm , https://forums.autodesk.com/t5/powermill-forum/stl-tolerance-setting/td-p/9942005
- Practical Machinist: mid-range CAM (PowerMill, Mastercam, Surfcam, VisualMill, OneCNC, Gibbs) computes toolpaths on tessellations of surfaces/solids anyway. https://www.practicalmachinist.com/forum/threads/running-toolpath-on-an-stl.136396/
- hyperMILL (OPEN MIND): hyperMILL Mesh module (formerly hyperCAD-S Mesh) prepares meshes for CAM: decimate, mesh from faces, split; 2022.1 adds "Faces from Mesh" (a plane per triangle) and mesh analysis. Imports .stl/.stla/.stlb; reseller lists 3MF, OBJ, point clouds. https://www.openmind-tech.com/en/cad/versions/polygon-meshes/ , https://www.3way.si/wp-content/uploads/2025/01/hyperMILL-brosura-2022.1.pdf , https://www.johnhart.com.au/131-software/hypermill
- OPEN MIND vendor statement against mesh-based paths (True Shape / HPM): toolpaths "calculated directly on the CAD surfaces rather than on tessellated geometry"; mesh-based toolpaths "can reproduce minor facets or geometric deviations from the mesh on the finished surface" even at small tolerances. https://www.openmind-tech.com/en/cam/high-precision-machining-hpm/
- ModuleWorks (component supplier to many CAM vendors): 3-axis mesh-based machining toolpaths (2013.4); 2026.08 "SuperFinish on NURBS Surfaces combines mesh-based calculation reliability with NURBS-based toolpath accuracy" (vendor claim, about 2026-08). https://www.moduleworks.com/moduleworks-is-a-leading-provider-of-cadcam/ , https://www.etmm-online.com/module-works-2026-08-gpu-simulation-toolpath-optimisation-a-cd3672a0e9ac739f80b43632be236777/
- 3-axis vs 5-axis on meshes: Mastercam documents multiaxis mesh toolpaths (X5). For hyperMILL no source confirmed 5-axis finishing directly on mesh. Research: 5-axis adaptive curvilinear toolpaths for STL surfaces, compared against Mastercam and NX9 (ResearchGate, about 2016). https://www.researchgate.net/publication/301741449_Five-axis_machining_of_STL_surfaces_by_adaptive_curvilinear_toolpaths ; robust toolpaths on poor-quality scan meshes (arXiv 2004.09136, 2020). https://arxiv.org/pdf/2004.09136

### 2.2 Stock models (dexel / voxel)

- ModuleWorks simulation uses a discrete dexel model; GPU version up to 10x faster than CPU (vendor benchmark, Aachen). 2024.12 adds dexel-based "graph-based triangulation" for faster, more watertight stock mesh export. https://www.moduleworks.com/precision-meets-performance-with-gpu-accelerated-simulation/ , https://www.moduleworks.com/moduleworks-releases-2024-12-software-components/
- Fusion adaptive clearing uses a voxel or octree in-process stock (third-party article). https://industrialmonitordirect.com/blogs/knowledgebase/fixing-fusion-360-cam-toolpath-failure-ram-exhaustion-and-disk
- No published accuracy figure for dexel stock resolution was found.

### 2.3 Toolpaths directly from implicit fields, SDFs or voxels (research)

- Dutta, Zhang, Liu, Chen, Wang, "Implicit Neural Field-Based Process Planning for Multi-Axis Manufacturing", arXiv 2511.17578 (v1 2025-11-15, v2 2026-04-20): layers and toolpaths represented as sinusoidal neural implicit fields, joint differentiable optimisation with collision avoidance; demonstrated on additive and subtractive examples; abstract gives no accuracy numbers or physical tests. (fetched) https://arxiv.org/abs/2511.17578
- Geodesic distance field process planning for 5-axis machining (He, Li, Li, Tang, J. Manuf. Sci. Eng. 143, 2021), cited in the above.
- Spline NC paths generated from an implicit point-set surface (PSS); CL points from the implicit surface, method works for any representation with a closest-distance query. https://www.researchgate.net/publication/277594971_Tool_Path_Generation_for_Complex_Surface_Machining_Using_Point_Cloud_Data
- Level-set offset surfaces for NC toolpaths (contour-parallel paths from equidistant level sets). https://www.researchgate.net/publication/226470205_Generation_of_offset_surface_for_tool_path_in_NC_machining_through_level_set_methods
- Additive only: OpenVCAD SDF slicing via marching squares, Wade et al. 2025, arXiv 2505.08093. https://arxiv.org/pdf/2505.08093
- NOT FOUND: any commercial CAM product that takes an implicit/SDF model as machining geometry.

### 2.4 How machine shops treat mesh vs STEP

- Shops and bureaus say send STEP, not STL; STL lacks exact surfaces, faces and edges; an STL quote request triggers a request for STEP. https://www.pcbway.com/blog/CNC_Machining/Why_the_stl_file_can_not_be_used_for_CNC_dca4b7d9.html , https://hlhrapid.com/blog/cad-file-formats-to-use-for-cnc-machining/ , https://www.simplemachining.com/blogs/cad-file-formats-explained-step-stl-iges-and-what-to-send-a-manufacturer
- Converting STL to STEP gives "a STEP wrapper around a faceted blob"; advice: send the STL honestly rather than a converted STEP. https://frankworks.ca/blog/step-file-ap242-guide-cnc-drawing-checklist (2026)
- Shops want STEP plus a 2D PDF drawing for tolerances unless AP242 with PMI. Same source and https://www.sincere-machining.com/what-files-are-needed-for-a-cnc-machining-quote/

---

## 3. Mesh or field to B-rep (commercial)

### Autodesk Fusion "Convert Mesh"
- Official help: Faceted (each triangle a face), Prismatic ("Merges groups of faces in prismatic features into singular faces"; face groups infer features), Organic (T-spline/organic body; refinement Low/Medium/High/Precise or by face count). No numeric deviation reported by the command. https://help.autodesk.com/view/fusion360/ENU/?guid=MESH-CONVERT-TO-SOLID (snapshot about 3 years old)
- Autodesk blog (about 2026-05): Prismatic for mechanical parts, commercial licences only; "Accurate" face-group generation with a Boundary tolerance parameter; face groups recognise planar, cylindrical, radial, fillet features. https://www.autodesk.com/products/fusion-360/blog/generate-face-groups-fusion/ (403 when fetched; content via search summary)
- Third-party (not Autodesk): Prismatic fits planes, cylinders, cones, spheres and tori. https://www.varsitytutors.com/practice/subjects/autodesk-fusion-360/lessons/mesh-to-solid-conversion . UNVERIFIED by Autodesk documentation whether output faces are analytic cylinder entities or NURBS.
- Open-source clone: stlToSolid writes STEP with plane/cylinder/cone/sphere/torus faces, unfitted regions stay faceted. https://github.com/Crypto69/stlToSolid

### Geomagic Design X (Hexagon)
- Auto Segment classifies mesh regions by curvature; wizards extract planes, cylinders, revolves, extrudes, fillets; "Accuracy Analyzer" gives real-time deviation between the CAD body and the mesh; two approaches: "Reflect design intent (design parameter extraction) or Generate an exact duplicate (surface fitting)" (2014 user guide). https://www.engineering.pitt.edu/contentassets/52314f399aba40fa86709314a569641c/geomagicdesignx2014userguide.pdf
- Reseller claims: fits primitives "within a tolerance you control" and reports deviation; pushes history-based parametric models to SOLIDWORKS, NX, Solid Edge, Inventor, Creo; tiers Go/Plus/Pro differ on Auto Surface and Exact Surfacing (TriMech page internally inconsistent). https://mfg.trimech.com/product/geomagic-design-x/ , https://3dmakerworld.com/products/geomagic-design-x-reverse-engineering-software , https://www.goengineer.com/blog/geomagic-design-x-features-and-tutorials

### QuickSurface
- Best-fit primitives: planes, cylinders, cones, spheres; constraints (parallel, coincident axes, perpendicular, offset); constraining for design intent "compromises with the quality of approximation of the scan points"; Analysis colour deviation map per primitive; Deviation Analyser on the whole model; STEP/IGES export or live parametric tree in SOLIDWORKS. https://www.quicksurface.com/how-to-make-3d-reverse-engineering-with-quicksurface/ , https://www.mesh2surface.com/media/D3D_OCT19_QuickSurface.pdf (DEVELOP3D review, Oct 2019)
- No source states the numeric deviation or that exported STEP faces are analytic rather than B-spline. UNVERIFIED.

### SOLIDWORKS ScanTo3D
- Surface Wizard: Automatic (organic shapes) and Guided (paint regions, then extract "analytic surfaces, such as planes, cylinders, etc."); Guided suits mostly analytic, machined parts, not complex parts. https://help.solidworks.com/2023/english/SolidWorks/scanto3d/HIDD_SURFWIZ_SOLSURF_CREATION.htm , https://help.solidworks.com/2022/english/SolidWorks/scanto3d/c_two_methods_for_converting_scan_data.htm
- Deviation Analysis appears in the semi-manual workflow, not in the Surface Wizard. https://www.mlc-cad.com/resources/solidworks/utilize-scanned-data-to-create-solid-models-using-solidworks-scanto3d/

### Siemens NX
- Convergent modelling keeps facet bodies faceted; Siemens GTAC (NX 12 era, quoted on Eng-Tips) said a convergent body "does not become analytic"; mixing facet and analytic ends faceted. https://www.eng-tips.com/threads/nx-12-0-0-27-extruded-surface-creates-a-convergent-body.445541/
- Fit Surface on facet regions creates "analytical and freeform surfaces"; Rapid Surfacing as an alternative. https://community.sw.siemens.com/s/question/0D54O000061xGpBSAU , https://saratech.com/2018/12/reverse-engineering-in-nx-12/
- No deviation readout details found for current NX releases.

### nTop implicit to STEP
- Route 1: Mesh from Implicit Body -> clean -> Quadrangulate Mesh -> optional Refine (Catmull-Clark) -> CAD Body from Quad Mesh -> STEP or Parasolid. Fewer, feature-aligned NURBS faces; only reliable on simple shapes such as topology-optimised parts; lattices fail. https://support.ntop.com/hc/en-us/articles/7323484481555-How-to-export-geometry-as-a-CAD-part-STEP-or-Parasolid , https://learn.ntop.com/courses/102-guide-to-meshing/lessons/conversion-to-cad-with-quadrangulate-mesh/
- Route 2: CAD Body from Implicit Body (Beta, nTop 3.28+): NURBS not neatly arranged, works on complex parts; "Preserved CAD" input tries to sew original CAD faces back in; example STEP over 600 MB vs .ntop under 30 kB. Same sources.
- So a cylinder that came from the implicit comes out as NURBS patches unless an original CAD face is preserved. No deviation figure published by nTop in the pages found.

---

## 4. Research: recovering exact CAD from meshes, point clouds, implicits

| Work | Input | Output faces | Exact primitives? | Tolerance guarantee? |
|---|---|---|---|---|
| Schnabel, Wahl, Klein 2007, Efficient RANSAC, CGF 26(2):214-226 | point cloud | plane, sphere, cylinder, cone, torus proxies | yes, analytic shape proxies; not a B-rep | inlier threshold epsilon only; no topology | 
| Point2CAD (Liu et al., CVPR 2024) | segmented point cloud (x,y,z,surface id) | primitives plus INR freeform; analytic intersection for edges/corners | primitives yes (types not listed on pages fetched) | none stated |
| ComplexGen (Guo et al., SIGGRAPH 2022) | point cloud | chain complex of corners, curves, patches with types | type list unverified | global optimisation for topology; no tolerance |
| NVD-Net / Split-and-Fit (Liu et al., SIGGRAPH 2024) | point cloud, mesh or UDF -> 256^3 voxel grid | Voronoi cells, one analytic primitive per cell, analytic intersections | yes per cell (secondary source) | none stated |
| NH-Rep (Guo et al., SIGGRAPH Asia 2022) | points sampled from a B-rep | neural Boolean tree of implicits (B-rep -> implicit, opposite direction) | no; learned leaves | none; empirical |
| BR-DF (Zhang et al., arXiv 2511.14870, 2025-11-18) | generated SDF + per-face UDFs | faceted B-rep via extended marching cubes | no ("strictly speaking a faceted B-Rep") | watertight, "never fail", 100% success; no geometric tolerance |
| Mesh2Brep (Shen et al., IEEE TVCG 31(10), online 2025-01-06) | mesh | robust primitive fitting plus tangency/collinearity constraints | yes (classical fitting) | not in abstract |
| DualBrep (Liu, Jayaraman et al., arXiv 2606.31579, 2026-06-30, SIGGRAPH 2026) | mesh, point cloud, image | trimmed B-spline per face, STEP via OpenCASCADE | no, B-spline (16x16 grids) | median gap 5.8e-3 in unit box on 100 test shapes; BRepCheck tolerance 0.1 |
| CADFit (arXiv 2605.01171, 2026-05) | watertight mesh | CadQuery program: sketch lines/arcs/circles, extrude, revolve, fillet, chamfer | circles via algebraic fit with RMS constraint, spline fallback | kernel validity only; convergence bound under assumptions; DeepCAD IoU 0.960, CD 0.016 |

Sources:
- Schnabel 2007: https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8659.2007.01016.x ; CGAL implementation detects plane, sphere, cylinder, cone, torus: https://www.cgal.org/2015/03/15/cgal-wipRansac/
- Point2CAD: https://arxiv.org/pdf/2312.04962 , https://openaccess.thecvf.com/content/CVPR2024/html/Liu_Point2CAD_Reverse_Engineering_CAD_Models_from_3D_Point_Clouds_CVPR_2024_paper.html , code outputs "unclipped", "clipped", "topo" artifacts with no format named (fetched) https://github.com/prs-eth/point2cad . Authors: ML's main role is segmentation; method is sensitive to segmentation quality.
- ComplexGen: https://arxiv.org/abs/2205.14573 , https://haopan.github.io/complexgen.html
- NVD-Net: https://history.siggraph.org/learning/split-and-fit-learning-b-reps-via-structure-aware-voronoi-partitioning/ ; pipeline detail from secondary source https://www.emergentmind.com/topics/generative-modeling-of-b-reps
- NH-Rep: https://arxiv.org/abs/2209.10191 , https://github.com/guohaoxiang/NH-Rep ; tested on about 10,000 manifold B-reps including NURBS; manifold B-reps only.
- BR-DF: (fetched) https://arxiv.org/abs/2511.14870
- Mesh2Brep: https://doi.org/10.1109/tvcg.2025.3525844
- DualBrep: (fetched) https://arxiv.org/abs/2606.31579 , (fetched) https://github.com/AutodeskAILab/DualBrep (Apache-2.0; training code not released)
- CADFit: (fetched) https://arxiv.org/html/2605.01171v1 ; 240-900 s per shape; no lofts or sweeps.
- Other 2026 CAD-program reconstruction: CADENA arXiv 2608.00799 https://arxiv.org/html/2608.00799v1 ; CADReasoner arXiv 2603.29847 https://arxiv.org/html/2603.29847v1 ; Point2Primitive arXiv 2505.02043 https://arxiv.org/pdf/2505.02043 ; curated list https://github.com/Bigger-and-Stronger/awesome-brep-reconstruction
- Survey of unstructured data to CSG: arXiv 2305.01220 https://arxiv.org/pdf/2305.01220

None of the research methods found states a guaranteed physical tolerance between output B-rep and input.

---

## 5. Standards

STEP (ISO 10303):
- AP242 editions: Ed1 2014, Ed2 2020, Ed3 2022 (ISO 10303-242:2022), Ed4 ISO 10303-242:2025 published 2025-08 per ISO; a new revision (ISO/CD 10303-242) approved 2025-11-10, CD registered 2026-02-23. https://www.iso.org/standard/84300.html , https://www.iso.org/standard/93277.html . SIS lists approval 2025-08-25; a reseller says 2025-01-01. https://www.sis.se/en/produkter/manufacturing-engineering/industrial-automation-systems/industrial-process-measurement-and-control/ss-iso-10303-2422025/
- Ed4 abstract scope includes additive manufacturing part design; AFNeT lists Isogeometry and Polyhedral B-Rep extensions; neither mentions lattice or implicit geometry. (fetched) https://atlas.afnet.fr/en/iso-step-ap242-edition-4/
- AM module (from Ed2): build orientation and placement, minimum build volume, support geometry identification; supports were exchanged as tessellated geometry in prototype tests. https://www.ap242.org/additive-manufacturing.html (certificate error on fetch; via search summary)
- Lattices were listed as a planned extension for AP242 Ed3 and later in an Airbus GPDIS 2019 talk (unit cells replicated in a volume, or treated like heterogeneous materials). https://gpdisonline.com/wp-content/uploads/2019/09/Airbus-Delaunay-UseOfStepAp242Ed1ToSuppor-MBD-Open.pdf . No evidence it was adopted.
- Trivariate B-spline volumes exist in Part 42 but were not in ARM Part 1801; EU TERRIFIC/CAxMan proposed extensions (arXiv 1803.05756, 2018). https://arxiv.org/pdf/1803.05756
- Part 42 does define CSG: csg_solid, boolean_result, half_space_solid (semi-infinite solid bounded by a surface, optional box), and primitives (cone, cylinder, sphere, torus, block, wedge, ellipsoid, tetrahedron, pyramid). https://steptools.com/stds/smrl/data/resource_docs/geometric_and_topological_representation/sys/6_schema.htm . Application module ISO/TS 10303-1068 (CSG 3D). Which APs and translators support it was NOT found.
- STEP AP242 can carry tessellated and point-cloud data (Ed2). https://www.ap242.org/edition-2.html

AMF (ISO/ASTM 52915):
- Current: ISO/ASTM 52915:2020, AMF Version 1.2. Curved triangles via vertex normals or edge tangents with second-degree Hermite interpolation; subdivision depth exactly five (1024 flat triangles per curved triangle) in 2016/2020 editions; max 25% out of plane, no inflections. https://cdn.standards.iteh.ai/samples/74640/6ab4418c993c4a0d8e711fa896995478/ISO-ASTM-52915-2020.pdf , https://pdfcoffee.com/bs-iso-astm-52915-2016pdf-pdf-free.html
- `<composite>` material proportions can be formulas of x, y, z (graded, periodic substructure) or texture lookups; this is functional material, not implicit geometry. https://holdsworthwiki.robarts.ca/files/astm_iso_standards/iso-astm_52915_2013-e-additive_manufacturing_file_format_amf_version1.1.pdf
- No lattice or implicit-geometry element found in AMF.

3MF beam lattice: Beam Lattice Extension v1.2.0 (02/27/2025) carries lattices as beams (nodes and radii), separate from implicit. https://3mf.io/spec/

---

## Not found or uncertain

- Which shipping slicer, printer or build processor reads 3MF Volumetric levelset data. None documented.
- Whether nTop, EOSPRINT or Magics read or write 3MF implicit (vs the proprietary .implicit).
- How EOSPRINT and Materialise Build Processors sample an implicit, and the resulting accuracy; no independent accuracy test.
- Whether Fusion Prismatic writes analytic cylinder entities or NURBS; Autodesk docs silent.
- Whether QuickSurface, Design X, ScanTo3D, NX Fit Surface export analytic STEP entities (cylindrical_surface) vs B-spline; numeric deviations not published.
- nTop CAD Body fitting deviation.
- Current official list of Fusion CAM operations that accept meshes.
- 5-axis finishing directly on meshes in hyperMILL and PowerMill.
- Commercial CAM taking an implicit/SDF directly: none found.
- Point2CAD and ComplexGen exact primitive type lists (not on fetched pages).
- Whether AP242 Ed4 conformance includes Part 42 CSG/half-space solids; no lattice or implicit entity in AP242 found.
