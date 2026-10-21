# Research F: large CAD assemblies (display, mates, patterns, interference) and SDF use

Compiled 2026-10-21. Sources were read through web search snippets and a few direct page fetches.
"Vendor claim" marks a statement made by the product's maker or reseller in marketing or docs, not independently measured.
"Accessed" means the page carried no visible date, so only the access date is known.
Gaps are listed at the end.

## 1. Display of large assemblies

### SOLIDWORKS (Dassault Systemes)

- Three open modes: Resolved, Lightweight, Large Design Review (LDR).
  Source: SOLIDWORKS blog "How to Open (Large) Assemblies Even Faster", 2023-01-18, https://blogs.solidworks.com/products/solidworks/how-to-open-large-assemblies-even-faster
- Resolved loads features, solid bodies, surfaces, mates and display data.
  Lightweight loads everything except the feature history: body data and mates are loaded, so components can still be added and mates edited.
  LDR loads only display (graphics) data. The post says assemblies open "within seconds instead of minutes" (vendor claim, no part counts or timings given).
  Source: same post.
- A reseller says lightweight components hold "both the tessellation data and the parasolid data".
  Source: Hawk Ridge Systems support article (undated), https://support.hawkridgesys.com/hc/en-us/articles/203486053-SOLIDWORKS-File-Opening-Modes-Resolved-Lightweight-Large-Design-Review-Detailing-Quick-View
- LDR on open: reads the top-level file, loads saved positions and the last saved graphics data, then computes display triangles. It does not validate or update subcomponents or solve mates. An option ("Automatic check and update all components") checks for out-of-date graphics data.
  Reseller example: an assembly that took nearly 8 minutes resolved opened "in seconds" in LDR (reseller claim, no part count given).
  Sources: Javelin, 2017-03, https://www.javelin-tech.com/blog/2017/03/solidworks-large-design-review-mode/ ; Hawk Ridge, above.
- LDR introduced in SOLIDWORKS 2012. Later additions: "Edit Assembly" in LDR allows inserting components, adding and editing mates, and assembly-level patterns (not with flexible subassemblies); 2020 added linear, circular, driven and sketch-driven component patterns in top-level LDR assemblies and mating to component reference geometry; 2022 allowed opening a subassembly resolved from an LDR parent.
  Source: search summary of SOLIDWORKS reseller posts (TriMech, Javelin, Engineers Rule), https://trimech.com/using-large-design-review-solidworks-large-assembly-best-practices/ , https://www.engineersrule.com/the-only-two-tips-you-need-for-solidworks-assemblies/ (dates on individual claims not verified).
- Since SOLIDWORKS 2023 there is an option "Automatically optimize resolved mode, hide lightweight mode" that uses lightweight loading selectively.
  Large Assembly Settings (Tools > Options > System Options > Assemblies) switch to lightweight above a user-set component count threshold.
  Source: SOLIDWORKS blog 2023-01-18 (above), Hawk Ridge (above).
- Official SolidPractices guide: rebuild order is load referenced parts, update out-of-date ones, solve all mates and position components, update in-context and assembly features, then redraw. "Lightweight Mode loads only the body and the mate data". Large Assembly Mode loads lightweight by default.
  Source: SolidPractices "Best Practices for Managing Large Assemblies" PDF (date not shown in snippet), https://my.solidworks.com/docs/BP/TS_BP_SWX_LargeAssemblies.pdf

### Siemens NX, Teamcenter visualisation, JT

- NX fact sheet: NX Active Mockup uses JT so designers "can load thousands of components from different CAD systems in seconds", and "JT enables you to load assemblies with up to a million parts with exceptional interactive performance" (vendor claim, no hardware or test conditions given).
  Source: NX for Mechanical Design fact sheet (hosted by ATA Engineering, file uploaded 2024-10; fact sheet's own date unknown), https://www.ata-e.com/wp-content/uploads/2024/10/NXDesignFactSheets.pdf
- Teamcenter Visualization Mockup "leverages the lightweight JT visualization format" for loading, navigating and analysing very complex multi-CAD assemblies, and includes clearance and interference checks. Visualization Professional fact sheet: CAD users often use Teamcenter for large assembly viewing because their native CAD cannot handle the data (vendor claim).
  Sources: Teamcenter 2008-era fact sheets, https://media.plm.automation.siemens.com/teamcenter/Teamcenter2008_Web_English/PublishFolder/collateral/LV_Mockup_FS.pdf , https://www.plm.automation.siemens.com/fr_fr/Images/4591_tcm68-4790.pdf
- Third-party NX guide: reference sets (MATE, SIMPLE, DRAFT, LIGHTWEIGHT), partial loading, lightweight representations that "use lightweight JT + outer-boundary representation, reducing memory by up to 90%" (third-party claim), and Open by Proximity, which fully loads only components inside a box or sphere and leaves the rest lightweight.
  Source: whole-spec.com (undated), https://whole-spec.com/en/cad-tips/nx/nxopen-teamcenter-troubleshoot/
- Third-party JT viewer 3DViewStation claims to load "a complete system consisting of several million parts" (vendor claim).
  Source: https://www.3dviewstation.com/im-looking-for/nx-jt-viewer-visualization-dmu.html (accessed 2026-10-21)
- JT format (ISO 14306): each part carries a mandatory multi-LOD tessellation; exact B-rep is optional (JT native B-rep, Parasolid XT B-rep, or STEP B-rep from ISO 14306:2017). Segments include logical scene graph, shape LOD segments, optional B-rep, PMI and metadata. Late-loaded property atoms are part of the spec; the format is streamable and supports single or multi-file partitioning.
  ISO PAS 14306 accepted 2009-09-18; ISO 14306:2012 published December 2012; ISO 14306:2017 added STEP B-rep; ISO 14306-1/-2:2024 and -3:2025 listed by one source.
  Sources: Library of Congress format description, https://loc.gov/preservation/digital/formats/fdd/fdd000562.shtml ; CADInterop, https://www.cadinterop.com/en/formats/neutral-format/jt.html ; ISO 14306:2017 sample TOC, https://standards.iteh.ai/catalog/standards/iso/18fd19fd-36e6-452c-8494-aebda2332a93/iso-14306-2017

### Onshape (PTC)

- Rendering in the browser uses tessellation chosen automatically per Part Studio complexity. Official help: "Small components in assemblies are culled graphically to maintain framerate. For assemblies with fewer than 5,000 bodies, the threshold is 2 pixel screen size. For larger assemblies, the threshold is 10 pixels." "Culling is graphical only; the BOM does not change." When zooming or panning complex assemblies, "selecting an instance forces its tessellation to be refined first".
  Source: Onshape Help "Performance Considerations" (fetched 2026-10-21), https://cad.onshape.com/help/Content/Home/performance_considerations.htm
- Onshape employee on the forum: "there is no limit on the number of parts or memory (on our servers at least, your browser may need more memory)", limits depend on the graphics card. An Onshape architect described a slow case where finely tessellated screws instanced many times gave "100M facets / silhouette edges"; workaround was deleting thread faces.
  Source: Onshape forum discussion 258 (2015-era thread), https://forum.onshape.com/discussion/258/large-assembly-performance
- User reports: failures at about 5,019 parts and 7,000 mates; a patterned test of 13,842 part instances (1,538 per subassembly) ran at 25 to 30 fps on a GTX 970 (2020 post); another user reports lag above about 10,000 generated parts.
  Sources: https://forum.onshape.com/discussion/13233/large-assembly-modeling-what-is-the-maximum-reasonable-size-for-an-assembly ; https://forum.onshape.com/discussion/31026/performance-issues-when-handling-large-amount-of-parts
- Onshape Q&A: assemblies with many shared components have components "brought in just once".
  Source: upFront.eZine 1093, 2021-04-26, https://upfrontezine.substack.com/p/upfrontezine-1093-q-and-a-from-onshapes

### PTC Creo

- Simplified representations control what is loaded into memory (can omit feature history). Shrinkwrap collects external surfaces into a single non-associative part; methods Surface Subset, Faceted Solid, Merged Solid. PTC says shrinkwrap can "reduce disk and memory usage by more than 90 percent" (vendor claim).
  Sources: PTC blog "3 Ways to Keep Large Assemblies from Slowing You Down" (undated), https://www.ptc.com/en/blogs/cad/3-ways-to-keep-large-assemblies-from-slowing-you-down ; PDSVISION, https://pdsvision.com/blog/advanced-assembly-functionality-large-assembly-management-in-creo/
- Not found: specifics of Creo's "graphics representation" (what mesh/LOD is stored) or Creo part-count figures.

### CATIA (Dassault Systemes)

- V5 cache: the first time a component is inserted it is tessellated into a CGR file stored in a local cache; later loads read the CGR rather than the CATPart. A read-only "released cache" on the network can be shared. Visualization mode shows only the tessellated skin; geometry is not available.
  Source: CATIA V5 Cache Management docs (V5 era, about R14, circa 2004), http://catiadoc.free.fr/online/bascupst_C2/bascupst0100.htm ; Design Tech Academy, 2018-12, https://www.designtechacademy.com/2018/12/how-to-use-cgr-mode-for-large-assembly.html
- Default tessellation sag ("3D fixed accuracy") is 0.2 mm; sag can be changed per object with Tools > Modify SAG for clash work.
  Source: CATIA V5 DMU Space Analysis guide R14, https://catiadesign.org/_doc/docv5r14/DMU%20Space%20Analysis.pdf
- 3DEXPERIENCE "lightweight opening mode": index-based information, progressive viewpoint-dependent loading (cloud only with VIK licence), only visualization data kept in memory, optimized and simplified 3D representations. Claimed to accelerate opening of large assemblies "by up to 100 times" (vendor community claim, no part count).
  Source: 3DSwym CATIA community post (about 2026-01), https://3dswym.3dexperience.3ds.com/post/catia-user-community/do-you-know-the-new-lightweight-opening-mode_HNPXiThGS5mSWW0fNWjfTg
- EDS Technologies (reseller): associative assembly patterns slow "significantly" above about 1,000 parts.
  Source: https://edstechnologies.com/blog/boosting-performance-in-large-assembly-patterns-drafting-views-with-3dexperience-catia/ (undated)

### HOOPS Visualize (Tech Soft 3D)

- Fixed-framerate mode draws in priority order by screen size and distance to viewer, interrupts to hold the target rate, and continues drawing undrawn items when the camera stops. Small objects are sacrificed first.
  Static model mode; "AttributeSpatial" adds spatial segments and is "more appropriate for very large models".
  Culling: backplane, view frustum, extent (screen-size) culling; occlusion culling on qualifying OpenGL hardware. Legacy 3DF default framerate 0.1 s (10 fps) with 150-pixel threshold; "Target Framerate" mode supports LODs.
  Sources: HOOPS Visualize HPS 2024 docs, https://docs.techsoft3d.com/hps/2024/prog_guide/0703_performance_considerations.html ; 3DF docs, https://docs.techsoft3d.com/hoops/visualize-3df/prog_guide/mvo/05_maintaining_constant_framerate.html
- No part-count benchmarks found in HOOPS docs.

### Parasolid facet bodies (Siemens)

- Convergent Modeling (from Parasolid v26 per a third party) combines classic B-rep and facet (triangle) geometry in one body. v31.1: when importing facet data "the creation of the mesh may be delayed until it is needed". v32.0: mixed models with a subset of classic surfaces (analytic planes, cylinders, spheres), stated use case "restore precise, classic mating surfaces" to topology-optimized facet models. v33.0 (2021-01): intersection of meshes with all classic surface types. v39.0 (about 2026-08): B-rep, facet and implicit geometry in the same model; mesh-to-B-rep fitting from mesh vertices.
  Sources: Parasolid v32.0 highlights, https://blogs.sw.siemens.com/plm-components/parasolid-v32-0-release-highlights/ ; Architosh 2021-01, https://architosh.com/2021/01/siemens-expends-convergent-modeling-new-parasolid-v33-0/ ; Parasolid v39.0, https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/ ; CADInterop, https://www.cadinterop.com/en/formats/neutral-format/parasolid.html
- Not found: Parasolid's own "lightweight body" display representation or facet-cache mechanism for viewing.

### NVIDIA Omniverse and OpenUSD instancing

- OpenUSD scenegraph instancing: "instance" prims share a "prototype"; descendants of instanceable prims are read-only (local opinions ignored). Point instancing is a vectorized form for large counts with per-instance position, orientation, scale.
  Source: OpenUSD docs, https://openusd.org/dev/api/_usd__page__scenegraph_instancing.html ; Omniverse modularity guide, https://docs.omniverse.nvidia.com/usd/latest/learn-openusd/independent/modularity-guide/instancing.html
- Omniverse CAD converters (HOOPS, JT, DGN) set instanceable = true on product occurrences that reference a prototype; the JT converter only does this for JtkINSTANCE nodes. Stated benefits: geometry and materials stored once, faster traversal and rendering, smaller files.
  Source: https://docs.omniverse.nvidia.com/kit/docs/omni.kit.converter.hoops/latest/Instancing.html ; https://docs.omniverse.nvidia.com/kit/docs/omni.kit.converter.jt/509.0.0/Instancing.html
- No published instance-count or frame-time benchmarks found.

### Unreal Engine Nanite and Datasmith

- Datasmith imports STEP, CATIA V5, NX, SOLIDWORKS, Inventor; tessellates B-rep at import with per-mesh tessellation overrides; keeps hierarchy; deliberately avoids merging the scene into one mesh.
  Source: Epic docs, https://dev.epicgames.com/documentation/unreal-engine/importing-cad-files-into-unreal-engine-using-datasmith
- Nanite Assemblies (UE 5.7): one assembly holds up to 65k part instances; part clusters are encoded into the final mesh's hierarchy without duplicating geometry, and per-instance transforms are resolved during cluster culling. Epic demo: largest tree's disk size dropped from 3.5 GB to about 29 MB. All documented examples are foliage, authored through USD.
  Source: Epic docs, https://dev.epicgames.com/documentation/unreal-engine/nanite-assemblies ; https://dev.epicgames.com/documentation/en-us/unreal-engine/nanite-foliage
- Twinmotion 2025.2 (2025-09/10) added Nanite for imported CAD/Datasmith data; Epic says it handles "hundreds of millions or even billions of polygons" in real time (vendor claim); static meshes only.
  Sources: Digital Production 2025-09-30, https://digitalproduction.com/2025/09/30/nanite-marches-in-twinmotion-2025-2-update/ ; Epic Twinmotion docs, https://dev.epicgames.com/documentation/twinmotion/optimizing-geometry-with-nanite-in-twinmotion
- Academic data point (not CAD): 1-billion-triangle terrain viewed in real time with Nanite; "Import and conversion to Nanite crashes when meshes are larger than 40 million triangles".
  Source: Uruk-Warka digital twin, International Journal of Digital Earth, 2024, https://www.tandfonline.com/doi/full/10.1080/17538947.2024.2324964
- Third party (CADInterop) advises pre-converting B-rep to optimized mesh with LOD for assemblies over 10,000 parts.
  Source: https://www.cadinterop.com/en/formats/mesh/unreal.html

### Other rendering data points

- Transparency rendering of an assembly with over 16 million triangles in about 10,000 parts.
  Source: "Interactive transparency rendering for large CAD models", IEEE TVCG 2005, https://doi.org/10.1109/tvcg.2005.82
- Large BIM models: hardware instancing on top of occlusion culling gave additional 1.25x to 1.7x speedups.
  Source: Chalmers thesis/paper, https://publications.lib.chalmers.se/records/fulltext/173349/local_173349.pdf

## 2. SDF and implicit rendering at scale

### Dreams (Media Molecule, PS4)

- SIGGRAPH 2015 talk (Alex Evans): engine runs almost entirely on PS4 compute "(no triangles!)"; scenes are CSG edit lists evaluated on the fly to high-resolution SDFs, from which dense multi-resolution point clouds are generated. SDF stored in 8x8x8 fp16 bricks, incrementally updated, with pre-filtered mips for LOD. Download-size constraint: SDFs are not stored, so a fast evaluator is needed. An abandoned "refinement renderer" was 4 to 10 times too slow for PS4 and gigavoxel brick memory was a problem.
  Sources: Advances in Real-Time Rendering 2015, https://advances.realtimerendering.com/s2015/ ; slides PDF, http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf ; Umbra Ignite video post, https://www.mediamolecule.com/blog/article/alex_at_umbra_ignite_2015_learning_from_failure_video
- After 2015 (Evans, quoted on forums, not verified at source): the atomic splat renderer left microscopic holes, so the "solidness" of the brick engine was brought back; sculpts are edit lists evaluated by async compute into a mega-texture of 8x8x8 bricks. A forum user (not Evans) describes the later renderer as rasterising bounding cubes and ray tracing the SDF in a shader.
  Sources: Beyond3D thread, https://forum.beyond3d.com/threads/signed-distance-field-rendering-pros-and-cons-as-used-in-ps4-title-dreams-spawn.57006/ ; Hacker News, https://news.ycombinator.com/item?id=16839166
- Released 2020-02-14. No published object counts or frame times found.

### Claybook (Second Order, GDC 2018)

- World and characters are SDFs; target 60 fps with ray-traced visuals and GPU simulation on consoles. World SDF 1024x1024x512, 8-bit signed, 586 MB with 5 mips, distance range [-4, +4] voxels (1/32 voxel precision). Brushes are pre-baked volumes 32^3 to 128^3; "Runtime performance doesn't depend on brush count".
  Source: Aaltonen, GDC 2018 slides, https://media.gdcvault.com/gdc2018/presentations/Aaltonen_Sebastian_GPU_Based_Clay.pdf
- Per-pass millisecond timings not captured in this search.

### Unreal Engine mesh distance fields (SDF used at scale for lighting, not primary display)

- Each static mesh gets an SDF volume texture, max 128^3 (8 MB). A camera-following Global Distance Field merges per-object SDFs into clipmaps; only changed regions update. Lumen software ray tracing traces per-mesh SDFs for the first 2 m of a ray, then the global field. CVar help warns that scenes with high instance density (overlapping meshes) have high tracing costs.
  Sources: Epic docs, https://dev.epicgames.com/documentation/en-us/unreal-engine/mesh-distance-fields-in-unreal-engine ; https://dev.epicgames.com/documentation/en-us/unreal-engine/lumen-technical-details-in-unreal-engine

### NanoVDB and sparse volumes

- NanoVDB: pointer-free linearized read-only GPU version of OpenVDB, portable to CUDA, OpenCL, GL, WebGL, DX12, OptiX, HLSL, GLSL.
  Source: Museth, SIGGRAPH 2021 Talks, DOI 10.1145/3450623.3464653, https://history.siggraph.org/learning/nanovdb-a-gpu-friendly-and-portable-vdb-data-structure-for-real-time-rendering-and-simulation-by-museth/
- Path tracer on NanoVDB: 1920x1080, RTX 3090, 4 spp with denoising in 61 ms (volumes, not CAD surfaces).
  Source: NVIDIA, "Interactive Path Tracing and Reconstruction of Sparse Volumes", 2021-03, https://research.nvidia.com/publication/2021-03_interactive-path-tracing-and-reconstruction-sparse-volumes
- fVDB uses hierarchical DDA for GPU ray marching on NanoVDB-style trees.
  Source: arXiv 2407.01781 (2024), https://arxiv.org/html/2407.01781v1

### Implicit engineering tools

- nTop (nTopology) 3.0: GPU acceleration for real-time visualization of implicit bodies; claims "10x to 100x performance increase" (vendor claim). For export, nTop's "High" complexity tier is over 10k beams/faces or over 1k TPMS cells, and it recommends a "Simplified Body Representation" for those.
  Sources: nTop 3.0 announcement, https://resources.ntop.com/resources/product-updates/ntopology-3-0-release-announcement/ ; nTop export guide, https://support.ntop.com/hc/en-us/articles/360048784874-How-to-export-geometry-Overview-of-export-workflows
- Research: a lattice with 14,063,027 struts handled in streaming mode; "The solid implicit model is rendered by direct ray-tracing".
  Source: "Memory-Efficient Modeling and Slicing of Large-Scale Adaptive Lattice Structures", arXiv 2101.05031 (2021), https://arxiv.org/pdf/2101.05031
- Blog on SDFs for CAD: ray marching is cheap because each step can advance by f(p); Fornjot chose B-rep because common operations (booleans) do not preserve the exact distance property, only a lower bound.
  Source: James Stanley, "Thoughts on Signed Distance Functions for CAD" (undated), https://incoherency.co.uk/blog/stories/sdf-thoughts.html

### CAD or engineering viewer rendering whole assemblies from SDFs

- Not found. No commercial or academic CAD assembly viewer that renders thousands of B-rep parts from SDFs turned up in these searches.

## 3. Mates and assembly constraints

### Geometric references

- SOLIDWORKS treats a subassembly as a rigid body (except flexible subassemblies); "it is possible to locate even the largest subassembly with only three mates". Flexible subassemblies promote their mates into the top-level solve.
  Sources: SolidPractices PDF (above); Hawk Ridge (2026), https://hawkridgesys.com/blog/tips-to-eliminate-solidworks-large-assembly-performance-issues
- Onshape: mates do not reference geometry directly. Each mate aligns two mate connectors, each a full local coordinate system; mate type sets which of the six DOF remain free. Implicit connectors are inferred at face centroids, edge midpoints, vertices, centres of circular openings, centres of cuts, and virtual sharps of conic faces; "changing the geometry will change the location of the Mate".
  Sources: Onshape Help, https://cad.onshape.com/help/Content/mateconnector.htm ; https://cad.onshape.com/help/Content/Assembly/mates.htm ; Onshape blog, https://www.onshape.com/en/blog/kinematics-assemblies-mates-simulation
- Onshape: rigid subassembly mates are not solved in higher-level assemblies.
  Source: Onshape Help Performance Considerations (above).

### Solver: Siemens D-Cubed 3D DCM

- 3D DCM solves assembly constraints, kinematics, and 3D sketching; Siemens calls it "the most widely adopted 3D geometric constraint solver for CAD/CAM/CAE applications" (vendor claim).
  Geometry supported: "points, lines, planes, circles, ellipses, cylinders, tori and parametric curves/surfaces", and constraining "free forms curves and surfaces". Constraints: distance, angle, radius, parallel, perpendicular, tangent, concentric, symmetric, normal, midpoint, pattern; dimensions can be bounded to ranges.
  Source: Siemens 3D DCM page (fetched 2026-10-21), https://www.siemens.com/en-us/products/plm-components/d-cubed/3d-dcm/
- A part is a "set" (rigid body); an unconstrained part has 6 DOF. The application supplies geometry and evaluators to DCM through the "Frustum" callback interface, and DCM moves geometry back.
  Sources: search summary of Siemens docs; USPTO patent 6963824, https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6963824
- 2013 Procedia CIRP paper describes 3D DCM as "a fast non sequential solver" (simultaneous rather than one-by-one).
  Source: Procedia CIRP 10 (2013) 169, https://www.sciencedirect.com/science/article/pii/S2212827113005775
- Version 59: automatic partitioning of unrelated entities, solved independently in one evaluate call, "significant performance improvements for larger models that can be partitioned" (vendor claim, no numbers). Version 58: accuracy far from the origin.
  Sources: https://blogs.sw.siemens.com/plm-components/d-cubed-3d-components-version-59-release/ ; https://blogs.sw.siemens.com/plm-components/d-cubed-3d-dcm-version-57-0/
- Licensees on Siemens' pages include Autodesk, DataSolid, Fujitsu, IronCAD, ISD, Onshape, PTC, Shapr3D, ZWCAD (and Top Systems in a blog list). Siemens also uses it in NX and Solid Edge (engineering.com). Wikipedia names SOLIDWORKS; Siemens' own lists do not include Dassault. Licensee lists do not say whether a product uses 3D DCM for mates or only 2D DCM for sketches.
  Sources: https://www.siemens.com/en-us/products/plm-components/d-cubed/ ; https://mobile.engineering.com/amp/17309.html ; https://en.wikipedia.org/wiki/Geometric_constraint_solving

### How many mates

- No published mate-count limit or scaling figure found for 3D DCM, SOLIDWORKS, or Onshape.
- Anecdotes: Onshape user failures at about 5,000 parts and 7,000 mates (forum, above); CAD Booster blog claims advanced mates are about 10 times slower than standard mates (blog claim), https://cadbooster.com/best-practices-assemblies/ ; Javelin advises few fixed references because more fixed items slow solving, https://www.javelin-tech.com/blog/2023/07/3-tips-for-mating-best-practices-in-solidworks-2023/

### Mates on non-analytic faces (mesh, freeform)

- SOLIDWORKS: reference axes, planes and coordinate systems can be built on graphics mesh bodies and mesh BREP bodies by selecting facets (planar reference), facet fins (linear edge reference), and facet vertices (point reference). For cylinders, "Surface From Mesh" fits a cylinder from selected facets, then an axis on it is used for a concentric mate; a fitted plane is used for coincident or width mates. For scanned or sculpted meshes, Paint Select picks many facets for the fit.
  Sources: SOLIDWORKS 2020 What's New, https://help.solidworks.com/2020/english/WhatsNew/t_add_reference_axes_planes.htm ; SOLIDWORKS blog "Lesson 4: Mesh Body Modeling Tools" (2025-05), https://blogs.solidworks.com/products/solidworks/lesson-4-solidworks-mesh-body-modeling-tools/
- SOLIDWORKS "Segment Imported Mesh Body" groups facets into planar or cylindrical faces; works best on CAD exported to mesh and back, scans often too irregular. Importing a mesh as standard BREP makes each facet a face, which "can decrease performance for large models".
  Source: SOLIDWORKS 2025 Help, https://help.solidworks.com/2025/english/solidworks/sldworks/c_segmented_imported_mesh_body_tool.htm ; Hawk Ridge, https://support.hawkridgesys.com/hc/en-us/articles/203485653-Working-with-Mesh-Files-in-SOLIDWORKS
- SOLIDWORKS 2018 (per Engineers Rule) allowed inserting and "precisely" mating components containing graphics bodies; supported mate types not stated.
  Source: https://www.engineersrule.com/what-a-fine-mesh-were-in/
- Parasolid v32 Convergent Modeling: restore classic planes/cylinders/spheres onto facet models as mating surfaces (vendor use case, above).
- 3D DCM accepts parametric and freeform curves and surfaces (vendor page, above); no detail on how tangency or distance to freeform surfaces is solved or how robust it is.

## 4. Patterns and instancing

- SOLIDWORKS Geometry Pattern option copies faces and edges of the seed instead of re-solving the feature; ignores end conditions. Without it, rebuild time scales roughly with instance count (example: 0.1 s seed x 8 instances ~ 0.8 s). Blog measurements: 0.50 s to 0.23 s and 0.81 s to 0.26 s (small test models).
  Sources: SOLIDWORKS blog 2013-09, https://blogs.solidworks.com/solidworksblog/2013/09/solidworks-geometry-pattern-a-simple-option-to-make-patterns-regenerate-faster.html ; SOLIDWORKS 2016 Help, https://help.solidworks.com/2016/English/solidworks/sldworks/c_Geometry_Pattern_Concept.htm ; Javelin 2022-03, https://www.javelin-tech.com/blog/2022/03/solidworks-geometry-pattern-option/
- User measurement: linear pattern rebuild time grows faster than linearly; trendline projects 11,289 s at 1,320 instances. Body pattern plus Combine was 6x faster than a feature pattern at 255 instances but slower beyond about 1,080 instances.
  Source: SOLIDWORKS forum thread 201949 (date not captured), https://forum.solidworks.com/thread/201949
- SOLIDWORKS blog case: switching feature patterns to body pattern plus Combine cut rebuild from 108 s to 12 s; without Combine 0.75 s.
  Source: https://blogs.solidworks.com/products/solidworks/increasing-performance-in-solidworks-patterns/
- SOLIDWORKS: mates are solved before assembly-level patterns, so mating to a patterned component forces recalculation (reseller claim).
  Source: search summary of CAD Booster / Peter CAD, https://petercad.com/solidworks-best-practices/best-practices-assemblies/
- Onshape: copying or patterning inside a Part Studio creates separate parts (BOM shows quantity 1 each); inserting or patterning in an Assembly creates instances of one part (BOM rolls up), and Onshape says instancing improves regeneration and document performance. Assembly patterns are not rigid; patterned components can be moved individually.
  Sources: Onshape tech tip (about 2026-10), https://www.onshape.com/en/resource-center/tech-tips/multi-part-studio-vs-assembly ; https://cad.onshape.com/help/Content/Assembly/assembly_linear_pattern.htm
- Onshape forum: replicate or pattern fasteners rather than inserting them one by one; adding mates slows with hundreds of replicated fasteners.
  Source: https://forum.onshape.com/discussion/19639/advice-on-assembly-structure-or-methodology-for-populating-instances-in-large-quantities
- Instancing on the display side: OpenUSD prototypes, Omniverse converter instancing, Nanite Assemblies up to 65k parts per assembly (section 1).
- Not found: NX pattern component storage (occurrence vs copy); CATIA pattern storage beyond the reseller slowdown note; cost figures for component patterns of thousands of instances.

## 5. Interference and clearance detection

### Vendor tools

- NX Clearance Analysis: "Calculate Using" exact or lightweight geometry; "Lightweight setting is the fastest and least accurate, and the Exact setting is the slowest and most accurate". Results can differ in both directions, "most likely to occur with curved surfaces that use a coarse faceting tolerance". "The penetration depth is always calculated using lightweight geometry." Types: hard, soft (within clearance zone), touching, containment.
  Source: NX Help (hosted by Univ. of Rochester ME204; NX version not shown), http://www2.me.rochester.edu/courses/ME204/nx_help/en_US/tdocExt/content/6/xid677726.xml and .../anly_asmclr_clr_prop_basic_op.xml (host did not resolve on direct fetch 2026-10-21; content from search snippets).
- NX user reports: exact mode reports false "hard" interferences that are really touching; STEP imports with cylinders approximated by polyhedra gave interferences nearly everywhere.
  Sources: https://www.eng-tips.com/threads/clearance-analysis-false-hard-interferences.402568/ ; https://community.sw.siemens.com/s/question/0D54O000061wxCnSAI/clearance-analysis-nx10
- SOLIDWORKS 2025: interference detection in LDR. Official help: "In Large Design Review mode, the volume of interference is not available and calculations for interference detection are approximate. For accurate results, resolve the components and recalculate". "False positives ... can occur when curved geometries are coincident or close enough to cause facets to collide." Treat coincidence as interference and size filters are not available in LDR. (That LDR uses tessellated data is implied by the facet wording, not stated outright.)
  Source: SOLIDWORKS 2025 Help, https://help.solidworks.com/2025/English/SWConnected/swdotworks/HIDD_DVE_INTERFERENCE.htm
- SOLIDWORKS 2024 SP3 added interference between surface bodies.
  Source: https://help.solidworks.com/2024/english/WhatsNew/c_wn2024_assemblies_interference_detection.htm
- CATIA V5 DMU Space Analysis: clash results depend on tessellation sag (default 0.2 mm); finer sag via Modify SAG in visualization mode with the cache on. On CGR-only data, clash, sectioning and compare are supported; distance/band analysis only partly.
  Sources: DMU Space Analysis R14 guide (above); "Working with CGRs in DMU", http://catiadoc.free.fr/online/cfyugdmu_C2/cfyugdmu0011.htm
- Teamcenter Vis Mockup: clearance analysis "in a matter of minutes" (blog, no assembly size).
  Source: PLM TechTalk 2019-11-14, https://plmtechtalk.com/2019/11/14/teamcenter-visualization-mockup-clearance-analysis/
- 3DViewStation: clash and clearance results "in seconds, even on assemblies of several thousand parts" (vendor claim).
  Source: https://www.cadinterop.com/en/our-products/3dviewstation/421-the-simplest-solution-to-detect-collision-in-a-cad-model.html
- Navisworks Clash Detective (BIM): Hard, Hard (Conservative), Clearance, Duplicate tests on exported triangle meshes with tolerances; practitioner guides say large models take 30 to 120 s per test (unsourced practitioner figure). Users report Hard misses clashes near mesh limits that Conservative catches.
  Sources: Autodesk Community thread (2021), https://forums.autodesk.com/t5/navisworks-forum/navisworks-clash-distances-and-hard-vs-conservative-issues/td-p/10306328 ; guide, https://www.atelier7arch.com/2026/04/navisworks-clash-detective-complete.html

### Methods in literature

- OBBTree (Gottschalk, Lin, Manocha, SIGGRAPH 1996): oriented bounding box hierarchy over triangles; cost model T = N x B + P x C (BV tests times BV cost plus primitive pairs times primitive cost). Names interference detection as a bottleneck for CAD virtual prototyping.
  Source: UNC TR96-013, https://www.cs.unc.edu/techreports/96-013.pdf
- GPU image-space interference (CInDeR, 2003) via frame-buffer virtual ray casting; related GPU work reports 95 ms for all collisions in a half-million-triangle scene and about 4 fps at 800K triangles (game/animation scenes, not CAD).
  Sources: https://www.researchgate.net/publication/2875907_CInDeR_Collision_and_Interference_Detection_in_Real-time_using_Graphics_Hardware ; https://www.researchgate.net/publication/228646794_INTERACTIVE_COLLISION_DETECTION_FOR_DEFORMABLE_AND_GPU_OBJECTS
- IBM patent US 5,497,453: graphics-hardware interference detection, motivated by exact solid interference being too slow for interactive sessions on industrial models.
  Source: https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/5497453

### Distance fields and voxels for interference

- Boeing Voxmap PointShell (VPS), McNeely, Puterbaugh, Troy, SIGGRAPH 1999: static scene as voxel map, moving object as surface point shell with normals; sustains 1000 Hz haptic rate; used in a Boeing CAD application for assembly and maintenance verification of complex aircraft systems. Approximate; assumes no deep penetration.
  2006 follow-up: order-of-magnitude spatial accuracy improvement by adding distance fields, temporal coherence and culling, trading away constant-time performance.
  Sources: https://www.researchgate.net/publication/234820066_Advances_in_voxel-based_6-DOF_haptic_rendering ; https://digital.lib.washington.edu/researchworks/bitstreams/54910b02-279c-481b-878b-20e30327a227/download ; Boeing figure, https://www.researchgate.net/figure/Boeing-computer-aided-design-CAD-application-for-assembly-and-maintenance-verification_fig4_226035748
- DLR improvements to VPS: fast voxmap/pointshell generation, 1 kHz collision response "with arbitrarily complex scenarios" (authors' claim).
  Source: https://www.semanticscholar.org/paper/Improvements-of-the-Voxmap-PointShell-Algorithm-of-Sagardia-Hulin/7ee9d51b409debd9e3e114d62a05bdd7cc434995
- 6-DOF haptic system to verify assembly procedures on CAD models built on VPS.
  Source: arXiv 1909.12714 (2019), https://arxiv.org/pdf/1909.12714
- "Assemble Them All" (SIGGRAPH Asia 2022): physics-based assembly-by-disassembly planning; each collision shape carries an SDF for contact with complex geometry such as fine screw threads; dataset of "thousands of physically valid industrial assemblies"; solves 80-part assemblies "within several minutes".
  Source: arXiv 2211.03977, https://arxiv.org/abs/2211.03977 (SDF detail from search snippet of the PDF)
- Convex body collision via SDF (Computer-Aided Design vol. 170, 2024): minimum distance and penetration depth via optimization (ellipsoid algorithm for convex); performance similar to GJK and MPR.
  Source: https://www.sciencedirect.com/science/article/pii/S0010448524000125
- Real-time collision detection between general SDFs (CAGD 2024): splits SDFs into parts with bounding volumes, uses interval arithmetic and gradients to find penetration depth, contact points and normals; accepts continuous and discrete SDFs.
  Source: http://www.cad.zju.edu.cn/home/jin/papers/Real_Time_CD_between_SDFs.pdf
- Distance fields for rigid/deformable collision (GraphiCon 2003): authors say "more exact than previous approaches using voxelisation".
  Source: https://www.graphicon.ru/html/2003/Proceedings/Technical/paper495.pdf
- Macklin et al. 2020, "Local optimization for robust signed distance field collision", Proc. ACM CGIT 3 (cited by the CAGD paper; not read directly).

## Not found or uncertain

- Mate-count limits or scaling numbers for any solver, including 3D DCM.
- Whether SOLIDWORKS uses 3D DCM for mates (only Wikipedia says D-Cubed is in SOLIDWORKS).
- How 3D DCM handles constraints on freeform or mesh faces in practice (accuracy, convergence).
- Any CAD or engineering viewer that renders full assemblies from SDFs; no part counts or frame times for SDF assembly rendering.
- Dreams object counts or frame times; Claybook per-pass ms timings.
- Independent benchmarks for the "up to a million parts" JT claim, Omniverse instance counts, Nanite on CAD assemblies with part counts.
- Creo graphics representation contents; NX and CATIA pattern storage (occurrence vs copy); component-pattern cost at thousands of instances.
- Peer-reviewed comparisons of exact B-rep versus tessellated interference checks with speed and accuracy on large assemblies.
- Whether any commercial CAD interference tool uses distance fields (only research and Boeing VPS found).
