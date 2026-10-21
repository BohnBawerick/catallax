# Research A: commercial CAD systems and kernels that mix B-rep with implicit / SDF / voxel geometry

Compiled 2026-10-21. Every fact carries a source URL and date where known.
"Vendor claim" marks marketing or self-reported performance statements.
Where a source did not answer a question, the gap is marked GAP.


## 1. Siemens Parasolid (kernel)

### Timeline of mixed-representation features

- v28 era (2016): Convergent Modeling introduced - facet (mesh) B-rep surfaces in the same body as classic B-rep. Press release "Siemens' Parasolid introduces Convergent Modeling on mixed models": https://news.siemens.com/en-us/parasolid-convergent-modeling-mixed-models/ (date not captured; GAP).
- v32.1: lattices "supported in Parasolid v32.1 as a new geometry type to address the need to accurately represent lattice materials in parts intended for additive manufacturing" - these were rod-and-ball (strut) lattices, not implicit. https://blogs.sw.siemens.com/plm-components/parasolid-v32-1-release-highlights/
- v33.1 (press release 2021-07-07): "new lattice modeling functions help users to create structures that benefit from a high strength-to-mass ratio"; Convergent Modeling expanded with "blending, tapering, offsetting, hollowing and thickening mixed models", letting users "integrate their precisely engineered (B-rep) designs with organic (facet-based) shapes". The press release never mentions implicit geometry. https://news.siemens.com/en-us/siemens-geometric-modeling-parasolid-release/ . v33.1 highlights blog also lists "Native Lattice Support": lattices attached to parts as construction geometry, clipping to a surface, meshes from lattices, lattices transmitted/received in files. https://blogs.sw.siemens.com/plm-components/parasolid-v33-1-release-highlights/
- v34.1: more lattice features (rods with different radii, combining lattices). https://blogs.sw.siemens.com/plm-components/parasolid-v34-1-release-highlights/
- v36.x, v37.x: GAP - no release-highlight page found that describes implicit features for these versions. A third-party site describes v34-v37 (2023 onward) as "Convergent Modeling as core technology" with no implicit mention. https://www.cadinterop.com/en/formats/neutral-format/parasolid.html
- v38.0 (blog dated 2025-08-04): first official wording found that uses "implicit": "Cellular bodies containing implicit pattern lattices" with slicing and enlarging, mass properties and local density of the sub-topology; "Support for the visualization of standalone implicit pattern lattices". Embedded rod-and-ball lattices get containment and faceting control. The post does not define "implicit" beyond "implicit pattern lattices". https://blogs.sw.siemens.com/plm-components/whats-new-in-parasolid-v-38-0/
- v39.0 (blog dated 2026-08-25; Wikipedia gives release 2026-07-31): Siemens says Parasolid can "support boundary representations (B-reps), facets, and implicit geometry within the same 3D model". Implicit work named in the post is lattice-only: enhanced clash for "cellular bodies containing embedded lattices", option to "exclude lattice and frame contributions from mass property calculations", slicing that generates "geometrically matched polylines in overlap regions". No implicit Booleans, no general implicit functions, no implicit-to-B-rep conversion described. The one conversion feature is mesh-to-B-rep fitting: "Surfaces can now be fitted using mesh vertices alone" for scan-to-CAD. https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/ ; https://en.wikipedia.org/wiki/Parasolid

### What "implicit" covers in Parasolid (answer to the brief's question)

Best source: CIMdata commentary on the Siemens PLM Components Innovation Conference, Cambridge UK, 2026-04-28/29. https://www.cimdata.com/en/component/docman/doc_download/4175-siemens-plm-components-2026-commentary
- "Siemens is expanding Parasolid's Convergent Modeling framework to natively support and mix implicit geometry seamlessly without forcing full format conversions."
- Parasolid now enables "operations directly on cellular bodies containing implicit pattern lattices, ensuring that downstream inquiries like slicing or mass property calculations can occur without sacrificing accuracy or design intent."
- Key limit: "Currently lattices and transient objects leverage the new technology, and implicit objects will become persistent in future releases." So as of April 2026, implicit geometry in Parasolid is lattice-scoped plus transient (non-saved) implicit objects; general persistent implicit bodies are stated as future work.
- George Allen (nTop Technical Fellow) at the same event: converting implicit bodies into B-reps "required by many downstream manufacturing and drafting applications remains a significant challenge - often considered to be a necessary evil that can result in distorted surfaces and reliability issues."
- PTC Onshape (Paul Chastell) demonstrated "convergent parts where boundary representation features and AI generated implicit geometries coexist", using FeatureScript and an MCP server, on Parasolid.
- Conference customer/partner presenters listed: Altair (Siemens), AMC Bridge, Cenit, Cimatron (Sandvik), Encube, Hexagon, InfinitForm, Kreodx, Materialise, nTop, Onshape (PTC), Quicksurface, Synera, Tech Soft 3D, Toolpath.com. (Presence on this list shows use of Siemens PLM Components, not necessarily Parasolid for implicit work.)

Summary for Parasolid:
1. Master: B-rep remains the primary persistent model; facets are a peer representation inside convergent bodies; implicit content is (as of v38/v39) implicit pattern lattices inside cellular bodies plus transient implicit objects.
2. Sync: lattices are carried inside the body (not converted); slicing, mass properties, clash work directly on them.
3. Users: AM lattice parts via host apps (NX, Solid Edge, Onshape, etc.).
4. Limits: general persistent implicit bodies stated as future; no published implicit-to-exact-B-rep conversion.
5. Formats: Parasolid X_T/X_B (lattices "transmitted and received" since v33.1); slicing output as polylines.

## 2. Siemens NX (now branded Designcenter NX) - Implicit Modeling

- Introduced in NX 1953 (release post 2020-12-11): "new optional module for the creation of equation driven structures", for structures "not feasible using traditional CAD methods"; filed under the Convergent Modeling heading. https://blogs.sw.siemens.com/designcenter/new-nx-december-release-details/
- Siemens Additive blog (2021-01-25, Aaron Frankel): "general equation builder with pre-defined examples" and out-of-the-box TPMS including gyroids; vendor claim that its "shelling, Boolean, and blend operations" are "extremely robust". Applications named: aerospace heat exchangers, medical implant fills. Post says nothing about outputs or B-rep. https://blogs.sw.siemens.com/additive/aa-e10/
- GfxSpeak (2021-02-11, quoting Siemens): dedicated task environment "with a library of pre-defined examples. The new tool can output convergent geometry for downstream reuse." https://gfxspeak.com/2021/02/11/processors-incorporates-paradigms/
- Older NX lattice tool: "the lattice command creates a faceted body, and the Unite command combines it with the part into a single convergent body". https://community.sw.siemens.com/s/article/introducing-lattice-in-nx (date GAP)
- NX 2021 (DEVELOP3D review): Voronoi lattices, lattices in multibody parts, lattice connected to part bulk; lattice-to-bulk smoothing not automatic at that time. https://develop3d.com/cad/siemens-nx-2021-q1-review-cad-modelling/
- NX 2506 (June 2025) release notes (Scribd copy): "GPU based preview for Lattice and Implicit is only supported for fully shaded style." Official "What's New in NX 2506" lists "Lattice Designer enhancements"; no separate implicit entry seen in the excerpt. https://www.scribd.com/document/890015554/NX-2506-Release-Notes ; https://docs.sw.siemens.com/en-US/doc/209349590/PL20241101461013487.whatsnew/xid1851820
- Convergent bodies and exchange (user forum, Siemens community): older answer (about 2018) "convergent bodies are not b-rep bodies, they are a new type of object called a Parasolid Mesh"; newer answer (about 2023): convergent bodies can go out "in Parasolid binary files not text files", and "Step AP242 also supports Convergent Bodies" - but the asker still got import errors in Ansys Electronics and Feko. https://community.sw.siemens.com/s/question/0D54O000061xKO5SAM/exporting-convergent-bodies-as-parasolid-or-stp-files ; Eng-Tips user: convergent tools are nice "only if you stay within NX". https://www.eng-tips.com/threads/export-convergent-body-as-a-step-iges-parasolid-etc.510654/
- Siemens forum (about 2017): "STEP only supports b-rep type bodies of solids and surfaces, not tessellated bodies like facet body commands" (pre-AP242-tessellation answer). https://community.sw.siemens.com/s/question/0D54O000061x7bnSAA/convert-facet-body

Summary for NX:
1. Master: implicit definition (equation/field) inside the implicit task environment; result handed downstream as convergent (facet) geometry. Exact mechanism GAP (no official doc page reached).
2. Sync: implicit features live in the NX feature tree of the task environment; output convergent body for downstream. Associativity of output to the implicit recipe: GAP.
3. Users: TPMS heat exchangers, implant fills, lattice-filled AM parts.
4. Limits: output is faceted (convergent), not exact analytic B-rep; convergent body exchange outside NX is unreliable per users. No source found claiming exact STEP B-rep with true cylinders from implicit bodies.
5. Formats: Parasolid binary (convergent), STEP AP242 (convergent/tessellated), JT; AM build prep inside NX AM.

## 3. Altair Inspire - Implicit Modeling (Altair now part of Siemens per CIMdata 2026 list "Altair (Siemens)")

- Convert to Implicit: inputs "BRep (NURBS), such as Parasolid and STEP", PolyNURBS, polymesh (.stl/.obj/.3mf), and topology optimization results. "The output is a Signed Distance Field (SDF)"; surfaces/curves give an unsigned distance field (UDF). Visualization Quality Low..Very High controls sharpness ("Medium quality shows rounded edges, while Very High preserves sharper features"). Page does not state whether the implicit body stays linked to the source B-rep. https://help.altair.com/inspire/en_us/topics/implicit/convert_t.htm (current help, 2026)
- FAQ: B-Rep and PolyNURBS "converted on the fly"; "prefers a GPU based system"; NVIDIA RTX 16 GB suggested; export STL, OBJ, 3MF; slice export "Export implicit parts as slice files (.cli and .3mf) to slice the native implicit model without meshing"; Python API with "75+ documented demos". https://help.altair.com/inspire/en_us/topics/implicit/implicit_faq.htm
- Field-driven: surface lattice shell thickness, strut diameter, unit cell size, stochastic valence all field-drivable; fields from OptiStruct/Fluids results. https://help.altair.com/inspire/en_us/topics/implicit/fields_thinking_c.htm ; https://help.altair.com/inspire/en_us/topics/implicit/surface_t.htm
- Back to CAD: "The route back to parametric CAD is via the PolyNURBS tool", fitting subdivision surfaces; point-edge sets page: "While we don't recommend it for complex models, you can fit a PolyNURBS surface to an implicit body". https://help.altair.com/inspire/en_us/topics/implicit/point_edge_sets_c.htm
- Altair Community moderator: to export STEP/Parasolid it "must be expressed with NURBS"; not recommended for lattices because the patch count makes the model "unworkable"; there is "no standard file format for implicit functions". https://community.altair.com/discussion/62855/how-to-save-or-export-an-implicit-model-in-desired-cad-extension-such-as-stp-file ; https://community.altair.com/discussion/26306/what-is-the-best-way-to-convert-an-implicit-body-to-a-solid-for-export (dates approx. 2023-2024, from search snippets)
- Inspire general export: Parasolid 25.0 to 37.1, STEP AP203. https://help.altair.com/inspire/en_us/topics/shared/unity_basics/file_formats_st_r.htm

Summary for Inspire:
1. Master: implicit SDF is the working representation in the implicit environment; B-rep inputs converted on the fly.
2. Sync: implicit modeling history in Inspire (construction-history style); B-rep source to implicit conversion; live link to source B-rep GAP.
3. Users: lattice-filled AM parts, topology-optimised parts.
4. Limits: no direct STEP/Parasolid export of implicit bodies; only PolyNURBS (subdivision) fit, which is approximate and impractical for lattices. No analytic cylinders.
5. Formats: STL, OBJ, 3MF meshes; CLI and 3MF slices direct from implicit; B-rep only via PolyNURBS fit.

## 4. nTop (formerly nTopology)

Kernel:
- Own implicit (F-rep / SDF) kernel. Third-party: "nTopology has their own geometry kernel. They don't rely on Parasolid, ACIS..." https://revolutioninsimulation.org/tag/ntopology/
- nTop 5.0 (released 2024-06-24): new implicit modeling kernel "built from the ground-up" (vendor claim), "enables greater precision and accelerates common design operations" (vendor claim); support article lists better precision for thin TPMS walls and "micron-scale tolerance for machined parts" (vendor claim). https://develop3d.com/cad/ntop-5-adds-interoperability-for-major-systems/ (2024-06-25); https://support.ntop.com/hc/en-us/articles/26062971882131-nTop-5-0-New-Implicit-Modeling-Kernel (403 on fetch; content from search snippet)
- nTop's data model "is a signed distance field-based implicit model - which runs natively in GPUs" (NVIDIA partnership coverage, 2024/2025; includes NVentures investment). https://www.3dprintingindustry.com/news/ntop-and-nvidia-enhance-computational-design-in-new-partnership-232786/ ; https://www.fabbaloo.com/news/ntop-partners-with-nvidia-to-integrate-gpu-powered-design-in-additive-manufacturing

Parasolid inside nTop:
- nTopology 3.39 release notes (about 2023): "We upgraded our Parasolid kernel to v35 for enhanced CAD import capability." So Parasolid is embedded for B-rep import/export, not as the modeling kernel. https://support.ntop.com/hc/en-us/articles/12882714639507-nTopology-3-39-What-s-New
- File types (article applies to nTop 5.45): import Parasolid .x_t/.x_b 9.0-38.0.163, STEP AP203/214/242, NX .prt 11 to NX 2412; "Parasolid formats tend to yield the most dependable results". CAD export Parasolid and STEP 242; native "nTop Implicit File, *.implicit". https://support.ntop.com/hc/en-us/articles/360048249514-What-file-types-can-nTop-import-and-export (403 on fetch; from search snippets)
- Import workflow: CAD imported then converted to implicit body ("less labour-intensive ... to run computations on these implicit bodies"). https://develop3d.com/cad/ntopology-3-0-review/

Export to CAD (B-rep):
- Two routes. (a) Mesh -> quad mesh -> refine -> "CAD Body from Quad Mesh" -> Export Part (STEP or x_t); for simple models such as topology-optimisation results. (b) "CAD Body from Implicit Body (Beta)" (since nTop 3.28) generates "a smooth (G1-continuous) solid part directly from an Implicit Body"; designed for small/medium lattices, rough limits "around 10,000 beams or faces for graph lattices and up to around 1,000 unit cells for TPMS". nTop recommends Parasolid for NX/SolidWorks. Example: lattice cube .step "over 600MB" vs .ntop "less than 30kB". https://support.ntop.com/hc/en-us/articles/7323484481555-How-to-export-geometry-as-a-CAD-part-STEP-or-Parasolid ; https://support.ntop.com/hc/en-us/articles/7323109441683-nTopology-3-28-What-s-New
- Output surfaces are fitted (G1 smooth patches from quad mesh or implicit body), not analytic recognition. No source found saying nTop recovers true planes/cylinders from implicit geometry. "Preserved CAD" input exists on export (from search snippet; exact semantics GAP - may let original imported B-rep faces pass through).
- Simplified body export: "simple CAD envelope with the mass properties of implicit designs" for downstream CAD. https://support.ntop.com/hc/en-us/articles/7170769142291-How-to-export-a-Simplified-Body-Representation-of-Implicit-as-a-CAD-part
- Platform page (vendor): exports "STEP, Parasolid, STL, and 3MF", "direct mesh export for analysis tools and Convergent Body for modern CAD systems". https://www.ntop.com/platform/ (accessed 2026-10-21)
- Users report: complex lattice export "seemed impossible"; exported geometry "deviated significantly". https://www.researchgate.net/post/What_is_the_best_way_to_export_complex_lattice_from_nTopology_software_into_Ansys_Worckbench (older)
- To NX: support article "How to Export Curves from nTop to NX for CAD reconstruction" - a rebuild-in-NX workaround. https://support.ntop.com/hc/en-us/articles/53668810987667-How-to-Export-Curves-from-nTop-to-NX-for-CAD-reconstruction

Export to printers:
- Implicit Interop: .implicit file; vendor claims up to 99% smaller than meshes, up to 500x faster generation, "accurate and lossless". EOS plug-in for EOSPRINT 2.14 (announced 2023-04-21); 3MF Consortium standardisation work. Siemens Energy heat exchanger printed on EOS M 290 from <1 MB implicit file. https://www.tctmagazine.com/ntopology-and-eos-bring-implicit-interop-to-3d-printing/ ; https://www.additivemanufacturing.media/news/eos-ntop-create-file-capability-to-limit-design-data-bottlenecks
- nTop Core (developer library, 2023; updated to new kernel with nTop 5): partners read/query .implicit natively. nTop 5 partners: Materialise (Magics, NxG Build Processor), Autodesk (Fusion "nTop Connector" add-in), Hexagon (Cradle CFD scSTREAM), Intact Solutions, cloudfluid. https://www.metal-am.com/ntop-announces-ntop-5-release-and-adds-five-integration-partners/ ; https://develop3d.com/cad/ntop-5-adds-interoperability-for-major-systems/
- Materialise: nTop kernel "optional add-on to all new Next-Gen Build Processors" (TCT interview). https://www.tctmagazine.com/ushering-in-a-new-era-ntop-ceo-brad-rothenberg-nvidia-materialise/
- SimScale (June 2025) imports nTop implicit geometry natively via immersed boundary. https://www.simscale.com/press/simscale-and-ntop-integration-announcement/

Export to CNC:
- No nTop-native CAM or toolpath output found. Route found: Autodesk Fusion integration (2023-11-03) where "CAM features expose nTop functions to CNC users" (article wording). https://www.tctmagazine.com/autodesk-integrate-ntop-implicit-interop-into-fusion-cad/ . Otherwise CNC goes through exported STEP/Parasolid (fitted) B-rep.

Summary for nTop:
1. Master: implicit field (SDF) in nTop's own kernel; imported B-rep converted to implicit.
2. Sync: notebook/block graph re-evaluated from recipe; B-rep and meshes are derived at export.
3. Users: TPMS heat exchangers, lattices, implants, cold plates, mostly AM.
4. Limits: B-rep export only for small/medium complexity (stated limits), fitted G1 surfaces, huge STEP files; no evidence of analytic face output.
5. Formats: .implicit (Implicit Interop), STL/3MF/OBJ mesh, STEP 242 and Parasolid, Convergent Body; printing via EOS/Materialise plug-ins.

## 5. Autodesk Fusion - Volumetric modeling / Volumetric Lattice / Volumetric API

- API docs (Fusion October 2026 API Help): the API "lets you define and edit implicit geometry programmatically" (namespace adsk.volume); a model is "a pair of graphs" whose output channels combine into "one signed distance field" which can be converted to mesh geometry; graph editing has "no equivalent user interface in Fusion". https://help.autodesk.com/cloudhelp/ENU/Fusion-360-API/files/VolumetricModeling_UM.htm
- From search snippets of the same doc: cell graph has lattice cell (whole volume) and texture cell (surface only); a BoundarySDF channel limits the lattice to the boundary body; "Implicit models lack explicit edges and faces. The exact geometry is determined at evaluation time, influenced by resolution." Mesh via VolumetricModelToMeshFeature: "Marching Cubes" or a slower method that "preserves the sharp features on the boundary of the boundary body".
- Volumetric Lattice feature (UI) generates cells "as an appearance"; must "Convert Volumetric Model to Mesh" before STEP/STL export carries it. https://www.autodesk.com/support/technical/article/caas/sfdcarticles/sfdcarticles/Exported-Volumetric-Lattice-part-from-Fusion-360-doesn-t-contain-Lattice-element.html ; https://forums.autodesk.com/t5/fusion-support-forum/no-way-to-export-volumetric-lattice/td-p/11574198
- Help: "only apply one volumetric lattice per component"; latticed bodies "can only be used in the Additive setups"; otherwise convert to mesh body. https://help.autodesk.com/view/fusion360/ENU/?contextId=SLD-VOLUMETRIC-LATTICE
- AU 2022 handout: conversion is not parametric - "editing the volumetric lattice will NOT update the mesh body". https://static.au-uw2-prd.autodesk.com/Class_Handout_CP502032_ClassHandout-CP502032-Lichtman-AU2022_v05.pdf
- Users: to get SAT/STEP solids, mesh then convert using "faceted" method; heavy files; 100k-1M+ elements. https://forums.autodesk.com/t5/fusion-support-forum/how-to-convert-volumetric-lattice-to-solid-body/td-p/12351038
- nTop Implicit Interop integration announced 2023-11-03 (release planned early 2024); nTop Connector add-in in Fusion App Store (2024). https://www.tctmagazine.com/autodesk-integrate-ntop-implicit-interop-into-fusion-cad/

Summary for Fusion:
1. Master: B-rep body is master; volumetric model is an SDF layered on it (BoundarySDF from the B-rep body).
2. Sync: volumetric features re-evaluate with the B-rep in the timeline; mesh conversion is a one-way, non-associative snapshot.
3. Users: lattice-filled printed parts in Additive setups.
4. Limits: lattice is display-only until meshed; one per component; no implicit-to-exact-B-rep; STEP only as faceted solid.
5. Formats: mesh (STL/3MF/OBJ) after conversion; additive build output from Additive setups; STEP of faceted solids.

## 6. Dyndrite

- Launched AMUG 2019-03-31: Accelerated Geometry Kernel (AGK), vendor claim "world's first fully GPU-native geometry kernel"; "Hybrid" kernel for "splines (NURBs), surface tessellations, volumetric data, tetrahedra, and voxels". https://www.dyndrite.com/press-releases/dyndrite-corporation-reveals-new-accelerated-geometry-kernel-the-worlds-first-fully-gpu-native ; https://www.voxelmatters.com/dyndrite-first-gpu-native-geometry-kernel-additive-toolkit/
- Later called Accelerated Computation Engine (ACE): "hybrid geometry kernel capable of hosting any geometry type without demoting to lower quality triangles" (vendor). Users can switch representations "without destroying the original underlying data format" (vendor). C++ and Python APIs.
- Toolpathing API: maintain "native CAD-based spline and B-Rep data throughout the entire workflow"; CAD colour metadata drives print strategies; GPU voxel engine for 3D geometric queries. https://www.dyndrite.com/technology/raster-toolpathing-api
- LPBF Pro: Native CAD Import with Metadata, supports, slicing, direct-to-print for Aconity3D, EOS, Renishaw, SLM etc.; CLI/CLI+, STL/3MF, OVF/ILT; scripted "build recipe". https://www.dyndrite.com/lpbf-product
- Summary: B-rep (native CAD) master for build prep; voxel/implicit used for queries and lattice/toolpath compute. Not a design CAD; output is machine build files, not STEP. Implicit-to-B-rep: not applicable / GAP. Independent benchmarks: GAP.

## 7. PTC (Creo, Onshape)

- Creo+ "IWP - New Formula-Driven Lattice Cell Type" (TPMS) under Lattice > Formula Driven. https://support.ptc.com/help/creo/creo_plus_early_visibility/usascii/whats_new/core_lattice_iwp.html (date GAP). PTC calls these "formula driven" not "implicit". Export behaviour of formula-driven lattices to STEP: GAP.
- Onshape: April 2026 demo at Siemens PLM Components of FeatureScript + MCP server producing "convergent parts where boundary representation features and AI generated implicit geometries coexist", on Parasolid. CIMdata commentary (see section 1). Shipping status: GAP.

## 8. Dassault Systemes (CATIA, SOLIDWORKS, 3DEXPERIENCE), Spatial

- CATIA Lattice Designer (3DEXPERIENCE): GoEngineer: "uniquely combines the strengths of our traditional Boundary Representation (BREP) CAD tools with new and advanced Signed Distance Field (SDF or Implicit Modeling) design tools"; recognises "Volumes" as a separate entity type, which restricts the features that work on them. https://www.goengineer.com/blog/3dexperience-lattice-designer-technical-deep-dive
- Dassault community wiki: "robust and scalable complex lattice design without limitations on model size, visualization or number of cells" (vendor), conformal bar lattices, TPMS "leveraging gpgpu compute capabilities for instant preview", heat-exchanger tools; 2025 GA: TPMS material distribution tuning. https://3dswym.3dexperience.3ds.com/wiki/catia-user-community/catia-lattice-designer_KTE65-ffSyCYRcKYv1ms7A ; https://3dswym.3dexperience.3ds.com/wiki/catia-user-community/2025-ga-catia-lattice-designer_Jrhj3_FcToqxBfhwYa00SA
- MDPI review (2025/26): CATIA Lattice Designer supports TPMS "with associative CAD updates". https://www.mdpi.com/2504-477X/10/1/32 (403 on fetch; from snippet)
- SOLIDWORKS desktop implicit features: GAP (GoEngineer says SOLIDWORKS "only knows 'solids' and 'surfaces'").
- Spatial ACIS/CGM: release posts 2025 1.0, 2025 1.0.1 (2025-04-30), 2026 1.0 (about 2025-11) mention mesh prep, polyhedral bodies, thickening polyhedral sheets for AM (2024 1.0.1) - no implicit, SDF or lattice features found. https://blog.spatial.com/news/2026-1-0 ; https://blog.spatial.com/news/2025-1-0 ; https://blog.spatial.com/news/2024-1-0-1

## 9. Ansys (Discovery, SpaceClaim)

- SpaceClaim lattice/infill is facet-based: gyroid created "as faceted geometry (like an STL), which can then be solidified into a CAD model (STEP/IGES)", hard to manage due to face count. https://www.fastwayengineering.com/blog/gyroid-lattice-heat-exchanger-a-fluid-flow-study ; lattice tools in the Facets menu (SpaceClaim 19.2, 2018-12-04). https://www.padtinc.com/2018/12/04/how-to-use-lattice-optimization-in-ansys-mechanical-and-ansys-spaceclaim-19-2/
- No implicit/SDF representation found in Discovery/SpaceClaim. Implicit (nTop) gyroids reach Ansys as STL/OBJ. https://www.mdpi.com/2673-4591/160/1/7 ; https://innovationspace.ansys.com/forum/forums/topic/any-body-know-about-ntopology-software-and-how-we-can-import-geometry-from-ntop-to-ansys-workbench/

## 10. Hexagon, Materialise, Carbon, Hyperganic, Synera

- Hexagon: Cradle CFD scSTREAM integrates nTop Core for direct implicit transfer (vendor claims "reduces manual processes by 80%", "5-10x faster"). https://develop3d.com/cad/ntop-5-adds-interoperability-for-major-systems/ . Hexagon presented at Siemens PLM Components 2026 (CIMdata). Own implicit CAD: GAP.
- Materialise: Magics and NxG Build Processor take nTop's implicit kernel (above). 3-matic "extends lattice capabilities to the mesh level" (MDPI review snippet) - mesh-based. Materialise presented at Siemens PLM Components 2026.
- Carbon Design Engine: "uses the same kernel as Carbon's powerful internal advanced lattice design tools"; conformal lattices, five lattice types (vendor). https://www.carbon3d.com/resources/blog/carbon-advanced-lattices . Listed as 3MF importer. https://3mf.io/compatibility-matrix/ . Representation and B-rep export: GAP.
- Hyperganic: engine "fundamentally based on voxels"; Core 1 direct output to printers incl. G-code; BASF Forward AM Ultrasim 3D Lattice Engine uses Hyperganic's voxel engine; HyDesign exports 3MF, STL, PNG and CLI slices. No STEP/B-rep export found. https://jlk.ae/2022/06/15/3-cores-now-its-your-turn/ ; https://3dprint.com/299689/basf-forward-am-to-debut-ultrasim-3d-lattice-engine-software-at-rapid-tct/ ; https://www.voxelmatters.com/hyperganic-launches-hydesign-for-3d-printed-lattices/
- Synera: presented at Siemens PLM Components 2026 (CIMdata list). Implicit features and kernel: GAP (search returned nothing on Synera).

## 11. Cross-cutting observations stated by sources

- No source in this sweep claims that any commercial system converts general implicit geometry into an exact STEP B-rep with analytic planes, cylinders and circles. Every documented route produces either mesh, faceted/convergent bodies, or fitted NURBS/subdivision (nTop G1 fit, Inspire PolyNURBS fit, Fusion faceted solid).
- Parasolid's stated direction (CIMdata, April 2026) is to avoid conversion by keeping implicit content inside the kernel next to B-rep and facets ("deferring discretization for as long as possible"); persistent implicit objects are announced as future.
- Exact analytic B-rep survives only where the B-rep stays master and the implicit part is an attached lattice/infill (Parasolid cellular bodies, Fusion BoundarySDF, CATIA volumes) or where the original B-rep is passed through on export (nTop "Preserved CAD" - semantics unconfirmed).

## Gaps and uncertain items

- Parasolid v36.x, v37.x release notes on implicit: not found.
- Parasolid public API names for implicit (PK_*): not public.
- NX implicit output exact options and associativity: no official doc page reached.
- nTop "Preserved CAD" export input semantics; nTop Convergent Body export mechanics.
- nTop native CNC/CAM: none found.
- Creo formula-driven lattice export, Onshape implicit shipping status, SOLIDWORKS desktop, Synera, Hexagon own CAD, Carbon representation: not found.
- support.ntop.com and mdpi.com pages returned 403; their content above comes from search-result snippets, not full-page reads.
- PicoGK / LEAP 71 and Oqton/others: not researched within budget.
