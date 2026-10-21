# Research E: implicit-mode features in implicit/SDF CAD tools and how they mix with B-rep

Compiled 2026-10-21. About 40 tool calls.

Source-quality notes:
- "Vendor claim" marks statements made by the vendor or a reseller about its own product.
- support.ntop.com, forums.autodesk.com, eng-tips.com, javelin-tech.com and mdpi.com returned HTTP 403 to direct fetch. Facts tagged [search-index] come from the search engine's rendering of those pages, not from a page I read directly. Treat them as one step less reliable than a direct read.
- Dates are the page's own date where shown, otherwise the release or version the page names, otherwise "date unknown".

---

## 1. nTop (nTopology)

### 1.1 Implicit features and how users drive them
- Interface is a "notebook" of blocks (node-graph style, linear notebook). nTop staff describe the modelling approach as "a variation of constructive solid geometry (CSG) with Signed-Distance Fields (SDFs)". LinkedIn post by Hans Gruber (nTop), roughly mid-2025. https://www.linkedin.com/posts/hans-gruber-29997336_why-implicit-modeling-is-reshaping-engineering-activity-7333367691329597441-R8S7
- Ramp block: input scalar field can be an implicit body, a point or axis, or a simulation result such as von Mises stress or displacement; values clamp below In min / above In max; Continuity options Geometric, Tangential, Curvature. nTop 5.0 caution: do not give Ramp a zero-length range, use If-Else instead. [search-index] https://support.ntop.com/hc/en-us/articles/360041676813-How-do-I-use-the-Ramp-block (date unknown, 5.x era)
- Remap Field: Fields > Remap; Scalar Field plus Coordinate Space inputs; produces an infinite field that must be Boolean-intersected with a bounding box. nTop 5.43 added a new default overload taking a vector field, plus ONNX neural-network surrogate import and a "Track to Block" tool that jumps from viewport to notebook. [search-index] https://support.ntop.com/hc/en-us/articles/50080923413907-nTop-5-43-What-s-New (date unknown) and https://support.ntop.com/hc/en-us/articles/1500004072802-How-to-use-Remap-Field-to-scale-or-translate-an-object
- Simulation results to field: Simulation tab > Exchange; build a von Mises Point Map, then Field from Point Map; displacement gives a vector field, use the Length/X/Y/Z chip to get a scalar for Ramp. [search-index] https://support.ntop.com/hc/en-us/articles/360054149874-How-to-use-simulation-results-to-create-a-Point-Map-or-Field
- Training course "Intermediate Field-Driven Design" covers Boolean operations with fields, capping TPMS with ramps, remap fields with Booleans for planar, cylindrical and spherical woven patterns. https://learn.ntop.com/courses/recorded-training-intermediate-field-driven-design/
- nTop 5.0: scalar field and implicit body values "may differ" from earlier versions because of higher precision; any block taking or outputting Scalar Field, Vector Field or Implicit Body may give different results. [search-index] https://support.ntop.com/hc/en-us/articles/27281187796883-nTop-5-0-Outline-of-expected-differences
- Not confirmed in this pass: exact names of smooth Boolean / blend blocks, conformal lattice blocks, shell/thicken blocks (training pages imply they exist; no direct read).

### 1.2 B-rep operations on or after an implicit body
- nTop has no in-tool B-rep editing path documented in what I read. Its own FAQ positions nTop for early-stage generation and says teams "export to the CAD, CAE, and PLM tools they already trust for detailed design, final validation, and documentation". Vendor claim of "minimal handoff friction". [search-index] https://www.ntop.com/platform/
- Direct implicit import into Autodesk Fusion ("nTop for Fusion" add-in), showcased Formnext 2023: implicit designs import into Fusion without meshing first. Vendor claim. https://www.autodesk.com/products/fusion-360/blog/ntop-and-autodesk-fusion-join-forces-direct-implicit-model-import-will-be-showcased-at-formnext-2023/
- SimScale integration to import nTop implicit geometry, June 2025; Siemens Energy named as a user. https://engtechnica.com/simscale-integrates-with-ntop-to-import-implicit-geometry/

### 1.3 Outputs
- Mesh: Export Mesh block (STL, OBJ, PLY), Export 3MF block. "In most cases you must first convert the implicit body into a mesh", and the mesh needs refining first. https://learn.ntop.com/courses/101-ntop-essentials/lessons/exporting-parts/
- CAD: Export Part block (.STEP or .x_t), after converting to a CAD body; overload must match input type. Parasolid recommended for SolidWorks and NX so no translation is needed. Same page; [search-index] https://support.ntop.com/hc/en-us/articles/7323484481555-How-to-export-geometry-as-a-CAD-part-STEP-or-Parasolid
- Slices: Manufacturing tab offers multiple slice formats (names not given on the page). https://learn.ntop.com/courses/101-ntop-essentials/lessons/exporting-parts/
- Mass properties on the implicit itself: Mass Properties block (Modeling > Evaluation), density input, Relative Error default 0.01. Surface Area block accepts implicit body, mesh or profile. https://learn.ntop.com/courses/101-ntop-essentials/lessons/evaluating-mass-volume-and-weight-savings/ ; [search-index] https://support.ntop.com/hc/en-us/articles/4402931752979-How-to-calculate-the-mass-or-volume-of-a-part ; https://support.ntop.com/hc/en-us/articles/5469662981139-How-to-calculate-surface-area
- FE meshing from implicit bodies exists (5.0 meshing updates; a 2025 paper "Best practices for generating finite element meshes from implicit bodies in nTop"). Not read in detail. https://support.ntop.com/hc/en-us/articles/30049593128211-nTop-5-0-Meshing-updates-for-Simulation-Optimization ; https://www.researchgate.net/publication/398926710_BEST_PRACTICES_FOR_GENERATING_FINITE_ELEMENT_MESHES_FROM_IMPLICIT_BODIES_IN_NTOP

### 1.4 Conversion costs and warnings
- Two implicit-to-CAD routes [search-index of the STEP/Parasolid article and 3.28 notes]:
  - Route 1 (simple models: topology-optimisation results, unfilled shells, coarse ribs): Mesh from Implicit Body > Quadrangulate Mesh > Refine Mesh > CAD Body from Quad Mesh > Export Part. NURBS align with features, fewer faces, easier to edit in CAD; "only works reliably on fairly simple models".
  - Route 2 (lattices and medium complexity): CAD Body from Implicit Body (Beta in 3.28, needs Beta blocks enabled). NURBS "not necessarily nicely arranged" but more likely to succeed. Converts "within a specified Tolerance"; larger tolerance gives fewer faces and faster generation.
  - Suggested limits: graph lattices up to about 10,000 beams or faces; TPMS up to about 1,000 unit cells. Above that, use Simplified Body Representation.
  - Example: a lattice cube exported to STEP was over 600 MB while the .ntop file was under 30 kB. nTop says export "does not overcome some inherited shortcomings of traditional CAD file formats".
  - Tutorial mesh step uses a 1 mm tolerance; Block Details shows closed / manifold / self-intersection checks.
  - Sources: https://support.ntop.com/hc/en-us/articles/7323484481555-How-to-export-geometry-as-a-CAD-part-STEP-or-Parasolid ; https://support.ntop.com/hc/en-us/articles/7323109441683-nTopology-3-28-What-s-New (3.28, 2022) ; https://www.ntop.com/resources/blog/whats-new-july-2022/ ("Lattices in... CAD?", July 2022; returned 404 on direct fetch in Oct 2026)
- Quadrangulate route, direct read: "this should not be used with lattices or TPMS structures"; "We do not recommend it for use with a complex lattice structure". No numeric accuracy or time figures given. https://learn.ntop.com/courses/102-guide-to-meshing/lessons/conversion-to-cad-with-quadrangulate-mesh/
- Simplified Body Representation (proxy for CAD): v3.15 (beta) assigns a modified density to the original CAD design space so CAD shows correct mass, optional Exact COG; v3.27 adds Match Inertia, which turns the envelope into surface bodies plus multiple solid spheres with computed densities. Output Parasolid. Export guide recommends it for graph lattices over 10k beams or TPMS over 1k cells, for visualisation, clearances, weight studies, assemblies, drawings, archiving. In nTop's NX example the bracket "visually doesn't show the lattice" but mass properties match. [search-index] https://support.ntop.com/hc/en-us/articles/7170769142291-How-to-export-a-Simplified-Body-Representation-of-Implicit-as-a-CAD-part ; https://support.ntop.com/hc/en-us/articles/4410154496915-nTopology-3-15-What-s-New ; https://support.ntop.com/hc/en-us/articles/6886452460819-nTopology-3-27-What-s-New ; https://support.ntop.com/hc/en-us/articles/360048784874-How-to-export-geometry-Overview-of-export-workflows

### 1.5 Friction reports
- ResearchGate thread: one user said converting complex structures to STEP "seemed impossible, even we use fillet"; another said a student's exported geometry deviated significantly from the nTop lattice, and the student imported a mesh into ANSYS instead. Date unknown. [search-index] https://www.researchgate.net/post/What_is_the_best_way_to_export_complex_lattice_from_nTopology_software_into_Ansys_Worckbench
- Siemens PLM Components 2026 conference: George Allen (nTop) said converting implicit bodies into the B-reps that downstream manufacturing and drafting applications require "remains a significant challenge". https://www.demystifyingplm.com/siemens-plm-components-2026-conference-report ; https://www.cimdata.com/en/component/docman/doc_download/4175-siemens-plm-components-2026-commentary

---

## 2. Altair Inspire Implicit Modeling (Altair now part of Siemens)

### 2.1 Features and controls
- Ribbon tools: primitives (cuboid, cylinder, sphere), surface, planar, strut and stochastic lattices, point clouds, Field, Convert to Implicit Geometry, Move, Combine, Subtract, Intersect, Offset, Invert, Morph, Smooth, Fillet, Remap; plus View Field and Visualization Quality. https://2024.help.altair.com/2024/inspire/en_us/topics/implicit/fields_thinking_c.htm (Inspire 2024)
- Strut lattice cell types: BCC, face-centred primitive, FCC, cubic primitive, Isotruss, Octagon, Hex Truss, Kelvin, Fluorite, Diamond, Truncated Cube; fills Parasolid, STL, PolyNURBS or implicit bodies. https://help.altair.com/inspire/en_us/topics/implicit/strut_t.htm
- Surface (TPMS) lattice density by slider/value, variable, or field. https://help.altair.com/inspire/en_us/topics/implicit/surface_t.htm
- Planar lattice cells (hexagon, diamond, rectangle, triangle, isogrid, circle) with "Conform to" that builds conformal curved spaces over U, V, W. https://help.altair.com/inspire/en_us/topics/implicit/planar_t.htm
- Stochastic lattice: beams filtered by field value. https://2024.help.altair.com/2024/inspire/en_us/topics/implicit/stochastic_t.htm
- Lattices can get an outer shell (inward, outward, both) or be combined with another body; transition options such as "Sharp".
- Field-driven: most sizes (strut diameter, wall thickness, cell size, shell thickness) take a value, a variable, or a field. Example: strut thickness 1 mm to 4 mm by distance to a plane. Fields can be made in context next to a parameter. Point clouds drive fields via inverse-distance weighting or nearest neighbour. Fields can come straight from OptiStruct or fluids results ("Create Field" at bottom of Analysis Explorer). Stated uses include varying fillet radii by stress level. https://2024.help.altair.com/2024/inspire/en_us/topics/implicit/fields_thinking_c.htm
- Python API hwx.inspire.implicit (2024) exposes primitives, lattices, point clouds, fields, offset, Booleans, invert, smooth, morph, fillet. https://help.altair.com/inspire/inspireapi/PackageReference/Inspire_/Implicit_/Implicit.html
- Convert to Implicit accepts B-rep (Parasolid, STEP), PolyNURBS, polymesh (STL, OBJ, 3MF), topology-optimisation results. https://help.altair.com/inspire/en_us/topics/implicit/convert_t.htm
- Reseller: implicit came into Inspire with the 2023 release, carrying over Altair Sulis capability. https://www.trueinsight.io/altair/altair-inspire

### 2.2 B-rep operations on implicit bodies
- Altair: "Implicit Modeling has no such features and they cannot be readily selected, highlighted, or edited"; "It is therefore difficult to select points on a surface, lines that separate features, and isolated faces." The whole surface is one isosurface. So edge fillets in the B-rep sense do not apply; Inspire's Fillet is an implicit operation. https://2024.help.altair.com/2024/inspire/en_us/topics/implicit/fields_thinking_c.htm
- Altair warns implicit combine/smooth/morph/offset may not preserve distance: "Simply combining multiple SDFs does not necessarily result in a SDF after the operation"; blends, fillets and offsets of distance-like fields may not give the exact thickness or radius wanted. Same page.

### 2.3 Outputs
- Mesh export: "easy to mesh and can be output as .stl, .obj, and .3mf". Same page.
- Structural analysis: official tutorial creates lattice volume plus overlap volume, Boolean-combines into "a single mesh part", then runs analysis. A community article says after Combine "the lattice gets a mesh... no longer just an implicit model", and the final run switches to SimSolid (meshless). https://help.altair.com/inspire/en_us/topics/inspire/tutorials/tut_implicit_struct_analysis_c.htm ; https://community.sw.siemens.com/s/article/37119-inspire-workflow-structural-analysis-topology-optimization-implicit-modeling
- Reseller claim: SimSolid means no meshing or simplification needed, including for implicit structures. https://www.trueinsight.io/blog/comparing-inspire-simlab-simsolid
- STEP: only via PolyNURBS fit (below). Drawings, interference checks: not found.

### 2.4 Conversion costs and warnings
- Grid sampling set by Visualization Quality; coarse grids round sharp edges; finish "in the highest visualization quality that your computer can tolerate"; complex models give very high triangle counts. https://2024.help.altair.com/2024/inspire/en_us/topics/implicit/fields_thinking_c.htm
- Back to CAD: PolyNURBS tool "can be used to fit subdivision surfaces to the Implicit model"; "you will be trading off computing time and accuracy during the fitting process"; "Implicit geometry cannot always convert back to other formats as easily"; "Parametric CAD is a poor choice for a very complex model." Same page.
- PolyNURBS Fit settings (2021.1): Number of PolyNURBS Faces default 2500, Curvature default 50%, Auto Creasing default No; Shrinkwrap first for better results. https://2021.help.altair.com/2021.1/inspire/en_us/topics/shared/geometry/polynurb_fit_st_t.htm

### 2.5 Friction reports
- Altair Community (about 2023): Altair staff said Fit is normally best; dense parts need more patches; for very dense or complex parts Fit "may not work" or quality may be insufficient; alternative is a 1:1 polygon-to-NURBS conversion in Inspire Studio that keeps the exact shape but also struggles with complex parts. [search-index; direct fetch returned title only] https://community.altair.com/discussion/26306/what-is-the-best-way-to-convert-an-implicit-body-to-a-solid-for-export

---

## 3. Siemens NX Implicit Modeling (plus Lattice, Convergent Modeling)

### 3.1 Features and controls
- Introduced NX 1953 as an optional module for "equation-driven structures that can't be designed through traditional CAD methods"; general equation builder with predefined examples; built-in TPMS incl. gyroid; shelling, Boolean, blend operations called "extremely robust" (vendor claim). Siemens blog, 2021-01-25. https://blogs.sw.siemens.com/additive/aa-e10/
- Forum summary: TPMS gyroid, Schwarz, Neovius, Schoen, Diamond; implicit primitives; Unite, Subtract, Intersect, Shell/Offset. [search-index] https://www.eng-tips.com/threads/nx-implicit-modeling.478264/
- Runs in an "Implicit Modeling task environment". June 2024 release (blog 2024-07-10): Set Resolution ("set the voxel size for any selected implicit bodies"); Unit Cell (a triply periodic one-body unit used by Body Lattice outside the task environment); Unite gains Continuous, Circular, Angular blend types. Lattice Designer: Voronoi pore distributions Uniform / Variable / Gaussian, Filter Lattice by thickness or ratio, Connect Dangling Rods. https://blogs.sw.siemens.com/designcenter/new-nx-jun-2024-generative-design/
- Another June release note: rod and ball diameters "controlled associatively by fields", conformal body lattice type, organic blending factor. https://blogs.sw.siemens.com/designcenter/new-june-nx-release-details/ (year not stated)
- NX 2021 added algorithmic modelling (logic editor) and Voronoi lattices. https://develop3d.com/cad/siemens-nx-2021-q1-review-cad-modelling/

### 3.2 B-rep interplay
- Lattice command (NX 11.0.2+) output is a faceted body; Unite with the part makes "a single convergent body" editable with Convergent Modeling. https://community.sw.siemens.com/s/article/introducing-lattice-in-nx
- Convergent Modeling on Parasolid holds mesh and B-rep in one body; operations can act on both, e.g. a blend where mesh and solid meet; Replace Face of Mesh; example of exact B-rep mold with mesh gyroid lattice. https://blogs.sw.siemens.com/additive/siemensconvergentvision/
- Mixing rule from forum: if B-rep and convergent bodies are combined either way, "the result will end up being faceted". Convert Facet Body "creates a facet model not a b-rep body"; "STEP export only support b-rep type bodies ... not tessellated bodies". NX 11/12 era. [search-index] https://community.sw.siemens.com/s/question/0D54O000061x0weSAA/facet-body-to-solid-model ; https://community.sw.siemens.com/s/question/0D54O000061x7bnSAA/convert-facet-body
- 2026: CIMdata says Siemens is extending Parasolid Convergent Modeling to natively support and mix implicit geometry "without forcing full format conversions", and lists Parasolid as a B-rep, facet and implicit kernel. Vendor-direction claim. https://www.cimdata.com/en/component/docman/doc_download/4175-siemens-plm-components-2026-commentary ; https://www.demystifyingplm.com/siemens-plm-components-2026-conference-report

### 3.3 Outputs, conversion costs
- Not found: how an implicit body leaves the task environment (convergent facet body or other), STEP of implicit results, drawings, mass properties, accuracy figures beyond the voxel-size setting.

---

## 4. Autodesk Fusion (Volumetric Lattice, volumetric modelling API)

### 4.1 Features and controls
- Volumetric Lattice (Design Extension). Cells: Gyroid, Cross, X-Cell, 2.5D X-Cell, Schwarz P, Schwarz D, Custom (from selected geometry). Proportions uniform/nonuniform, Scale, Move/Rotate. Solidity tab: Uniform or Gradient Along Path plus solidity %. Offset tab: pick bodies/faces/points/lines to define a region with its own solidity and a Blend Distance. https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-CREATE-VOLUMETRIC-LATTICE.htm
- Limits: one volumetric lattice node per component; one solid body or closed mesh per lattice. Same page.
- November 2024: Volumetric Lattice and Create Mesh from Lattice became parametric timeline features, editable in the parameters table and usable in configurations. https://www.autodesk.com/products/fusion-360/blog/november-2024-product-update-whats-new/
- June 2025 preview of a volumetric modelling API: implicit (SDF) definitions via DAGs of nodes over scalars, vectors, colours; PrimaryGraph (overall shape) and CellGraph (internal pattern); channels such as BoundarySDF. https://3dprintingindustry.com/news/autodesk-previews-volumetric-modeling-api-for-fusion-240899/

### 4.2 B-rep operations
- Lattice is held as volume properties on the body. Autodesk support: exporting directly "you lose the lattice"; it is not in STL; users report STEP exports return the original solid. https://www.autodesk.com/support/technical/article/caas/sfdcarticles/sfdcarticles/Exported-Volumetric-Lattice-part-from-Fusion-360-doesn-t-contain-Lattice-element.html
- Help: bodies with volumetric lattices "can only be used in the Additive setups"; for other workspaces convert to a mesh body. [search-index] https://help.autodesk.com/view/fusion360/ENU/?contextId=SLD-VOLUMETRIC-LATTICE
- Forum: lattice does not interact with B-rep, T-splines or mesh, so it cannot be subtracted from B-rep; joining needs both converted to mesh first. [search-index] https://forums.autodesk.com/t5/fusion-support-forum/unable-to-use-cut-operation-with-volumetric-lattice/td-p/12230001

### 4.3 Outputs and conversion costs
- Create Mesh (right-click body) with High/Medium/Low/Custom refinement; Element Size sets the grid; result "can be very dense", suggest Repair/Remesh/Reduce. https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-CREATE-VOLUMETRIC-LATTICE.htm
- Additive: slicing lattice bodies "much slower" than slicing a mesh; supports and orientation use only bounding geometry, so may add unneeded supports or ignore lattice overhangs. Same page.
- No direct smooth STEP/SAT route; only workaround is mesh then faceted B-rep, which a forum expert says is very slow to make and to use. Low-accuracy meshes 100k to 1M elements; higher quickly exceeds 1M; low resolution makes the gyroid lose smoothness. [search-index] https://forums.autodesk.com/t5/fusion-support-forum/how-to-convert-volumetric-lattice-to-solid-body/td-p/12351038 ; https://forums.autodesk.com/t5/fusion-support-forum/no-way-to-export-volumetric-lattice/td-p/11574198 ; https://forums.autodesk.com/t5/fusion-design-validate-document/a-way-to-export-volumetric-lattices-as-body/td-p/12631629

### 4.4 Friction
- Users wanting the lattice in Abaqus or STEP-only tools report no good native route and a heavy double-meshing load. [search-index] same forum threads.

---

## 5. CATIA / 3DEXPERIENCE Lattice Designer

### 5.1 Features and controls
- Combines B-rep with SDF tools; bars-and-points lattices (rectangular or circular arrays), TPMS (Gyroid, Neovius, Fischer-Koch, Diamond), conformal TPMS with varying cell sizes and thickness. Guided Lattice Design Assistant: design set, model, cavity, lattice definition, generation. Cavities via Volume Extrude, Partition, Multiple Extract. Default BCC 5 mm cell, 0.5 mm bar radius. Visual Script Designer gives node-based custom specs. GoEngineer (reseller), Mike Britton, 2025-04-07. https://www.goengineer.com/blog/3dexperience-lattice-designer-technical-deep-dive
- 2025 GA (wiki 2025-01-27): Variable distribution, cell size or thickness varied along a reference curve (Relative Cell Size / Relative Thickness); cavities from intersecting surfaces with auto-extrapolation; replication modes Rectangular, Circular, Conformal between faces, Conformal through sections, Custom. https://3dswym.3dexperience.3ds.com/wiki/catia-user-community/2025-ga-catia-lattice-designer_Jrhj3_FcToqxBfhwYa00SA
- R2026x FD02: "Custom" cell type from user-defined equations; asymmetric offset of TPMS for heat exchangers; GPGPU instant preview (vendor claim). [search-index] https://3dswym.3dexperience.3ds.com/wiki/catia-user-community/catia-lattice-designer_KTE65-ffSyCYRcKYv1ms7A

### 5.2 B-rep interplay
- Reseller claim: lattice creation is non-destructive, "There is no degradation to a mesh and maintains full associativity to your original designed part"; changing part dimensions updates the lattice. But final detailed-design output is polyhedral after decimation. https://www.goengineer.com/blog/3dexperience-lattice-designer-technical-deep-dive
- Not found: whether B-rep fillets or holes can be applied to the generated lattice.

### 5.3 Outputs
- STL, AMF, 3MF exports "guaranteed to be watertight" (vendor claim). Automatic FE model generation from the lattice for structural validation (needs a SIMULIA role); usable in Fluid Dynamics Engineer for CFD; feeds DELMIA print prep. [search-index] https://3dswym.3dexperience.3ds.com/wiki/catia-user-community/catia-lattice-designer_KTE65-ffSyCYRcKYv1ms7A
- Drawings, print prep, simulation, CAM downstream on platform; concept weight-saving estimate (example 51.62%). https://www.goengineer.com/blog/3dexperience-lattice-designer-technical-deep-dive

### 5.4 Conversion costs
- Decimation sag controls roughness when going to polyhedral; higher smoothing means more deviation; finer settings raise processing time. No numbers. Same GoEngineer page.

---

## 6. LEAP 71 PicoGK (and ShapeKernel)

- Open-source voxel kernel on OpenVDB (Apache-2.0), with ShapeKernel as an abstraction layer. Code-only: C# (community Python binding PycoGK). https://github.com/leap71/PicoGK/discussions/3 ; https://github.com/ghedo44/PycoGK/blob/main/README.md
- Implicits enter via SDF callbacks rendered into voxels (Voxels_RenderImplicit); v1.7.7.5 added IBoundedImplicit interface. https://github.com/leap71/PicoGK/releases ; https://deepwiki.com/leap71/PicoGKRuntime
- Outputs: STL via Mesh from Voxels then SaveToStlFile (Lin Kayser: none of the ShapeKernel examples export STL, you add it). https://github.com/leap71/PicoGK/discussions/63
- .CLI slice export for LPBF (EOS etc.) from v1.7.0, Marching Squares vectorisation; v1.7.1 fixed a Y-flip bug that mirrored printed parts. https://github.com/leap71/PicoGK/discussions/41 ; https://github.com/leap71/PicoGK/releases
- VDB native I/O. LEAP 71 downloads ship as STL plus .CLI only. https://leap71.com/downloads/
- B-rep: no STEP or B-rep export found in any source checked; no B-rep operations exist in the kernel.

## 7. Metafold

- Kernel stores shapes as implicit signed-distance shape functions, not meshes. Shapes are DAGs of operator nodes; every graph starts with GenerateSamplePoints and ends with Redistance and Threshold. Python SDK, jobs such as evaluate_metrics (relative density) and export_triangle_mesh. https://docs.metafold3d.com/core-concepts
- SampleSurfaceLattice (gyroid etc.), SampleLattice (beam graph with nodes, edges, section_radius); Shell-and-Fill workflow (fill with lattice, shell, union). https://docs.metafold3d.com/creating-a-surface-lattice ; https://docs.metafold3d.com/shell-and-fill
- Exports STL, 3MF, PNG slice stacks (older 2022 coverage). Fusion plugin imports implicit/volumetric models into Fusion. https://www.fabbaloo.com/news/rapidly-design-lightweight-parts-with-metafolds-3d-modeling-tool ; https://apps.autodesk.com/en/Publisher/PublisherHomepage?ID=VV4J5YX8LTFQAK2Q
- Now markets "Geometric Intelligence" platform, local or cloud kernel deployment. https://www.metafold3d.com/
- B-rep: nothing found on STEP export or B-rep operations.

## 8. Hyperganic

- Core 3 public release 2022-06-15; voxel kernel, much in assembly; users write C# in Visual Studio; simulation via Finite Cell Method after acquiring DirectFEM (April 2022). https://www.tctmagazine.com/hyperganic-public-rollout-ai-powered-design-software/ ; https://3dprintingindustry.com/news/hyperganic-launches-its-ai-based-algorithmic-design-software-to-the-public-210883/
- PitchBook lists Hyperganic "Out of Business", last deal "Liquidation", no date. https://pitchbook.com/profiles/company/222979-51
- An MDPI review mentions a browser-based "HyDesign" launched 2024; status unverified (MDPI page 403). [search-index] https://www.mdpi.com/2504-477X/10/1/32

## 9. Carbon Design Engine

- Cloud lattice generator; input part STL, output latticed STL. Conformal lattices (no hanging nodes, vendor claim), multi-zone with blended zones (Pro). Tools: Adaptivity (recommended max cell size from thickness), Add Shape (zones), Analysis (adaptivity feedback, not FEA), Apply Base, Automation via API, Boolean Operations, Combine. https://learn.carbon3d.com/software/tools-and-features ; https://www.carbon3d.com/products/carbon-design-engine ; https://3dprint.com/279381/carbon-design-engine-generates-conformal-lattices-for-3d-printing/ (2021)
- Opened to all buyers late 2021. https://www.carbon3d.com/news/press-releases/carbon-makes-design-engine-software-available-to-all
- No B-rep path, no STEP, no 3MF, no FEA found.

## 10. Newcomers 2024-2026

- Spherene SphereneNXT, browser-based, launched 2026-06-30: ADMS (adaptive density minimal surface), TPMS infills, field-driven designs with continuously varying cell size, density, thickness, orientation; upload CAD or mesh incl. multi-body assemblies. Also runs inside nTop, Fusion, Grasshopper. https://www.metal-am.com/spherene-launches-browser-based-spherenenxt-for-internal-structure-design/ ; https://www.tctmagazine.com/spherene-launches-browser-based-platform-for-3d-printed-internal-structure-design/
- Spherene joined 3MF Consortium (May 2026); the 3MF Volumetric Extension supports voxel and implicit (SDF) representations. https://3mf.io/news/2026/05/spherene-joins-the-3mf-consortium-to-bring-intelligent-internal-structures-to-additive-manufacturing/
- implicit-FreeCAD (experimental fork): field graphs, B-rep signed-distance snapshots, mixed field composition via an OCCT extension; README says it "does not yet make implicit geometry interchangeable with every Part/B-Rep operation"; a tolerance-controlled B-spline bridge is unfinished. https://github.com/kevinsmia1939/implicit-FreeCAD
- implicitModeler (GitHub), code-first SDF modeller. https://github.com/jacksonwhite19/implicitModeler
- Dyndrite hybrid GPU kernel (NURBS, tessellations, volumetric, voxels); article date unclear. https://www.tctmagazine.com/dyndrite-launches-accelerated-geometry-kernel-additive/

## 11. Cross-tool and academic

- Zhao and Zou, "Translating TPMS models to STEP files", arXiv 2405.07946 (May 2024, rev. Feb 2025): F-rep vs B-rep mismatch has meant "limited use of commercial CAD/CAM/CAE software packages" for TPMS; translation must bound deviation (2 epsilon) and keep C2 continuity; meshing loses sharp features. https://arxiv.org/abs/2405.07946
- Explainer: feature trees, drawing views, tolerance annotations and NC toolpaths all live in B-rep. https://zanerobotics.substack.com/p/implicit-modeling-in-3d-math-mechanics

## Gaps (not found or uncertain)

- nTop: exact smooth Boolean / blend / conformal / thicken block names; whether CAD Body from Implicit Body has left Beta in 5.x; any timing numbers.
- NX: what body type leaves the Implicit Modeling task environment; STEP, drawings, mass properties for implicit results; accuracy costs beyond voxel size.
- CATIA: whether B-rep fillets/holes apply to the lattice; STEP of lattices; the R2026x FD02 items came only via search index.
- Fusion: whether the 2025 volumetric API left preview; any 2026 changes.
- Carbon, PicoGK, Metafold: no STEP or B-rep route found (absence, not a confirmed refusal).
- Hyperganic: shutdown date; HyDesign status.
- First-hand user reviews (G2, Reddit) on drawings or GD&T for implicit parts: none found.
- Interference or clearance checks on raw implicit bodies: no tool documented one, apart from nTop's Simplified Body envelope used for clearances in CAD.
