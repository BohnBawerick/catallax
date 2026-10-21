# Research B: pure implicit / SDF / voxel geometry kernels for engineering parts

Compiled 2026-10-21. Every fact carries a source URL and a date where known.
"[vendor claim]" or "[author claim]" marks marketing or self-reported statements.
"[inferred]" marks my own reading of a source, not a direct statement.
Repo metadata (created, last push, stars, licence) comes from the GitHub API, read 2026-10-21.

---

## 1. PicoGK (LEAP 71)

### What it is
- Open-source geometry kernel by LEAP 71 (Dubai), written by Lin Kayser, announced at GITEX Global 2023-10-18.
  https://leap71.com/2023/10/18/introducing-picogk-an-open-source-kernel-for-computational-engineering/
- Licence Apache-2.0 for both PicoGK and PicoGK Runtime ("allows both commercial and non-commercial use").
  https://leap71.com/picogk/ ; https://github.com/leap71/PicoGK
- Language: runtime in C++, higher-level API and LEAP 71 ShapeKernel in C#; viewer on GLFW; "PicoGK incorporates other Free and Open-Source projects, notably OpenVDB."
  https://leap71.com/picogk/
- Representation: voxel field on OpenVDB sparse volumes. Kayser's book: PicoGK "uses voxels as geometry representation to avoid the pitfalls of building a fragile stack of complex floating point math." Also described as "a thin-yet-powerful layer of functionality on top of OpenVDB."
  https://picogk.org/coding-for-engineers/8-first-steps-in-picogk.html (chapter undated; footer (c) 2024-2026)
- Design stance: "RISC" reduced instruction set, a move away from Hyperganic Core, which built latticing, infills and warping into the core. Lowest tier lets you "paint matter into space". Kayser, 2023-12-06.
  https://jlk.ae/2023/12/06/the-power-of-reduced-instruction-sets/
- Repo: created 2023-10-09, last push 2026-08-27, about 1.1k stars.
  Releases: v1.7.7.5 (2025-08-05), v2.0.0 (2026-04-21, large rewrite of memory management, new viewer, slice visualisation, overhang-angle visualisation "important for 3D printing", moved to OpenVDB transforms "for better compatibility with .VDB exchange"), v2.1.0 (2026-05-25, cylindrical/spherical coordinate types), v2.2.0 (2026-06-05, `Overhang` type, OpenVDB updated to v13), v2.3.0 (2026-08-05, `Rad` angle type).
  https://github.com/leap71/PicoGK/releases
- Kayser post (2026, undated in result) says he spent about 8 months rewriting core internals; PicoGK ships as a NuGet package; he frames LLMs writing engineering code as a goal. [author claim]
  https://x.com/linkayser/status/2046927799882321943
- Who uses it: LEAP 71 itself (all Noyron work is built on it); the RP/CEM propulsion model built on it is used by The Exploration Company, while RP/CEM stays proprietary.
  https://develop3d.com/cad/leap-71-picogk-open-source/ (2023) ; https://www.tctmagazine.com/leap-71-picogk-open-source-computational-engineering-framework/ (2023)
- Third-party: taucad/picovoxel, a WebAssembly + TypeScript port of PicoGK (repo created 2026-07-17, Apache-2.0, explicitly "not endorsed by, or supported by LEAP 71"). It exposes voxel CSG, offsets, fillets, shells, lattices, STL/GLB/VDB I/O, and a slicing module producing contours, SVG and CLI files.
  https://github.com/taucad/picovoxel

### Real parts made (LEAP 71 / Noyron, which is built on PicoGK)
Note: LEAP 71 press releases name Noyron, not PicoGK. The PicoGK link comes from LEAP 71 saying all its work is built on PicoGK (https://x.com/leap_71/status/2103485378577797570 ; Kayser post above).
- 2024 (June/July): TKL-5, 5 kN kerolox thruster, hot-fired twice (3.5 s and 12 s) with Airborne Engineering and University of Sheffield; reported undamaged.
  https://www.tctmagazine.com/leap-71-world-first-hot-fire-test-computationally-engineered-rocket-engine/
- 2024-11 (Formnext): 5 kN aerospike printed by Aconity3D on AconityMIDI+ LPBF in CuCrZr, 210 mm tall, under 62 h, 6.9 kg powder; depowdered by Solukon, heat treated by Fraunhofer ILT. Hot-fired 2024-12-18 at Westcott, UK.
  https://leap71.com/2024/11/21/leap-71-and-aconity3d-unveil-3d-printed-aerospike-rocket-engine-at-formnext-2024/ ; https://leap71.com/2024/12/23/leap-71-hot-fires-advanced-aerospike-rocket-engine-designed-by-computational-ai/ ; https://www.tctmagazine.com/leap71-successfully-hot-fire-tests-aerospike-engine/ (2025-01-02)
- 2024-11 (Formnext): with Eplus3D, "world's largest additively manufactured rocket thruster" [vendor claim], printed in one 354 h run on EP-M650H, about 43.5 kg, "no post-processing applied".
  https://www.tctmagazine.com/leap-71-presents-two-large-scale-single-piece-rocket-engines-computational-engineering/
- 2025-11-06: 600 mm part for XRB-2E6 methalox engine printed on Nikon SLM NXG 600E. [vendor claim: "one of the largest and most complex"]
  https://leap71.com/2025/11/06/nikon-slm-leap71/
- 2025-11-12 (Formnext): 1.5 m tall hypersonic precooler (heat exchanger) concept printed on Farsoon FS811M-U-8.
  https://leap71.com/2025/11/12/farsoon-leap71-formnext25/ ; https://3dprint.com/322203/pushing-hypersonic-limits-leap-71-farsoon-reveal-ai-designed-precooler-at-formnext/
- Summer 2025: multi-material heat exchanger, multi-material aerospike (~300 mm), and an electric motor concept (~300 mm) made with Fraunhofer IGCV's dual-metal PBF-LB process on a Nikon SLM machine (Metal AM magazine; page returned HTTP 403 to fetch, detail from search snippet only).
  https://www.metal-am.com/articles/leap-71-why-engineering-must-move-beyond-cad-to-realise-the-promise-of-ai-and-additive-manufacturing/
- 2025-12-11: two 20 kN methalox engines (bell and aerospike), CuCrZr, printed by Aconity3D, hot-fired; bell reached steady state at nominal chamber pressure with ">93% combustion efficiency"; aerospike reached 50 bar but completed only one burn. "Under three weeks from specification to first flame" and "designed entirely without human intervention" [vendor claims]. 200 kN and 2,000 kN designs in manufacturing validation, 2026 hot-fire target.
  https://leap71.com/2025/12/11/leap-71-hot-fires-two-orbital-class-methalox-engines-designed-autonomously-by-noyron/
- 2026-02-28: stackable liquid-liquid monolithic heat exchanger "printable on an 800x800mm metal 3D printer" (X post, design shown, print status not stated).
  https://x.com/leap_71/status/2027743142871519726

### Path to manufacturing and accuracy
- All shown parts are metal laser powder bed fusion (LPBF / PBF-LB). No machined-from-stock part found.
- Machining: for the 5 kN aerospike, "Except for the cutting of the threads, no post-machining is required." TCT Magazine (around late 2025).
  https://www.tctmagazine.com/leap-71-presents-two-large-scale-single-piece-rocket-engines-computational-engineering/
- Export formats in PicoGK source: `IO/MeshIo.cs` (STL; internal units mm; STL `UNIT=` header written, with the code comment that "many STL importers out there, ignore the UNIT info"), `IO/Cli.cs` (ASCII CLI Common Layer Interface writer `WriteSlicesToCliFile`, units configurable), `IO/OpenVdbFile.cs` (.vdb), `IO/ImageIo.cs`, `IO/VoxelsIo.cs`.
  https://github.com/leap71/PicoGK/tree/main/IO ; https://github.com/leap71/PicoGK/blob/main/PicoGK_MeshIo.cs
- ShapeKernel getting-started: by default PicoGK exports STLs; example changes voxel size "from 0.5mm to 0.1mm to increase the resolution".
  https://github.com/leap71/LEAP71_ShapeKernel/blob/main/Documentation/README-GettingStarted.md
- Voxel size is a single global parameter (`Library.Go(fVoxelSizeMM, ...)`, mm). Book: "Every voxel has a size, details below that size cannot be expressed"; "the size of a voxel is ultimately a measure for the precision of your output"; going from 0.2 to 0.4 mm gives "about a 6-8 times speed increase".
  https://picogk.org/coding-for-engineers/8-first-steps-in-picogk.html
- Third-party rule of thumb (Tau): voxel size about one third of thinnest wall; time scales 1/h^3, memory and triangles 1/h^2.
  https://github.com/taucad/tau/pull/764
- Voxel size actually used for the flight-relevant engines: NOT FOUND.

### Stated limits
- No STEP or B-rep export found in the PicoGK repo IO folder or docs [inferred from file list].
- No statement found from LEAP 71 on exact dimensions, true cylinders, fillet radius, drawings, or tolerances.
- PyVCAD paper (2025) says PicoGK "does not support multi-material heterogeneous design".
  https://matterassembly.org/assets/pdf/preprints/2025_Wade_pyvcad.pdf
- Kayser's broader position: "CAD and the visual process that predates it had a good run" [author claim].
  https://x.com/linkayser/status/2046927799882321943

---

## 2. Hyperganic Core

- Proprietary voxel geometry kernel for additive manufacturing (Munich, founded 2015, CEO Lin Kayser at the time). "Much of the platform is written directly in Assembly"; GPU rendering of "every particle" [vendor claims].
  https://3dprintingindustry.com/news/hyperganic-launches-its-ai-based-algorithmic-design-software-to-the-public-210883/ (2022)
- Funding: $7.8M first VC round, 2021-02-17, alongside Core 2.0 and Print Framework 2.0.
  https://3dprint.com/280450/hyperganic-raises-7-8m-in-first-venture-round-for-ai-driven-3d-printing/
- Core 3 public release June 2022.
  https://www.tctmagazine.com/hyperganic-public-rollout-ai-powered-design-software/
- Acquired DirectFEM (TUM spinout) to add simulation.
  https://3dadept.com/hyperganic-to-create-a-one-stop-solution-from-engineering-to-manufacturing-with-the-acquisition-of-directfem/
- Parts: 2022 EOS partnership, full aerospike engine printed on EOS LPBF; earlier AI-designed combustion chamber.
  https://www.3dnatives.com/en/eos-and-hyperganic-3d-print-more-efficient-space-propulsion-components-120520224/ ; https://www.additivemanufacturing.media/news/eos-hyperganic-collaborate-on-design-performance-of-space-propulsion-components
- 2024: launched HyDesign, a cloud lattice design app (X post 2024-07-24); joined IAM3DHUB March 2024.
  https://x.com/hyperganic ; https://www.iam3dhub.org/hyperganic-new-technological-partner-on-board-at-iam3dhub/
- Status: PitchBook lists "Out of Business", latest deal type "Liquidation". Singapore entity Hyperganic Pte. Ltd. "Dissolved - Creditors Voluntary Winding Up" (last annual return 2024-09-19). Crunchbase still says active (likely stale).
  https://pitchbook.com/profiles/company/222979-51 ; https://recordowl.com/company/hyperganic-pte-ltd ; https://www.crunchbase.com/organization/hyperganic
- Output formats and accuracy: NOT FOUND.

---

## 3. nTop implicit side (pure-implicit workflow only)

- Implicit bodies are native; implicit data moves in and out only as nTop's own `.implicit` file; to export otherwise "you must first convert the implicit body into a mesh" (support doc, lists nTop 5.45).
  https://support.ntop.com/hc/en-us/articles/360048249514-What-file-types-can-nTop-import-and-export
- Direct slicing without a mesh: Slice Body block (layer height, feature size, frame), marching-squares variant, Export Slices; black-and-white image-stack slicing for DLP and binder jet; build-volume blocks for EOS M and Renishaw.
  https://support.ntop.com/hc/en-us/articles/360042790233-How-to-create-and-export-slices ; https://support.ntop.com/hc/en-us/articles/1500007187802-How-to-Use-B-W-Image-Slicing
- 3MF export carries outer mesh plus beam lattice (surface lattices not supported in 3MF route).
  https://support.ntop.com/hc/en-us/articles/360054069993-How-to-export-a-mesh-and-lattice-to-3MF-format
- "Implicit Interop": around Formnext 2022, Siemens Energy heat exchanger sliced directly from an nTop implicit in EOSPRINT because the STL or CLI would have been too large; nTop file "<1Mb"; "you can even get more accuracy than when slicing an already discritized mesh" [vendor claim, CEO Bradley Rothenberg].
  https://cy.linkedin.com/posts/bradleyrothenberg_implicitinterop-implicitinterop-activity-6997282905710329856-TFfj ; https://www.linkedin.com/posts/ntopology_implicit-interop-in-ntop-export-without-activity-7079173770141208577-dHJM (2023)
- Limits stated by nTop: getting STEP/Parasolid requires conversion (quad-mesh route, "CAD Body from Implicit Body" within a tolerance, or "Simplified Body Representation" for million-face lattices); sharp edges round unless refined or sharpened; nTop says users do early design in nTop then export to their CAD/CAE/PLM tools "for detailed design, final validation, and documentation".
  https://support.ntop.com/hc/en-us/articles/7323484481555-How-to-export-geometry-as-a-CAD-part-STEP-or-Parasolid ; https://support.ntop.com/hc/en-us/articles/7323109441683-nTopology-3-28-What-s-New ; https://www.ntop.com/platform/
- (B-rep side covered by another helper.)

---

## 4. libfive (Matt Keeter)

- F-rep (functional representation) solid modelling framework. C++ core with C API; Python and Guile Scheme bindings; Studio GUI (OpenSCAD-like live coding plus direct push/pull editing that changes script variables).
  https://github.com/libfive/libfive
- Licence: libfive, stdlib and Python bindings MPL-2.0; libfive-guile and Studio GPL-2.0-or-later.
- Meshing: "watertight, manifold, hierarchical, feature-preserving triangle meshes" [author claim]. Output is triangle mesh (STL).
- Repo: created 2015-12-12, last push 2026-10-16, ~1.7k stars, no tagged releases.
- Engineering parts made with it: NOT FOUND.
- STEP, drawings, tolerances: not mentioned in README.

## 5. Fidget and Halfspace (Matt Keeter, Rust)

- Fidget: "blazing fast implicit surface evaluation"; tape simplification from interval evaluation; hand-written aarch64/x86_64 JIT; meshing by its own Manifold Dual Contouring implementation; canonical bytecode format; wasm32 support (interpreter only).
  https://github.com/mkeeter/fidget
- Licence MPL-2.0; created 2022-06-24; latest release v0.5.1 (2026-09-27), v0.5.0 (2026-08-03); ~500 stars.
- Author-stated limits: "strong Lego-kit-without-a-manual energy"; "if you're doing serious work with the library, expect to fork it"; "personal-scale experimental project"; "its meshing is much less battle-tested than libfive"; only a debug GUI.
- Halfspace: "experimental IDE for doing solid modeling with distance fields", Rust, MPL-2.0, created 2025-04-06, last push 2026-10-19; scripted via fidget-rhai; web and native.
  https://github.com/mkeeter/halfspace ; https://mattkeeter.com/projects/halfspace
- Engineering parts: NOT FOUND.

## 6. ImplicitCAD (Haskell)

- Programmatic CAD in Haskell with an OpenSCAD-compatible parser (extopenscad); CSG plus rounded unions and rounded extrusions; "directly provides some GCode generation" (2D).
  https://github.com/Haskell-Things/ImplicitCAD
- Licence AGPL-3.0; repo created 2011-10-21; last push 2026-04-13; latest release v0.4.1.0 (2023-12-18).
- Outputs seen in README: SVG, PNG, 3D meshes (STL), 2D gcode.
- Engineering parts: NOT FOUND. STEP, tolerances: not mentioned.

## 7. Curv

- Language "for making art using mathematics" (repo description). Programs compile to GPU fragment shaders. Exports STL, OBJ and X3D meshes "for 3D printing"; X3D carries full colour; meshes claimed "defect free: watertight, manifold" [author claim].
  https://github.com/curv3d/curv
- Licence Apache-2.0; created 2017-05-15; last push 2023-09-27 (no activity in about three years).
- Engineering use: none found; aim stated as art.

## 8. SDFX (Go)

- "A simple CAD package written in Go"; 2D and 3D SDFs defined in Go code; rendered to STL/3MF "to be viewed and/or 3d printed".
  https://github.com/deadsy/sdfx
- Licence MIT; created 2017-01-09; last push 2026-08-29; ~630 stars; no tags.
- Claim: "SDFs can easily do filleting and chamfering (hard to do with OpenSCAD)" [author claim]. Hackaday echo, 2022-11-28: https://hackaday.com/2022/11/28/taking-distance-based-cad-to-the-next-level/
- Gallery of mechanical-type parts: pottery wheel casting pattern and core box, cylinder head, M-square casting pattern, involute gear, camshaft, Geneva mechanism, nuts and bolts. Whether these were cast, printed or machined: not stated.
- Roadmap lists known gaps: fix distance function for non-linear extrusions "E.g. screws"; generate smaller/better STL files.
  https://github.com/deadsy/sdfx/blob/master/docs/ROADMAP.md
- Related Go: soypat/gsdf (GPU-accelerated, BSD-3-Clause, created 2024-08-11). https://github.com/soypat/gsdf
- Related Python: fogleman/sdf (MIT, 2021) https://github.com/fogleman/sdf ; kylebeck/sdf fork (GPU via PyTorch/Taichi, dual contouring with QEF for sharp edges, aims at "complex mechanical assemblies" [author claim]) https://github.com/kylebeck/sdf

## 9. OpenVCAD / pyvcad (MacCurdy lab, CU Boulder)

- Implicit multi-material modelling: at any XYZ it returns signed distance (geometry) and volume fractions (material). Scales to "hundreds of billions of voxels" [author claim].
  https://matterassembly.org/MacCurdyLab.github.io/openvcad
- Paper: Wade et al., "OpenVCAD: An open source volumetric multi-material geometry compiler", Additive Manufacturing 79 (2024) 103912, doi:10.1016/j.addma.2023.103912; includes an open-source inkjet slicer.
  https://www.sciencedirect.com/science/article/abs/pii/S2214860423005250
- pyvcad paper: arXiv 2509.15562 (2025-09-19), ACM SCF 2025. Python API; outputs PNG stacks (inkjet), surface meshes for FFF/standard slicers, FEA meshes with adaptive sizing. Example: multi-material bicycle seat, graded lattices.
  https://arxiv.org/abs/2509.15562 ; https://pypi.org/project/OpenVCAD/
- Gradient-informed slicing for FFF toolpaths: arXiv 2505.08093 / doi:10.1016/j.addma.2025.104963.
  https://arxiv.org/pdf/2505.08093
- The pyvcad paper says OpenVDB and NanoVDB "are primarily tailored for rendering and lack specialized functionalities needed to convert sparse voxel structures into manufacturable formats suitable for additive manufacturing".
  https://matterassembly.org/assets/pdf/preprints/2025_Wade_pyvcad.pdf
- Licence: NOT CONFIRMED (PyPI package exists; lab's own repo not found; a 2026 fork has NOASSERTION).
- Printer voxel resolution used: NOT FOUND in abstract.

## 10. ShapeJS / AbFab3D (Shapeways)

- JavaScript scripts over the AbFab3D voxel library (Java). Each voxel stores a normalized, limited-range signed distance, 8-bit by default (256 levels). Default voxel size 0.1 mm. Exports STL, X3D (colour), experimental SVX voxel format. Back end in CUDA on Shapeways GPUs (2020).
  https://shapejs.shapeways.com/v2/docs ; https://en.wikipedia.org/wiki/ShapeJS ; https://stackoverflow.blog/2020/02/06/how-shapeways-software-enables-3d-printing-at-scale/
- AbFab3D repo: last push 2022-09-23. https://github.com/AbFab3D/AbFab3D
- Use was consumer and jewellery-style printed goods; one user reports a 3-minute model timeout. https://www.instructables.com/Ember-Printer-Creating-and-Printing-Models-Using-S/
- Shapeways filed Chapter 7 on 2024-07-02; brand relaunched 2024-12-03 without marketplace and shops. ShapeJS status after that: NOT FOUND.
  https://hackaday.com/2024/07/04/shapeways-files-for-bankruptcy/ ; https://en.wikipedia.org/wiki/Shapeways

## 11. OpenVDB / NanoVDB as a modelling kernel

- OpenVDB: Apache-2.0, ASWF, created 2013, active (push 2026-10-21). https://github.com/AcademySoftwareFoundation/openvdb
- Engineering modelling uses found:
  - PicoGK (above) is the clearest case of a full kernel on OpenVDB.
  - Dendro: volumetric modelling plug-in for Grasshopper built on OpenVDB (MPL-2.0, created 2018, push 2026-08-01). Used in a J. Computational Design and Engineering paper to model heterogeneous lattices; compared in a 2024 Applied Sciences review of lattice CAD tools.
    https://github.com/ryein/dendro ; https://academic.oup.com/jcde/article/9/5/1565/6656371 ; https://doi.org/10.3390/app14083306
  - Autodesk patent US 11947334 uses the OpenVDB method to extract the final level set in lattice topology optimisation.
    https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11947334
- Counterpoint: OpenVCAD authors (2025) say OpenVDB/NanoVDB are rendering-focused and lack AM output conversion (quoted above).
- NanoVDB used as an engineering modelling kernel: NOT FOUND.

## 12. 2024-2026 startups and research projects (SDF-first CAD for mechanical work)

- Metafold 3D (Toronto, founded 2020): cloud, API-first implicit geometry kernel; shape functions give signed distance; Python and JS SDKs, REST API; focus on lattices and metamaterials for AM; seed $1.78M (around 2023). Docs admit "implicit shapes often fail to capture features with sharp edges". Newer positioning is "Geometric Intelligence" (CAD feature detection for quoting).
  https://docs.metafold3d.com/core-concepts ; https://www.metafold3d.com/ ; https://develop3d.com/cad/metafold-secures-funding-geometry-engine/
- Gradient Control Laboratories (Blake Courter, Luke Church): research on Unit Gradient Fields (UGFs) as a better-conditioned alternative to SDFs for engineering; announced Omega, "the world's first pure geometry as code development system" [author claim], 2026-02-07, self-described "at its earliest stages". nTop funds GCL and keeps access to its research. Author: "We are aware that we have little evidence for our extraordinary claims."
  https://www.gradientcontrol.com/ ; https://www.blakecourter.com/2026/02/07/geometry-as-code.html ; https://www.blakecourter.com/2023/06/03/foreword.html
- Tau (taucad): open-source "AI-native CAD platform", browser-based, multi-kernel (Replicad, Manifold, Zoo KCL, JSCAD, OpenCascade, build123d, PicoGK, PicoVoxel). Repo created 2025-01-04, Apache-2.0. PicoVoxel port created 2026-07-17. Not SDF-only.
  https://github.com/taucad/tau ; https://github.com/taucad/picovoxel
- Matt Keeter's Halfspace (2025) and Fidget v0.5 (2026): research-grade, see section 5.
- Aurorin CAD (YC W2026): new kernel, but B-rep, not SDF. Listed here only to record that it is not a match.
  https://www.startuphub.ai/ai-news/claudes-corner/2026/claudes-corner-aurorin-cad-yc-w2026
- Academic B-rep/implicit bridging (not kernels for production): NH-Rep (2022) https://arxiv.org/pdf/2209.10191 ; Patch-Grid https://arxiv.org/pdf/2308.13934 ; DualBrep (2026, SDF+UDF, neural rebuild can leave "small non-watertight artifacts") https://arxiv.org/pdf/2606.31579 ; "RBF Your SDF" (2026) https://arxiv.org/pdf/2609.06209
- Voxel resolution guidance for machining features (ML context): at least 4 voxels to detect a hole, 16+ to represent its size and curvature. arXiv 2601.06334 (2026). https://arxiv.org/pdf/2601.06334
- Funded 2024-2026 startup that sells a pure SDF CAD kernel for machined mechanical parts: NOT FOUND after four searches.

---

## Cross-cutting: what the sources say about limits

| Limit | What sources say |
|---|---|
| Exact dimensions | PicoGK: precision bounded by voxel size; no tolerance figures published. nTop: conversion to CAD body is within a tolerance. No pure-SDF tool found publishing dimensional accuracy figures. |
| True circles/cylinders | No tool found that keeps an analytic cylinder through export; all export meshes or slices. |
| Fillet of a given radius | SDFX, ImplicitCAD, picovoxel advertise fillets/rounded unions; no source verifies a measured radius after meshing. Courter: blends "always succeed" [author claim]. |
| Sharp edges | Metafold docs: implicits "often fail to capture features with sharp edges". nTop needs sharpening/refinement. libfive and kylebeck/sdf claim feature-preserving or QEF dual contouring. |
| STEP export | None of the pure tools found export STEP natively. nTop does via approximate conversion. Tau converts formats but through other kernels. |
| Drawings | No pure-SDF tool found producing 2D engineering drawings. nTop says documentation happens downstream. |
| CNC machining | No pure-SDF tool found generating milling toolpaths. ImplicitCAD has 2D gcode. LEAP 71 parts: only thread cutting post-print on the 5 kN aerospike. |

## Gaps
- Voxel size LEAP 71 used for its hot-fired engines.
- Whether LEAP 71 machines sealing faces or flanges on other engines.
- Hyperganic Core output formats and accuracy.
- Any engineering part made with libfive, Fidget, ImplicitCAD or Curv.
- OpenVCAD licence and inkjet voxel resolution.
- ShapeJS status after 2024.
- Metal AM "beyond CAD" article (HTTP 403); details taken from search snippets only.
