## 2026-10-26: qualified B-rep print output

Pieces: mesh consumer, features, deterministic verification interface and verification driver.

Change: `source: print` selects the B-rep route for ordinary recipes and the field route for result dependencies containing lattice, gyroid, smooth-blend or mesh-import definitions. It never retries another route after a refusal. `source: field` with a print wrapper explicitly retains the field route. The default accuracy remains 0.01 mm, the cell cap 2,000,000 and the vertex cap 1,000,000. No modelling tolerance or work allowance changes.

Why: ordinary finite faces and trims need certified tessellation, not uniform field cells. A 4 x 4 x 1 mm box needs twelve planar facets.

Migration: `catallax-print/1` accepts optional `route: "brep" | "field"` and `certified: true`. The latter opts into `evaluate_certified` and consumes its geometry-bound uniform model certificate. It cannot be combined with the field route. Ordinary unmarked exact results have zero model error. Marked geometry without a complete certificate refuses, rather than treating a mark as a whole-result bound. Consumers that require native-field qualification must request `route: "field"` or use `source: field`. Existing field print construction and its separate acceptance gate are unchanged.

### Qualification

The B-rep adapter calls `mesh_body` with shared physical-edge vertices. It sums the maximum certified face bound and maximum finite-trim edge bound, conservatively including endpoint replacement and binary64 rounding. The mesher receives half the remaining requested accuracy for each bound; final composition must still leave room for encoding. Each binary32 vertex displacement has an exact rational L1 bound, which dominates Euclidean distance. Barycentric interpolation extends that bound in both directions over every corresponding facet. The complete certificate is `W_file <= U_model + V_tessellation + V_encoding <= A`.

Before publication, the encoded mesh must have oppositely paired edges and one connected cyclic link per vertex. The adapter checks the entire affine vertex motion from binary64 to binary32. A bounding-volume tree proposes pairs. Exact rational Bernstein coefficients certify noncollapse, retained facet orientation and separating planes over the whole motion interval. Each facet pair must remain disjoint except for its common indexed simplex. This excludes self-intersections and merged sheets, including contacts between components or cavities. A missing separation proof is a refusal, not a sampled success. Producer material senses are retained, with positive signed volume for each outer shell and negative volume for each cavity. Exact ray parity checks cavity ownership and excludes nested overlapping material components. No shell is flipped independently.

The adapter does not accept a preview or cached display as print evidence. It does not close the native-field print gate or implement adaptive field cells and trade-off estimates.

### File and manifest contract

Every print result has a `printing` record naming its route, requested `tolerance_mm`, effective `max_cells` and `max_vertices`, `cell_cap_applicable`, units and coordinate encoding. B-rep `cells` is zero and its cell cap is not applicable. Refused outputs have no achieved `deviation_bound_mm` or file. Successful outputs record model, tessellation and encoding contributions, the achieved bound, counts and separate completeness, topology and orientation statuses. `files` records SHA-256 and byte length for both artifacts; `printing.sha256` is the STL digest. Verify the digest against actual file bytes before trusting the accompanying certificate. Hashes bind metadata and geometry together, not the authority of an untrusted producer.

The binary STL header is `catallax-print-stl/1`. All multibyte fields are little endian. Offsets are zero-based; reserved bytes must be zero.

| Bytes | Meaning |
| --- | --- |
| 0..8 | ASCII `CXPRINT1` |
| 8 | Route, 1 B-rep or 2 field |
| 9 | Units, 1 millimetres |
| 10 | Coordinates, 1 IEEE-754 binary32 little endian |
| 11 | Cell-cap applicability, 0 not applicable or 1 applicable |
| 12..20 | Exact requested accuracy numerator, signed 64-bit integer |
| 20..28 | Exact requested accuracy denominator, unsigned 64-bit positive integer |
| 28..36 | Configured cell cap, unsigned 64-bit integer |
| 36..44 | Configured vertex cap, unsigned 64-bit integer |
| 44..52 | Model bound, finite nonnegative binary64 |
| 52..60 | Tessellation bound, finite nonnegative binary64 |
| 60..68 | Coordinate-rounding bound, finite nonnegative binary64 |
| 68..76 | Achieved whole-file bound, finite nonnegative binary64 |
| 76..80 | Reserved zero bytes |

The requested rational is stored losslessly. Bounds are rounded upward to exactly represented dyadic values, then summed and rounded upward again. The manifest uses the exact rational values of these stored binary64 bounds. The header's achieved bound must dominate the sum and must not exceed the exact requested rational. Accuracy values outside the integer encoding range refuse with `print.metadata-range`; they are never truncated. `print.metadata-accuracy` refuses if conservative metadata rounding exhausts the accuracy. OBJ's leading `# catallax-print/1` JSON comment carries the same full settings and certificate. Its decimal coordinates round-trip to the binary32 values stored in STL.

STL itself standardizes neither units nor accuracy metadata. A generic reader may ignore this header and loses the guarantee. A conforming consumer must recognize this version, reject unknown units, encoding, route, flags or nonzero reserved bytes, validate the rational and bound arithmetic, and preserve the header or use the accompanying metadata-bearing OBJ. Facet winding is authoritative; advisory STL normals do not replace topology and orientation checks. Deterministic comparisons include the new header bytes, manifest, SHA-256 values and OBJ comments.

Specific new refusals are `cli.print-route`, `print.invalid-options`, `print.brep-field-definition`, `print.brep-model-bound`, `print.brep-topology`, `print.brep-accuracy`, `print.brep-encoding-range`, `print.brep-encoding-topology`, `print.brep-orientation`, `print.metadata-range` and `print.metadata-accuracy`. Build and mesher refusals retain their original typed codes. No refusal publishes partial geometry or silently loosens controls.
