# Field STEP import for the current analytic subset

`Plan::compile_with_files` and `Plan::compile_with_files_checked` accept the existing immutable `ImportFiles` map, keyed by the SHA-256 in each recipe file reference. Ordinary compile methods pass an empty map. An import without its bytes now returns `field.import-missing-source`.

STEP admission, saved repair decisions and rigid placement use the same import route as the B-rep. The field retains finite analytic definitions after admission and discards the import model. Point membership, distance, normals, named-face trimming, sampling and display run on those retained definitions. This shares source admission; it does not use B-rep point containment or B-rep distances to answer field queries.

The admitted field domains are closed polygonal planar parts, including nonconvex parts, and one outward cylindrical or 45-degree conical wall between two perpendicular circular disks. The round wall and both disks must share their complete rims. Other curved trims return `field.import-domain`, rather than replacing the part by supporting halfspaces or a convex hull. Face origins retain the file's `ADVANCED_FACE` entity keys, including every origin of a merged planar face.

Planar source edges add arrangement delimiters with no material contribution. When coplanar faces merge geometrically, named queries retain the original closed source polygons. A name cannot expand across an adjacent face's trim. These polygons follow subsequent rigid transforms with the field.

Planar extraction also refuses with `field.import-domain` if it produces no positive material component or loses any admitted face name. An inward shell cannot silently become a proven-empty imported part. The refusal reports the attained and expected name counts.

Spline faces return `field.import-spline` and identify M10. Sphere, torus and other unsupported source definitions return `field.import-unsupported`. The shared importer still limits cones to its existing exact 45-degree degree-unit representation. This change does not extend shared STEP admission to general cones, spherical faces, toroidal faces or `VERTEX_LOOP`.

`mass::volume` integrates one retained planar or axial field cell. It returns a certified enclosure under the existing quantity Accuracy and precision/work limits. Other field definitions return `field.volume-definition`; their existing `mass::properties` query remains available.

`Revision::check_with_files` requests the existing two-way boundary comparison using the same source bytes in both kernels. Imported parts additionally compare 125 exact containment samples on a fixed grid extending outside the part's bounds. They also compare independently integrated volume enclosures. `Options::import_volume_accuracy` defaults to the existing quantity Accuracy: relative `1e-12`, absolute `1e-30` cubic millimetres. Each integration receives one quarter of that width. Agreement requires the worst difference over both complete enclosures to fit the requested bound. Unresolved classification, volume or boundary checks remain unresolved. A membership or volume disagreement cannot produce `Agree`. Samples do not certify agreement at every point.

The verification CLI accepts `catallax-import-recipe/1` wrappers for the recipe, field and compare sources. A wrapper contains exactly `format`, `recipe` and `files`. `files` maps SHA-256 strings to complete UTF-8 STEP text. Every supplied hash is checked. File names remain labels; the CLI opens no paths from the recipe. Portable question replay retains the source bytes and original work allowance.

For import wrappers, comparison manifests add `comparison.import_agreement`, with exact sample positions and classifications, both volume enclosures, and the attained and allowed volume difference. Kernel refusals retain their code and detail. Ordinary recipe manifests keep their existing shape.

No modelling tolerance, saved replay answer or work allowance changes. Raw evidence retained privately.
