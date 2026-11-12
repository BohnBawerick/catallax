## 2026-11-12: the field builds imported spheres and tori with trims, as bodies of revolution

Pieces: field, compare

Change: `Plan::compile_with_files` builds an imported STEP part with trimmed spherical or toroidal faces when the part is a body of revolution: every face is a plane across one common axis or a cylinder, cone, sphere or ring torus about it, and every boundary is a whole parallel circle, or a pole or apex on the axis. The part is read back as its meridian profile and kept as the field's full revolution of that profile, with each face under the name or names of the imported face. Hemispheres, domed cylinders, lenses, spherical zones, torus halves and caps, and cylinders with a rounded edge build this way. The axis may point anywhere, provided its direction and a perpendicular have rational unit vectors. `field::mass::volume` gives the exact volume of such a part, as `pi` times a rational number plus tube-turn angles for torus arcs.

Why: M21 asks for imported parts made of planes, cylinders, cones, spheres and tori in the field. Trimmed spheres and tori were the remaining analytic gap and returned `field.import-domain`.

Details:

- Orientation is checked, not assumed. Each face's outward normal and the way its boundary circle runs fix the direction of the profile with material on its left. Every circle must be left by one face and entered by one, the axis must be entered at the top of each interval and left at the bottom, and each profile loop must run counter-clockwise at even nesting depth and clockwise at odd depth. A reversed face or an inward-facing shell returns `field.import-domain` and is never published as a filled solid.
- Still refused, with `field.import-domain` and the reason: a sphere or torus cut about two axes, a boundary split into arcs, a boundary that is not a whole circle, an irrational axis, circle radius or placement, two faces meeting at one point of the axis, a horn or spindle torus, and a complete sphere or torus beside other faces. Spline faces keep `field.import-spline`.
- Cylinders, cones and planes alone keep their existing round and planar cells; only a part with a sphere or torus face takes the new route. Complete spheres and tori keep their primitives.
- The revolved faces keep the roles they were built with for sampling and trimming, so candidate sampling and face queries work under the imported names.
- M16's distance, clearance and measurement queries are not extended; they answer as they do for any revolved sweep.

No existing answer changes. Every part these rules admit was refused before. `field::mass::volume` keeps its results for cells and primitives; its precision loop is shared. No tolerance, saved replay answer or work allowance changes. Raw evidence retained privately.

Migration: Callers that expected `field.import-domain` for a hemisphere, half torus or other body of revolution now receive a built field. The refusal text for the remaining cases now begins with the face count and `faces with spherical or toroidal trims`, followed by the reason.
