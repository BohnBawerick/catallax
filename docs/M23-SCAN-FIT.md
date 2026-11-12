# M23: classic fitting of scans to a clean part

`catallax_field::reconstruction::scan::fit_scan` takes a `scan_import` step (RECIPE.md section 4.19.1) and its file bytes, and returns candidate surfaces, typed refusals for whatever it cannot explain or place, and, when the block family below fits, a clean part as an ordinary recipe built in both kernels with both deviation numbers. No recipe or feature tree is input. Learned segmentation is not used. Free-form regions and the written verdict on learned segmentation such as DualBrep wait for M10 and M11.

## What is proved and what is not

Every fitted surface is a candidate: a binary64 least-squares fit to the samples it was given. The part is the recipe the fit writes, and only that recipe's values are exact: the B-rep and the field evaluate it under their existing definitions and proofs, and the B-rep passes structural validation before it is returned. Nothing claims that the part equals the scanned object. Both deviation numbers are unproved binary64 maxima over stated samples, not Hausdorff bounds. The only proven quantity in the deviation is the display mesh's chord bound from `catallax-mesh`.

## Method

1. **Samples.** The scan's placed points: triangle corners, then isolated samples, in the reference's order. Between 100 and 1,000,000 samples; outside that the fit refuses with `reconstruction.too-few-samples` or `reconstruction.too-many-samples`.
2. **Neighbourhoods.** A k-d tree with index tie-breaks gives each sample its 16 nearest samples. Their covariance gives a normal, a flatness score and a local residual. The noise deviation is the median local residual, corrected for the three fitted parameters, unless the caller states it. The 10 nearest neighbours, symmetrised, form the region graph.
3. **Region growing.** From the flattest unassigned sample, a plane grows over graph neighbours within 4 noise deviations of it whose normals are within 30 degrees, refitted each time the region doubles. Every sample ends in some region; curved surfaces end as many narrow strips.
4. **Merging.** Regions are merged, largest first, with each neighbour in order of shared graph edges, when one surface explains both: the simplest of plane, cylinder and cone whose root-mean-square residual on each side is within 1.5 noise deviations, with at most one sample in a hundred beyond the growth bound. A cylinder or cone must bend visibly over its samples (its sagitta across their angular span exceeds the growth bound), and a cone must taper visibly (its radius changes by more than the growth bound along them), so noise never upgrades a plane to a cylinder or a cylinder to a cone. A region less than a twentieth the size of its neighbour joins only on the neighbour's own surface. Failed pairs are retried only after one side changes kind or doubles.
5. **Assignment.** Regions of at least 30 samples keep their surfaces. Surfaces then grow over unassigned neighbours within the growth bound, boundary samples move to the surface with the smallest residual, surfaces are refitted, and surfaces that one surface explains are joined, connected or not; a smaller surface whose samples a larger one explains nine in ten gives those samples up. A refitted cone that no longer tapers visibly becomes a cylinder. One recovery pass then fits each connected group of 30 or more unexplained samples as a whole, which catches a narrow round whose growth strips were too thin to merge. Every round updates all samples from the previous round's state.
6. **Assembly.** Plane normals within 1 degree are grouped; the three mutually perpendicular groups with the most samples set the stock frame. A plane across a cylinder or cone axis whose samples lie within its radius is that surface's cap. The outermost other planes of each group are the stock faces. A cylinder or cone parallel to a stock axis whose samples go more than half way round it is a boss when its samples lie beyond a stock face and it has a cap, and a bore when they lie within the stock, through or blind. A cylinder parallel to a stock edge whose axis lies one radius inside both faces, within the larger of 4 noise deviations and a tenth of the radius, and whose samples face the edge, is a convex round; its radius is then refitted with the axis held tangent to both faces, as the recipe's fillet defines it. Rounds whose radii agree within the larger of 4 noise deviations and 2 per cent share one fillet step. Unplaced surfaces smaller than one per cent of the samples are released and the samples assigned again before the final assembly.
7. **Recipe.** `catallax-recipe/4`: the scan step, a block, each boss united, each bore subtracted, then the fillets. Lengths are rounded to 1 µm and quaternion components to 10⁻⁶; tool placements are exact rational offsets from the stock's corner, so a boss stands exactly on its face. Through cuts run 1 mm past the faces they open. A kernel refusal is returned as the part's refusal, with the attempted recipe kept.
8. **Deviation.** The B-rep is meshed with a proven chord bound of 0.001 mm. Scan to part: every sample, plus interior points of scan triangles at four barycentric subdivisions, to that mesh. Part to scan: the mesh's corners, its edges at 0.25 mm and its interior at about one point per 0.0625 mm², on an additive recurrence, to the scan reference. Each side reports its sample count, maximum, root-mean-square and 95th percentile; the maps keep every value, and each part sample's B-rep face.

Normal estimation runs in bounded parallel batches and writes results in input order. Searches, fitting, assignment and deviation queries charge the caller in input order and check its deadline. Answers and work counts are identical at 1, 2 and 4 workers. Direct calls use the decider policy's work cap; an enclosing operation allowance remains authoritative.

## Refusals

- `reconstruction.unexplained-region`: samples that no plane, cylinder or cone explains, with their total, the number of connected regions and the largest. Each region's samples and centroid are returned.
- `reconstruction.unplaced-surface`: a fitted surface the block family cannot place, with its sample count and the reason.
- `reconstruction.no-stock`: no three perpendicular plane groups with two faces each, or the admitted groups do not fit the orthogonalized frame.
- `reconstruction.deviation-too-large`: either direction would retain more than 1,000,000 samples, the mesh index would exceed that size, or sampling counts or spacing are nonfinite. Counts and work are checked before allocating sample vectors.
- `reconstruction.nonfinite-sample`: a placed scan point is outside the finite binary64 fitting range.
- `numbers.operation-budget` and `numbers.deadline`: the caller's work or time limit stopped fitting.
- `reconstruction.recipe`, `reconstruction.invalid-part`, `reconstruction.nonfinite-deviation`, or the kernel's own refusal when the written recipe does not build.

## Seeded controls

Point clouds of isolated OBJ vertices, written to 1 nm by [the sampler](../crates/catallax-field/examples/scan_fit/sampler.rs) from each part's analytic patches, uniform by area, moved along the true normal by noise of deviation 0.01 mm (the sum of twelve uniforms from a SplitMix64 stream). Numbers in mm, from `cargo run --release -p catallax-field --example scan_fit -- <part> <count> 0.01 <seed>` on one 24-thread machine.

| Control | Samples, seed | Surfaces | Part | Scan to part max / RMS | Part to scan max / RMS |
| --- | ---: | --- | --- | ---: | ---: |
| 30 × 20 × 8 plate, bore r 4 | 20,000, 1 | 6 planes, 1 cylinder | block 30.002 × 20 × 8, bore r 4 | 0.0375 / 0.0103 | 0.578 / 0.184 |
| same, cone bore r 3 to 5 | 20,000, 2 | 6 planes, 1 cone (slope 0.2497) | block 29.999 × 20 × 8, cone tool r 2.751 to 5.248 over 10 | 0.0396 / 0.0100 | 0.552 / 0.184 |
| 40 × 30 × 10 block, boss r 8, four top edges rounded r 1 | 30,000, 3 | 7 planes, 5 cylinders | block 40.002 × 29.999 × 10, boss r 8 h 11.999, one fillet r 1.008 on four edges | 0.0429 / 0.0103 | 0.727 / 0.218 |

These timings were recorded before the work-accounting fixes. Each fit took about 1 to 3 s. Scan to part stays within about four noise deviations. Part to scan is dominated by the sampling: a part point is as far from the scan as the nearest sample, and the largest empty patches in 20,000 to 30,000 random samples are about 0.5 to 0.8 mm across. A hole in the scan shows in that direction only: removing a 3 mm-radius patch of the plate's top raises part to scan above 2.5 mm while scan to part is unchanged. A spherical dome on the plate is refused as an unexplained region with its count. Rounds fitted freely over a quarter arc come out 1 to 3 per cent small at this density; held tangent to their faces they come out within 1.5 per cent. Seeds 1 to 6 of each control give the same surfaces and the same family of part. On the radius-1 rounds at about 7 samples per mm², noise of 0.005 and 0.01 mm keeps all four rounds; at 0.02 mm one round is lost to the unexplained count, and at 0.03 mm the rounds are refused.

## Limits

The family is a block stock with bosses, bores and convex rounds on stock edges. Rounds where a boss meets the stock (tori), chamfers, pockets, steps, and features off the stock axes are refused as unexplained or unplaced. Rounds of different radii that meet at a corner become separate fillet steps, which the kernels may refuse. A bore or boss scanned less than half way round is refused as unplaced. Small rounds need noise well below their radius at the scan's density. Free-form regions need M10 and M11. The deviation numbers are sampled and unproved.

## The P3 round trip

PLAN M23 asks that 30,000 noisy points taken from P3 come back as P3's planes and cylinders, with a deviation map. [The sampler](../crates/catallax-field/examples/scan_fit/sampler.rs) writes P3 (`parts/P3-filleted-boss.json`) patch by patch: the six block faces, the top face outside the junction, the boss wall and top, the four radius-1 top-edge rounds mitred at the corners, and the radius-2 concave junction round, a quarter torus. A regression checks that 5,000 noise-free samples lie within 0.0001 mm of P3's own B-rep, meshed by `catallax-mesh` with a proven chord bound of 0.0001 mm. The fitter never sees P3's recipe.

Seed 3, 30,000 samples, noise 0.01 mm (estimated 0.0108 mm). The historical runtime before the work-accounting fixes was 2.2 s on 24 threads:

| P3 surface | Truth | Fitted |
| --- | --- | --- |
| Block faces x, y, z | 0 and 40, 0 and 30, 0 and 10 | -0.0013 and 40.0002, 0.0006 and 29.9994, -0.0002 and 10.0010 |
| Boss top (cap) | z 22 | 21.9987 |
| Boss wall | r 8 about (20, 15), axis z | r 8.0005 about (20.0003, 15.0002), axis within 0.009 degrees of z |
| Top-edge rounds | r 1, four edges | r 1.0059 to 1.0127, tangent to their faces |
| Junction round | torus, r 2 | refused: 885 samples in one connected region, all of them torus samples |

The recipe written is a block 40 × 30 × 10, a boss of radius 8 and height 12 standing on its top, and one fillet of radius 1.008 on the four top edges, built in both kernels. The junction stays sharp: the classic fit has no torus, so the part carries no round there and says so. Of the 1,205 torus samples, the boss wall and the top face explain the 320 within four noise deviations of them, near the round's tangent lines.

Both deviation numbers locate that difference. Scan to part: maximum 0.584 mm, RMS 0.054 mm, 95th percentile 0.025 mm; every sample away from the torus is within 0.06 mm, and the torus samples reach the expected 2 (1 - 1/√2) = 0.586 mm. Part to scan: maximum 1.027 mm, RMS 0.243 mm, 95th percentile 0.423 mm over 845,673 part samples; the sharp junction lies 2 (√2 - 1) = 0.828 mm from the round, and elsewhere the figure measures the spacing of 30,000 random samples. The deviation map, the scan and a text summary are raw evidence retained privately; `cargo run --release -p catallax-field --example scan_fit -- p3 30000 0.01 3 --map <dir>` reproduces them on the same platform: two rows of three orthographic views in a PNG, scan to part on every scan sample and part to scan on every 15th part sample, and the same values as coloured PLY points.

Seeds 1 to 16 all give the six block faces, the cap, the boss cylinder and the four rounds, with 820 to 941 torus samples refused and no unplaced surface. The regression `crates/catallax-field/tests/scan_p3_round_trip.rs` pins seed 3 with the bounds above (faces and boss within 0.005 mm, rounds within 0.02 mm, scan to part within six noise deviations off the torus) and checks identical sample coordinates, fitted surfaces and residuals, assignments, unexplained regions, refusal records, recipe text, both complete deviation maps and statistics, decision statistics and work counts at 1, 2 and 4 workers.

**What this says about learned segmentation.** On P3 the classic route recovers every plane and cylinder. What it cannot explain is the junction torus, a blend between a plane and a cylinder: that needs a torus or rolling-ball fit, not a different segmentation. This is evidence for the verdict, not the verdict; the written verdict on learned segmentation such as DualBrep waits for M10 and M11, with free-form regions.
