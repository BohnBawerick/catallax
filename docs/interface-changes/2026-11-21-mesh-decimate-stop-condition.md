## 2026-11-21: mesh decimation states the stop condition its queue gives

Pieces: catallax-recipe, catallax-import, catallax-features
Change: corrects the stop condition stated with [the mesh decimation operation](2026-11-21-mesh-boolean-and-decimation.md). `mesh_decimate` and `catallax_import::mesh::decimate::decimate` stop when at most `triangles` remain or their queue of collapses is empty, not when no allowed collapse is left. Each edge is queued at one placement, a refused collapse leaves the queue, and a refused edge comes back only when a collapse moves one of its ends ([RECIPE.md section 4.22](../RECIPE.md#422-mesh-decimation)). Results, work and refusals are unchanged.
Why: the earlier text promised more than the queue checks, and a thin torus of 144 triangles asked for one stops at 40 with allowed collapses left. SMLib's `SmPolyDecimate::DoQuadricDecimation` works the same way.
Migration: none in code. A caller who needs fewer triangles after a stop above the count can add a second `mesh_decimate` step, which starts a new queue.
