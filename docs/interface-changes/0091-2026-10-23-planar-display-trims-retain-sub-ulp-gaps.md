## 2026-10-23: planar display trims retain sub-ulp gaps

Pieces: display mesh
Change: planar faces now triangulate their trims with certified `Real` orientation and containment predicates before rounding display positions. Hole bridges and ears retain every sampled boundary vertex and shared edge index. Curved faces keep the periodic CDT path.
Why: a valid flange union leaves a gap smaller than one binary64 ulp between a cap's outer loop and a hole. Projecting rounded positions first made those constraints overlap and split without orientation labels. The exact body and Boolean are unchanged.
Migration: no API, model, STEP or chord-tolerance changes. A planar triangle's bound is the maximum certified endpoint displacement, since its exact support is affine. Sub-ulp strips may have coincident rounded positions but retain distinct indices and closed incidence. The minimized two-block union and multiple-bore regressions cover these boundaries.

