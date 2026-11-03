## 2026-11-03: cheaper display and measure for a ring or loop fillet

Pieces: shapes, mesh, measure
Change: a rolling saddle's display grid now gives each fitted span its own column width, from that span's own curvature bounds, so a loop's mesh has fewer triangles at the same chord bound. The mesh no longer evaluates a rail's exact points when no planar face uses that rail, and a saddle's rails never meet a plane. A saddle's volume, area and centroid integrals now halve only the panels that need it, instead of every panel at each depth. Every chord bound and every measure enclosure keeps its meaning.
Why: a loop's display and measure cost several times a through ring's. Most of the display time went on exact rail points that no face read, and most of the measure time on panels that already met their target.
Migration: none. A saddle band's triangle count changes, and so may a loop measure's enclosure, within its requested width.
