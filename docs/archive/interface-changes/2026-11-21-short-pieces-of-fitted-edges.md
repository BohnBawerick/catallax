Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-21: short pieces of fitted edges in hidden-line views

Pieces: application (`catallax_app::views`).
Change: `hidden_line_view` no longer starts a fitted edge's ray from a point moved farther across the edge than the piece's own ends. A piece whose ends lie nearer its sample point than twice the fit's tolerance (or rule T1's floor, if larger) takes the visibility of the piece before it, and the answer carries a new `VIEW_SHORT_PIECE_MARK` (`app.view-short-piece`) mark whose bound is the longest such piece. A fitted edge with no longer piece refuses with `app.view-sample`.
Why: post-merge review of the hidden-line views of intersection-branch edges. A short piece's first ray could start past an outline next to it and take the visibility on the outline's far side.
Migration: code that lists the marks of a hidden-line view should accept `app.view-short-piece`. A view's strokes change only where such a piece used to be labelled apart from the piece before it.
