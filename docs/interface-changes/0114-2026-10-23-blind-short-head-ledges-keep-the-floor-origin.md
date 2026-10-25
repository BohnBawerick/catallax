## 2026-10-23: blind short-head ledges keep the floor origin

Pieces: recipe, features, field
Change: a blind hole with a head can structurally reference `floor` even with a cone tip. The planar annulus at blind depth between an unfinished head and the narrower cone tip carries that origin. Its point set is unchanged.
Why: a legal short countersink or counterbore can leave this ledge; the earlier flat-tip-only structural rule gave it no usable name.
Migration: no new role or alias. Existing through-hole floor references remain invalid.

