## 2026-11-03: name the unfinished equal-cylinder crossing corner

Pieces: features, application, verification

Change: a failed cylinder/cylinder fillet now reports the existing `fillet.unsupported-corner` code when an exact proof identifies a collapsed equal-radius crossing pole on its selected input arc. Previously that case returned `fillet.unsupported-contact`. The detector lives in `fillet_equal_crossing.rs`; the shared builder adds only a named failure-path hook. No construction, successful body, surface type, function signature or recipe schema changes.

Why: for two equal cylinders on crossing axes, with equal positive signed offset radii, the two contact feet coincide at each pole. Their contact directions have a zero-area spherical cone. `SphereSector::new` requires three independent directions, so installing a finite spherical corner there would change the rolling-ball definition. The refused geometry needs explicit collapsed-boundary topology rather than a forced sphere patch. The existing SphereSector implementation is unchanged.

The proof requires exact equality of source radii and exact crossing of nonparallel axes, agreeing face senses and an actual selected-arc endpoint at a pole, or a complete closed branch. It excludes an inward offset that has already collapsed. A trimmed arc that stops before either pole and contacts on opposite ball sides keep their earlier refusal. Displaced-axis and unequal-radius regular fillets retain their successful construction paths. Modelling-tolerance coincidence alone cannot authorize this classification; missing proof preserves the original refusal.

Migration: callers that matched `fillet.unsupported-contact` for the proved equal-crossing pole case now receive `fillet.unsupported-corner` and no radius limit. This updates the current behavior described in [the strict contact-limit change](2026-11-03-fillet-contact-limit-equality.md), which retained the generic code for this case. The result still refuses; it publishes no new face, mesh, STEP or 100 nm face certificate. This is not proof that the recipe is impossible, and it does not enable singular fillet geometry.

SolidDraft must update its Catallax dependency and route `fillet.unsupported-corner` through its unsupported-feature UI, without presenting a radius-too-large limit or enabling this geometry. Existing successful recipes need no migration.
