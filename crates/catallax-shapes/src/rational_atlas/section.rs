//! Complete analytic sections over every chart of a joined atlas face.

use super::*;
use crate::rational_patch::{Cover, search_cover, vanishing_side_row};
use crate::{PatchSectionBranch, SectionSearch, SectionSurface};

/// A complete checked section of one logical atlas face. Every component is a
/// cell-cover branch over the face's charts, joined across chart seams.
#[derive(Clone, Debug)]
pub struct AtlasSection {
    branches: Vec<PatchSectionBranch>,
    search: SectionSearch,
}

impl AtlasSection {
    /// Exact connected components, open between physical borders or closed.
    pub fn branches(&self) -> &[PatchSectionBranch] {
        &self.branches
    }
    /// Cell counts and decision work over all charts.
    pub fn search(&self) -> &SectionSearch {
        &self.search
    }
}

impl RationalAtlas {
    /// Complete finite sections by a plane, cylinder or sphere. Each chart is
    /// searched in the source parameters; arc ends join exactly across chart
    /// seams and periodic rows. Open ends lie on physical face borders.
    pub fn sections(
        &self,
        surface: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<AtlasSection, Refusal> {
        surface.validate(dec)?;
        let charts: Arc<[RationalPatch]> =
            self.charts().iter().map(|c| c.patch().clone()).collect();
        for chart in charts.iter() {
            if vanishing_side_row(chart, surface, dec)? {
                return Err(Refusal::new(
                    "intersect.atlas-section-row",
                    "the section contains a whole chart boundary row",
                ));
            }
        }
        let sides = [
            SplineBoundary::UMin,
            SplineBoundary::UMax,
            SplineBoundary::VMin,
            SplineBoundary::VMax,
        ];
        let cover = Cover {
            physical: self
                .charts()
                .iter()
                .map(|c| sides.map(|side| c.neighbor(side).is_none()))
                .collect(),
            charts,
            wrap_u: self
                .has_periodic_join()
                .then(|| self.source().domain()[0].map(Clone::clone)),
        };
        let (branches, search) = search_cover(&cover, surface, dec)?;
        Ok(AtlasSection { branches, search })
    }
}
