//! Complete analytic sections over every chart of a joined atlas face.

use super::*;
use crate::rational_patch::{Cover, search_cover, vanishing_side_row};
use crate::{BSplineCurve, PatchSectionBranch, Point3, SectionSearch, SectionSurface};

mod axis;

/// A complete checked section of one logical atlas face, retaining rational
/// components or cell-cover branches joined across chart seams.
#[derive(Clone, Debug)]
pub struct AtlasSection {
    branches: Vec<PatchSectionBranch>,
    points: Vec<Point3>,
    curves: Vec<(BSplineCurve, bool)>,
    search: SectionSearch,
}

impl AtlasSection {
    /// Finite rational components and their exact odd/even crossing parity.
    pub fn curves(&self) -> impl Iterator<Item = (&BSplineCurve, bool)> {
        self.curves
            .iter()
            .map(|(curve, crossing)| (curve, *crossing))
    }
    /// Isolated contacts retained by the complete cell search.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }
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
    /// Complete finite sections by an admitted analytic surface. Each chart is
    /// searched in the source parameters; arc ends join exactly across chart
    /// seams and periodic rows. Open ends lie on physical face borders.
    pub fn sections(
        &self,
        surface: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<AtlasSection, Refusal> {
        self.section_search(surface, false, dec)
    }

    /// Complete sections whose retained cells also have strict power-form
    /// transverse derivative bounds for Taylor tracing. Refines the same
    /// section search without changing ordinary query admission or allowances.
    pub fn sections_for_tracing(
        &self,
        surface: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<AtlasSection, Refusal> {
        self.section_search(surface, true, dec)
    }

    fn section_search(
        &self,
        surface: &SectionSurface,
        trace_ready: bool,
        dec: &mut Decider,
    ) -> Result<AtlasSection, Refusal> {
        surface.validate(dec)?;
        let start = dec.stats().work;
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
        if matches!(surface, SectionSurface::Plane(_))
            && let Some(curves) = axis::section(self, surface, dec)?
        {
            return Ok(AtlasSection {
                branches: Vec::new(),
                points: Vec::new(),
                curves,
                search: SectionSearch {
                    work: dec.stats().work - start,
                    ..SectionSearch::default()
                },
            });
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
            wrap: std::array::from_fn(|axis| {
                self.has_periodic_join(axis)
                    .then(|| self.source().domain()[axis].map(Clone::clone))
            }),
        };
        let section = search_cover(&cover, surface, trace_ready, dec)?;
        Ok(AtlasSection {
            branches: section.branches,
            points: section.points,
            curves: Vec::new(),
            search: section.search,
        })
    }
}
