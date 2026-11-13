//! Retained chart covers for the existing cell-pair search.

use super::*;
use crate::{RationalAtlas, SplineBoundary};

#[derive(Clone, Copy)]
pub(super) enum Source<'a> {
    Patch(&'a RationalPatch),
    Atlas(&'a RationalAtlas),
}

#[derive(Clone, Debug)]
pub(super) enum RetainedSource {
    Patch(RationalPatch),
    Atlas(RationalAtlas),
}

impl RetainedSource {
    pub(super) fn span(&self, chart: usize, span: usize) -> &Span {
        match self {
            Self::Patch(patch) => &patch.data.spans[span],
            Self::Atlas(atlas) => &atlas.charts()[chart].patch().data.spans[span],
        }
    }
}

impl<'a> Source<'a> {
    fn charts(self) -> impl Iterator<Item = &'a RationalPatch> {
        let count = match self {
            Self::Patch(_) => 1,
            Self::Atlas(atlas) => atlas.charts().len(),
        };
        (0..count).map(move |chart| match self {
            Self::Patch(patch) => patch,
            Self::Atlas(atlas) => atlas.charts()[chart].patch(),
        })
    }

    pub(super) fn spans(self) -> impl Iterator<Item = (usize, usize, &'a Span)> {
        self.charts().enumerate().flat_map(|(chart, patch)| {
            patch
                .data
                .spans
                .iter()
                .enumerate()
                .map(move |(index, span)| (chart, index, span))
        })
    }

    pub(super) fn retained(self) -> RetainedSource {
        match self {
            Self::Patch(patch) => RetainedSource::Patch(patch.clone()),
            Self::Atlas(atlas) => RetainedSource::Atlas(atlas.clone()),
        }
    }

    pub(super) fn physical(self, chart: usize, span: &Span) -> [bool; 4] {
        let sides = [
            SplineBoundary::UMin,
            SplineBoundary::UMax,
            SplineBoundary::VMin,
            SplineBoundary::VMax,
        ];
        let (patch, boundary) = match self {
            Self::Patch(patch) => (patch, [true; 4]),
            Self::Atlas(atlas) => {
                let chart = &atlas.charts()[chart];
                (
                    chart.patch(),
                    sides.map(|side| chart.neighbor(side).is_none()),
                )
            }
        };
        std::array::from_fn(|side| {
            boundary[side] && span.domain[side / 2][side % 2] == patch.domain()[side / 2][side % 2]
        })
    }

    pub(super) fn periodic(self) -> [Option<[Rational; 2]>; 2] {
        match self {
            Self::Patch(_) => [None, None],
            Self::Atlas(atlas) => std::array::from_fn(|axis| {
                atlas
                    .has_periodic_join(axis)
                    .then(|| atlas.source().domain()[axis].map(Clone::clone))
            }),
        }
    }
}

impl Context {
    /// Opposite periodic row values have one canonical stratum. The atlas
    /// constructor has already proved position and orientation along that row.
    pub(super) fn canonical<'a>(
        &'a self,
        axis: usize,
        value: &'a Rational,
        span: usize,
    ) -> &'a Rational {
        if let Some([lo, hi]) = &self.sides[axis / 2].periodic[axis % 2]
            && value == hi
            && !self.physical(axis, value, span)
        {
            lo
        } else {
            value
        }
    }
}

impl RationalAtlas {
    /// Complete meeting of two retained logical faces, with exact joins across
    /// their chart boundaries and identified periodic rows. The same per-cell
    /// field, normal-cone and root certificates as `RationalPatch::patch_meeting`
    /// decide every arc; the atlas supplies only its proved cover and ownership.
    pub fn patch_meeting(
        &self,
        other: &Self,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        meeting(
            [Source::Atlas(self), Source::Atlas(other)],
            [None, None],
            dec,
        )
    }

    /// Complete meeting with explicit reuse of exactly equal rational span
    /// polynomials and unrounded Bernstein nets. Every original chart and its
    /// ownership remain in the unchanged cell-pair search. The legacy method
    /// retains its original preparation and work.
    pub fn patch_meeting_with_cached_spans(
        &self,
        other: &Self,
        dec: &mut Decider,
    ) -> Result<PatchPairIntersection, Refusal> {
        prepared::meeting_cached([Source::Atlas(self), Source::Atlas(other)], dec)
    }
}
