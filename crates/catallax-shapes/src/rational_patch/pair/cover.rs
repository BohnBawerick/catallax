//! Retained chart covers for the existing cell-pair search.

use super::original::OriginalAtlas;
use super::*;
use crate::{RationalAtlas, SplineBoundary};

#[derive(Clone, Copy)]
pub(super) enum Source<'a> {
    Patch(&'a RationalPatch),
    Atlas(&'a RationalAtlas),
    Original(&'a Arc<OriginalAtlas>),
}

#[derive(Clone, Debug)]
pub(super) enum RetainedSource {
    Patch(RationalPatch),
    Atlas(RationalAtlas),
    Original(Arc<OriginalAtlas>),
}

impl RetainedSource {
    fn source(&self) -> Source<'_> {
        match self {
            Self::Patch(patch) => Source::Patch(patch),
            Self::Atlas(atlas) => Source::Atlas(atlas),
            Self::Original(source) => Source::Original(source),
        }
    }

    pub(super) fn same_definition(&self, other: &Self) -> bool {
        if let (Self::Original(a), Self::Original(b)) = (self, other) {
            if !a.atlas.source().same_definition(b.atlas.source()) {
                return false;
            }
        } else if matches!(self, Self::Original(_)) || matches!(other, Self::Original(_)) {
            return false;
        }
        let (a, b) = (self.source(), other.source());
        if a.periodic() != b.periodic() {
            return false;
        }
        let (mut a, mut b) = (a.charts(), b.charts());
        loop {
            match (a.next(), b.next()) {
                (Some(a), Some(b)) => {
                    if a.data.domain != b.data.domain || !a.carrier().same_definition(b.carrier()) {
                        return false;
                    }
                }
                (None, None) => return true,
                _ => return false,
            }
        }
    }

    pub(super) fn query_charts(
        &self,
        unit: usize,
    ) -> impl Iterator<Item = (usize, &RationalPatch)> + '_ {
        let (single, cover) = match self {
            Self::Patch(patch) => (Some((0, patch)), None),
            Self::Atlas(atlas) => (Some((unit, atlas.charts()[unit].patch())), None),
            Self::Original(source) => (
                None,
                Some((&source.atlas, source.units[unit].charts.as_slice())),
            ),
        };
        single
            .into_iter()
            .chain(cover.into_iter().flat_map(|(atlas, indices)| {
                indices
                    .iter()
                    .map(move |&index| (index, atlas.charts()[index].patch()))
            }))
    }

    pub(super) fn span(&self, chart: usize, span: usize) -> &Span {
        match self {
            Self::Patch(patch) => &patch.data.spans[span],
            Self::Atlas(atlas) => &atlas.charts()[chart].patch().data.spans[span],
            Self::Original(source) => {
                assert_eq!(span, 0, "an original source unit contains one tensor span");
                &source.units[chart].span
            }
        }
    }
}

impl<'a> Source<'a> {
    fn charts(self) -> impl Iterator<Item = &'a RationalPatch> {
        let count = match self {
            Self::Patch(_) => 1,
            Self::Atlas(atlas) => atlas.charts().len(),
            Self::Original(source) => source.atlas.charts().len(),
        };
        (0..count).map(move |chart| match self {
            Self::Patch(patch) => patch,
            Self::Atlas(atlas) => atlas.charts()[chart].patch(),
            Self::Original(source) => source.atlas.charts()[chart].patch(),
        })
    }

    pub(super) fn spans(self) -> impl Iterator<Item = (usize, usize, &'a Span)> {
        let original = match self {
            Self::Original(source) => Some(source),
            _ => None,
        };
        let charted = if original.is_some() { None } else { Some(self) };
        charted
            .into_iter()
            .flat_map(|source| {
                source.charts().enumerate().flat_map(|(chart, patch)| {
                    patch
                        .data
                        .spans
                        .iter()
                        .enumerate()
                        .map(move |(index, span)| (chart, index, span))
                })
            })
            .chain(original.into_iter().flat_map(|source| {
                source
                    .units
                    .iter()
                    .enumerate()
                    .map(|(index, unit)| (index, 0, &unit.span))
            }))
    }

    pub(super) fn retained(self) -> RetainedSource {
        match self {
            Self::Patch(patch) => RetainedSource::Patch(patch.clone()),
            Self::Atlas(atlas) => RetainedSource::Atlas(atlas.clone()),
            Self::Original(source) => RetainedSource::Original(source.clone()),
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
            Self::Original(source) => return source.units[chart].physical,
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
        let atlas = match self {
            Self::Patch(_) => return [None, None],
            Self::Atlas(atlas) => atlas,
            Self::Original(source) => &source.atlas,
        };
        std::array::from_fn(|axis| {
            atlas
                .has_periodic_join(axis)
                .then(|| atlas.source().domain()[axis].map(Clone::clone))
        })
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
