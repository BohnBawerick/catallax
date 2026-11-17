use super::*;
use crate::geometry_source::{Geometry, MappingError, Surface as SourceSurface};
use catallax_numbers::NumericalCheck;
use catallax_shapes::{BSplineSurface, RationalAtlas};

mod extrusions;
mod offsets;

pub(super) struct MappedSurface {
    pub carrier: Surface,
    pub source_checks: Vec<NumericalCheck>,
    // Source identifications guide trim reconstruction; they do not certify
    // exact periodic joins of a subsequently approximated offset carrier.
    pub source_periodic: [bool; 2],
    pub ruled_v: bool,
}

pub(super) fn mapping(error: MappingError) -> ImportError {
    match error {
        MappingError::Source(error) => error.into(),
        MappingError::Geometry(error) => error.into(),
    }
}

fn checked_spline(
    carrier: BSplineSurface,
    periodic: [bool; 2],
    dec: &mut Decider,
) -> Result<Surface, ImportError> {
    let mut components = RationalAtlas::partition(carrier, periodic, dec)?;
    if components.len() != 1 {
        return Err(Refusal::new(
            "import.family-chart-components",
            "the complete source face requires more than one G1 atlas component; no component is discarded",
        ).into());
    }
    Ok(Surface::BSpline(
        components.pop().expect("one checked component").into(),
    ))
}

pub(super) fn surface(
    face: &crate::Face,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<MappedSurface>, ImportError> {
    match &face.surface {
        crate::Surface::BSpline(source) => {
            let periodic = source.closed.map(|closed| closed == Some(true));
            Ok(Some(MappedSurface {
                carrier: checked_spline(source.geometry.to_surface(dec)?, periodic, dec)?,
                source_checks: Vec::new(),
                source_periodic: periodic,
                ruled_v: false,
            }))
        }
        crate::Surface::Retained(retained) => {
            let source = retained.source();
            let id = retained.entity();
            let scale = retained.length_mm();
            match source.geometry(id, dec)? {
                Geometry::Surface(SourceSurface::Extrusion(_, _)) => {
                    extrusions::surface(source, id, scale, face, dec).map(Some)
                }
                Geometry::Surface(SourceSurface::Offset(_, _, _)) => {
                    offsets::surface(source, id, scale, face, options, dec).map(Some)
                }
                _ => Ok(None),
            }
        }
        _ => Ok(None),
    }
}
