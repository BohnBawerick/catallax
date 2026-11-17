use super::*;
use crate::geometry_source::{EntityId, ParsedStep};
use catallax_shapes::Torus;

fn error(detail: &str) -> ImportError {
    Refusal::new("import.offset-basis-unverified", detail).into()
}

pub(super) fn surface(
    source: &ParsedStep,
    id: EntityId,
    scale: &Rational,
    face: &crate::Face,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<MappedSurface, ImportError> {
    use crate::planar_arithmetic as a;
    let Geometry::Surface(SourceSurface::Offset(parent, distance, _)) = source.geometry(id, dec)?
    else {
        return Err(error("expected OFFSET_SURFACE"));
    };
    let distance = a::arithmetic(&distance, scale, '*', dec)?;
    let base = source.geometry(parent, dec)?;
    if matches!(base, Geometry::Surface(SourceSurface::Offset(_, _, _))) {
        return Err(error(
            "nested offsets need a checked normal-error composition; distances are not simply added",
        ));
    }
    if matches!(
        base,
        Geometry::Surface(
            SourceSurface::Plane(_) | SourceSurface::Cylinder(_, _) | SourceSurface::Sphere(_, _)
        )
    ) {
        // OFFSET_SURFACE's final logical field describes self-intersection;
        // it does not reverse the basis parameter normal.
        return Ok(MappedSurface {
            carrier: source.native_surface(id, scale, dec).map_err(mapping)?,
            source_checks: Vec::new(),
            source_periodic: [false; 2],
            ruled_v: false,
        });
    }
    let mut basis = match base {
        Geometry::Surface(SourceSurface::Extrusion(_, _)) => {
            extrusions::surface(source, parent, scale, face, dec)?
        }
        Geometry::Surface(SourceSurface::BSpline(_)) => {
            let geometry = source.spline_surface(parent, scale, dec).map_err(mapping)?;
            MappedSurface {
                carrier: checked_spline(geometry.to_surface(dec)?, [false; 2], dec)?,
                source_checks: Vec::new(),
                source_periodic: [false; 2],
                ruled_v: false,
            }
        }
        Geometry::Surface(SourceSurface::Torus(_, major, minor)) => {
            if major.sign() != Sign::Positive || minor.sign() != Sign::Positive || minor >= major {
                return Err(error(
                    "offset torus requires the original positive-radius ring sheet",
                ));
            }
            let Surface::Torus(torus) =
                source.native_surface(parent, scale, dec).map_err(mapping)?
            else {
                return Err(error("offset torus basis has no complete ring support"));
            };
            return Ok(MappedSurface {
                carrier: Surface::Torus(Torus::new(
                    torus.center,
                    torus.axis,
                    torus.major,
                    torus.minor + Real::from(distance),
                    torus.ref_dir,
                    dec,
                )?),
                source_checks: Vec::new(),
                source_periodic: [true; 2],
                ruled_v: false,
            });
        }
        _ => {
            return Err(error(
                "offset basis requires an admitted plane, cylinder, sphere, ring torus, extrusion or finite spline carrier",
            ));
        }
    };
    if distance.is_zero() {
        return Ok(basis);
    }
    if let Surface::Plane(plane) = &mut basis.carrier {
        let normal = plane
            .normal
            .scaled(&(&Real::one() / &plane.normal.norm_squared().sqrt()));
        plane.origin = plane.origin.plus(&normal.scaled(&Real::from(distance)));
        return Ok(basis);
    }
    let Surface::BSpline(support) = &basis.carrier else {
        return Err(error(
            "normal interpolation requires the checked finite spline basis",
        ));
    };
    let fit_tolerance = a::arithmetic(&options.tolerance_mm, &Rational::from_i64(4), '/', dec)?;
    let (candidate, mut checks) =
        n::surface::normal_offset_surface(support.carrier(), &distance, &fit_tolerance, dec)?;
    for check in &mut checks {
        if let NumericalCheck::Samples { scope, .. } = check {
            *scope = format!("source offset face {}: {scope}", face.name);
        }
    }
    basis.source_checks.append(&mut checks);
    basis.carrier = checked_spline(candidate.to_surface(dec)?, [false; 2], dec)?;
    Ok(basis)
}
