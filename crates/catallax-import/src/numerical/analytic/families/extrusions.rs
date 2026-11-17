use super::*;
use crate::geometry_source::{Curve as SourceCurve, EntityId, ParsedStep};
use catallax_shapes::{BSplineCurve, Plane, Point3};

fn error(detail: &str) -> ImportError {
    Refusal::new("import.extrusion-unverified", detail).into()
}

fn conic(
    source: &ParsedStep,
    placement: EntityId,
    radii: [Rational; 2],
    scale: &Rational,
    dec: &mut Decider,
) -> Result<BSplineCurve, ImportError> {
    if radii.iter().any(|radius| radius.sign() != Sign::Positive) {
        return Err(error("extrusion profile radii must be positive"));
    }
    let frame = source
        .axis_placement(placement, scale, dec)
        .map_err(mapping)?;
    let x = frame.axes[0].scaled(&Real::from(&radii[0] * scale));
    let y = frame.axes[1].scaled(&Real::from(&radii[1] * scale));
    let nx = x.scaled(&Real::from_i64(-1));
    let ny = y.scaled(&Real::from_i64(-1));
    let first = frame.origin.plus(&x);
    let controls = vec![
        first.clone(),
        first.plus(&y),
        frame.origin.plus(&y),
        frame.origin.plus(&nx).plus(&y),
        frame.origin.plus(&nx),
        frame.origin.plus(&nx).plus(&ny),
        frame.origin.plus(&ny),
        first.plus(&ny),
        first,
    ];
    // Unequal endpoint weights keep the complete conic rational without a
    // rounded sqrt(2). Each quarter is the exact half-angle parametrization.
    Ok(BSplineCurve::new(
        2,
        [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
            .into_iter()
            .map(Rational::from_i64)
            .collect(),
        controls,
        [1, 1, 2, 1, 1, 1, 2, 1, 1]
            .into_iter()
            .map(Real::from_i64)
            .collect(),
        dec,
    )?)
}

fn profile(
    source: &ParsedStep,
    id: EntityId,
    scale: &Rational,
    dec: &mut Decider,
) -> Result<Curve, ImportError> {
    Ok(match source.geometry(id, dec)? {
        Geometry::Curve(SourceCurve::Circle(placement, radius)) => Curve::BSpline(conic(
            source,
            placement,
            [radius.clone(), radius],
            scale,
            dec,
        )?),
        Geometry::Curve(SourceCurve::Ellipse(placement, a, b)) => {
            Curve::BSpline(conic(source, placement, [a, b], scale, dec)?)
        }
        Geometry::Curve(SourceCurve::BSpline(_)) => Curve::BSpline(
            source
                .spline_curve(id, &std::array::from_fn(|_| scale.clone()), dec)
                .map_err(mapping)?
                .to_curve(dec)?,
        ),
        Geometry::Curve(SourceCurve::Line(_, _)) => {
            source.native_curve(id, scale, dec).map_err(mapping)?
        }
        _ => {
            return Err(error(
                "extrusion requires a retained line, conic or finite spline profile; wrappers are not discarded",
            ));
        }
    })
}

fn include(range: &mut Option<[Rational; 2]>, lo: Rational, hi: Rational) {
    if let Some([current_lo, current_hi]) = range {
        if lo < *current_lo {
            *current_lo = lo;
        }
        if hi > *current_hi {
            *current_hi = hi;
        }
    } else {
        *range = Some([lo, hi]);
    }
}

fn projection(p: &Point3, direction: &Vec3, dec: &mut Decider) -> Result<Ball, ImportError> {
    let mut value = Ball::from_rational(&Rational::zero(), BITS);
    for (coordinate, ratio) in [&p.x, &p.y, &p.z].into_iter().zip(direction) {
        value =
            value + dec.enclose_budgeted_at(coordinate, BITS)? * Ball::from_rational(ratio, BITS);
    }
    Ok(value)
}

pub(super) fn surface(
    source: &ParsedStep,
    id: EntityId,
    scale: &Rational,
    face: &crate::Face,
    dec: &mut Decider,
) -> Result<MappedSurface, ImportError> {
    use crate::planar_arithmetic as a;
    let Geometry::Surface(SourceSurface::Extrusion(curve, axis)) = source.geometry(id, dec)? else {
        return Err(error("expected SURFACE_OF_LINEAR_EXTRUSION"));
    };
    let Geometry::Vector(direction, magnitude) = source.geometry(axis, dec)? else {
        return Err(error("extrusion requires a VECTOR"));
    };
    let Geometry::Direction(ratios) = source.geometry(direction, dec)? else {
        return Err(error("extrusion VECTOR requires a DIRECTION"));
    };
    let mut direction: Vec3 = ratios
        .try_into()
        .map_err(|_| error("extrusion direction must be three-dimensional"))?;
    if magnitude.is_zero() || direction.iter().all(Rational::is_zero) {
        return Err(error("extrusion vector must not be zero"));
    }
    if magnitude.sign() == Sign::Negative {
        direction = direction.map(|ratio| -ratio);
    }
    let profile = profile(source, curve, scale, dec)?;
    let Curve::BSpline(profile) = profile else {
        let Curve::Line(line) = profile else {
            unreachable!()
        };
        let normal = line.dir.cross(&vector(&direction));
        return Ok(MappedSurface {
            carrier: Surface::Plane(Plane::new(line.origin, normal, line.dir, dec)?),
            source_checks: Vec::new(),
            source_periodic: [false; 2],
            ruled_v: true,
        });
    };
    let mut profile_bounds = None;
    for point in profile.controls() {
        dec.charge_work(1)?;
        let value = projection(point, &direction, dec)?;
        include(&mut profile_bounds, value.lo(), value.hi());
    }
    let [profile_lo, profile_hi] = profile_bounds.expect("validated spline has controls");
    let mut boundary = None;
    for bound in &face.loops {
        for edge in &bound.edges {
            let [lo, hi] = crate::swept::curve_range(edge, &direction, dec)?;
            include(&mut boundary, lo, hi);
        }
        if let Some(vertex) = &bound.vertex {
            let value = a::dot(&vertex.point, &direction, dec)?;
            include(&mut boundary, value.clone(), value);
        }
    }
    let [boundary_lo, boundary_hi] =
        boundary.ok_or_else(|| error("extrusion face has no finite boundary"))?;
    let norm = a::dot(&direction, &direction, dec)?;
    let lo = a::arithmetic(
        &a::arithmetic(&boundary_lo, &profile_hi, '-', dec)?,
        &norm,
        '/',
        dec,
    )?;
    let hi = a::arithmetic(
        &a::arithmetic(&boundary_hi, &profile_lo, '-', dec)?,
        &norm,
        '/',
        dec,
    )?;
    if lo >= hi {
        return Err(error(
            "extrusion face has no positive finite sweep interval",
        ));
    }
    let count = profile
        .control_count()
        .checked_mul(2)
        .ok_or(ImportRefusal::ResourceLimit)?;
    dec.charge_work(u64::try_from(count).map_err(|_| ImportRefusal::ResourceLimit)?)?;
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    controls
        .try_reserve_exact(count)
        .map_err(|_| ImportRefusal::ResourceLimit)?;
    weights
        .try_reserve_exact(count)
        .map_err(|_| ImportRefusal::ResourceLimit)?;
    let axis = vector(&direction);
    let lower = axis.scaled(&Real::from(lo.clone()));
    let upper = axis.scaled(&Real::from(hi.clone()));
    for (point, weight) in profile.controls().zip(profile.weights()) {
        controls.push(point.plus(&lower));
        controls.push(point.plus(&upper));
        weights.push(weight.clone());
        weights.push(weight.clone());
    }
    let periodic = [profile.is_closed() == Some(true), false];
    let carrier = BSplineSurface::new(
        [profile.degree(), 1],
        [
            profile.knots().to_vec(),
            vec![lo.clone(), lo, hi.clone(), hi],
        ],
        [profile.control_count(), 2],
        controls,
        weights,
        dec,
    )?;
    Ok(MappedSurface {
        carrier: checked_spline(carrier, periodic, dec)?,
        source_checks: Vec::new(),
        source_periodic: periodic,
        ruled_v: true,
    })
}
