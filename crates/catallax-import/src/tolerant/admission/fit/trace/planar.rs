use super::super::super::ProjectedTrace;
use super::super::*;
use crate::geometry_source::{Curve as SourceCurve, Geometry, Surface as SourceSurface};

struct Spline<'a> {
    degree: usize,
    knots: &'a [f64],
    weights: &'a [f64],
    controls: &'a [[f64; 2]],
}
impl Spline<'_> {
    fn value(&self, parameter: f64, work: &mut Vec<[f64; 3]>) -> Option<[f64; 2]> {
        let n = self.controls.len();
        let span = self
            .knots
            .partition_point(|knot| *knot <= parameter)
            .saturating_sub(1)
            .clamp(self.degree, n - 1);
        work.clear();
        for index in span - self.degree..=span {
            let w = self.weights[index];
            work.push([self.controls[index][0] * w, self.controls[index][1] * w, w]);
        }
        for level in 1..=self.degree {
            for j in (level..=self.degree).rev() {
                let index = span - self.degree + j;
                let width = self.knots[index + self.degree + 1 - level] - self.knots[index];
                let alpha = if width == 0.0 {
                    0.0
                } else {
                    (parameter - self.knots[index]) / width
                };
                work[j] = std::array::from_fn(|axis| {
                    (1.0 - alpha) * work[j - 1][axis] + alpha * work[j][axis]
                });
            }
        }
        let value = work[self.degree];
        let uv = [value[0] / value[2], value[1] / value[2]];
        (value[2] > 0.0 && uv.iter().all(|value| value.is_finite())).then_some(uv)
    }
}

pub(in super::super) fn project(
    context: &Context<'_>,
    surface: &NativeSurface,
    plane: Option<&PlaneFrame>,
    edge: &EdgeSamples,
    samples: &[UvSample],
    precision: f64,
    dec: &mut Decider,
) -> Result<Option<ProjectedTrace>, Error> {
    let Some(plane) = plane else {
        return Ok(None);
    };
    if !matches!(
        context.parsed.geometry(surface.source(), dec)?,
        Geometry::Surface(SourceSurface::Plane(_))
    ) {
        return Ok(None);
    }
    let Some(mut source) = edge.source_curve else {
        return Ok(None);
    };
    let mut basis = None;
    for _ in 0..128 {
        match context.parsed.geometry(source, dec)? {
            Geometry::Curve(SourceCurve::SurfaceCurve(parent, ..)) => source = parent,
            Geometry::Curve(SourceCurve::BSpline(_)) => {
                basis = Some(source);
                break;
            }
            Geometry::Curve(SourceCurve::Circle(..) | SourceCurve::Ellipse(..)) => {
                return conic(plane, surface, edge, samples, precision, dec);
            }
            _ => return Ok(None),
        }
    }
    let Some(source) = basis else {
        return Ok(None);
    };
    let curve = context
        .parsed
        .spline_curve(source, &std::array::from_fn(|_| context.scale.clone()), dec)
        .map_err(crate::tolerant::native::mapped)?;
    let mut controls = Vec::with_capacity(curve.controls().len());
    for control in curve.controls() {
        let projected = plane.project(point(control)?);
        controls.push(projected.location.parameters);
    }
    let largest = curve.weights().iter().max().expect("positive weights");
    let weights: Vec<_> = curve
        .weights()
        .iter()
        .map(|weight| (weight / largest).to_f64_nearest())
        .collect();
    let knots: Vec<_> = curve.knots().iter().map(Rational::to_f64_nearest).collect();
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight <= 0.0)
        || knots.iter().any(|knot| !knot.is_finite())
    {
        return Ok(None);
    }
    let spline = Spline {
        degree: curve.degree(),
        knots: &knots,
        weights: &weights,
        controls: &controls,
    };
    let mut work = Vec::with_capacity(curve.degree() + 1);
    for sample in samples {
        dec.charge_work(1)?;
        let Some(uv) = spline.value(sample.curve.parameters[0], &mut work) else {
            return Ok(None);
        };
        let actual = surface.sample(surface.location(uv)?, dec)?;
        let expected = surface.sample(sample.surface, dec)?;
        let error = norm(subtract(actual.point_mm, expected.point_mm))
            + actual.uncertainty_mm
            + expected.uncertainty_mm;
        if error > precision.max(64.0 * f64::EPSILON * norm(expected.point_mm).max(1.0)) {
            return Ok(None);
        }
    }
    Ok(Some(ProjectedTrace::Spline {
        source_curve: source,
        controls,
    }))
}

fn conic(
    plane: &PlaneFrame,
    surface: &NativeSurface,
    edge: &EdgeSamples,
    samples: &[UvSample],
    precision: f64,
    dec: &mut Decider,
) -> Result<Option<ProjectedTrace>, Error> {
    let PhysicalCurve::Native(curve) = &edge.carrier else {
        return Ok(None);
    };
    let period = curve.angular_period();
    let mut points = [[0.0; 2]; 3];
    for (index, point) in points.iter_mut().enumerate() {
        let (_, sample) = edge
            .carrier
            .sample_at(period * index as f64 / 4.0, 0, dec)?;
        *point = plane.project(sample).location.parameters;
    }
    let centre = std::array::from_fn(|axis| points[0][axis] * 0.5 + points[2][axis] * 0.5);
    let vectors = [points[0], points[1]]
        .map(|point| std::array::from_fn::<_, 2, _>(|axis| point[axis] - centre[axis]));
    let radii = vectors.map(|[x, y]| libm::hypot(x, y));
    if radii
        .iter()
        .any(|radius| !radius.is_finite() || *radius <= 0.0)
    {
        return Ok(None);
    }
    let axes =
        std::array::from_fn::<_, 2, _>(|index| vectors[index].map(|value| value / radii[index]));
    if (axes[0][0] * axes[1][0] + axes[0][1] * axes[1][1]).abs() > 128.0 * f64::EPSILON {
        return Ok(None);
    }
    for sample in samples {
        let angle = sample.curve.parameters[0] / period * std::f64::consts::TAU;
        let (sin, cos) = libm::sincos(angle);
        let uv = std::array::from_fn(|axis| {
            centre[axis] + vectors[0][axis] * cos + vectors[1][axis] * sin
        });
        let actual = surface.sample(surface.location(uv)?, dec)?;
        let expected = surface.sample(sample.surface, dec)?;
        let residual = norm(subtract(actual.point_mm, expected.point_mm))
            + actual.uncertainty_mm
            + expected.uncertainty_mm;
        if residual > precision.max(128.0 * f64::EPSILON * norm(expected.point_mm).max(1.0)) {
            return Ok(None);
        }
    }
    Ok(Some(ProjectedTrace::Conic {
        centre,
        axes,
        radii,
    }))
}
