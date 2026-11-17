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

pub(super) fn project(
    context: &Context<'_>,
    surface: &NativeSurface,
    edge: &EdgeSamples,
    samples: &[UvSample],
    precision: f64,
    dec: &mut Decider,
) -> Result<Option<ProjectedTrace>, Error> {
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
        let projected = surface.project(point(control)?.point_mm, None, dec)?;
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
    Ok(Some(ProjectedTrace {
        source_curve: source,
        controls,
    }))
}
