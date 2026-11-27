use std::{collections::BTreeMap, f64::consts::TAU};

use catallax_numbers::{Rational, Refusal};
use catallax_shapes::{Curve, Surface};
use catallax_topology::{
    Sense,
    numerical::{NumericalBody, NumericalEdge},
};

use super::{
    DISPLAY_DESCRIPTION_VERSION, DisplayDescription, DisplayEdge, DisplayEdgeUse, DisplayFace,
    DisplayLoop, DisplayOptions, DisplayTrimSample, Uv, Xyz, edges::sample_span, evaluate::*,
    failed, spline::SplineEdge,
};

/// Read an admitted numerical B-rep without rechecking its geometry or fits.
/// Face and edge handles are the definition's stable keys, not exact-model IDs.
pub fn describe_numerical_body(
    body: &NumericalBody,
    options: &DisplayOptions,
) -> Result<DisplayDescription<String, String>, Refusal> {
    super::validate_options(options)?;
    let def = body.definition();
    let read = |p: &[Rational; 3]| p.each_ref().map(Rational::to_f64_nearest);
    let origin = def.vertices.first().map_or([0.0; 3], |v| read(&v.position));
    let vertices = def
        .vertices
        .iter()
        .map(|v| xyz(sub(read(&v.position), origin)))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| failed("numerical vertex exceeds binary64"))?;
    let mut result = DisplayDescription {
        version: DISPLAY_DESCRIPTION_VERSION,
        origin,
        edges: Vec::with_capacity(def.edges.len()),
        faces: Vec::with_capacity(def.faces.len()),
    };
    let mut parameters = Vec::with_capacity(def.edges.len());
    let mut remaining = options.max_edge_points;
    for e in &def.edges {
        let mut ts = Vec::new();
        let mut points = Vec::new();
        sample(
            e,
            e.vertices.map(|i| vertices[i]),
            origin,
            options.chord_tolerance,
            &mut |t, p| {
                if remaining == 0 {
                    return None;
                }
                remaining -= 1;
                points.push(xyz(p)?);
                ts.push(t);
                Some(())
            },
        )
        .ok_or_else(|| failed(format!("numerical edge {} could not be sampled", e.key)))?;
        parameters.push(ts);
        result.edges.push(DisplayEdge {
            edge: e.key.clone(),
            points,
        });
    }
    let trims = body
        .native_trims()
        .iter()
        .map(|b| ((b.face.as_str(), b.boundary, b.use_index), &b.curve))
        .collect::<BTreeMap<_, _>>();
    for face in &def.faces {
        let surface = describe(&face.surface, origin)
            .ok_or_else(|| failed(format!("numerical face {} has no float carrier", face.key)))?;
        let mut chart = Projector::new(&surface, Some(&face.surface))
            .ok_or_else(|| failed(format!("numerical face {} has no finite chart", face.key)))?;
        let units = match face.surface {
            Surface::Cone(_) | Surface::Cylinder(_) => [1.0 / TAU, 1.0],
            Surface::Sphere(_) | Surface::Torus(_) => [1.0 / TAU; 2],
            _ => [1.0; 2],
        };
        let into_chart = |uv: Uv| std::array::from_fn(|i| uv[i] * units[i]);
        let mut trim_bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
        let mut loops = Vec::with_capacity(face.loops.len());
        for (boundary, row) in face.loops.iter().enumerate() {
            let mut uses = Vec::with_capacity(row.len());
            for (use_index, trim) in row.iter().enumerate() {
                let native = trims
                    .get(&(face.key.as_str(), boundary, use_index))
                    .map(|curve| {
                        SplineEdge::new(curve, [0.0; 3]).ok_or_else(|| {
                            failed(format!(
                                "numerical face {} has an unrepresentable pcurve",
                                face.key
                            ))
                        })
                    })
                    .transpose()?;
                let mut eval_trim = native.as_ref().map(SplineEdge::evaluator);
                let ends = trim
                    .uv
                    .each_ref()
                    .map(|p| into_chart(p.each_ref().map(Rational::to_f64_nearest)));
                let edge = trim.edge;
                let ts = &parameters[edge];
                let points = &result.edges[edge].points;
                let reverse = trim.sense == Sense::Reversed;
                let [start, end] = if reverse {
                    [ts[ts.len() - 1], ts[0]]
                } else {
                    [ts[0], ts[ts.len() - 1]]
                };
                let projected = matches!(face.surface, Surface::Plane(_))
                    && !matches!(def.edges[edge].curve, Curve::Line(_));
                let mut samples = Vec::with_capacity(ts.len());
                for i in 0..ts.len() {
                    let point = if reverse { ts.len() - 1 - i } else { i };
                    let uv = if i == 0 {
                        ends[0]
                    } else if i == ts.len() - 1 {
                        ends[1]
                    } else if let Some(eval) = &mut eval_trim {
                        let [u, v, _] = eval(ts[point]);
                        into_chart([u, v])
                    } else if projected {
                        chart.project(points[point], None)
                    } else {
                        let t = (ts[point] - start) / (end - start);
                        std::array::from_fn(|i| ends[0][i] + t * (ends[1][i] - ends[0][i]))
                    };
                    let [u, v] = uv.map(finite);
                    let nonfinite =
                        || failed(format!("numerical face {} has nonfinite trim UV", face.key));
                    let uv = [u.ok_or_else(nonfinite)?, v.ok_or_else(nonfinite)?];
                    for i in 0..2 {
                        trim_bounds[i][0] = trim_bounds[i][0].min(uv[i]);
                        trim_bounds[i][1] = trim_bounds[i][1].max(uv[i]);
                    }
                    samples.push(DisplayTrimSample { point, uv });
                }
                uses.push(DisplayEdgeUse { edge, samples });
            }
            loops.push(DisplayLoop { uses, vertex: None });
        }
        let bounds = std::array::from_fn(|i| chart.domain[i].unwrap_or(trim_bounds[i]));
        if bounds
            .iter()
            .any(|[lo, hi]| !lo.is_finite() || !hi.is_finite() || lo >= hi)
        {
            return Err(failed(format!(
                "numerical face {} has an empty or unbounded display chart",
                face.key
            )));
        }
        let periods = chart.periods;
        result.faces.push(DisplayFace {
            face: face.key.clone(),
            reversed: face.sense == Sense::Reversed,
            surface,
            bounds,
            periods,
            loops,
        });
    }
    Ok(result)
}

fn xyz(p: Xyz) -> Option<Xyz> {
    let [x, y, z] = p.map(finite);
    Some([x?, y?, z?])
}

fn sample(
    edge: &NumericalEdge,
    ends: [Xyz; 2],
    origin: Xyz,
    tolerance: f64,
    out: &mut impl FnMut(f64, Xyz) -> Option<()>,
) -> Option<()> {
    let [start, end] = edge.parameters.each_ref().map(Rational::to_f64_nearest);
    finite(start)?;
    finite(end)?;
    if start == end {
        return None;
    }
    out(start, ends[0])?;
    match &edge.curve {
        Curve::Line(_) => out(end, ends[1]),
        Curve::Circle(c) => {
            let f = frame(&c.center, &c.normal, &c.ref_dir, origin)?;
            let radius = positive(libm::sqrt(scalar(&c.radius_sq)))?;
            let span = end - start;
            let step = libm::sqrt(8.0 * tolerance / radius).min(TAU / 8.0);
            let count = libm::ceil(span.abs() / step) as usize;
            for i in 1..count {
                let t = start + span * i as f64 / count as f64;
                let (s, c) = libm::sincos(t);
                out(
                    t,
                    add(
                        f.origin,
                        add(scale(f.x, radius * c), scale(f.y, radius * s)),
                    ),
                )?;
            }
            out(end, ends[1])
        }
        Curve::BSpline(s) => {
            let curve = SplineEdge::new(s, origin)?;
            let mut at = curve.evaluator();
            let mut a = start;
            let mut p = ends[0];
            let [lo, hi] = [start.min(end), start.max(end)];
            let knots = curve.breaks();
            let knots =
                &knots[knots.partition_point(|&k| k <= lo)..knots.partition_point(|&k| k < hi)];
            for b in (0..knots.len())
                .map(|i| knots[if start < end { i } else { knots.len() - 1 - i }])
                .chain([end])
            {
                if a == b {
                    continue;
                }
                let q = if b == end { ends[1] } else { at(b) };
                sample_span(&mut at, [a, b], [p, q], tolerance, out)?;
                (a, p) = (b, q);
            }
            Some(())
        }
        _ => None,
    }
}
