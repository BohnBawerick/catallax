use super::super::super::{TolerantEdge, TolerantUse};
use super::*;

pub(super) fn collapsed(
    surface: &NativeSurface,
    ends: [[f64; 2]; 2],
    target: NativeSample,
    precision: f64,
    dec: &mut Decider,
) -> Result<Option<f64>, Error> {
    let mut substantial = false;
    for axis in 0..2 {
        let span = surface.periods()[axis].or_else(|| {
            let mut low = f64::INFINITY;
            let mut high = f64::NEG_INFINITY;
            for chart in 0..surface.chart_count() {
                let limits = &surface.domain(chart)[axis];
                low = low.min(limits[0].as_ref()?.to_f64_nearest());
                high = high.max(limits[1].as_ref()?.to_f64_nearest());
            }
            Some(high - low)
        });
        substantial |= span.is_some_and(|span| {
            span > 0.0 && (ends[1][axis] - ends[0][axis]).abs() >= span * 0.000001
        });
    }
    if !substantial {
        return Ok(None);
    }
    let mut maximum: f64 = 0.0;
    for station in 0..=STATIONS {
        dec.charge_work(1)?;
        let fraction = station as f64 / STATIONS as f64;
        let uv =
            std::array::from_fn(|axis| ends[0][axis] + fraction * (ends[1][axis] - ends[0][axis]));
        let sample = match surface
            .location(uv)
            .and_then(|location| surface.sample(location, dec))
        {
            Ok(sample) => sample,
            Err(error) if super::super::super::native::control(&error) => return Err(error),
            Err(_) => return Ok(None),
        };
        maximum = maximum.max(
            norm(subtract(sample.point_mm, target.point_mm))
                + sample.uncertainty_mm
                + target.uncertainty_mm,
        );
        if maximum > precision {
            return Ok(None);
        }
    }
    Ok(Some(maximum))
}

pub(super) fn align_endpoints(
    surface: &NativeSurface,
    samples: &mut [UvSample],
    targets: [NativeSample; 2],
    precision: f64,
    dec: &mut Decider,
) -> Result<bool, Error> {
    let mut changed = false;
    for (endpoint, (index, neighbor)) in [(0, 1), (samples.len() - 1, samples.len() - 2)]
        .into_iter()
        .enumerate()
    {
        for axis in 0..2 {
            let before = samples[index].surface.parameters;
            let mut after = before;
            after[axis] = samples[neighbor].surface.parameters[axis];
            if collapsed(surface, [before, after], targets[endpoint], precision, dec)?.is_some() {
                samples[index].surface = surface.location(after)?;
                changed = true;
            }
        }
    }
    Ok(changed)
}

pub(super) fn insert(
    context: &mut Context<'_>,
    body: &mut TolerantBody,
    evidence: &mut TolerantEvidence,
    dec: &mut Decider,
) -> Result<usize, Error> {
    let mut inserted = 0;
    for shell in 0..body.shells.len() {
        for face in 0..body.shells[shell].faces.len() {
            let name = body.shells[shell].faces[face].name.clone();
            let Some(entity) = id(&name) else {
                continue;
            };
            let Topology::Face {
                surface: Some(entity),
                ..
            } = context.parsed.topology(entity, dec)?
            else {
                continue;
            };
            let surface = context.surface(entity, dec)?;
            for boundary in 0..body.shells[shell].faces[face].loops.len() {
                let count = body.shells[shell].faces[face].loops[boundary].edges.len();
                if count == 0 {
                    continue;
                }
                let trims: Vec<_> = evidence
                    .native_trims
                    .iter()
                    .enumerate()
                    .filter_map(|(index, trim)| {
                        (trim.body == body.name && trim.face == name && trim.boundary == boundary)
                            .then_some(index)
                    })
                    .collect();
                if trims.len() != count {
                    return Err(invalid(
                        "import.tolerant-pole-trims",
                        "one trace is required for every boundary use",
                    ));
                }
                let mut uses = Vec::with_capacity(count);
                for usage in 0..count {
                    let original = &body.shells[shell].faces[face].loops[boundary].edges[usage];
                    let vertex = body.edges[original.edge].vertices[usize::from(original.forward)];
                    evidence.native_trims[trims[usage]].use_index = uses.len();
                    uses.push(original.clone());
                    let target = point(&body.vertices[vertex].source.point)?;
                    let first = evidence.native_trims[trims[usage]]
                        .samples
                        .last()
                        .expect("finite trace")
                        .surface;
                    let next = evidence.native_trims[trims[(usage + 1) % count]].samples[0].surface;
                    let Some(maximum) = collapsed(
                        &surface,
                        [first.parameters, next.parameters],
                        target,
                        body.options.initial_mm.to_f64_nearest(),
                        dec,
                    )?
                    else {
                        continue;
                    };
                    let edge_name = format!("{name}/boundary/{boundary}/pole/{usage}");
                    let index = body.edges.len();
                    let controls = vec![body.vertices[vertex].source.point.clone(); 2];
                    let geometry = catallax_shapes::spline_tools::RationalSplineCurve::new(
                        1,
                        vec![
                            Rational::zero(),
                            Rational::zero(),
                            Rational::one(),
                            Rational::one(),
                        ],
                        controls,
                        vec![Rational::one(); 2],
                        dec,
                    )?;
                    body.edge_indices.insert(edge_name.clone(), index);
                    body.edges.push(TolerantEdge {
                        name: edge_name.clone(),
                        curve: Curve::BSpline(crate::StepSplineCurve {
                            geometry,
                            closed: Some(true),
                            self_intersect: None,
                        }),
                        curve_same_sense: true,
                        vertices: [vertex; 2],
                        tolerance_mm: body.options.initial_mm.clone(),
                        source_range: None,
                        degenerate: true,
                    });
                    evidence.native_trims.push(NativeTrim {
                        body: body.name.clone(),
                        face: name.clone(),
                        boundary,
                        use_index: uses.len(),
                        edge: edge_name.clone(),
                        forward: true,
                        source_pcurve: None,
                        reconstructed: true,
                        projected: None,
                        samples: [first, next]
                            .into_iter()
                            .enumerate()
                            .map(|(endpoint, surface)| UvSample {
                                fraction: endpoint as f64,
                                curve: NativeLocation {
                                    chart: 0,
                                    parameters: [endpoint as f64, 0.0],
                                },
                                surface,
                                residual_mm: maximum,
                                converged: true,
                            })
                            .collect(),
                    });
                    evidence.tolerances.push(EntityTolerance {
                        body: body.name.clone(),
                        kind: "edge",
                        entity: edge_name.clone(),
                        tolerance_mm: body.options.initial_mm.to_f64_nearest(),
                        maximum_sample_residual_mm: maximum,
                        samples: STATIONS + 1,
                    });
                    evidence.repairs.push(
                        super::super::super::healing::HealingChange::DegenerateEdgeInserted {
                            body: body.name.clone(),
                            face: name.clone(),
                            boundary,
                            edge: edge_name,
                            vertex: body.vertices[vertex].source.name.clone(),
                            maximum_residual_mm: format!("{maximum:.17e}"),
                        },
                    );
                    uses.push(TolerantUse {
                        edge: index,
                        forward: true,
                    });
                    inserted += 1;
                }
                body.shells[shell].faces[face].loops[boundary].edges = uses;
            }
        }
    }
    Ok(inserted)
}
