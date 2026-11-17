use super::*;

pub(super) struct State<'a> {
    pub vertices: &'a mut [NativeSample],
    pub edges: &'a mut [Option<EdgeSamples>],
    pub vertex_tolerances: &'a mut [Accumulator],
    pub edge_tolerances: &'a mut [Accumulator],
}

fn segment(a: [[f64; 3]; 2], b: [[f64; 3]; 2]) -> Option<[f64; 2]> {
    let u = subtract(a[1], a[0]);
    let v = subtract(b[1], b[0]);
    let w = subtract(b[0], a[0]);
    let (aa, ab, bb) = (dot(u, u), dot(u, v), dot(v, v));
    let determinant = aa * bb - ab * ab;
    if !determinant.is_finite() || determinant <= 128.0 * f64::EPSILON * aa * bb {
        return None;
    }
    let aw = dot(u, w);
    let bw = dot(v, w);
    let p = [
        (aw * bb - bw * ab) / determinant,
        (aw * ab - bw * aa) / determinant,
    ];
    p.iter()
        .all(|v| v.is_finite() && *v >= 0.0 && *v <= 1.0)
        .then_some(p)
}

fn evaluate(edge: &EdgeSamples, parameter: f64, dec: &mut Decider) -> Result<NativeSample, Error> {
    Ok(edge
        .carrier
        .sample_at(parameter, edge.points[0].0.chart, dec)?
        .1)
}

fn intersection(
    a: &EdgeSamples,
    b: &EdgeSamples,
    mut p: [f64; 2],
    precision: f64,
    dec: &mut Decider,
) -> Result<Option<([f64; 2], NativeSample)>, Error> {
    let limits = [a, b].map(|edge| {
        let first = edge.points[0].0.parameters[0];
        let last = edge.points.last().expect("finite edge").0.parameters[0];
        [first.min(last), first.max(last)]
    });
    for _ in 0..16 {
        dec.charge_work(1)?;
        let x = evaluate(a, p[0], dec)?;
        let y = evaluate(b, p[1], dec)?;
        let gap = norm(subtract(x.point_mm, y.point_mm)) + x.uncertainty_mm + y.uncertainty_mm;
        if gap <= precision {
            return Ok(Some((
                p,
                NativeSample {
                    point_mm: std::array::from_fn(|axis| {
                        x.point_mm[axis] * 0.5 + y.point_mm[axis] * 0.5
                    }),
                    uncertainty_mm: gap,
                },
            )));
        }
        let mut derivatives = [[0.0; 3]; 2];
        for (axis, edge) in [a, b].into_iter().enumerate() {
            let delta = ((limits[axis][1] - limits[axis][0]) * 0.00001)
                .max(64.0 * f64::EPSILON * p[axis].abs().max(1.0));
            let low = (p[axis] - delta).max(limits[axis][0]);
            let high = (p[axis] + delta).min(limits[axis][1]);
            if high <= low {
                return Ok(None);
            }
            let low_point = evaluate(edge, low, dec)?;
            let high_point = evaluate(edge, high, dec)?;
            derivatives[axis] =
                subtract(high_point.point_mm, low_point.point_mm).map(|v| v / (high - low));
        }
        let delta = subtract(y.point_mm, x.point_mm);
        let [u, v] = derivatives;
        let (aa, ab, bb) = (dot(u, u), dot(u, v), dot(v, v));
        let determinant = aa * bb - ab * ab;
        if !determinant.is_finite() || determinant <= 128.0 * f64::EPSILON * aa * bb {
            return Ok(None);
        }
        let (aw, bw) = (dot(u, delta), dot(v, delta));
        let next = [
            p[0] + (aw * bb - bw * ab) / determinant,
            p[1] + (aw * ab - bw * aa) / determinant,
        ];
        let next = std::array::from_fn(|axis| next[axis].clamp(limits[axis][0], limits[axis][1]));
        if next.map(f64::to_bits) == p.map(f64::to_bits) {
            return Ok(None);
        }
        p = next;
    }
    Ok(None)
}

pub(super) fn cut(
    context: &mut Context<'_>,
    body: &mut TolerantBody,
    state: State<'_>,
    evidence: &mut TolerantEvidence,
    dec: &mut Decider,
) -> Result<(), Error> {
    let precision = body.options.initial_mm.to_f64_nearest();
    for shell in 0..body.shells.len() {
        for face in 0..body.shells[shell].faces.len() {
            for boundary in 0..body.shells[shell].faces[face].loops.len() {
                let count = body.shells[shell].faces[face].loops[boundary].edges.len();
                for usage in 0..count {
                    let left = body.shells[shell].faces[face].loops[boundary].edges[usage].clone();
                    let right = body.shells[shell].faces[face].loops[boundary].edges
                        [(usage + 1) % count]
                        .clone();
                    if left.edge == right.edge {
                        continue;
                    }
                    let vertex = body.edges[left.edge].vertices[usize::from(left.forward)];
                    if body.edges[right.edge].vertices[usize::from(!right.forward)] != vertex {
                        continue;
                    }
                    let a = state.edges[left.edge].as_ref().expect("active edge");
                    let b = state.edges[right.edge].as_ref().expect("active edge");
                    let endpoint_parameters = [a, b].map(|edge| {
                        [
                            edge.points[0].0.parameters[0],
                            edge.points.last().expect("finite edge").0.parameters[0],
                        ]
                    });
                    let mut crossing: Option<([f64; 2], NativeSample)> = None;
                    for x in a.points.windows(2) {
                        for y in b.points.windows(2) {
                            dec.charge_work(1)?;
                            let Some(fractions) = segment(
                                [x[0].1.point_mm, x[1].1.point_mm],
                                [y[0].1.point_mm, y[1].1.point_mm],
                            ) else {
                                continue;
                            };
                            let seed = [
                                x[0].0.parameters[0]
                                    + fractions[0] * (x[1].0.parameters[0] - x[0].0.parameters[0]),
                                y[0].0.parameters[0]
                                    + fractions[1] * (y[1].0.parameters[0] - y[0].0.parameters[0]),
                            ];
                            let Some(found) = intersection(a, b, seed, precision, dec)? else {
                                continue;
                            };
                            let ordinary_endpoint = (0..2).any(|first| {
                                (0..2).any(|second| {
                                    body.edges[left.edge].vertices[first]
                                        == body.edges[right.edge].vertices[second]
                                        && [first, second].into_iter().enumerate().all(
                                            |(axis, endpoint)| {
                                                let value = endpoint_parameters[axis][endpoint];
                                                let width = (endpoint_parameters[axis][1]
                                                    - endpoint_parameters[axis][0])
                                                    .abs();
                                                (found.0[axis] - value).abs()
                                                    <= 64.0
                                                        * f64::EPSILON
                                                        * width
                                                            .max(value.abs())
                                                            .max(found.0[axis].abs())
                                            },
                                        )
                                })
                            });
                            if ordinary_endpoint {
                                continue;
                            }
                            let displacement =
                                norm(subtract(found.1.point_mm, state.vertices[vertex].point_mm))
                                    + found.1.uncertainty_mm
                                    + state.vertices[vertex].uncertainty_mm;
                            if displacement <= precision {
                                continue;
                            }
                            if displacement > context.repair_cap {
                                return Err(invalid(
                                    "import.tolerant-crossing-cap",
                                    "cutting an adjacent wire crossing exceeds the stated repair width",
                                ));
                            }
                            if let Some(previous) = crossing
                                && norm(subtract(previous.1.point_mm, found.1.point_mm)) > precision
                            {
                                return Err(invalid(
                                    "import.tolerant-ambiguous-crossing",
                                    "multiple adjacent wire intersections require a material-side decision",
                                ));
                            }
                            crossing = Some(found);
                        }
                    }
                    let Some((parameters, point)) = crossing else {
                        continue;
                    };
                    let mut ranges = [a, b].map(|edge| {
                        [
                            edge.points[0].0.parameters[0],
                            edge.points.last().expect("finite").0.parameters[0],
                        ]
                    });
                    ranges[0][usize::from(left.forward)] = parameters[0];
                    ranges[1][usize::from(!right.forward)] = parameters[1];
                    body.vertices[vertex].source.point = point
                        .point_mm
                        .map(|v| Rational::from_f64_exact(v).expect("finite intersection"));
                    state.vertices[vertex] = point;
                    for (edge, range) in [left.edge, right.edge].into_iter().zip(ranges) {
                        body.edges[edge].source_range =
                            Some(range.map(|v| Rational::from_f64_exact(v).expect("finite trim")));
                    }
                    evidence.repairs.push(
                        super::super::super::healing::HealingChange::WireCrossingCut {
                            body: body.name.clone(),
                            face: body.shells[shell].faces[face].name.clone(),
                            boundary,
                            edges: [
                                body.edges[left.edge].name.clone(),
                                body.edges[right.edge].name.clone(),
                            ],
                            point_mm: body.vertices[vertex]
                                .source
                                .point
                                .each_ref()
                                .map(ToString::to_string),
                            parameters: parameters.map(|v| format!("{v:.17e}")),
                        },
                    );
                    for index in 0..body.edges.len() {
                        if state.edges[index].is_some()
                            && body.edges[index].vertices.contains(&vertex)
                        {
                            let (samples, tolerance) = finite_edge(
                                context,
                                body,
                                index,
                                state.vertices,
                                state.vertex_tolerances,
                                evidence,
                                dec,
                            )?;
                            state.edges[index] = Some(samples);
                            state.edge_tolerances[index] = tolerance;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent_crossing_cuts_real_finite_ranges_and_refuses_outside_the_cap() {
        let source = b"ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
#1=CARTESIAN_POINT('',(-1.,0.,0.));#2=CARTESIAN_POINT('',(0.01,-0.01,0.));
#3=CARTESIAN_POINT('',(0.,1.,0.));#4=CARTESIAN_POINT('',(0.,0.,0.));
#5=DIRECTION('',(1.,0.,0.));#6=DIRECTION('',(0.,1.,0.));#7=DIRECTION('',(-1.,-1.,0.));
#8=DIRECTION('',(0.,0.,1.));#9=VECTOR('',#5,1.);#10=VECTOR('',#6,1.);#11=VECTOR('',#7,1.4142135623730951);
#12=LINE('',#1,#9);#13=LINE('',#4,#10);#14=LINE('',#3,#11);
#15=VERTEX_POINT('',#1);#16=VERTEX_POINT('',#2);#17=VERTEX_POINT('',#3);
#18=EDGE_CURVE('',#15,#16,#12,.T.);#19=EDGE_CURVE('',#16,#17,#13,.T.);#20=EDGE_CURVE('',#17,#15,#14,.T.);
#21=ORIENTED_EDGE('',*,*,#18,.T.);#22=ORIENTED_EDGE('',*,*,#19,.T.);#23=ORIENTED_EDGE('',*,*,#20,.T.);
#24=EDGE_LOOP('',(#21,#22,#23));#25=FACE_OUTER_BOUND('',#24,.T.);
#26=AXIS2_PLACEMENT_3D('',#4,#8,#5);#27=PLANE('',#26);#28=ADVANCED_FACE('',(#25),#27,.T.);
#29=OPEN_SHELL('',(#28));#30=SHELL_BASED_SURFACE_MODEL('',(#29));
ENDSEC;END-ISO-10303-21;";
        let mut dec = Decider::default();
        let parsed = ParsedStep::read(source, &mut dec).unwrap();
        let scale = Rational::one();
        for cap in [0.02, 0.005] {
            let retained = TolerantBody::retained(
                "#30",
                false,
                parsed.retained_body(30, &scale),
                &Default::default(),
            );
            let mut body = retained.materialize(&mut dec).unwrap().into_owned();
            for edge in &mut body.edges {
                let Curve::Retained(geometry) = &edge.curve else {
                    unreachable!()
                };
                edge.curve = geometry.recognized_line(&mut dec).unwrap().unwrap();
            }
            let original = body.clone();
            let mut context =
                Context::new(&parsed, &scale, AngleUnit::Radians, BTreeMap::new(), cap);
            let mut evidence = TolerantEvidence {
                method: "crossing-regression",
                source_sha256: parsed.source_sha256().into(),
                authority: "sampled",
                tolerances: vec![],
                sampled_checks: vec![],
                repairs: vec![],
                topology_checks: vec![],
                topology_measurements: vec![],
                warnings: vec![],
                native_trims: vec![],
                native_domains: vec![],
            };
            let mut vertices = body
                .vertices
                .iter()
                .map(|vertex| point(&vertex.source.point).unwrap())
                .collect::<Vec<_>>();
            let mut vertex_tolerances = vec![Accumulator::new(0.00001); vertices.len()];
            let mut edges = Vec::new();
            let mut edge_tolerances = Vec::new();
            for index in 0..body.edges.len() {
                let (samples, tolerance) = finite_edge(
                    &mut context,
                    &body,
                    index,
                    &vertices,
                    &mut vertex_tolerances,
                    &mut evidence,
                    &mut dec,
                )
                .unwrap();
                edges.push(Some(samples));
                edge_tolerances.push(tolerance);
            }
            let result = cut(
                &mut context,
                &mut body,
                State {
                    vertices: &mut vertices,
                    edges: &mut edges,
                    vertex_tolerances: &mut vertex_tolerances,
                    edge_tolerances: &mut edge_tolerances,
                },
                &mut evidence,
                &mut dec,
            );
            if cap < 0.01 {
                assert!(
                    matches!(result, Err(Error::Numerical(error)) if error.code == "import.tolerant-crossing-cap")
                );
                assert_eq!(body, original);
            } else {
                result.unwrap();
                assert!(norm(point(&body.vertices[1].source.point).unwrap().point_mm) < 1e-12);
                assert!(
                    (body.edges[0].source_range.as_ref().unwrap()[1].to_f64_nearest() - 1.0).abs()
                        < 1e-12
                );
                assert!(
                    body.edges[1].source_range.as_ref().unwrap()[0]
                        .to_f64_nearest()
                        .abs()
                        < 1e-12
                );
                assert!(matches!(
                    evidence.repairs.as_slice(),
                    [super::super::super::super::healing::HealingChange::WireCrossingCut { .. }]
                ));
                assert_ne!(
                    body.vertices[1].source.point,
                    original.vertices[1].source.point
                );
            }
        }
    }
}
