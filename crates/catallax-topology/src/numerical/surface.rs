//! Definition-held spline sheets. Trims are exact parameter rectangles; extension
//! and normal thickening retain their numerical checks and generating sheet.

use super::*;
use catallax_shapes::spline_tools::{RationalPoint, RationalSplineCurve, RationalSplineSurface};
use catallax_shapes::{BSplineJet, BSplineSurface, RationalPatch};
use catallax_trace::{spline as trace, surface_extend, vector::Vector2};

/// An editable finite spline sheet, admitted as a surface body without material.
#[derive(Clone, Debug)]
pub struct SplineSurfaceBody {
    source: RationalSplineSurface,
    body: NumericalBody,
}

fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("surface-body.construction", detail)
}

fn reserve<T>(count: usize, dec: &mut Decider) -> Result<Vec<T>, Refusal> {
    dec.charge_work(u64::try_from(count).map_err(|_| failure("storage count overflows"))?)?;
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|e| failure(e.to_string()))?;
    Ok(out)
}

fn product(a: usize, b: usize) -> Result<usize, Refusal> {
    a.checked_mul(b)
        .ok_or_else(|| failure("surface construction dimensions overflow"))
}

fn domain(source: &RationalSplineSurface) -> [[Rational; 2]; 2] {
    source.domain().map(|axis| axis.map(Clone::clone))
}

fn corners(source: &RationalSplineSurface) -> [[Rational; 2]; 4] {
    let [[a, b], [c, d]] = domain(source);
    [
        [a.clone(), c.clone()],
        [b.clone(), c],
        [b, d.clone()],
        [a, d],
    ]
}

fn positions(source: &RationalSplineSurface) -> [RationalPoint; 4] {
    let [u, v] = source.control_count();
    [0, (u - 1) * v, u * v - 1, v - 1].map(|i| source.controls()[i].clone())
}

fn boundary(
    source: &RationalSplineSurface,
    side: usize,
    dec: &mut Decider,
) -> Result<RationalSplineCurve, Refusal> {
    let [u, v] = source.control_count();
    let axis = usize::from(side == 1 || side == 3);
    let indices: Vec<_> = match side {
        0 => (0..u).map(|i| i * v).collect(),
        1 => (0..v).map(|j| (u - 1) * v + j).collect(),
        2 => (0..u).map(|i| i * v + v - 1).collect(),
        3 => (0..v).collect(),
        _ => unreachable!(),
    };
    dec.charge_work((indices.len() * 4) as u64)?;
    let curve = RationalSplineCurve::new(
        source.degrees()[axis],
        source.knots()[axis].to_vec(),
        indices
            .iter()
            .map(|i| source.controls()[*i].clone())
            .collect(),
        indices
            .iter()
            .map(|i| source.weights()[*i].clone())
            .collect(),
        dec,
    )?;
    let curve = if side >= 2 {
        curve.reversed(dec)?
    } else {
        curve
    };
    curve.reparameterized([&Rational::zero(), &Rational::one()], dec)
}

fn support(source: &RationalSplineSurface, dec: &mut Decider) -> Result<Surface, Refusal> {
    Ok(Surface::BSpline(
        RationalPatch::new(source.to_surface(dec)?, dec)?.into(),
    ))
}

fn definition(
    name: &str,
    source: &RationalSplineSurface,
    dec: &mut Decider,
) -> Result<NumericalDefinition, Refusal> {
    let uv = corners(source);
    let vertices = positions(source)
        .into_iter()
        .enumerate()
        .map(|(i, position)| NumericalVertex {
            key: format!("{name}/corner/{i}"),
            position,
        })
        .collect();
    let mut edges = reserve(4, dec)?;
    let mut trims = reserve(4, dec)?;
    for i in 0..4 {
        edges.push(NumericalEdge {
            key: format!("{name}/boundary/{i}"),
            curve: Curve::BSpline(boundary(source, i, dec)?.to_curve(dec)?),
            vertices: [i, (i + 1) % 4],
            parameters: [Rational::zero(), Rational::one()],
        });
        trims.push(NumericalTrim {
            edge: i,
            sense: Sense::Same,
            uv: [uv[i].clone(), uv[(i + 1) % 4].clone()],
        });
    }
    Ok(NumericalDefinition {
        name: name.into(),
        vertices,
        edges,
        faces: vec![NumericalFace {
            key: name.into(),
            surface: support(source, dec)?,
            sense: Sense::Same,
            loops: vec![trims],
        }],
    })
}

fn jet(
    source: &BSplineSurface,
    uv: &[Rational; 2],
    workspace: &mut [Vec<[Ball; 4]>; 2],
    dec: &mut Decider,
) -> Result<BSplineJet, Refusal> {
    dec.charge_work(1)?;
    let count = source.control_count();
    let span = std::array::from_fn(|a| {
        if &uv[a] == source.domain()[a][1] {
            count[a] - 1
        } else {
            source.knots()[a].partition_point(|k| k <= &uv[a]) - 1
        }
    });
    source.jet_on_span(
        span,
        &uv.clone().map(|q| [q.clone(), q]),
        BITS,
        dec,
        workspace,
    )
}

fn offset_point(
    value: &BSplineJet,
    thickness: &Rational,
    dec: &mut Decider,
) -> Result<[Ball; 3], Refusal> {
    dec.charge_work(24)?;
    let [u, v] = &value.first;
    let n = [
        &u[1] * &v[2] - &u[2] * &v[1],
        &u[2] * &v[0] - &u[0] * &v[2],
        &u[0] * &v[1] - &u[1] * &v[0],
    ];
    let length = (n[0].square() + n[1].square() + n[2].square())
        .sqrt()
        .ok_or_else(|| failure("surface normal magnitude is unresolved"))?;
    let scale = Ball::from_rational(thickness, BITS)
        .checked_div(&length)
        .ok_or_else(|| failure("surface normal may vanish"))?;
    Ok(std::array::from_fn(|a| &value.value[a] + &(&n[a] * &scale)))
}

fn parameter(lo: &Rational, hi: &Rational, i: usize, subdivisions: usize) -> Rational {
    lo + &(&(hi - lo) * &Rational::from_ratio(i as i64, subdivisions as i64).unwrap())
}

fn append_checks(body: &mut NumericalBody, extra: Vec<NumericalCheck>) -> Result<(), Refusal> {
    let mut checks = body.evidence.checks().to_vec();
    checks.extend(extra);
    body.evidence = NumericalEvidence::new(
        body.options.tolerance_mm.clone(),
        checks,
        body.evidence.repairs().to_vec(),
    )?;
    Ok(())
}

impl SplineSurfaceBody {
    /// Builds one bounded four-sided spline sheet, with physical boundary curves.
    /// The generating source remains exact; admission is numerical, not a solid proof.
    pub fn new(
        name: &str,
        source: RationalSplineSurface,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        guarded(|| {
            let source = source.into_clamped(dec)?;
            let body = build_surface(definition(name, &source, dec)?, options, dec)?;
            Ok(Self { source, body })
        })
    }

    /// The retained finite spline definition, including exact knots and weights.
    pub fn surface(&self) -> &RationalSplineSurface {
        &self.source
    }

    /// The surface-kind B-rep and its scoped checks, suitable for numerical STEP.
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }

    /// Restricts the sheet to an exact parameter rectangle and rebuilds its trims.
    /// Outside or empty rectangles refuse without altering the input sheet.
    pub fn trim(&self, rectangle: [[Rational; 2]; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        guarded(|| {
            let source = self
                .source
                .extracted(rectangle.each_ref().map(|a| a.each_ref()), dec)?;
            let mut next = Self::new(
                &self.body.definition.name,
                source,
                self.body.options.clone(),
                dec,
            )?;
            next.body.inherit(&self.body);
            Ok(next)
        })
    }

    /// Continues one chart boundary tangentially to an external parameter value.
    /// Reuses the existing SMLib G1 ruled-extension candidate. Independent samples
    /// check the entire old rectangle's grid, not only the newly moved corners.
    /// No whole-domain approximation or tangent-continuity proof is claimed.
    pub fn extend(
        &self,
        axis: usize,
        target: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        guarded(|| {
            if axis > 1 {
                return Err(failure("extension axis must be zero or one"));
            }
            let old_domain = domain(&self.source);
            if target >= old_domain[axis][0] && target <= old_domain[axis][1] {
                return Err(failure(
                    "extension parameter must be outside the active interval",
                ));
            }
            let candidate = extend_candidate(&self.source, axis, &target, dec)?;
            let old = self.source.to_surface(dec)?;
            let new = candidate.to_surface(dec)?;
            let mut old_work = old.span_workspace();
            let mut new_work = new.span_workspace();
            let mut samples = reserve(81, dec)?;
            for i in 0..=8 {
                for j in 0..=8 {
                    let uv = [
                        parameter(&old_domain[0][0], &old_domain[0][1], i, 8),
                        parameter(&old_domain[1][0], &old_domain[1][1], j, 8),
                    ];
                    let before = jet(&old, &uv, &mut old_work, dec)?;
                    let after = jet(&new, &uv, &mut new_work, dec)?;
                    samples.push(NumericalSample {
                        parameters: uv.to_vec(),
                        distance_upper_mm: distance(&before.value, &after.value),
                    });
                }
            }
            let mut next = Self::new(
                &self.body.definition.name,
                candidate,
                self.body.options.clone(),
                dec,
            )?;
            append_checks(
                &mut next.body,
                vec![NumericalCheck::Samples {
                    name: "surface.extension-retained-source-samples".into(),
                    scope: self.body.definition.name.clone(),
                    samples,
                }],
            )?;
            next.body.inherit(&self.body);
            Ok(next)
        })
    }

    /// Thickens along the unit normal by a nonzero signed distance. The source
    /// cap remains exact. Cubic interpolation uses every source knot span;
    /// independent per-span checks include both one-sided knot limits.
    /// Exact ruled rims close the shell without requiring equal cap bases.
    /// Checks are sampled, not a whole-body deviation or injectivity proof.
    pub fn thicken(
        &self,
        thickness: Rational,
        dec: &mut Decider,
    ) -> Result<NumericalBody, Refusal> {
        guarded(|| {
            if thickness.is_zero() {
                return Err(failure("thickness must be nonzero"));
            }
            let (offset, checks) = fit_normal_offset(
                &self.source,
                &thickness,
                &self.body.options.tolerance_mm,
                dec,
            )?;
            let front = thickness.sign() == Sign::Positive;
            let (lower, upper) = if front {
                (&self.source, &offset)
            } else {
                (&offset, &self.source)
            };
            let mut result = build(
                solid_definition(&self.body.definition.name, lower, upper, front, dec)?,
                self.body.options.clone(),
                dec,
            )?;
            append_checks(&mut result, checks)?;
            result.inherit(&self.body);
            Ok(result)
        })
    }
}

fn extend_candidate(
    source: &RationalSplineSurface,
    axis: usize,
    target: &Rational,
    dec: &mut Decider,
) -> Result<RationalSplineSurface, Refusal> {
    let domains = domain(source);
    let before = target < &domains[axis][0];
    let mut extent = domains[axis].clone();
    let degree = source.degrees()[axis];
    let mut count = source.control_count()[axis];
    while target < &extent[0] || target > &extent[1] {
        dec.charge_work(4)?;
        count = count
            .checked_add(degree)
            .ok_or_else(|| failure("extension count overflows"))?;
        let size = product(count, source.control_count()[1 - axis])?;
        let cost = product(product(size, degree + 1)?, product(degree + 1, 8)?)?;
        dec.charge_work(u64::try_from(cost).map_err(|_| failure("extension work overflows"))?)?;
        let step = (&extent[1] - &extent[0]) * Rational::from_ratio(9, 10).unwrap();
        if before {
            extent[0] = &extent[0] - &step;
        } else {
            extent[1] = &extent[1] + &step;
        }
    }
    let original_floats = domains
        .each_ref()
        .map(|d| d.each_ref().map(|q| q.to_f64_nearest()));
    let knots = source
        .knots()
        .map(|k| k.iter().map(Rational::to_f64_nearest).collect());
    let mut controls = reserve(source.controls().len(), dec)?;
    for (p, w) in source.controls().iter().zip(source.weights()) {
        dec.charge_work(4)?;
        controls.push(trace::ControlPoint {
            coordinates: [
                (&p[0] * w).to_f64_nearest(),
                (&p[1] * w).to_f64_nearest(),
                (&p[2] * w).to_f64_nearest(),
                w.to_f64_nearest(),
            ],
        });
    }
    let raw = trace::BSplineSurface::new(source.degrees(), knots, source.control_count(), controls)
        .map_err(|e| failure(format!("extension input: {e:?}")))?;
    let mut uv = [
        (original_floats[0][0] + original_floats[0][1]) * 0.5,
        (original_floats[1][0] + original_floats[1][1]) * 0.5,
    ];
    let share = (target - &domains[axis][0]) / (&domains[axis][1] - &domains[axis][0]);
    let mapped = original_floats[axis][0]
        + share.to_f64_nearest() * (original_floats[axis][1] - original_floats[axis][0]);
    uv[axis] = if before {
        mapped.next_down()
    } else {
        mapped.next_up()
    };
    let mut extension = None;
    surface_extend::evaluate_unbounded(&raw, &mut extension, Vector2::new(uv[0], uv[1]))
        .map_err(|e| failure(format!("tangent extension: {e:?}")))?;
    let extended = extension
        .ok_or_else(|| failure("periodic extension did not create a finite candidate"))?
        .surface()
        .map_err(|e| failure(format!("extension result: {e:?}")))?;
    let rational = |x| {
        Rational::from_f64_exact(x).ok_or_else(|| failure("extension generated a nonfinite value"))
    };
    let mut output_knots = [Vec::new(), Vec::new()];
    for a in 0..2 {
        let flo = rational(original_floats[a][0])?;
        let fhi = rational(original_floats[a][1])?;
        let scale = (&domains[a][1] - &domains[a][0]) / (&fhi - &flo);
        output_knots[a] = reserve(extended.knots()[a].len(), dec)?;
        for x in extended.knots()[a] {
            dec.charge_work(3)?;
            output_knots[a].push(&domains[a][0] + &(&(&rational(*x)? - &flo) * &scale));
        }
    }
    let mut points = reserve(extended.control_points().len(), dec)?;
    let mut weights = reserve(extended.control_points().len(), dec)?;
    for cp in extended.control_points() {
        dec.charge_work(4)?;
        let (p, w) = cp.get(trace::PointForm::Euclidean);
        points.push([rational(p.x)?, rational(p.y)?, rational(p.z)?]);
        weights.push(rational(w)?);
    }
    let candidate = RationalSplineSurface::new(
        extended.degrees(),
        output_knots,
        extended.counts(),
        points,
        weights,
        dec,
    )?;
    let mut requested = domains;
    requested[axis][usize::from(!before)] = target.clone();
    candidate.extracted(requested.each_ref().map(|a| a.each_ref()), dec)
}

struct OffsetInterval {
    source_span: usize,
    bounds: [Rational; 2],
}

fn offset_intervals(
    source: &RationalSplineSurface,
    axis: usize,
    divisions: usize,
    dec: &mut Decider,
) -> Result<Vec<OffsetInterval>, Refusal> {
    let knots = source.knots()[axis];
    let first = source.degrees()[axis];
    let last = source.control_count()[axis];
    dec.charge_work(((last - first) as u64).saturating_mul(2))?;
    let count = (first..last).filter(|i| knots[*i] != knots[*i + 1]).count();
    let mut intervals = reserve(product(count, divisions)?, dec)?;
    for span in first..last {
        if knots[span] == knots[span + 1] {
            continue;
        }
        for part in 0..divisions {
            intervals.push(OffsetInterval {
                source_span: span,
                bounds: [
                    parameter(&knots[span], &knots[span + 1], part, divisions),
                    parameter(&knots[span], &knots[span + 1], part + 1, divisions),
                ],
            });
        }
    }
    Ok(intervals)
}

fn offset_knots(intervals: &[OffsetInterval], dec: &mut Decider) -> Result<Vec<Rational>, Refusal> {
    let first = intervals
        .first()
        .ok_or_else(|| failure("source has no active knot span"))?;
    let count = intervals
        .len()
        .checked_add(7)
        .ok_or_else(|| failure("offset knot count overflows"))?;
    let mut knots = reserve(count, dec)?;
    knots.extend(std::iter::repeat_n(first.bounds[0].clone(), 4));
    knots.extend(
        intervals
            .iter()
            .skip(1)
            .map(|interval| interval.bounds[0].clone()),
    );
    knots.extend(std::iter::repeat_n(
        intervals.last().expect("active intervals exist").bounds[1].clone(),
        4,
    ));
    Ok(knots)
}

fn offset_nodes(knots: &[Rational], dec: &mut Decider) -> Result<Vec<Rational>, Refusal> {
    let count = knots.len() - 4;
    let mut nodes = reserve(count, dec)?;
    let three = Rational::from_i64(3);
    for i in 0..count {
        dec.charge_work(3)?;
        let sum = &(&knots[i + 1] + &knots[i + 2]) + &knots[i + 3];
        nodes.push(&sum / &three);
    }
    Ok(nodes)
}

fn fit_normal_offset(
    source: &RationalSplineSurface,
    thickness: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(RationalSplineSurface, Vec<NumericalCheck>), Refusal> {
    let carrier = source.to_surface(dec)?;
    let mut workspace = carrier.span_workspace();
    let mut divisions = source.degrees().map(|degree| degree.div_ceil(3).max(1));
    loop {
        let intervals = [
            offset_intervals(source, 0, divisions[0], dec)?,
            offset_intervals(source, 1, divisions[1], dec)?,
        ];
        let knots = [
            offset_knots(&intervals[0], dec)?,
            offset_knots(&intervals[1], dec)?,
        ];
        let parameters = [offset_nodes(&knots[0], dec)?, offset_nodes(&knots[1], dec)?];
        let counts = parameters.each_ref().map(|nodes| nodes.len());
        let size = product(counts[0], counts[1])?;
        let mut samples = reserve(size, dec)?;
        for u in &parameters[0] {
            for v in &parameters[1] {
                let value = jet(&carrier, &[u.clone(), v.clone()], &mut workspace, dec)?;
                samples.push(
                    offset_point(&value, thickness, dec)?
                        .each_ref()
                        .map(Ball::mid),
                );
            }
        }
        let offset = RationalSplineSurface::through_grid_on_knots(
            [&parameters[0], &parameters[1]],
            [3, 3],
            knots,
            &samples,
            dec,
        )?;
        let candidate = offset.to_surface(dec)?;
        let mut candidate_work = candidate.span_workspace();
        let mut checks = reserve(product(intervals[0].len(), intervals[1].len())?, dec)?;
        let mut maximum = Rational::zero();
        for (ui, u) in intervals[0].iter().enumerate() {
            for (vi, v) in intervals[1].iter().enumerate() {
                let mut evidence = reserve(49, dec)?;
                for i in 0..=6 {
                    for j in 0..=6 {
                        let uv = [
                            parameter(&u.bounds[0], &u.bounds[1], i, 6),
                            parameter(&v.bounds[0], &v.bounds[1], j, 6),
                        ];
                        dec.charge_work(1)?;
                        let original = carrier.jet_on_span(
                            [u.source_span, v.source_span],
                            &uv.clone().map(|q| [q.clone(), q]),
                            BITS,
                            dec,
                            &mut workspace,
                        )?;
                        let actual = jet(&candidate, &uv, &mut candidate_work, dec)?;
                        let error =
                            distance(&offset_point(&original, thickness, dec)?, &actual.value);
                        maximum = maximum.max(error.clone());
                        evidence.push(NumericalSample {
                            parameters: uv.to_vec(),
                            distance_upper_mm: error,
                        });
                    }
                }
                checks.push(NumericalCheck::Samples {
                    name: "surface.thicken-normal-offset-samples".into(),
                    scope: format!(
                        "source spans {},{}; interpolation cell {ui},{vi}",
                        u.source_span, v.source_span
                    ),
                    samples: evidence,
                });
            }
        }
        if maximum <= *tolerance {
            return Ok((offset, checks));
        }
        for divisions in &mut divisions {
            *divisions = divisions
                .checked_mul(2)
                .ok_or_else(|| failure("offset refinement dimensions overflow"))?;
        }
    }
}

fn solid_definition(
    name: &str,
    lower: &RationalSplineSurface,
    upper: &RationalSplineSurface,
    source_on_lower: bool,
    dec: &mut Decider,
) -> Result<NumericalDefinition, Refusal> {
    let role = |layer| {
        if (layer == 0) == source_on_lower {
            "source"
        } else {
            "offset"
        }
    };
    let mut vertices = reserve(8, dec)?;
    for (layer, sheet) in [lower, upper].into_iter().enumerate() {
        for (i, position) in positions(sheet).into_iter().enumerate() {
            vertices.push(NumericalVertex {
                key: format!("{name}/{}/corner/{i}", role(layer)),
                position,
            });
        }
    }
    let mut edges = reserve(12, dec)?;
    let mut boundaries = [Vec::new(), Vec::new()];
    for (layer, sheet) in [lower, upper].into_iter().enumerate() {
        for i in 0..4 {
            let curve = boundary(sheet, i, dec)?;
            edges.push(NumericalEdge {
                key: format!("{name}/{}/boundary/{i}", role(layer)),
                curve: Curve::BSpline(curve.to_curve(dec)?),
                vertices: [layer * 4 + i, layer * 4 + (i + 1) % 4],
                parameters: [Rational::zero(), Rational::one()],
            });
            boundaries[layer].push(curve);
        }
    }
    for i in 0..4 {
        let a = &vertices[i].position;
        let b = &vertices[i + 4].position;
        let curve = RationalSplineCurve::new(
            1,
            vec![
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::one(),
            ],
            vec![a.clone(), b.clone()],
            vec![Rational::one(); 2],
            dec,
        )?;
        edges.push(NumericalEdge {
            key: format!("{name}/rim/{i}"),
            curve: Curve::BSpline(curve.to_curve(dec)?),
            vertices: [i, i + 4],
            parameters: [Rational::zero(), Rational::one()],
        });
    }
    let mut faces = reserve(6, dec)?;
    for (layer, sheet) in [lower, upper].into_iter().enumerate() {
        let uv = corners(sheet);
        let order = if layer == 0 {
            [3, 2, 1, 0]
        } else {
            [0, 1, 2, 3]
        };
        let trims = order
            .into_iter()
            .map(|i| NumericalTrim {
                edge: layer * 4 + i,
                sense: if layer == 0 {
                    Sense::Reversed
                } else {
                    Sense::Same
                },
                uv: if layer == 0 {
                    [uv[(i + 1) % 4].clone(), uv[i].clone()]
                } else {
                    [uv[i].clone(), uv[(i + 1) % 4].clone()]
                },
            })
            .collect();
        faces.push(NumericalFace {
            key: if role(layer) == "source" {
                name.into()
            } else {
                format!("{name}/offset")
            },
            surface: support(sheet, dec)?,
            sense: if layer == 0 {
                Sense::Reversed
            } else {
                Sense::Same
            },
            loops: vec![trims],
        });
    }
    let zero = Rational::zero();
    let one = Rational::one();
    let rim_uv = [
        [zero.clone(), zero.clone()],
        [one.clone(), zero.clone()],
        [one.clone(), one.clone()],
        [zero.clone(), one.clone()],
    ];
    let [lower_boundaries, upper_boundaries] = boundaries;
    for (i, (bottom, top)) in lower_boundaries
        .into_iter()
        .zip(upper_boundaries)
        .enumerate()
    {
        let rim = RationalSplineSurface::ruled(&bottom, &top, dec)?;
        let walk = [
            (i, Sense::Same),
            (8 + (i + 1) % 4, Sense::Same),
            (4 + i, Sense::Reversed),
            (8 + i, Sense::Reversed),
        ];
        let trims = walk
            .into_iter()
            .enumerate()
            .map(|(j, (edge, sense))| NumericalTrim {
                edge,
                sense,
                uv: [rim_uv[j].clone(), rim_uv[(j + 1) % 4].clone()],
            })
            .collect();
        faces.push(NumericalFace {
            key: format!("{name}/wall/{i}"),
            surface: support(&rim, dec)?,
            sense: Sense::Same,
            loops: vec![trims],
        });
    }
    Ok(NumericalDefinition {
        name: format!("{name}/thicken"),
        vertices,
        edges,
        faces,
    })
}
