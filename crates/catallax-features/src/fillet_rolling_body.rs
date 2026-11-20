//! Topology publication after rolling contacts meet real neighbouring faces.

use super::rolling::RolledStrip;
use super::*;
use catallax_shapes::NumericalSurfaceCut;

pub(super) fn occurrence(
    d: &NumericalDefinition,
    face: usize,
    edge: usize,
) -> Result<(usize, usize), Refusal> {
    let mut found = None;
    for (boundary, row) in d.faces[face].loops.iter().enumerate() {
        for (use_index, trim) in row.iter().enumerate() {
            if trim.edge == edge && found.replace((boundary, use_index)).is_some() {
                return Err(bad("a support repeats the selected edge in its boundary"));
            }
        }
    }
    found.ok_or_else(|| bad("a rolling support lost its selected boundary occurrence"))
}

pub(super) fn material_sense(
    strip: &RolledStrip,
    tau: &Rational,
    dec: &mut Decider,
) -> Result<Sense, Refusal> {
    let half = Rational::from_ratio(1, 2).unwrap();
    let jet = strip
        .fillet
        .surface()
        .source()
        .jet_box([tau, tau], [&half, &half], 128, dec)?;
    let centre = strip
        .fillet
        .spine()
        .centre_curve()
        .evaluate_box([tau, tau], 128, dec)?;
    let normal: [Ball; 3] = std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&jet.first[0][j] * &jet.first[1][k]) - &(&jet.first[0][k] * &jet.first[1][j])
    });
    let direction = (0..3).fold(Ball::from_rational(&q(0), 128), |sum, i| {
        &sum + &(&normal[i] * &(&jet.value[i] - &centre[i]))
    });
    let outward = if direction.lo().sign() == Sign::Positive {
        true
    } else if direction.hi().sign() == Sign::Negative {
        false
    } else {
        return Err(bad(
            "the trimmed rolling surface has an unresolved material normal",
        ));
    };
    Ok(if outward == strip.convex {
        Sense::Same
    } else {
        Sense::Reversed
    })
}

/// Clips a strip to its two incident neighbouring faces and publishes the edit.
/// Junction sections are supplied separately by the network publication layer.
pub(crate) fn cap_strip(
    source: &NumericalBody,
    strip: &RolledStrip,
    key: &str,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    let original = source.definition();
    let edge_key = format!("{:?}", strip.source.edge);
    let selected = original
        .edges
        .iter()
        .position(|edge| edge.key == edge_key)
        .ok_or_else(|| bad("the rolling strip has no retained source edge"))?;
    let input = &original.edges[selected];
    let face_index = |id| {
        let face_key = format!("{id:?}");
        original
            .faces
            .iter()
            .position(|face| face.key == face_key)
            .ok_or_else(|| bad("a rolling support has no retained source face"))
    };
    let sides = [
        face_index(strip.source.faces[0])?,
        face_index(strip.source.faces[1])?,
    ];
    let occurrences = [
        occurrence(original, sides[0], selected)?,
        occurrence(original, sides[1], selected)?,
    ];
    let senses = std::array::from_fn::<_, 2, _>(|side| {
        original.faces[sides[side]].loops[occurrences[side].0][occurrences[side].1].sense
    });
    let caps = [
        end_face(original, input.vertices[0], sides)?,
        end_face(original, input.vertices[1], sides)?,
    ];
    let cut = |end: usize, dec: &mut Decider| -> Result<NumericalSurfaceCut, Refusal> {
        let face = &original.faces[caps[end]];
        let xyz = &original.vertices[input.vertices[end]].position;
        let point = Point3::new(
            xyz[0].clone().into(),
            xyz[1].clone().into(),
            xyz[2].clone().into(),
        );
        let seed =
            crate::numerical_embedding::surface_parameters(&face.surface, &point, tolerance, dec)?;
        strip
            .fillet
            .surface_cut(
                face.surface.clone(),
                [q(end as i64), seed[0].clone(), seed[1].clone()],
                dec,
            )
            .map_err(|error| at_stage(&format!("neighbour {} intersection", face.key), error))
    };
    let cuts = [cut(0, dec)?, cut(1, dec)?];
    let mut d = original.clone();
    // As in `fillet_numerical::build_impl`: earlier derived copies are
    // measured again.
    d.reset_copy_tags();
    d.name = format!("{} / {key}", original.name);
    let mut native = native_rows(source)?;
    let mut parameters: Parameters = BTreeMap::new();
    let mut vertices = [[0; 2]; 2];
    for side in 0..2 {
        let rail = strip
            .fillet
            .rail(side)
            .ok_or_else(|| bad("a rolling rail is missing"))?;
        if cuts[0].tau_ends()[side] >= cuts[1].tau_ends()[side] {
            return Err(bad(
                "neighbour contacts consume or reverse the rolling rail interval",
            ));
        }
        for end in 0..2 {
            let tau = &cuts[end].tau_ends()[side];
            let point = rail
                .evaluate_box([tau, tau], 128, dec)?
                .map(|value| value.mid());
            let vertex = d.vertices.len();
            vertices[side][end] = vertex;
            d.vertices.push(NumericalVertex {
                key: format!(
                    "{key}/{edge_key}/rail/{}/end/{end}",
                    original.faces[sides[side]].key
                ),
                position: point,
            });
            let border = neighbour(original, sides[side], selected, input.vertices[end])?;
            let t = project(
                &original.edges[border],
                &d.vertices[vertex].position,
                tolerance,
                dec,
            )
            .map_err(|error| {
                at_stage(
                    &format!("border {} contact", original.edges[border].key),
                    error,
                )
            })?;
            let slot = usize::from(original.edges[border].vertices[1] == input.vertices[end]);
            d.edges[border].vertices[slot] = vertex;
            d.edges[border].parameters[slot] = t.clone();
            for (face_index, face) in d.faces.iter_mut().enumerate() {
                for (boundary, row) in face.loops.iter_mut().enumerate() {
                    for (use_index, trim) in row.iter_mut().enumerate() {
                        if trim.edge != border {
                            continue;
                        }
                        let uv = match &native[face_index][boundary][use_index] {
                            Some(curve) => uv_at(curve, &t, dec)?,
                            None => plane_uv(&face.surface, &d.vertices[vertex].position, dec)?,
                        };
                        let endpoint = if trim.sense == Sense::Same {
                            slot
                        } else {
                            1 - slot
                        };
                        trim.uv[endpoint] = uv.clone();
                        parameters.insert((face_index, vertex), uv);
                    }
                }
            }
        }
    }
    let mut rails = [0; 2];
    for side in 0..2 {
        let interval = [
            cuts[0].tau_ends()[side].clone(),
            cuts[1].tau_ends()[side].clone(),
        ];
        let endpoints = vertices[side].map(|vertex| parameters[&(sides[side], vertex)].clone());
        let (pcurve, _) = constrained_rail(
            strip
                .fillet
                .parameter_rail(side)
                .ok_or_else(|| bad("a rolling contact pcurve is missing"))?,
            &interval,
            &endpoints,
            None,
            dec,
        )?;
        let (boundary, use_index) = occurrences[side];
        native[sides[side]][boundary][use_index] = Some(pcurve);
        rails[side] = d.edges.len();
        d.edges.push(NumericalEdge {
            key: format!("{key}/{edge_key}/rail/{}", original.faces[sides[side]].key),
            curve: Curve::BSpline(strip.fillet.rail(side).unwrap().clone()),
            vertices: vertices[side],
            parameters: interval,
        });
        let uv = walked(&d.edges[rails[side]], senses[side])
            .map(|vertex| parameters[&(sides[side], vertex)].clone());
        d.faces[sides[side]].loops[boundary][use_index] = new_trim(rails[side], senses[side], uv);
    }
    let mut cut_edges = [0; 2];
    for end in 0..2 {
        let face = caps[end];
        cut_edges[end] = d.edges.len();
        d.edges.push(NumericalEdge {
            key: format!("{key}/{edge_key}/end/{end}/{}", original.faces[face].key),
            curve: Curve::BSpline(cuts[end].curve().clone()),
            vertices: [vertices[0][end], vertices[1][end]],
            parameters: [q(0), q(1)],
        });
        let endpoints = [vertices[0][end], vertices[1][end]];
        let uv = endpoints.map(|vertex| parameters[&(face, vertex)].clone());
        let (pcurve, _) =
            constrained_rail(&cuts[end].parameters()[1], &[q(0), q(1)], &uv, None, dec)?;
        let old_loops = std::mem::take(&mut d.faces[face].loops);
        let old_native = std::mem::take(&mut native[face]);
        for (boundary, bindings) in old_loops.into_iter().zip(old_native) {
            let first = boundary
                .first()
                .ok_or_else(|| bad("a neighbour boundary is empty"))?;
            let first_uv = first.uv[0].clone();
            let first_vertex = walked(&d.edges[first.edge], first.sense)[0];
            let mut replacement = Vec::with_capacity(boundary.len() + 1);
            let mut replacement_native = Vec::with_capacity(boundary.len() + 1);
            let mut uses = boundary.into_iter().zip(bindings).peekable();
            while let Some((trim, curve)) = uses.next() {
                let a = walked(&d.edges[trim.edge], trim.sense)[1];
                let b = uses.peek().map_or(first_vertex, |(next, _)| {
                    walked(&d.edges[next.edge], next.sense)[0]
                });
                replacement.push(trim);
                replacement_native.push(curve);
                if a != b {
                    let sense = if [a, b] == endpoints {
                        Sense::Same
                    } else if [b, a] == endpoints {
                        Sense::Reversed
                    } else {
                        return Err(bad(
                            "a neighbour loop has a gap unrelated to this rolling strip",
                        ));
                    };
                    let next_uv = uses
                        .peek()
                        .map_or_else(|| first_uv.clone(), |(next, _)| next.uv[0].clone());
                    replacement.push(new_trim(
                        cut_edges[end],
                        sense,
                        [replacement.last().unwrap().uv[1].clone(), next_uv],
                    ));
                    replacement_native.push(Some(pcurve.clone()));
                }
            }
            d.faces[face].loops.push(replacement);
            native[face].push(replacement_native);
        }
    }
    let face_index = d.faces.len();
    for side in 0..2 {
        for end in 0..2 {
            parameters.insert(
                (face_index, vertices[side][end]),
                [cuts[end].tau_ends()[side].clone(), q(side as i64)],
            );
        }
    }
    let first = (rails[0], flip(senses[0]));
    let reaches = usize::from(first.1 == Sense::Same);
    let path = [
        first,
        (cut_edges[reaches], Sense::Same),
        (rails[1], flip(senses[1])),
        (cut_edges[1 - reaches], Sense::Reversed),
    ];
    let mut blend_native = Vec::with_capacity(4);
    let mut boundary = Vec::with_capacity(4);
    for (edge, sense) in path {
        let pcurve = if edge == rails[0] {
            parameter_line(strip.fillet.rail(0).unwrap().domain(), 0, dec)?
        } else if edge == rails[1] {
            parameter_line(strip.fillet.rail(1).unwrap().domain(), 1, dec)?
        } else {
            cuts[usize::from(edge == cut_edges[1])].parameters()[0].clone()
        };
        blend_native.push(Some(pcurve));
        let mut trim = new_trim(
            edge,
            sense,
            walked(&d.edges[edge], sense).map(|vertex| parameters[&(face_index, vertex)].clone()),
        );
        if rails.contains(&edge) {
            trim.tag = blend_row();
        }
        boundary.push(trim);
    }
    let face_key = format!("{key}/blend/{edge_key}");
    let tau = (&cuts[0].tau_ends()[0] + &cuts[1].tau_ends()[0]).mul_pow2(-1);
    d.faces.push(NumericalFace {
        key: face_key.clone(),
        surface: Surface::BSpline(SplineSupport::RationalAtlas(strip.fillet.surface().clone())),
        sense: material_sense(strip, &tau, dec)?,
        loops: vec![boundary],
    });
    native.push(vec![blend_native]);
    let mut sampled = source.sampled_carriers().clone();
    sampled.insert(face_key);
    publish(source, d, native, sampled, dec)
}

pub(super) fn publish(
    source: &NumericalBody,
    mut d: NumericalDefinition,
    native: NativeRows,
    sampled: BTreeSet<String>,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    let mut bindings = Vec::new();
    for (face, rows) in d.faces.iter().zip(native) {
        for (boundary, (row, curves)) in face.loops.iter().zip(rows).enumerate() {
            for (use_index, (trim, curve)) in row.iter().zip(curves).enumerate() {
                if let Some(curve) = curve {
                    bindings.push(NumericalNativeTrim {
                        // The u/v curve of a use is made by the same rule
                        // as its endpoints.
                        tag: trim.tag.clone(),
                        face: face.key.clone(),
                        boundary,
                        use_index,
                        edge: d.edges[trim.edge].key.clone(),
                        curve,
                    });
                }
            }
        }
    }
    compact(&mut d);
    source
        .rebuild_with_sampled_carriers(
            d,
            bindings,
            sampled,
            // The class cap of the body the strip is cut into (rule T2):
            // the fit aimed at `tolerance`, and admission measures again.
            source.options().clone(),
            dec,
        )
        .map_err(|error| at_stage("real-domain publication", error))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};
    use catallax_shapes::Pose;
    use catallax_topology::Model;
    use catallax_topology::numerical::NumericalOptions;

    fn planar_area(body: &NumericalBody, face: &NumericalFace, dec: &Decider) -> f64 {
        let mut integral = 0.0;
        let definition = body.definition();
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, trim) in row.iter().enumerate() {
                let binding = body.native_trims().iter().find(|binding| {
                    binding.face == face.key
                        && binding.boundary == boundary
                        && binding.use_index == use_index
                });
                let Some(binding) = binding else {
                    let uv = &trim.uv;
                    integral += uv[0][0].to_f64_nearest() * uv[1][1].to_f64_nearest()
                        - uv[1][0].to_f64_nearest() * uv[0][1].to_f64_nearest();
                    continue;
                };
                let edge = &definition.edges[trim.edge];
                let mut range = edge.parameters.clone();
                if trim.sense == Sense::Reversed {
                    range.reverse();
                }
                let direction = if range[0] < range[1] { 1.0 } else { -1.0 };
                range.sort();
                for (_, span) in binding.curve.spans() {
                    let lo = span[0].clone().max(range[0].clone());
                    let hi = span[1].clone().min(range[1].clone());
                    if lo >= hi {
                        continue;
                    }
                    let (lo, hi) = (lo.to_f64_nearest(), hi.to_f64_nearest());
                    let middle = (lo + hi) * 0.5;
                    let half = (hi - lo) * 0.5;
                    for (node, weight) in [
                        (-(3.0f64 / 5.0).sqrt(), 5.0 / 9.0),
                        (0.0, 8.0 / 9.0),
                        ((3.0f64 / 5.0).sqrt(), 5.0 / 9.0),
                    ] {
                        let t = Rational::from_f64_exact(middle + half * node).unwrap();
                        let jet = binding.curve.jet_box([&t, &t], 64, dec).unwrap();
                        let value = jet.value.map(|v| v.mid().to_f64_nearest());
                        let first = jet.first.map(|v| v.mid().to_f64_nearest());
                        integral +=
                            direction * half * weight * (value[0] * first[1] - value[1] * first[0]);
                    }
                }
            }
        }
        integral.abs() * 0.5
    }

    #[test]
    fn rolling_caps_replace_the_sharp_corner_in_every_original_plane_domain() {
        let (body, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let mut model = Model::default();
            let body = catallax_topology::block(&mut model, [20, 30, 40].map(q), &Pose::identity())
                .unwrap();
            let edge = model
                .edges_of_body(body)
                .into_iter()
                .find(|id| {
                    let edge = model.edge(*id);
                    [edge.start, edge.end].into_iter().all(|vertex| {
                        let p = model.geom().position(model.vertex(vertex.unwrap()).point);
                        p.x.as_rational().is_some_and(|v| v.is_zero())
                            && p.y.as_rational().is_some_and(|v| v.is_zero())
                    })
                })
                .unwrap();
            let tolerance = Rational::parse("0.0000001").unwrap();
            let (definition, native) = crate::numerical_embedding::definition_from_exact_model(
                &model,
                body,
                "one rolling edge".into(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let source = catallax_topology::numerical::build_with_native_trims(
                definition,
                native,
                NumericalOptions {
                    tolerance_mm: tolerance.clone(),
                    vertex_tolerance_mm: None,
                },
                &mut dec,
            )
            .unwrap();
            let chains = super::super::rolling::construct(
                &model,
                body,
                &[edge],
                &q(1),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            cap_strip(&source, &chains[0].strips[0], "round", &tolerance, &mut dec).unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        let dec = Decider::new(BudgetPolicy {
            max_work: 100_000_000,
            ..Default::default()
        });
        let cap_area = 600.0 - (1.0 - std::f64::consts::FRAC_PI_4);
        for (axis, coordinate, expected) in [
            (0, 0.0, 1160.0),
            (0, 20.0, 1200.0),
            (1, 0.0, 760.0),
            (1, 30.0, 800.0),
            (2, 0.0, cap_area),
            (2, 40.0, cap_area),
        ] {
            let face = body
                .definition()
                .faces
                .iter()
                .find(|face| {
                    let Surface::Plane(plane) = &face.surface else {
                        return false;
                    };
                    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
                    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
                    normal[axis].as_rational().is_some_and(|v| !v.is_zero())
                        && origin[axis]
                            .as_rational()
                            .is_some_and(|v| (v.to_f64_nearest() - coordinate).abs() < 1e-12)
                })
                .unwrap();
            let area = planar_area(&body, face, &dec);
            assert!(
                (area - expected).abs() <= 0.00001,
                "axis={axis}, coordinate={coordinate}: area={area}, expected={expected}"
            );
        }
        println!(
            "one rolling edge, real-domain caps and publication work={}",
            work.used
        );
    }

    /// Property test of the producer of blend-row tags on one rolled edge
    /// (decision Q8, section 4.4): a block of random size, one of its edges
    /// rounded at a random radius, and the full check on both rails.
    #[test]
    fn the_rails_of_a_rolled_edge_are_rows_of_their_blend_face() {
        let mut random = super::super::tag_audit::Stream::seeded("rolled edge");
        for _ in 0..2 {
            let size = [(); 3].map(|()| q(random.between(8, 40)));
            // A radius between 1/4 and 3 mm, in quarters.
            let radius = q(random.between(1, 12)).mul_pow2(-2);
            let (body, work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 100_000_000,
                    ..Default::default()
                });
                let mut model = Model::default();
                let body =
                    catallax_topology::block(&mut model, size.clone(), &Pose::identity()).unwrap();
                let edges = model.edges_of_body(body);
                let edge = edges[random.between(0, edges.len() as i64 - 1) as usize];
                let tolerance = Rational::parse("0.0000001").unwrap();
                let (definition, native) = crate::numerical_embedding::definition_from_exact_model(
                    &model,
                    body,
                    "one rolling edge".into(),
                    &tolerance,
                    &mut dec,
                )
                .unwrap();
                let source = catallax_topology::numerical::build_with_native_trims(
                    definition,
                    native,
                    NumericalOptions {
                        tolerance_mm: tolerance.clone(),
                        vertex_tolerance_mm: None,
                    },
                    &mut dec,
                )
                .unwrap();
                let chains = super::super::rolling::construct(
                    &model,
                    body,
                    &[edge],
                    &radius,
                    &tolerance,
                    &mut dec,
                )
                .unwrap();
                cap_strip(&source, &chains[0].strips[0], "round", &tolerance, &mut dec).unwrap()
            });
            assert!(!work.exhausted, "{work:?}");
            println!(
                "size={:?} radius={radius}",
                size.each_ref().map(ToString::to_string)
            );
            super::super::tag_audit::assert_blend_rows(&body, 1);
        }
    }
}
