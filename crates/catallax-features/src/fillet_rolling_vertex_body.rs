//! Shared topology for rolling strips whose ends meet three-support vertices.

use super::body::{material_sense, occurrence, publish};
use super::rolling::RolledChain;
use super::vertex::{VertexBlend, boundaries};
use super::*;
use catallax_topology::{FaceId, Model};

pub(crate) fn build(
    source: &NumericalBody,
    model: &Model,
    chains: &[RolledChain],
    corners: &[VertexBlend],
    key: &str,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    let original = source.definition();
    let mut offsets = Vec::with_capacity(chains.len() + 1);
    offsets.push(0usize);
    for chain in chains {
        offsets.push(offsets.last().unwrap() + chain.strips.len());
    }
    let strip_count = *offsets.last().unwrap();
    let corner_face_start = original.faces.len() + strip_count;
    let face_keys: BTreeMap<_, _> = original
        .faces
        .iter()
        .enumerate()
        .map(|(i, face)| (face.key.as_str(), i))
        .collect();
    let edge_keys: BTreeMap<_, _> = original
        .edges
        .iter()
        .enumerate()
        .map(|(i, edge)| (edge.key.as_str(), i))
        .collect();
    let mut face_indices = BTreeMap::<FaceId, usize>::new();
    for corner in corners {
        for face in corner.faces {
            if let std::collections::btree_map::Entry::Vacant(entry) = face_indices.entry(face) {
                let name = format!("{face:?}");
                let index = face_keys
                    .get(name.as_str())
                    .copied()
                    .ok_or_else(|| bad("a vertex support has no retained source face"))?;
                entry.insert(index);
            }
        }
    }
    let mut end_sections = BTreeMap::new();
    let mut geometry = Vec::with_capacity(corners.len());
    for (corner_index, vertex) in corners.iter().enumerate() {
        let section = |index: usize, dec: &mut Decider| {
            let section = &vertex.sections[index];
            boundaries::construct(
                &chains[section.strip.chain],
                vertex,
                section,
                tolerance,
                dec,
            )
        };
        geometry.push([section(0, dec)?, section(1, dec)?, section(2, dec)?]);
        for (index, section) in vertex.sections.iter().enumerate() {
            let address = (
                section.strip.chain,
                section.strip.strip,
                section.strip.source_end,
            );
            if end_sections
                .insert(address, (corner_index, index))
                .is_some()
            {
                return Err(bad("a rolling end has more than one shared vertex section"));
            }
        }
    }
    for (chain_index, chain) in chains.iter().enumerate() {
        for (strip_index, strip) in chain.strips.iter().enumerate() {
            let edge = model.edge(strip.source.edge);
            for (end, vertex) in [edge.start, edge.end].into_iter().enumerate() {
                let Some(vertex) = vertex else {
                    return Err(bad("a vertex network contains an unbounded or ring end"));
                };
                let (corner, _) = end_sections
                    .get(&(chain_index, strip_index, end))
                    .copied()
                    .ok_or_else(|| {
                        bad(format!(
                            "rolling end {vertex:?} does not have a three-support vertex section"
                        ))
                    })?;
                if corners[corner].source != vertex {
                    return Err(bad("a shared section names a different source vertex"));
                }
            }
        }
    }
    let mut d = original.clone();
    // As in `fillet_numerical::build_impl`: earlier derived copies are
    // measured again.
    d.reset_copy_tags();
    d.name = format!("{} / {key}", original.name);
    let mut native = native_rows(source)?;
    let mut parameters: Parameters = BTreeMap::new();
    let mut vertices = Vec::with_capacity(corners.len());
    for (corner_index, corner) in corners.iter().enumerate() {
        let mut canonical: [Option<&boundaries::BoundaryContact>; 3] = [None; 3];
        for boundary in &geometry[corner_index] {
            for contact in &boundary.contacts {
                if let Some(previous) = canonical[contact.contact]
                    && (previous.xyz != contact.xyz || previous.sphere_uv != contact.sphere_uv)
                {
                    return Err(bad(
                        "shared sections disagree on a canonical vertex contact",
                    ));
                }
                canonical[contact.contact] = Some(contact);
            }
        }
        let mut ids = [0; 3];
        for contact in 0..3 {
            let canonical = canonical[contact]
                .ok_or_else(|| bad("a vertex contact has no section boundary"))?;
            let id = d.vertices.len();
            ids[contact] = id;
            d.vertices.push(NumericalVertex {
                key: format!(
                    "{key}/vertex/{:?}/contact/{:?}",
                    corner.source, corner.faces[contact]
                ),
                position: canonical.xyz.clone(),
            });
            parameters.insert(
                (face_indices[&corner.faces[contact]], id),
                corner.ball.parameters()[contact].clone(),
            );
            parameters.insert(
                (corner_face_start + corner_index, id),
                canonical.sphere_uv.clone(),
            );
        }
        vertices.push(ids);
    }
    let mut section_edges = Vec::with_capacity(corners.len());
    for (corner_index, corner) in corners.iter().enumerate() {
        let mut edges = [0; 3];
        for (index, boundary) in geometry[corner_index].iter().enumerate() {
            let section = &corner.sections[index];
            let strip = &chains[section.strip.chain].strips[section.strip.strip];
            edges[index] = d.edges.len();
            d.edges.push(NumericalEdge {
                key: format!(
                    "{key}/vertex/{:?}/section/{:?}",
                    corner.source, strip.source.edge
                ),
                curve: Curve::BSpline(boundary.curve.clone()),
                vertices: boundary
                    .contacts
                    .each_ref()
                    .map(|contact| vertices[corner_index][contact.contact]),
                parameters: [q(0), q(1)],
            });
        }
        section_edges.push(edges);
    }
    let mut sampled = source.sampled_carriers().clone();
    for (chain_index, chain) in chains.iter().enumerate() {
        for (strip_index, strip) in chain.strips.iter().enumerate() {
            let edge_key = format!("{:?}", strip.source.edge);
            let selected = edge_keys
                .get(edge_key.as_str())
                .copied()
                .ok_or_else(|| bad("a rolling network edge has no retained source occurrence"))?;
            let sides = strip.source.faces.map(|face| {
                face_indices
                    .get(&face)
                    .copied()
                    .ok_or_else(|| bad("a rolling support is absent from its shared vertex fan"))
            });
            let [a, b] = sides;
            let sides = [a?, b?];
            let occurrences = [
                occurrence(original, sides[0], selected)?,
                occurrence(original, sides[1], selected)?,
            ];
            let senses = std::array::from_fn::<_, 2, _>(|side| {
                original.faces[sides[side]].loops[occurrences[side].0][occurrences[side].1].sense
            });
            let ends = [
                end_sections[&(chain_index, strip_index, 0)],
                end_sections[&(chain_index, strip_index, 1)],
            ];
            let tau = ends.map(|(corner, section)| corners[corner].sections[section].tau.clone());
            if tau[0] >= tau[1] {
                return Err(bad(
                    "vertex contacts consume or reverse a rolling strip interval",
                ));
            }
            let face = original.faces.len() + offsets[chain_index] + strip_index;
            let mut rails = [0; 2];
            for side in 0..2 {
                let mut ids = [0; 2];
                for (end, &(corner, section)) in ends.iter().enumerate() {
                    let contact = geometry[corner][section].contacts[side].contact;
                    if corners[corner].faces[contact] != strip.source.faces[side] {
                        return Err(bad(
                            "a vertex section contact names a different rolling support",
                        ));
                    }
                    ids[end] = vertices[corner][contact];
                    parameters.insert((face, ids[end]), [tau[end].clone(), q(side as i64)]);
                }
                let uv = ids.map(|vertex| parameters[&(sides[side], vertex)].clone());
                let pcurve = strip
                    .fillet
                    .parameter_rail(side)
                    .ok_or_else(|| bad("a rolling contact pcurve is missing"))?;
                let (pcurve, _) = constrained_rail(pcurve, &tau, &uv, None, dec)?;
                let (boundary, use_index) = occurrences[side];
                native[sides[side]][boundary][use_index] = Some(pcurve);
                rails[side] = d.edges.len();
                d.edges.push(NumericalEdge {
                    key: format!("{key}/{edge_key}/rail/{}", original.faces[sides[side]].key),
                    curve: Curve::BSpline(
                        strip
                            .fillet
                            .rail(side)
                            .ok_or_else(|| bad("a rolling rail is missing"))?
                            .clone(),
                    ),
                    vertices: ids,
                    parameters: tau.clone(),
                });
                let uv = walked(&d.edges[rails[side]], senses[side])
                    .map(|vertex| parameters[&(sides[side], vertex)].clone());
                d.faces[sides[side]].loops[boundary][use_index] =
                    new_trim(rails[side], senses[side], uv);
            }
            let cuts = ends.map(|(corner, section)| section_edges[corner][section]);
            let first = (rails[0], flip(senses[0]));
            let reaches = usize::from(first.1 == Sense::Same);
            let path = [
                first,
                (cuts[reaches], Sense::Same),
                (rails[1], flip(senses[1])),
                (cuts[1 - reaches], Sense::Reversed),
            ];
            let mut loop_ = Vec::with_capacity(4);
            let mut curves = Vec::with_capacity(4);
            for (edge, sense) in path {
                let pcurve = if edge == rails[0] {
                    parameter_line(strip.fillet.rail(0).unwrap().domain(), 0, dec)?
                } else if edge == rails[1] {
                    parameter_line(strip.fillet.rail(1).unwrap().domain(), 1, dec)?
                } else {
                    let (corner, section) = ends[usize::from(edge == cuts[1])];
                    geometry[corner][section].parameters[0].clone()
                };
                curves.push(Some(pcurve));
                let mut trim = new_trim(
                    edge,
                    sense,
                    walked(&d.edges[edge], sense).map(|vertex| parameters[&(face, vertex)].clone()),
                );
                if rails.contains(&edge) {
                    trim.tag = blend_row();
                }
                loop_.push(trim);
            }
            let face_key = format!("{key}/blend/{edge_key}");
            sampled.insert(face_key.clone());
            d.faces.push(NumericalFace {
                key: face_key,
                surface: Surface::BSpline(SplineSupport::RationalAtlas(
                    strip.fillet.surface().clone(),
                )),
                sense: material_sense(strip, &(&tau[0] + &tau[1]).mul_pow2(-1), dec)?,
                loops: vec![loop_],
            });
            native.push(vec![curves]);
        }
    }
    for (corner_index, corner) in corners.iter().enumerate() {
        let face = corner_face_start + corner_index;
        let mut loop_ = Vec::with_capacity(3);
        let mut curves = Vec::with_capacity(3);
        for (index, geometry) in geometry[corner_index].iter().enumerate() {
            let edge = section_edges[corner_index][index];
            let sense = if geometry.reversed {
                Sense::Reversed
            } else {
                Sense::Same
            };
            loop_.push(new_trim(
                edge,
                sense,
                walked(&d.edges[edge], sense).map(|vertex| parameters[&(face, vertex)].clone()),
            ));
            curves.push(Some(geometry.parameters[1].clone()));
        }
        d.faces.push(NumericalFace {
            key: format!("{key}/corner/{:?}", corner.source),
            surface: corner.surface.clone(),
            sense: corner.sense,
            loops: vec![loop_],
        });
        native.push(vec![curves]);
    }
    publish(source, d, native, sampled, dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};
    use catallax_shapes::Pose;
    use catallax_topology::numerical::{NumericalOptions, build_with_native_trims};

    #[test]
    fn twelve_rolling_edges_publish_one_body_with_the_trimmed_source_domains() {
        let (body, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let mut model = Model::default();
            let body = catallax_topology::block(&mut model, [20, 30, 40].map(q), &Pose::identity())
                .unwrap();
            let radius = q(1);
            let tolerance = Rational::parse("0.0000001").unwrap();
            let selected = model.edges_of_body(body);
            let chains = super::super::rolling::construct(
                &model, body, &selected, &radius, &tolerance, &mut dec,
            )
            .unwrap();
            let corners =
                super::super::vertex::construct(&model, &chains, &radius, &tolerance, &mut dec)
                    .unwrap();
            let (definition, native) = crate::numerical_embedding::definition_from_exact_model(
                &model,
                body,
                "rolling vertex network".into(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let source = build_with_native_trims(
                definition,
                native,
                NumericalOptions {
                    tolerance_mm: tolerance.clone(),
                    vertex_tolerance_mm: None,
                },
                &mut dec,
            )
            .unwrap();
            build(
                &source, &model, &chains, &corners, "round", &tolerance, &mut dec,
            )
            .unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        for (axis, coordinate, expected) in [
            (0, 0, 1064),
            (0, 20, 1064),
            (1, 0, 684),
            (1, 30, 684),
            (2, 0, 504),
            (2, 40, 504),
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
                        && origin[axis].as_rational() == Some(q(coordinate))
                })
                .unwrap();
            let twice_area = face.loops.iter().flatten().fold(q(0), |sum, trim| {
                sum + &trim.uv[0][0] * &trim.uv[1][1] - &trim.uv[1][0] * &trim.uv[0][1]
            });
            assert!(
                (twice_area.abs().mul_pow2(-1) - q(expected)).abs()
                    <= Rational::parse("0.00001").unwrap(),
                "axis={axis}, coordinate={coordinate}: trimmed original face area is wrong"
            );
        }
        println!(
            "twelve rolling edges, shared corner boundaries and closed-body publication work={}",
            work.used
        );
    }

    /// Property test of the producer of blend-row tags on a rolled network
    /// (decision Q8, section 4.4): a block of random size with all twelve
    /// edges rounded at a random radius, and the full check on the two rails
    /// of each of the twelve blend faces.
    #[test]
    fn the_rails_of_a_rolled_network_are_rows_of_their_blend_faces() {
        let mut random = super::super::tag_audit::Stream::seeded("rolled network");
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
            let tolerance = Rational::parse("0.0000001").unwrap();
            let selected = model.edges_of_body(body);
            let chains = super::super::rolling::construct(
                &model, body, &selected, &radius, &tolerance, &mut dec,
            )
            .unwrap();
            let corners =
                super::super::vertex::construct(&model, &chains, &radius, &tolerance, &mut dec)
                    .unwrap();
            let (definition, native) = crate::numerical_embedding::definition_from_exact_model(
                &model,
                body,
                "rolling vertex network".into(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let source = build_with_native_trims(
                definition,
                native,
                NumericalOptions {
                    tolerance_mm: tolerance.clone(),
                    vertex_tolerance_mm: None,
                },
                &mut dec,
            )
            .unwrap();
            build(
                &source, &model, &chains, &corners, "round", &tolerance, &mut dec,
            )
            .unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        println!(
            "size={:?} radius={radius}",
            size.each_ref().map(ToString::to_string)
        );
        super::super::tag_audit::assert_blend_rows(&body, 12);
    }
}
