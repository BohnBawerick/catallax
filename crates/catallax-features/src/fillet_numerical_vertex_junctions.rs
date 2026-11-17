//! Vertex-local shared curves and native charts, before topology publication.
//! Tangent pairs follow SMLib's checked shared-section case. Agreement between
//! fitted carriers is sampled numerical evidence, not exact coincidence.

#[path = "fillet_numerical_vertex_assembly.rs"]
pub(crate) mod assembly;

#[path = "fillet_numerical_vertex_cliff.rs"]
pub(crate) mod cliff;
#[path = "fillet_numerical_vertex_ends.rs"]
pub(crate) mod ends;

use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalSample, Rational, Real, Refusal, Sign,
};
use catallax_shapes::{BSplineCurve, Point3, Surface};
use catallax_topology::{FaceId, Model, Sense, VertexId};

use super::super::rolling::{RolledChain, RolledStrip};
use super::{BITS, StripEnd, boundaries, cross, dot, unit, upper_distance};

#[derive(Debug)]
pub(crate) struct VertexJunction {
    pub source: VertexId,
    pub contacts: Vec<VertexContact>,
    pub boundaries: Vec<JunctionBoundary>,
    pub terminations: Vec<StripTermination>,
    pub patches: Vec<JunctionPatch>,
    pub checks: Vec<NumericalCheck>,
}

#[derive(Debug)]
pub(crate) struct VertexContact {
    pub position: [Rational; 3],
    /// Keep each rail's native branch, including distinct periodic lifts on
    /// the same face. A termination identifies its own entries below.
    pub parameters: Vec<(FaceId, [Rational; 2])>,
}

#[derive(Debug)]
pub(crate) struct JunctionBoundary {
    pub curve: BSplineCurve,
    pub contacts: [usize; 2],
    /// Both pcurves use the shared world curve's parameter 0..1.
    pub sides: [BoundarySide; 2],
}

#[derive(Debug)]
pub(crate) struct BoundarySide {
    pub face: JunctionFace,
    pub parameters: BSplineCurve,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum JunctionFace {
    Strip(StripEnd),
    Source(FaceId),
    Patch(usize),
}

#[derive(Debug)]
pub(crate) struct StripTermination {
    pub strip: StripEnd,
    /// Longitudinal parameters and contact indices in strip rail order.
    pub tau: [Rational; 2],
    pub contacts: [usize; 2],
    /// Per-rail index into that contact's native-parameter entries.
    pub support_parameters: [usize; 2],
    /// Ordered from rail 0 to rail 1, not a material face-loop orientation.
    pub boundaries: Vec<(usize, Sense)>,
}

#[derive(Debug)]
pub(crate) struct JunctionPatch {
    pub surface: Surface,
    pub sense: Sense,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-vertex-junction", detail)
}

fn strip<'a>(
    model: &Model,
    chains: &'a [RolledChain],
    vertex: VertexId,
    end: StripEnd,
) -> Result<&'a RolledStrip, Refusal> {
    let source = chains
        .get(end.chain)
        .and_then(|chain| chain.strips.get(end.strip))
        .ok_or_else(|| bad("junction end does not identify a rolling strip"))?;
    let edge = model.edge(source.source.edge);
    if [edge.start, edge.end].get(end.source_end) != Some(&Some(vertex)) {
        return Err(bad("junction end does not meet its source vertex"));
    }
    Ok(source)
}

fn length(vector: &[Ball; 3]) -> Result<Rational, Refusal> {
    vector
        .iter()
        .fold(
            Ball::from_rational(&Rational::zero(), BITS),
            |sum, value| &sum + &value.square(),
        )
        .sqrt()
        .map(|value| value.hi())
        .ok_or_else(|| bad("junction direction or distance cannot be enclosed"))
}

fn gap(a: &[Ball; 3], b: &[Ball; 3]) -> Result<Rational, Refusal> {
    length(&std::array::from_fn(|i| &a[i] - &b[i]))
}

fn source_geometry(
    model: &Model,
    source: &RolledStrip,
    vertex: VertexId,
    end: StripEnd,
    dec: &mut Decider,
) -> Result<([Ball; 3], [[Ball; 3]; 2]), Refusal> {
    let edge = model.edge(source.source.edge);
    let point = model.vertex(vertex).point;
    let positive = (end.source_end == 0) == (edge.sense == Sense::Same);
    let tangents = catallax_intersect::tangent_at(model.geom(), edge.curve, point, dec)?;
    let tangent = &tangents[usize::from(positive)];
    let tangent = unit([
        dec.enclose_budgeted_at(&tangent.x, BITS)?,
        dec.enclose_budgeted_at(&tangent.y, BITS)?,
        dec.enclose_budgeted_at(&tangent.z, BITS)?,
    ])?;
    let tangent = if positive {
        tangent
    } else {
        tangent.map(|v| -v)
    };
    let point = model.geom().position(point);
    let normal = |side: usize, dec: &mut Decider| -> Result<[Ball; 3], Refusal> {
        let face = model.face(source.source.faces[side]);
        let n = model.geom().surface(face.surface).normal_at(&point, dec)?;
        let n = unit([
            dec.enclose_budgeted_at(&n.x, BITS)?,
            dec.enclose_budgeted_at(&n.y, BITS)?,
            dec.enclose_budgeted_at(&n.z, BITS)?,
        ])?;
        Ok(if face.sense == Sense::Same {
            n
        } else {
            n.map(|v| -v)
        })
    };
    Ok((tangent, [normal(0, dec)?, normal(1, dec)?]))
}

fn section_parameters(
    source: &RolledStrip,
    tau: &Rational,
    reversed: bool,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let mut domain = source.fillet.surface().source().domain()[1];
    if reversed {
        domain.reverse();
    }
    BSplineCurve::new(
        1,
        [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        domain
            .map(|v| Point3::new(tau.clone().into(), v.clone().into(), Real::zero()))
            .to_vec(),
        vec![Real::one(); 2],
        dec,
    )
}

/// Retains one shared section of two source-tangent strips. Sharp corners must
/// use their own construction, including unequal reach and exposed neighbors.
pub(crate) fn tangent_pair(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    ends: [StripEnd; 2],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexJunction, Refusal> {
    if radius.sign() != Sign::Positive || tolerance.sign() != Sign::Positive {
        return Err(bad("junction radius and tolerance must be positive"));
    }
    let sources = [
        strip(model, chains, vertex, ends[0])?,
        strip(model, chains, vertex, ends[1])?,
    ];
    if sources[0].source.edge == sources[1].source.edge || sources[0].convex != sources[1].convex {
        return Err(bad(
            "a tangent pair needs distinct strips of the same convexity",
        ));
    }
    let common = sources[0]
        .source
        .faces
        .iter()
        .enumerate()
        .find_map(|(i, face)| {
            sources[1]
                .source
                .faces
                .iter()
                .position(|other| other == face)
                .map(|j| (i, j))
        })
        .ok_or_else(|| bad("tangent strips have no common original support"))?;
    let swap = usize::from(common.0 != common.1);
    let geometry = [
        source_geometry(model, sources[0], vertex, ends[0], dec)?,
        source_geometry(model, sources[1], vertex, ends[1], dec)?,
    ];
    let mut source_error = &length(&std::array::from_fn(|i| {
        &geometry[0].0[i] + &geometry[1].0[i]
    }))? * radius;
    for side in 0..2 {
        source_error =
            source_error.max(&gap(&geometry[0].1[side], &geometry[1].1[side ^ swap])? * radius);
    }
    if source_error > *tolerance {
        return Err(bad(
            "source edges or support normals do not form a tangent pair at the requested radius",
        ));
    }
    let tau = ends.map(|end| Rational::from_i64(end.source_end as i64));
    let surfaces = sources.map(|source| source.fillet.surface().source());
    let curve = boundaries::isocurve(surfaces[0], &tau[0], dec)?;
    let parameters = [
        section_parameters(sources[0], &tau[0], false, dec)?,
        section_parameters(sources[1], &tau[1], swap != 0, dec)?,
    ];
    let centre = sources[0]
        .fillet
        .spine()
        .centre_curve()
        .evaluate_box([&tau[0], &tau[0]], BITS, dec)?
        .map(|value| value.mid());
    let other_centre =
        sources[1]
            .fillet
            .spine()
            .centre_curve()
            .evaluate_box([&tau[1], &tau[1]], BITS, dec)?;
    let centre_error = upper_distance(&other_centre, &centre);
    if centre_error > *tolerance {
        return Err(bad(
            "fitted tangent strips have different ball centers at their join",
        ));
    }
    let mut contacts = Vec::with_capacity(2);
    let mut support_parameters = [[0; 2]; 2];
    let mut endpoint_error = [Rational::zero(), Rational::zero()];
    for contact in 0..2 {
        let s = Rational::from_i64(contact as i64);
        let position = curve
            .evaluate_box([&s, &s], BITS, dec)?
            .map(|value| value.mid());
        let mut native = Vec::with_capacity(2);
        for strip_index in 0..2 {
            let side = contact ^ if strip_index == 0 { 0 } else { swap };
            let source = sources[strip_index];
            let t = &tau[strip_index];
            let rail = source
                .fillet
                .rail(side)
                .ok_or_else(|| bad("tangent contact has no world rail"))?;
            endpoint_error[contact] = endpoint_error[contact].clone().max(upper_distance(
                &rail.evaluate_box([t, t], BITS, dec)?,
                &position,
            ));
            let uv = source
                .fillet
                .parameter_rail(side)
                .ok_or_else(|| bad("tangent contact has no native support rail"))?
                .evaluate_box([t, t], BITS, dec)?;
            let uv = [uv[0].mid(), uv[1].mid()];
            let enclosed = uv.each_ref().map(|value| Ball::from_rational(value, BITS));
            let face = source.source.faces[side];
            let world = model.geom().surface(model.face(face).surface).evaluate(
                &enclosed[0],
                &enclosed[1],
                dec,
            )?;
            endpoint_error[contact] = endpoint_error[contact]
                .clone()
                .max(upper_distance(&world, &position));
            support_parameters[strip_index][side] = native.len();
            native.push((face, uv));
        }
        if endpoint_error[contact] > *tolerance {
            return Err(bad(
                "tangent contact misses an original support or its fitted rail",
            ));
        }
        contacts.push(VertexContact {
            position,
            parameters: native,
        });
    }
    let mut samples = Vec::new();
    for span in curve.knots().windows(2).filter(|span| span[0] < span[1]) {
        for index in 0..=8 {
            dec.charge_work(1)?;
            let s = &span[0] + &((&span[1] - &span[0]) * Rational::from_ratio(index, 8).unwrap());
            let world = curve.evaluate_box([&s, &s], BITS, dec)?;
            let mut error = centre_error.clone();
            let radial: [Ball; 3] =
                std::array::from_fn(|i| &world[i] - &Ball::from_rational(&centre[i], BITS));
            let distance = dot(&radial, &radial)
                .sqrt()
                .ok_or_else(|| bad("tangent section radius cannot be enclosed"))?;
            error = error.max(
                (&distance.lo() - radius)
                    .abs()
                    .max((&distance.hi() - radius).abs()),
            );
            let normal =
                |side: usize, dec: &mut Decider| -> Result<([Ball; 3], Rational), Refusal> {
                    let uv = parameters[side].evaluate_box([&s, &s], BITS, dec)?;
                    let point = surfaces[side].evaluate_box(
                        [&uv[0].lo(), &uv[0].hi()],
                        [&uv[1].lo(), &uv[1].hi()],
                        BITS,
                        dec,
                    )?;
                    let error = gap(&world, &point)?;
                    let [du, dv] = surfaces[side].derivatives(&uv[0], &uv[1], dec)?;
                    Ok((unit(cross(&du, &dv))?, error))
                };
            let (first, first_error) = normal(0, dec)?;
            let (second, second_error) = normal(1, dec)?;
            let mut normals = [first, second];
            error = error.max(first_error).max(second_error);
            let alignment = dot(&normals[0], &normals[1]);
            if alignment.hi().sign() == Sign::Negative {
                for value in &mut normals[1] {
                    *value = -&*value;
                }
            } else if alignment.lo().sign() != Sign::Positive {
                return Err(bad("tangent section normal orientation is unresolved"));
            }
            error = error.max(&gap(&normals[0], &normals[1])? * radius);
            if s.is_zero() {
                error = error.max(endpoint_error[0].clone());
            }
            if s == Rational::one() {
                error = error.max(endpoint_error[1].clone());
            }
            if error > *tolerance {
                return Err(bad(
                    "fitted tangent sections miss their shared curve, radius or tangent tolerance",
                ));
            }
            samples.push(NumericalSample {
                parameters: vec![s],
                distance_upper_mm: error,
            });
        }
    }
    let [first, second] = parameters;
    Ok(VertexJunction {
        source: vertex,
        contacts,
        boundaries: vec![JunctionBoundary {
            curve,
            contacts: [0, 1],
            sides: [
                BoundarySide { face: JunctionFace::Strip(ends[0]), parameters: first },
                BoundarySide { face: JunctionFace::Strip(ends[1]), parameters: second },
            ],
        }],
        terminations: ends.into_iter().enumerate().map(|(index, end)| StripTermination {
            strip: end,
            tau: [tau[index].clone(), tau[index].clone()],
            contacts: if index == 0 || swap == 0 { [0, 1] } else { [1, 0] },
            support_parameters: support_parameters[index],
            boundaries: vec![(0, if index == 0 || swap == 0 { Sense::Same } else { Sense::Reversed })],
        }).collect(),
        patches: Vec::new(),
        checks: vec![
            NumericalCheck::Samples {
                name: "fillet.vertex.source-tangent-pair-at-radius".into(),
                scope: format!("{vertex:?}; enclosed source endpoint directions scaled by radius"),
                samples: vec![NumericalSample { parameters: vec![radius.clone()], distance_upper_mm: source_error }],
            },
            NumericalCheck::Samples {
                name: "fillet.vertex.shared-tangent-section".into(),
                scope: "world and paired native curves, radius and tangent gaps; finite sample locations only".into(),
                samples,
            },
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, current_operation_work, with_operation_budget};
    use catallax_recipe::{Recipe, Role};

    #[test]
    fn owner_plane_cylinder_join_retains_one_world_curve_and_both_native_charts() {
        let ((rolling_work, junction_work), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let recipe =
                Recipe::parse(include_str!("../tests/fixtures/owner-curved-pocket.json")).unwrap();
            let built = crate::evaluate(&recipe, &mut dec).unwrap();
            let model = &built.model;
            let selected: Vec<_> = model.edges_of_body(built.result).into_iter().filter(|edge| {
                model.edge_name(built.result, *edge).is_some_and(|name| {
                    name.between().iter().any(|origin| origin.step.as_str() == "blank")
                        && name.between().iter().any(|origin| {
                            origin.step.as_str() == "pocket-tool"
                                && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "far")
                        })
                })
            }).collect();
            assert_eq!(selected.len(), 2);
            let radius = Rational::one();
            let tolerance = Rational::parse("0.0000001").unwrap();
            let before = current_operation_work().used;
            let chains = super::super::super::rolling::construct(
                model,
                built.result,
                &selected,
                &radius,
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let after_rolling = current_operation_work().used;
            assert_eq!(chains.len(), 1);
            assert_eq!(chains[0].strips.len(), 2);
            let first = model.edge(chains[0].strips[0].source.edge);
            let second = model.edge(chains[0].strips[1].source.edge);
            let vertex = [first.start, first.end]
                .into_iter()
                .flatten()
                .find(|vertex| second.start == Some(*vertex) || second.end == Some(*vertex))
                .unwrap();
            let ends: [StripEnd; 2] = std::array::from_fn(|index| StripEnd {
                chain: 0,
                strip: index,
                source_end: usize::from([first, second][index].start != Some(vertex)),
            });
            let mut junction_work = 0;
            for order in [ends, [ends[1], ends[0]]] {
                let before = current_operation_work().used;
                let junction =
                    tangent_pair(model, &chains, vertex, order, &radius, &tolerance, &mut dec)
                        .unwrap();
                junction_work += current_operation_work().used - before;
                assert_eq!(junction.source, vertex);
                assert!(
                    junction.patches.is_empty(),
                    "a tangent seam creates no corner face"
                );
                assert_eq!(junction.boundaries.len(), 1);
                let boundary = &junction.boundaries[0];
                for index in 0..=31 {
                    let s = Rational::from_ratio(index, 31).unwrap();
                    let world = boundary.curve.evaluate_box([&s, &s], BITS, &dec).unwrap();
                    for side in &boundary.sides {
                        let JunctionFace::Strip(end) = side.face else {
                            panic!("tangent seam must be shared by two strips")
                        };
                        let uv = side.parameters.evaluate_box([&s, &s], BITS, &dec).unwrap();
                        let surface = chains[end.chain].strips[end.strip]
                            .fillet
                            .surface()
                            .source();
                        let composed = surface
                            .evaluate_box(
                                [&uv[0].lo(), &uv[0].hi()],
                                [&uv[1].lo(), &uv[1].hi()],
                                BITS,
                                &dec,
                            )
                            .unwrap();
                        assert!(
                            gap(&world, &composed).unwrap() <= tolerance,
                            "native seam composition at {s}"
                        );
                    }
                }
                for termination in &junction.terminations {
                    let source = &chains[termination.strip.chain].strips[termination.strip.strip];
                    for side in 0..2 {
                        let contact = &junction.contacts[termination.contacts[side]];
                        let (face, uv) = &contact.parameters[termination.support_parameters[side]];
                        assert_eq!(*face, source.source.faces[side]);
                        let uv = uv.each_ref().map(|value| Ball::from_rational(value, BITS));
                        let support = model
                            .geom()
                            .surface(model.face(*face).surface)
                            .evaluate(&uv[0], &uv[1], &mut dec)
                            .unwrap();
                        assert!(upper_distance(&support, &contact.position) <= tolerance);
                        let t = &termination.tau[side];
                        let rail = source
                            .fillet
                            .rail(side)
                            .unwrap()
                            .evaluate_box([t, t], BITS, &dec)
                            .unwrap();
                        assert!(upper_distance(&rail, &contact.position) <= tolerance);
                    }
                }
            }
            (after_rolling - before, junction_work)
        });
        assert!(!work.exhausted);
        println!(
            "owner surface-tangent junction rolling_work={rolling_work} junction_work={junction_work} total_work={}; geometry only, no body publication",
            work.used
        );
    }
}

fn tangent_neighbors(chains: &[RolledChain], a: StripEnd, b: StripEnd) -> bool {
    if a.chain != b.chain || a.strip == b.strip {
        return false;
    }
    let chain = &chains[a.chain];
    let entry =
        |end: StripEnd| usize::from(chain.strips[end.strip].source.sense == Sense::Reversed);
    let follows = |from: StripEnd, to: StripEnd| {
        (from.strip + 1 == to.strip
            || chain.closed && from.strip + 1 == chain.strips.len() && to.strip == 0)
            && from.source_end == 1 - entry(from)
            && to.source_end == entry(to)
    };
    follows(a, b) || follows(b, a)
}

/// Construct candidates at finite original vertices. Closed one-edge contours
/// have no vertex fan here; their retained seam belongs to contour construction.
pub(crate) fn construct(
    model: &Model,
    body: catallax_topology::BodyId,
    chains: &[RolledChain],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Vec<VertexJunction>, Refusal> {
    if radius.sign() != Sign::Positive || tolerance.sign() != Sign::Positive {
        return Err(bad("vertex radius and tolerance must be positive"));
    }
    let count: usize = chains.iter().map(|chain| chain.strips.len()).sum();
    let mut ends = Vec::with_capacity(count.saturating_mul(2));
    for (chain, group) in chains.iter().enumerate() {
        for (index, source) in group.strips.iter().enumerate() {
            let edge = model.edge(source.source.edge);
            if edge.start == edge.end {
                continue;
            }
            for (source_end, vertex) in [edge.start, edge.end].into_iter().enumerate() {
                if let Some(vertex) = vertex {
                    ends.push((
                        vertex,
                        StripEnd {
                            chain,
                            strip: index,
                            source_end,
                        },
                    ));
                }
            }
        }
    }
    ends.sort_unstable_by_key(|&(vertex, end)| (vertex, end.chain, end.strip, end.source_end));
    let groups = ends.chunk_by(|a, b| a.0 == b.0).count();
    let mut result = Vec::with_capacity(groups);
    for group in ends.chunk_by(|a, b| a.0 == b.0) {
        dec.charge_work(1)?;
        let vertex = group[0].0;
        let candidate = match group {
            [(_, end)] => ends::face_end(model, body, chains, vertex, *end, tolerance, dec)?,
            [(_, first), (_, second)] => {
                let pair = [*first, *second];
                if tangent_neighbors(chains, pair[0], pair[1]) {
                    tangent_pair(model, chains, vertex, pair, radius, tolerance, dec)?
                } else {
                    let tau = cliff::common_contact(model, chains, vertex, pair, tolerance, dec)?;
                    let section = cliff::shared_section(
                        model, body, chains, vertex, pair, tau, tolerance, dec,
                    )?;
                    assembly::cliff(model, chains, vertex, section, tolerance, dec)?
                }
            }
            [(_, first), (_, second), (_, third)] => assembly::ball(
                model,
                chains,
                vertex,
                &[*first, *second, *third],
                radius,
                tolerance,
                dec,
            )?,
            _ => {
                return Err(bad(
                    "this vertex producer needs one, two or three incident selected strip ends",
                ));
            }
        };
        result.push(candidate);
    }
    Ok(result)
}
