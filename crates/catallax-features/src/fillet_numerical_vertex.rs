//! Numerical common-ball closures of three incident rolling strips.
//! The owner publishes these candidates only after actual face-domain clipping.

#[path = "fillet_numerical_vertex_boundaries.rs"]
pub(crate) mod boundaries;

use catallax_numbers::{Ball, Decider, NumericalCheck, NumericalSample, Rational, Refusal, Sign};
use catallax_shapes::{
    NumericalFilletCorner, NumericalFilletSupport, RPoint3, RVector3, Sphere, Surface,
};
use catallax_topology::{FaceId, Model, Sense, VertexId};

use super::rolling::{RolledChain, RolledStrip};

const BITS: u32 = 128;

#[derive(Clone, Copy, Debug)]
pub(crate) struct StripEnd {
    pub chain: usize,
    pub strip: usize,
    pub source_end: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct VertexSection {
    pub strip: StripEnd,
    /// Original fitted-strip U coordinate; section V runs from zero to one.
    pub tau: Rational,
    /// Contact indices in `VertexBlend::ball`, in strip side order.
    pub contacts: [usize; 2],
    /// Whether the corner loop walks the section from V=1 to V=0.
    pub reversed: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct VertexBlend {
    pub source: VertexId,
    pub faces: [FaceId; 3],
    pub ball: NumericalFilletCorner,
    pub surface: Surface,
    pub sense: Sense,
    pub sections: [VertexSection; 3],
    pub checks: Vec<NumericalCheck>,
}

struct FaceSeed {
    face: FaceId,
    parameters: [Rational; 2],
    offset: Rational,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-vertex", detail)
}

fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&a[j] * &b[k]) - &(&a[k] * &b[j])
    })
}

fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, i| {
        &sum + &(&a[i] * &b[i])
    })
}

fn unit(a: [Ball; 3]) -> Result<[Ball; 3], Refusal> {
    let square = dot(&a, &a);
    if square.lo().sign() != Sign::Positive {
        return Err(bad("vertex contact has a singular or unresolved direction"));
    }
    let inverse = square
        .sqrt()
        .and_then(|length| length.recip())
        .ok_or_else(|| bad("vertex contact direction cannot be normalized"))?;
    Ok(a.map(|value| &value * &inverse))
}

fn upper_distance(a: &[Ball; 3], b: &[Rational; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        sum + (&a[i].lo() - &b[i]).abs().max((&a[i].hi() - &b[i]).abs())
    })
}

fn strip<'a>(chains: &'a [RolledChain], end: &StripEnd) -> &'a RolledStrip {
    &chains[end.chain].strips[end.strip]
}

fn face_seed(
    model: &Model,
    strip: &RolledStrip,
    side: usize,
    near: &[Rational; 3],
    radius: &Rational,
    dec: &mut Decider,
) -> Result<FaceSeed, Refusal> {
    let station = strip
        .fillet
        .spine()
        .nodes()
        .iter()
        .min_by_key(|node| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&node.centre[i] - &near[i]).abs()
            })
        })
        .ok_or_else(|| bad("incident strip has no numerical stations"))?;
    let parameters = station.parameters[side].clone();
    let at = parameters
        .each_ref()
        .map(|value| Ball::from_rational(value, BITS));
    let surface = model
        .geom()
        .surface(model.face(strip.source.faces[side]).surface);
    let point = surface.evaluate(&at[0], &at[1], dec)?;
    let [du, dv] = surface.derivatives(&at[0], &at[1], dec)?;
    let normal = unit(cross(&du, &dv))?;
    let delta = std::array::from_fn(|i| &Ball::from_rational(&station.centre[i], BITS) - &point[i]);
    let distance = dot(&delta, &normal);
    let offset = if distance.lo().sign() == Sign::Positive {
        radius.clone()
    } else if distance.hi().sign() == Sign::Negative {
        -radius
    } else {
        return Err(bad(
            "incident strip does not separate its signed support offset",
        ));
    };
    dec.charge_work(32)?;
    Ok(FaceSeed {
        face: strip.source.faces[side],
        parameters,
        offset,
    })
}

fn section_tau(
    strip: &RolledStrip,
    centre: &[Rational; 3],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let spine = strip.fillet.spine();
    let initial = spine
        .nodes()
        .iter()
        .min_by_key(|node| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&node.centre[i] - &centre[i]).abs()
            })
        })
        .ok_or_else(|| bad("incident strip has no numerical stations"))?;
    let curve = spine.centre_curve();
    let domain = curve.domain();
    let bounds = domain.map(Rational::to_f64_nearest);
    let target = centre.each_ref().map(Rational::to_f64_nearest);
    let mut parameter = initial.tau.to_f64_nearest();
    let limit = tolerance.mul_pow2(-1);
    for _ in 0..32 {
        dec.charge_work(1)?;
        let tau = Rational::from_f64_exact(parameter)
            .ok_or_else(|| bad("vertex section parameter is not finite"))?
            .max(domain[0].clone())
            .min(domain[1].clone());
        let jet = curve.jet_box([&tau, &tau], BITS, dec)?;
        if upper_distance(&jet.value, centre) <= limit {
            return Ok(tau);
        }
        let value = jet
            .value
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let tangent = jet
            .first
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let speed: f64 = tangent.iter().map(|value| value * value).sum();
        let step = (0..3)
            .map(|i| (value[i] - target[i]) * tangent[i])
            .sum::<f64>()
            / speed;
        if !step.is_finite() || !speed.is_finite() || speed <= 0.0 {
            return Err(bad(
                "vertex section has an unresolved or singular centre tangent",
            ));
        }
        let next = (parameter - step).clamp(bounds[0], bounds[1]);
        if next.to_bits() == parameter.to_bits() {
            return Err(bad("fitted incident spine misses the common vertex ball"));
        }
        parameter = next;
    }
    Err(bad(
        "vertex section solve did not meet its declared numerical tolerance",
    ))
}

fn carrier(ball: &NumericalFilletCorner, dec: &mut Decider) -> Result<Surface, Refusal> {
    let radial: [[Rational; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| &ball.contacts()[i][j] - &ball.centre()[j])
    });
    // A mixed-sign combination lies outside either side of the contact cone,
    // keeping both definition poles away from a regular three-contact sector.
    let axis: [Rational; 3] = std::array::from_fn(|i| &radial[0][i] - &radial[1][i]);
    let mean: [Rational; 3] =
        std::array::from_fn(|i| &(&radial[0][i] + &radial[1][i]) + &radial[2][i]);
    let square = axis
        .iter()
        .fold(Rational::zero(), |sum, value| sum + value * value);
    let along = (0..3).fold(Rational::zero(), |sum, i| sum + &axis[i] * &mean[i]);
    let reference: [Rational; 3] =
        std::array::from_fn(|i| &(&mean[i] * &square) - &(&axis[i] * &along));
    Ok(Surface::Sphere(Sphere::new(
        RPoint3::new(
            ball.centre()[0].clone(),
            ball.centre()[1].clone(),
            ball.centre()[2].clone(),
        )
        .into(),
        ball.radius().clone().into(),
        RVector3::new(axis[0].clone(), axis[1].clone(), axis[2].clone()).into(),
        RVector3::new(
            reference[0].clone(),
            reference[1].clone(),
            reference[2].clone(),
        )
        .into(),
        dec,
    )?))
}

fn check_section(
    source: &RolledStrip,
    section: &VertexSection,
    ball: &NumericalFilletCorner,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Vec<NumericalSample>, Refusal> {
    let surface = source.fillet.surface().source();
    let tau = Ball::from_rational(&section.tau, BITS);
    let radius = Ball::from_rational(ball.radius(), BITS);
    let mut checks = Vec::with_capacity(9);
    for index in 0..=8 {
        dec.charge_work(1)?;
        let share = Rational::from_ratio(index, 8).unwrap();
        let v = Ball::from_rational(&share, BITS);
        let point =
            surface.evaluate_box([&section.tau, &section.tau], [&share, &share], BITS, dec)?;
        let radial =
            std::array::from_fn(|i| &point[i] - &Ball::from_rational(&ball.centre()[i], BITS));
        let length = dot(&radial, &radial)
            .sqrt()
            .ok_or_else(|| bad("vertex radius cannot be enclosed"))?;
        let mut error = (&length.lo() - ball.radius())
            .abs()
            .max((&length.hi() - ball.radius()).abs());
        let radial = unit(radial)?;
        let [du, dv] = surface.derivatives(&tau, &v, dec)?;
        let mut normal = unit(cross(&du, &dv))?;
        let alignment = dot(&normal, &radial);
        if alignment.hi().sign() == Sign::Negative {
            normal = normal.map(|value| -value);
        } else if alignment.lo().sign() != Sign::Positive {
            return Err(bad(
                "incident strip tangency cannot be oriented against the common ball",
            ));
        }
        let gap: [Ball; 3] = std::array::from_fn(|i| &normal[i] - &radial[i]);
        let gap_square = gap.iter().fold(
            Ball::from_rational(&Rational::zero(), BITS),
            |sum, value| sum + value.square(),
        );
        let angular_gap = gap_square
            .sqrt()
            .ok_or_else(|| bad("vertex tangent gap cannot be enclosed"))?;
        let angular_mm = &angular_gap * &radius;
        error = error.max(angular_mm.hi());
        if index == 0 || index == 8 {
            error = error.max(upper_distance(
                &point,
                &ball.contacts()[section.contacts[usize::from(index != 0)]],
            ));
        }
        if error > *tolerance {
            return Err(bad(
                "incident section misses the vertex ball, source contact or tangent within the declared tolerance",
            ));
        }
        checks.push(NumericalSample {
            parameters: vec![section.tau.clone(), share],
            distance_upper_mm: error,
        });
    }
    Ok(checks)
}

fn corner(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    ends: &[StripEnd],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexBlend, Refusal> {
    let convex = strip(chains, &ends[0]).convex;
    if ends.iter().any(|end| strip(chains, end).convex != convex) {
        return Err(bad(
            "incident strips require incompatible material sides of one rolling ball",
        ));
    }
    let point = model.geom().position(model.vertex(vertex).point);
    let near = [
        dec.enclose_budgeted_at(&point.x, BITS)?.mid(),
        dec.enclose_budgeted_at(&point.y, BITS)?.mid(),
        dec.enclose_budgeted_at(&point.z, BITS)?.mid(),
    ];
    let mut seeds = Vec::<FaceSeed>::with_capacity(3);
    let mut contacts = [[0usize; 2]; 3];
    for (index, end) in ends.iter().enumerate() {
        let source = strip(chains, end);
        for side in 0..2 {
            let face = source.source.faces[side];
            let found = if let Some(found) = seeds.iter().position(|seed| seed.face == face) {
                found
            } else {
                if seeds.len() == 3 {
                    return Err(bad(
                        "three incident strips do not bound one three-support vertex fan",
                    ));
                }
                seeds.push(face_seed(model, source, side, &near, radius, dec)?);
                seeds.len() - 1
            };
            contacts[index][side] = found;
        }
    }
    if seeds.len() != 3
        || (0..3).any(|face| contacts.iter().flatten().filter(|&&i| i == face).count() != 2)
    {
        return Err(bad(
            "incident section contacts do not form one three-sided vertex loop",
        ));
    }
    let faces = std::array::from_fn(|i| seeds[i].face);
    let parameters = std::array::from_fn(|i| seeds[i].parameters.clone());
    let supports: [NumericalFilletSupport; 3] = seeds
        .into_iter()
        .map(|seed| {
            NumericalFilletSupport::new(
                model.geom().surface(model.face(seed.face).surface).clone(),
                seed.offset,
                dec,
            )
        })
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| bad("vertex support inventory changed"))?;
    let ball = NumericalFilletCorner::solve(&supports, parameters, tolerance.clone(), dec)?;
    let mut sections = Vec::with_capacity(3);
    let mut checks = Vec::with_capacity(4);
    checks.push(NumericalCheck::Samples {
        name: "fillet.vertex.common-ball".into(),
        scope: format!("{vertex:?}"),
        samples: vec![NumericalSample {
            parameters: ball.parameters().iter().flatten().cloned().collect(),
            distance_upper_mm: ball.sampled_error_mm().clone(),
        }],
    });
    for (index, end) in ends.iter().enumerate() {
        let source = strip(chains, end);
        let section = VertexSection {
            strip: *end,
            tau: section_tau(source, ball.centre(), tolerance, dec)?,
            contacts: contacts[index],
            reversed: false,
        };
        checks.push(NumericalCheck::Samples {
            name: "fillet.vertex.section-contact-and-tangency".into(),
            scope: format!("{vertex:?}/{:?}", source.source.edge),
            samples: check_section(source, &section, &ball, tolerance, dec)?,
        });
        sections.push(section);
    }
    // Orient the section cycle through its shared support contacts.
    let start = sections[0].contacts[0];
    let mut current = sections[0].contacts[1];
    for index in 1..3 {
        let next = (index..3)
            .find(|&j| sections[j].contacts.contains(&current))
            .ok_or_else(|| bad("vertex sections do not close through their source supports"))?;
        sections.swap(index, next);
        sections[index].reversed = sections[index].contacts[1] == current;
        current = sections[index].contacts[usize::from(!sections[index].reversed)];
    }
    if current != start {
        return Err(bad("vertex section cycle is open"));
    }
    let contact_order = [
        sections[0].contacts[usize::from(sections[0].reversed)],
        sections[1].contacts[usize::from(sections[1].reversed)],
        sections[2].contacts[usize::from(sections[2].reversed)],
    ];
    let rays: [[Rational; 3]; 3] = contact_order
        .map(|contact| std::array::from_fn(|i| &ball.contacts()[contact][i] - &ball.centre()[i]));
    let orientation = (0..3).fold(Rational::zero(), |sum, i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        sum + &rays[0][i] * &(&(&rays[1][j] * &rays[2][k]) - &(&rays[1][k] * &rays[2][j]))
    });
    if orientation.is_zero() {
        return Err(bad(
            "vertex contacts do not orient a finite spherical sector",
        ));
    }
    if (orientation.sign() == Sign::Positive) != convex {
        sections.reverse();
        for section in &mut sections {
            section.reversed = !section.reversed;
        }
    }
    Ok(VertexBlend {
        source: vertex,
        faces,
        surface: carrier(&ball, dec)?,
        sense: if convex { Sense::Same } else { Sense::Reversed },
        ball,
        sections: sections
            .try_into()
            .map_err(|_| bad("vertex section inventory changed"))?,
        checks,
    })
}

/// Builds numerical vertex candidates. Single-strip ends and tangent joins stay
/// with the trim owner; no closed-body admission is claimed by this producer.
pub(crate) fn construct(
    model: &Model,
    chains: &[RolledChain],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Vec<VertexBlend>, Refusal> {
    if radius.sign() != Sign::Positive || tolerance.sign() != Sign::Positive {
        return Err(bad(
            "vertex radius and numerical tolerance must be positive",
        ));
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
    let mut result = Vec::with_capacity(ends.len() / 3);
    for group in ends.chunk_by(|a, b| a.0 == b.0) {
        dec.charge_work(1)?;
        if group.len() == 3 {
            let sections = [group[0].1, group[1].1, group[2].1];
            result.push(corner(
                model, chains, group[0].0, &sections, radius, tolerance, dec,
            )?);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, current_operation_work, with_operation_budget};
    use catallax_shapes::Pose;
    use std::collections::{BTreeMap, BTreeSet};

    #[test]
    fn twelve_box_strips_share_eight_checked_vertex_balls_at_each_radius() {
        for requested in ["0.5", "1", "2", "5"] {
            let ((strip_work, vertex_work), work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 100_000_000,
                    ..Default::default()
                });
                let radius = Rational::parse(requested).unwrap();
                let tolerance = Rational::parse("0.0000001").unwrap();
                let mut model = Model::default();
                let size = [20, 30, 40].map(Rational::from_i64);
                let body =
                    catallax_topology::block(&mut model, size.clone(), &Pose::identity()).unwrap();
                let selected = model.edges_of_body(body);
                let before_strips = current_operation_work().used;
                let chains = super::super::rolling::construct(
                    &model, body, &selected, &radius, &tolerance, &mut dec,
                )
                .unwrap();
                let after_strips = current_operation_work().used;
                let vertices = construct(&model, &chains, &radius, &tolerance, &mut dec).unwrap();
                let after_vertices = current_operation_work().used;
                assert_eq!(vertices.len(), 8);
                assert_eq!(
                    vertices
                        .iter()
                        .map(|vertex| vertex.source)
                        .collect::<BTreeSet<_>>(),
                    model.vertices_of_body(body).into_iter().collect()
                );
                let mut ends = BTreeMap::<_, BTreeSet<_>>::new();
                for vertex in &vertices {
                    let source = model.geom().position(model.vertex(vertex.source).point);
                    let source = [&source.x, &source.y, &source.z]
                        .map(|value| dec.enclose_budgeted_at(value, BITS).unwrap().mid());
                    for axis in 0..3 {
                        let expected = if source[axis].is_zero() {
                            radius.clone()
                        } else {
                            &size[axis] - &radius
                        };
                        assert!(
                            (&vertex.ball.centre()[axis] - &expected).abs() <= tolerance,
                            "radius={requested}, vertex={:?}, axis={axis}",
                            vertex.source
                        );
                    }
                    let Surface::Sphere(surface) = &vertex.surface else {
                        panic!("common-ball closure must retain its native sphere");
                    };
                    assert_eq!(surface.radius.rational_ref(), Some(&radius));
                    assert_eq!(vertex.sense, Sense::Same);
                    assert_eq!(vertex.faces.into_iter().collect::<BTreeSet<_>>().len(), 3);
                    for section in &vertex.sections {
                        let source = strip(&chains, &section.strip);
                        let curve = source.fillet.spine().centre_curve();
                        let centre = curve
                            .evaluate_box([&section.tau, &section.tau], BITS, &dec)
                            .unwrap();
                        assert!(upper_distance(&centre, vertex.ball.centre()) <= tolerance);
                        ends.entry(source.source.edge)
                            .or_default()
                            .insert(section.strip.source_end);
                    }
                    for check in &vertex.checks {
                        if let NumericalCheck::Samples { samples, .. } = check {
                            assert!(
                                samples
                                    .iter()
                                    .all(|sample| sample.distance_upper_mm <= tolerance)
                            );
                        }
                    }
                }
                assert_eq!(ends.len(), 12);
                assert!(ends.values().all(|ends| ends.iter().copied().eq([0, 1])));
                (after_strips - before_strips, after_vertices - after_strips)
            });
            assert!(!work.exhausted, "radius={requested}: {work:?}");
            println!(
                "vertex candidates radius={requested} strip_work={strip_work} vertex_work={vertex_work} total_work={}; not closed-body acceptance",
                work.used
            );
        }
    }
}
