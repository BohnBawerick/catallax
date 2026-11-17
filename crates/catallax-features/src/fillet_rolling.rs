//! One guide-normal rolling-ball construction over regular support adapters.
//! Surface recipes do not select the algorithm or qualify its numerical result.

use catallax_numbers::{Ball, Decider, NumericalCheck, Rational, Refusal, Sign};
use catallax_shapes::{
    NumericalFilletSpine, NumericalFilletSupport, NumericalSplineFillet, Point3, SplineSupport,
    Surface,
};
use catallax_topology::{BodyId, EdgeId, Model, Sense};

use crate::{
    fillet_build::chains::{self, ChainEdge},
    numerical_embedding,
};

pub(crate) struct RolledStrip {
    pub source: ChainEdge,
    pub fillet: NumericalSplineFillet,
    pub convex: bool,
}

pub(crate) struct RolledChain {
    pub strips: Vec<RolledStrip>,
    pub closed: bool,
    pub checks: Vec<NumericalCheck>,
}

struct SideSeed {
    endpoints: [[Rational; 2]; 2],
    outward: [Ball; 3],
    alignment: Sign,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-contact", detail)
}

fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn separated_sign(value: &Ball, detail: &str) -> Result<Sign, Refusal> {
    if value.lo().sign() == Sign::Positive {
        Ok(Sign::Positive)
    } else if value.hi().sign() == Sign::Negative {
        Ok(Sign::Negative)
    } else {
        Err(bad(detail))
    }
}

// These are initial guesses, not projected-point or incidence certificates.
// A chart seam may have no unique inverse. A small Cartesian displacement can
// choose a neighbouring chart; the contact walk then checks the actual support.
fn parameters(
    surface: &Surface,
    point: &Point3,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let mut p = [
        dec.enclose_budgeted_at(&point.x, 128)?,
        dec.enclose_budgeted_at(&point.y, 128)?,
        dec.enclose_budgeted_at(&point.z, 128)?,
    ];
    if let Some(uv) = surface.parameters_of(&p, dec)? {
        return Ok(uv.each_ref().map(Ball::mid));
    }
    let displacement = Ball::from_rational(&tolerance.mul_pow2(-6), 128);
    for axis in 0..3 {
        let shifted = &p[axis] + &displacement;
        let original = std::mem::replace(&mut p[axis], shifted);
        let proposal = surface.parameters_of(&p, dec)?;
        p[axis] = original;
        if let Some(uv) = proposal {
            return Ok(uv.each_ref().map(Ball::mid));
        }
    }
    Err(bad(
        "the selected edge endpoint has no regular support-chart seed",
    ))
}

fn support(
    surface: Surface,
    distance: Rational,
    endpoints: &[[Rational; 2]; 2],
    margin: &Rational,
    dec: &mut Decider,
) -> Result<NumericalFilletSupport, Refusal> {
    let surface = match surface {
        Surface::BSpline(SplineSupport::RationalAtlas(mut atlas)) if !margin.is_zero() => {
            let domain = super::domain(&atlas);
            for axis in 0..2 {
                let delta = (&endpoints[1][axis] - &endpoints[0][axis]).abs();
                if delta.is_zero() {
                    continue;
                }
                let width = &domain[axis][1] - &domain[axis][0];
                let fraction = &(margin * &delta) / &width;
                atlas = NumericalFilletSpine::extended_support(&atlas, axis, &fraction, dec)?;
            }
            Surface::BSpline(SplineSupport::RationalAtlas(atlas))
        }
        other => other,
    };
    NumericalFilletSupport::new(surface, distance, dec)
}

/// Builds candidate strips and their two contact rails. Real face-domain trims
/// and junction blends must be constructed before this becomes an admitted body.
/// No exact carrier-intersection certificate gates this construction.
pub(crate) fn construct(
    model: &Model,
    body: BodyId,
    selected: &[EdgeId],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Vec<RolledChain>, Refusal> {
    if radius.sign() != Sign::Positive || tolerance.sign() != Sign::Positive {
        return Err(bad(
            "rolling-ball radius and numerical tolerance must be positive",
        ));
    }
    let mut output = Vec::new();
    for chain in chains::group(model, body, selected, radius, tolerance, dec)? {
        let mut result = RolledChain {
            strips: Vec::with_capacity(chain.edges.len()),
            closed: chain.closed,
            checks: chain.checks,
        };
        for source in chain.edges {
            dec.charge_work(1)?;
            let guide = numerical_embedding::fillet_guide(model, source.edge, tolerance, dec)?;
            let surfaces = source
                .faces
                .map(|face| model.geom().surface(model.face(face).surface).clone());
            let tangent = guide
                .curve
                .jet_box([&Rational::zero(), &Rational::zero()], 128, dec)?
                .first;
            let end_tangent = guide
                .curve
                .jet_box([&Rational::one(), &Rational::one()], 128, dec)?
                .first;
            let speed = dot(&tangent, &tangent)
                .sqrt()
                .ok_or_else(|| bad("guide start speed cannot be enclosed"))?
                .lo()
                .min(
                    dot(&end_tangent, &end_tangent)
                        .sqrt()
                        .ok_or_else(|| bad("guide end speed cannot be enclosed"))?
                        .lo(),
                );
            if speed.sign() != Sign::Positive {
                return Err(bad("the selected edge has a singular guide endpoint"));
            }
            let is_ring = model.edge(source.edge).start == model.edge(source.edge).end;
            let margin = if is_ring {
                Rational::zero()
            } else {
                radius.mul_pow2(1) / speed
            };
            let seed = |side: usize, dec: &mut Decider| -> Result<SideSeed, Refusal> {
                let endpoints = [
                    parameters(&surfaces[side], &guide.ends[0], tolerance, dec)?,
                    parameters(&surfaces[side], &guide.ends[1], tolerance, dec)?,
                ];
                let normal = surfaces[side].normal_at(&guide.ends[0], dec)?;
                let values = [
                    dec.enclose_budgeted_at(&normal.x, 128)?,
                    dec.enclose_budgeted_at(&normal.y, 128)?,
                    dec.enclose_budgeted_at(&normal.z, 128)?,
                ];
                let outward = if model.face(source.faces[side]).sense == Sense::Same {
                    values
                } else {
                    values.map(|v| -v)
                };
                let uv = endpoints[0].each_ref().map(|q| Ball::from_rational(q, 128));
                let [du, dv] = surfaces[side].derivatives(&uv[0], &uv[1], dec)?;
                let alignment = separated_sign(
                    &dot(&cross(&du, &dv), &outward),
                    "support chart orientation cannot be separated from the source face normal",
                )?;
                Ok(SideSeed {
                    endpoints,
                    outward,
                    alignment,
                })
            };
            let seeds = [seed(0, dec)?, seed(1, dec)?];
            let mut use_sense = None;
            for boundary in &model.face(source.faces[0]).loops {
                if let catallax_topology::Loop::Coedges(uses) = model.loop_(*boundary) {
                    for coedge in uses {
                        let coedge = model.coedge(*coedge);
                        if coedge.edge == source.edge {
                            use_sense = Some(coedge.sense);
                        }
                    }
                }
            }
            let use_sense = use_sense
                .ok_or_else(|| bad("the selected edge has no first-support occurrence"))?;
            let orientation = dot(&cross(&seeds[0].outward, &seeds[1].outward), &tangent);
            let sign = separated_sign(
                &orientation,
                "the selected edge has a smooth, singular or unresolved dihedral",
            )?;
            let convex = (sign == Sign::Positive) == (use_sense == Sense::Same);
            let offsets = seeds.each_ref().map(|seed| {
                if convex == (seed.alignment == Sign::Positive) {
                    -radius
                } else {
                    radius.clone()
                }
            });
            let [first, second] = surfaces;
            let [first_offset, second_offset] = offsets;
            let supports = [
                support(first, first_offset, &seeds[0].endpoints, &margin, dec)?,
                support(second, second_offset, &seeds[1].endpoints, &margin, dec)?,
            ];
            let [[first_seed, _], [second_seed, _]] = seeds.map(|seed| seed.endpoints);
            let spine = NumericalFilletSpine::walk(
                supports,
                guide.curve,
                [-&margin, &Rational::one() + &margin],
                [first_seed, second_seed],
                tolerance.clone(),
                dec,
            )
            .map_err(|error| super::at_stage("general rolling-ball contact walk", error))?;
            let fillet = NumericalSplineFillet::fit(spine, dec)
                .map_err(|error| super::at_stage("general strip and contact rails", error))?;
            result.checks.push(guide.check);
            result.strips.push(RolledStrip {
                source,
                fillet,
                convex,
            });
        }
        output.push(result);
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};

    fn signed_plane_distance(
        model: &Model,
        face: catallax_topology::FaceId,
        point: &[Ball; 3],
        dec: &mut Decider,
    ) -> Ball {
        let definition = model.face(face);
        let Surface::Plane(plane) = model.geom().surface(definition.surface) else {
            panic!("box support is a plane");
        };
        let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z]
            .map(|value| dec.enclose_budgeted_at(value, 128).unwrap());
        let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z]
            .map(|value| dec.enclose_budgeted_at(value, 128).unwrap());
        let delta = std::array::from_fn(|i| &point[i] - &origin[i]);
        let value = dot(&delta, &normal);
        if definition.sense == Sense::Same {
            value
        } else {
            -value
        }
    }

    fn within(value: &Ball, expected: &Rational, tolerance: &Rational) {
        assert!((&value.lo() - expected).abs() <= *tolerance, "{value:?}");
        assert!((&value.hi() - expected).abs() <= *tolerance, "{value:?}");
    }

    #[test]
    fn rolling_box_contacts_keep_the_material_side_for_every_edge_orientation() {
        for radius in ["0.5", "1", "2", "5"] {
            let ((), work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 100_000_000,
                    ..Default::default()
                });
                let radius = Rational::parse(radius).unwrap();
                let tolerance = Rational::parse("0.0000001").unwrap();
                let mut model = Model::default();
                let body = catallax_topology::block(
                    &mut model,
                    [20, 30, 40].map(Rational::from_i64),
                    &catallax_shapes::Pose::identity(),
                )
                .unwrap();
                let selected = model.edges_of_body(body);
                let chains =
                    construct(&model, body, &selected, &radius, &tolerance, &mut dec).unwrap();
                assert_eq!(
                    chains.iter().map(|chain| chain.strips.len()).sum::<usize>(),
                    12
                );
                let half = Rational::from_ratio(1, 2).unwrap();
                for chain in chains {
                    for strip in chain.strips {
                        let centre = strip
                            .fillet
                            .spine()
                            .centre_curve()
                            .evaluate_box([&half, &half], 128, &dec)
                            .unwrap();
                        for side in 0..2 {
                            within(
                                &signed_plane_distance(
                                    &model,
                                    strip.source.faces[side],
                                    &centre,
                                    &mut dec,
                                ),
                                &(-&radius),
                                &tolerance,
                            );
                            let foot = strip
                                .fillet
                                .rail(side)
                                .unwrap()
                                .evaluate_box([&half, &half], 128, &dec)
                                .unwrap();
                            within(
                                &signed_plane_distance(
                                    &model,
                                    strip.source.faces[side],
                                    &foot,
                                    &mut dec,
                                ),
                                &Rational::zero(),
                                &tolerance,
                            );
                        }
                        for across in [Rational::zero(), half.clone(), Rational::one()] {
                            let point = strip
                                .fillet
                                .surface()
                                .source()
                                .evaluate_box([&half, &half], [&across, &across], 128, &dec)
                                .unwrap();
                            let delta: [Ball; 3] = std::array::from_fn(|i| &point[i] - &centre[i]);
                            let distance = dot(&delta, &delta).sqrt().unwrap();
                            within(&distance, &radius, &tolerance);
                        }
                    }
                }
            });
            println!(
                "untrimmed box strips radius={radius} work={work}; not a closed-body acceptance"
            );
        }
    }
}
