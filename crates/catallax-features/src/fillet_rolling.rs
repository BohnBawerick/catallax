//! One guide-normal rolling-ball construction over regular support adapters.
//! Surface recipes do not select the algorithm or qualify its numerical result.

use catallax_numbers::{Ball, Decider, NumericalCheck, Rational, Refusal, Sign};
use catallax_shapes::{
    NumericalFilletSpine, NumericalFilletSupport, NumericalSplineFillet, SplineSupport, Surface,
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
                    numerical_embedding::surface_parameters(
                        &surfaces[side],
                        &guide.ends[0],
                        tolerance,
                        dec,
                    )?,
                    numerical_embedding::surface_parameters(
                        &surfaces[side],
                        &guide.ends[1],
                        tolerance,
                        dec,
                    )?,
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
            let fillet = NumericalSplineFillet::fit(spine, dec).map_err(|error| {
                super::at_stage(
                    &format!(
                        "general strip and contact rails for edge {:?}, supports {:?}, name {:?}, guide margin {:.6e}",
                        source.edge, source.faces, model.edge_name(body, source.edge),
                        margin.to_f64_nearest(),
                    ),
                    error,
                )
            })?;
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
    use catallax_shapes::Point3;

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
    fn curved_owner_strip_keeps_its_radius_across_both_guide_extensions() {
        let tolerance = Rational::parse("0.0000001").unwrap();
        let (fillet, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let recipe = catallax_recipe::Recipe::parse(include_str!(
                "../tests/fixtures/owner-curved-pocket.json"
            ))
            .unwrap();
            let source = crate::evaluate(&recipe, &mut dec).unwrap();
            let selection: catallax_recipe::EdgeSelection = serde_json::from_str(
                r#"[
                {"between":[
                    {"step":"blank","role":{"side":{"segment":"right-round"}}},
                    {"step":"pocket-tool","role":{"side":{"segment":"near"}}}
                ]}
            ]"#,
            )
            .unwrap();
            let selected = crate::fillet_build::select(
                &source.model,
                source.result,
                &selection,
                true,
                &[],
                &mut dec,
            )
            .unwrap();
            construct(
                &source.model,
                source.result,
                &selected,
                &Rational::one(),
                &tolerance,
                &mut dec,
            )
            .unwrap()
            .into_iter()
            .flat_map(|chain| chain.strips)
            .find(|strip| strip.source.edge == selected[0])
            .unwrap()
            .fillet
        });
        println!(
            "curved owner source and one regular rolling strip work={}",
            work.used
        );
        assert!(!work.exhausted, "{work:?}");

        let dec = Decider::new(BudgetPolicy {
            max_work: 100_000_000,
            ..Default::default()
        });
        let carrier = fillet.surface().source();
        let domain = fillet.spine().centre_curve().domain();
        let width = domain[1] - domain[0];
        let x = Ball::from_rational(&Rational::from_i64(15), 128);
        let z = Ball::from_rational(&Rational::from_i64(35), 128);
        let centre_y = Rational::parse("-50.335").unwrap();
        let radius = Rational::from_i64(24);
        let check_tolerance = tolerance.mul_pow2(1);
        for i in 0..=17 {
            let share = Rational::from_ratio(i, 17).unwrap();
            let tau = domain[0] + &(&width * &share);
            let centre = fillet
                .spine()
                .centre_curve()
                .evaluate_box([&tau, &tau], 128, &dec)
                .unwrap();
            within(&centre[1], &centre_y, &check_tolerance);
            let dx = &centre[0] - &x;
            let dz = &centre[2] - &z;
            within(
                &(dx.square() + dz.square()).sqrt().unwrap(),
                &radius,
                &check_tolerance,
            );
            for j in 0..=7 {
                let v = Rational::from_ratio(j, 7).unwrap();
                let point = carrier
                    .evaluate_box([&tau, &tau], [&v, &v], 128, &dec)
                    .unwrap();
                let square = (0..3)
                    .fold(Ball::from_rational(&Rational::zero(), 128), |sum, axis| {
                        &sum + &(&point[axis] - &centre[axis]).square()
                    });
                within(&square.sqrt().unwrap(), &Rational::one(), &check_tolerance);
            }
        }
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
                "untrimmed box strips radius={radius} work={}; not a closed-body acceptance",
                work.used
            );
        }
    }

    #[test]
    fn nonplanar_neighbour_cut_retains_both_native_curves() {
        let ((), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let radius = Rational::one();
            let tolerance = Rational::parse("0.0000001").unwrap();
            let mut model = Model::default();
            let body = catallax_topology::block(
                &mut model,
                [20, 30, 40].map(Rational::from_i64),
                &catallax_shapes::Pose::identity(),
            )
            .unwrap();
            let edge = model
                .edges_of_body(body)
                .into_iter()
                .find(|edge| {
                    let edge = model.edge(*edge);
                    [edge.start, edge.end].into_iter().all(|vertex| {
                        let point = model.geom().position(model.vertex(vertex.unwrap()).point);
                        point.x.as_rational().is_some_and(|q| q.is_zero())
                            && point.y.as_rational().is_some_and(|q| q.is_zero())
                    })
                })
                .unwrap();
            let ends = [model.edge(edge).start, model.edge(edge).end].map(|vertex| {
                model
                    .geom()
                    .position(model.vertex(vertex.unwrap()).point)
                    .z
                    .as_rational()
                    .unwrap()
            });
            let mut chains =
                construct(&model, body, &[edge], &radius, &tolerance, &mut dec).unwrap();
            let strip = chains.remove(0).strips.remove(0);

            // Quadratic tensor patch: x,y in [-1,10], z=10+x²/100+y²/50.
            // The independent oracle below uses that equation, not the patch.
            let xy = [
                Rational::from_i64(-1),
                Rational::from_ratio(9, 2).unwrap(),
                Rational::from_i64(10),
            ];
            let squared = [1, -10, 100].map(Rational::from_i64);
            let a = Rational::from_ratio(1, 100).unwrap();
            let b = Rational::from_ratio(1, 50).unwrap();
            let mut controls = Vec::with_capacity(9);
            for i in 0..3 {
                for j in 0..3 {
                    let height = Rational::from_i64(10) + &a * &squared[i] + &b * &squared[j];
                    controls.push(Point3::new(
                        xy[i].clone().into(),
                        xy[j].clone().into(),
                        height.into(),
                    ));
                }
            }
            let surface = catallax_shapes::BSplineSurface::new(
                [2, 2],
                std::array::from_fn(|_| [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec()),
                [3, 3],
                controls,
                vec![catallax_numbers::Real::one(); 9],
                &mut dec,
            )
            .unwrap();
            let [atlas] = <[_; 1]>::try_from(
                catallax_shapes::RationalAtlas::partition(surface, [false, false], &mut dec)
                    .unwrap(),
            )
            .unwrap();
            let guess = (&Rational::from_i64(10) - &ends[0]) / (&ends[1] - &ends[0]);
            let foot = strip
                .fillet
                .rail(0)
                .unwrap()
                .evaluate_box([&guess, &guess], 128, &dec)
                .unwrap();
            let eleven = Rational::from_i64(11);
            let seed = [
                guess,
                &(&foot[0].mid() + &Rational::one()) / &eleven,
                &(&foot[1].mid() + &Rational::one()) / &eleven,
            ];
            let cut = strip
                .fillet
                .surface_cut(
                    Surface::BSpline(SplineSupport::RationalAtlas(atlas)),
                    seed,
                    &mut dec,
                )
                .unwrap();
            let one = Ball::from_rational(&Rational::one(), 128);
            let ten = Ball::from_rational(&Rational::from_i64(10), 128);
            let eleven = Ball::from_rational(&eleven, 128);
            let a = Ball::from_rational(&a, 128);
            let b = Ball::from_rational(&b, 128);
            let start = Ball::from_rational(&ends[0], 128);
            let length = Ball::from_rational(&(&ends[1] - &ends[0]), 128);
            for j in 0..=32 {
                let t = Rational::from_ratio(j, 32).unwrap();
                let point = cut.curve().evaluate_box([&t, &t], 128, &dec).unwrap();
                let uv = cut
                    .parameters()
                    .each_ref()
                    .map(|curve| curve.evaluate_box([&t, &t], 128, &dec).unwrap());
                let dx = &point[0] - &one;
                let dy = &point[1] - &one;
                within(
                    &(dx.square() + dy.square()).sqrt().unwrap(),
                    &radius,
                    &tolerance,
                );
                let height = &(&ten + &(&a * &point[0].square())) + &(&b * &point[1].square());
                within(&(&point[2] - &height), &Rational::zero(), &tolerance);
                let z = &start + &(&length * &uv[0][0]);
                within(&(&point[2] - &z), &Rational::zero(), &tolerance);
                for axis in 0..2 {
                    let coordinate = &(&eleven * &uv[1][axis]) - &one;
                    within(&(&point[axis] - &coordinate), &Rational::zero(), &tolerance);
                }
            }
        });
        println!("general nonplanar neighbour cut work={}", work.used);
    }
}
