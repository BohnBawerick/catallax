//! A plane removes the two neighbor parameters from the intersection solve.
//! Cubic Hermite world/native curves retain SMLib-style composition samples.
//! The fitted curves are numerical; no between-sample bound is claimed.

use catallax_numbers::{Ball, Decider, NumericalCheck, NumericalSample, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, Surface};

use super::super::{BITS, bad, cross, dot, unit, upper_distance};
use super::EndCurves;

struct Node {
    at: Rational,
    values: [[Rational; 3]; 3],
    slopes: [[Rational; 3]; 3],
}

struct Fit<'a> {
    strip: &'a BSplineSurface,
    neighbor: &'a Surface,
    origin: [Ball; 3],
    axes: [[Ball; 3]; 2],
    normal: [Ball; 3],
    tolerance: &'a Rational,
    dec: &'a mut Decider,
    knots: Vec<Rational>,
    controls: [Vec<Point3>; 3],
    samples: Vec<NumericalSample>,
}

impl Fit<'_> {
    fn node(&mut self, at: Rational, mut parameter: Rational) -> Result<Node, Refusal> {
        let domain = self.strip.domain()[0];
        let limit = self.tolerance.mul_pow2(-6);
        for _ in 0..32 {
            self.dec.charge_work(1)?;
            let jet = self
                .strip
                .jet_box([&parameter, &parameter], [&at, &at], BITS, self.dec)?;
            let offset = std::array::from_fn(|i| &jet.value[i] - &self.origin[i]);
            let equation = dot(&self.normal, &offset);
            let du = dot(&self.normal, &jet.first[0]);
            let dv = dot(&self.normal, &jet.first[1]);
            let inverse = du
                .recip()
                .ok_or_else(|| bad("plane end is not transverse to the retained strip"))?;
            if equation.lo().abs().max(equation.hi().abs()) <= limit {
                let slope = -(&dv * &inverse);
                let world_slope: [Ball; 3] =
                    std::array::from_fn(|i| &(&jet.first[0][i] * &slope) + &jet.first[1][i]);
                let uv = self.axes.each_ref().map(|axis| dot(axis, &offset).mid());
                let uv_slope = self
                    .axes
                    .each_ref()
                    .map(|axis| dot(axis, &world_slope).mid());
                let [neighbor_u, neighbor_v] = uv;
                let [neighbor_du, neighbor_dv] = uv_slope;
                return Ok(Node {
                    at: at.clone(),
                    values: [
                        jet.value.map(|value| value.mid()),
                        [parameter, at, Rational::zero()],
                        [neighbor_u, neighbor_v, Rational::zero()],
                    ],
                    slopes: [
                        world_slope.map(|value| value.mid()),
                        [slope.mid(), Rational::one(), Rational::zero()],
                        [neighbor_du, neighbor_dv, Rational::zero()],
                    ],
                });
            }
            let step = equation.mid().to_f64_nearest() / du.mid().to_f64_nearest();
            let proposal = parameter.to_f64_nearest() - step;
            let next = Rational::from_f64_exact(proposal)
                .ok_or_else(|| bad("plane-end scalar solve produced a non-finite proposal"))?
                .max(domain[0].clone())
                .min(domain[1].clone());
            if next == parameter {
                return Err(bad(
                    "plane-end scalar solve stalled in the finite strip domain",
                ));
            }
            parameter = next;
        }
        Err(bad(
            "plane-end scalar solve did not meet its numerical tolerance",
        ))
    }

    fn interval(&mut self, left: &Node, right: &Node, depth: u32) -> Result<(), Refusal> {
        let third = (&right.at - &left.at) / Rational::from_i64(3);
        let controls: [[[Rational; 3]; 4]; 3] = std::array::from_fn(|curve| {
            [
                left.values[curve].clone(),
                std::array::from_fn(|i| {
                    &left.values[curve][i] + &(&third * &left.slopes[curve][i])
                }),
                std::array::from_fn(|i| {
                    &right.values[curve][i] - &(&third * &right.slopes[curve][i])
                }),
                right.values[curve].clone(),
            ]
        });
        let mut samples: [Option<NumericalSample>; 9] = std::array::from_fn(|_| None);
        let mut accepted = true;
        let limit = self.tolerance.mul_pow2(-2);
        let domain = self.strip.domain()[0];
        for index in 0_i64..=8 {
            self.dec.charge_work(1)?;
            let complement = 8 - index;
            let weights = [
                complement.pow(3),
                3 * complement.pow(2) * index,
                3 * complement * index.pow(2),
                index.pow(3),
            ]
            .map(Rational::from_i64);
            let values: [[Rational; 3]; 3] = std::array::from_fn(|curve| {
                std::array::from_fn(|i| {
                    (0..4)
                        .fold(Rational::zero(), |sum, j| {
                            &sum + &(&weights[j] * &controls[curve][j][i])
                        })
                        .mul_pow2(-9)
                })
            });
            if values[1][0] < *domain[0] || values[1][0] > *domain[1] {
                accepted = false;
                break;
            }
            let on_strip = self.strip.evaluate_box(
                [&values[1][0], &values[1][0]],
                [&values[1][1], &values[1][1]],
                BITS,
                self.dec,
            )?;
            let uv = [&values[2][0], &values[2][1]].map(|value| Ball::from_rational(value, BITS));
            let on_neighbor = self.neighbor.evaluate(&uv[0], &uv[1], self.dec)?;
            let error =
                upper_distance(&on_strip, &values[0]).max(upper_distance(&on_neighbor, &values[0]));
            if error > limit {
                accepted = false;
            }
            let share = Rational::from_ratio(index, 8).unwrap();
            samples[index as usize] = Some(NumericalSample {
                parameters: vec![&left.at + &((&right.at - &left.at) * share)],
                distance_upper_mm: error,
            });
        }
        if !accepted {
            if depth == 24 {
                return Err(bad(
                    "plane end cannot meet its sampled world/native tolerance within the subdivision limit",
                ));
            }
            let at = (&left.at + &right.at).mul_pow2(-1);
            let seed = (&left.values[1][0] + &right.values[1][0]).mul_pow2(-1);
            let middle = self.node(at, seed)?;
            self.interval(left, &middle, depth + 1)?;
            return self.interval(&middle, right, depth + 1);
        }
        if self.controls[0].len() >= 3 * 16_384 + 1 {
            return Err(bad(
                "plane end exceeds the retained numerical station limit",
            ));
        }
        let skip = usize::from(!self.controls[0].is_empty());
        for (output, row) in self.controls.iter_mut().zip(controls) {
            output.extend(
                row.into_iter()
                    .skip(skip)
                    .map(|[x, y, z]| Point3::new(x.into(), y.into(), z.into())),
            );
        }
        self.knots.extend(std::iter::repeat_n(
            right.at.clone(),
            if right.at == Rational::one() { 4 } else { 3 },
        ));
        self.samples.extend(samples.into_iter().flatten());
        Ok(())
    }
}

pub(super) fn construct(
    strip: &BSplineSurface,
    neighbor: &Surface,
    seed: Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<EndCurves, Refusal> {
    let Surface::Plane(plane) = neighbor else {
        return Err(bad("scalar plane end requires an actual plane neighbor"));
    };
    let zero = Ball::from_rational(&Rational::zero(), BITS);
    let axes = neighbor.derivatives(&zero, &zero, dec)?;
    let normal = unit(cross(&axes[0], &axes[1]))?;
    let origin = [
        dec.enclose_budgeted_at(&plane.origin.x, BITS)?,
        dec.enclose_budgeted_at(&plane.origin.y, BITS)?,
        dec.enclose_budgeted_at(&plane.origin.z, BITS)?,
    ];
    let mut fit = Fit {
        strip,
        neighbor,
        origin,
        axes,
        normal,
        tolerance,
        dec,
        knots: vec![Rational::zero(); 4],
        controls: std::array::from_fn(|_| Vec::new()),
        samples: Vec::new(),
    };
    let left = fit.node(Rational::zero(), seed)?;
    let right = fit.node(Rational::one(), left.values[1][0].clone())?;
    let tau = [left.values[1][0].clone(), right.values[1][0].clone()];
    fit.interval(&left, &right, 0)?;
    let Fit {
        knots,
        controls,
        samples,
        dec,
        ..
    } = fit;
    let [world, strip, neighbor] = controls;
    let build = |knots: Vec<Rational>, controls: Vec<Point3>, dec: &mut Decider| {
        let count = controls.len();
        BSplineCurve::new(3, knots, controls, vec![Real::one(); count], dec)
    };
    let curve = build(knots.clone(), world, dec)?;
    let strip = build(knots.clone(), strip, dec)?;
    let neighbor = build(knots, neighbor, dec)?;
    Ok(EndCurves {
        curve,
        parameters: [strip, neighbor],
        tau,
        check: NumericalCheck::Samples {
            name: "fillet.vertex.plane-end-world-and-native-curves".into(),
            scope: "scalar plane constraint; nine composition samples per accepted cubic span, not a Hausdorff bound".into(),
            samples,
        },
    })
}
