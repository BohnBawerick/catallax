//! One fixed-strip section and its paired native curves. Sphere agreement is
//! sampled numerical evidence; the strip isocurve itself is an exact extraction.

use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalSample, Rational, Real, Refusal, Sign,
};
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, Surface, Vector3};

use super::super::rolling::RolledChain;
use super::{BITS, VertexBlend, VertexSection, dot, unit};

#[derive(Clone, Debug)]
pub(crate) struct BoundaryContact {
    /// Index into the common ball's three contacts, not a topology ID.
    pub contact: usize,
    pub xyz: [Rational; 3],
    pub sphere_uv: [Rational; 2],
}

#[derive(Clone, Debug)]
pub(crate) struct VertexBoundary {
    /// Shared by the strip and cap, with parameter 0..1 running strip V min..max.
    pub curve: BSplineCurve,
    /// `[fillet UV, sphere UV]`, both at the shared world-curve parameter.
    pub parameters: [BSplineCurve; 2],
    /// Canonical common-ball contacts in the world curve's forward order.
    pub contacts: [BoundaryContact; 2],
    /// The cap loop walks all three curves backwards when this is true.
    pub reversed: bool,
    pub sampled_error_mm: Rational,
    pub check: NumericalCheck,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-vertex-boundary", detail)
}

fn point(uv: &[Rational; 2]) -> Point3 {
    Point3::new(
        Real::from_rational(uv[0].clone()),
        Real::from_rational(uv[1].clone()),
        Real::zero(),
    )
}

fn gap(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        let delta = &a[i] - &b[i];
        sum + delta.lo().abs().max(delta.hi().abs())
    })
}

fn enclosed(point: &[Rational; 3]) -> [Ball; 3] {
    point.each_ref().map(|x| Ball::from_rational(x, BITS))
}

fn vector(v: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&v.x, BITS)?,
        dec.enclose_budgeted_at(&v.y, BITS)?,
        dec.enclose_budgeted_at(&v.z, BITS)?,
    ])
}

/// Contract the homogeneous U control rows with exact de Boor coefficients.
pub(crate) fn isocurve(
    surface: &BSplineSurface,
    tau: &Rational,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let [nu, nv] = surface.control_count();
    let [degree, across] = surface.degrees();
    let [knots, vknots] = surface.knots();
    let domain = surface.domain();
    if tau < domain[0][0] || tau > domain[0][1] {
        return Err(bad("fixed section leaves the fitted strip domain"));
    }
    let span = knots
        .partition_point(|knot| knot <= tau)
        .saturating_sub(1)
        .min(nu - 1);
    let mut row = Vec::with_capacity(degree + 1);
    let mut controls = Vec::with_capacity(nv);
    let mut weights = Vec::with_capacity(nv);
    for v in 0..nv {
        row.clear();
        for u in span - degree..=span {
            let i = u * nv + v;
            let p = &surface.controls()[i];
            let w = &surface.weights()[i];
            row.push([&p.x * w, &p.y * w, &p.z * w, w.clone()]);
        }
        for r in 1..=degree {
            for j in (r..=degree).rev() {
                dec.charge_work(1)?;
                let i = span - degree + j;
                let denominator = &knots[i + degree + 1 - r] - &knots[i];
                if denominator.sign() != Sign::Positive {
                    return Err(bad("fixed section has a zero-support knot interval"));
                }
                let alpha = Real::from_rational((tau - &knots[i]) / denominator);
                let complement = &Real::one() - &alpha;
                let (before, current) = row.split_at_mut(j);
                for k in 0..4 {
                    current[0][k] = &(&complement * &before[j - 1][k]) + &(&alpha * &current[0][k]);
                }
            }
        }
        let [x, y, z, w] = row
            .pop()
            .ok_or_else(|| bad("fixed section has no control row"))?;
        controls.push(Point3::new(&x / &w, &y / &w, &z / &w));
        weights.push(w);
    }
    let width = domain[1][1] - domain[1][0];
    let knots = vknots
        .iter()
        .map(|knot| &(knot - domain[1][0]) / &width)
        .collect();
    BSplineCurve::new(across, knots, controls, weights, dec)
}

struct Chart {
    origin: [Ball; 3],
    axes: [[Ball; 3]; 3],
    angle: Rational,
    cosine: Ball,
    sine: Ball,
}

impl Chart {
    fn new(vertex: &VertexBlend, dec: &mut Decider) -> Result<Self, Refusal> {
        let Surface::Sphere(sphere) = &vertex.surface else {
            return Err(bad("common-ball boundary needs its spherical carrier"));
        };
        let origin = [
            dec.enclose_budgeted_at(&sphere.center.x, BITS)?,
            dec.enclose_budgeted_at(&sphere.center.y, BITS)?,
            dec.enclose_budgeted_at(&sphere.center.z, BITS)?,
        ];
        let axes = [
            unit(vector(&sphere.ref_dir, dec)?)?,
            unit(vector(&sphere.axis.cross(&sphere.ref_dir), dec)?)?,
            unit(vector(&sphere.axis, dec)?)?,
        ];
        let mut angles = [0.0; 3];
        for (i, contact) in vertex.ball.contacts().iter().enumerate() {
            let delta =
                std::array::from_fn(|j| &Ball::from_rational(&contact[j], BITS) - &origin[j]);
            angles[i] = dot(&delta, &axes[1])
                .mid()
                .to_f64_nearest()
                .atan2(dot(&delta, &axes[0]).mid().to_f64_nearest());
        }
        angles.sort_by(f64::total_cmp);
        let mut largest = 0;
        let mut width = f64::NEG_INFINITY;
        for i in 0..3 {
            let next = angles[(i + 1) % 3] + if i == 2 { std::f64::consts::TAU } else { 0.0 };
            if next - angles[i] > width {
                width = next - angles[i];
                largest = i;
            }
        }
        // The largest contact gap proposes a seam outside this corner. Every
        // inverse and subsequent world/native composition is still checked.
        let start = angles[(largest + 1) % 3];
        let end = angles[largest]
            + if largest == 2 {
                0.0
            } else {
                std::f64::consts::TAU
            };
        let angle = Rational::from_f64_exact((start + end) * 0.5)
            .ok_or_else(|| bad("corner chart proposal is not finite"))?;
        let rotation = Ball::from_rational(&angle, BITS);
        Ok(Self {
            origin,
            axes,
            angle,
            cosine: rotation.cos(),
            sine: rotation.sin(),
        })
    }

    fn coordinates(&self, value: &[Ball; 3]) -> [Ball; 3] {
        let x = dot(value, &self.axes[0]);
        let y = dot(value, &self.axes[1]);
        [
            &(&x * &self.cosine) + &(&y * &self.sine),
            &(&y * &self.cosine) - &(&x * &self.sine),
            dot(value, &self.axes[2]),
        ]
    }

    fn uv(&self, position: &[Ball; 3]) -> Result<[Ball; 2], Refusal> {
        let delta = std::array::from_fn(|i| &position[i] - &self.origin[i]);
        let [x, y, z] = self.coordinates(&delta);
        let rho = (&x.square() + &y.square())
            .sqrt()
            .ok_or_else(|| bad("corner chart radial norm is unresolved"))?;
        if rho.lo().sign() != Sign::Positive {
            return Err(bad("corner boundary reaches a spherical chart pole"));
        }
        let u = Ball::atan2(&y, &x)
            .ok_or_else(|| bad("corner boundary crosses the proposed chart seam"))?;
        let v = Ball::atan2(&z, &rho).ok_or_else(|| bad("corner latitude is unresolved"))?;
        Ok([&u + &Ball::from_rational(&self.angle, BITS), v])
    }

    fn node(&self, curve: &BSplineCurve, t: Rational, dec: &Decider) -> Result<Node, Refusal> {
        let jet = curve.jet_box([&t, &t], BITS, dec)?;
        let uv = self.uv(&jet.value)?.map(|value| value.mid());
        let delta = std::array::from_fn(|i| &jet.value[i] - &self.origin[i]);
        let [x, y, z] = self.coordinates(&delta);
        let [dx, dy, dz] = self.coordinates(&jet.first);
        let rho_sq = &x.square() + &y.square();
        let rho = rho_sq
            .sqrt()
            .ok_or_else(|| bad("corner chart derivative has no radial norm"))?;
        let du = (&(&x * &dy) - &(&y * &dx))
            .checked_div(&rho_sq)
            .ok_or_else(|| bad("corner longitude derivative reaches a pole"))?;
        let drho = (&(&x * &dx) + &(&y * &dy))
            .checked_div(&rho)
            .ok_or_else(|| bad("corner latitude derivative reaches a pole"))?;
        let dv = (&(&rho * &dz) - &(&z * &drho))
            .checked_div(&(&rho_sq + &z.square()))
            .ok_or_else(|| bad("corner latitude derivative is singular"))?;
        Ok(Node {
            t,
            uv,
            derivative: [du.mid(), dv.mid()],
        })
    }
}

struct Node {
    t: Rational,
    uv: [Rational; 2],
    derivative: [Rational; 2],
}

struct Fit<'a> {
    curve: &'a BSplineCurve,
    fillet: &'a BSplineSurface,
    sphere: &'a Surface,
    chart: &'a Chart,
    tau: &'a Rational,
    tolerance: &'a Rational,
    dec: &'a mut Decider,
    knots: Vec<Rational>,
    controls: Vec<Point3>,
    samples: Vec<NumericalSample>,
    upper: Rational,
}

impl Fit<'_> {
    fn interval(&mut self, left: &Node, right: &Node, depth: usize) -> Result<(), Refusal> {
        let width = &right.t - &left.t;
        let third = &width / &Rational::from_i64(3);
        let controls: [[Rational; 2]; 4] = [
            left.uv.clone(),
            std::array::from_fn(|i| &left.uv[i] + &(&third * &left.derivative[i])),
            std::array::from_fn(|i| &right.uv[i] - &(&third * &right.derivative[i])),
            right.uv.clone(),
        ];
        let mut errors: [Rational; 9] = std::array::from_fn(|_| Rational::zero());
        let mut accepted = true;
        for (i, error) in errors.iter_mut().enumerate() {
            self.dec.charge_work(1)?;
            let share = Rational::from_ratio(i as i64, 8).expect("constant denominator");
            let t = &left.t + &(&width * &share);
            let mut values = controls.clone();
            for level in (1..4).rev() {
                for j in 0..level {
                    for k in 0..2 {
                        values[j][k] = &(&values[j][k] * &(&Rational::one() - &share))
                            + &(&values[j + 1][k] * &share);
                    }
                }
            }
            let uv = values[0].each_ref().map(|x| Ball::from_rational(x, BITS));
            let world = self.curve.evaluate_box([&t, &t], BITS, self.dec)?;
            let on_sphere = self.sphere.evaluate(&uv[0], &uv[1], self.dec)?;
            let domain = self.fillet.domain()[1];
            let v = domain[0] + &(&(domain[1] - domain[0]) * &t);
            let on_strip =
                self.fillet
                    .evaluate_box([self.tau, self.tau], [&v, &v], BITS, self.dec)?;
            *error = gap(&world, &on_sphere).max(gap(&world, &on_strip));
            if *error > *self.tolerance {
                accepted = false;
                break;
            }
        }
        if accepted {
            self.controls.extend(controls[1..].iter().map(point));
            if right.t < Rational::one() {
                self.knots.extend(std::iter::repeat_n(right.t.clone(), 3));
            }
            for (i, error) in errors.into_iter().enumerate() {
                if i == 0 && !self.samples.is_empty() {
                    continue;
                }
                self.upper = self.upper.clone().max(error.clone());
                self.samples.push(NumericalSample {
                    parameters: vec![
                        &left.t
                            + &(&width
                                * &Rational::from_ratio(i as i64, 8)
                                    .expect("constant denominator")),
                    ],
                    distance_upper_mm: error,
                });
            }
            return Ok(());
        }
        if depth == 12 {
            return Err(bad(
                "sphere native-curve fit does not meet its sampled composition tolerance",
            ));
        }
        let middle = self
            .chart
            .node(self.curve, (&left.t + &right.t).mul_pow2(-1), self.dec)?;
        self.interval(left, &middle, depth + 1)?;
        self.interval(&middle, right, depth + 1)
    }
}

/// Build one section from its containing chain. All sections of the same ball
/// choose the same chart and canonical contacts, independent of call order.
pub(crate) fn construct(
    chain: &RolledChain,
    vertex: &VertexBlend,
    section: &VertexSection,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexBoundary, Refusal> {
    if tolerance.sign() != Sign::Positive {
        return Err(bad("boundary tolerance must be positive"));
    }
    let strip = chain
        .strips
        .get(section.strip.strip)
        .ok_or_else(|| bad("fixed section does not identify a strip in this chain"))?;
    let fillet = strip.fillet.surface().source();
    let curve = isocurve(fillet, &section.tau, dec)?;
    let chart = Chart::new(vertex, dec)?;
    let contact = |index: usize| -> Result<BoundaryContact, Refusal> {
        let xyz = vertex
            .ball
            .contacts()
            .get(index)
            .ok_or_else(|| bad("fixed section has no common-ball contact"))?
            .clone();
        let uv = chart.uv(&enclosed(&xyz))?;
        let on_sphere = vertex.surface.evaluate(&uv[0], &uv[1], dec)?;
        if gap(&enclosed(&xyz), &on_sphere) > *tolerance {
            return Err(bad(
                "canonical contact does not meet its sphere composition tolerance",
            ));
        }
        Ok(BoundaryContact {
            contact: index,
            xyz,
            sphere_uv: uv.map(|x| x.mid()),
        })
    };
    let contacts = [contact(section.contacts[0])?, contact(section.contacts[1])?];
    let mut left = chart.node(&curve, Rational::zero(), dec)?;
    let mut right = chart.node(&curve, Rational::one(), dec)?;
    // Canonical endpoints keep the three cap UV boundaries joined. Their
    // agreement with the unmodified strip endpoints remains numerical.
    left.uv = contacts[0].sphere_uv.clone();
    right.uv = contacts[1].sphere_uv.clone();
    let mut endpoint_errors = [Rational::zero(), Rational::zero()];
    for (i, (t, contact)) in [(&left.t, &contacts[0]), (&right.t, &contacts[1])]
        .into_iter()
        .enumerate()
    {
        let world = curve.evaluate_box([t, t], BITS, dec)?;
        endpoint_errors[i] = gap(&world, &enclosed(&contact.xyz));
        if endpoint_errors[i] > *tolerance {
            return Err(bad(
                "fixed section endpoint misses its canonical common-ball contact",
            ));
        }
    }
    let mut fit = Fit {
        curve: &curve,
        fillet,
        sphere: &vertex.surface,
        chart: &chart,
        tau: &section.tau,
        tolerance,
        dec,
        knots: vec![Rational::zero(); 4],
        controls: vec![point(&left.uv)],
        samples: Vec::new(),
        upper: Rational::zero(),
    };
    fit.interval(&left, &right, 0)?;
    for (i, error) in endpoint_errors.into_iter().enumerate() {
        fit.upper = fit.upper.max(error.clone());
        fit.samples.push(NumericalSample {
            parameters: vec![Rational::from_i64(i as i64)],
            distance_upper_mm: error,
        });
    }
    let Fit {
        mut knots,
        controls,
        samples,
        upper,
        dec,
        ..
    } = fit;
    knots.extend(std::iter::repeat_n(Rational::one(), 4));
    let count = controls.len();
    let sphere_uv = BSplineCurve::new(3, knots, controls, vec![Real::one(); count], dec)?;
    let domain = fillet.domain()[1];
    let strip_uv = BSplineCurve::new(
        1,
        [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        domain
            .map(|v| point(&[section.tau.clone(), v.clone()]))
            .to_vec(),
        vec![Real::one(); 2],
        dec,
    )?;
    Ok(VertexBoundary {
        curve,
        parameters: [strip_uv, sphere_uv],
        contacts,
        reversed: section.reversed,
        sampled_error_mm: upper,
        check: NumericalCheck::Samples {
            name: "fillet.vertex.world-and-paired-native-curves".into(),
            scope: "fixed-tau section; finite sample locations only".into(),
            samples,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, current_operation_work, with_operation_budget};
    use catallax_shapes::{Pose, Sphere};
    use catallax_topology::{Model, Sense};
    use std::collections::BTreeMap;

    #[test]
    fn box_vertex_boundaries_share_parameters_contacts_and_seam_safe_orientation() {
        for requested in ["0.5", "1", "2", "5"] {
            let ((strip_work, corner_work, boundary_work), work) =
                with_operation_budget(100_000_000, || {
                    let mut dec = Decider::new(BudgetPolicy {
                        max_work: 100_000_000,
                        ..Default::default()
                    });
                    let radius = Rational::parse(requested).unwrap();
                    let tolerance = Rational::parse("0.0000001").unwrap();
                    let mut model = Model::default();
                    let body = catallax_topology::block(
                        &mut model,
                        [20, 30, 40].map(Rational::from_i64),
                        &Pose::identity(),
                    )
                    .unwrap();
                    let before = current_operation_work().used;
                    let chains = crate::fillet::numerical::rolling::construct(
                        &model,
                        body,
                        &model.edges_of_body(body),
                        &radius,
                        &tolerance,
                        &mut dec,
                    )
                    .unwrap();
                    let after_strips = current_operation_work().used;
                    let mut vertices =
                        super::super::construct(&model, &chains, &radius, &tolerance, &mut dec)
                            .unwrap();
                    let after_corners = current_operation_work().used;
                    let mut boundary_work = 0;
                    let mut crossed_native_seam = false;
                    for (i, vertex) in vertices.iter_mut().enumerate() {
                        if i % 2 == 1 {
                            let Surface::Sphere(sphere) = &vertex.surface else {
                                unreachable!()
                            };
                            vertex.surface = Surface::Sphere(
                                Sphere::new(
                                    sphere.center.clone(),
                                    sphere.radius.clone(),
                                    sphere.axis.clone(),
                                    Vector3::new(
                                        -&sphere.ref_dir.x,
                                        -&sphere.ref_dir.y,
                                        -&sphere.ref_dir.z,
                                    ),
                                    &mut dec,
                                )
                                .unwrap(),
                            );
                        }
                        let mut contacts = BTreeMap::new();
                        let mut loop_ends = Vec::new();
                        let mut area = Rational::zero();
                        for section in &vertex.sections {
                            let chain = &chains[section.strip.chain];
                            let start = current_operation_work().used;
                            let boundary =
                                construct(chain, vertex, section, &tolerance, &mut dec).unwrap();
                            boundary_work += current_operation_work().used - start;
                            assert_eq!(boundary.reversed, section.reversed);
                            assert!(boundary.sampled_error_mm <= tolerance);
                            let NumericalCheck::Samples { samples, .. } = &boundary.check else {
                                panic!("boundary composition needs numerical evidence")
                            };
                            assert!(
                                samples
                                    .iter()
                                    .all(|sample| sample.distance_upper_mm <= tolerance)
                            );
                            for contact in &boundary.contacts {
                                let definition = (contact.xyz.clone(), contact.sphere_uv.clone());
                                if let Some(previous) =
                                    contacts.insert(contact.contact, definition.clone())
                                {
                                    assert_eq!(
                                        previous, definition,
                                        "all sections must use one canonical contact and chart"
                                    );
                                }
                            }
                            let [a, b] = if boundary.reversed {
                                [&boundary.contacts[1], &boundary.contacts[0]]
                            } else {
                                [&boundary.contacts[0], &boundary.contacts[1]]
                            };
                            loop_ends.push((a.contact, b.contact));
                            area = area + &a.sphere_uv[0] * &b.sphere_uv[1]
                                - &b.sphere_uv[0] * &a.sphere_uv[1];
                            for j in 0..=17 {
                                let t = Rational::from_ratio(j, 17).unwrap();
                                let world =
                                    boundary.curve.evaluate_box([&t, &t], BITS, &dec).unwrap();
                                let uv = boundary.parameters[0]
                                    .evaluate_box([&t, &t], BITS, &dec)
                                    .unwrap();
                                let native = chain.strips[section.strip.strip]
                                    .fillet
                                    .surface()
                                    .source()
                                    .evaluate(&uv[0], &uv[1], &dec)
                                    .unwrap();
                                assert!(gap(&world, &native) <= tolerance);
                                assert!(uv[0].lo() <= section.tau && uv[0].hi() >= section.tau);
                                assert!(uv[1].lo() <= t && uv[1].hi() >= t);
                                let uv = boundary.parameters[1]
                                    .evaluate_box([&t, &t], BITS, &dec)
                                    .unwrap();
                                let native = vertex.surface.evaluate(&uv[0], &uv[1], &dec).unwrap();
                                assert!(
                                    gap(&world, &native) <= tolerance,
                                    "radius={requested}, t={t}"
                                );
                                crossed_native_seam |= uv[0].lo() > Ball::pi(BITS).hi();
                            }
                            for (i, contact) in boundary.contacts.iter().enumerate() {
                                let t = Rational::from_i64(i as i64);
                                let endpoint =
                                    boundary.curve.evaluate_box([&t, &t], BITS, &dec).unwrap();
                                assert!(gap(&endpoint, &enclosed(&contact.xyz)) <= tolerance);
                                let uv = boundary.parameters[1]
                                    .evaluate_box([&t, &t], BITS, &dec)
                                    .unwrap();
                                for (value, expected) in uv[..2].iter().zip(&contact.sphere_uv) {
                                    assert!(value.lo() <= *expected && value.hi() >= *expected);
                                }
                            }
                        }
                        for j in 0..loop_ends.len() {
                            assert_eq!(loop_ends[j].1, loop_ends[(j + 1) % loop_ends.len()].0);
                        }
                        assert_eq!(
                            area.sign(),
                            if vertex.sense == Sense::Same {
                                Sign::Positive
                            } else {
                                Sign::Negative
                            }
                        );
                    }
                    assert!(
                        crossed_native_seam,
                        "rotated carriers must exercise the unwrapped sphere chart"
                    );
                    (
                        after_strips - before,
                        after_corners - after_strips,
                        boundary_work,
                    )
                });
            assert!(!work.exhausted, "radius={requested}: {work:?}");
            println!(
                "vertex boundaries radius={requested} strip_work={strip_work} corner_work={corner_work} boundary_work={boundary_work} total_control_work={}; not closed-body acceptance",
                work.used
            );
        }
    }
}
