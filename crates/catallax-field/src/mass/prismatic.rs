//! Analytic CSG moments of axis boxes and finite right circular cylinders.
//!
//! Indicator algebra distributes the Boolean tape into signed intersections.
//! Bounds remove only zero-volume intersections. Every remaining intersection
//! must contain parallel cylinders, so it is a finite extrusion of a rectangle
//! clipped by disks. Green's theorem integrates its complete line/arc boundary.

use super::{
    Accuracy, Bounds, Contribution, Options, PAIRS, Properties, Totals, axis_block, failure, graph,
    q, reached, valid,
};
use crate::{Instruction, Plan, geometry::Leaf, metric::Kind};
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_recipe::BooleanKind;
use catallax_shapes::{RPoint3, Vector3};
use std::collections::BTreeMap;

mod section;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Circle {
    axis: usize,
    radius: Rational,
    centre: [Rational; 2],
}
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Shape {
    bounds: Bounds,
    circles: Vec<Circle>,
}
type Polynomial = BTreeMap<Shape, Rational>;

struct Counter {
    used: usize,
    limit: usize,
}
impl Counter {
    fn step(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        if self.used == self.limit {
            return Err(failure(
                "field.mass-budget",
                format!(
                    "{} CSG intersection/boundary cells exhaust limit {}",
                    self.used, self.limit
                ),
            ));
        }
        self.used += 1;
        dec.charge_work(1)
    }
}
fn coverage(detail: impl Into<String>) -> Refusal {
    failure("field.mass-prismatic-coverage", detail)
}
fn axes(axis: usize) -> [usize; 2] {
    [(axis + 1) % 3, (axis + 2) % 3]
}

fn leaf(leaf: &Leaf, dec: &mut Decider) -> Result<Shape, Refusal> {
    dec.charge_work(1)?;
    match leaf {
        Leaf::Primitive(p) => {
            if let Some(bounds) = axis_block(p) {
                return Ok(Shape {
                    bounds,
                    circles: Vec::new(),
                });
            }
            if p.kind != Kind::Round
                || p.parameters[0] != p.parameters[1]
                || p.parameters[0] <= q(0)
                || p.parameters[2] <= q(0)
            {
                return Err(coverage(
                    "requires an axis box or a finite right circular cylinder",
                ));
            }
            let world = p.inverse_placement().inverse();
            let origin = world
                .apply_rational_point(&RPoint3::origin())
                .ok_or_else(|| coverage("cylinder origin is not rational"))?;
            let direction = world
                .apply_vector(&Vector3::new(0.into(), 0.into(), 1.into()))
                .as_rational()
                .ok_or_else(|| coverage("cylinder axis is not rational"))?;
            let direction = [direction.x, direction.y, direction.z];
            let mut nonzero = direction.iter().enumerate().filter(|(_, v)| !v.is_zero());
            let Some((axis, sign)) = nonzero.next() else {
                return Err(coverage("cylinder axis is zero"));
            };
            if nonzero.next().is_some() || sign.abs() != q(1) {
                return Err(coverage("cylinder axis must be a coordinate axis"));
            }
            let origin = [origin.x, origin.y, origin.z];
            let radius = p.parameters[0].clone();
            let end = &origin[axis] + &(sign * &p.parameters[2]);
            let mut bounds = std::array::from_fn(|i| [&origin[i] - &radius, &origin[i] + &radius]);
            bounds[axis] = [
                origin[axis].clone().min(end.clone()),
                origin[axis].clone().max(end),
            ];
            let [u, v] = axes(axis);
            Ok(Shape {
                bounds,
                circles: vec![Circle {
                    axis,
                    centre: [origin[u].clone(), origin[v].clone()],
                    radius,
                }],
            })
        }
        Leaf::Feature(feature) => graph::prism(feature, dec)?
            .map(|bounds| Shape {
                bounds,
                circles: Vec::new(),
            })
            .ok_or_else(|| coverage("feature must prove it is a bounded axis box")),
        Leaf::Cell(cell) => cell
            .convex_planes()
            .and_then(graph::axis_planes)
            .map(|bounds| Shape {
                bounds,
                circles: Vec::new(),
            })
            .ok_or_else(|| coverage("cell must prove it is a bounded axis box")),
        _ => Err(coverage("live leaf is not a covered box or cylinder")),
    }
}

impl Shape {
    fn intersection(&self, other: &Self) -> Option<Self> {
        let bounds: Bounds = std::array::from_fn(|i| {
            [
                self.bounds[i][0].clone().max(other.bounds[i][0].clone()),
                self.bounds[i][1].clone().min(other.bounds[i][1].clone()),
            ]
        });
        if bounds.iter().any(|v| v[0] >= v[1]) {
            return None;
        }
        let mut circles = self.circles.clone();
        circles.extend(other.circles.iter().cloned());
        circles.sort();
        circles.dedup();
        // Equal disks centred on a line: the two extreme centres imply every
        // intermediate constraint by convexity of squared distance in the centre.
        let mut normalized = Vec::new();
        let mut start = 0;
        while start < circles.len() {
            let mut end = start + 1;
            while end < circles.len()
                && circles[end].axis == circles[start].axis
                && circles[end].radius == circles[start].radius
            {
                end += 1;
            }
            let group = &circles[start..end];
            if group.len() > 2 {
                let d = [
                    &group[1].centre[0] - &group[0].centre[0],
                    &group[1].centre[1] - &group[0].centre[1],
                ];
                let collinear = group.iter().all(|c| {
                    (&(&c.centre[0] - &group[0].centre[0]) * &d[1])
                        == (&(&c.centre[1] - &group[0].centre[1]) * &d[0])
                });
                if collinear {
                    let projection = |c: &Circle| &(&c.centre[0] * &d[0]) + &(&c.centre[1] * &d[1]);
                    normalized.push(group.iter().min_by_key(|c| projection(c)).unwrap().clone());
                    normalized.push(group.iter().max_by_key(|c| projection(c)).unwrap().clone());
                    start = end;
                    continue;
                }
            }
            normalized.extend_from_slice(group);
            start = end;
        }
        normalized.sort();
        Some(Self {
            bounds,
            circles: normalized,
        })
    }
}

fn add(polynomial: &mut Polynomial, shape: Shape, coefficient: Rational) {
    let value = polynomial.entry(shape.clone()).or_insert_with(|| q(0));
    *value = &*value + &coefficient;
    if value.is_zero() {
        polynomial.remove(&shape);
    }
}
fn product(
    a: &Polynomial,
    b: &Polynomial,
    counter: &mut Counter,
    dec: &mut Decider,
) -> Result<Polynomial, Refusal> {
    let mut result = Polynomial::new();
    for (a, ca) in a {
        for (b, cb) in b {
            counter.step(dec)?;
            if let Some(shape) = a.intersection(b) {
                add(&mut result, shape, ca * cb);
            }
        }
    }
    Ok(result)
}
fn polynomial(
    plan: &Plan,
    counter: &mut Counter,
    dec: &mut Decider,
) -> Result<Polynomial, Refusal> {
    let mut leaves = Vec::with_capacity(plan.leaves.len());
    for (i, definition) in plan.leaves.iter().enumerate() {
        leaves.push(if plan.leaf_registers[i].is_some() {
            Some(leaf(definition, dec)?)
        } else {
            None
        });
    }
    let mut tape: Vec<Polynomial> = Vec::with_capacity(plan.instructions.len());
    for instruction in &plan.instructions {
        counter.step(dec)?;
        let result = match *instruction {
            Instruction::Empty => Polynomial::new(),
            Instruction::Leaf(i) => {
                let shape = leaves[i].as_ref().expect("live field instruction");
                if shape.bounds.iter().any(|v| v[0] >= v[1]) {
                    Polynomial::new()
                } else {
                    BTreeMap::from([(shape.clone(), q(1))])
                }
            }
            Instruction::Boolean { kind, a, b } => {
                let intersection = product(&tape[a], &tape[b], counter, dec)?;
                match kind {
                    BooleanKind::Intersect => intersection,
                    BooleanKind::Union | BooleanKind::Subtract => {
                        let mut result = tape[a].clone();
                        if kind == BooleanKind::Union {
                            for (shape, value) in &tape[b] {
                                add(&mut result, shape.clone(), value.clone());
                            }
                        }
                        for (shape, value) in intersection {
                            add(&mut result, shape, -value);
                        }
                        result
                    }
                }
            }
        };
        tape.push(result);
    }
    Ok(tape.swap_remove(plan.root))
}

fn power_integral(limits: &[Rational; 2], power: usize) -> Rational {
    let n = (power + 1) as i32;
    (limits[1].pow(n).unwrap() - limits[0].pow(n).unwrap())
        .checked_div(&q(n as i64))
        .unwrap()
}
fn moment_exponents(index: usize) -> [usize; 3] {
    let mut result = [0; 3];
    if (1..=3).contains(&index) {
        result[index - 1] = 1;
    }
    if index >= 4 {
        let (i, j) = PAIRS[index - 4];
        result[i] += 1;
        result[j] += 1;
    }
    result
}
fn box_moments(shape: &Shape) -> [Rational; 10] {
    std::array::from_fn(|index| {
        let exponents = moment_exponents(index);
        (0..3).fold(q(1), |value, i| {
            value * power_integral(&shape.bounds[i], exponents[i])
        })
    })
}
fn moments(
    shape: &Shape,
    bits: u32,
    counter: &mut Counter,
    dec: &mut Decider,
) -> Result<[Ball; 10], Refusal> {
    if shape.circles.is_empty() {
        return Ok(box_moments(shape)
            .each_ref()
            .map(|value| Ball::from_rational(value, bits)));
    }
    let axis = shape.circles[0].axis;
    if shape.circles.iter().any(|c| c.axis != axis) {
        return Err(coverage(format!(
            "intersection has nonparallel cylinder axes; enclosing box {:?}",
            shape.bounds
        )));
    }
    let [u, v] = axes(axis);
    let radial = section::moments(
        &shape.circles,
        &[shape.bounds[u].clone(), shape.bounds[v].clone()],
        bits,
        counter,
        dec,
    )?;
    const POWERS: [[usize; 2]; 6] = [[0, 0], [1, 0], [0, 1], [2, 0], [0, 2], [1, 1]];
    Ok(std::array::from_fn(|index| {
        let exponents = moment_exponents(index);
        let pair = [exponents[u], exponents[v]];
        let index = POWERS.iter().position(|p| *p == pair).unwrap();
        &radial[index]
            * &Ball::from_rational(&power_integral(&shape.bounds[axis], exponents[axis]), bits)
    }))
}

/// Certified moments of Boolean combinations of axis boxes and finite right
/// cylinders. Every nonzero intersection term must have parallel cylinders;
/// empty box intersections are proved zero before this coverage check.
/// The original field tape supplies the set algebra. Complete analytic arc and
/// line boundaries supply all moments, independently of B-rep and print meshes.
/// This opt-in entry leaves the existing `properties` route unchanged.
pub fn properties(
    plan: &Plan,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Result<Properties, Refusal> {
    if options.max_cells == 0
        || ![&accuracy.volume, &accuracy.centroid, &accuracy.inertia]
            .into_iter()
            .all(valid)
    {
        return Err(failure(
            "field.mass-options",
            "nonnegative accuracies and a positive cell limit are required",
        ));
    }
    let mut counter = Counter {
        used: 0,
        limit: options.max_cells,
    };
    let polynomial = polynomial(plan, &mut counter, dec)?;
    if polynomial.is_empty() {
        return Ok(Properties {
            volume: super::exact(q(0)),
            centroid: None,
            inertia: None,
            estimate: None,
            cells: counter.used,
        });
    }
    if polynomial.keys().all(|shape| shape.circles.is_empty()) {
        let mut sums: [Rational; 10] = std::array::from_fn(|_| q(0));
        for (shape, coefficient) in &polynomial {
            counter.step(dec)?;
            let moments = box_moments(shape);
            for i in 0..10 {
                sums[i] = &sums[i] + &(&moments[i] * coefficient);
            }
        }
        let totals = Totals {
            lo: sums.clone(),
            hi: sums.clone(),
            estimate: sums,
        };
        return totals
            .properties(&[q(0), q(0), q(0)], counter.used)
            .ok_or_else(|| {
                failure(
                    "field.mass-volume",
                    format!("exact box CSG has inconsistent volume {}", totals.lo[0]),
                )
            });
    }
    let maximum = dec.policy().max_precision_bits;
    if maximum == 0 {
        return Err(failure(
            "numbers.precision-limit",
            "no prismatic moment precision available",
        ));
    }
    let mut bits = maximum.min(128);
    loop {
        let mut sums: [Ball; 10] = std::array::from_fn(|_| Ball::from_rational(&q(0), bits));
        let result = (|| -> Result<_, Refusal> {
            for (shape, coefficient) in &polynomial {
                counter.step(dec)?;
                let moments = moments(shape, bits, &mut counter, dec)?;
                let coefficient = Ball::from_rational(coefficient, bits);
                for i in 0..10 {
                    sums[i] = &sums[i] + &(&moments[i] * &coefficient);
                }
            }
            let contribution = Contribution {
                lo: sums.each_ref().map(Ball::lo),
                hi: sums.each_ref().map(Ball::hi),
                estimate: sums.each_ref().map(|v| (&v.lo() + &v.hi()).mul_pow2(-1)),
            };
            let mut totals = Totals::new();
            totals.update(&contribution, false);
            Ok(totals.properties(&[q(0), q(0), q(0)], counter.used))
        })();
        let mut refinement = None;
        match result {
            Ok(Some(result))
                if reached(&result.volume, &accuracy.volume)
                    && result
                        .centroid
                        .as_ref()
                        .is_none_or(|c| c.iter().all(|v| reached(v, &accuracy.centroid)))
                    && result.inertia.as_ref().is_none_or(|c| {
                        c.iter().flatten().all(|v| reached(v, &accuracy.inertia))
                    }) =>
            {
                return Ok(result);
            }
            Err(r) if r.code != "field.mass-prismatic-refine" => return Err(r),
            Err(r) => refinement = Some(r),
            _ => {}
        }
        if bits == maximum {
            if let Some(refinement) = refinement {
                let mut refusal = refinement.with_code("field.mass-accuracy");
                refusal.detail = format!(
                    "{} CSG/boundary cells at {bits} bits; {}",
                    counter.used, refusal.detail
                );
                return Err(refusal);
            }
            return Err(failure(
                "field.mass-accuracy",
                format!(
                    "{} CSG/boundary cells at {bits} bits; reached volume [{},{}] mm3; requested property widths or positive volume not certified",
                    counter.used,
                    sums[0].lo().to_decimal_bounds(10).0,
                    sums[0].hi().to_decimal_bounds(10).1
                ),
            ));
        }
        bits = bits.saturating_mul(2).min(maximum);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Cell;
    use catallax_recipe::Recipe;
    use catallax_shapes::{Plane, RVector3, Surface};

    #[test]
    fn contradictory_axis_halfspaces_never_subtract_negative_material() {
        let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"empty cell subtraction","steps":[{"id":"a","op":"block","size":["2","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"b","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"c","op":"boolean","kind":"subtract","target":"a","tool":"b"}],"result":"c"}"#).unwrap();
        for upper_x in [0, 1] {
            let mut plan = Plan::compile(&recipe).unwrap();
            let bounds = [[q(1), q(upper_x)], [q(0), q(1)], [q(0), q(1)]];
            let mut planes = Vec::new();
            for axis in 0..3 {
                for side in 0..2 {
                    let mut origin = [q(0), q(0), q(0)];
                    origin[axis] = bounds[axis][side].clone();
                    let normal: [Rational; 3] = std::array::from_fn(|i| {
                        if i == axis {
                            q(if side == 0 { -1 } else { 1 })
                        } else {
                            q(0)
                        }
                    });
                    let reference: [Rational; 3] =
                        std::array::from_fn(|i| q(if i == (axis + 1) % 3 { 1 } else { 0 }));
                    planes.push(Surface::Plane(Plane {
                        origin: RPoint3::new(
                            origin[0].clone(),
                            origin[1].clone(),
                            origin[2].clone(),
                        )
                        .to_real(),
                        normal: RVector3::new(
                            normal[0].clone(),
                            normal[1].clone(),
                            normal[2].clone(),
                        )
                        .to_real(),
                        ref_dir: RVector3::new(
                            reference[0].clone(),
                            reference[1].clone(),
                            reference[2].clone(),
                        )
                        .to_real(),
                    }));
                }
            }
            *plan.leaves.last_mut().unwrap() = Leaf::Cell(Cell::new(planes, None).unwrap());
            let result = properties(
                &plan,
                &Accuracy::default(),
                &Options::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert_eq!(result.volume.lo, q(8));
            assert_eq!(result.volume.hi, q(8));
        }
    }
}
