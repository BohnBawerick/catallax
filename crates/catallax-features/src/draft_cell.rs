//! Bounded convex cells of an exact oriented plane arrangement.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Plane, Point3, Vector3};
use std::collections::BTreeSet;

use crate::draft::sign as exact_sign;

/// One outward support and its complete surviving origin set.
#[derive(Clone, Debug)]
pub struct PlanePatch {
    /// Material is on the nonpositive side.
    pub plane: Plane,
    /// Origins of coincident same-sense source patches.
    pub origins: Vec<Origin>,
}

/// A closed, bounded, full-dimensional intersection of oriented halfspaces.
#[derive(Clone, Debug)]
pub struct PlanarCell {
    /// Exact support definitions, including redundant supports.
    pub supports: Vec<PlanePatch>,
    /// Exact triple-support intersections on the boundary.
    pub vertices: Vec<Point3>,
    /// Outward boundary cycles, indexed like supports. Empty means swallowed.
    pub cycles: Vec<Vec<usize>>,
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}

/// Orders a convex face from its first vertex, counterclockwise about the
/// outward normal. Each boundary edge lies on one other support, so the order
/// needs no exact zero proofs. `None` when the supports do not determine it.
fn adjacent_walk(
    support: usize,
    cycle: &[usize],
    incident: &[BTreeSet<usize>],
    vertices: &[Point3],
    patch: &PlanePatch,
    dec: &mut Decider,
) -> Result<Option<Vec<usize>>, Refusal> {
    let neighbours = |a: usize| {
        cycle
            .iter()
            .copied()
            .filter(|b| {
                *b != a
                    && incident[a]
                        .iter()
                        .any(|s| *s != support && incident[*b].contains(s))
            })
            .collect::<Vec<_>>()
    };
    if cycle.iter().any(|v| neighbours(*v).len() != 2) {
        return Ok(None);
    }
    let start = cycle[0];
    let around = neighbours(start);
    let turn = vertices[around[0]]
        .minus(&vertices[start])
        .cross(&vertices[around[1]].minus(&vertices[start]))
        .dot(&patch.plane.normal);
    let mut next = match exact_sign(&turn, dec)? {
        Sign::Positive => around[0],
        Sign::Negative => around[1],
        Sign::Zero => return Ok(None),
    };
    let mut walk = vec![start];
    let mut previous = start;
    while next != start {
        if walk.len() == cycle.len() {
            return Ok(None);
        }
        walk.push(next);
        let [a, b] = <[usize; 2]>::try_from(neighbours(next)).expect("checked above");
        (previous, next) = (next, if a == previous { b } else { a });
    }
    Ok((walk.len() == cycle.len()).then_some(walk))
}

impl PlanarCell {
    /// Constructs the arrangement without assuming the input adjacency survives.
    /// No bounding cutter or approximation contributes a cap.
    pub fn new(supports: Vec<PlanePatch>, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let mut unique: Vec<PlanePatch> = Vec::with_capacity(supports.len());
        for mut support in supports {
            let mut coincident = None;
            for (i, previous) in unique.iter().enumerate() {
                if exact_sign(
                    &previous
                        .plane
                        .normal
                        .cross(&support.plane.normal)
                        .norm_squared(),
                    dec,
                )? == Sign::Zero
                    && exact_sign(&value(&previous.plane, &support.plane.origin), dec)?
                        == Sign::Zero
                {
                    if exact_sign(&previous.plane.normal.dot(&support.plane.normal), dec)?
                        == Sign::Negative
                    {
                        return Ok(None);
                    }
                    coincident = Some(i);
                    break;
                }
            }
            if let Some(i) = coincident {
                unique[i].origins.append(&mut support.origins);
                unique[i].origins.sort();
                unique[i].origins.dedup();
            } else {
                unique.push(support);
            }
        }
        let supports = unique;
        let mut vertices: Vec<Point3> = Vec::new();
        // Supports known to contain each vertex. Their exact zeros need no proof,
        // which an oblique algebraic frame could not afford for every triple.
        let mut incident: Vec<BTreeSet<usize>> = Vec::new();
        for i in 0..supports.len() {
            for j in i + 1..supports.len() {
                for k in j + 1..supports.len() {
                    dec.check_deadline()?;
                    let [a, b, c] = [&supports[i].plane, &supports[j].plane, &supports[k].plane];
                    let bc = b.normal.cross(&c.normal);
                    let determinant = a.normal.dot(&bc);
                    if exact_sign(&determinant, dec)? == Sign::Zero {
                        continue;
                    }
                    let offset =
                        |plane: &Plane| plane.origin.minus(&Point3::origin()).dot(&plane.normal);
                    let vector = bc
                        .scaled(&offset(a))
                        .plus(&c.normal.cross(&a.normal).scaled(&offset(b)))
                        .plus(&a.normal.cross(&b.normal).scaled(&offset(c)))
                        .scaled(&(&Real::one() / &determinant));
                    let point = Point3::origin().plus(&vector);
                    let mut inside = true;
                    for (index, support) in supports.iter().enumerate() {
                        if index != i
                            && index != j
                            && index != k
                            && exact_sign(&value(&support.plane, &point), dec)? == Sign::Positive
                        {
                            inside = false;
                            break;
                        }
                    }
                    if !inside {
                        continue;
                    }
                    let mut duplicate = false;
                    for (index, previous) in vertices.iter().enumerate() {
                        if exact_sign(&point.minus(previous).norm_squared(), dec)? == Sign::Zero {
                            incident[index].extend([i, j, k]);
                            duplicate = true;
                            break;
                        }
                    }
                    if !duplicate {
                        vertices.push(point);
                        incident.push(BTreeSet::from([i, j, k]));
                    }
                }
            }
        }
        if vertices.len() < 4 {
            return Ok(None);
        }
        let mut volume = false;
        let first = &vertices[0];
        'basis: for i in 1..vertices.len() {
            for j in i + 1..vertices.len() {
                let cross = vertices[i].minus(first).cross(&vertices[j].minus(first));
                let shared = |k: usize| {
                    incident[0]
                        .iter()
                        .any(|support| [i, j, k].iter().all(|v| incident[*v].contains(support)))
                };
                for (k, vertex) in vertices.iter().enumerate().skip(j + 1) {
                    if !shared(k)
                        && exact_sign(&cross.dot(&vertex.minus(first)), dec)? != Sign::Zero
                    {
                        volume = true;
                        break 'basis;
                    }
                }
            }
        }
        if !volume {
            return Ok(None);
        }
        // A nontrivial recession cone in R3 has an extreme ray at two support
        // planes, unless their normals fail to span R3. Full-dimensional bounded
        // input gives rank three; the explicit axes also cover rank loss.
        let axes = [
            Vector3::new(1.into(), 0.into(), 0.into()),
            Vector3::new(0.into(), 1.into(), 0.into()),
            Vector3::new(0.into(), 0.into(), 1.into()),
        ];
        let mut directions = axes.to_vec();
        for (i, a) in supports.iter().enumerate() {
            for b in &supports[i + 1..] {
                directions.push(a.plane.normal.cross(&b.plane.normal));
            }
            for axis in &axes {
                directions.push(a.plane.normal.cross(axis));
            }
        }
        for direction in directions {
            if exact_sign(&direction.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            for orientation in [-1, 1] {
                let direction = direction.scaled(&Real::from_i64(orientation));
                let mut recession = true;
                for support in &supports {
                    if exact_sign(&support.plane.normal.dot(&direction), dec)? == Sign::Positive {
                        recession = false;
                        break;
                    }
                }
                if recession {
                    return Err(Refusal::new(
                        "draft.nonclosing-boundary",
                        "drafted halfspaces have an unbounded material direction",
                    ));
                }
            }
        }
        let mut cycles = Vec::with_capacity(supports.len());
        for (index, support) in supports.iter().enumerate() {
            let mut cycle = Vec::new();
            for (i, point) in vertices.iter().enumerate() {
                if incident[i].contains(&index)
                    || exact_sign(&value(&support.plane, point), dec)? == Sign::Zero
                {
                    cycle.push(i);
                }
            }
            if cycle.len() < 3 {
                cycles.push(Vec::new());
                continue;
            }
            let mut area = false;
            let first = &vertices[cycle[0]];
            for i in 1..cycle.len() {
                for j in i + 1..cycle.len() {
                    if exact_sign(
                        &vertices[cycle[i]]
                            .minus(first)
                            .cross(&vertices[cycle[j]].minus(first))
                            .norm_squared(),
                        dec,
                    )? == Sign::Positive
                    {
                        area = true;
                        break;
                    }
                }
                if area {
                    break;
                }
            }
            if !area {
                cycles.push(Vec::new());
                continue;
            }
            if let Some(walk) = adjacent_walk(index, &cycle, &incident, &vertices, support, dec)? {
                cycles.push(walk);
                continue;
            }
            let sum = cycle.iter().fold(Vector3::zero(), |sum, i| {
                sum.plus(&vertices[*i].minus(&Point3::origin()))
            });
            let centre = Point3::origin()
                .plus(&sum.scaled(&(&Real::one() / &Real::from_i64(cycle.len() as i64))));
            let u = vertices[cycle[0]].minus(&centre);
            let v = support.plane.normal.cross(&u);
            let coordinates = |i: usize| {
                let d = vertices[i].minus(&centre);
                [d.dot(&u), d.dot(&v)]
            };
            let half = |p: &[Real; 2], dec: &mut Decider| -> Result<bool, Refusal> {
                Ok(match exact_sign(&p[1], dec)? {
                    Sign::Positive => true,
                    Sign::Negative => false,
                    Sign::Zero => exact_sign(&p[0], dec)? != Sign::Negative,
                })
            };
            // Insertion sort propagates exact predicate failures instead of
            // hiding them in a comparator or choosing a floating-point order.
            for i in 1..cycle.len() {
                let mut j = i;
                while j > 0 {
                    let a = coordinates(cycle[j - 1]);
                    let b = coordinates(cycle[j]);
                    let (ha, hb) = (half(&a, dec)?, half(&b, dec)?);
                    let out_of_order = if ha != hb {
                        !ha
                    } else {
                        exact_sign(&(&a[0] * &b[1] - &a[1] * &b[0]), dec)? == Sign::Negative
                    };
                    if !out_of_order {
                        break;
                    }
                    cycle.swap(j - 1, j);
                    j -= 1;
                }
            }
            cycles.push(cycle);
        }
        Ok(Some(Self {
            supports,
            vertices,
            cycles,
        }))
    }

    /// Strict saved-intent membership. Boundary picks never choose a component.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        for support in &self.supports {
            if exact_sign(&value(&support.plane, point), dec)? != Sign::Negative {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// A relatively open neutral-plane section, including a coincident facet
    /// but excluding a sole edge or vertex contact.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        // A support on the neutral plane already places its vertices there.
        let mut contained = BTreeSet::new();
        for (support, cycle) in self.supports.iter().zip(&self.cycles) {
            if exact_sign(
                &support.plane.normal.cross(&neutral.normal).norm_squared(),
                dec,
            )? == Sign::Zero
                && exact_sign(&value(neutral, &support.plane.origin), dec)? == Sign::Zero
            {
                contained.extend(cycle.iter().copied());
            }
        }
        let mut negative = false;
        let mut positive = false;
        let mut zero = Vec::new();
        for (i, point) in self.vertices.iter().enumerate() {
            if contained.contains(&i) {
                zero.push(i);
                continue;
            }
            match exact_sign(&value(neutral, point), dec)? {
                Sign::Negative => negative = true,
                Sign::Positive => positive = true,
                Sign::Zero => zero.push(i),
            }
        }
        if negative && positive {
            return Ok(true);
        }
        if zero.len() < 3 {
            return Ok(false);
        }
        let p = &self.vertices[zero[0]];
        for i in 1..zero.len() {
            for j in i + 1..zero.len() {
                if exact_sign(
                    &self.vertices[zero[i]]
                        .minus(p)
                        .cross(&self.vertices[zero[j]].minus(p))
                        .norm_squared(),
                    dec,
                )? == Sign::Positive
                {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}
