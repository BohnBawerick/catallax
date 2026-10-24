//! Constant-cubic fields reduce to analytic columns on nearest-grid sectors.
//!
//! On an octant sector 0 <= u <= v <= a/2, the cubic leaf is the union
//! of vertical cylinders u²+v² <= r² and horizontal intervals
//! |z-z_grid| <= sqrt(r²-u²). Boxes and parallel cylinders add only
//! constant or circular column endpoints. Certified endpoint ordering makes
//! every inner integral polynomial. Only the outer coordinate needs cubature.
//! Taylor's theorem bounds its remainder; unresolved endpoint order retains
//! a material-density enclosure instead of assuming a topology.

use super::{
    Accuracy, Bounds, Cell, Contribution, Options, Properties, Totals, axis_block, failure, q,
    reached,
};
use crate::{Instruction, Plan, geometry, metric::Kind};
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use std::{
    cmp::Ordering,
    collections::{BTreeSet, BinaryHeap},
};

const PRECISION: u32 = 128;
const TERMS: usize = 10;
const EXPONENTS: [[usize; 3]; 10] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
];
fn ball(value: &Rational) -> Ball {
    Ball::from_rational(value, PRECISION)
}
fn integer(value: i64) -> Ball {
    ball(&q(value))
}
fn cover(lo: &Rational, hi: &Rational) -> Ball {
    ball(lo).hull(&ball(hi))
}
fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

/// Coefficients f^(k)/k! at a point or over an interval, not a sampled fit.
#[derive(Clone)]
struct Jet([Ball; TERMS]);
impl Jet {
    fn constant(value: Ball) -> Self {
        Self(std::array::from_fn(|i| {
            if i == 0 { value.clone() } else { integer(0) }
        }))
    }
    fn variable(value: Ball) -> Self {
        let mut jet = Self::constant(value);
        jet.0[1] = integer(1);
        jet
    }
    fn add(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] + &other.0[i]))
    }
    fn sub(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] - &other.0[i]))
    }
    fn scale(&self, scale: &Rational) -> Self {
        let scale = ball(scale);
        Self(std::array::from_fn(|i| &self.0[i] * &scale))
    }
    fn mul(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|n| {
            let mut sum = integer(0);
            for k in 0..=n {
                if !self.0[k].is_exact_zero() && !other.0[n - k].is_exact_zero() {
                    sum = &sum + &(&self.0[k] * &other.0[n - k]);
                }
            }
            sum
        }))
    }
    fn sqrt(&self) -> Option<Self> {
        let mut result = Self::constant(self.0[0].sqrt()?);
        if self.0[1..].iter().all(Ball::is_exact_zero) {
            return Some(result);
        }
        if self.0[0].lo() <= q(0) {
            return None;
        }
        let denominator = &integer(2) * &result.0[0];
        for n in 1..TERMS {
            let mut sum = integer(0);
            for k in 1..n {
                sum = &sum + &(&result.0[k] * &result.0[n - k]);
            }
            result.0[n] = (&self.0[n] - &sum).checked_div(&denominator)?;
        }
        Some(result)
    }
}

#[derive(Clone, PartialEq, Eq)]
enum Curve {
    Constant(Rational),
    Variable,
    Root {
        offset: Rational,
        shift: Rational,
        radius: Rational,
        sign: i64,
    },
}
impl Curve {
    fn evaluate(&self, u: &Jet) -> Option<Jet> {
        match self {
            Self::Constant(c) => Some(Jet::constant(ball(c))),
            Self::Variable => Some(u.clone()),
            Self::Root {
                offset,
                shift,
                radius,
                sign,
            } => {
                let shifted = u.sub(&Jet::constant(ball(shift)));
                let radical = Jet::constant(ball(&(radius * radius))).sub(&shifted.mul(&shifted));
                Some(Jet::constant(ball(offset)).add(&radical.sqrt()?.scale(&q(*sign))))
            }
        }
    }
}
fn cut(curves: &mut Vec<Curve>, curve: Curve) -> usize {
    if let Some(i) = curves.iter().position(|old| *old == curve) {
        i
    } else {
        curves.push(curve);
        curves.len() - 1
    }
}
fn constant(curves: &mut Vec<Curve>, value: &Rational) -> usize {
    cut(curves, Curve::Constant(value.clone()))
}

#[derive(Clone)]
enum Definition {
    Empty,
    Box(Bounds),
    Cylinder {
        center: [Rational; 2],
        radius: Rational,
        z: [Rational; 2],
    },
    Lattice,
}
struct Model<'a> {
    plan: &'a Plan,
    definitions: Vec<Definition>,
    grid: [Rational; 3],
    size: Rational,
    radius: Rational,
    bounds: Bounds,
    anchor: [Rational; 3],
    z_nodes: Vec<Rational>,
}
fn cylinder(primitive: &geometry::Primitive) -> Option<Definition> {
    if primitive.kind != Kind::Round || primitive.parameters[0] != primitive.parameters[1] {
        return None;
    }
    let map = primitive.inverse_placement();
    let m = map.rational_matrix()?;
    if !m[2][0].is_zero() || !m[2][1].is_zero() || m[2][2].abs() != q(1) {
        return None;
    }
    let translation = map.translation.as_rational()?;
    let t = [translation.x, translation.y, translation.z];
    let center: [Rational; 2] =
        std::array::from_fn(|i| -(0..3).fold(q(0), |sum, j| &sum + &(m[j][i] * &t[j])));
    let z = if *m[2][2] == q(1) {
        [-&t[2], &primitive.parameters[2] - &t[2]]
    } else {
        [&t[2] - &primitive.parameters[2], t[2].clone()]
    };
    Some(Definition::Cylinder {
        center,
        radius: primitive.parameters[0].clone(),
        z,
    })
}
impl<'a> Model<'a> {
    fn new(plan: &'a Plan, bounds: &Bounds) -> Option<Self> {
        let mut grid = None;
        let mut definitions = Vec::with_capacity(plan.leaves.len());
        for (i, leaf) in plan.leaves.iter().enumerate() {
            if plan.leaf_registers[i].is_none() {
                definitions.push(Definition::Empty);
                continue;
            }
            definitions.push(match leaf {
                geometry::Leaf::Primitive(p) => {
                    if let Some(bounds) = axis_block(p) {
                        Definition::Box(bounds)
                    } else {
                        cylinder(p)?
                    }
                }
                geometry::Leaf::Lattice(lattice) => {
                    let candidate = lattice.axis_cubic_grid()?;
                    if grid.as_ref().is_some_and(|old| *old != candidate) {
                        return None;
                    }
                    grid = Some(candidate);
                    Definition::Lattice
                }
                _ => return None,
            });
        }
        let (grid, size, radius) = grid?;
        // Larger radii couple adjacent grid lines. Keep the general certified
        // integrator for that geometry instead of dropping overlapping cells.
        if radius.mul_pow2(1) >= size {
            return None;
        }
        Some(Self {
            plan,
            definitions,
            grid,
            size,
            radius,
            bounds: bounds.clone(),
            anchor: std::array::from_fn(|i| (&bounds[i][0] + &bounds[i][1]).mul_pow2(-1)),
            z_nodes: Vec::new(),
        })
    }
    fn truth(&self, leaves: &[bool], registers: &mut [bool]) -> bool {
        for (i, instruction) in self.plan.instructions.iter().enumerate() {
            registers[i] = match *instruction {
                Instruction::Empty => false,
                Instruction::Leaf(leaf) => leaves[leaf],
                Instruction::Boolean { kind, a, b } => match kind {
                    catallax_recipe::BooleanKind::Union => registers[a] || registers[b],
                    catallax_recipe::BooleanKind::Intersect => registers[a] && registers[b],
                    catallax_recipe::BooleanKind::Subtract => registers[a] && !registers[b],
                },
            };
        }
        registers[self.plan.root]
    }
    fn nodes(&self, axis: usize, max_nodes: usize) -> Option<Vec<Rational>> {
        let half = self.size.mul_pow2(-1);
        let first = (&(&(&self.bounds[axis][0] - &self.grid[axis]) - &half) / &self.size)
            .ceil()
            .to_i64()?;
        let last = (&(&(&self.bounds[axis][1] - &self.grid[axis]) + &half) / &self.size)
            .floor()
            .to_i64()?;
        if last.saturating_sub(first) as u64 >= max_nodes as u64 {
            return None;
        }
        Some(
            (first..=last)
                .map(|k| &self.grid[axis] + &(&self.size * &q(k)))
                .collect(),
        )
    }
}
#[derive(Clone)]
struct Sector {
    node: [Rational; 2],
    outer: usize,
    signs: [i64; 2],
}
impl Sector {
    fn inner(&self) -> usize {
        1 - self.outer
    }
    fn position(&self, axis: usize, distance: &Rational) -> Rational {
        &self.node[axis] + &(&q(self.signs[axis]) * distance)
    }
    fn coordinate(&self, axis: usize, position: &Rational) -> Rational {
        &q(self.signs[axis]) * &(position - &self.node[axis])
    }
}
#[derive(Clone)]
enum Region {
    Empty,
    Product {
        y: [usize; 2],
        z: [usize; 2],
    },
    Lattice {
        vertical: usize,
        horizontal: Vec<[usize; 2]>,
    },
}
struct OrderingProof {
    curves: Vec<Curve>,
    ranks: Vec<usize>,
}
fn ordered(curves: &[Curve], lower: &Curve, upper: &Curve, u: &Ball) -> Option<OrderingProof> {
    let variable = Jet::constant(u.clone());
    let lo = lower.evaluate(&variable)?.0[0].clone();
    let hi = upper.evaluate(&variable)?.0[0].clone();
    let mut clipped = vec![lower.clone(), upper.clone()];
    let mut aliases = Vec::with_capacity(curves.len());
    for curve in curves {
        let value = curve.evaluate(&variable)?.0[0].clone();
        let alias = if value.hi() <= lo.lo() {
            0
        } else if value.lo() >= hi.hi() {
            1
        } else {
            cut(&mut clipped, curve.clone())
        };
        aliases.push(alias);
    }
    let values: Vec<Ball> = clipped
        .iter()
        .map(|c| c.evaluate(&variable).map(|j| j.0[0].clone()))
        .collect::<Option<_>>()?;
    let mut indices: Vec<usize> = (0..clipped.len()).collect();
    indices.sort_by_key(|&i| values[i].mid());
    for pair in indices.windows(2) {
        if values[pair[0]].hi() > values[pair[1]].lo() {
            return None;
        }
    }
    if indices.first() != Some(&0) || indices.last() != Some(&1) {
        return None;
    }
    let mut ranks = vec![0; clipped.len()];
    for (rank, &i) in indices.iter().enumerate() {
        ranks[i] = rank;
    }
    Some(OrderingProof {
        curves: indices.into_iter().map(|i| clipped[i].clone()).collect(),
        ranks: aliases.into_iter().map(|i| ranks[i]).collect(),
    })
}
struct Band {
    y: [usize; 2],
    z: Vec<[usize; 2]>,
}
struct Atlas {
    y: Vec<Curve>,
    z: Vec<Curve>,
    bands: Vec<Band>,
}
impl Atlas {
    fn new(
        model: &Model<'_>,
        sector: &Sector,
        lo: &Rational,
        hi: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let u = cover(lo, hi);
        let outer = sector.outer;
        let inner = sector.inner();
        let half = model.size.mul_pow2(-1);
        let mut y = Vec::new();
        let mut z = Vec::new();
        let mut regions = Vec::with_capacity(model.definitions.len());
        for definition in &model.definitions {
            dec.check_deadline()?;
            regions.push(match definition {
                Definition::Empty => Region::Empty,
                Definition::Box(bounds) => {
                    let a = sector.position(outer, lo);
                    let b = sector.position(outer, hi);
                    let low = a.clone().min(b.clone());
                    let high = a.max(b);
                    if high <= bounds[outer][0] || low >= bounds[outer][1] {
                        Region::Empty
                    } else if low < bounds[outer][0] || high > bounds[outer][1] {
                        return Ok(None);
                    } else {
                        let a = sector.coordinate(inner, &bounds[inner][0]);
                        let b = sector.coordinate(inner, &bounds[inner][1]);
                        let yc = [
                            constant(&mut y, &a.clone().min(b.clone())),
                            constant(&mut y, &a.max(b)),
                        ];
                        let zc = [
                            constant(&mut z, &bounds[2][0]),
                            constant(&mut z, &bounds[2][1]),
                        ];
                        Region::Product { y: yc, z: zc }
                    }
                }
                Definition::Cylinder {
                    center,
                    radius,
                    z: extent,
                } => {
                    let shift = sector.coordinate(outer, &center[outer]);
                    let delta = &u - &ball(&shift);
                    let radicand = &ball(&(radius * radius)) - &delta.square();
                    if radicand.hi() <= q(0) {
                        Region::Empty
                    } else if radicand.lo() < q(0) {
                        return Ok(None);
                    } else {
                        let offset = sector.coordinate(inner, &center[inner]);
                        let yc = [-1, 1].map(|sign| {
                            cut(
                                &mut y,
                                Curve::Root {
                                    offset: offset.clone(),
                                    shift: shift.clone(),
                                    radius: radius.clone(),
                                    sign,
                                },
                            )
                        });
                        let zc = [constant(&mut z, &extent[0]), constant(&mut z, &extent[1])];
                        Region::Product { y: yc, z: zc }
                    }
                }
                Definition::Lattice => {
                    if *lo >= model.radius {
                        Region::Empty
                    } else if *hi > model.radius {
                        return Ok(None);
                    } else {
                        let vertical = cut(
                            &mut y,
                            Curve::Root {
                                offset: q(0),
                                shift: q(0),
                                radius: model.radius.clone(),
                                sign: 1,
                            },
                        );
                        let mut horizontal = Vec::new();
                        for node in &model.z_nodes {
                            horizontal.push([-1, 1].map(|sign| {
                                cut(
                                    &mut z,
                                    Curve::Root {
                                        offset: node.clone(),
                                        shift: q(0),
                                        radius: model.radius.clone(),
                                        sign,
                                    },
                                )
                            }));
                        }
                        Region::Lattice {
                            vertical,
                            horizontal,
                        }
                    }
                }
            });
        }
        let Some(yp) = ordered(&y, &Curve::Variable, &Curve::Constant(half), &u) else {
            return Ok(None);
        };
        let Some(zp) = ordered(
            &z,
            &Curve::Constant(model.bounds[2][0].clone()),
            &Curve::Constant(model.bounds[2][1].clone()),
            &u,
        ) else {
            return Ok(None);
        };
        for region in &mut regions {
            match region {
                Region::Product { y, z } => {
                    *y = y.map(|i| yp.ranks[i]);
                    *z = z.map(|i| zp.ranks[i]);
                }
                Region::Lattice {
                    vertical,
                    horizontal,
                } => {
                    *vertical = yp.ranks[*vertical];
                    for z in horizontal {
                        *z = z.map(|i| zp.ranks[i]);
                    }
                }
                Region::Empty => {}
            }
        }
        let mut leaves = vec![false; model.definitions.len()];
        let mut registers = vec![false; model.plan.instructions.len()];
        let mut admitted = vec![Vec::new(); yp.curves.len() - 1];
        for (yi, bands) in admitted.iter_mut().enumerate() {
            for zi in 0..zp.curves.len() - 1 {
                dec.check_deadline()?;
                for (i, region) in regions.iter().enumerate() {
                    leaves[i] = match region {
                        Region::Empty => false,
                        Region::Product { y, z } => {
                            y[0] <= yi && yi < y[1] && z[0] <= zi && zi < z[1]
                        }
                        Region::Lattice {
                            vertical,
                            horizontal,
                        } => yi < *vertical || horizontal.iter().any(|z| z[0] <= zi && zi < z[1]),
                    };
                }
                if model.truth(&leaves, &mut registers) {
                    bands.push(zi);
                }
            }
        }
        // Cancel interior subdivision walls before interval differentiation.
        // Evaluating F(b)-F(a)+F(c)-F(b) as intervals would retain a false
        // derivative uncertainty from b, even when b is not a field boundary.
        let mut bands: Vec<Band> = Vec::new();
        for (yi, columns) in admitted.into_iter().enumerate() {
            let mut z: Vec<[usize; 2]> = Vec::new();
            for zi in columns {
                if let Some(last) = z.last_mut()
                    && last[1] == zi
                {
                    last[1] = zi + 1;
                } else {
                    z.push([zi, zi + 1]);
                }
            }
            if z.is_empty() {
                continue;
            }
            if let Some(last) = bands.last_mut()
                && last.y[1] == yi
                && last.z == z
            {
                last.y[1] = yi + 1;
            } else {
                bands.push(Band { y: [yi, yi + 1], z });
            }
        }
        Ok(Some(Self {
            y: yp.curves,
            z: zp.curves,
            bands,
        }))
    }
    fn evaluate(&self, model: &Model<'_>, sector: &Sector, u: &Jet) -> Option<[Jet; 10]> {
        if self.bands.is_empty() {
            return Some(std::array::from_fn(|_| Jet::constant(integer(0))));
        }
        let inner = sector.inner();
        let offset = &sector.node[inner] - &model.anchor[inner];
        let sign = q(sector.signs[inner]);
        let mut yvalues: Vec<Option<[Jet; 3]>> = vec![None; self.y.len()];
        let mut zvalues: Vec<Option<[Jet; 3]>> = vec![None; self.z.len()];
        for band in &self.bands {
            for i in band.y {
                if yvalues[i].is_none() {
                    let value = self.y[i]
                        .evaluate(u)?
                        .scale(&sign)
                        .add(&Jet::constant(ball(&offset)));
                    yvalues[i] = Some(antiderivatives(value));
                }
            }
            for endpoints in &band.z {
                for &i in endpoints {
                    if zvalues[i].is_none() {
                        let value = self.z[i]
                            .evaluate(u)?
                            .sub(&Jet::constant(ball(&model.anchor[2])));
                        zvalues[i] = Some(antiderivatives(value));
                    }
                }
            }
        }
        let mut yz: [[Jet; 3]; 3] =
            std::array::from_fn(|_| std::array::from_fn(|_| Jet::constant(integer(0))));
        for band in &self.bands {
            let lower = yvalues[band.y[0]].as_ref().expect("admitted endpoint");
            let upper = yvalues[band.y[1]].as_ref().expect("admitted endpoint");
            let yints: [Jet; 3] = std::array::from_fn(|i| upper[i].sub(&lower[i]).scale(&sign));
            let mut zsum: [Jet; 3] = std::array::from_fn(|_| Jet::constant(integer(0)));
            for endpoints in &band.z {
                let lower = zvalues[endpoints[0]].as_ref().expect("admitted endpoint");
                let upper = zvalues[endpoints[1]].as_ref().expect("admitted endpoint");
                for (k, sum) in zsum.iter_mut().enumerate() {
                    *sum = sum.add(&upper[k].sub(&lower[k]));
                }
            }
            for (a, row) in yz.iter_mut().enumerate() {
                for (b, value) in row.iter_mut().enumerate().take(3 - a) {
                    *value = value.add(&yints[a].mul(&zsum[b]));
                }
            }
        }
        let x = Jet::constant(ball(
            &(&sector.node[sector.outer] - &model.anchor[sector.outer]),
        ))
        .add(&u.scale(&q(sector.signs[sector.outer])));
        let xp = [Jet::constant(integer(1)), x.clone(), x.mul(&x)];
        Some(std::array::from_fn(|i| {
            let ex = EXPONENTS[i];
            xp[ex[sector.outer]].mul(&yz[ex[inner]][ex[2]])
        }))
    }
}

fn antiderivatives(value: Jet) -> [Jet; 3] {
    let square = value.mul(&value);
    let cube = square.mul(&value);
    [
        value,
        square.scale(&q(1).mul_pow2(-1)),
        cube.scale(&q(1).checked_div(&q(3)).unwrap()),
    ]
}

struct Pending {
    sector: usize,
    lo: Rational,
    hi: Rational,
    contribution: Contribution,
    priority: Rational,
    serial: usize,
}
impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Pending {
    fn cmp(&self, other: &Self) -> Ordering {
        self.priority
            .cmp(&other.priority)
            .then_with(|| other.serial.cmp(&self.serial))
    }
}
struct Integration<'a> {
    model: Model<'a>,
    sectors: Vec<Sector>,
    totals: Totals,
    pending: BinaryHeap<Pending>,
    cells: usize,
    max_cells: usize,
}
impl Integration<'_> {
    fn admit(
        &mut self,
        sector: usize,
        lo: Rational,
        hi: Rational,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if self.cells == self.max_cells {
            return Err(failure(
                "field.mass-budget",
                "analytic column cell limit reached before requested Accuracy",
            ));
        }
        dec.check_deadline()?;
        let index = self.cells;
        self.cells += 1;
        let geometry = &self.sectors[sector];
        let mut bounds = self.model.bounds.clone();
        let a = geometry.position(geometry.outer, &lo);
        let b = geometry.position(geometry.outer, &hi);
        bounds[geometry.outer] = [a.clone().min(b.clone()), a.max(b)];
        let a = geometry.position(geometry.inner(), &lo);
        let b = geometry.position(geometry.inner(), &self.model.size.mul_pow2(-1));
        bounds[geometry.inner()] = [a.clone().min(b.clone()), a.max(b)];
        let volume = (0..3).fold(q(1), |v, i| &v * &(&bounds[i][1] - &bounds[i][0]));
        let density = Contribution::new(
            &Cell {
                bounds,
                volume,
                priority: q(0),
                serial: 0,
            },
            &self.model.anchor,
            false,
        );
        let contribution = if let Some(atlas) = Atlas::new(&self.model, geometry, &lo, &hi, dec)? {
            let interval = cover(&lo, &hi);
            let mid = (&lo + &hi).mul_pow2(-1);
            let radius = (&hi - &lo).mul_pow2(-1);
            let center = atlas.evaluate(&self.model, geometry, &Jet::variable(ball(&mid)));
            let derivatives =
                atlas.evaluate(&self.model, geometry, &Jet::variable(interval.clone()));
            let values: [Ball; 10] = if let (Some(center), Some(derivatives)) =
                (center, derivatives)
            {
                std::array::from_fn(|i| {
                    let mut integral = integer(0);
                    let mut power = radius.clone();
                    for k in 0..TERMS - 1 {
                        if k % 2 == 0 {
                            let weight = power.mul_pow2(1).checked_div(&q((k + 1) as i64)).unwrap();
                            integral = &integral + &(&center[i].0[k] * &ball(&weight));
                        }
                        power = &power * &radius;
                    }
                    let error = (&magnitude(&derivatives[i].0[TERMS - 1]) * &power)
                        .mul_pow2(1)
                        .checked_div(&q(TERMS as i64))
                        .unwrap();
                    &integral + &cover(&(-&error), &error)
                })
            } else {
                let range = atlas
                    .evaluate(&self.model, geometry, &Jet::constant(interval))
                    .expect("certified column interval range");
                std::array::from_fn(|i| &range[i].0[0] * &ball(&(&hi - &lo)))
            };
            // A Taylor remainder can be much wider than the material box
            // near a singular endpoint. Intersect independent certificates
            // before selecting refinement priority, so it cannot starve
            // the rest of the domain with an artificial unbounded error.
            let lo: [Rational; 10] =
                std::array::from_fn(|i| values[i].lo().max(density.lo[i].clone()));
            let hi: [Rational; 10] =
                std::array::from_fn(|i| values[i].hi().min(density.hi[i].clone()));
            if (0..10).any(|i| lo[i] > hi[i]) {
                return Err(failure(
                    "field.mass-accuracy",
                    "independent column moment certificates do not overlap",
                ));
            }
            let estimate = std::array::from_fn(|i| (&lo[i] + &hi[i]).mul_pow2(-1));
            Contribution { lo, hi, estimate }
        } else {
            density
        };
        let priority = contribution
            .lo
            .iter()
            .zip(&contribution.hi)
            .map(|(lo, hi)| hi - lo)
            .max()
            .unwrap();
        self.totals.update(&contribution, false);
        if !priority.is_zero() {
            self.pending.push(Pending {
                sector,
                lo,
                hi,
                contribution,
                priority,
                serial: index,
            });
        }
        Ok(())
    }
    fn properties(&self) -> Option<Properties> {
        self.totals.properties(&self.model.anchor, self.cells)
    }
    fn refusal(&self) -> Refusal {
        let detail=self.properties().map_or_else(||format!("{} analytic intervals; positive field volume not yet certified",self.cells),|p|format!("{} analytic intervals; volume width {} mm3, maximum inertia width {:?} mm5; requested Accuracy not reached",self.cells,p.volume.width(),p.inertia.as_ref().map(|v|v.iter().flatten().map(catallax_measure::Enclosure::width).max().unwrap().to_string())));
        failure("field.mass-budget", detail)
    }
}

pub(super) fn properties(
    plan: &Plan,
    bounds: &Bounds,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Option<Result<Properties, Refusal>> {
    let model = Model::new(plan, bounds)?;
    Some(integrate(model, accuracy, options, dec))
}
fn integrate(
    mut model: Model<'_>,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Result<Properties, Refusal> {
    let xs = model.nodes(0, options.max_cells).ok_or_else(|| {
        failure(
            "field.mass-budget",
            "lattice column index exceeds the cell limit",
        )
    })?;
    let ys = model.nodes(1, options.max_cells).ok_or_else(|| {
        failure(
            "field.mass-budget",
            "lattice column index exceeds the cell limit",
        )
    })?;
    if xs.len().saturating_mul(ys.len()).saturating_mul(8) > options.max_cells {
        return Err(failure(
            "field.mass-budget",
            "initial lattice sectors exceed the cell limit",
        ));
    }
    model.z_nodes = model.nodes(2, options.max_cells).ok_or_else(|| {
        failure(
            "field.mass-budget",
            "lattice column index exceeds the cell limit",
        )
    })?;
    let mut integration = Integration {
        model,
        sectors: Vec::new(),
        totals: Totals::new(),
        pending: BinaryHeap::new(),
        cells: 0,
        max_cells: options.max_cells,
    };
    let half = integration.model.size.mul_pow2(-1);
    for x in xs {
        for y in &ys {
            for outer in 0..2 {
                for sx in [-1, 1] {
                    for sy in [-1, 1] {
                        let sector = Sector {
                            node: [x.clone(), y.clone()],
                            outer,
                            signs: [sx, sy],
                        };
                        let mut cuts =
                            BTreeSet::from([q(0), half.clone(), integration.model.radius.clone()]);
                        for definition in &integration.model.definitions {
                            match definition {
                                Definition::Box(bounds) => {
                                    for endpoint in &bounds[outer] {
                                        cuts.insert(sector.coordinate(outer, endpoint));
                                    }
                                }
                                Definition::Cylinder { center, radius, .. } => {
                                    let shift = sector.coordinate(outer, &center[outer]);
                                    cuts.insert(&shift - radius);
                                    cuts.insert(&shift + radius);
                                }
                                _ => {}
                            }
                        }
                        let cuts: Vec<Rational> = cuts
                            .into_iter()
                            .filter(|v| *v >= q(0) && *v <= half)
                            .collect();
                        let index = integration.sectors.len();
                        integration.sectors.push(sector);
                        for pair in cuts.windows(2) {
                            integration.admit(index, pair[0].clone(), pair[1].clone(), dec)?;
                        }
                    }
                }
            }
        }
    }
    loop {
        if let Some(result) = integration.properties()
            && reached(&result.volume, &accuracy.volume)
            && result
                .centroid
                .as_ref()
                .is_none_or(|c| c.iter().all(|v| reached(v, &accuracy.centroid)))
            && result
                .inertia
                .as_ref()
                .is_none_or(|m| m.iter().flatten().all(|v| reached(v, &accuracy.inertia)))
        {
            return Ok(result);
        }
        for _ in 0..128 {
            if integration.cells.saturating_add(2) > options.max_cells {
                return Err(integration.refusal());
            }
            let Some(parent) = integration.pending.pop() else {
                return Err(integration.refusal());
            };
            integration.totals.update(&parent.contribution, true);
            let mid = (&parent.lo + &parent.hi).mul_pow2(-1);
            integration.admit(parent.sector, parent.lo, mid.clone(), dec)?;
            integration.admit(parent.sector, mid, parent.hi, dec)?;
        }
    }
}
