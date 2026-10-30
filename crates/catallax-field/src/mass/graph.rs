//! Moments of an axis box clipped by one smooth field leaf, through monotone graphs.
//!
//! The region is an intersection of axis-aligned blocks or prisms. The leaf is a
//! gyroid, or an axis cubic lattice whose radius is an axis grade. Blocks of the
//! region are split until the leaf is classified, or until one smooth expression
//! `F` defines the leaf over the block and one partial derivative `F_d` keeps its
//! sign there. The material of such a block lies on one side of a graph
//! `d = ζ(e)` over the other two coordinates `e`, so each column's inner integral
//! is an exact polynomial of `ζ`.
//!
//! Over a cell whose root stays strictly inside the column, the outer integral is
//! the midpoint rule. Taylor's theorem bounds its error by
//! `A (M₁₁ w₁²/24 + M₁₂ w₁ w₂/16 + M₂₂ w₂²/24)` for cell area `A`, widths `w` and
//! bounds `M` on the integrand's second derivatives over the cell. Those come from
//! implicit derivatives of `ζ`, with `F`'s derivatives enclosed over the cell and
//! the root's range. A cell where the root leaves its column keeps exact moments of
//! its certain material plus a material-density enclosure of the rest, as does a
//! block that never becomes a graph. Nothing is sampled away.
//!
//! Sums are exact fixed-point integers, so removing a refined parent's
//! contribution never widens the totals.

use super::{
    Accuracy, Bounds, Cell, Contribution, Options, Properties, Totals, axis_block, failure, q,
    reached,
};
use crate::{
    Instruction, Plan, geometry,
    grade::Piece as GradePiece,
    gyroid::Gyroid,
    lattice::{BoxState, Lattice, LinePiece},
};
use catallax_numbers::{Decider, Integer, Interval, Rational, Refusal, Sign};
use std::{cmp::Ordering, collections::BinaryHeap};

/// Fixed-point fraction bits of the moment sums.
const SCALE: i64 = 48;
/// Exponents of the ten coordinate moments, in the order of [`Contribution`].
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

type Moments = [[i128; 2]; 10];

fn point(v: f64) -> Interval {
    Interval::point(v).unwrap()
}
fn enclose(bounds: &[Rational; 2]) -> Interval {
    Interval::new(
        Interval::from_rational(&bounds[0]).lo(),
        Interval::from_rational(&bounds[1]).hi(),
    )
    .unwrap()
}
fn magnitude(v: Interval) -> f64 {
    v.lo().abs().max(v.hi().abs())
}
/// `xⁿ` for `n <= 3`, monotone where it can be.
fn power(x: Interval, n: usize) -> Interval {
    match n {
        0 => point(1.0),
        1 => x,
        2 => x.square(),
        _ => {
            let lo = Interval::point(x.lo()).unwrap();
            let hi = Interval::point(x.hi()).unwrap();
            let a = lo.square() * lo;
            let b = hi.square() * hi;
            Interval::new(a.lo(), b.hi()).unwrap()
        }
    }
}
/// Outward fixed-point endpoints. Scaling by a power of two is exact, so only the
/// final floor and ceiling round, and both round outward.
fn fixed(value: Interval) -> Option<[i128; 2]> {
    let scale = (1u64 << SCALE) as f64;
    let lo = (value.lo() * scale).floor();
    let hi = (value.hi() * scale).ceil();
    // Far below i128's range even after summing millions of cells.
    let limit = 1e30;
    (lo.is_finite() && hi.is_finite() && lo.abs() < limit && hi.abs() < limit)
        .then_some([lo as i128, hi as i128])
}
fn rational(value: i128) -> Rational {
    let magnitude = value.unsigned_abs();
    let high = Rational::from_integer(Integer::from_u64((magnitude >> 64) as u64)).mul_pow2(64);
    let low = Rational::from_integer(Integer::from_u64(magnitude as u64));
    let result = (&high + &low).mul_pow2(-SCALE);
    if value < 0 { -&result } else { result }
}
fn rational_moments(contribution: &Contribution) -> Option<Moments> {
    let mut moments = [[0; 2]; 10];
    for (i, entry) in moments.iter_mut().enumerate() {
        *entry = fixed(Interval::new(
            Interval::from_rational(&contribution.lo[i]).lo(),
            Interval::from_rational(&contribution.hi[i]).hi(),
        )?)?;
    }
    Some(moments)
}
fn volume(bounds: &Bounds) -> Rational {
    (0..3).fold(q(1), |v, i| &v * &(&bounds[i][1] - &bounds[i][0]))
}
fn box_moments(bounds: &Bounds, anchor: &[Rational; 3], inside: bool) -> Contribution {
    Contribution::new(
        &Cell {
            bounds: bounds.clone(),
            volume: volume(bounds),
            priority: q(0),
            serial: 0,
        },
        anchor,
        inside,
    )
}

/// The smooth leaf.
#[derive(Clone, Copy)]
enum Smooth<'a> {
    Gyroid(&'a Gyroid),
    Lattice(&'a Lattice),
}
/// One smooth expression of the leaf over a block.
#[derive(Clone, Debug)]
enum Piece {
    Gyroid(GradePiece),
    Line(LinePiece),
}
impl Smooth<'_> {
    fn state(self, world: &[Interval; 3], dec: &mut Decider) -> Result<BoxState, Refusal> {
        match self {
            Self::Gyroid(gyroid) => gyroid.box_state(world, dec),
            Self::Lattice(lattice) => lattice.box_state(world, dec),
        }
    }
    fn piece(self, bounds: &Bounds) -> Option<Piece> {
        match self {
            Self::Gyroid(gyroid) => gyroid.level_piece(bounds).map(Piece::Gyroid),
            Self::Lattice(lattice) => lattice.line_piece(bounds).map(Piece::Line),
        }
    }
    fn value(self, piece: &Piece, world: &[Interval; 3]) -> Option<Interval> {
        match (self, piece) {
            (Self::Gyroid(gyroid), Piece::Gyroid(_)) => gyroid.residual(world),
            (Self::Lattice(lattice), Piece::Line(line)) => lattice.line_value(line, world),
            _ => None,
        }
    }
    fn derivatives(
        self,
        piece: &Piece,
        world: &[Interval; 3],
    ) -> Option<([Interval; 3], [[Interval; 3]; 3])> {
        match (self, piece) {
            (Self::Gyroid(gyroid), Piece::Gyroid(level)) => gyroid.derivatives(world, *level),
            (Self::Lattice(lattice), Piece::Line(line)) => lattice.line_derivatives(line, world),
            _ => None,
        }
    }
    fn planes(self) -> Option<(usize, [Rational; 2])> {
        match self {
            Self::Gyroid(gyroid) => gyroid.level_planes(),
            Self::Lattice(lattice) => lattice.grade_planes(),
        }
    }
}

#[derive(Clone)]
enum Kind {
    Block,
    Column {
        piece: Piece,
        axis: usize,
        rising: bool,
    },
}
struct Item {
    bounds: Bounds,
    kind: Kind,
    moments: Moments,
    priority: i128,
    split: usize,
    serial: u64,
}
impl PartialEq for Item {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
}
impl Eq for Item {}
impl PartialOrd for Item {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Item {
    fn cmp(&self, other: &Self) -> Ordering {
        self.priority
            .cmp(&other.priority)
            .then_with(|| other.serial.cmp(&self.serial))
    }
}

/// What one column cell holds.
enum Column {
    Empty,
    Full,
    /// Enclosed moments and the cross axis whose refinement shrinks them most.
    Graph(Box<[Interval; 10]>, usize),
    /// The column part certain to be material, and the uncertain slab beside it.
    Kink(Box<[Bounds; 2]>),
}

struct Model<'a> {
    smooth: Smooth<'a>,
    region: Bounds,
    anchor: [Rational; 3],
    centre: [Interval; 3],
}

/// The axis bounds of a convex prism whose faces are all axis planes.
fn prism(feature: &crate::features::Feature, dec: &mut Decider) -> Result<Option<Bounds>, Refusal> {
    let Some(planes) = feature.convex_planes(dec)? else {
        return Ok(None);
    };
    let mut lower: [Option<Rational>; 3] = [None, None, None];
    let mut upper: [Option<Rational>; 3] = [None, None, None];
    for plane in planes {
        let (Some(normal), Some(origin)) = (plane.normal.as_rational(), plane.origin.as_rational())
        else {
            return Ok(None);
        };
        let normal = [normal.x, normal.y, normal.z];
        let origin = [origin.x, origin.y, origin.z];
        let mut axes = (0..3).filter(|&i| !normal[i].is_zero());
        let (Some(axis), None) = (axes.next(), axes.next()) else {
            return Ok(None);
        };
        // Outward normals: a positive component bounds the coordinate from above.
        let side = if normal[axis].sign() == Sign::Positive {
            &mut upper[axis]
        } else {
            &mut lower[axis]
        };
        let tighter = match side.take() {
            Some(old) if normal[axis].sign() == Sign::Positive => old.min(origin[axis].clone()),
            Some(old) => old.max(origin[axis].clone()),
            None => origin[axis].clone(),
        };
        *side = Some(tighter);
    }
    let mut bounds: Bounds = std::array::from_fn(|_| [q(0), q(0)]);
    for i in 0..3 {
        let (Some(lo), Some(hi)) = (lower[i].take(), upper[i].take()) else {
            return Ok(None);
        };
        bounds[i] = [lo, hi];
    }
    Ok(Some(bounds))
}

impl<'a> Model<'a> {
    /// The region and the smooth leaf of a supported plan, or `None`.
    fn new(plan: &'a Plan, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        // The result must be an intersection of leaves.
        let mut leaves = Vec::new();
        let mut pending = vec![plan.root];
        while let Some(node) = pending.pop() {
            match plan.instructions[node] {
                Instruction::Leaf(leaf) => leaves.push(leaf),
                Instruction::Boolean {
                    kind: catallax_recipe::BooleanKind::Intersect,
                    a,
                    b,
                } => pending.extend([a, b]),
                _ => return Ok(None),
            }
        }
        leaves.sort_unstable();
        leaves.dedup();
        let mut smooth = None;
        let mut region: Option<Bounds> = None;
        for leaf in leaves {
            let bounds = match &plan.leaves[leaf] {
                geometry::Leaf::Gyroid(gyroid) => {
                    if smooth.replace(Smooth::Gyroid(gyroid)).is_some() {
                        return Ok(None);
                    }
                    continue;
                }
                geometry::Leaf::Lattice(lattice) if lattice.axis_graded() => {
                    if smooth.replace(Smooth::Lattice(lattice)).is_some() {
                        return Ok(None);
                    }
                    continue;
                }
                geometry::Leaf::Primitive(primitive) => match axis_block(primitive) {
                    Some(bounds) => bounds,
                    None => return Ok(None),
                },
                geometry::Leaf::Feature(feature) => match prism(feature, dec)? {
                    Some(bounds) => bounds,
                    None => return Ok(None),
                },
                _ => return Ok(None),
            };
            region = Some(match region {
                None => bounds,
                Some(old) => std::array::from_fn(|i| {
                    [
                        old[i][0].clone().max(bounds[i][0].clone()),
                        old[i][1].clone().min(bounds[i][1].clone()),
                    ]
                }),
            });
        }
        let (Some(smooth), Some(region)) = (smooth, region) else {
            return Ok(None);
        };
        let anchor: [Rational; 3] =
            std::array::from_fn(|i| (&region[i][0] + &region[i][1]).mul_pow2(-1));
        Ok(Some(Self {
            smooth,
            centre: anchor.each_ref().map(Interval::from_rational),
            region,
            anchor,
        }))
    }

    fn world(bounds: &Bounds) -> [Interval; 3] {
        bounds.each_ref().map(enclose)
    }

    /// A column axis along which `F` keeps its sign of change over the block.
    fn monotone(&self, piece: &Piece, world: &[Interval; 3]) -> Option<(usize, bool)> {
        let (gradient, _) = self.smooth.derivatives(piece, world)?;
        let (axis, slope) = gradient
            .iter()
            .enumerate()
            .filter(|(_, g)| g.lo() > 0.0 || g.hi() < 0.0)
            .max_by(|(_, a), (_, b)| {
                a.lo()
                    .abs()
                    .min(a.hi().abs())
                    .total_cmp(&b.lo().abs().min(b.hi().abs()))
            })?;
        Some((axis, slope.lo() > 0.0))
    }

    fn column(
        &self,
        bounds: &Bounds,
        piece: &Piece,
        d: usize,
        rising: bool,
        dec: &mut Decider,
    ) -> Result<Option<Column>, Refusal> {
        let world = Self::world(bounds);
        let at = |t: Interval| {
            let mut face = world;
            face[d] = t;
            face
        };
        let [d0, d1] = bounds[d].each_ref().map(Interval::from_rational);
        let (Some(low), Some(high)) = (
            self.smooth.value(piece, &at(d0)),
            self.smooth.value(piece, &at(d1)),
        ) else {
            return Ok(None);
        };
        // Material is where F <= 0: below the root when F rises along d.
        if (rising && low.lo() > 0.0) || (!rising && high.lo() > 0.0) {
            return Ok(Some(Column::Empty));
        }
        if (rising && high.hi() < 0.0) || (!rising && low.hi() < 0.0) {
            return Ok(Some(Column::Full));
        }
        let smooth = if rising {
            low.hi() < 0.0 && high.lo() > 0.0
        } else {
            low.lo() > 0.0 && high.hi() < 0.0
        };
        if !smooth {
            return self.kink(bounds, piece, d, rising, dec);
        }
        let [e1, e2] = if d == 1 {
            [0, 2]
        } else {
            [(d + 1) % 3, (d + 2) % 3]
        };
        let centre_of = |i: usize| (&bounds[i][0] + &bounds[i][1]).mul_pow2(-1);
        let mut centre = world;
        centre[e1] = Interval::from_rational(&centre_of(e1));
        centre[e2] = Interval::from_rational(&centre_of(e2));
        let Some(root) = self.root(piece, centre, d, [d0.lo(), d1.hi()], dec)? else {
            return Ok(None);
        };
        let Some((gradient, _)) = self.smooth.derivatives(piece, &world) else {
            return Ok(None);
        };
        let widths = [e1, e2].map(|i| Interval::from_rational(&(&bounds[i][1] - &bounds[i][0])));
        let fd = gradient[d];
        if fd.lo() <= 0.0 && fd.hi() >= 0.0 {
            return Ok(None);
        }
        let mut spread = point(0.0);
        for (axis, width) in [e1, e2].into_iter().zip(widths) {
            spread = spread + point(magnitude(gradient[axis] / fd)) * width * point(0.5);
        }
        let spread = spread.hi();
        let zeta = Interval::new(
            (root.lo() - spread).max(d0.lo()),
            (root.hi() + spread).min(d1.hi()),
        )
        .unwrap();
        let mut tight = world;
        tight[d] = zeta;
        let Some((g, h)) = self.smooth.derivatives(piece, &tight) else {
            return Ok(None);
        };
        let fd = g[d];
        if fd.lo() <= 0.0 && fd.hi() >= 0.0 {
            return Ok(None);
        }
        let z = [e1, e2].map(|i| -(g[i] / fd));
        let second = |i: usize, j: usize| {
            let (a, b) = ([e1, e2][i], [e1, e2][j]);
            -((h[a][b] + h[a][d] * z[j] + h[b][d] * z[i] + h[d][d] * z[i] * z[j]) / fd)
        };
        let zz = [[second(0, 0), second(0, 1)], [second(1, 0), second(1, 1)]];
        let relative = |value: Interval, axis: usize| value - self.centre[axis];
        let u = [e1, e2].map(|i| relative(world[i], i));
        let uc = [e1, e2].map(|i| relative(centre[i], i));
        let zr = relative(zeta, d);
        let zc = relative(root, d);
        let base = relative(if rising { d0 } else { d1 }, d);
        let area = widths[0] * widths[1];
        let mut values = [point(0.0); 10];
        let mut errors = [0.0_f64; 2];
        for (m, exponents) in EXPONENTS.iter().enumerate() {
            let [a, b] = [exponents[e1], exponents[e2]];
            let gamma = exponents[d];
            let sign = point(if rising { 1.0 } else { -1.0 });
            let h0 = |zeta: Interval| {
                sign * (power(zeta, gamma + 1) - power(base, gamma + 1)) / point((gamma + 1) as f64)
            };
            let h1 = sign * power(zr, gamma);
            let h2 = if gamma == 0 {
                point(0.0)
            } else {
                sign * point(gamma as f64) * power(zr, gamma - 1)
            };
            let p = |x: [Interval; 2], [i, j]: [usize; 2]| {
                let factor = |n: usize, k: usize, x: Interval| {
                    if k > n {
                        point(0.0)
                    } else {
                        let falling = (0..k).fold(1.0, |f, s| f * (n - s) as f64);
                        point(falling) * power(x, n - k)
                    }
                };
                factor(a, i, x[0]) * factor(b, j, x[1])
            };
            let hz = h0(zr);
            let value = p(uc, [0, 0]) * h0(zc);
            let d11 = p(u, [2, 0]) * hz
                + point(2.0) * p(u, [1, 0]) * h1 * z[0]
                + p(u, [0, 0]) * (h2 * z[0].square() + h1 * zz[0][0]);
            let d12 = p(u, [1, 1]) * hz
                + p(u, [1, 0]) * h1 * z[1]
                + p(u, [0, 1]) * h1 * z[0]
                + p(u, [0, 0]) * (h2 * z[0] * z[1] + h1 * zz[0][1]);
            let d22 = p(u, [0, 2]) * hz
                + point(2.0) * p(u, [0, 1]) * h1 * z[1]
                + p(u, [0, 0]) * (h2 * z[1].square() + h1 * zz[1][1]);
            let [w1, w2] = widths;
            let e11 = point(magnitude(d11)) * w1.square() / point(24.0);
            let e12 = point(magnitude(d12)) * w1 * w2 / point(16.0);
            let e22 = point(magnitude(d22)) * w2.square() / point(24.0);
            let error = (area * (e11 + e12 + e22)).hi();
            if !error.is_finite() {
                return Ok(None);
            }
            errors[0] = errors[0].max((area * e11).hi());
            errors[1] = errors[1].max((area * e22).hi());
            values[m] = area * value + Interval::new(-error, error).unwrap();
        }
        let split = if errors[0] >= errors[1] { e1 } else { e2 };
        Ok(Some(Column::Graph(Box::new(values), split)))
    }

    /// The root of `F` along `d` at a cell centre, by interval Newton steps.
    fn root(
        &self,
        piece: &Piece,
        centre: [Interval; 3],
        d: usize,
        [lo, hi]: [f64; 2],
        dec: &mut Decider,
    ) -> Result<Option<Interval>, Refusal> {
        let mut t = Interval::new(lo, hi).unwrap();
        for _ in 0..64 {
            dec.check_deadline()?;
            let mid = t.lo() * 0.5 + t.hi() * 0.5;
            let mut at = centre;
            at[d] = point(mid);
            let Some(value) = self.smooth.value(piece, &at) else {
                return Ok(None);
            };
            at[d] = t;
            let Some((gradient, _)) = self.smooth.derivatives(piece, &at) else {
                return Ok(None);
            };
            let slope = gradient[d];
            if slope.lo() <= 0.0 && slope.hi() >= 0.0 {
                return Ok(None);
            }
            let step = point(mid) - value / slope;
            let lo = t.lo().max(step.lo());
            let hi = t.hi().min(step.hi());
            if lo > hi {
                return Ok(None);
            }
            let next = Interval::new(lo, hi).unwrap();
            let settled = next.hi() - next.lo() >= 0.5 * (t.hi() - t.lo());
            t = next;
            if settled {
                break;
            }
        }
        Ok(Some(t))
    }

    /// Brackets the root over a cell whose root leaves the column somewhere.
    fn kink(
        &self,
        bounds: &Bounds,
        piece: &Piece,
        d: usize,
        rising: bool,
        dec: &mut Decider,
    ) -> Result<Option<Column>, Refusal> {
        let world = Self::world(bounds);
        // Whether every point of the cell at `t` is strictly on the material side
        // (`below`) or strictly on the void side of the root.
        let side = |t: f64, material: bool, dec: &mut Decider| -> Result<bool, Refusal> {
            dec.check_deadline()?;
            let mut face = world;
            face[d] = point(t);
            Ok(self
                .smooth
                .value(piece, &face)
                .is_some_and(|v| if material { v.hi() < 0.0 } else { v.lo() > 0.0 }))
        };
        let [d0, d1] = [bounds[d][0].to_f64_nearest(), bounds[d][1].to_f64_nearest()];
        // The material end of the column is d0 when F rises, d1 when it falls.
        let (start, end) = if rising { (d0, d1) } else { (d1, d0) };
        let mut certain = start;
        let mut beyond = end;
        let mut probe = [start, end];
        for _ in 0..12 {
            let mid = probe[0] * 0.5 + probe[1] * 0.5;
            if side(mid, true, dec)? {
                probe[0] = mid;
                certain = mid;
            } else {
                probe[1] = mid;
            }
        }
        let mut probe = [certain, end];
        for _ in 0..12 {
            let mid = probe[0] * 0.5 + probe[1] * 0.5;
            if side(mid, false, dec)? {
                probe[1] = mid;
                beyond = mid;
            } else {
                probe[0] = mid;
            }
        }
        let (Some(certain), Some(beyond)) = (
            Rational::from_f64_exact(certain),
            Rational::from_f64_exact(beyond),
        ) else {
            return Ok(None);
        };
        let clamp = |v: Rational| v.max(bounds[d][0].clone()).min(bounds[d][1].clone());
        let (certain, beyond) = (clamp(certain), clamp(beyond));
        let mut sure = bounds.clone();
        let mut open = bounds.clone();
        if rising {
            sure[d][1] = certain.clone();
            open[d] = [certain.clone(), beyond.clone().max(certain)];
        } else {
            sure[d][0] = certain.clone();
            open[d] = [beyond.clone().min(certain.clone()), certain];
        }
        Ok(Some(Column::Kink(Box::new([sure, open]))))
    }
}

struct Integration<'a> {
    model: Model<'a>,
    lo: [i128; 10],
    hi: [i128; 10],
    pending: BinaryHeap<Item>,
    cells: usize,
    serial: u64,
    max_cells: usize,
}

fn add_moments(total: &mut Integration<'_>, moments: &Moments, remove: bool) {
    for ((lo, hi), [low, high]) in total.lo.iter_mut().zip(&mut total.hi).zip(moments) {
        if remove {
            *lo -= low;
            *hi -= high;
        } else {
            *lo += low;
            *hi += high;
        }
    }
}

impl Integration<'_> {
    fn range_failure() -> Refusal {
        failure(
            "field.mass-budget",
            "a field moment exceeds the fixed-point summation range",
        )
    }

    fn admit(
        &mut self,
        bounds: Bounds,
        column: Option<(Piece, usize, bool)>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if self.cells == self.max_cells {
            return Err(failure(
                "field.mass-budget",
                "field moment cell limit reached before the requested Accuracy",
            ));
        }
        self.cells += 1;
        dec.check_deadline()?;
        let anchor = &self.model.anchor;
        let world = Model::world(&bounds);
        let state = self.model.smooth.state(&world, dec)?;
        if state == BoxState::Outside {
            return Ok(());
        }
        if state == BoxState::Inside {
            let moments = rational_moments(&box_moments(&bounds, anchor, true))
                .ok_or_else(Self::range_failure)?;
            add_moments(self, &moments, false);
            return Ok(());
        }
        let column = match column {
            Some(column) => Some(column),
            None => self.model.smooth.piece(&bounds).and_then(|piece| {
                let (axis, rising) = self.model.monotone(&piece, &world)?;
                Some((piece, axis, rising))
            }),
        };
        let density = box_moments(&bounds, anchor, false);
        let longest = (0..3)
            .max_by(|&a, &b| (&bounds[a][1] - &bounds[a][0]).cmp(&(&bounds[b][1] - &bounds[b][0])))
            .unwrap();
        let (moments, kind, split) = match column {
            None => (
                rational_moments(&density).ok_or_else(Self::range_failure)?,
                Kind::Block,
                longest,
            ),
            Some((piece, axis, rising)) => {
                let cross = |bounds: &Bounds| {
                    [(axis + 1) % 3, (axis + 2) % 3]
                        .into_iter()
                        .max_by(|&a, &b| {
                            (&bounds[a][1] - &bounds[a][0]).cmp(&(&bounds[b][1] - &bounds[b][0]))
                        })
                        .unwrap()
                };
                let kind = Kind::Column {
                    piece: piece.clone(),
                    axis,
                    rising,
                };
                match self.model.column(&bounds, &piece, axis, rising, dec)? {
                    Some(Column::Empty) => return Ok(()),
                    Some(Column::Full) => {
                        let moments = rational_moments(&box_moments(&bounds, anchor, true))
                            .ok_or_else(Self::range_failure)?;
                        add_moments(self, &moments, false);
                        return Ok(());
                    }
                    Some(Column::Graph(values, split)) => {
                        let density = rational_moments(&density).ok_or_else(Self::range_failure)?;
                        let mut moments = [[0; 2]; 10];
                        for i in 0..10 {
                            let graph = fixed(values[i]).ok_or_else(Self::range_failure)?;
                            moments[i] = [graph[0].max(density[i][0]), graph[1].min(density[i][1])];
                            if moments[i][0] > moments[i][1] {
                                return Err(failure(
                                    "field.mass-accuracy",
                                    "independent graph and density moment certificates do not overlap",
                                ));
                            }
                        }
                        (moments, kind, split)
                    }
                    Some(Column::Kink(parts)) => {
                        let [certain, uncertain] = *parts;
                        let sure = rational_moments(&box_moments(&certain, anchor, true))
                            .ok_or_else(Self::range_failure)?;
                        let open = rational_moments(&box_moments(&uncertain, anchor, false))
                            .ok_or_else(Self::range_failure)?;
                        let mut moments = [[0; 2]; 10];
                        for i in 0..10 {
                            moments[i] = [sure[i][0] + open[i][0], sure[i][1] + open[i][1]];
                        }
                        (moments, kind, cross(&bounds))
                    }
                    None => (
                        rational_moments(&density).ok_or_else(Self::range_failure)?,
                        kind,
                        cross(&bounds),
                    ),
                }
            }
        };
        add_moments(self, &moments, false);
        let priority = moments.iter().map(|m| m[1] - m[0]).max().unwrap();
        if priority > 0 {
            self.pending.push(Item {
                bounds,
                kind,
                moments,
                priority,
                split,
                serial: self.serial,
            });
            self.serial += 1;
        }
        Ok(())
    }

    fn totals(&self) -> Totals {
        Totals {
            lo: self.lo.map(rational),
            hi: self.hi.map(rational),
            estimate: std::array::from_fn(|i| rational(self.lo[i] / 2 + self.hi[i] / 2)),
        }
    }

    fn properties(&self) -> Option<Properties> {
        self.totals().properties(&self.model.anchor, self.cells)
    }

    fn budget_failure(&self) -> Refusal {
        let detail = if let Some(properties) = self.properties() {
            let centroid = properties.centroid.as_ref().map(|c| {
                c.iter()
                    .map(catallax_measure::Enclosure::width)
                    .max()
                    .unwrap()
            });
            let inertia = properties.inertia.as_ref().map(|m| {
                m.iter()
                    .flatten()
                    .map(catallax_measure::Enclosure::width)
                    .max()
                    .unwrap()
            });
            format!(
                "{} graph cells; reached volume width {} mm3, maximum centroid width {:?} mm, maximum inertia width {:?} mm5; requested Accuracy not reached",
                self.cells,
                properties.volume.width(),
                centroid.map(|q| q.to_string()),
                inertia.map(|q| q.to_string())
            )
        } else {
            format!(
                "{} graph cells; positive field volume is not yet certified",
                self.cells
            )
        };
        failure("field.mass-budget", detail)
    }
}

/// Moments of a supported plan, or `None` for every other plan.
pub(super) fn properties(
    plan: &Plan,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Option<Result<Properties, Refusal>> {
    let model = match Model::new(plan, dec) {
        Ok(Some(model)) => model,
        Ok(None) => return None,
        Err(error) => return Some(Err(error)),
    };
    Some(integrate(model, accuracy, options, dec))
}

fn integrate(
    model: Model<'_>,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Result<Properties, Refusal> {
    if model.region.iter().any(|b| b[1] <= b[0]) {
        return Ok(Properties {
            volume: catallax_measure::Enclosure { lo: q(0), hi: q(0) },
            centroid: None,
            inertia: None,
            estimate: None,
            cells: 0,
        });
    }
    // A graded leaf is affine on each side of its anchor planes; start the
    // partition there so no cell straddles a clamp.
    let mut blocks = vec![model.region.clone()];
    if let Some((axis, planes)) = model.smooth.planes() {
        for plane in planes {
            blocks = blocks
                .into_iter()
                .flat_map(|block| {
                    if plane > block[axis][0] && plane < block[axis][1] {
                        let mut lower = block.clone();
                        let mut upper = block;
                        lower[axis][1] = plane.clone();
                        upper[axis][0] = plane.clone();
                        vec![lower, upper]
                    } else {
                        vec![block]
                    }
                })
                .collect();
        }
    }
    let mut integration = Integration {
        model,
        lo: [0; 10],
        hi: [0; 10],
        pending: BinaryHeap::new(),
        cells: 0,
        serial: 0,
        max_cells: options.max_cells,
    };
    for block in blocks {
        integration.admit(block, None, dec)?;
    }
    loop {
        if let Some(result) = integration.properties() {
            let centroid_ok = result
                .centroid
                .as_ref()
                .is_none_or(|c| c.iter().all(|v| reached(v, &accuracy.centroid)));
            let inertia_ok = result
                .inertia
                .as_ref()
                .is_none_or(|m| m.iter().flatten().all(|v| reached(v, &accuracy.inertia)));
            if reached(&result.volume, &accuracy.volume) && centroid_ok && inertia_ok {
                return Ok(result);
            }
        }
        if integration.pending.is_empty() {
            return Err(failure(
                "field.mass-accuracy",
                "graph moment certificates did not meet the requested Accuracy",
            ));
        }
        for _ in 0..128 {
            if !integration.pending.is_empty()
                && integration.cells.saturating_add(2) > integration.max_cells
            {
                return Err(integration.budget_failure());
            }
            let Some(item) = integration.pending.pop() else {
                break;
            };
            add_moments(&mut integration, &item.moments, true);
            let axis = item.split;
            let middle = (&item.bounds[axis][0] + &item.bounds[axis][1]).mul_pow2(-1);
            let mut lower = item.bounds.clone();
            let mut upper = item.bounds;
            lower[axis][1] = middle.clone();
            upper[axis][0] = middle;
            let column = match item.kind {
                Kind::Block => None,
                Kind::Column {
                    piece,
                    axis,
                    rising,
                } => Some((piece, axis, rising)),
            };
            integration.admit(lower, column.clone(), dec)?;
            integration.admit(upper, column, dec)?;
        }
    }
}
