//! Material correspondence for constant axis-aligned cubic infill.
//!
//! Plantinga and Vegter subdivide boxes that can meet an implicit surface.
//! At hard Boolean joins use Raffalli's separating-direction flow (Theorem 2),
//! retaining all possible active gradients on each side of a crease. Deform
//! each squared coordinate into its continuous affine interpolant before
//! applying the same Boolean expression.
//! https://www.math.rug.nl/~gert/documents/2004/pv-sgp-2004.pdf
//! https://raffalli.eu/pdfs/isotopic.pdf

use super::{Options, bounds, q, refusal};
use crate::{Instruction, Plan, geometry::Leaf, metric::Kind};
use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_recipe::BooleanKind;

mod polygons;

type Point = [Rational; 3];
type Box3 = [[Rational; 2]; 3];
type Range = [Interval; 3];

enum Support {
    Block(Box3),
    Cylinder {
        centre: Point,
        radius: Rational,
        height: Rational,
    },
    Cubic {
        centre: Point,
        radius: Rational,
        period: Rational,
    },
}

impl Support {
    fn new(leaf: &Leaf) -> Option<Self> {
        match leaf {
            Leaf::Primitive(p) => {
                let inverse = p.inverse_placement();
                let matrix = inverse.rational_matrix()?;
                if (0..3).any(|i| (0..3).any(|j| *matrix[i][j] != q(i64::from(i == j)))) {
                    return None;
                }
                let centre = inverse.translation.rational_components()?.map(|x| -x);
                match p.kind {
                    Kind::Block => Some(Self::Block(std::array::from_fn(|i| {
                        [centre[i].clone(), &centre[i] + &p.parameters[i]]
                    }))),
                    Kind::Round if p.parameters[0] == p.parameters[1] => Some(Self::Cylinder {
                        centre,
                        radius: p.parameters[0].clone(),
                        height: p.parameters[2].clone(),
                    }),
                    _ => None,
                }
            }
            Leaf::Lattice(lattice) => {
                let (centre, period, radius) = lattice.axis_cubic_grid()?;
                (&radius * &q(2) < period).then_some(Self::Cubic {
                    centre,
                    radius,
                    period,
                })
            }
            _ => None,
        }
    }
}

struct Model<'a> {
    plan: &'a Plan,
    supports: Vec<Option<Support>>,
    pitch: Rational,
    domain: Box3,
    displacement: Rational,
    encoding_reserve: Rational,
}

pub(super) struct Certified {
    pub points: Vec<Point>,
    pub triangles: Vec<[usize; 3]>,
    pub displacement: Rational,
    pub cells: usize,
}

pub(super) fn candidate(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
) -> Option<Result<Certified, Refusal>> {
    if options.tolerance <= q(0) || options.max_cells == 0 || options.max_vertices < 4 {
        return None;
    }
    let mut supports = Vec::with_capacity(plan.leaves.len());
    for (idx, leaf) in plan.leaves.iter().enumerate() {
        let support = if plan.leaf_registers[idx].is_some() {
            Some(Support::new(leaf)?)
        } else {
            None
        };
        supports.push(support);
    }
    let pitch = supports
        .iter()
        .flatten()
        .find_map(|support| match support {
            Support::Cubic { period, radius, .. } if options.tolerance < *radius => Some(period),
            _ => None,
        })?
        .mul_pow2(-1);
    Some(build(plan, supports, pitch, options, dec))
}

fn build(
    plan: &Plan,
    supports: Vec<Option<Support>>,
    mut pitch: Rational,
    options: &Options,
    dec: &mut Decider,
) -> Result<Certified, Refusal> {
    // Common pitch must divide each half period so a periodic crease is a knot.
    // Admit one period per request here. The initial pitch spends at most half
    // the requested displacement on an isolated radial support (|grad| = 2r).
    let half = pitch.clone();
    let mut radius = None::<Rational>;
    for support in supports.iter().flatten() {
        let r = match support {
            Support::Cubic { period, radius, .. } if period.mul_pow2(-1) == half => radius,
            Support::Cubic { .. } => {
                return Err(refusal(
                    "field.print-detail-topology",
                    "cubic lattice print needs a common period",
                ));
            }
            Support::Cylinder { radius, .. } => radius,
            Support::Block(_) => continue,
        };
        radius = Some(radius.map_or_else(|| r.clone(), |old| old.min(r.clone())));
    }
    let diameter = &radius.expect("cubic support is present") * &q(2);
    let initial_error = &diameter * &options.tolerance;
    while &pitch * &pitch > initial_error {
        dec.charge_work(1)?;
        pitch = pitch.mul_pow2(-1);
    }

    let domain = bounds(plan, dec)?.map(|[a, b]| [&a - &q(1), &b + &q(1)]);
    // One binary32 spacing at the maximum absolute coordinate, per axis.
    // The extra spacing covers conversion through binary64.
    let mut encoding_reserve = q(0);
    for limits in &domain {
        let extent = limits[0].abs().max(limits[1].abs());
        let upper = (Interval::from_rational(&extent).hi() as f32).next_up();
        let spacing = f64::from(upper.next_up()) - f64::from(upper);
        encoding_reserve = &encoding_reserve
            + &Rational::from_f64_exact(spacing).ok_or_else(|| {
                refusal(
                    "field.print-coordinate",
                    "lattice placement exceeds binary32",
                )
            })?;
    }
    let displacement = &options.tolerance - &encoding_reserve;
    // Leave room to round each bound upwards when composing the STL header.
    let limit = displacement.to_f64_nearest().next_down();
    let displacement = Rational::from_f64_exact(limit)
        .filter(|v| *v > q(0))
        .ok_or_else(|| {
            refusal(
                "field.print-coordinate",
                "encoding consumes the requested accuracy",
            )
        })?;
    let model = Model {
        plan,
        supports,
        pitch,
        domain,
        displacement,
        encoding_reserve,
    };
    let mut budget = Budget { cells: 0, options };
    model.prove(&mut budget, dec)?;
    let (points, triangles) = polygons::mesh(&model, &mut budget, dec)?;
    Ok(Certified {
        points,
        triangles,
        displacement: model.displacement,
        cells: budget.cells,
    })
}

struct Budget<'a> {
    cells: usize,
    options: &'a Options,
}
impl Budget<'_> {
    fn visit(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        self.visit_many(1, dec)
    }

    fn visit_many(&mut self, count: usize, dec: &mut Decider) -> Result<(), Refusal> {
        if count > self.options.max_cells - self.cells {
            return Err(refusal(
                "field.print-budget",
                "lattice print cell limit reached",
            ));
        }
        self.cells += count;
        dec.charge_work(count as u64)
    }
}

fn ip(x: f64) -> Interval {
    Interval::point(x).expect("finite scalar interval")
}
fn negative(x: Interval) -> Interval {
    Interval::new(-x.hi(), -x.lo()).expect("interval negation")
}
fn enclosure([a, b]: &[Rational; 2]) -> Interval {
    Interval::new(
        Interval::from_rational(a).lo(),
        Interval::from_rational(b).hi(),
    )
    .expect("ordered bounds")
}

// Filip, Magedson and Markot's interpolation bound for Hessian [2, 0, 2].
// Replace this interface with catallax_shapes::tessellation::chord_bound_interval
// when the display mesher lands.
fn radial_interpolation_error(width: Interval) -> Interval {
    width.square() * ip(0.5)
}

#[derive(Clone, Copy)]
struct Branch {
    gradient: Range,
    motion: Interval,
}
impl Branch {
    fn reversed(self) -> Self {
        Self {
            gradient: self.gradient.map(negative),
            motion: self.motion,
        }
    }
}

struct Jet {
    value: Interval,
    active: Vec<Branch>,
}
impl Jet {
    fn new(value: Interval, gradient: Range, motion: Interval) -> Self {
        Self {
            value,
            active: vec![Branch { gradient, motion }],
        }
    }
    fn plane(value: Interval, axis: usize, sense: f64) -> Self {
        let mut gradient = [ip(0.0); 3];
        gradient[axis] = ip(sense);
        Self::new(value, gradient, ip(0.0))
    }
}

fn extreme_value(maximum: bool, a: Interval, b: Interval) -> Interval {
    let pick = if maximum { f64::max } else { f64::min };
    Interval::new(pick(a.lo(), b.lo()), pick(a.hi(), b.hi())).unwrap()
}

fn can_be_active(maximum: bool, value: Interval, operand: Interval) -> bool {
    if maximum {
        operand.hi() >= value.lo()
    } else {
        operand.lo() <= value.hi()
    }
}

fn extreme(maximum: bool, parts: &[Jet]) -> Jet {
    let value = parts
        .iter()
        .skip(1)
        .fold(parts[0].value, |a, b| extreme_value(maximum, a, b.value));
    let active = parts
        .iter()
        .filter(|j| can_be_active(maximum, value, j.value))
        .flat_map(|j| j.active.iter().copied())
        .collect();
    Jet { value, active }
}

fn operand_value(tape: &[Jet], index: usize, reversed: bool) -> Interval {
    let value = tape[index].value;
    if reversed { negative(value) } else { value }
}

fn join_blocks(a: &Box3, b: &Box3) -> Option<Box3> {
    let axis = (0..3).find(|&axis| a[axis] != b[axis])?;
    if (0..3).any(|i| i != axis && a[i] != b[i])
        || (a[axis][1] != b[axis][0] && b[axis][1] != a[axis][0])
    {
        return None;
    }
    let mut joined = a.clone();
    joined[axis] = [
        a[axis][0].clone().min(b[axis][0].clone()),
        a[axis][1].clone().max(b[axis][1].clone()),
    ];
    Some(joined)
}

// Each interval contains squared distance and every one-sided derivative
// through the interpolation. Floor indices are candidates only: the outward
// interval is used to include both branches when it overlaps a seam.
fn periodic(x: Interval, period: Interval, pitch: Interval) -> (Interval, Interval) {
    let turns = x / period + ip(0.5);
    let lo = turns.lo().floor();
    let hi = turns.hi().floor();
    let half = period * ip(0.5);
    let seam_gradient = Interval::new(-period.hi(), period.hi()).unwrap();
    let broad = || {
        (
            Interval::new(0.0, half.square().hi()).unwrap(),
            seam_gradient,
        )
    };
    if !lo.is_finite() || !hi.is_finite() || lo.abs().max(hi.abs()) >= (1u64 << 52) as f64 {
        return broad();
    }
    let at = |k| x - period * ip(k);
    if lo.total_cmp(&hi).is_eq() {
        let delta = at(lo);
        if turns.contains(lo) {
            return (delta.square(), seam_gradient);
        }
        let ends = (delta * ip(2.0))
            .hull(&(delta * ip(2.0) + Interval::new(-pitch.hi(), pitch.hi()).unwrap()));
        return (delta.square(), ends);
    }
    if (hi - lo).total_cmp(&1.0).is_eq() {
        let left = at(lo);
        let right = at(hi);
        let floor_left = if left.lo() > 0.0 {
            ip(left.lo()).square().lo()
        } else {
            0.0
        };
        let floor_right = if right.hi() < 0.0 {
            ip(right.hi()).square().lo()
        } else {
            0.0
        };
        return (
            Interval::new(floor_left.min(floor_right), half.square().hi()).unwrap(),
            seam_gradient,
        );
    }
    broad()
}

enum Prepared {
    Block([Range; 2]),
    Cylinder {
        centre: Range,
        radius2: Interval,
        height: Interval,
    },
    Cubic {
        centre: Range,
        radius2: Interval,
        period: Interval,
    },
}
impl Prepared {
    fn new(s: &Support) -> Self {
        let i = Interval::from_rational;
        match s {
            Support::Block(b) => Self::Block(std::array::from_fn(|end| {
                std::array::from_fn(|axis| i(&b[axis][end]))
            })),
            Support::Cylinder {
                centre,
                radius,
                height,
            } => Self::Cylinder {
                centre: centre.each_ref().map(i),
                radius2: i(radius).square(),
                height: i(height),
            },
            Support::Cubic {
                centre,
                radius,
                period,
            } => Self::Cubic {
                centre: centre.each_ref().map(i),
                radius2: i(radius).square(),
                period: i(period),
            },
        }
    }
    fn jet(&self, world: &Range, pitch: Interval, error: Interval) -> Jet {
        let radial = |values: &[(Interval, Interval); 3], a: usize, b: usize, radius2| {
            let value = values[a].0 + values[b].0 - radius2;
            let value = Interval::new(value.lo(), (value + error).hi()).unwrap();
            let mut gradient = [ip(0.0); 3];
            gradient[a] = values[a].1;
            gradient[b] = values[b].1;
            Jet::new(value, gradient, error)
        };
        match *self {
            Self::Block([low, high]) => {
                let planes: [Jet; 6] = std::array::from_fn(|i| {
                    let k = i / 2;
                    if i % 2 == 0 {
                        Jet::plane(low[k] - world[k], k, -1.0)
                    } else {
                        Jet::plane(world[k] - high[k], k, 1.0)
                    }
                });
                extreme(true, &planes)
            }
            Self::Cylinder {
                centre,
                radius2,
                height,
            } => {
                let local: Range = std::array::from_fn(|i| world[i] - centre[i]);
                let slack = Interval::new(-pitch.hi(), pitch.hi()).unwrap();
                let jets = local.map(|x| (x.square(), (x * ip(2.0)) + slack));
                extreme(
                    true,
                    &[
                        Jet::plane(negative(local[2]), 2, -1.0),
                        Jet::plane(local[2] - height, 2, 1.0),
                        radial(&jets, 0, 1, radius2),
                    ],
                )
            }
            Self::Cubic {
                centre,
                radius2,
                period,
            } => {
                let jets = std::array::from_fn(|i| periodic(world[i] - centre[i], period, pitch));
                extreme(
                    false,
                    &[
                        radial(&jets, 0, 1, radius2),
                        radial(&jets, 0, 2, radius2),
                        radial(&jets, 1, 2, radius2),
                    ],
                )
            }
        }
    }
}

fn separates(branches: &[Branch], distance: Interval, dec: &mut Decider) -> Result<bool, Refusal> {
    let mut fits = |dir: [f64; 3]| {
        dec.charge_work(branches.len() as u64)?;
        let v = dir.map(ip);
        let norm2 = v.iter().fold(ip(0.0), |s, a| s + a.square());
        Ok::<_, Refusal>(branches.iter().all(|branch| {
            let slope = (0..3).fold(ip(0.0), |s, i| s + branch.gradient[i] * v[i]);
            slope.lo() > 0.0
                && (slope * distance).square().lo() > (branch.motion.square() * norm2).hi()
        }))
    };
    let midpoint = |b: &Branch| b.gradient.map(|x| x.mid());
    let mut combined = [0.0; 3];
    for b in branches {
        let dir = midpoint(b);
        let norm = crate::metric::software_sqrt(dir.iter().map(|v| v * v).sum());
        if norm > 0.0 {
            for i in 0..3 {
                combined[i] += dir[i] / norm;
            }
        }
    }
    if fits(combined)? {
        return Ok(true);
    }
    for b in branches {
        if fits(midpoint(b))? {
            return Ok(true);
        }
    }
    for x in -1..=1 {
        for y in -1..=1 {
            for z in -1..=1 {
                if fits([f64::from(x), f64::from(y), f64::from(z)])? {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

impl Model<'_> {
    fn joined_blocks(&self) -> Vec<Option<Prepared>> {
        let mut boxes: Vec<Option<Box3>> = Vec::with_capacity(self.plan.instructions.len());
        let mut prepared = Vec::with_capacity(self.plan.instructions.len());
        for instruction in &self.plan.instructions {
            let joined = match *instruction {
                Instruction::Boolean {
                    kind: BooleanKind::Union,
                    a,
                    b,
                } => boxes[a]
                    .as_ref()
                    .zip(boxes[b].as_ref())
                    .and_then(|(a, b)| join_blocks(a, b)),
                _ => None,
            };
            prepared.push(joined.clone().map(|b| Prepared::new(&Support::Block(b))));
            boxes.push(match *instruction {
                Instruction::Leaf(i) => match &self.supports[i] {
                    Some(Support::Block(b)) => Some(b.clone()),
                    _ => None,
                },
                _ => joined,
            });
        }
        prepared
    }

    fn prove(&self, budget: &mut Budget<'_>, dec: &mut Decider) -> Result<(), Refusal> {
        let prepared: Vec<_> = self
            .supports
            .iter()
            .map(|s| s.as_ref().map(Prepared::new))
            .collect();
        let pitch = Interval::from_rational(&self.pitch);
        let error = radial_interpolation_error(pitch);
        let distance = Interval::from_rational(&self.displacement);
        let mut pending = vec![self.domain.each_ref().map(enclosure)];
        let instructions = &self.plan.instructions;
        let tape_work = instructions.len() as u64;
        dec.charge_work(tape_work)?;
        // A face-to-face block union is exactly a block throughout the homotopy.
        // Its shared face is interior, even where the original min field is zero.
        let joined = self.joined_blocks();
        let mut tape: Vec<Jet> = Vec::with_capacity(instructions.len());
        let mut needed = vec![[false; 2]; instructions.len()];
        let mut branches = Vec::new();
        while let Some(cell) = pending.pop() {
            budget.visit(dec)?;
            dec.charge_work(tape_work)?;
            tape.clear();
            for (index, instruction) in instructions.iter().enumerate() {
                if let Some(block) = &joined[index] {
                    tape.push(block.jet(&cell, pitch, error));
                    continue;
                }
                let jet = match *instruction {
                    Instruction::Empty => Jet::new(ip(1.0), [ip(0.0); 3], ip(0.0)),
                    Instruction::Leaf(i) => prepared[i]
                        .as_ref()
                        .expect("live support")
                        .jet(&cell, pitch, error),
                    Instruction::Boolean { kind, a, b } => Jet {
                        value: extreme_value(
                            kind != BooleanKind::Union,
                            tape[a].value,
                            operand_value(&tape, b, kind == BooleanKind::Subtract),
                        ),
                        active: Vec::new(),
                    },
                };
                tape.push(jet);
            }
            if !tape[self.plan.root].value.contains_zero() {
                continue;
            }
            // Walk the active DAG once per sign. Shared operands may reach
            // the root along exponentially many paths but have only two signs.
            dec.charge_work(tape_work)?;
            needed.fill([false; 2]);
            needed[self.plan.root][0] = true;
            branches.clear();
            for (index, instruction) in instructions.iter().enumerate().rev() {
                let [forward, reversed] = needed[index];
                if !forward && !reversed {
                    continue;
                }
                let jet = &tape[index];
                if let Instruction::Boolean { kind, a, b } = *instruction
                    && joined[index].is_none()
                {
                    for (operand, flip) in [(a, false), (b, kind == BooleanKind::Subtract)] {
                        if can_be_active(
                            kind != BooleanKind::Union,
                            jet.value,
                            operand_value(&tape, operand, flip),
                        ) {
                            needed[operand][usize::from(flip)] |= forward;
                            needed[operand][usize::from(!flip)] |= reversed;
                        }
                    }
                } else {
                    if forward {
                        branches.extend(jet.active.iter().copied());
                    }
                    if reversed {
                        branches.extend(jet.active.iter().map(|branch| branch.reversed()));
                    }
                }
            }
            if separates(&branches, distance, dec)? {
                continue;
            }
            let middle = cell.map(|axis| axis.mid());
            let spacing = Interval::from_rational(&self.encoding_reserve).hi();
            // Refinement below the output spacing cannot improve a file-space proof.
            if cell
                .iter()
                .zip(middle)
                .any(|(axis, mid)| !(axis.lo() < mid && mid < axis.hi()) || axis.width() <= spacing)
            {
                return Err(refusal(
                    "field.print-detail-topology",
                    format!("no common direction through the lattice homotopy at {cell:?}"),
                ));
            }
            for octant in (0..8).rev() {
                pending.push(std::array::from_fn(|i| {
                    let [lo, hi] = if octant & (1 << i) == 0 {
                        [cell[i].lo(), middle[i]]
                    } else {
                        [middle[i], cell[i].hi()]
                    };
                    Interval::new(lo, hi).unwrap()
                }));
            }
        }
        // The exact box carriers and shrinking radial carriers preserve the
        // original interval bounds of each positive operand throughout time.
        // Boolean containment therefore keeps every zero inside the padded box.
        // A smooth partition of local directions gives a common flow. Every
        // active branch rises along that flow faster than its time derivative
        // divided by the displacement budget. The zero of the Lipschitz field
        // moves uniquely on each flow line, retaining material on the same side.
        Ok(())
    }
}
