//! Certified unit-density moments of the field, independent of print meshes.
//!
//! Analytic columns use exact inner antiderivatives and certified Taylor
//! remainders. Other definitions use exact admitted-box moments and unresolved
//! densities in [0,1]. Neither method samples away an unresolved region.

use crate::{Plan, geometry, lattice::BoxState, metric::Kind, printing::Classifier};
use catallax_features::stand_in::{MassTarget, StandInBodies};
use catallax_measure::{Accuracy as QuantityAccuracy, Enclosure};
use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign};
use catallax_shapes::Surface;
use std::{
    cmp::Ordering,
    collections::{BTreeSet, BinaryHeap},
};

mod columns;

type Bounds = [[Rational; 2]; 3];
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Requested widths, in the units of each measured quantity.
#[derive(Clone, Debug, Default)]
pub struct Accuracy {
    /// Volume accuracy in cubic millimetres.
    pub volume: QuantityAccuracy,
    /// Centre-of-gravity accuracy in millimetres.
    pub centroid: QuantityAccuracy,
    /// Centroidal inertia accuracy in millimetres to the fifth power.
    pub inertia: QuantityAccuracy,
}
/// Deterministic integration resource control.
#[derive(Clone, Debug)]
pub struct Options {
    /// Maximum classified boxes or analytic intervals, including initial partitions.
    pub max_cells: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            max_cells: 2_000_000,
        }
    }
}
/// Certified field properties and a rational estimate inside those enclosures.
#[derive(Clone, Debug)]
pub struct Properties {
    /// True field volume, not mesh volume.
    pub volume: Enclosure,
    /// True field centre of gravity; absent exactly for zero volume.
    pub centroid: Option<[Enclosure; 3]>,
    /// True field centroidal inertia; absent exactly for zero volume.
    pub inertia: Option<[[Enclosure; 3]; 3]>,
    /// Rational coordinate-moment proposal enclosed by the integration bounds.
    /// This estimate is independent of any print mesh.
    pub estimate: Option<MassTarget>,
    /// Total classified boxes or analytic intervals.
    pub cells: usize,
}
impl Properties {
    /// Whether a surrogate's exact properties lie in every field enclosure.
    pub fn contains(&self, target: &MassTarget) -> bool {
        self.volume.contains(&target.volume)
            && self
                .centroid
                .as_ref()
                .is_some_and(|c| (0..3).all(|i| c[i].contains(&target.centroid[i])))
            && self.inertia.as_ref().is_some_and(|tensor| {
                (0..3).all(|i| (0..3).all(|j| tensor[i][j].contains(&target.inertia[i][j])))
            })
    }
    /// Whether every exchanged surrogate property stays inside the field bounds.
    pub fn encloses(&self, target: &SurrogateProperties) -> bool {
        fn contains(outer: &Enclosure, inner: &Enclosure) -> bool {
            outer.lo <= inner.lo && inner.hi <= outer.hi
        }
        contains(&self.volume, &target.volume)
            && self
                .centroid
                .as_ref()
                .is_some_and(|c| (0..3).all(|i| contains(&c[i], &target.centroid[i])))
            && self.inertia.as_ref().is_some_and(|inertia| {
                (0..3).all(|i| (0..3).all(|j| contains(&inertia[i][j], &target.inertia[i][j])))
            })
    }
}

/// Unit-density proxy moments after a bounded convex-body exchange deviation.
#[derive(Clone, Debug)]
pub struct SurrogateProperties {
    /// Proxy volume enclosure.
    pub volume: Enclosure,
    /// Proxy centre-of-gravity enclosures.
    pub centroid: [Enclosure; 3],
    /// Proxy centroidal inertia enclosures.
    pub inertia: [[Enclosure; 3]; 3],
}

/// Encloses the exchange of the solver's convex boxes, independently of a mesh.
///
/// Both convex bodies fit a cube of side `2(radius + deviation)`. Monotonicity
/// of the Steiner coefficients bounds each outer parallel layer by the cube's
/// layer. Twice that bound covers the symmetric difference. Integrating its
/// coordinate bounds covers volume, first and second moments; division by a
/// positive volume lower bound then covers the translated centroid and inertia.
pub fn surrogate_properties(
    surrogate: &StandInBodies,
    deviation: &Rational,
) -> Result<SurrogateProperties, Refusal> {
    if deviation.sign() == Sign::Negative {
        return Err(failure(
            "field.surrogate-deviation",
            "exchange deviation must be nonnegative",
        ));
    }
    let radius = &surrogate.support_radius + deviation;
    let diameter = radius.mul_pow2(1);
    let square = deviation * deviation;
    let layer =
        &(&(&q(6) * &(&diameter * &diameter)) * deviation) + &(&(&q(12) * &diameter) * &square);
    let layer = &layer + &(&q(8) * &(&square * deviation));
    let difference = &layer * &q(2 * surrogate.bodies.len() as i64);
    let mass = &surrogate.achieved.volume;
    let lower = mass - &difference;
    if lower.sign() != Sign::Positive {
        return Err(failure(
            "field.surrogate-deviation",
            "exchange uncertainty does not preserve a positive proxy volume",
        ));
    }
    let upper = mass + &difference;
    let centroid_error = (&difference * &radius)
        .checked_div(&lower)
        .expect("positive volume lower bound");
    let second_error = &difference * &(&radius * &radius);
    let shifted_error = &upper * &(&centroid_error * &centroid_error);
    let inertia_error = (&second_error + &shifted_error).mul_pow2(1);
    Ok(SurrogateProperties {
        volume: Enclosure {
            lo: lower,
            hi: upper,
        },
        centroid: std::array::from_fn(|i| Enclosure {
            lo: &surrogate.achieved.centroid[i] - &centroid_error,
            hi: &surrogate.achieved.centroid[i] + &centroid_error,
        }),
        inertia: std::array::from_fn(|i| {
            std::array::from_fn(|j| Enclosure {
                lo: &surrogate.achieved.inertia[i][j] - &inertia_error,
                hi: &surrogate.achieved.inertia[i][j] + &inertia_error,
            })
        }),
    })
}

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn failure(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}
fn exact(value: Rational) -> Enclosure {
    Enclosure {
        lo: value.clone(),
        hi: value,
    }
}
fn add(a: &Enclosure, b: &Enclosure) -> Enclosure {
    Enclosure {
        lo: &a.lo + &b.lo,
        hi: &a.hi + &b.hi,
    }
}
fn subtract(a: &Enclosure, b: &Enclosure) -> Enclosure {
    Enclosure {
        lo: &a.lo - &b.hi,
        hi: &a.hi - &b.lo,
    }
}
fn multiply(a: &Enclosure, b: &Enclosure) -> Enclosure {
    let products = [&a.lo * &b.lo, &a.lo * &b.hi, &a.hi * &b.lo, &a.hi * &b.hi];
    Enclosure {
        lo: products.iter().min().unwrap().clone(),
        hi: products.iter().max().unwrap().clone(),
    }
}
fn divide(a: &Enclosure, b: &Enclosure) -> Enclosure {
    let reciprocal = Enclosure {
        lo: q(1).checked_div(&b.hi).unwrap(),
        hi: q(1).checked_div(&b.lo).unwrap(),
    };
    multiply(a, &reciprocal)
}
fn reached(value: &Enclosure, accuracy: &QuantityAccuracy) -> bool {
    let magnitude = if value.lo.sign() == value.hi.sign() {
        value.lo.abs().min(value.hi.abs())
    } else {
        q(0)
    };
    value.width() <= (&accuracy.relative * &magnitude).max(accuracy.absolute.clone())
}
fn valid(accuracy: &QuantityAccuracy) -> bool {
    accuracy.relative.sign() != Sign::Negative && accuracy.absolute.sign() != Sign::Negative
}

/// Checks that replacing lattice material by its region is an outer envelope.
/// Subtraction reverses containment, so negative lattice occurrences need a
/// separate lower-envelope construction and are not covered by that exchange.
pub fn check_region_envelope(plan: &Plan) -> Result<(), Refusal> {
    use crate::Instruction;
    let mut signs = vec![0u8; plan.instructions.len()];
    signs[plan.root] = 1;
    // Compiled tape operands precede their Boolean instruction.
    for node in (0..plan.instructions.len()).rev() {
        let sign = signs[node];
        if sign == 0 {
            continue;
        }
        match plan.instructions[node] {
            Instruction::Leaf(leaf)
                if sign & 2 != 0 && matches!(plan.leaves[leaf], geometry::Leaf::Lattice(_)) =>
            {
                return Err(failure(
                    "field.stand-in-envelope-coverage",
                    "region substitution is not an outer-envelope proof for a negative lattice operand",
                ));
            }
            Instruction::Boolean { kind, a, b } => {
                signs[a] |= sign;
                signs[b] |= if kind == catallax_recipe::BooleanKind::Subtract {
                    ((sign & 1) << 1) | ((sign & 2) >> 1)
                } else {
                    sign
                };
            }
            _ => {}
        }
    }
    Ok(())
}

#[derive(Clone)]
struct Cell {
    bounds: Bounds,
    volume: Rational,
    priority: Rational,
    serial: u64,
}
impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.priority
            .cmp(&other.priority)
            .then_with(|| other.serial.cmp(&self.serial))
    }
}
struct Contribution {
    lo: [Rational; 10],
    hi: [Rational; 10],
    estimate: [Rational; 10],
}
impl Contribution {
    fn new(cell: &Cell, anchor: &[Rational; 3], inside: bool) -> Self {
        let local: Bounds = std::array::from_fn(|i| {
            [
                &cell.bounds[i][0] - &anchor[i],
                &cell.bounds[i][1] - &anchor[i],
            ]
        });
        let centers: [Rational; 3] =
            std::array::from_fn(|i| (&local[i][0] + &local[i][1]).mul_pow2(-1));
        let widths: [Rational; 3] = std::array::from_fn(|i| &local[i][1] - &local[i][0]);
        let moments: [Rational; 10] = std::array::from_fn(|index| match index {
            0 => cell.volume.clone(),
            1..=3 => &cell.volume * &centers[index - 1],
            _ => {
                let (i, j) = PAIRS[index - 4];
                let mut product = &centers[i] * &centers[j];
                if i == j {
                    product = product + (&widths[i] * &widths[i]).checked_div(&q(12)).unwrap();
                }
                &cell.volume * &product
            }
        });
        if inside {
            return Self {
                lo: moments.clone(),
                hi: moments.clone(),
                estimate: moments,
            };
        }
        let ranges: [Enclosure; 10] = std::array::from_fn(|index| {
            let range = match index {
                0 => {
                    return Enclosure {
                        lo: q(0),
                        hi: cell.volume.clone(),
                    };
                }
                1..=3 => Enclosure {
                    lo: local[index - 1][0].clone(),
                    hi: local[index - 1][1].clone(),
                },
                _ => {
                    let (i, j) = PAIRS[index - 4];
                    if i == j {
                        let squares = [&local[i][0] * &local[i][0], &local[i][1] * &local[i][1]];
                        Enclosure {
                            lo: if local[i][0] <= q(0) && local[i][1] >= q(0) {
                                q(0)
                            } else {
                                squares[0].clone().min(squares[1].clone())
                            },
                            hi: squares[0].clone().max(squares[1].clone()),
                        }
                    } else {
                        multiply(
                            &Enclosure {
                                lo: local[i][0].clone(),
                                hi: local[i][1].clone(),
                            },
                            &Enclosure {
                                lo: local[j][0].clone(),
                                hi: local[j][1].clone(),
                            },
                        )
                    }
                }
            };
            Enclosure {
                lo: (&cell.volume * &range.lo).min(q(0)),
                hi: (&cell.volume * &range.hi).max(q(0)),
            }
        });
        Self {
            lo: ranges.each_ref().map(|r| r.lo.clone()),
            hi: ranges.each_ref().map(|r| r.hi.clone()),
            estimate: moments.map(|v| v.mul_pow2(-1)),
        }
    }
}
struct Totals {
    lo: [Rational; 10],
    hi: [Rational; 10],
    estimate: [Rational; 10],
}
impl Totals {
    fn new() -> Self {
        Self {
            lo: std::array::from_fn(|_| q(0)),
            hi: std::array::from_fn(|_| q(0)),
            estimate: std::array::from_fn(|_| q(0)),
        }
    }
    fn update(&mut self, value: &Contribution, remove: bool) {
        for (i, (lo, hi)) in self.lo.iter_mut().zip(&mut self.hi).enumerate() {
            if remove {
                *lo = &*lo - &value.lo[i];
                *hi = &*hi - &value.hi[i];
                self.estimate[i] = &self.estimate[i] - &value.estimate[i];
            } else {
                *lo = &*lo + &value.lo[i];
                *hi = &*hi + &value.hi[i];
                self.estimate[i] = &self.estimate[i] + &value.estimate[i];
            }
        }
    }
    fn properties(&self, anchor: &[Rational; 3], cells: usize) -> Option<Properties> {
        let volume = Enclosure {
            lo: self.lo[0].clone(),
            hi: self.hi[0].clone(),
        };
        if volume.hi.is_zero() {
            return Some(Properties {
                volume,
                centroid: None,
                inertia: None,
                estimate: None,
                cells,
            });
        }
        if volume.lo <= q(0) {
            return None;
        }
        let means: [Enclosure; 3] = std::array::from_fn(|i| {
            divide(
                &Enclosure {
                    lo: self.lo[i + 1].clone(),
                    hi: self.hi[i + 1].clone(),
                },
                &volume,
            )
        });
        let centroid = std::array::from_fn(|i| add(&means[i], &exact(anchor[i].clone())));
        let covariance: [[Enclosure; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let index = 4 + PAIRS
                    .iter()
                    .position(|&(a, b)| (a == i && b == j) || (a == j && b == i))
                    .unwrap();
                subtract(
                    &Enclosure {
                        lo: self.lo[index].clone(),
                        hi: self.hi[index].clone(),
                    },
                    &multiply(&volume, &multiply(&means[i], &means[j])),
                )
            })
        });
        let inertia = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    add(
                        &covariance[(i + 1) % 3][(i + 1) % 3],
                        &covariance[(i + 2) % 3][(i + 2) % 3],
                    )
                } else {
                    Enclosure {
                        lo: -&covariance[i][j].hi,
                        hi: -&covariance[i][j].lo,
                    }
                }
            })
        });
        let estimated_volume = &self.estimate[0];
        let mean: [Rational; 3] =
            std::array::from_fn(|i| self.estimate[i + 1].checked_div(estimated_volume).unwrap());
        let estimated_covariance: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let index = 4 + PAIRS
                    .iter()
                    .position(|&(a, b)| (a == i && b == j) || (a == j && b == i))
                    .unwrap();
                &self.estimate[index] - &(estimated_volume * &(&mean[i] * &mean[j]))
            })
        });
        let estimate = MassTarget {
            volume: estimated_volume.clone(),
            centroid: std::array::from_fn(|i| &mean[i] + &anchor[i]),
            inertia: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    if i == j {
                        &estimated_covariance[(i + 1) % 3][(i + 1) % 3]
                            + &estimated_covariance[(i + 2) % 3][(i + 2) % 3]
                    } else {
                        -&estimated_covariance[i][j]
                    }
                })
            }),
        };
        Some(Properties {
            volume,
            centroid: Some(centroid),
            inertia: Some(inertia),
            estimate: Some(estimate),
            cells,
        })
    }
}
fn axis_block(primitive: &geometry::Primitive) -> Option<Bounds> {
    if primitive.kind != Kind::Block {
        return None;
    }
    let map = primitive.inverse_placement();
    let matrix = map.rational_matrix()?;
    let translation = map.translation.as_rational()?;
    let offset = [translation.x, translation.y, translation.z];
    let mut result: Bounds = std::array::from_fn(|_| [q(0), q(0)]);
    let mut used = [false; 3];
    for (i, row) in matrix.iter().enumerate() {
        let axis = row.iter().position(|n| !n.is_zero())?;
        if used[axis]
            || row
                .iter()
                .enumerate()
                .any(|(j, n)| j != axis && !n.is_zero())
        {
            return None;
        }
        used[axis] = true;
        result[axis] = if *row[axis] == q(1) {
            [-&offset[i], &primitive.parameters[i] - &offset[i]]
        } else if *row[axis] == q(-1) {
            [&offset[i] - &primitive.parameters[i], offset[i].clone()]
        } else {
            return None;
        };
    }
    Some(result)
}
fn closed_box(cell: &Bounds, solid: &Bounds) -> Option<BoxState> {
    if (0..3).any(|i| cell[i][1] <= solid[i][0] || cell[i][0] >= solid[i][1]) {
        Some(BoxState::Outside)
    } else if (0..3).all(|i| cell[i][0] >= solid[i][0] && cell[i][1] <= solid[i][1]) {
        Some(BoxState::Inside)
    } else {
        None
    }
}
struct Integration<'a> {
    classifier: Classifier<'a>,
    exact_boxes: Vec<Option<Bounds>>,
    anchor: [Rational; 3],
    totals: Totals,
    pending: BinaryHeap<Cell>,
    cells: usize,
    serial: u64,
    max_cells: usize,
}
impl Integration<'_> {
    fn admit(&mut self, bounds: Bounds, dec: &mut Decider) -> Result<(), Refusal> {
        if self.cells == self.max_cells {
            return Err(failure(
                "field.mass-budget",
                "field moment cell limit reached before the requested Accuracy",
            ));
        }
        self.cells += 1;
        let cover = std::array::from_fn(|i| {
            Interval::new(
                Interval::from_rational(&bounds[i][0]).lo(),
                Interval::from_rational(&bounds[i][1]).hi(),
            )
            .unwrap()
        });
        let boxes = &self.exact_boxes;
        let state = self.classifier.classify_with(&cover, dec, |i| {
            boxes[i]
                .as_ref()
                .and_then(|solid| closed_box(&bounds, solid))
        })?;
        if state == BoxState::Outside {
            return Ok(());
        }
        let volume = (0..3).fold(q(1), |v, i| &v * &(&bounds[i][1] - &bounds[i][0]));
        let mut cell = Cell {
            bounds,
            volume,
            priority: q(0),
            serial: self.serial,
        };
        self.serial += 1;
        let contribution = Contribution::new(&cell, &self.anchor, state == BoxState::Inside);
        self.totals.update(&contribution, false);
        if state == BoxState::Unknown {
            // Refine the largest unresolved coordinate-moment enclosure,
            // rather than treating equal-volume boxes as equally important.
            cell.priority = contribution
                .lo
                .iter()
                .zip(&contribution.hi)
                .map(|(lo, hi)| hi - lo)
                .max()
                .expect("ten coordinate moments");
            self.pending.push(cell);
        }
        Ok(())
    }
    fn budget_failure(&self) -> Refusal {
        let detail = if let Some(properties) = self.totals.properties(&self.anchor, self.cells) {
            let centroid = properties
                .centroid
                .as_ref()
                .map(|c| c.iter().map(Enclosure::width).max().unwrap());
            let inertia = properties
                .inertia
                .as_ref()
                .map(|m| m.iter().flatten().map(Enclosure::width).max().unwrap());
            format!(
                "{} classified cells; reached volume width {} mm3, maximum centroid width {:?} mm, maximum inertia width {:?} mm5; requested Accuracy not reached",
                self.cells,
                properties.volume.width(),
                centroid.map(|q| q.to_string()),
                inertia.map(|q| q.to_string())
            )
        } else {
            format!(
                "{} classified cells; positive field volume is not yet certified",
                self.cells
            )
        };
        failure("field.mass-budget", detail)
    }
}
/// Integrates the original field to the requested property widths.
/// Planar subdivision walls have zero volume; exact closed-box admission avoids
/// charging their arbitrary floating-point neighbourhood to mass uncertainty.
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
            "nonnegative property Accuracy values and a positive cell limit are required",
        ));
    }
    if plan.is_proven_empty() {
        return Ok(Properties {
            volume: exact(q(0)),
            centroid: None,
            inertia: None,
            estimate: None,
            cells: 0,
        });
    }
    let bounds = crate::printing::bounds(plan, dec)?;
    if bounds.iter().any(|b| b[1] <= b[0]) {
        return Ok(Properties {
            volume: exact(q(0)),
            centroid: None,
            inertia: None,
            estimate: None,
            cells: 0,
        });
    }
    if let Some(properties) = columns::properties(plan, &bounds, accuracy, options, dec) {
        return properties;
    }
    let anchor = std::array::from_fn(|i| (&bounds[i][0] + &bounds[i][1]).mul_pow2(-1));
    let mut axes: [BTreeSet<Rational>; 3] =
        std::array::from_fn(|i| BTreeSet::from([bounds[i][0].clone(), bounds[i][1].clone()]));
    for face in &plan.faces {
        if plan.leaf_registers[face.leaf].is_none() {
            continue;
        }
        if let Surface::Plane(plane) = &face.surface
            && let (Some(normal), Some(origin)) =
                (plane.normal.as_rational(), plane.origin.as_rational())
        {
            let normal = [normal.x, normal.y, normal.z];
            let origin = [origin.x, origin.y, origin.z];
            if normal.iter().filter(|n| !n.is_zero()).count() == 1 {
                let axis = normal.iter().position(|n| !n.is_zero()).unwrap();
                let coordinate = &origin[axis];
                if coordinate > &bounds[axis][0] && coordinate < &bounds[axis][1] {
                    axes[axis].insert(coordinate.clone());
                }
            }
        }
    }
    let axes = axes.map(|a| a.into_iter().collect::<Vec<_>>());
    let exact_boxes = plan
        .leaves
        .iter()
        .map(|leaf| match leaf {
            geometry::Leaf::Primitive(p) => axis_block(p),
            _ => None,
        })
        .collect();
    let mut integration = Integration {
        classifier: Classifier::new(plan, dec)?,
        exact_boxes,
        anchor,
        totals: Totals::new(),
        pending: BinaryHeap::new(),
        cells: 0,
        serial: 0,
        max_cells: options.max_cells,
    };
    for x in axes[0].windows(2) {
        for y in axes[1].windows(2) {
            for z in axes[2].windows(2) {
                integration.admit(
                    [
                        [x[0].clone(), x[1].clone()],
                        [y[0].clone(), y[1].clone()],
                        [z[0].clone(), z[1].clone()],
                    ],
                    dec,
                )?;
            }
        }
    }
    loop {
        if let Some(result) = integration
            .totals
            .properties(&integration.anchor, integration.cells)
        {
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
                "exact field moment propagation did not meet the requested Accuracy",
            ));
        }
        for _ in 0..128 {
            if !integration.pending.is_empty()
                && integration.cells.saturating_add(2) > integration.max_cells
            {
                // Keep the parent's complete enclosure until both children
                // can be admitted, including on a resource-limit refusal.
                return Err(integration.budget_failure());
            }
            let Some(cell) = integration.pending.pop() else {
                break;
            };
            integration
                .totals
                .update(&Contribution::new(&cell, &integration.anchor, false), true);
            let mut axis = 0;
            let mut longest = &cell.bounds[0][1] - &cell.bounds[0][0];
            for i in 1..3 {
                let width = &cell.bounds[i][1] - &cell.bounds[i][0];
                if width > longest {
                    axis = i;
                    longest = width;
                }
            }
            let middle = (&cell.bounds[axis][0] + &cell.bounds[axis][1]).mul_pow2(-1);
            let mut lower = cell.bounds.clone();
            let mut upper = cell.bounds;
            lower[axis][1] = middle.clone();
            upper[axis][0] = middle;
            integration.admit(lower, dec)?;
            integration.admit(upper, dec)?;
        }
    }
}
