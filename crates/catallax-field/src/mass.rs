//! Certified unit-density moments of the field, independent of print meshes.
//!
//! Exact box moments cover admitted material. Unresolved boxes admit densities
//! in [0,1]; their coordinate moments give outward rational bounds. Refinement
//! replaces those contributions, never samples away an unresolved region.

use crate::{Plan, geometry, lattice::BoxState, metric::Kind, printing::Classifier};
use catallax_features::stand_in::MassTarget;
use catallax_measure::{Accuracy as QuantityAccuracy, Enclosure};
use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign};
use catallax_shapes::Surface;
use std::{
    cmp::Ordering,
    collections::{BTreeSet, BinaryHeap},
};

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
    /// Maximum classified boxes, including the initial plane partition.
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
    /// Moments of the admitted material plus half of each unresolved box.
    /// This positive-density estimate is not measured from a print mesh.
    pub estimate: Option<MassTarget>,
    /// Total classified boxes.
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

#[derive(Clone)]
struct Cell {
    bounds: Bounds,
    volume: Rational,
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
        self.volume
            .cmp(&other.volume)
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
        self.cells += 1;
        if self.cells > self.max_cells {
            return Err(failure(
                "field.mass-budget",
                "field moment cell limit reached before the requested Accuracy",
            ));
        }
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
        let cell = Cell {
            bounds,
            volume,
            serial: self.serial,
        };
        self.serial += 1;
        self.totals.update(
            &Contribution::new(&cell, &self.anchor, state == BoxState::Inside),
            false,
        );
        if state == BoxState::Unknown {
            self.pending.push(cell);
        }
        Ok(())
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
    let anchor = std::array::from_fn(|i| (&bounds[i][0] + &bounds[i][1]).mul_pow2(-1));
    let mut axes: [BTreeSet<Rational>; 3] =
        std::array::from_fn(|i| BTreeSet::from([bounds[i][0].clone(), bounds[i][1].clone()]));
    for face in &plan.faces {
        if plan.leaf_registers[face.leaf].is_none() {
            continue;
        }
        if let Surface::Plane(plane) = &face.surface {
            if let (Some(normal), Some(origin)) =
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
