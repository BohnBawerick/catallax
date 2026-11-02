//! Bezier patches of a rational sheet over a `u, v` chart with quadratic `v`:
//! subdivision, the display grid and the complete line solver. A sheet's own
//! law supplies only its rows and the exact equations of a span.
//!
//! ## Display grid
//!
//! Each cell is one subdivided rational Bezier patch P / W. Against the bilinear
//! map B through its corners, S - B = (P - B W) / W. Both P and B W are polynomials
//! of one degree higher in each direction, so the Bernstein coefficients of the
//! numerator over the smallest weight bound the distance at matching parameters.
//! That bound follows the patch's real curvature, whereas comparing rational
//! controls with their Greville points adds a term from the weights' own slope.
//! The quad's split then adds its twist.
//!
//! ## Line solver
//!
//! Over the chart, two planes through the line give a continuous piecewise
//! polynomial map G whose zeros are the line's hits. Bernstein hulls discard
//! each cell that keeps G from zero, so every hit stays in a surviving cell.
//! Surviving cells join into rectangles. Where a preconditioner Y gives
//! |I - Y J| < 1 over a rectangle's interval Jacobian, G is injective there:
//! the mean of J along any segment, knots included, stays in that hull. Such a
//! rectangle owns at most one hit, and its Jacobian keeps one sign, so the line
//! crosses the carrier there. A rectangle with no such Y refuses.
//!
//! The mean-value image of a rectangle then decides its hit. Clear of every
//! knot and row, the hit gets a certified root of the span's own equations. On
//! a knot or row it is an exact point of that section, or the image narrows off
//! it. Nothing is accepted from a nearby enclosure alone.

use crate::rolling_envelope::{plus, times};
use crate::{BSplineSurface, Point3, Vector3, math};
use catallax_numbers::{
    Ball, BoxRoots, Decider, Outcome, PolySystem, Rational, Real, RealPoly, Refusal, RootBox, Sign,
};
use std::sync::Arc;

// Proposal bounds, not tolerances: running out refuses and proves nothing.
const MAX_DEPTH: usize = 32;
const MAX_CELLS: usize = 256;
const MAX_NARROWING: usize = 16;

fn unresolved(detail: impl Into<String>) -> Refusal {
    Refusal::new("numbers.geometric-gap-unproven", detail)
}
fn b(value: &Rational, bits: u32) -> Ball {
    Ball::from_rational(value, bits)
}
fn q(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}
fn meets(a: &Ball, c: &Ball) -> bool {
    a.lo() <= c.hi() && c.lo() <= a.hi()
}
fn dot(a: &[Ball; 3], c: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &c[0]) + &(&a[1] * &c[1])) + &(&a[2] * &c[2])
}

// Enclosed frame of one placed copy.
pub(crate) struct LocalFrame {
    pub(crate) origin: [Ball; 3],
    pub(crate) axes: [[Ball; 3]; 3],
}

impl LocalFrame {
    pub(crate) fn new(
        origin: &Point3,
        axes: [&Vector3; 3],
        bits: u32,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let vector = |v: &Vector3| -> Result<[Ball; 3], Refusal> {
            Ok([
                dec.enclose_at(&v.x, bits)?,
                dec.enclose_at(&v.y, bits)?,
                dec.enclose_at(&v.z, bits)?,
            ])
        };
        Ok(Self {
            origin: math::point(origin, bits, dec)?,
            axes: [vector(axes[0])?, vector(axes[1])?, vector(axes[2])?],
        })
    }
    // A query's frame: charged to the caller and bound by its precision policy.
    pub(crate) fn budgeted(
        origin: &Point3,
        axes: [&Vector3; 3],
        bits: u32,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            origin: math::budgeted_point(origin, bits, dec)?,
            axes: [
                math::budgeted_vector(axes[0], bits, dec)?,
                math::budgeted_vector(axes[1], bits, dec)?,
                math::budgeted_vector(axes[2], bits, dec)?,
            ],
        })
    }
    pub(crate) fn local(&self, point: &[Ball; 3]) -> [Ball; 3] {
        let d: [Ball; 3] = std::array::from_fn(|i| &point[i] - &self.origin[i]);
        self.axes.each_ref().map(|axis| dot(&d, axis))
    }
    pub(crate) fn world(&self, local: &[Ball; 3]) -> [Ball; 3] {
        std::array::from_fn(|i| {
            (0..3).fold(self.origin[i].clone(), |sum, axis| {
                sum + &self.axes[axis][i] * &local[axis]
            })
        })
    }
}

/// One tensor Bezier patch of a sheet, with homogeneous enclosed controls.
#[derive(Clone)]
pub(crate) struct Patch {
    pub(crate) degree: usize,
    pub(crate) cell: [[Rational; 2]; 2],
    pub(crate) net: Vec<[Ball; 4]>,
}
impl Patch {
    pub(crate) fn all(
        carrier: &BSplineSurface,
        controls: &[[Ball; 4]],
        degree: usize,
    ) -> Vec<Self> {
        carrier
            .spans(0)
            .enumerate()
            .map(|(index, span)| Self {
                degree,
                cell: [
                    [span[0].clone(), span[1].clone()],
                    [Rational::zero(), Rational::one()],
                ],
                net: controls[index * degree * 3..(index * degree + degree + 1) * 3].to_vec(),
            })
            .collect()
    }
    // The same patches with their stored enclosures rounded outward to `bits`.
    pub(crate) fn at(
        carrier: &BSplineSurface,
        controls: &[[Ball; 4]],
        stored: u32,
        degree: usize,
        bits: u32,
    ) -> Vec<Self> {
        let mut patches = Self::all(carrier, controls, degree);
        if bits < stored {
            for patch in &mut patches {
                for h in &mut patch.net {
                    *h = std::array::from_fn(|i| h[i].with_precision(bits));
                }
            }
        }
        patches
    }
    pub(crate) fn split(self) -> [Self; 2] {
        // A span's u extent is about a tenth of v's; compare in a metric that
        // roughly equalises the two chart speeds.
        let du = &self.cell[0][1] - &self.cell[0][0];
        let dv = (&self.cell[1][1] - &self.cell[1][0]) / q(8);
        let axis = usize::from(dv > du);
        let half = (&self.cell[axis][0] + &self.cell[axis][1]) / q(2);
        let mut low = self.cell.clone();
        let mut high = self.cell;
        low[axis][1] = half.clone();
        high[axis][0] = half;
        let [left, right] = split_net(&self.net, [self.degree, 2], axis);
        [
            Self {
                degree: self.degree,
                cell: low,
                net: left,
            },
            Self {
                degree: self.degree,
                cell: high,
                net: right,
            },
        ]
    }
    pub(crate) fn cartesian(&self) -> Result<Vec<[Ball; 3]>, Refusal> {
        self.net
            .iter()
            .map(|h| {
                let w = &h[3];
                if w.lo() <= Rational::zero() {
                    return Err(unresolved("subdivided sheet weight is not positive"));
                }
                Ok([0, 1, 2].map(|i| h[i].checked_div(w).expect("positive weight")))
            })
            .collect()
    }
    pub(crate) fn hull(&self) -> Result<[Ball; 3], Refusal> {
        let points = self.cartesian()?;
        Ok(std::array::from_fn(|i| {
            points[1..]
                .iter()
                .fold(points[0][i].clone(), |hull, p| hull.hull(&p[i]))
        }))
    }
}

impl Patch {
    // The sub-patch over a rectangle inside this cell, by two De Casteljau cuts
    // along each axis. A degenerate rectangle gives the net of one point.
    pub(crate) fn clip(&self, rect: &Rect, bits: u32) -> Self {
        let mut net = self.net.clone();
        for (axis, [from, to]) in rect.iter().enumerate() {
            let [lo, hi] = &self.cell[axis];
            let cut = |net: &[[Ball; 4]], ratio: Rational| {
                let right = b(&ratio, bits);
                let left = b(&(Rational::one() - ratio), bits);
                divide_net(
                    net,
                    [self.degree, 2],
                    axis,
                    |a: &[Ball; 4], c: &[Ball; 4]| {
                        std::array::from_fn(|i| &(&a[i] * &left) + &(&c[i] * &right))
                    },
                )
            };
            if from > lo {
                let [_, right] = cut(&net, (from - lo) / (hi - lo));
                net = right;
            }
            if to < hi && hi > from {
                let [left, _] = cut(&net, (to - from) / (hi - from));
                net = left;
            }
        }
        Self {
            degree: self.degree,
            cell: rect.clone(),
            net,
        }
    }
}

// De Casteljau halving of a tensor Bernstein net along one axis.
pub(crate) fn split_net<T: Clone + Halve>(
    net: &[T],
    degrees: [usize; 2],
    axis: usize,
) -> [Vec<T>; 2] {
    divide_net(net, degrees, axis, T::halve)
}

// De Casteljau division along one axis, with `mix` the step between neighbours.
pub(crate) fn divide_net<T: Clone>(
    net: &[T],
    degrees: [usize; 2],
    axis: usize,
    mix: impl Fn(&T, &T) -> T,
) -> [Vec<T>; 2] {
    let [m, n] = degrees;
    let index = |i: usize, j: usize| i * (n + 1) + j;
    let mut left = net.to_vec();
    let mut right = net.to_vec();
    let (lines, length) = if axis == 0 {
        (n + 1, m + 1)
    } else {
        (m + 1, n + 1)
    };
    for line in 0..lines {
        let at = |k: usize| {
            if axis == 0 {
                index(k, line)
            } else {
                index(line, k)
            }
        };
        let mut row: Vec<T> = (0..length).map(|k| net[at(k)].clone()).collect();
        let mut low = Vec::with_capacity(length);
        let mut high = vec![row[length - 1].clone(); length];
        for level in 0..length {
            low.push(row[0].clone());
            high[length - 1 - level] = row[length - 1 - level].clone();
            row = (0..row.len() - 1)
                .map(|k| mix(&row[k], &row[k + 1]))
                .collect();
            if row.is_empty() {
                break;
            }
        }
        for k in 0..length {
            left[at(k)] = low[k].clone();
            right[at(k)] = high[k].clone();
        }
    }
    [left, right]
}

pub(crate) trait Halve {
    fn halve(&self, other: &Self) -> Self;
}
impl Halve for Ball {
    fn halve(&self, other: &Self) -> Self {
        &(self + other)
            * &b(
                &Rational::from_ratio(1, 2).expect("nonzero"),
                self.precision(),
            )
    }
}
impl Halve for [Ball; 4] {
    fn halve(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i].halve(&other[i]))
    }
}

/// Grid over the chart, u-major. Boundary rows are the carrier's own rows.
#[derive(Clone, Debug)]
pub struct EnvelopeGrid {
    /// Grid lines along `u`.
    pub columns: usize,
    /// Grid lines along `v`.
    pub rows: usize,
    /// Enclosed world points, `columns * rows` of them.
    pub points: Vec<[Ball; 3]>,
    /// Two-way bound between the carrier and the grid's triangles, in mm.
    pub chord: Rational,
}

fn l1(a: &[Ball; 3], c: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + magnitude(&(&a[i] - &c[i])))
}

pub(crate) fn along(cell: Vec<Patch>, axis: usize) -> Vec<Patch> {
    cell.into_iter()
        .flat_map(|patch| {
            let half = (&patch.cell[axis][0] + &patch.cell[axis][1]) / q(2);
            let mut low = patch.cell.clone();
            let mut high = patch.cell.clone();
            low[axis][1] = half.clone();
            high[axis][0] = half;
            let [left, right] = split_net(&patch.net, [patch.degree, 2], axis);
            [
                Patch {
                    degree: patch.degree,
                    cell: low,
                    net: left,
                },
                Patch {
                    degree: patch.degree,
                    cell: high,
                    net: right,
                },
            ]
        })
        .collect()
}

impl Patch {
    fn corners(&self) -> Result<[[Ball; 3]; 4], Refusal> {
        let degree = self.degree;
        let at = |i: usize, j: usize| -> Result<[Ball; 3], Refusal> {
            let h = &self.net[i * 3 + j];
            let divide = |value: &Ball| {
                value
                    .checked_div(&h[3])
                    .ok_or_else(|| unresolved("subdivided sheet weight is not positive"))
            };
            Ok([divide(&h[0])?, divide(&h[1])?, divide(&h[2])?])
        };
        Ok([at(0, 0)?, at(degree, 0)?, at(0, 2)?, at(degree, 2)?])
    }
    fn chord(&self) -> Result<Rational, Refusal> {
        let degree = self.degree;
        let [a, c, d, e] = self.corners()?;
        let bits = a[0].precision();
        let zero: [Ball; 3] = std::array::from_fn(|_| b(&Rational::zero(), bits));
        // (1 - t) and t raise a Bernstein coefficient i of degree n to i and i + 1.
        let lifted = |degree: usize, step: usize, i: usize| {
            (if step == 0 { degree + 1 - i } else { i + 1 }) as i64
        };
        let corners = [[&a, &d], [&c, &e]];
        let mut numerator = Rational::zero();
        for k in 0..=degree + 1 {
            for l in 0..=3usize {
                let mut sum = zero.clone();
                for (step_u, along_v) in corners.iter().enumerate() {
                    let Some(i) = k.checked_sub(step_u).filter(|i| *i <= degree) else {
                        continue;
                    };
                    for (step_v, corner) in along_v.iter().enumerate() {
                        let Some(j) = l.checked_sub(step_v).filter(|j| *j <= 2) else {
                            continue;
                        };
                        let h = &self.net[i * 3 + j];
                        let factor = b(
                            &Rational::from_ratio(
                                lifted(degree, step_u, i) * lifted(2, step_v, j),
                                (degree as i64 + 1) * 3,
                            )
                            .expect("degree"),
                            bits,
                        );
                        for (m, value) in sum.iter_mut().enumerate() {
                            *value = &*value + &(&(&h[m] - &(&corner[m] * &h[3])) * &factor);
                        }
                    }
                }
                numerator = numerator.max(l1(&sum, &zero));
            }
        }
        let low = self.net.iter().fold(None::<Rational>, |low, h| {
            let w = h[3].lo();
            Some(low.map_or(w.clone(), |value| value.min(w)))
        });
        let low = low.expect("controls");
        if low <= Rational::zero() {
            return Err(unresolved("subdivided sheet weight is not positive"));
        }
        let twist: [Ball; 3] = std::array::from_fn(|k| {
            &(&(&(&a[k] - &c[k]) - &d[k]) + &e[k]) * &b(&Rational::from_ratio(1, 4).unwrap(), bits)
        });
        Ok(numerator / low + l1(&twist, &zero))
    }
}

/// Subdivides `patches`, one per span along `u`, until the chord is at most
/// `tolerance`. Doubling keeps every knot on a grid line. A tolerance out of reach
/// after 2^10 cells per span direction refuses instead of widening the report.
pub(crate) fn sheet_grid(
    patches: Vec<Patch>,
    tolerance: &Rational,
    max_points: usize,
    dec: &Decider,
) -> Result<EnvelopeGrid, Refusal> {
    // Two cuts across start each meridian below a sixteenth of a turn per cell.
    let mut patches = along(along(patches, 1), 1);
    let spans = patches.len() / 4;
    let (mut cuts_u, mut cuts_v) = (0u32, 2u32);
    let (patches, chord) = loop {
        dec.check_deadline()?;
        let mut worst = Rational::zero();
        let mut worst_shape = [Rational::zero(), Rational::zero()];
        for patch in &patches {
            dec.check_deadline()?;
            let error = patch.chord()?;
            if error > worst {
                worst = error;
                let [a, c, d, _] = patch.corners()?;
                worst_shape = [l1(&a, &c), l1(&a, &d)];
            }
        }
        if worst <= *tolerance {
            break (patches, worst);
        }
        let cells = (spans << cuts_u) << cuts_v;
        if cuts_u >= 10 || cuts_v >= 12 || cells.saturating_mul(4) > max_points {
            return Err(Refusal::new(
                "mesh.resource-limit",
                "sheet grid cannot reach the chord tolerance within its limits",
            ));
        }
        // Refine across the longer side of the worst cell.
        if worst_shape[0] >= worst_shape[1] {
            patches = along(patches, 0);
            cuts_u += 1;
        } else {
            patches = along(patches, 1);
            cuts_v += 1;
        }
    };
    let per_span = 1usize << cuts_v;
    let rows = per_span + 1;
    let columns = patches.len() / per_span + 1;
    let mut points = Vec::with_capacity(columns * rows);
    let mut ordered: Vec<&Patch> = patches.iter().collect();
    ordered.sort_by(|a, c| {
        a.cell[0][0]
            .cmp(&c.cell[0][0])
            .then_with(|| a.cell[1][0].cmp(&c.cell[1][0]))
    });
    for column in 0..columns {
        let source = column.min(columns - 2);
        for row in 0..rows {
            let patch = ordered[source * per_span + row.min(per_span - 1)];
            let [a, c, d, e] = patch.corners()?;
            let corner = match (column > source, row >= per_span) {
                (false, false) => a,
                (true, false) => c,
                (false, true) => d,
                (true, true) => e,
            };
            points.push(corner);
        }
    }
    Ok(EnvelopeGrid {
        columns,
        rows,
        points,
        chord,
    })
}

/// Encloses every chart preimage of a box in the patches' frame by pruned
/// subdivision. Existence still needs a separate incidence proof.
pub(crate) fn cover(
    mut patches: Vec<Patch>,
    point: &[Ball; 3],
    dec: &Decider,
) -> Result<Option<[[Rational; 2]; 2]>, Refusal> {
    for _ in 0..24 {
        dec.check_deadline()?;
        let mut next = Vec::with_capacity(patches.len() * 2);
        for patch in patches {
            let hull = patch.hull()?;
            if (0..3).all(|i| meets(&hull[i], &point[i])) {
                next.extend(patch.split());
            }
        }
        if next.is_empty() {
            return Ok(None);
        }
        if next.len() > 4096 {
            return Err(unresolved("sheet preimage cover does not localise"));
        }
        patches = next;
    }
    let mut hull: Option<[[Rational; 2]; 2]> = None;
    for patch in patches {
        let cell = patch.cell;
        hull = Some(match hull {
            None => cell,
            Some(h) => [
                [
                    h[0][0].clone().min(cell[0][0].clone()),
                    h[0][1].clone().max(cell[0][1].clone()),
                ],
                [
                    h[1][0].clone().min(cell[1][0].clone()),
                    h[1][1].clone().max(cell[1][1].clone()),
                ],
            ],
        });
    }
    Ok(hull)
}

pub(crate) type Rect = [[Rational; 2]; 2];

fn unresolved_line(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

fn meets_rect(a: &Rect, c: &Rect) -> bool {
    (0..2).all(|axis| a[axis][0] <= c[axis][1] && c[axis][0] <= a[axis][1])
}
pub(crate) fn common(a: &Rect, c: &Rect) -> Option<Rect> {
    meets_rect(a, c).then(|| {
        std::array::from_fn(|axis| {
            [
                a[axis][0].clone().max(c[axis][0].clone()),
                a[axis][1].clone().min(c[axis][1].clone()),
            ]
        })
    })
}
fn join_rect(a: &Rect, c: &Rect) -> Rect {
    std::array::from_fn(|axis| {
        [
            a[axis][0].clone().min(c[axis][0].clone()),
            a[axis][1].clone().max(c[axis][1].clone()),
        ]
    })
}
fn width(rect: &Rect) -> Rational {
    (&rect[0][1] - &rect[0][0]).max(&rect[1][1] - &rect[1][0])
}

/// What the line solver needs from a sheet beyond its Bezier patches.
pub(crate) trait LineSheet: Sized {
    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal>;
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal>;
    /// The line's origin and direction in the frame of the patches.
    fn local_line(&self, origin: &Point3, direction: &Vector3) -> [[Real; 3]; 2];
    /// Patches in the frame, one per span along `u`, enclosed to `bits`.
    fn spans(&self, bits: u32) -> Vec<Patch>;
    fn chart(&self) -> Rect;
    /// Span ends along `u`.
    fn knots(&self) -> Vec<Rational>;
    /// Candidate line parameters on the row `v = 1` when `high`, else `v = 0`.
    /// [`LineSheet::on_boundary`] then keeps those on the row.
    fn row_contacts(
        &self,
        high: bool,
        line: &[[Real; 3]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal>;
    fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal>;
    /// The three homogeneous controls of the section at knot `index`, in the
    /// frame, from low row to high row.
    fn knot_section(&self, index: usize) -> [[Real; 4]; 3];
    /// The exact equations of one span over `rect`, with the first variable the
    /// line parameter, or `None` when they have no rational coefficients.
    fn span_system(
        &self,
        query: &Query<'_, Self>,
        index: usize,
        rect: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal>;
}

// What a contact may lie on: the knot section of a node, or a long boundary row.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Section {
    Knot(usize),
    High,
    Low,
}

// A rectangle owning at most one hit: `cells` holds any hit, and `domain`, which
// holds `cells` with a margin, is where the mean-value image `image` was proved.
struct Piece {
    cells: Rect,
    domain: Rect,
    image: Rect,
}

struct Contact {
    parameter: Real,
    local: [Ball; 3],
}

pub(crate) struct Query<'a, S: LineSheet> {
    sheet: &'a S,
    pub(crate) bits: u32,
    pub(crate) origin: &'a Point3,
    pub(crate) direction: &'a Vector3,
    // The same line in the definition frame, exactly and enclosed.
    pub(crate) exact: [[Real; 3]; 2],
    pub(crate) enclosed: [[Ball; 3]; 2],
    // Two independent exact normals of the line, enclosed, with their offsets.
    pub(crate) normals: [[Real; 3]; 2],
    enclosed_normals: [[Ball; 3]; 2],
    offsets: [Ball; 2],
    pub(crate) spans: Vec<Patch>,
}

/// Every point where a whole line meets the closed finite sheet, in line order.
/// Each is an exact transverse crossing and appears once, including at a knot
/// and on a boundary row. A tangency, or a contact this cannot separate, refuses.
/// It never reports a hit or a miss from an enclosure.
pub(crate) fn line_hits<S: LineSheet>(
    sheet: &S,
    origin: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Vec<Point3>, Refusal> {
    dec.check_deadline()?;
    math::validate_point(origin, dec)?;
    math::validate_vector(direction, dec)?;
    let query = Query::new(sheet, origin, direction, dec)?;
    let pieces = query.pieces(dec)?;
    let mut hits: Vec<Real> = Vec::with_capacity(pieces.len());
    for index in 0..pieces.len() {
        if let Some(parameter) = query.decide(&pieces, index, dec)? {
            let at = hits
                .iter()
                .try_fold(0, |at, hit| -> Result<usize, Refusal> {
                    Ok(if math::sign(&(&parameter - hit), dec)? == Sign::Positive {
                        at + 1
                    } else {
                        at
                    })
                })?;
            hits.insert(at, parameter);
        }
    }
    Ok(hits.iter().map(|hit| query.world(hit)).collect())
}

impl<'a, S: LineSheet> Query<'a, S> {
    fn new(
        sheet: &'a S,
        origin: &'a Point3,
        direction: &'a Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let bits = sheet.query_bits(dec)?;
        let frame = sheet.budgeted_frame(bits, dec)?;
        let exact = sheet.local_line(origin, direction);
        let along = math::budgeted_vector(direction, bits, dec)?;
        let enclosed = [
            frame.local(&math::budgeted_point(origin, bits, dec)?),
            frame.axes.each_ref().map(|axis| dot(&along, axis)),
        ];
        // Crossing with the axis the direction leans on least keeps both normals long.
        let least = (0..3)
            .min_by_key(|&i| magnitude(&enclosed[1][i]))
            .expect("three axes");
        let unit: [Real; 3] = std::array::from_fn(|i| {
            if i == least {
                Real::one()
            } else {
                Real::zero()
            }
        });
        let first = cross(&exact[1], &unit);
        let second = cross(&exact[1], &first);
        let normals = [first, second];
        let mut enclosed_normals = Vec::with_capacity(2);
        let mut offsets = Vec::with_capacity(2);
        for normal in &normals {
            let normal: [Ball; 3] = [
                dec.enclose_budgeted_at(&normal[0], bits)?,
                dec.enclose_budgeted_at(&normal[1], bits)?,
                dec.enclose_budgeted_at(&normal[2], bits)?,
            ];
            offsets.push(dot(&normal, &enclosed[0]));
            enclosed_normals.push(normal);
        }
        Ok(Self {
            sheet,
            bits,
            origin,
            direction,
            exact,
            enclosed,
            normals,
            enclosed_normals: enclosed_normals.try_into().expect("two normals"),
            offsets: offsets.try_into().expect("two normals"),
            spans: sheet.spans(bits),
        })
    }

    pub(crate) fn world(&self, parameter: &Real) -> Point3 {
        self.origin.plus(&self.direction.scaled(parameter))
    }

    // G on the controls of a patch. Positive weights keep each hit a zero of G.
    fn equations(&self, patch: &Patch) -> Vec<[Ball; 2]> {
        patch
            .net
            .iter()
            .map(|h| {
                std::array::from_fn(|k| {
                    let along = dot(
                        &[h[0].clone(), h[1].clone(), h[2].clone()],
                        &self.enclosed_normals[k],
                    );
                    &along - &(&self.offsets[k] * &h[3])
                })
            })
            .collect()
    }

    // Rectangles that together hold every hit, each with its mean-value image.
    fn pieces(&self, dec: &mut Decider) -> Result<Vec<Piece>, Refusal> {
        let mut cells = self.spans.clone();
        for depth in 0..MAX_DEPTH {
            let mut kept = Vec::with_capacity(cells.len());
            for cell in cells {
                dec.check_deadline()?;
                let values = self.equations(&cell);
                let apart = (0..2).any(|k| {
                    values.iter().all(|value| value[k].lo() > Rational::zero())
                        || values.iter().all(|value| value[k].hi() < Rational::zero())
                });
                if !apart {
                    kept.push(cell);
                }
            }
            cells = kept;
            if cells.is_empty() {
                return Ok(Vec::new());
            }
            if cells.len() > MAX_CELLS {
                break;
            }
            if depth >= 2
                && let Some(pieces) = self.resolved(&cells, dec)?
            {
                return Ok(pieces);
            }
            cells = cells.into_iter().flat_map(Patch::split).collect();
        }
        Err(unresolved_line(
            "a line contact on this finite spline is not a separated transverse crossing",
        ))
    }

    fn resolved(&self, cells: &[Patch], dec: &mut Decider) -> Result<Option<Vec<Piece>>, Refusal> {
        let chart = self.sheet.chart();
        // Each cell's margin is its own size, so a hit sits inside its rectangle.
        let mut parts: Vec<(Rect, Rect)> = cells
            .iter()
            .map(|cell| {
                let wider = std::array::from_fn(|axis| {
                    let [lo, hi] = &cell.cell[axis];
                    let margin = hi - lo;
                    [
                        (lo - &margin).max(chart[axis][0].clone()),
                        (hi + &margin).min(chart[axis][1].clone()),
                    ]
                });
                (cell.cell.clone(), wider)
            })
            .collect();
        let mut index = 0;
        while index < parts.len() {
            let mut other = index + 1;
            let mut merged = false;
            while other < parts.len() {
                if meets_rect(&parts[index].1, &parts[other].1) {
                    let (cells, domain) = parts.swap_remove(other);
                    parts[index] = (
                        join_rect(&parts[index].0, &cells),
                        join_rect(&parts[index].1, &domain),
                    );
                    merged = true;
                } else {
                    other += 1;
                }
            }
            if !merged {
                index += 1;
            }
        }
        let mut pieces = Vec::with_capacity(parts.len());
        for (cells, domain) in parts {
            dec.check_deadline()?;
            let Some(image) = self.image(&domain) else {
                return Ok(None);
            };
            pieces.push(Piece {
                cells,
                domain,
                image,
            });
        }
        Ok(Some(pieces))
    }

    // m - Y G(m) + (I - Y J)(R - m), or `None` without a contracting Y on R.
    // J's hull covers every span of R, so the bound holds across knots.
    fn image(&self, rect: &Rect) -> Option<Rect> {
        let bits = self.bits;
        let degree = self.spans[0].degree;
        let mut hull: [[Option<Ball>; 2]; 2] = [[None, None], [None, None]];
        let mut centre = None;
        let middle: [Rational; 2] =
            std::array::from_fn(|axis| (&rect[axis][0] + &rect[axis][1]).mul_pow2(-1));
        for span in &self.spans {
            let Some(part) = common(&span.cell, rect) else {
                continue;
            };
            if span.cell[0][0] <= middle[0] && middle[0] <= span.cell[0][1] && centre.is_none() {
                let point =
                    std::array::from_fn(|axis| [middle[axis].clone(), middle[axis].clone()]);
                centre = Some(self.equations(&span.clip(&point, bits))[0].clone());
            }
            if part[0][0] == part[0][1] || part[1][0] == part[1][1] {
                continue;
            }
            let values = self.equations(&span.clip(&part, bits));
            let steps = [
                b(&(q(degree as i64) / (&part[0][1] - &part[0][0])), bits),
                b(&(q(2) / (&part[1][1] - &part[1][0])), bits),
            ];
            let at = |i: usize, j: usize| &values[i * 3 + j];
            for i in 0..=degree {
                for j in 0..3 {
                    for (k, row) in hull.iter_mut().enumerate() {
                        let mut include = |slope: Ball, axis: usize| {
                            let entry = &mut row[axis];
                            *entry = Some(match entry.take() {
                                Some(wide) => wide.hull(&slope),
                                None => slope,
                            });
                        };
                        if i < degree {
                            include(&(&at(i + 1, j)[k] - &at(i, j)[k]) * &steps[0], 0);
                        }
                        if j < 2 {
                            include(&(&at(i, j + 1)[k] - &at(i, j)[k]) * &steps[1], 1);
                        }
                    }
                }
            }
        }
        let [[Some(a), Some(c)], [Some(d), Some(e)]] = hull else {
            return None;
        };
        let jacobian = [[a, c], [d, e]];
        let value: [Ball; 2] = centre?;
        let mid = [
            [jacobian[0][0].mid(), jacobian[0][1].mid()],
            [jacobian[1][0].mid(), jacobian[1][1].mid()],
        ];
        let determinant = &mid[0][0] * &mid[1][1] - &mid[0][1] * &mid[1][0];
        let inverse = determinant.recip()?;
        let y = [
            [
                b(&(&mid[1][1] * &inverse), bits),
                b(&-(&mid[0][1] * &inverse), bits),
            ],
            [
                b(&-(&mid[1][0] * &inverse), bits),
                b(&(&mid[0][0] * &inverse), bits),
            ],
        ];
        let spread: [Ball; 2] = std::array::from_fn(|axis| {
            b(&rect[axis][0], bits).hull(&b(&rect[axis][1], bits)) - b(&middle[axis], bits)
        });
        let mut image = rect.clone();
        for row in 0..2 {
            let mut total = Rational::zero();
            let mut value_at =
                b(&middle[row], bits) - (&(&y[row][0] * &value[0]) + &(&y[row][1] * &value[1]));
            for column in 0..2 {
                let identity = if row == column { 1 } else { 0 };
                let residual = b(&q(identity), bits)
                    - (&(&y[row][0] * &jacobian[0][column]) + &(&y[row][1] * &jacobian[1][column]));
                total = total + magnitude(&residual);
                value_at = value_at + &residual * &spread[column];
            }
            if total >= Rational::one() {
                return None;
            }
            image[row] = [value_at.lo(), value_at.hi()];
        }
        Some(image)
    }

    fn touched(&self, image: &Rect) -> Vec<Section> {
        let mut touched: Vec<Section> = self
            .sheet
            .knots()
            .iter()
            .enumerate()
            .filter(|(_, knot)| image[0][0] <= **knot && **knot <= image[0][1])
            .map(|(index, _)| Section::Knot(index))
            .collect();
        if image[1][0] <= Rational::zero() {
            touched.push(Section::Low);
        }
        if image[1][1] >= Rational::one() {
            touched.push(Section::High);
        }
        touched
    }

    // The line parameter of the hit of one piece, if it has one.
    fn decide(
        &self,
        pieces: &[Piece],
        index: usize,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let piece = &pieces[index];
        let mut domain = piece.domain.clone();
        let mut image = piece.image.clone();
        let mut tried = Vec::new();
        for _ in 0..MAX_NARROWING {
            dec.check_deadline()?;
            let Some(core) = common(&image, &piece.cells) else {
                return Ok(None);
            };
            let touched = self.touched(&core);
            if touched.is_empty()
                && width(&core) < Rational::one().mul_pow2(-i64::from(self.bits / 4))
            {
                return self.interior(piece, &core, dec);
            }
            for section in touched {
                if tried.contains(&section) {
                    continue;
                }
                tried.push(section);
                if let Some(contact) = self.owned_contact(section, pieces, index, dec)? {
                    return Ok(Some(contact));
                }
            }
            let Some(next) = common(&image, &domain) else {
                return Ok(None);
            };
            let Some(narrower) = self.image(&next) else {
                break;
            };
            if next == domain && narrower == image {
                break;
            }
            domain = next;
            image = narrower;
        }
        Err(unresolved_line(
            "a line contact on this finite spline is not separated from a chart knot or boundary",
        ))
    }

    // The contact on a section that belongs to this piece. Every hit lies in some
    // piece, so a contact clear of all the others is the piece's unique hit.
    fn owned_contact(
        &self,
        section: Section,
        pieces: &[Piece],
        index: usize,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        for parameter in self.contacts(section, dec)? {
            let contact = self.contact(parameter, dec)?;
            let mut owners = Vec::new();
            for (other, piece) in pieces.iter().enumerate() {
                if self.may_hold(&piece.domain, &contact)? {
                    owners.push(other);
                }
            }
            match owners.as_slice() {
                [owner] if *owner == index => return Ok(Some(contact.parameter)),
                [_] => {}
                [] => {
                    return Err(unresolved_line(
                        "an exact line contact has no chart cell on this finite spline",
                    ));
                }
                _ if owners.contains(&index) => {
                    return Err(unresolved_line(
                        "line contacts on this finite spline are not separated",
                    ));
                }
                _ => {}
            }
        }
        Ok(None)
    }

    fn contact(&self, parameter: Real, dec: &mut Decider) -> Result<Contact, Refusal> {
        let along = dec.enclose_budgeted_at(&parameter, self.bits)?;
        Ok(Contact {
            parameter,
            local: std::array::from_fn(|i| &self.enclosed[0][i] + &(&along * &self.enclosed[1][i])),
        })
    }

    fn may_hold(&self, rect: &Rect, contact: &Contact) -> Result<bool, Refusal> {
        for span in &self.spans {
            let Some(part) = common(&span.cell, rect) else {
                continue;
            };
            let hull = span.clip(&part, self.bits).hull()?;
            if (0..3).all(|i| meets(&hull[i], &contact.local[i])) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    // Exact line parameters where the line meets one section of the chart.
    fn contacts(&self, section: Section, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        match section {
            Section::Knot(index) => self.knot_contacts(index, dec),
            Section::Low | Section::High => {
                let candidates =
                    self.sheet
                        .row_contacts(section == Section::High, &self.exact, dec)?;
                self.on_boundary(candidates, dec)
            }
        }
    }

    fn on_boundary(&self, candidates: Vec<Real>, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let mut kept = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            if self.sheet.on_boundary(&self.world(&candidate), dec)? {
                kept.push(candidate);
            }
        }
        Ok(kept)
    }

    // The section at a knot is one rational quadratic in v with exact controls,
    // so its contacts are the common roots of the two equations of the line.
    fn knot_contacts(&self, index: usize, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let section = self.sheet.knot_section(index);
        let [origin, direction] = &self.exact;
        let along = |normal: &[Real; 3], row: &[Real; 4]| -> Real {
            let offset = dot_real(normal, origin);
            &dot_real(normal, &[row[0].clone(), row[1].clone(), row[2].clone()])
                - &(&offset * &row[3])
        };
        let equation = |normal: &[Real; 3]| {
            let [c0, c1, c2] = [
                along(normal, &section[0]),
                along(normal, &section[1]),
                along(normal, &section[2]),
            ];
            let two = Real::from_i64(2);
            let linear = &two * &(&c1 - &c0);
            let quadratic = &(&c0 - &(&two * &c1)) + &c2;
            RealPoly::from_terms(
                1,
                vec![(vec![0], c0), (vec![1], linear), (vec![2], quadratic)],
            )
        };
        let shared = equation(&self.normals[0]).gcd(&equation(&self.normals[1]), dec)?;
        if shared.is_zero(dec)? {
            return Err(unresolved_line(
                "a line along a knot section needs an overlap solver",
            ));
        }
        let mut parameters = Vec::new();
        for v in shared.real_roots(dec)? {
            if math::sign(&v, dec)? == Sign::Negative
                || math::sign(&(&v - &Real::one()), dec)? == Sign::Positive
            {
                continue;
            }
            let rest = &Real::one() - &v;
            let basis = [
                rest.square(),
                &Real::from_i64(2) * &(&rest * &v),
                v.square(),
            ];
            let weighted = |k: usize| {
                (0..3).fold(Real::zero(), |sum, j| {
                    plus(&sum, &times(&section[j][k], &basis[j]))
                })
            };
            let weight = weighted(3);
            let relative: [Real; 3] =
                std::array::from_fn(|k| &(&weighted(k) / &weight) - &origin[k]);
            parameters.push(&dot_real(&relative, direction) / &dot_real(direction, direction));
        }
        Ok(parameters)
    }

    // A hit clear of every knot and row: one root of the span's exact equations.
    fn interior(
        &self,
        piece: &Piece,
        core: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let knots = self.sheet.knots();
        let index = knots
            .windows(2)
            .position(|pair| pair[0] <= core[0][0] && core[0][1] <= pair[1])
            .expect("a core clear of every knot lies in one span");
        let span = [&knots[index], &knots[index + 1]];
        // Room for the certificate, stopping short of the rows, knots and domain.
        let limits = [
            [
                span[0].clone().max(piece.domain[0][0].clone()),
                span[1].clone().min(piece.domain[0][1].clone()),
            ],
            [
                Rational::zero().max(piece.domain[1][0].clone()),
                Rational::one().min(piece.domain[1][1].clone()),
            ],
        ];
        let room: [Rational; 2] = std::array::from_fn(|axis| {
            let reach = (&core[axis][0] - &limits[axis][0]).min(&limits[axis][1] - &core[axis][1]);
            (reach / q(2))
                .min((&core[axis][1] - &core[axis][0]) * q(4) + Rational::one().mul_pow2(-24))
        });
        if room.iter().any(|room| room.sign() != Sign::Positive) {
            return Err(unresolved_line(
                "a line contact leaves no certificate room on this finite spline",
            ));
        }
        let box_: Rect = std::array::from_fn(|axis| {
            [&core[axis][0] - &room[axis], &core[axis][1] + &room[axis]]
        });
        let Some((system, candidate)) = self.sheet.span_system(self, index, &box_, dec)? else {
            return Err(unresolved_line(
                "this finite spline crossing has no rational equations",
            ));
        };
        match dec.certify_root_at(&system, &candidate, self.bits) {
            Outcome::Proven(BoxRoots::One(root)) => Ok(Some(root.coordinate(0))),
            Outcome::Proven(BoxRoots::None) => Ok(None),
            Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
            Outcome::ProvenCoincident(_) => {
                Err(unresolved_line("a root certificate reported coincidence"))
            }
            Outcome::Refused(refusal) => Err(refusal),
        }
    }
}

pub(crate) fn cross(a: &[Real; 3], c: &[Real; 3]) -> [Real; 3] {
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &times(&a[j], &c[k]) - &times(&a[k], &c[j])
    })
}
pub(crate) fn dot_real(a: &[Real; 3], c: &[Real; 3]) -> Real {
    (0..3).fold(Real::zero(), |sum, i| plus(&sum, &times(&a[i], &c[i])))
}
pub(crate) fn rational_dot(a: &[Rational; 3], c: &[Rational; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + &a[i] * &c[i])
}
