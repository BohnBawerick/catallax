//! The lattice leaf of RECIPE.md section 4.15: `{ p : β(p) <= ρ(p) }`, where
//! `β` is the distance to the nearest beam of the periodic beam set and `ρ` is
//! the radius grade read at the query point.
//!
//! The beam set is infinite. Evaluation reduces the query point into its unit
//! cell and searches the neighbouring cells, which is exact because every
//! preset's covering radius is below one cell and a custom cell's search span
//! is derived from its overhang. Nothing is enumerated per region, so the leaf
//! is correct everywhere; the plan intersects it with the region body.
//!
//! The leaf's scalar is `(β - ρ) / L`, with `L` a rational upper bound on the
//! Lipschitz constant of `β - ρ`. Its magnitude never exceeds the distance to
//! the lattice boundary, which is what the plan's residual bounds require. For
//! a constant radius `L = 1` and the value is the exact distance outside.

use crate::grade::GradeDef;
use crate::metric::{PRECISION, Scalar};
use crate::{Location, geometry};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_recipe::{Grade, LatticeCell};
use catallax_shapes::{Isometry, RPoint3, Rotation, Vector3};

#[derive(Clone, Debug)]
struct Beam {
    start: [Rational; 3],
    direction: [Rational; 3],
    length_sq: Rational,
    interval_start: [Interval; 3],
    interval_direction: [Interval; 3],
    interval_length_sq: Interval,
    interval_bounds: [Interval; 3],
    fast_start: [f64; 3],
    fast_direction: [f64; 3],
    fast_length_sq: f64,
    fast_lo: [f64; 3],
    fast_hi: [f64; 3],
}

/// Whether a box is proven inside, proven outside, or unresolved.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BoxState {
    /// Strictly inside on the closed box.
    Inside,
    /// Strictly outside on the closed box.
    Outside,
    /// Not settled by interval arithmetic at this box size.
    Unknown,
}

/// One grid line that is a nearest beam everywhere in a box, with the grade's piece
/// there. On that box the leaf is `β - ρ`, with `β` the distance to the line.
#[derive(Clone, Debug)]
pub(crate) struct LinePiece {
    /// The axis the line is parallel to.
    axis: usize,
    /// The line's coordinates on the other two axes; the entry on `axis` is unused.
    line: [Interval; 3],
    /// The grade's gradient on the box.
    slope: [Interval; 3],
}

/// A compiled lattice leaf.
#[derive(Clone, Debug)]
pub(crate) struct Lattice {
    rotation: [[Rational; 3]; 3],
    inverse_rotation: [[Rational; 3]; 3],
    translation: [Rational; 3],
    fast_inverse: [[f64; 3]; 3],
    fast_translation: [f64; 3],
    size: Rational,
    fast_size: f64,
    beams: Vec<Beam>,
    cubic: bool,
    span: i64,
    order: Vec<i64>,
    overhang_cells: i64,
    grade: GradeDef,
    radius_range: [Rational; 2],
    lipschitz: Rational,
    fast_lipschitz: f64,
}

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn dot<T: Scalar>(a: &[T; 3], b: &[T; 3]) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}
fn rational_dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    &(&a[0] * &b[0]) + &(&(&a[1] * &b[1]) + &(&a[2] * &b[2]))
}
fn f(r: &Rational) -> f64 {
    r.to_f64_nearest()
}
fn matrix(rotation: &Rotation) -> Result<[[Rational; 3]; 3], Refusal> {
    let m = rotation.rational_matrix().ok_or_else(|| {
        Refusal::new(
            "field.lattice-frame",
            "a lattice frame needs a rational rotation",
        )
    })?;
    Ok(std::array::from_fn(|i| {
        std::array::from_fn(|j| m[i][j].clone())
    }))
}
fn transpose(a: &[[Rational; 3]; 3]) -> [[Rational; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i].clone()))
}
fn apply(m: &[[Rational; 3]; 3], v: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| rational_dot(&m[i], v))
}
fn apply_fast(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (m[i][0] * v[0] + m[i][1] * v[1]) + m[i][2] * v[2])
}
fn fast_matrix(m: &[[Rational; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| f(&m[i][j])))
}
fn interval_apply(m: &[[Rational; 3]; 3], v: &[Interval; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| {
        (0..3).fold(Interval::point(0.0).unwrap(), |s, j| {
            s + Interval::from_rational(&m[i][j]) * v[j]
        })
    })
}
fn imin(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().min(b.lo()), a.hi().min(b.hi())).unwrap()
}
fn imax(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().max(b.lo()), a.hi().max(b.hi())).unwrap()
}
fn clamp01(t: Interval) -> Interval {
    imin(
        imax(t, Interval::point(0.0).unwrap()),
        Interval::point(1.0).unwrap(),
    )
}

fn periodic_distance(coordinate: Interval) -> Interval {
    let half = Interval::point(0.5).unwrap();
    if !coordinate.lo().is_finite()
        || !coordinate.hi().is_finite()
        || coordinate.hi() - coordinate.lo() >= 1.0
        || coordinate.lo().abs().max(coordinate.hi().abs()) > 4_503_599_627_370_496.0
    {
        return Interval::new(0.0, 0.5).unwrap();
    }
    let low = coordinate.lo().floor();
    let high = coordinate.hi().floor();
    if low as i64 == high as i64 {
        let fraction = coordinate - Interval::point(low).unwrap();
        let distance = imin(fraction, Interval::point(1.0).unwrap() - fraction);
        imax(Interval::point(0.0).unwrap(), imin(distance, half))
    } else {
        let delta = coordinate - Interval::point(high).unwrap();
        Interval::new(0.0, delta.lo().abs().max(delta.hi().abs()).min(0.5)).unwrap()
    }
}

impl Lattice {
    /// Compiles the leaf. Impossible parameters refuse with `lattice.*` codes.
    pub(crate) fn new(
        cell: &LatticeCell,
        size: &Rational,
        frame: &catallax_recipe::Pose,
        grade: &Grade,
    ) -> Result<Self, Refusal> {
        let template = catallax_features::lattice::validate(cell, size, grade)?;
        let radius_range = catallax_features::lattice::radius_range(grade)?;
        let [w, x, y, z] = frame.rotation.each_ref().map(|n| n.value().clone());
        let rotation = matrix(&Rotation::from_quaternion(w, x, y, z)?)?;
        let inverse_rotation = transpose(&rotation);
        let translation = frame.translation.0.each_ref().map(|n| n.value().clone());
        let beams = template
            .beams
            .iter()
            .map(|b| {
                let start = b.start.each_ref().map(|c| c * size);
                let end = b.end.each_ref().map(|c| c * size);
                let direction: [Rational; 3] = std::array::from_fn(|i| &end[i] - &start[i]);
                let lo: [Rational; 3] =
                    std::array::from_fn(|i| start[i].clone().min(end[i].clone()));
                let hi: [Rational; 3] =
                    std::array::from_fn(|i| start[i].clone().max(end[i].clone()));
                let length_sq = rational_dot(&direction, &direction);
                Beam {
                    fast_start: start.each_ref().map(f),
                    fast_direction: direction.each_ref().map(f),
                    fast_length_sq: f(&length_sq),
                    fast_lo: lo.each_ref().map(|r| f(r).next_down()),
                    fast_hi: hi.each_ref().map(|r| f(r).next_up()),
                    interval_start: start.each_ref().map(Interval::from_rational),
                    interval_direction: direction.each_ref().map(Interval::from_rational),
                    interval_length_sq: Interval::from_rational(&length_sq),
                    interval_bounds: std::array::from_fn(|i| {
                        Interval::new(
                            Interval::from_rational(&lo[i]).lo(),
                            Interval::from_rational(&hi[i]).hi(),
                        )
                        .unwrap()
                    }),
                    length_sq,
                    start,
                    direction,
                }
            })
            .collect();
        let overhang_cells = template.overhang.ceil().to_i64().ok_or_else(|| {
            Refusal::new("lattice.impossible", "custom cell overhang is too large")
        })?;
        // Presets have covering radius below one cell: the neighbouring 27 cells
        // always hold the nearest beam. A custom cell's beams lie within its cube
        // grown by the overhang, so every cell point is within sqrt(3)(1+2m) cells
        // of one of its own beams; cells farther than that plus the overhang
        // cannot hold a nearer one.
        let span = if template.preset {
            1
        } else {
            overhang_cells
                .checked_mul(5)
                .and_then(|n| n.checked_add(3))
                .ok_or_else(|| {
                    Refusal::new(
                        "lattice.impossible",
                        "custom cell overhang exceeds the search index range",
                    )
                })?
        };
        let grade = GradeDef::compile(grade, "lattice.impossible")?
            .expect("validate refuses distance grades");
        let lipschitz = &q(1) + &grade.slope_bound();
        Ok(Self {
            cubic: matches!(cell, LatticeCell::Cubic),
            fast_inverse: fast_matrix(&inverse_rotation),
            fast_translation: translation.each_ref().map(f),
            rotation,
            inverse_rotation,
            translation,
            fast_size: f(size),
            size: size.clone(),
            beams,
            order: (0..=span)
                .flat_map(|d| if d == 0 { vec![0] } else { vec![-d, d] })
                .collect(),
            span,
            overhang_cells,
            grade,
            radius_range,
            fast_lipschitz: f(&lipschitz),
            lipschitz,
        })
    }

    /// Whether the cell has no beams at all, so the lattice is the empty set.
    pub(crate) fn is_empty(&self) -> bool {
        self.beams.is_empty()
    }

    /// Whether the cell is cubic and its frame a signed axis permutation.
    fn axis_cubic(&self) -> bool {
        self.cubic
            && self.rotation.iter().all(|row| {
                row.iter().filter(|v| !v.is_zero()).count() == 1
                    && row.iter().any(|v| v.abs() == q(1))
            })
    }

    /// Parameters of a constant cubic grid aligned with the model axes.
    pub(crate) fn axis_cubic_grid(&self) -> Option<([Rational; 3], Rational, Rational)> {
        if !self.axis_cubic() {
            return None;
        }
        let GradeDef::Constant(radius) = &self.grade else {
            return None;
        };
        Some((self.translation.clone(), self.size.clone(), radius.clone()))
    }

    /// The smallest radius of a cubic grid aligned with the model axes, constant or
    /// graded.
    pub(crate) fn axis_cubic_radius(&self) -> Option<&Rational> {
        self.axis_cubic().then_some(&self.radius_range[0])
    }

    /// Whether the radius is a linear grade whose anchors differ along one part axis.
    pub(crate) fn axis_graded(&self) -> bool {
        self.grade.axis().is_some()
    }

    /// The axis and positions of an axis grade's anchor planes.
    pub(crate) fn grade_planes(&self) -> Option<(usize, [Rational; 2])> {
        self.grade.planes()
    }

    /// The grid line nearest every point of a closed rational box of an axis cubic
    /// grid, decided exactly, when it is the same line throughout, the box keeps away
    /// from it and the grade is one piece there. The cubic beams are the grid lines
    /// `x_i ∈ t_i + aℤ` parallel to each axis, whatever the signed axis permutation
    /// of the frame. The lines parallel to `axis` are nearest wherever the distance
    /// to the grid along `axis` is at least the distance along each other axis.
    pub(crate) fn line_piece(&self, bounds: &[[Rational; 2]; 3]) -> Option<LinePiece> {
        if !self.axis_cubic() {
            return None;
        }
        let piece = self.grade.exact_piece(
            &std::array::from_fn(|i| bounds[i][0].clone()),
            &std::array::from_fn(|i| bounds[i][1].clone()),
        )?;
        let half = self.size.mul_pow2(-1);
        let mut line: [Rational; 3] = std::array::from_fn(|_| q(0));
        let mut near: [[Rational; 2]; 3] = std::array::from_fn(|_| [q(0), q(0)]);
        for i in 0..3 {
            let index = (&(&(&bounds[i][0] - &self.translation[i]) / &self.size)
                + &Rational::from_ratio(1, 2).unwrap())
                .floor();
            let at = &self.translation[i] + &(&Rational::from_integer(index) * &self.size);
            if bounds[i][1] > &at + &half {
                return None;
            }
            let a = (&bounds[i][0] - &at).abs();
            let b = (&bounds[i][1] - &at).abs();
            let lower = if bounds[i][0] <= at && at <= bounds[i][1] {
                q(0)
            } else {
                a.clone().min(b.clone())
            };
            near[i] = [lower, a.max(b)];
            line[i] = at;
        }
        let axis = (0..3).max_by(|&a, &b| near[a][0].cmp(&near[b][0]))?;
        if (0..3).any(|g| g != axis && near[g][1] > near[axis][0])
            || (0..3).all(|i| i == axis || near[i][0].is_zero())
        {
            return None;
        }
        Some(LinePiece {
            axis,
            line: line.each_ref().map(Interval::from_rational),
            slope: self.grade.piece_gradient(piece),
        })
    }

    /// `β - ρ` over a box inside a line piece's box.
    pub(crate) fn line_value(&self, piece: &LinePiece, world: &[Interval; 3]) -> Option<Interval> {
        let [j, k] = [(piece.axis + 1) % 3, (piece.axis + 2) % 3];
        let beta =
            ((world[j] - piece.line[j]).square() + (world[k] - piece.line[k]).square()).sqrt()?;
        Some(beta - self.grade.interval(world))
    }

    /// `∇(β - ρ)` and its Hessian over a box inside a line piece's box. The radius is
    /// affine there, so the Hessian is that of `β`.
    pub(crate) fn line_derivatives(
        &self,
        piece: &LinePiece,
        world: &[Interval; 3],
    ) -> Option<([Interval; 3], [[Interval; 3]; 3])> {
        let [j, k] = [(piece.axis + 1) % 3, (piece.axis + 2) % 3];
        let xj = world[j] - piece.line[j];
        let xk = world[k] - piece.line[k];
        let beta = (xj.square() + xk.square()).sqrt()?;
        if beta.lo() <= 0.0 {
            return None;
        }
        let cube = beta * beta * beta;
        let zero = Interval::point(0.0).unwrap();
        let mut gradient = piece.slope.map(|v| zero - v);
        gradient[j] = gradient[j] + xj / beta;
        gradient[k] = gradient[k] + xk / beta;
        let mut hessian = [[zero; 3]; 3];
        hessian[j][j] = xk.square() / cube;
        hessian[k][k] = xj.square() / cube;
        hessian[j][k] = zero - xj * xk / cube;
        hessian[k][j] = hessian[j][k];
        Some((gradient, hessian))
    }

    fn local(&self, p: &[Rational; 3]) -> [Rational; 3] {
        let shifted: [Rational; 3] = std::array::from_fn(|i| &p[i] - &self.translation[i]);
        apply(&self.inverse_rotation, &shifted)
    }
    fn fast_local(&self, p: [f64; 3]) -> [f64; 3] {
        let shifted = std::array::from_fn(|i| p[i] - self.fast_translation[i]);
        apply_fast(&self.fast_inverse, &shifted)
    }
    fn radius(&self, p: &[Rational; 3]) -> Rational {
        self.grade.value(p)
    }
    fn fast_radius(&self, p: [f64; 3]) -> f64 {
        self.grade.fast(p)
    }
    fn interval_radius(&self, p: &[Interval; 3]) -> Interval {
        self.grade.interval(p)
    }

    /// Exact squared distance to the nearest beam, with every nearest beam's
    /// foot vector `local - foot`.
    fn nearest(
        &self,
        local: &[Rational; 3],
        collect_feet: bool,
    ) -> Option<(Rational, Vec<[Rational; 3]>)> {
        if self.cubic && !collect_feet {
            let distances: [Rational; 3] = std::array::from_fn(|i| {
                let cell = Rational::from_integer((&local[i] / &self.size).floor());
                let remainder = &local[i] - &(&cell * &self.size);
                let distance = remainder.clone().min(&self.size - &remainder);
                &distance * &distance
            });
            let squared = (&distances[0] + &distances[1])
                .min(&distances[0] + &distances[2])
                .min(&distances[1] + &distances[2]);
            return Some((squared, Vec::new()));
        }
        let base: [Rational; 3] =
            std::array::from_fn(|i| Rational::from_integer((&local[i] / &self.size).floor()));
        let mut best: Option<Rational> = None;
        let mut feet = Vec::new();
        for kx in -self.span..=self.span {
            for ky in -self.span..=self.span {
                for kz in -self.span..=self.span {
                    let shift: [Rational; 3] =
                        std::array::from_fn(|i| &self.size * &(&base[i] + &q([kx, ky, kz][i])));
                    for beam in &self.beams {
                        let w: [Rational; 3] =
                            std::array::from_fn(|i| &(&local[i] - &beam.start[i]) - &shift[i]);
                        let t = rational_dot(&w, &beam.direction)
                            .checked_div(&beam.length_sq)
                            .expect("validated beam")
                            .max(q(0))
                            .min(q(1));
                        let foot: [Rational; 3] =
                            std::array::from_fn(|i| &w[i] - &(&t * &beam.direction[i]));
                        let d = rational_dot(&foot, &foot);
                        match &best {
                            Some(b) if d > *b => continue,
                            Some(b) if d == *b => {
                                if collect_feet && !feet.contains(&foot) {
                                    feet.push(foot);
                                }
                            }
                            _ => {
                                best = Some(d);
                                feet.clear();
                                if collect_feet {
                                    feet.push(foot);
                                }
                            }
                        }
                    }
                }
            }
        }
        best.map(|b| (b, feet))
    }

    /// Exact classification of a rational world point.
    pub(crate) fn location(&self, p: &[Rational; 3]) -> Location {
        let Some((distance_sq, _)) = self.nearest(&self.local(p), false) else {
            return Location::Outside;
        };
        let radius = self.radius(p);
        let radius_sq = &radius * &radius;
        match distance_sq.cmp(&radius_sq) {
            std::cmp::Ordering::Less => Location::Inside,
            std::cmp::Ordering::Equal => Location::Boundary,
            std::cmp::Ordering::Greater => Location::Outside,
        }
    }

    /// Enclosure of `(β - ρ) / L`, negative inside.
    pub(crate) fn proven_distance(&self, p: &RPoint3, location: Location) -> Ball {
        if location == Location::Boundary {
            return Ball::from_rational(&Rational::zero(), PRECISION);
        }
        let coordinates = [p.x.clone(), p.y.clone(), p.z.clone()];
        let (distance_sq, _) = self
            .nearest(&self.local(&coordinates), false)
            .expect("a compiled lattice leaf has beams");
        let beta = Ball::from_rational(&distance_sq, PRECISION)
            .sqrt()
            .expect("non-negative squared distance");
        let value = &beta - &Ball::from_rational(&self.radius(&coordinates), PRECISION);
        value
            .checked_div(&Ball::from_rational(&self.lipschitz, PRECISION))
            .expect("positive Lipschitz bound")
    }

    /// Tentative `(β - ρ) / L` in binary64, with the software square root.
    pub(crate) fn fast_distance(&self, p: [f64; 3]) -> f64 {
        let local = self.fast_local(p);
        if self.cubic {
            let distances: [f64; 3] = std::array::from_fn(|i| {
                let remainder = local[i] - (local[i] / self.fast_size).floor() * self.fast_size;
                let distance = remainder.min(self.fast_size - remainder);
                distance * distance
            });
            let squared = (distances[0] + distances[1])
                .min(distances[0] + distances[2])
                .min(distances[1] + distances[2]);
            return (Scalar::sqrt(&squared) - self.fast_radius(p)) / self.fast_lipschitz;
        }
        let mut best = f64::INFINITY;
        let base: [i64; 3] = std::array::from_fn(|i| (local[i] / self.fast_size).floor() as i64);
        for &dx in &self.order {
            for &dy in &self.order {
                for &dz in &self.order {
                    let shift = [
                        (base[0] + dx) as f64 * self.fast_size,
                        (base[1] + dy) as f64 * self.fast_size,
                        (base[2] + dz) as f64 * self.fast_size,
                    ];
                    for beam in &self.beams {
                        let mut box_sq = 0.0;
                        for i in 0..3 {
                            let c = local[i] - shift[i];
                            let gap = (beam.fast_lo[i] - c).max(c - beam.fast_hi[i]).max(0.0);
                            box_sq += gap * gap;
                        }
                        if box_sq >= best {
                            continue;
                        }
                        let w: [f64; 3] =
                            std::array::from_fn(|i| local[i] - shift[i] - beam.fast_start[i]);
                        let t =
                            (dot(&w, &beam.fast_direction) / beam.fast_length_sq).clamp(0.0, 1.0);
                        let foot: [f64; 3] =
                            std::array::from_fn(|i| w[i] - t * beam.fast_direction[i]);
                        let d = dot(&foot, &foot);
                        if d < best {
                            best = d;
                        }
                    }
                }
            }
        }
        if !best.is_finite() {
            return f64::INFINITY;
        }
        (Scalar::sqrt(&best) - self.fast_radius(p)) / self.fast_lipschitz
    }

    /// Positive-infinitesimal ray classification from a rational point.
    pub(crate) fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let location = self.location(p);
        if location != Location::Boundary {
            return Ok(location);
        }
        let local = self.local(p);
        let (_, feet) = self
            .nearest(&local, true)
            .expect("boundary point has a nearest beam");
        // Along d, β changes at the smallest rate (foot · d)/|foot| over the
        // nearest beams, and every nearest foot has length ρ(p) > 0. The sign of
        // the rate of β - ρ is therefore the sign of min(foot · d) - ρ·(dρ/ds).
        let mut rate: Option<Real> = None;
        for foot in &feet {
            let [x, y, z] = apply(&self.rotation, foot);
            let candidate =
                &(&Real::from(x) * &d.x) + &(&(&Real::from(y) * &d.y) + &(&Real::from(z) * &d.z));
            rate = Some(match rate {
                None => candidate,
                Some(previous) => match dec.compare(&candidate, &previous) {
                    catallax_numbers::Outcome::Proven(catallax_numbers::StrictOrder::Less) => {
                        candidate
                    }
                    catallax_numbers::Outcome::Proven(_)
                    | catallax_numbers::Outcome::ProvenCoincident(_) => previous,
                    catallax_numbers::Outcome::NotExact(n) => {
                        return Err(Refusal::new(n.code, n.detail));
                    }
                    catallax_numbers::Outcome::Refused(r) => return Err(r),
                },
            });
        }
        let mut rate = rate.expect("at least one nearest beam");
        if let Some(change) = self.grade.rate(p, d, dec)? {
            rate = &rate - &(&change * &Real::from(self.radius(p)));
        }
        Ok(match geometry::sign(&rate, dec)? {
            Sign::Negative => Location::Inside,
            Sign::Positive => Location::Outside,
            Sign::Zero => Location::Indeterminate,
        })
    }

    /// Composes a rational rigid motion in front of the frame. A reflection is
    /// accepted: the frame matrix stays orthogonal, so its transpose is still the
    /// inverse and the local beam set is reused unchanged.
    pub(crate) fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let m = map.rational_matrix().ok_or_else(|| {
            Refusal::new(
                "field.lattice-frame",
                "lattice copies need a rational rigid motion",
            )
        })?;
        let m: [[Rational; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].clone()));
        let t = map.translation.as_rational().ok_or_else(|| {
            Refusal::new(
                "field.lattice-frame",
                "lattice copies need a rational rigid motion",
            )
        })?;
        let t = [t.x, t.y, t.z];
        // world' = M (R l + t0) + t = (M R) l + (M t0 + t)
        let rotation: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).fold(q(0), |s, k| s + &m[i][k] * &self.rotation[k][j]))
        });
        let shifted = apply(&m, &self.translation);
        self.translation = std::array::from_fn(|i| &shifted[i] + &t[i]);
        self.inverse_rotation = transpose(&rotation);
        self.fast_inverse = fast_matrix(&self.inverse_rotation);
        self.fast_translation = self.translation.each_ref().map(f);
        self.rotation = rotation;
        self.grade.transform(&m, &t);
        Ok(())
    }

    /// Interval classification of a world box. Candidate beams come from the
    /// cells the box can reach once grown by the largest radius.
    pub(crate) fn box_state(
        &self,
        world: &[Interval; 3],
        dec: &mut Decider,
    ) -> Result<BoxState, Refusal> {
        dec.check_deadline()?;
        if self.beams.is_empty() {
            return Ok(BoxState::Outside);
        }
        let cover = &self.size * &q(2 + 4 * self.overhang_cells);
        if self.radius_range[0] > cover {
            return Ok(BoxState::Inside);
        }
        let shifted: [Interval; 3] =
            std::array::from_fn(|i| world[i] - Interval::from_rational(&self.translation[i]));
        let local = interval_apply(&self.inverse_rotation, &shifted);
        let radius = self.interval_radius(world);
        let radius_sq = radius.square();
        let r_hi = radius.hi();
        let size = Interval::from_rational(&self.size);
        if size.lo() <= 0.0 || !r_hi.is_finite() {
            return Ok(BoxState::Unknown);
        }
        if self.cubic {
            // The cubic beam set is exactly the three families of grid lines.
            let distances = local.map(|coordinate| periodic_distance(coordinate / size) * size);
            let squared = imin(
                distances[0].square() + distances[1].square(),
                imin(
                    distances[0].square() + distances[2].square(),
                    distances[1].square() + distances[2].square(),
                ),
            );
            return Ok(if squared.hi() < radius_sq.lo() {
                BoxState::Inside
            } else if squared.lo() > radius_sq.hi() {
                BoxState::Outside
            } else {
                BoxState::Unknown
            });
        }
        let margin = self.overhang_cells as f64 + 1.0;
        let mut range = [[0i64; 2]; 3];
        for i in 0..3 {
            let cells = (local[i] + Interval::new(-r_hi, r_hi).unwrap()) / size;
            let lo = cells.lo().floor() - margin;
            let hi = cells.hi().ceil() + margin;
            // Conversion must not saturate an index and discard candidate beams.
            if !lo.is_finite() || !hi.is_finite() || lo <= i64::MIN as f64 || hi >= i64::MAX as f64
            {
                return Ok(BoxState::Unknown);
            }
            range[i] = [lo as i64, hi as i64];
        }
        let mut outside = true;
        for kx in range[0][0]..=range[0][1] {
            for ky in range[1][0]..=range[1][1] {
                for kz in range[2][0]..=range[2][1] {
                    let cell: [Interval; 3] = std::array::from_fn(|i| {
                        local[i] - size * Interval::from_rational(&q([kx, ky, kz][i]))
                    });
                    for beam in &self.beams {
                        dec.check_deadline()?;
                        let inflation = Interval::new(-r_hi, r_hi).unwrap();
                        if (0..3).any(|i| {
                            let beam_bound = beam.interval_bounds[i] + inflation;
                            cell[i].hi() < beam_bound.lo() || cell[i].lo() > beam_bound.hi()
                        }) {
                            continue;
                        }
                        let w: [Interval; 3] =
                            std::array::from_fn(|i| cell[i] - beam.interval_start[i]);
                        let t = (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
                            s + w[i] * beam.interval_direction[i]
                        }) / beam.interval_length_sq;
                        let t = clamp01(t);
                        let foot: [Interval; 3] =
                            std::array::from_fn(|i| w[i] - t * beam.interval_direction[i]);
                        let d_sq = foot[0].square() + foot[1].square() + foot[2].square();
                        if d_sq.hi() < radius_sq.lo() {
                            return Ok(BoxState::Inside);
                        }
                        if d_sq.lo() <= radius_sq.hi() {
                            outside = false;
                        }
                    }
                }
            }
        }
        Ok(if outside {
            BoxState::Outside
        } else {
            BoxState::Unknown
        })
    }
}
