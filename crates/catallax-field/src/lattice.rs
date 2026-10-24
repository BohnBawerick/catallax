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
    lo: [Rational; 3],
    hi: [Rational; 3],
    fast_start: [f64; 3],
    fast_direction: [f64; 3],
    fast_length_sq: f64,
    fast_lo: [f64; 3],
    fast_hi: [f64; 3],
}

#[derive(Clone, Debug)]
enum GradeDef {
    Constant(Rational),
    Linear {
        from: [Rational; 3],
        direction: [Rational; 3],
        length_sq: Rational,
        start: Rational,
        end: Rational,
    },
}

/// Whether a box is proven inside, proven outside, or unresolved.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum BoxState {
    /// Strictly inside on the closed box.
    Inside,
    /// Strictly outside on the closed box.
    Outside,
    /// Not settled by interval arithmetic at this box size.
    Unknown,
}

/// A compiled lattice leaf.
#[derive(Clone, Debug)]
pub(crate) struct Lattice {
    rotation: [[Rational; 3]; 3],
    inverse_rotation: [[Rational; 3]; 3],
    translation: [Rational; 3],
    fast_rotation: [[f64; 3]; 3],
    fast_inverse: [[f64; 3]; 3],
    fast_translation: [f64; 3],
    size: Rational,
    fast_size: f64,
    beams: Vec<Beam>,
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
fn point(p: &[Rational; 3]) -> RPoint3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone())
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
                Beam {
                    fast_start: start.each_ref().map(f),
                    fast_direction: direction.each_ref().map(f),
                    fast_length_sq: f(&rational_dot(&direction, &direction)),
                    fast_lo: lo.each_ref().map(|r| f(r).next_down()),
                    fast_hi: hi.each_ref().map(|r| f(r).next_up()),
                    length_sq: rational_dot(&direction, &direction),
                    start,
                    direction,
                    lo,
                    hi,
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
        let (grade, lipschitz) = match grade {
            Grade::Constant(r) => (GradeDef::Constant(r.value().clone()), q(1)),
            Grade::Linear {
                from,
                to,
                start,
                end,
            } => {
                let from = from.0.each_ref().map(|n| n.value().clone());
                let to = to.0.each_ref().map(|n| n.value().clone());
                let direction: [Rational; 3] = std::array::from_fn(|i| &to[i] - &from[i]);
                let length_sq = rational_dot(&direction, &direction);
                let length_lo = Ball::from_rational(&length_sq, PRECISION)
                    .sqrt()
                    .expect("positive squared length")
                    .lo();
                let slope = (end.value() - start.value()).abs();
                let lipschitz = &q(1)
                    + &slope.checked_div(&length_lo).ok_or_else(|| {
                        Refusal::new(
                            "lattice.impossible",
                            "a linear grade needs two distinct points",
                        )
                    })?;
                (
                    GradeDef::Linear {
                        from,
                        direction,
                        length_sq,
                        start: start.value().clone(),
                        end: end.value().clone(),
                    },
                    lipschitz,
                )
            }
            Grade::Distance { .. } => unreachable!("validate refuses distance grades"),
        };
        Ok(Self {
            fast_rotation: fast_matrix(&rotation),
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

    /// Smallest and largest radius of the grade.
    pub(crate) fn radius_range(&self) -> &[Rational; 2] {
        &self.radius_range
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
        match &self.grade {
            GradeDef::Constant(r) => r.clone(),
            GradeDef::Linear {
                from,
                direction,
                length_sq,
                start,
                end,
            } => {
                let rel: [Rational; 3] = std::array::from_fn(|i| &p[i] - &from[i]);
                let t = rational_dot(&rel, direction)
                    .checked_div(length_sq)
                    .expect("distinct grade points");
                let t = t.max(q(0)).min(q(1));
                start + &(&t * &(end - start))
            }
        }
    }
    fn fast_radius(&self, p: [f64; 3]) -> f64 {
        match &self.grade {
            GradeDef::Constant(r) => f(r),
            GradeDef::Linear {
                from,
                direction,
                length_sq,
                start,
                end,
            } => {
                let rel: [f64; 3] = std::array::from_fn(|i| p[i] - f(&from[i]));
                let t = dot(&rel, &direction.each_ref().map(f)) / f(length_sq);
                let t = t.max(0.0).min(1.0);
                f(start) + t * (f(end) - f(start))
            }
        }
    }
    fn interval_radius(&self, p: &[Interval; 3]) -> Interval {
        match &self.grade {
            GradeDef::Constant(r) => Interval::from_rational(r),
            GradeDef::Linear {
                from,
                direction,
                length_sq,
                start,
                end,
            } => {
                let rel: [Interval; 3] =
                    std::array::from_fn(|i| p[i] - Interval::from_rational(&from[i]));
                let t = (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
                    s + rel[i] * Interval::from_rational(&direction[i])
                }) / Interval::from_rational(length_sq);
                let t = clamp01(t);
                Interval::from_rational(start)
                    + t * (Interval::from_rational(end) - Interval::from_rational(start))
            }
        }
    }

    /// Exact squared distance to the nearest beam, with every nearest beam's
    /// foot vector `local - foot` and clamp parameter.
    fn nearest(
        &self,
        local: &[Rational; 3],
        collect_feet: bool,
    ) -> Option<(Rational, Vec<([Rational; 3], Rational)>)> {
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
                                if collect_feet && !feet.contains(&(foot.clone(), t.clone())) {
                                    feet.push((foot, t));
                                }
                            }
                            _ => {
                                best = Some(d);
                                feet.clear();
                                if collect_feet {
                                    feet.push((foot, t));
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
                        let t = (dot(&w, &beam.fast_direction) / beam.fast_length_sq)
                            .max(0.0)
                            .min(1.0);
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

    /// Outward (unnormalised) boundary direction at a boundary point, or
    /// `None` on a crease between equidistant beams or at a grade clamp.
    pub(crate) fn normal(&self, p: &[Rational; 3]) -> Option<Vector3> {
        let local = self.local(p);
        let (distance_sq, feet) = self.nearest(&local, true)?;
        if feet.len() != 1 || distance_sq.is_zero() {
            return None;
        }
        let foot = apply(&self.rotation, &feet[0].0);
        let GradeDef::Linear {
            from,
            direction,
            length_sq,
            start,
            end,
        } = &self.grade
        else {
            return Some(Vector3::new(
                Real::from(foot[0].clone()),
                Real::from(foot[1].clone()),
                Real::from(foot[2].clone()),
            ));
        };
        let rel: [Rational; 3] = std::array::from_fn(|i| &p[i] - &from[i]);
        let t = rational_dot(&rel, direction)
            .checked_div(length_sq)
            .expect("distinct grade points");
        if t <= q(0) || t >= q(1) {
            return None;
        }
        // ∇(β - ρ) = foot/|foot| - s·direction with s = (end - start)/|direction|².
        // The same direction is foot - |foot|·s·direction, which keeps one root.
        let scale = (end - start)
            .checked_div(length_sq)
            .expect("distinct grade points");
        let length = Real::from(distance_sq).sqrt();
        Some(Vector3::new(
            &Real::from(foot[0].clone()) - &(&length * &Real::from(&scale * &direction[0])),
            &Real::from(foot[1].clone()) - &(&length * &Real::from(&scale * &direction[1])),
            &Real::from(foot[2].clone()) - &(&length * &Real::from(&scale * &direction[2])),
        ))
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
        let world_feet: Vec<[Rational; 3]> = feet
            .iter()
            .map(|(foot, _)| apply(&self.rotation, foot))
            .collect();
        // Along d, β changes at the smallest rate (foot · d)/|foot| over the
        // nearest beams, and every nearest foot has length ρ(p) > 0. The sign of
        // the rate of β - ρ is therefore the sign of min(foot · d) - ρ·(dρ/ds).
        let mut rate: Option<Real> = None;
        for foot in &world_feet {
            let candidate = &(&Real::from(foot[0].clone()) * &d.x)
                + &(&(&Real::from(foot[1].clone()) * &d.y)
                    + &(&Real::from(foot[2].clone()) * &d.z));
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
        if let GradeDef::Linear {
            from,
            direction,
            length_sq,
            start,
            end,
        } = &self.grade
        {
            let rel: [Rational; 3] = std::array::from_fn(|i| &p[i] - &from[i]);
            let t = rational_dot(&rel, direction)
                .checked_div(length_sq)
                .expect("distinct grade points");
            let along = &(&Real::from(direction[0].clone()) * &d.x)
                + &(&(&Real::from(direction[1].clone()) * &d.y)
                    + &(&Real::from(direction[2].clone()) * &d.z));
            let along_sign = geometry::sign(&along, dec)?;
            let ramp = if t > q(0) && t < q(1) {
                true
            } else if t == q(0) {
                along_sign == Sign::Positive
            } else if t == q(1) {
                along_sign == Sign::Negative
            } else {
                false
            };
            if ramp {
                let slope = Real::from(
                    &(end - start)
                        .checked_div(length_sq)
                        .expect("distinct grade points")
                        * &self.radius(p),
                );
                rate = &rate - &(&slope * &along);
            }
        }
        Ok(match geometry::sign(&rate, dec)? {
            Sign::Negative => Location::Inside,
            Sign::Positive => Location::Outside,
            Sign::Zero => Location::Indeterminate,
        })
    }

    /// Composes a rational rigid motion in front of the frame.
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
        if map.reversing() {
            return Err(Refusal::new(
                "field.lattice-frame",
                "a mirrored lattice frame is not a rotation; unfinished coverage",
            ));
        }
        // world' = M (R l + t0) + t = (M R) l + (M t0 + t)
        let rotation: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).fold(q(0), |s, k| s + &m[i][k] * &self.rotation[k][j]))
        });
        let shifted = apply(&m, &self.translation);
        self.translation = std::array::from_fn(|i| &shifted[i] + &t[i]);
        self.inverse_rotation = transpose(&rotation);
        self.fast_rotation = fast_matrix(&rotation);
        self.fast_inverse = fast_matrix(&self.inverse_rotation);
        self.fast_translation = self.translation.each_ref().map(f);
        self.rotation = rotation;
        if let GradeDef::Linear {
            from, direction, ..
        } = &mut self.grade
        {
            let moved = apply(&m, from);
            *from = std::array::from_fn(|i| &moved[i] + &t[i]);
            *direction = apply(&m, direction);
        }
        Ok(())
    }

    /// Interval classification of a world box. Candidate beams come from the
    /// cells the box can reach once grown by the largest radius.
    pub(crate) fn box_state(&self, world: &[Interval; 3]) -> BoxState {
        if self.beams.is_empty() {
            return BoxState::Outside;
        }
        let shifted: [Interval; 3] =
            std::array::from_fn(|i| world[i] - Interval::from_rational(&self.translation[i]));
        let local = interval_apply(&self.inverse_rotation, &shifted);
        let radius = self.interval_radius(world);
        let radius_sq = radius.square();
        let r_hi = radius.hi();
        let a = self.fast_size;
        let margin = self.overhang_cells + 1;
        let range: [[i64; 2]; 3] = std::array::from_fn(|i| {
            [
                ((local[i].lo() - r_hi) / a).floor() as i64 - margin,
                ((local[i].hi() + r_hi) / a).floor() as i64 + margin,
            ]
        });
        let mut outside = true;
        for kx in range[0][0]..=range[0][1] {
            for ky in range[1][0]..=range[1][1] {
                for kz in range[2][0]..=range[2][1] {
                    let shift = [kx as f64 * a, ky as f64 * a, kz as f64 * a];
                    for beam in &self.beams {
                        let mut misses = false;
                        for i in 0..3 {
                            let lo = beam.fast_lo[i] + shift[i] - r_hi - a * 1e-9;
                            let hi = beam.fast_hi[i] + shift[i] + r_hi + a * 1e-9;
                            if local[i].hi() < lo || local[i].lo() > hi {
                                misses = true;
                                break;
                            }
                        }
                        if misses {
                            continue;
                        }
                        let w: [Interval; 3] = std::array::from_fn(|i| {
                            local[i]
                                - Interval::from_rational(&beam.start[i])
                                - Interval::from_rational(&(&self.size * &q([kx, ky, kz][i])))
                        });
                        let t = (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
                            s + w[i] * Interval::from_rational(&beam.direction[i])
                        }) / Interval::from_rational(&beam.length_sq);
                        let t = clamp01(t);
                        let foot: [Interval; 3] = std::array::from_fn(|i| {
                            w[i] - t * Interval::from_rational(&beam.direction[i])
                        });
                        let d_sq = foot[0].square() + foot[1].square() + foot[2].square();
                        if d_sq.hi() < radius_sq.lo() {
                            return BoxState::Inside;
                        }
                        if !(d_sq.lo() > radius_sq.hi()) {
                            outside = false;
                        }
                    }
                }
            }
        }
        if outside {
            BoxState::Outside
        } else {
            BoxState::Unknown
        }
    }

    /// Interval residual bound for the nearest-support attribution of display triangles.
    pub(crate) fn fast_residual(&self, p: [f64; 3]) -> f64 {
        self.fast_distance(p).abs() * self.fast_lipschitz
    }

    /// World point of a local rational point, for tests and witnesses.
    pub(crate) fn world(&self, local: &[Rational; 3]) -> RPoint3 {
        let rotated = apply(&self.rotation, local);
        point(&std::array::from_fn(|i| &rotated[i] + &self.translation[i]))
    }
}
