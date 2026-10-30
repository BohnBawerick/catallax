//! The network gyroid leaf of RECIPE.md section 4.16: `{ p : g(p) <= c(p) }`, where
//! `g = sin X cos Y + sin Y cos Z + sin Z cos X` and `X = 2πx/a`.
//!
//! The level is a constant or a supported linear grade (see `crate::grade`): its
//! anchors differ along one part axis, both endpoints lie in `[-1/4, 1/4]` and its
//! slope `K` is at most `1/10` per millimetre. The frame is the part frame.
//! Sheets, other grades and other frames refuse explicitly.
//!
//! The leaf's scalar is `(g - c) / L` with `L = 22/a + K`. Each partial derivative
//! of `g` is a difference of two products of sines and cosines scaled by `2π/a`, so
//! `|∇g| <= 2√3·2π/a < 22/a`, using `π < 22/7` and `√3 < 7/4`, and the clamped
//! grade changes by at most `K` per millimetre. The scalar's magnitude therefore
//! never exceeds the distance to the level set. It is a dimensionless residual
//! made into a bound, not a signed distance.
//!
//! A graded level is affine between its anchor planes, where `g - c` has the
//! Hessian of `g`. A box across an anchor plane uses the range of `g` less the
//! exact clamped range of `c`, never a second-order bound across the kink.
//!
//! Binary64 and interval evaluation use fixed polynomials, not a platform math
//! library. Their remainder terms are part of the interval result.

use crate::grade::{GradeDef, Piece, check_certified};
use crate::metric::PRECISION;
use crate::{Location, geometry};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign, sin_cos_degrees};
use catallax_recipe::{Grade, GyroidKind, Pose};
use catallax_shapes::{Isometry, RPoint3, Vector3};

use crate::lattice::BoxState;

/// Taylor coefficients of `sin u / u` in `u²`, through `u^16 / 17!`.
const SINE_TERMS: usize = 9;
/// Taylor coefficients of `cos u` in `u²`, through `u^16 / 16!`.
const COSINE_TERMS: usize = 9;
/// Largest `|t|` whose quarter-turn reduction keeps absolute precision.
const LARGEST_TURN: f64 = 1_099_511_627_776.0;

fn point(v: f64) -> Interval {
    Interval::point(v).unwrap()
}
/// `1, -1/((n+1)(n+2)), ...`, the alternating reciprocal factorials from `n`.
fn alternating_series(first: i64) -> Vec<Rational> {
    let mut value = Rational::one();
    let mut n = first;
    let mut series = Vec::with_capacity(SINE_TERMS);
    for term in 0..SINE_TERMS {
        series.push(value.clone());
        if term + 1 < SINE_TERMS {
            value = &(-&value) / &Rational::from_i64((n + 1) * (n + 2));
            n += 2;
        }
    }
    series
}

/// One value and its derivatives at a point, for a box's second-order bound.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Jet {
    /// `g`.
    pub value: Interval,
    /// `∇g`.
    pub gradient: [Interval; 3],
}

/// A compiled gyroid leaf.
#[derive(Clone, Debug)]
pub(crate) struct Gyroid {
    period: Rational,
    grade: GradeDef,
    lipschitz: Rational,
    turns: Interval,
    two_pi: Interval,
    wave: Interval,
    wave_squared: Interval,
    interval_level: Interval,
    sine: [Interval; SINE_TERMS],
    cosine: [Interval; COSINE_TERMS],
    tail_sine: f64,
    tail_cosine: f64,
    fast_turns: f64,
    fast_level: f64,
    fast_lipschitz: f64,
    fast_sine: [f64; SINE_TERMS],
    fast_cosine: [f64; COSINE_TERMS],
}

/// The certified level range and slope of a graded network level.
fn level_limits() -> ([Rational; 2], Rational) {
    let quarter = Rational::from_ratio(1, 4).unwrap();
    ([-&quarter, quarter], Rational::from_ratio(1, 10).unwrap())
}

impl Gyroid {
    /// Compiles the leaf. An impossible period or coincident grade anchors refuse
    /// with `gyroid.impossible`; deferred functions refuse with `gyroid.unsupported`.
    pub(crate) fn new(
        period: &Rational,
        level: &Grade,
        kind: GyroidKind,
        frame: &Pose,
    ) -> Result<Self, Refusal> {
        Self::certified(period, level, kind, frame).map_err(|error| {
            if error.code.starts_with("gyroid.grade-") {
                Refusal::new("gyroid.unsupported", error.detail)
            } else {
                error
            }
        })
    }

    /// As [`Gyroid::new`], but a linear level outside the supported function names
    /// its reason: `gyroid.grade-function-unsupported` for anchors that do not
    /// differ along one part axis, `gyroid.grade-range-unsupported` for an endpoint
    /// outside `[-1/4, 1/4]` or a slope above `1/10` per millimetre.
    pub(crate) fn certified(
        period: &Rational,
        level: &Grade,
        kind: GyroidKind,
        frame: &Pose,
    ) -> Result<Self, Refusal> {
        if period.sign() != Sign::Positive {
            return Err(Refusal::new(
                "gyroid.impossible",
                "a gyroid period must be positive",
            ));
        }
        if kind != GyroidKind::Network {
            return Err(Refusal::new(
                "gyroid.unsupported",
                "sheet gyroids are not supported yet; this field builds the network",
            ));
        }
        let Some(grade) = GradeDef::compile(level, "gyroid.impossible")? else {
            return Err(Refusal::new(
                "gyroid.unsupported",
                "a distance-graded gyroid level is not supported yet; use a constant or linear level",
            ));
        };
        let ([low, high], slope) = level_limits();
        check_certified(
            &grade,
            [
                "gyroid.grade-function-unsupported",
                "gyroid.grade-range-unsupported",
            ],
            [&low, &high],
            &slope,
        )?;
        let [w, x, y, z] = frame.rotation.each_ref().map(|n| n.value().clone());
        if !x.is_zero()
            || !y.is_zero()
            || !z.is_zero()
            || w.sign() != Sign::Positive
            || frame.translation.0.iter().any(|n| !n.value().is_zero())
        {
            return Err(Refusal::new(
                "gyroid.unsupported",
                "a gyroid frame must be the part frame",
            ));
        }
        // A graded level only ever lies in [-1/4, 1/4]; zero stands in where a
        // constant level is read.
        let level = grade.constant().cloned().unwrap_or_else(Rational::zero);
        let lipschitz = &(&Rational::from_i64(22) / period) + &grade.slope_bound();
        let turns = &Rational::one() / period;
        let pi = Ball::pi(PRECISION);
        let two_pi = (&pi + &pi).to_interval();
        let wave = (&(&pi + &pi) * &Ball::from_rational(&turns, PRECISION)).to_interval();
        let sine = alternating_series(1);
        let cosine = alternating_series(0);
        // `turn` admits |u| <= 4/5, so these bound the first omitted terms.
        let tail = |n: i64| {
            let mut bound = Rational::one();
            for i in 1..=n {
                bound = &bound * &Rational::from_ratio(4, 5 * i).unwrap();
            }
            Interval::from_rational(&bound).hi()
        };
        Ok(Self {
            lipschitz: lipschitz.clone(),
            turns: Interval::from_rational(&turns),
            two_pi,
            wave,
            wave_squared: wave.square(),
            interval_level: Interval::from_rational(&level),
            sine: std::array::from_fn(|i| Interval::from_rational(&sine[i])),
            cosine: std::array::from_fn(|i| Interval::from_rational(&cosine[i])),
            tail_sine: tail(19),
            tail_cosine: tail(18),
            fast_turns: turns.to_f64_nearest(),
            fast_level: level.to_f64_nearest(),
            fast_lipschitz: lipschitz.to_f64_nearest(),
            fast_sine: std::array::from_fn(|i| sine[i].to_f64_nearest()),
            fast_cosine: std::array::from_fn(|i| cosine[i].to_f64_nearest()),
            period: period.clone(),
            grade,
        })
    }

    /// An upper bound on the spectral norm of the Hessian of `g`. A row of the
    /// Hessian holds one diagonal entry `k²(sin X cos Y + sin Z cos X)`, at most
    /// `√2 k²` by Cauchy-Schwarz, and two off-diagonal entries of at most `k²`.
    pub(crate) fn curvature(&self) -> Interval {
        self.wave_squared * point(3.5)
    }

    /// An upper bound on `|∇g|`.
    pub(crate) fn lipschitz(&self) -> Interval {
        Interval::from_rational(&self.lipschitz)
    }

    /// The length of one period.
    pub(crate) fn period(&self) -> &Rational {
        &self.period
    }

    /// Whether the level lies at or above the largest value of `g`, so every point is
    /// in the network, or at or below the smallest, so its interior is empty.
    pub(crate) fn saturated(&self) -> Option<bool> {
        let level = self.grade.constant()?;
        let extreme = Rational::from_ratio(3, 2).unwrap();
        if *level >= extreme {
            Some(true)
        } else if *level <= -&extreme {
            Some(false)
        } else {
            None
        }
    }

    /// Whether the level is a non-constant grade.
    pub(crate) fn graded(&self) -> bool {
        self.grade.constant().is_none()
    }

    /// Whether the level is one affine function over a closed rational box, so the
    /// unclamped grade equals the level there.
    pub(crate) fn affine_over(&self, bounds: &[[Rational; 2]; 3]) -> bool {
        self.grade.affine_over(
            &std::array::from_fn(|i| bounds[i][0].clone()),
            &std::array::from_fn(|i| bounds[i][1].clone()),
        )
    }

    /// `sin` and `cos` of one turn coordinate, or `None` beyond reduction range.
    fn turn(&self, turns: Interval) -> Option<[Interval; 2]> {
        let magnitude = turns.lo().abs().max(turns.hi().abs());
        if magnitude >= LARGEST_TURN {
            return None;
        }
        let quarters = (turns.mid() * 4.0 + 0.5).floor();
        let u = (turns - point(quarters) * point(0.25)) * self.two_pi;
        if u.lo() < -0.8 || u.hi() > 0.8 {
            return None;
        }
        let u2 = u.square();
        let mut sine = self.sine[SINE_TERMS - 1];
        for c in self.sine[..SINE_TERMS - 1].iter().rev() {
            sine = sine * u2 + *c;
        }
        let mut cosine = self.cosine[COSINE_TERMS - 1];
        for c in self.cosine[..COSINE_TERMS - 1].iter().rev() {
            cosine = cosine * u2 + *c;
        }
        let sine = sine * u + Interval::new(-self.tail_sine, self.tail_sine).unwrap();
        let cosine = cosine + Interval::new(-self.tail_cosine, self.tail_cosine).unwrap();
        Some(match (quarters as i64).rem_euclid(4) {
            0 => [sine, cosine],
            1 => [cosine, -sine],
            2 => [-sine, -cosine],
            _ => [-cosine, sine],
        })
    }

    /// The grade's exact piece over a closed rational box, or `None` across an
    /// anchor plane. A constant level is one piece everywhere.
    pub(crate) fn level_piece(&self, bounds: &[[Rational; 2]; 3]) -> Option<Piece> {
        self.grade.exact_piece(
            &std::array::from_fn(|i| bounds[i][0].clone()),
            &std::array::from_fn(|i| bounds[i][1].clone()),
        )
    }

    /// The axis and positions of a graded level's anchor planes.
    pub(crate) fn level_planes(&self) -> Option<(usize, [Rational; 2])> {
        self.grade.planes()
    }

    /// Encloses `∇(g - c)` and the Hessian of `g` over a box on which the level is
    /// the given piece. Values at the centre are widened by the next derivative's
    /// bound times the half-widths: second derivatives of `g` are at most `2k²` on
    /// the diagonal and `k²` off it, and third derivatives at most `2k³`.
    pub(crate) fn derivatives(
        &self,
        world: &[Interval; 3],
        piece: Piece,
    ) -> Option<([Interval; 3], [[Interval; 3]; 3])> {
        let centre = world.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
        let radius = std::array::from_fn::<f64, 3, _>(|i| {
            (point(world[i].hi()) - point(centre[i]))
                .hi()
                .max((point(centre[i]) - point(world[i].lo())).hi())
        });
        let c = centre.map(point);
        let [sx, cx] = self.turn(c[0] * self.turns)?;
        let [sy, cy] = self.turn(c[1] * self.turns)?;
        let [sz, cz] = self.turn(c[2] * self.turns)?;
        let k = self.wave;
        let k2 = self.wave_squared;
        let k3 = (k2 * k).hi();
        let slope = self.grade.piece_gradient(piece);
        let r = radius.map(point);
        let reach = (r[0] + r[1] + r[2]).hi();
        let gradient = [
            (cx * cy - sz * sx) * k,
            (cy * cz - sx * sy) * k,
            (cz * cx - sy * sz) * k,
        ];
        let gradient = std::array::from_fn(|i| {
            let spread = (0..3).fold(point(0.0), |sum, j| {
                sum + k2 * point(if i == j { 2.0 } else { 1.0 }) * r[j]
            });
            let spread = spread.hi();
            gradient[i] - slope[i] + Interval::new(-spread, spread).unwrap()
        });
        let hessian = [
            [
                point(0.0) - (sx * cy + sz * cx) * k2,
                point(0.0) - cx * sy * k2,
                point(0.0) - sx * cz * k2,
            ],
            [
                point(0.0) - cx * sy * k2,
                point(0.0) - (sy * cz + sx * cy) * k2,
                point(0.0) - cy * sz * k2,
            ],
            [
                point(0.0) - sx * cz * k2,
                point(0.0) - cy * sz * k2,
                point(0.0) - (sz * cx + sy * cz) * k2,
            ],
        ];
        let spread = (point(2.0 * k3) * point(reach)).hi();
        let hessian = hessian.map(|row| row.map(|v| v + Interval::new(-spread, spread).unwrap()));
        Some((gradient, hessian))
    }

    /// Encloses `g` and `∇g` over a box small enough to be a rounded point.
    fn jet(&self, p: [Interval; 3]) -> Option<Jet> {
        let [sx, cx] = self.turn(p[0] * self.turns)?;
        let [sy, cy] = self.turn(p[1] * self.turns)?;
        let [sz, cz] = self.turn(p[2] * self.turns)?;
        Some(Jet {
            value: sx * cy + sy * cz + sz * cx,
            gradient: [
                (cx * cy - sz * sx) * self.wave,
                (cy * cz - sx * sy) * self.wave,
                (cz * cx - sy * sz) * self.wave,
            ],
        })
    }

    /// `g - ĉ` and its gradient over a box small enough to be a rounded point, where
    /// `ĉ` is the unclamped affine extension of the level. Between the anchor planes
    /// `ĉ` is the level; a constant level is its own extension.
    pub(crate) fn affine_jet(&self, p: [Interval; 3]) -> Option<Jet> {
        let jet = self.jet(p)?;
        let slope = self.grade.gradient();
        Some(Jet {
            value: jet.value - self.grade.affine(&p),
            gradient: std::array::from_fn(|i| jet.gradient[i] - slope[i]),
        })
    }

    /// A second-order expansion of `g - c` about the box centre, for a level that is
    /// `level` at the centre with constant gradient `slope` over the box. Off-diagonal
    /// second derivatives of `g` are at most `k²` and diagonal ones at most `2k²`,
    /// with `k = 2π/a`, so the quadratic term of Taylor's theorem is at most
    /// `k²(Σ rᵢ² + Σ_{i<j} rᵢ rⱼ)` for half-widths `r`.
    fn expansion(
        &self,
        world: &[Interval; 3],
        level: impl FnOnce([Interval; 3]) -> Interval,
        slope: [Interval; 3],
    ) -> Option<Interval> {
        if world
            .iter()
            .any(|b| !b.lo().is_finite() || !b.hi().is_finite())
        {
            return None;
        }
        let centre = world.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
        let radius = std::array::from_fn::<f64, 3, _>(|i| {
            (point(world[i].hi()) - point(centre[i]))
                .hi()
                .max((point(centre[i]) - point(world[i].lo())).hi())
        });
        let jet = self.jet(centre.map(point))?;
        let mut spread = point(0.0);
        for ((slope, grade), width) in jet.gradient.iter().zip(slope).zip(radius) {
            let slope = *slope - grade;
            spread = spread + point(slope.lo().abs().max(slope.hi().abs())) * point(width);
        }
        let r = radius.map(point);
        let curved =
            r[0].square() + r[1].square() + r[2].square() + r[0] * r[1] + r[1] * r[2] + r[2] * r[0];
        let bound = (spread + self.wave_squared * curved).hi();
        let value = jet.value - level(centre.map(point));
        Some(value + Interval::new(-bound, bound).unwrap())
    }

    /// The residual `g - c` over a closed box. Within one piece of the level this is
    /// the second-order expansion; across an anchor plane it is the expansion of `g`
    /// less the exact clamped range of the level over the box.
    pub(crate) fn residual(&self, world: &[Interval; 3]) -> Option<Interval> {
        let zero = [point(0.0); 3];
        if self.grade.constant().is_some() {
            return self.expansion(world, |_| self.interval_level, zero);
        }
        match self.grade.piece(world) {
            Some(piece) => self.expansion(
                world,
                |centre| match piece {
                    Piece::Ramp => self.grade.affine(&centre),
                    Piece::Start | Piece::End => self.grade.interval(&centre),
                },
                self.grade.piece_gradient(piece),
            ),
            None => Some(self.expansion(world, |_| point(0.0), zero)? - self.grade.interval(world)),
        }
    }

    /// The residual `g - ĉ` over a closed box, with `ĉ` the unclamped affine extension
    /// of the level. It is smooth everywhere and equals [`Gyroid::residual`] wherever
    /// the level is affine.
    pub(crate) fn affine_residual(&self, world: &[Interval; 3]) -> Option<Interval> {
        if self.grade.constant().is_some() {
            return self.expansion(world, |_| self.interval_level, [point(0.0); 3]);
        }
        self.expansion(
            world,
            |centre| self.grade.affine(&centre),
            self.grade.gradient(),
        )
    }

    /// The residual `g - ĉ` over a triangle: its centroid's value, the extreme of
    /// the centroid's gradient over the corners, and the Hessian bound over the
    /// corners' largest distance from the centroid. `ĉ` is as in
    /// [`Gyroid::affine_residual`].
    pub(crate) fn facet_residual(&self, corners: &[[Interval; 3]; 3]) -> Option<Interval> {
        let centre: [Interval; 3] =
            std::array::from_fn(|i| (corners[0][i] + corners[1][i] + corners[2][i]) / point(3.0));
        let jet = self.affine_jet(centre)?;
        let mut lo = 0.0_f64;
        let mut hi = 0.0_f64;
        let mut reach = 0.0_f64;
        for corner in corners {
            let d: [Interval; 3] = std::array::from_fn(|i| corner[i] - centre[i]);
            let rate = jet.gradient[0] * d[0] + jet.gradient[1] * d[1] + jet.gradient[2] * d[2];
            lo = lo.min(rate.lo());
            hi = hi.max(rate.hi());
            reach = reach.max((d[0].square() + d[1].square() + d[2].square()).hi());
        }
        let curved = (point(0.5) * self.curvature() * point(reach)).hi();
        Some(jet.value + Interval::new(lo - curved, hi + curved).unwrap())
    }

    /// Exact classification of a rational point.
    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        if let Some(value) = self.residual(&p.each_ref().map(Interval::from_rational)) {
            if value.lo() > 0.0 {
                return Ok(Location::Outside);
            }
            if value.hi() < 0.0 {
                return Ok(Location::Inside);
            }
        }
        let [sines, cosines] = self.exact_terms(p, dec)?;
        let value = &(&(&(&sines[0] * &cosines[1]) + &(&sines[1] * &cosines[2]))
            + &(&sines[2] * &cosines[0]))
            - &Real::from(self.grade.value(p));
        Ok(match geometry::sign(&value, dec)? {
            Sign::Negative => Location::Inside,
            Sign::Zero => Location::Boundary,
            Sign::Positive => Location::Outside,
        })
    }

    fn exact_terms(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<[[Real; 3]; 2], Refusal> {
        let mut sines: [Real; 3] = std::array::from_fn(|_| Real::zero());
        let mut cosines: [Real; 3] = std::array::from_fn(|_| Real::zero());
        for i in 0..3 {
            let degrees = &(&p[i] * &Rational::from_i64(360)) / &self.period;
            let [sine, cosine] = sin_cos_degrees(&degrees, dec)?;
            sines[i] = sine;
            cosines[i] = cosine;
        }
        Ok([sines, cosines])
    }

    /// Enclosure of `(g - c) / L`, negative inside.
    pub(crate) fn proven_distance(&self, p: &RPoint3, location: Location) -> Ball {
        if location == Location::Boundary {
            return Ball::from_rational(&Rational::zero(), PRECISION);
        }
        let pi = Ball::pi(PRECISION);
        let two_pi = &pi + &pi;
        let turn = |coordinate: &Rational| {
            let angle = &two_pi * &Ball::from_rational(&(coordinate / &self.period), PRECISION);
            [angle.sin(), angle.cos()]
        };
        let [sx, cx] = turn(&p.x);
        let [sy, cy] = turn(&p.y);
        let [sz, cz] = turn(&p.z);
        let level = self.grade.value(&[p.x.clone(), p.y.clone(), p.z.clone()]);
        let value = &(&(&(&sx * &cy) + &(&sy * &cz)) + &(&sz * &cx))
            - &Ball::from_rational(&level, PRECISION);
        value
            .checked_div(&Ball::from_rational(&self.lipschitz, PRECISION))
            .expect("positive Lipschitz bound")
    }

    /// Tentative `(g - c) / L` in binary64.
    pub(crate) fn fast_distance(&self, p: [f64; 3]) -> f64 {
        let turn = |coordinate: f64| {
            let turns = coordinate * self.fast_turns;
            let quarters = (turns * 4.0 + 0.5).floor();
            let u = (turns - quarters * 0.25) * core::f64::consts::TAU;
            let u2 = u * u;
            let sine = self.fast_sine.iter().rev().fold(0.0, |s, c| s * u2 + c) * u;
            let cosine = self.fast_cosine.iter().rev().fold(0.0, |s, c| s * u2 + c);
            match (quarters as i64).rem_euclid(4) {
                0 => [sine, cosine],
                1 => [cosine, -sine],
                2 => [-sine, -cosine],
                _ => [-cosine, sine],
            }
        };
        let [sx, cx] = turn(p[0]);
        let [sy, cy] = turn(p[1]);
        let [sz, cz] = turn(p[2]);
        let level = if self.graded() {
            self.grade.fast(p)
        } else {
            self.fast_level
        };
        (sx * cy + sy * cz + sz * cx - level) / self.fast_lipschitz
    }

    /// Positive-infinitesimal ray classification from a rational point.
    pub(crate) fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let location = self.location(p, dec)?;
        if location != Location::Boundary {
            return Ok(location);
        }
        let [s, c] = self.exact_terms(p, dec)?;
        // `∇g · d` without its positive factor `k = 2π/a`.
        let rate = &(&(&d.x * &(&(&c[0] * &c[1]) - &(&s[2] * &s[0])))
            + &(&d.y * &(&(&c[1] * &c[2]) - &(&s[0] * &s[1]))))
            + &(&d.z * &(&(&c[2] * &c[0]) - &(&s[1] * &s[2])));
        let along = geometry::sign(&rate, dec)?;
        let Some(change) = self.grade.rate(p, d, dec)? else {
            return Ok(match along {
                Sign::Negative => Location::Inside,
                Sign::Positive => Location::Outside,
                Sign::Zero => Location::Indeterminate,
            });
        };
        // The rate of `g - c` is `k (∇g · d) - ∇c · d`. Both terms are algebraic and
        // `k` is not, so it vanishes only when both terms do.
        let grade = geometry::sign(&change, dec)?;
        let sign = match (along, grade) {
            (Sign::Zero, Sign::Zero) => Sign::Zero,
            (along, Sign::Zero) => along,
            (Sign::Zero, Sign::Positive) => Sign::Negative,
            (Sign::Zero, Sign::Negative) => Sign::Positive,
            (along, grade) if along != grade => along,
            _ => self.compare_wave(&rate, &change, dec)?,
        };
        Ok(match sign {
            Sign::Negative => Location::Inside,
            Sign::Positive => Location::Outside,
            Sign::Zero => Location::Indeterminate,
        })
    }

    /// The sign of `k a - b` for nonzero algebraic `a` and `b`, with `k = 2π/a`. The
    /// two terms are never equal, so refinement separates them; a precision limit
    /// leaves the direction unresolved.
    fn compare_wave(&self, a: &Real, b: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
        for bits in [64u32, 128, 256, 512] {
            let width = Rational::one().mul_pow2(-i64::from(bits) / 2);
            let a = dec.enclose(a, &width)?;
            let b = dec.enclose(b, &width)?;
            let pi = Ball::pi(bits);
            let k = &(&pi + &pi) * &Ball::from_rational(&(&Rational::one() / &self.period), bits);
            let difference = &(&k * &a) - &b;
            if difference.lo() > Rational::zero() {
                return Ok(Sign::Positive);
            }
            if difference.hi() < Rational::zero() {
                return Ok(Sign::Negative);
            }
        }
        Ok(Sign::Zero)
    }

    /// A copy would leave the part frame, which this increment does not support.
    pub(crate) fn transform(&mut self, _map: &Isometry) -> Result<(), Refusal> {
        Err(Refusal::new(
            "gyroid.unsupported",
            "a gyroid frame must be the part frame",
        ))
    }

    /// Interval classification of a world box.
    pub(crate) fn box_state(
        &self,
        world: &[Interval; 3],
        dec: &mut Decider,
    ) -> Result<BoxState, Refusal> {
        dec.check_deadline()?;
        Ok(match self.residual(world) {
            Some(value) if value.lo() > 0.0 => BoxState::Outside,
            Some(value) if value.hi() < 0.0 => BoxState::Inside,
            _ => BoxState::Unknown,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_recipe::Number;

    fn q(value: &str) -> Rational {
        Rational::parse(value).unwrap()
    }
    fn gyroid(level: &str) -> Gyroid {
        let level = serde_json::from_str(&format!(r#"{{"constant":"{level}"}}"#)).unwrap();
        let frame =
            serde_json::from_str(r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#)
                .unwrap();
        Gyroid::new(
            Number::parse("10").unwrap().value(),
            &level,
            GyroidKind::Network,
            &frame,
        )
        .unwrap()
    }
    fn tenth_box(centre: [&str; 3]) -> [Interval; 3] {
        centre.map(|c| {
            let c = q(c);
            Interval::new(
                Interval::from_rational(&(&c - &q("1/10"))).lo(),
                Interval::from_rational(&(&c + &q("1/10"))).hi(),
            )
            .unwrap()
        })
    }

    /// K-GB, fixed by the independent V7c designer. Whole-cell exclusion rests on these
    /// dispositions: at level 1/2 the origin box is wholly material and the box about
    /// `(5/2, 0, 0)` wholly void, while at level 0 the origin box crosses the surface.
    #[test]
    fn gyroid_classifier_respects_containment_and_crossing() {
        let mut dec = Decider::default();
        assert_eq!(
            gyroid("1/2")
                .box_state(&tenth_box(["0", "0", "0"]), &mut dec)
                .unwrap(),
            BoxState::Inside
        );
        assert_eq!(
            gyroid("1/2")
                .box_state(&tenth_box(["5/2", "0", "0"]), &mut dec)
                .unwrap(),
            BoxState::Outside
        );
        assert_eq!(
            gyroid("0")
                .box_state(&tenth_box(["0", "0", "0"]), &mut dec)
                .unwrap(),
            BoxState::Unknown
        );
    }

    /// The polynomial enclosures hold the arbitrary-precision sine and cosine, across
    /// every quarter-turn reduction and at reduction boundaries.
    #[test]
    fn interval_jets_enclose_arbitrary_precision_values() {
        let field = gyroid("0");
        let pi = Ball::pi(PRECISION);
        let two_pi = &pi + &pi;
        for numerator in -1_000..=1_000 {
            let coordinate = Rational::from_ratio(numerator * 41, 256).unwrap();
            let p = [
                Interval::from_rational(&coordinate),
                Interval::from_rational(&(&coordinate + &q("7/3"))),
                Interval::from_rational(&(&coordinate * &q("-3/5"))),
            ];
            let exact = [
                coordinate.clone(),
                &coordinate + &q("7/3"),
                &coordinate * &q("-3/5"),
            ]
            .map(|x| &two_pi * &Ball::from_rational(&(&x / &q("10")), PRECISION));
            let value = &(&(&exact[0].sin() * &exact[1].cos())
                + &(&exact[1].sin() * &exact[2].cos()))
                + &(&exact[2].sin() * &exact[0].cos());
            let jet = field.jet(p).unwrap();
            let value = value.to_interval();
            assert!(jet.value.lo() <= value.lo() && jet.value.hi() >= value.hi());
            assert!(jet.value.width() < 1e-12);
        }
    }

    /// The box bound holds a dense sample of the box, where the second-order term
    /// and the sign of each first-order term all matter.
    #[test]
    fn box_residuals_enclose_sampled_values() {
        let field = gyroid("1/4");
        for corner in 0..200 {
            let lo = [
                corner as f64 * 0.137 - 9.0,
                corner as f64 * -0.091,
                2.0 - corner as f64 * 0.053,
            ];
            let world: [Interval; 3] =
                std::array::from_fn(|i| Interval::new(lo[i], lo[i] + 0.5).unwrap());
            let range = field.residual(&world).unwrap();
            for sample in 0..64 {
                let at: [f64; 3] = std::array::from_fn(|i| {
                    lo[i]
                        + 0.5 * f64::from((sample >> i) & 1) * (1.0 - f64::from(sample % 7) / 13.0)
                });
                let value = field.jet(at.map(point)).unwrap().value - field.interval_level;
                assert!(range.lo() <= value.hi() && range.hi() >= value.lo());
            }
        }
    }
}
