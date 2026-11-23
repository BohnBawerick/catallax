//! Clamped grades shared by lattice radii and gyroid levels.
//!
//! Linear anchors, centres and axes are fixed in part coordinates. Face-distance
//! grades resolve against the current region. The grade's endpoint values and slope
//! determine admission to terminal certification.

use crate::geometry;
use crate::metric::{PRECISION, Scalar};
use catallax_features::grade::{Linear, rational_dot};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_recipe::Grade;
use catallax_shapes::Vector3;

mod faces;
mod radius;
use radius::Radius;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn f(r: &Rational) -> f64 {
    r.to_f64_nearest()
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

/// On a rectangular body the distance to each trimmed face is its inward
/// plane distance. Outside the clip, the clamped planes extend the grade
/// continuously without changing the clipped solid.
#[derive(Clone, Debug)]
pub(crate) struct FaceDistance {
    planes: Vec<Linear>,
    slope: Rational,
}

/// A compiled grade.
#[derive(Clone, Debug)]
pub(crate) enum GradeDef {
    Constant(Rational),
    Linear(Box<Linear>),
    Faces(Box<FaceDistance>),
    Radius(Box<Radius>),
}

/// The piece of a clamped linear grade that holds over a whole box.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Piece {
    /// The constant start value, at or before the start anchor plane.
    Start,
    /// The linear ramp, between the anchor planes inclusive.
    Ramp,
    /// The constant end value, at or beyond the end anchor plane.
    End,
    /// The ramp on one selected rectangular face's side of the bisectors.
    Face(usize),
    /// A smooth ramp in the distance to a centre or axis.
    Curved,
}

impl GradeDef {
    pub(crate) fn from_faces(
        planes: impl IntoIterator<Item = (usize, bool, Rational)>,
        near: &Rational,
        far: &Rational,
        span: &Rational,
    ) -> Self {
        if near == far {
            return Self::Constant(near.clone());
        }
        let mut planes: Vec<_> = planes
            .into_iter()
            .map(|(axis, upper, position)| {
                let mut from = [q(0), q(0), q(0)];
                let mut direction = [q(0), q(0), q(0)];
                from[axis] = position;
                direction[axis] = if upper { -span } else { span.clone() };
                Linear {
                    from,
                    direction,
                    length_sq: span * span,
                    start: near.clone(),
                    end: far.clone(),
                }
            })
            .collect();
        if planes.len() == 1 {
            return Self::Linear(Box::new(planes.pop().expect("one selected face")));
        }
        Self::Faces(Box::new(FaceDistance {
            planes,
            slope: &(far - near).abs() / span,
        }))
    }

    /// Compiles a grade with no region dependency. Equal endpoint values give the
    /// constant. Face-distance grades are resolved separately against a region.
    pub(crate) fn compile(
        grade: &Grade,
        impossible: &'static str,
    ) -> Result<Option<Self>, Refusal> {
        Ok(Some(match grade {
            Grade::Constant(c) => Self::Constant(c.value().clone()),
            Grade::Linear {
                from,
                to,
                start,
                end,
            } => {
                if from == to {
                    return Err(Refusal::new(
                        impossible,
                        "a linear grade needs two distinct points",
                    ));
                }
                if start == end {
                    return Ok(Some(Self::Constant(start.value().clone())));
                }
                Self::Linear(Box::new(Linear::new(
                    from.0.each_ref().map(|n| n.value().clone()),
                    to.0.each_ref().map(|n| n.value().clone()),
                    start.value().clone(),
                    end.value().clone(),
                )))
            }
            Grade::Distance { .. } => return Ok(None),
            Grade::Point {
                centre,
                r0,
                r1,
                start,
                end,
            } => {
                catallax_features::lattice::validate_distance_radii(grade, impossible)?;
                Self::from_radius(
                    &centre.0,
                    None,
                    [r0.value(), r1.value()],
                    [start.value(), end.value()],
                )
            }
            Grade::Radial {
                origin,
                direction,
                r0,
                r1,
                start,
                end,
            } => {
                catallax_features::lattice::validate_distance_radii(grade, impossible)?;
                Self::from_radius(
                    &origin.0,
                    Some(&direction.0),
                    [r0.value(), r1.value()],
                    [start.value(), end.value()],
                )
            }
        }))
    }

    fn from_radius(
        origin: &[catallax_recipe::Number; 3],
        direction: Option<&[catallax_recipe::Number; 3]>,
        [r0, r1]: [&Rational; 2],
        [start, end]: [&Rational; 2],
    ) -> Self {
        if start == end {
            return Self::Constant(start.clone());
        }
        let values = |p: &[catallax_recipe::Number; 3]| p.each_ref().map(|n| n.value().clone());
        Self::Radius(Box::new(Radius::new(
            values(origin),
            direction.map(values),
            [r0.clone(), r1.clone()],
            [start.clone(), end.clone()],
        )))
    }

    /// The constant value, or `None` for a ramp.
    pub(crate) fn constant(&self) -> Option<&Rational> {
        match self {
            Self::Constant(c) => Some(c),
            Self::Linear(_) | Self::Faces(_) | Self::Radius(_) => None,
        }
    }

    /// The part axis the anchors differ along, when they differ along exactly one.
    pub(crate) fn axis(&self) -> Option<usize> {
        let Self::Linear(linear) = self else {
            return None;
        };
        let mut axes = (0..3).filter(|&i| !linear.direction[i].is_zero());
        let axis = axes.next()?;
        axes.next().is_none().then_some(axis)
    }

    /// A rational upper bound on the slope `|v1 - v0| / |g - f|`, exact when the
    /// anchors differ along one axis. Zero for a constant.
    pub(crate) fn slope_bound(&self) -> Rational {
        if let Self::Faces(faces) = self {
            return faces.slope.clone();
        }
        if let Self::Radius(radius) = self {
            return radius.ramp.slope.abs();
        }
        let Self::Linear(linear) = self else {
            return q(0);
        };
        let rise = (&linear.end - &linear.start).abs();
        let length = match self.axis() {
            Some(axis) => linear.direction[axis].abs(),
            None => Ball::from_rational(&linear.length_sq, PRECISION)
                .sqrt()
                .expect("positive squared length")
                .lo(),
        };
        rise.checked_div(&length)
            .expect("a positive lower bound on a nonzero length")
    }

    /// The exact value for grades whose value at a rational point is rational.
    pub(crate) fn rational_value(&self, p: &[Rational; 3]) -> Option<Rational> {
        Some(match self {
            Self::Constant(c) => c.clone(),
            Self::Linear(linear) => linear.interpolate(linear.parameter(p)),
            Self::Faces(faces) => faces.planes[0].interpolate(
                faces
                    .planes
                    .iter()
                    .map(|plane| plane.parameter(p))
                    .min()
                    .expect("selected faces"),
            ),
            Self::Radius(_) => return None,
        })
    }

    pub(crate) fn real(&self, p: &[Rational; 3]) -> Real {
        match self {
            Self::Radius(radius) => radius.real(p),
            _ => Real::from(self.rational_value(p).expect("rational grade")),
        }
    }

    pub(crate) fn ball(&self, p: &[Rational; 3]) -> Ball {
        match self {
            Self::Radius(radius) => radius.ball(p),
            _ => Ball::from_rational(&self.rational_value(p).expect("rational grade"), PRECISION),
        }
    }

    fn fast_parameter(linear: &Linear, p: [f64; 3]) -> f64 {
        let rel: [f64; 3] = std::array::from_fn(|i| p[i] - f(&linear.from[i]));
        let direction = linear.direction.each_ref().map(f);
        rel[0]
            .mul(&direction[0])
            .add(&rel[1].mul(&direction[1]))
            .add(&rel[2].mul(&direction[2]))
            / f(&linear.length_sq)
    }

    /// The tentative value in binary64.
    pub(crate) fn fast(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Constant(c) => f(c),
            Self::Linear(linear) => Self::fast_interpolate(linear, Self::fast_parameter(linear, p)),
            Self::Faces(faces) => Self::fast_interpolate(
                &faces.planes[0],
                faces
                    .planes
                    .iter()
                    .map(|plane| Self::fast_parameter(plane, p))
                    .fold(f64::INFINITY, f64::min),
            ),
            Self::Radius(radius) => radius.fast(p),
        }
    }

    fn fast_interpolate(linear: &Linear, t: f64) -> f64 {
        f(&linear.start) + t.clamp(0.0, 1.0) * (f(&linear.end) - f(&linear.start))
    }

    fn interval_parameter(linear: &Linear, p: &[Interval; 3]) -> Interval {
        let rel: [Interval; 3] =
            std::array::from_fn(|i| p[i] - Interval::from_rational(&linear.from[i]));
        (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
            s + rel[i] * Interval::from_rational(&linear.direction[i])
        }) / Interval::from_rational(&linear.length_sq)
    }

    /// An enclosure of every value over a box.
    pub(crate) fn interval(&self, p: &[Interval; 3]) -> Interval {
        match self {
            Self::Constant(c) => Interval::from_rational(c),
            Self::Linear(linear) => {
                Self::interval_interpolate(linear, clamp01(Self::interval_parameter(linear, p)))
            }
            Self::Faces(faces) => Self::interval_interpolate(
                &faces.planes[0],
                clamp01(
                    faces
                        .planes
                        .iter()
                        .map(|plane| Self::interval_parameter(plane, p))
                        .reduce(imin)
                        .expect("selected faces"),
                ),
            ),
            Self::Radius(radius) => radius.interval(p),
        }
    }

    fn interval_interpolate(linear: &Linear, t: Interval) -> Interval {
        Interval::from_rational(&linear.start)
            + t * (Interval::from_rational(&linear.end) - Interval::from_rational(&linear.start))
    }

    /// The piece that holds over the whole box, or `None` when the box meets an
    /// anchor plane's two sides. A constant is one ramp-free piece everywhere.
    pub(crate) fn piece(&self, p: &[Interval; 3]) -> Option<Piece> {
        if let Self::Radius(radius) = self {
            return radius.piece(p);
        }
        if let Self::Faces(faces) = self {
            return face_piece(
                &faces.planes,
                |plane| {
                    let t = Self::interval_parameter(plane, p);
                    [t.lo(), t.hi()]
                },
                [0.0, 1.0],
            );
        }
        let Self::Linear(linear) = self else {
            return Some(Piece::Start);
        };
        let t = Self::interval_parameter(linear, p);
        if t.hi() <= 0.0 {
            Some(Piece::Start)
        } else if t.lo() >= 1.0 {
            Some(Piece::End)
        } else if t.lo() >= 0.0 && t.hi() <= 1.0 {
            Some(Piece::Ramp)
        } else {
            None
        }
    }

    pub(crate) fn has_affine_extension(&self) -> bool {
        matches!(self, Self::Constant(_) | Self::Linear(_))
    }

    /// The unclamped affine extension of the ramp over a box: `v0 + t (v1 - v0)`.
    /// It equals the grade between the anchor planes. A constant returns itself.
    pub(crate) fn affine(&self, p: &[Interval; 3]) -> Interval {
        match self {
            Self::Constant(c) => Interval::from_rational(c),
            Self::Linear(linear) => {
                Interval::from_rational(&linear.start)
                    + Self::interval_parameter(linear, p)
                        * (Interval::from_rational(&linear.end)
                            - Interval::from_rational(&linear.start))
            }
            Self::Faces(_) | Self::Radius(_) => unreachable!("this grade has no affine extension"),
        }
    }

    /// The gradient of the affine extension, enclosed. Zero for a constant.
    pub(crate) fn gradient(&self) -> [Interval; 3] {
        match self {
            Self::Constant(_) => [Interval::point(0.0).unwrap(); 3],
            Self::Linear(linear) => Self::linear_gradient(linear),
            Self::Faces(_) | Self::Radius(_) => unreachable!("this grade has no affine extension"),
        }
    }

    fn linear_gradient(linear: &Linear) -> [Interval; 3] {
        let scale = (&linear.end - &linear.start)
            .checked_div(&linear.length_sq)
            .expect("distinct grade points");
        std::array::from_fn(|i| Interval::from_rational(&(&linear.direction[i] * &scale)))
    }

    pub(crate) fn piece_value(&self, piece: Piece, p: &[Interval; 3]) -> Interval {
        match piece {
            Piece::Ramp => self.affine(p),
            Piece::Start | Piece::End | Piece::Curved => self.interval(p),
            Piece::Face(i) => {
                let Self::Faces(faces) = self else {
                    unreachable!("face grade piece")
                };
                let plane = &faces.planes[i];
                Self::interval_interpolate(plane, Self::interval_parameter(plane, p))
            }
        }
    }

    /// The gradient of one piece, enclosed: the ramp's, or zero on a clamp.
    pub(crate) fn piece_gradient(&self, piece: Piece) -> [Interval; 3] {
        match piece {
            Piece::Ramp => self.gradient(),
            Piece::Start | Piece::End | Piece::Curved => [Interval::point(0.0).unwrap(); 3],
            Piece::Face(i) => {
                let Self::Faces(faces) = self else {
                    unreachable!("face grade piece")
                };
                Self::linear_gradient(&faces.planes[i])
            }
        }
    }

    pub(crate) fn curved_derivatives(
        &self,
        p: &[Interval; 3],
    ) -> Option<([Interval; 3], [[Interval; 3]; 3])> {
        let Self::Radius(radius) = self else {
            return None;
        };
        radius.derivatives(p)
    }

    /// The exact range of `t` over a closed rational box.
    fn parameter_range(linear: &Linear, lo: &[Rational; 3], hi: &[Rational; 3]) -> [Rational; 2] {
        // `t` is affine, so its extremes are at the corners selected by the signs of
        // the direction's components.
        let corner = |upper: bool| -> [Rational; 3] {
            std::array::from_fn(|i| {
                if (linear.direction[i].sign() == Sign::Positive) == upper {
                    hi[i].clone()
                } else {
                    lo[i].clone()
                }
            })
        };
        [corner(false), corner(true)].map(|p| linear.parameter(&p))
    }

    /// The piece that holds over a whole closed rational box, decided exactly, or
    /// `None` when the box meets both sides of an anchor plane.
    pub(crate) fn exact_piece(&self, lo: &[Rational; 3], hi: &[Rational; 3]) -> Option<Piece> {
        if let Self::Radius(radius) = self {
            return radius.piece(&interval_box(lo, hi));
        }
        if let Self::Faces(faces) = self {
            return face_piece(
                &faces.planes,
                |plane| Self::parameter_range(plane, lo, hi),
                [q(0), q(1)],
            );
        }
        let Self::Linear(linear) = self else {
            return Some(Piece::Start);
        };
        let [first, last] = Self::parameter_range(linear, lo, hi);
        if last <= q(0) {
            Some(Piece::Start)
        } else if first >= q(1) {
            Some(Piece::End)
        } else if first >= q(0) && last <= q(1) {
            Some(Piece::Ramp)
        } else {
            None
        }
    }

    /// The axis and positions of the two anchor planes of a grade whose anchors
    /// differ along one part axis.
    pub(crate) fn planes(&self) -> Option<(usize, [Rational; 2])> {
        let axis = self.axis()?;
        let Self::Linear(linear) = self else {
            return None;
        };
        let from = linear.from[axis].clone();
        let to = &from + &linear.direction[axis];
        Some((axis, [from, to]))
    }

    /// Bounds a grade over a closed rational box. Constants and linear grades
    /// return the exact extremes.
    pub(crate) fn range(&self, lo: &[Rational; 3], hi: &[Rational; 3]) -> [Rational; 2] {
        match self {
            Self::Constant(c) => [c.clone(), c.clone()],
            Self::Faces(faces) => {
                let [lo, hi] = faces
                    .planes
                    .iter()
                    .map(|plane| Self::parameter_range(plane, lo, hi))
                    .reduce(|[a, b], [c, d]| [a.min(c), b.min(d)])
                    .expect("selected faces");
                let [a, b] = [lo, hi].map(|t| faces.planes[0].interpolate(t));
                [a.clone().min(b.clone()), a.max(b)]
            }
            Self::Linear(linear) => {
                let [a, b] = Self::parameter_range(linear, lo, hi).map(|t| linear.interpolate(t));
                [a.clone().min(b.clone()), a.max(b)]
            }
            Self::Radius(radius) => {
                let value = radius.interval(&interval_box(lo, hi));
                let [a, b] = &radius.ramp.values;
                [
                    Rational::from_f64_exact(value.lo()).unwrap_or_else(|| a.min(b).clone()),
                    Rational::from_f64_exact(value.hi()).unwrap_or_else(|| a.max(b).clone()),
                ]
            }
        }
    }

    /// The rate of change along `d` in the positive-infinitesimal direction from `p`,
    /// or `None` where the grade is constant that way.
    pub(crate) fn rate(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        if let Self::Faces(faces) = self {
            return faces.rate(p, d, dec);
        }
        if let Self::Radius(radius) = self {
            return radius.rate(p, d, dec);
        }
        let Self::Linear(linear) = self else {
            return Ok(None);
        };
        let Linear {
            direction,
            length_sq,
            start,
            end,
            ..
        } = linear.as_ref();
        let t = linear.parameter(p);
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
        Ok(ramp.then(|| {
            &Real::from(
                (end - start)
                    .checked_div(length_sq)
                    .expect("distinct grade points"),
            ) * &along
        }))
    }

    /// Composes a rational rigid motion `x -> M x + t` into the anchors.
    pub(crate) fn transform(&mut self, m: &[[Rational; 3]; 3], t: &[Rational; 3]) {
        let planes = match self {
            Self::Constant(_) => return,
            Self::Linear(linear) => std::slice::from_mut(linear.as_mut()),
            Self::Faces(faces) => &mut faces.planes,
            Self::Radius(radius) => return radius.transform(m, t),
        };
        for linear in planes {
            let apply = |v: &[Rational; 3]| -> [Rational; 3] {
                std::array::from_fn(|i| rational_dot(&m[i], v))
            };
            let moved = apply(&linear.from);
            linear.from = std::array::from_fn(|i| &moved[i] + &t[i]);
            linear.direction = apply(&linear.direction);
        }
    }
}

fn interval_box(lo: &[Rational; 3], hi: &[Rational; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| Interval::from_rational(&lo[i]).hull(&Interval::from_rational(&hi[i])))
}

fn face_piece<T: PartialOrd>(
    planes: &[Linear],
    range: impl Fn(&Linear) -> [T; 2],
    [zero, one]: [T; 2],
) -> Option<Piece> {
    let ranges = planes.iter().map(&range);
    if ranges.clone().any(|[_, hi]| hi <= zero) {
        return Some(Piece::Start);
    }
    if ranges.clone().all(|[lo, _]| lo >= one) {
        return Some(Piece::End);
    }
    for (i, [lo, hi]) in ranges.clone().enumerate() {
        if lo >= zero
            && hi <= one
            && ranges
                .clone()
                .enumerate()
                .all(|(j, [low, _])| i == j || hi <= low)
        {
            return Some(Piece::Face(i));
        }
    }
    None
}

impl FaceDistance {
    fn rate(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let parameters: Vec<_> = self.planes.iter().map(|plane| plane.parameter(p)).collect();
        let first = parameters.iter().min().expect("selected faces");
        if first < &q(0) || first > &q(1) {
            return Ok(None);
        }
        let mut rate: Option<Real> = None;
        for (plane, parameter) in self.planes.iter().zip(&parameters) {
            if parameter != first {
                continue;
            }
            let [x, y, z] = plane.direction.clone().map(Real::from);
            let along = Vector3::new(x, y, z).dot(d);
            if let Some(previous) = &rate
                && geometry::sign(&(&along - previous), dec)? != Sign::Negative
            {
                continue;
            }
            rate = Some(along);
        }
        let rate = rate.expect("one closest face");
        if (first == &q(0) && geometry::sign(&rate, dec)? != Sign::Positive)
            || (first == &q(1) && geometry::sign(&rate, dec)? != Sign::Negative)
        {
            return Ok(None);
        }
        let plane = &self.planes[0];
        Ok(Some(
            &rate * &Real::from(&(&plane.end - &plane.start) / &plane.length_sq),
        ))
    }
}

/// Checks the endpoint and slope limits of one field family. Constants keep
/// their own route. Face distances check these limits before face resolution.
pub(crate) fn check_certified(
    grade: &GradeDef,
    range: &'static str,
    [low, high]: [&Rational; 2],
    max_slope: &Rational,
) -> Result<(), Refusal> {
    if let GradeDef::Radius(radius) = grade {
        return check_limits(
            radius.ramp.values.each_ref(),
            [low, high],
            || radius.ramp.slope.abs(),
            max_slope,
            range,
        );
    }
    let GradeDef::Linear(linear) = grade else {
        return Ok(());
    };
    if grade.axis().is_none() {
        check_endpoints([&linear.start, &linear.end], [low, high], range)?;
        let rise = &linear.end - &linear.start;
        let allowed_sq = &(max_slope * max_slope) * &linear.length_sq;
        if &rise * &rise > allowed_sq {
            return Err(slope_refusal(&grade.slope_bound(), max_slope, range));
        }
        return Ok(());
    }
    check_limits(
        [&linear.start, &linear.end],
        [low, high],
        || grade.slope_bound(),
        max_slope,
        range,
    )
}

fn check_limits(
    endpoints: [&Rational; 2],
    [low, high]: [&Rational; 2],
    slope: impl FnOnce() -> Rational,
    max_slope: &Rational,
    range: &'static str,
) -> Result<(), Refusal> {
    check_endpoints(endpoints, [low, high], range)?;
    let slope = slope();
    if &slope > max_slope {
        return Err(slope_refusal(&slope, max_slope, range));
    }
    Ok(())
}

fn slope_refusal(slope: &Rational, max_slope: &Rational, range: &'static str) -> Refusal {
    Refusal::new(
        range,
        format!("grade slope {slope} per mm exceeds the certified {max_slope} per mm"),
    )
}

fn check_endpoints(
    endpoints: [&Rational; 2],
    [low, high]: [&Rational; 2],
    range: &'static str,
) -> Result<(), Refusal> {
    for value in endpoints {
        if value < low || value > high {
            return Err(Refusal::new(
                range,
                format!("grade endpoint {value} is outside the certified range [{low}, {high}]"),
            ));
        }
    }
    Ok(())
}

pub(crate) fn check_distance(
    [near, far, span]: [&Rational; 3],
    [low, high]: [&Rational; 2],
    max_slope: &Rational,
    [impossible, range]: [&'static str; 2],
) -> Result<(), Refusal> {
    if span.sign() != Sign::Positive {
        return Err(Refusal::new(
            impossible,
            "a distance grade needs a positive span",
        ));
    }
    if near == far {
        return Ok(());
    }
    check_limits(
        [near, far],
        [low, high],
        || &(far - near).abs() / span,
        max_slope,
        range,
    )
}

/// Bounds a grade over a closed rational box. Constant and linear grades return
/// their exact clamped range. Point and radial grades return an enclosure.
///
/// Refusals: `field.grade-box` when a lower corner exceeds its upper corner,
/// `field.grade-unsupported` for a grade needing region faces,
/// `field.grade-impossible` for invalid grade parameters, and the decider's
/// work and deadline outcomes.
pub fn grade_range(
    grade: &Grade,
    lo: &[Rational; 3],
    hi: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    dec.check_deadline()?;
    if (0..3).any(|i| lo[i] > hi[i]) {
        return Err(Refusal::new(
            "field.grade-box",
            "each lower box corner must not exceed its upper corner",
        ));
    }
    let Some(compiled) = GradeDef::compile(grade, "field.grade-impossible")? else {
        return Err(Refusal::new(
            "field.grade-unsupported",
            "distance grades need trimmed-face distance; unfinished coverage",
        ));
    };
    Ok(compiled.range(lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grade(text: &str) -> Grade {
        serde_json::from_str(text).unwrap()
    }
    fn r(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn exact_range_clamps_across_an_anchor_plane() {
        let level = grade(
            r#"{"linear":{"from":["0","0","-1"],"to":["0","0","1"],"start":"-1/10","end":"1/10"}}"#,
        );
        let range = grade_range(
            &level,
            &[r("0"), r("0"), r("999/1000")],
            &[r("0"), r("0"), r("1001/1000")],
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(range, [r("999/10000"), r("1/10")]);
    }

    #[test]
    fn equal_endpoints_compile_to_the_constant() {
        let equal = grade(
            r#"{"linear":{"from":["0","0","20"],"to":["0","0","30"],"start":"1","end":"1"}}"#,
        );
        let compiled = GradeDef::compile(&equal, "lattice.impossible")
            .unwrap()
            .unwrap();
        assert_eq!(compiled.constant(), Some(&r("1")));
    }
}
