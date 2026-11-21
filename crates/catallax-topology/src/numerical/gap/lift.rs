//! Binary64 forms of the two copies of an edge, in a local frame.
//!
//! Rule T1 of the tolerance rules runs every operation in a local frame
//! reached by an exact rational transform. Here that transform is a shift:
//! one rational point near the edge is subtracted from every coordinate, and
//! one rational u/v point from the parameters of planes and cylinders, before
//! anything is rounded to binary64. A gap of 1e-5 mm at 1e12 mm from the
//! world origin is then measured at the scale of the edge.

use super::numeric::{
    TRIG_ROUNDINGS, UNIT_ROUNDOFF, basis_functions, find_span, roundings, sin_cos,
};
use super::{
    EdgeCopies, FaceCopyRef, GapCurve, GapSurface, SampledCurve, SampledSurface, unsupported,
};
use crate::numerical::BITS;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface, Curve, Point3, Surface, Vector3};

/// A spline curve with homogeneous controls; the last coordinate is the weight.
pub(super) struct Spline<const D: usize> {
    pub degree: usize,
    pub knots: Vec<f64>,
    /// Each exact knot minus its binary64 value.
    low: Vec<f64>,
    pub controls: Vec<[f64; D]>,
    /// Largest Cartesian coordinate magnitude of each control.
    reach: Vec<f64>,
}

/// A spline surface with homogeneous controls in u-major order.
pub(super) struct Net {
    pub degrees: [usize; 2],
    pub knots: [Vec<f64>; 2],
    low: [Vec<f64>; 2],
    pub count: [usize; 2],
    controls: Vec<[f64; 4]>,
    reach: Vec<f64>,
}

pub(super) enum Curve64 {
    /// `origin + t dir`.
    Line {
        origin: [f64; 3],
        dir: [f64; 3],
    },
    /// `center + cos t major + sin t minor`: a circle or an ellipse.
    Conic {
        center: [f64; 3],
        major: [f64; 3],
        minor: [f64; 3],
    },
    Spline(Spline<4>),
    /// A curve with no kernel form, read through its own evaluator in the
    /// frame that evaluator works in.
    Sampled(SampledCurve),
}

/// The kinds of surface of revolution and their two radii.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Revolved {
    Cylinder,
    Cone,
    Sphere,
    Torus,
    EllipticCylinder,
}

pub(super) enum Surface64 {
    /// `origin + u e1 + v e2`, with `origin` already at the u/v shift.
    Plane {
        origin: [f64; 3],
        frame: [[f64; 3]; 3],
    },
    /// The kernel's charts of the analytic curved surfaces.
    Revolved {
        kind: Revolved,
        origin: [f64; 3],
        frame: [[f64; 3]; 3],
        /// Radius; slope; radius; major and minor; major and minor semi-axis.
        sizes: [f64; 2],
    },
    Spline(Net),
    /// A surface with no kernel form, read through its own evaluator in the
    /// frame that evaluator works in.
    Sampled(SampledSurface),
}

pub(super) enum Face64 {
    Native(Spline<3>),
    /// Straight in u/v from the first edge parameter to the second.
    Linear {
        ends: [[f64; 2]; 2],
        /// Whether the copy keeps its v and moves its u at exactly the rate
        /// of the edge parameter, forwards or backwards. Decided on the
        /// exact values, before any rounding.
        turns_with_edge: bool,
    },
    /// The foot of the 3D point on a plane.
    Projected,
}

/// Both copies and the surface, ready to evaluate.
pub(super) struct Lifted {
    pub curve: Curve64,
    pub face: Face64,
    pub surface: Surface64,
    /// Edge parameters in traversal order.
    pub parameters: [f64; 2],
    /// Each exact edge parameter minus its binary64 value.
    low: [f64; 2],
    /// The exact u/v shift subtracted from the face's parameters.
    pub shift: [Rational; 2],
}

/// Reusable basis-function storage.
#[derive(Default)]
pub(super) struct Scratch {
    u: Vec<f64>,
    v: Vec<f64>,
    work: Vec<f64>,
}

/// A value with the size of its binary64 rounding error.
pub(super) struct Bounded<T> {
    pub value: T,
    /// An estimate of the absolute rounding error of each coordinate.
    pub error: f64,
}

fn rational(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    match value.rational_ref() {
        Some(exact) => Ok(exact.clone()),
        None => Ok(dec.enclose_budgeted_at(value, BITS)?.mid()),
    }
}

fn rational_point(p: &Point3, dec: &mut Decider) -> Result<[Rational; 3], Refusal> {
    Ok([
        rational(&p.x, dec)?,
        rational(&p.y, dec)?,
        rational(&p.z, dec)?,
    ])
}

/// Binary64 values of exact numbers, and each number minus its value.
fn split<'a>(exact: impl Iterator<Item = &'a Rational>) -> (Vec<f64>, Vec<f64>) {
    exact
        .map(|number| {
            let value = number.to_f64_nearest();
            let stored = Rational::from_f64_exact(value).expect("finite knot");
            (value, (number - &stored).to_f64_nearest())
        })
        .unzip()
}

fn ball(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}

fn rounded(value: &Ball) -> f64 {
    value.mid().to_f64_nearest()
}

/// A unit vector as 128-bit balls, as `math::frame` builds the kernel's frames.
fn unit(v: &[Rational; 3]) -> Result<[Ball; 3], Refusal> {
    let balls = v.each_ref().map(ball);
    let length = (balls[0].square() + balls[1].square() + balls[2].square())
        .sqrt()
        .expect("sum of squares");
    let mut out = Vec::with_capacity(3);
    for component in &balls {
        out.push(
            component
                .checked_div(&length)
                .ok_or_else(|| unsupported("a frame direction has zero length"))?,
        );
    }
    Ok(out.try_into().expect("three components"))
}

/// The kernel's frame `[unit(reference), unit(axis x reference), unit(axis)]`.
fn frame(
    axis: &Vector3,
    reference: &Vector3,
    dec: &mut Decider,
) -> Result<[[Ball; 3]; 3], Refusal> {
    let vector = |v: &Vector3, dec: &mut Decider| -> Result<[Rational; 3], Refusal> {
        Ok([
            rational(&v.x, dec)?,
            rational(&v.y, dec)?,
            rational(&v.z, dec)?,
        ])
    };
    let (a, r) = (vector(axis, dec)?, vector(reference, dec)?);
    let cross = [
        &a[1] * &r[2] - &a[2] * &r[1],
        &a[2] * &r[0] - &a[0] * &r[2],
        &a[0] * &r[1] - &a[1] * &r[0],
    ];
    Ok([unit(&r)?, unit(&cross)?, unit(&a)?])
}

fn rounded_frame(frame: &[[Ball; 3]; 3]) -> [[f64; 3]; 3] {
    frame.each_ref().map(|e| e.each_ref().map(rounded))
}

fn root(square: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    ball(&rational(square, dec)?)
        .sqrt()
        .ok_or_else(|| unsupported("a squared size is negative"))
}

/// `point - origin + sum of shares[k] * frame[k]`, rounded once at the end.
fn placed(
    point: &[Rational; 3],
    origin: &[Rational; 3],
    frame: &[[Ball; 3]; 3],
    shares: [&Ball; 3],
) -> [f64; 3] {
    std::array::from_fn(|axis| {
        let mut value = ball(&(&point[axis] - &origin[axis]));
        for k in 0..3 {
            value = &value + &(&frame[k][axis] * shares[k]);
        }
        rounded(&value)
    })
}

fn spline_curve<const D: usize>(
    source: &BSplineCurve,
    shift: &[Rational],
    dec: &mut Decider,
) -> Result<Spline<D>, Refusal> {
    let count = source.control_count();
    dec.charge_work((source.knots().len() + D * count) as u64)?;
    let mut controls = Vec::with_capacity(count);
    let mut reach = Vec::with_capacity(count);
    for i in 0..count {
        let weight = rational(source.weight(i), dec)?;
        let point = rational_point(source.control(i), dec)?;
        let mut control = [0.0; D];
        let mut size = 0.0_f64;
        for axis in 0..D - 1 {
            let local = &point[axis] - &shift[axis];
            size = size.max(local.to_f64_nearest().abs());
            control[axis] = (&local * &weight).to_f64_nearest();
        }
        control[D - 1] = weight.to_f64_nearest();
        controls.push(control);
        reach.push(size);
    }
    let (knots, low) = split(source.knots().iter());
    Ok(Spline {
        degree: source.degree(),
        knots,
        low,
        controls,
        reach,
    })
}

fn net(source: &BSplineSurface, origin: &[Rational; 3], dec: &mut Decider) -> Result<Net, Refusal> {
    let count = source.control_count();
    let knots = source.knots();
    dec.charge_work((knots[0].len() + knots[1].len() + 4 * count[0] * count[1]) as u64)?;
    let mut controls = Vec::with_capacity(count[0] * count[1]);
    let mut reach = Vec::with_capacity(controls.capacity());
    for (point, weight) in source.controls().iter().zip(source.weights()) {
        let weight = rational(weight, dec)?;
        let point = rational_point(point, dec)?;
        let mut control = [0.0; 4];
        let mut size = 0.0_f64;
        for axis in 0..3 {
            let local = &point[axis] - &origin[axis];
            size = size.max(local.to_f64_nearest().abs());
            control[axis] = (&local * &weight).to_f64_nearest();
        }
        control[3] = weight.to_f64_nearest();
        controls.push(control);
        reach.push(size);
    }
    let (u, v) = (split(knots[0].iter()), split(knots[1].iter()));
    Ok(Net {
        degrees: source.degrees(),
        knots: [u.0, v.0],
        low: [u.1, v.1],
        count,
        controls,
        reach,
    })
}

/// The face copy's u/v at the first edge parameter, as exact rationals. It is
/// the u/v shift of planes and cylinders.
fn face_start(copies: &EdgeCopies<'_>, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
    Ok(match copies.face_copy {
        FaceCopyRef::Native(curve) => {
            let first = curve.control(0);
            [rational(&first.x, dec)?, rational(&first.y, dec)?]
        }
        FaceCopyRef::Linear(ends) => ends[0].clone(),
        FaceCopyRef::Projected => [Rational::zero(), Rational::zero()],
    })
}

/// Rounds both copies and the surface to binary64 about a local origin.
pub(super) fn lift(copies: &EdgeCopies<'_>, dec: &mut Decider) -> Result<Lifted, Refusal> {
    let source = match (copies.curve, copies.surface) {
        (GapCurve::Sampled(curve), GapSurface::Sampled(surface)) => {
            return lift_sampled(copies, curve, surface, dec);
        }
        // An evaluator works in its owner's frame, and a kernel carrier is
        // shifted here to a frame of its own: the two cannot be read
        // against each other.
        (GapCurve::Sampled(_), _) | (_, GapSurface::Sampled(_)) => {
            return Err(unsupported(
                "a sampled carrier is read against a sampled carrier in the same frame",
            ));
        }
        (GapCurve::Kernel(source), _) => source,
    };
    let origin = match source {
        Curve::Line(line) => rational_point(&line.origin, dec)?,
        Curve::Circle(circle) => rational_point(&circle.center, dec)?,
        Curve::Ellipse(ellipse) => rational_point(&ellipse.center, dec)?,
        Curve::BSpline(curve) => rational_point(curve.control(0), dec)?,
        _ => return Err(unsupported("the 3D copy has no parameter evaluator")),
    };
    let zero = Ball::from_rational(&Rational::zero(), BITS);
    let curve = match source {
        Curve::Line(line) => {
            let dir = &line.dir;
            Curve64::Line {
                origin: [0.0; 3],
                dir: [
                    rational(&dir.x, dec)?.to_f64_nearest(),
                    rational(&dir.y, dec)?.to_f64_nearest(),
                    rational(&dir.z, dec)?.to_f64_nearest(),
                ],
            }
        }
        Curve::Circle(circle) => {
            let frame = frame(&circle.normal, &circle.ref_dir, dec)?;
            let radius = root(&circle.radius_sq, dec)?;
            Curve64::Conic {
                center: [0.0; 3],
                major: placed(&origin, &origin, &frame, [&radius, &zero, &zero]),
                minor: placed(&origin, &origin, &frame, [&zero, &radius, &zero]),
            }
        }
        Curve::Ellipse(ellipse) => {
            let frame = frame(&ellipse.normal, &ellipse.major_dir, dec)?;
            let (major, minor) = (root(&ellipse.major_sq, dec)?, root(&ellipse.minor_sq, dec)?);
            Curve64::Conic {
                center: [0.0; 3],
                major: placed(&origin, &origin, &frame, [&major, &zero, &zero]),
                minor: placed(&origin, &origin, &frame, [&zero, &minor, &zero]),
            }
        }
        Curve::BSpline(source) => Curve64::Spline(spline_curve(source, &origin, dec)?),
        _ => unreachable!("refused above"),
    };
    let start = face_start(copies, dec)?;
    let mut shift = [Rational::zero(), Rational::zero()];
    let analytic = |kind: Revolved,
                    point: &Point3,
                    axis: &Vector3,
                    reference: &Vector3,
                    sizes: [Ball; 2],
                    height: &Rational,
                    dec: &mut Decider|
     -> Result<Surface64, Refusal> {
        let frame = frame(axis, reference, dec)?;
        let point = rational_point(point, dec)?;
        Ok(Surface64::Revolved {
            kind,
            origin: placed(&point, &origin, &frame, [&zero, &zero, &ball(height)]),
            frame: rounded_frame(&frame),
            sizes: sizes.each_ref().map(rounded),
        })
    };
    let surface = match copies.surface {
        GapSurface::Carrier(carrier) => Surface64::Spline(net(carrier, &origin, dec)?),
        GapSurface::Kernel(Surface::BSpline(support)) => {
            Surface64::Spline(net(support.carrier(), &origin, dec)?)
        }
        GapSurface::Kernel(Surface::Plane(plane)) => {
            shift = start;
            let frame = frame(&plane.normal, &plane.ref_dir, dec)?;
            let point = rational_point(&plane.origin, dec)?;
            Surface64::Plane {
                origin: placed(
                    &point,
                    &origin,
                    &frame,
                    [&ball(&shift[0]), &ball(&shift[1]), &zero],
                ),
                frame: rounded_frame(&frame),
            }
        }
        GapSurface::Kernel(Surface::Cylinder(s)) => {
            shift[1] = start[1].clone();
            let radius = root(&s.radius_sq, dec)?;
            analytic(
                Revolved::Cylinder,
                &s.origin,
                &s.axis,
                &s.ref_dir,
                [radius, zero.clone()],
                &shift[1],
                dec,
            )?
        }
        GapSurface::Kernel(Surface::EllipticCylinder(s)) => {
            shift[1] = start[1].clone();
            let sizes = [root(&s.major_sq, dec)?, root(&s.minor_sq, dec)?];
            analytic(
                Revolved::EllipticCylinder,
                &s.origin,
                &s.axis,
                &s.major_dir,
                sizes,
                &shift[1],
                dec,
            )?
        }
        GapSurface::Kernel(Surface::Cone(s)) => {
            let slope = ball(&rational(&s.slope, dec)?);
            analytic(
                Revolved::Cone,
                &s.apex,
                &s.axis,
                &s.ref_dir,
                [slope, zero.clone()],
                &Rational::zero(),
                dec,
            )?
        }
        GapSurface::Kernel(Surface::Sphere(s)) => {
            let radius = ball(&rational(&s.radius, dec)?);
            analytic(
                Revolved::Sphere,
                &s.center,
                &s.axis,
                &s.ref_dir,
                [radius, zero.clone()],
                &Rational::zero(),
                dec,
            )?
        }
        GapSurface::Kernel(Surface::Torus(s)) => {
            let sizes = [
                ball(&rational(&s.major, dec)?),
                ball(&rational(&s.minor, dec)?),
            ];
            analytic(
                Revolved::Torus,
                &s.center,
                &s.axis,
                &s.ref_dir,
                sizes,
                &Rational::zero(),
                dec,
            )?
        }
        GapSurface::Kernel(Surface::TorusSheet(s)) => {
            let sizes = [
                ball(&rational(&s.major, dec)?),
                ball(&rational(&s.minor, dec)?),
            ];
            analytic(
                Revolved::Torus,
                &s.center,
                &s.axis,
                &s.ref_dir,
                sizes,
                &Rational::zero(),
                dec,
            )?
        }
        GapSurface::Kernel(_) => return Err(unsupported("the surface has no chart evaluator")),
        GapSurface::Sampled(_) => unreachable!("refused above"),
    };
    lifted(copies, curve, surface, shift, dec)
}

/// The face copy rounded about the u/v shift, and the edge's parameters.
fn lifted(
    copies: &EdgeCopies<'_>,
    curve: Curve64,
    surface: Surface64,
    shift: [Rational; 2],
    dec: &mut Decider,
) -> Result<Lifted, Refusal> {
    let face = match copies.face_copy {
        FaceCopyRef::Native(source) => Face64::Native(spline_curve(source, &shift, dec)?),
        FaceCopyRef::Linear(ends) => {
            let [first, last] = copies.parameters;
            let (turned, run) = (&ends[1][0] - &ends[0][0], last - first);
            Face64::Linear {
                ends: ends.each_ref().map(|uv| {
                    std::array::from_fn(|axis| (&uv[axis] - &shift[axis]).to_f64_nearest())
                }),
                turns_with_edge: ends[0][1] == ends[1][1]
                    && !run.is_zero()
                    && turned.abs() == run.abs(),
            }
        }
        FaceCopyRef::Projected => {
            if !matches!(surface, Surface64::Plane { .. }) {
                return Err(unsupported("only a plane has a projected face copy"));
            }
            Face64::Projected
        }
    };
    let (parameters, low) = split(copies.parameters.into_iter());
    Ok(Lifted {
        curve,
        face,
        surface,
        parameters: [parameters[0], parameters[1]],
        low: [low[0], low[1]],
        shift,
    })
}

/// Both copies of an edge use whose 3D curve and surface are read through
/// their own evaluators. Nothing is shifted: the evaluators share a frame of
/// their owner's choice, and the face copy keeps the surface's own u/v.
fn lift_sampled(
    copies: &EdgeCopies<'_>,
    curve: &SampledCurve,
    surface: &SampledSurface,
    dec: &mut Decider,
) -> Result<Lifted, Refusal> {
    lifted(
        copies,
        Curve64::Sampled(curve.clone()),
        Surface64::Sampled(surface.clone()),
        [Rational::zero(), Rational::zero()],
        dec,
    )
}

pub(super) fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl<const D: usize> Spline<D> {
    /// The Cartesian point at `t`. Past either end the end span's polynomial
    /// continues.
    /// Source: Piegl and Tiller, algorithms A3.1 and A4.1.
    pub(super) fn point(&self, t: f64, scratch: &mut Scratch) -> Bounded<[f64; D]> {
        let span = find_span(&self.knots, self.degree, self.controls.len(), t);
        basis_functions(
            &self.knots,
            &self.low,
            self.degree,
            span,
            t,
            &mut scratch.u,
            &mut scratch.work,
        );
        let first = span - self.degree;
        let mut sum = [0.0; D];
        let (mut reach, mut spread) = (0.0_f64, 0.0);
        for (i, share) in scratch.u.iter().enumerate() {
            let control = &self.controls[first + i];
            for axis in 0..D {
                sum[axis] += share * control[axis];
            }
            spread += share.abs() * control[D - 1];
            reach = reach.max(self.reach[first + i]);
        }
        let weight = sum[D - 1];
        for value in sum.iter_mut().take(D - 1) {
            *value /= weight;
        }
        Bounded {
            value: sum,
            error: UNIT_ROUNDOFF * roundings(self.degree) * reach * spread / weight.abs(),
        }
    }
}

impl Net {
    /// The surface point at `(u, v)`, continuing the end spans past the border.
    /// Source: Piegl and Tiller, algorithms A3.5 and A4.3.
    fn point(&self, u: f64, v: f64, scratch: &mut Scratch) -> Bounded<[f64; 3]> {
        let [du, dv] = self.degrees;
        let su = find_span(&self.knots[0], du, self.count[0], u);
        let sv = find_span(&self.knots[1], dv, self.count[1], v);
        basis_functions(
            &self.knots[0],
            &self.low[0],
            du,
            su,
            u,
            &mut scratch.u,
            &mut scratch.work,
        );
        basis_functions(
            &self.knots[1],
            &self.low[1],
            dv,
            sv,
            v,
            &mut scratch.v,
            &mut scratch.work,
        );
        let mut sum = [0.0; 4];
        let (mut reach, mut spread) = (0.0_f64, 0.0);
        for (i, share_u) in scratch.u.iter().enumerate() {
            let row = (su - du + i) * self.count[1] + (sv - dv);
            for (j, share_v) in scratch.v.iter().enumerate() {
                let share = share_u * share_v;
                let control = &self.controls[row + j];
                for axis in 0..4 {
                    sum[axis] += share * control[axis];
                }
                spread += share.abs() * control[3];
                reach = reach.max(self.reach[row + j]);
            }
        }
        Bounded {
            value: [sum[0] / sum[3], sum[1] / sum[3], sum[2] / sum[3]],
            error: UNIT_ROUNDOFF * (roundings(du) + roundings(dv)) * reach * spread / sum[3].abs(),
        }
    }

    /// Whether `(u, v)` lies outside the knot rectangle by more than rounding.
    /// The homogeneous controls whose basis functions are nonzero somewhere
    /// on the u/v rectangle `rect`: the local support of a B-spline (Piegl
    /// and Tiller, "The NURBS Book", section 2.2: local support), so their hull holds the surface there.
    pub(super) fn controls_over(&self, rect: [[f64; 2]; 2]) -> Vec<[f64; 4]> {
        let spans: [[usize; 2]; 2] = std::array::from_fn(|axis| {
            let [lo, hi] = rect[axis];
            let d = self.degrees[axis];
            let (knots, count) = (&self.knots[axis], self.count[axis]);
            [
                find_span(knots, d, count, lo) - d,
                find_span(knots, d, count, hi),
            ]
        });
        let mut out = Vec::new();
        for i in spans[0][0]..=spans[0][1] {
            for j in spans[1][0]..=spans[1][1] {
                out.push(self.controls[i * self.count[1] + j]);
            }
        }
        out
    }

    pub(super) fn outside(&self, uv: [f64; 2]) -> bool {
        (0..2).any(|axis| {
            let [lo, hi] = self.domain(axis);
            let slack = OUTSIDE_ROUNDINGS * UNIT_ROUNDOFF * lo.abs().max(hi.abs()).max(hi - lo);
            uv[axis] < lo - slack || uv[axis] > hi + slack
        })
    }

    pub(super) fn domain(&self, axis: usize) -> [f64; 2] {
        [
            self.knots[axis][self.degrees[axis]],
            self.knots[axis][self.count[axis]],
        ]
    }
}

/// Source: an operation count. A face copy meant to run along a border is
/// rounded when stored and again when evaluated at each of its two
/// coordinates' sums; sixteen units cover both without calling a real
/// excursion a rounding error.
const OUTSIDE_ROUNDINGS: f64 = 16.0;

/// Source: an operation count of the analytic charts below: two products and
/// two sums per frame term, three terms, over the rounded frame and sizes.
const CHART_ROUNDINGS: f64 = 12.0;

impl Surface64 {
    /// The point at `uv`.
    pub(super) fn point(&self, uv: [f64; 2], scratch: &mut Scratch) -> Bounded<[f64; 3]> {
        match self {
            Self::Spline(net) => net.point(uv[0], uv[1], scratch),
            Self::Sampled(sampled) => {
                let (value, error) = (sampled.point)(uv);
                Bounded { value, error }
            }
            Self::Plane { origin, frame } => {
                let value = std::array::from_fn(|axis| {
                    origin[axis] + uv[0] * frame[0][axis] + uv[1] * frame[1][axis]
                });
                let reach = origin
                    .iter()
                    .fold(uv[0].abs() + uv[1].abs(), |sum, c| sum.max(c.abs()));
                Bounded {
                    value,
                    error: UNIT_ROUNDOFF * CHART_ROUNDINGS * reach,
                }
            }
            Self::Revolved {
                kind,
                origin,
                frame,
                sizes,
            } => {
                let (sin_u, cos_u) = sin_cos(uv[0]);
                let v = uv[1];
                // The kernel's charts, `Surface::evaluate`.
                let [a, b, c] = match kind {
                    Revolved::Cylinder => [sizes[0] * cos_u, sizes[0] * sin_u, v],
                    Revolved::EllipticCylinder => [sizes[0] * cos_u, sizes[1] * sin_u, v],
                    Revolved::Cone => {
                        let radial = v * sizes[0];
                        [radial * cos_u, radial * sin_u, v]
                    }
                    Revolved::Sphere => {
                        let (sin_v, cos_v) = sin_cos(v);
                        let radial = sizes[0] * cos_v;
                        [radial * cos_u, radial * sin_u, sizes[0] * sin_v]
                    }
                    Revolved::Torus => {
                        let (sin_v, cos_v) = sin_cos(v);
                        let radial = sizes[0] + sizes[1] * cos_v;
                        [radial * cos_u, radial * sin_u, sizes[1] * sin_v]
                    }
                };
                let value = std::array::from_fn(|axis| {
                    origin[axis] + a * frame[0][axis] + b * frame[1][axis] + c * frame[2][axis]
                });
                let reach = origin
                    .iter()
                    .fold(a.abs() + b.abs() + c.abs(), |sum, x| sum.max(x.abs()));
                Bounded {
                    value,
                    error: UNIT_ROUNDOFF * (CHART_ROUNDINGS + 2.0 * TRIG_ROUNDINGS) * reach,
                }
            }
        }
    }

    /// The foot on the surface of `target`, from a seed: Gauss-Newton on the
    /// squared distance with forward-difference partials, as
    /// [`Lifted::allowance`] takes them (Dennis and Schnabel 1983, sections
    /// 5.4 and 10.1). `spend` is asked before each evaluation of the surface.
    /// Where it gives `None` that evaluation is not made, nor any after it,
    /// and the search gives `None`.
    pub(super) fn foot_spending(
        &self,
        target: &[f64; 3],
        seed: [f64; 2],
        scratch: &mut Scratch,
        spend: impl FnMut() -> Option<()>,
    ) -> Option<[f64; 2]> {
        nearest_within(|uv| self.point(uv, scratch).value, target, seed, spend)
    }

    /// [`Self::foot_spending`] under an allowance: each evaluation of the
    /// surface, with its stated work, is taken from `left` before it is
    /// made. `None` when `left` runs out.
    pub(super) fn foot_within(
        &self,
        target: &[f64; 3],
        seed: [f64; 2],
        scratch: &mut Scratch,
        left: &mut u64,
    ) -> Option<[f64; 2]> {
        self.foot_spending(target, seed, scratch, || self.spend(left))
    }
}

/// The two parameters at which `map` is nearest `target`, from a seed:
/// Gauss-Newton on the squared distance with forward-difference partials, as
/// [`Lifted::allowance`] takes them (Dennis and Schnabel 1983, sections 5.4
/// and 10.1). `spend` is asked before each evaluation of `map`. Where it
/// gives `None` that evaluation is not made, nor any after it, and the search
/// gives `None`.
pub(in crate::numerical) fn nearest_within(
    mut map: impl FnMut([f64; 2]) -> [f64; 3],
    target: &[f64; 3],
    seed: [f64; 2],
    mut spend: impl FnMut() -> Option<()>,
) -> Option<[f64; 2]> {
    let mut uv = seed;
    for _ in 0..FOOT_STEPS {
        spend()?;
        let base = map(uv);
        let r: [f64; 3] = std::array::from_fn(|k| target[k] - base[k]);
        let mut partial = [[0.0; 3]; 2];
        for (axis, column) in partial.iter_mut().enumerate() {
            let step = f64::EPSILON.sqrt() * uv[axis].abs().max(1.0);
            let mut moved = uv;
            moved[axis] += step;
            spend()?;
            let there = map(moved);
            *column = std::array::from_fn(|k| (there[k] - base[k]) / step);
        }
        let (a, b, c) = (
            dot(&partial[0], &partial[0]),
            dot(&partial[0], &partial[1]),
            dot(&partial[1], &partial[1]),
        );
        let (g0, g1) = (dot(&partial[0], &r), dot(&partial[1], &r));
        let det = a * c - b * b;
        if det == 0.0 || !det.is_finite() {
            break;
        }
        let d = [(c * g0 - b * g1) / det, (a * g1 - b * g0) / det];
        uv = [uv[0] + d[0], uv[1] + d[1]];
        if d[0].abs() + d[1].abs() <= f64::EPSILON * (uv[0].abs() + uv[1].abs()).max(1.0) {
            break;
        }
    }
    Some(uv)
}

/// Takes one evaluation from what is left of an allowance.
pub(super) fn spend(left: &mut u64) -> Option<()> {
    spend_units(left, 1)
}

/// Takes `units` from what is left of an allowance, before they are spent.
fn spend_units(left: &mut u64, units: u64) -> Option<()> {
    *left = left.checked_sub(units)?;
    Some(())
}

impl Surface64 {
    /// The work a sampled surface states for one evaluation, on top of the
    /// checker's own unit. A kernel surface has none.
    pub(super) fn work(&self) -> u64 {
        match self {
            Self::Sampled(sampled) => sampled.work,
            Self::Plane { .. } | Self::Revolved { .. } | Self::Spline(_) => 0,
        }
    }

    /// Takes one evaluation of this surface, with its stated work, from what
    /// is left of an allowance, before the evaluation is made.
    pub(super) fn spend(&self, left: &mut u64) -> Option<()> {
        spend_units(left, self.work().checked_add(1)?)
    }
}

impl Curve64 {
    /// The work a sampled curve states for one evaluation, on top of the
    /// checker's own unit. A kernel curve has none.
    pub(super) fn work(&self) -> u64 {
        match self {
            Self::Sampled(sampled) => sampled.work,
            Self::Line { .. } | Self::Conic { .. } | Self::Spline(_) => 0,
        }
    }

    /// Takes one evaluation of this curve, with its stated work, from what
    /// is left of an allowance, before the evaluation is made.
    pub(super) fn spend(&self, left: &mut u64) -> Option<()> {
        spend_units(left, self.work().checked_add(1)?)
    }

    pub(super) fn point(&self, t: f64, scratch: &mut Scratch) -> Bounded<[f64; 3]> {
        match self {
            Self::Line { origin, dir } => {
                let value = std::array::from_fn(|axis| origin[axis] + t * dir[axis]);
                let reach = (0..3).fold(0.0_f64, |sum, axis| sum.max((t * dir[axis]).abs()));
                Bounded {
                    value,
                    error: UNIT_ROUNDOFF * 2.0 * reach,
                }
            }
            Self::Conic {
                center,
                major,
                minor,
            } => {
                let (sine, cosine) = sin_cos(t);
                let value = std::array::from_fn(|axis| {
                    center[axis] + cosine * major[axis] + sine * minor[axis]
                });
                let reach = (0..3).fold(0.0_f64, |sum, axis| {
                    sum.max(major[axis].abs() + minor[axis].abs())
                });
                Bounded {
                    value,
                    error: UNIT_ROUNDOFF * (CHART_ROUNDINGS + TRIG_ROUNDINGS) * reach,
                }
            }
            Self::Spline(spline) => {
                let Bounded { value, error } = spline.point(t, scratch);
                Bounded {
                    value: [value[0], value[1], value[2]],
                    error,
                }
            }
            Self::Sampled(sampled) => {
                let (value, error) = (sampled.point)(t);
                Bounded { value, error }
            }
        }
    }
}

impl Lifted {
    /// Edge parameters in ascending order.
    pub(super) fn range(&self) -> [f64; 2] {
        let [a, b] = self.parameters;
        [a.min(b), a.max(b)]
    }

    /// The work of one residual evaluation: the checker's unit, and what the
    /// two carriers state for one evaluation each. It is one unit for kernel
    /// carriers.
    pub(super) fn sample_work(&self) -> u64 {
        self.curve
            .work()
            .saturating_add(self.surface.work())
            .saturating_add(1)
    }

    /// The work of [`Lifted::allowance`], the closing probes of a search:
    /// its evaluations, and what the two carriers state for the points it
    /// reads of each.
    pub(super) fn closing_work(&self) -> u64 {
        ALLOWANCE_EVALUATIONS
            .saturating_add(ALLOWANCE_CURVE_POINTS.saturating_mul(self.curve.work()))
            .saturating_add(ALLOWANCE_SURFACE_POINTS.saturating_mul(self.surface.work()))
    }

    /// The face copy's u/v at `t`, in the shifted parameters. `None` for a
    /// projected copy, which has no u/v of its own.
    pub(super) fn face_at(&self, t: f64, scratch: &mut Scratch) -> Option<Bounded<[f64; 2]>> {
        match &self.face {
            Face64::Native(spline) => {
                let Bounded { value, error } = spline.point(t, scratch);
                Some(Bounded {
                    value: [value[0], value[1]],
                    error,
                })
            }
            Face64::Linear { ends, .. } => {
                // Differences of the exact parameters, as for knots.
                let [a, b] = self.parameters;
                let share = ((t - a) - self.low[0]) / ((b - a) + (self.low[1] - self.low[0]));
                let value: [f64; 2] = std::array::from_fn(|axis| {
                    ends[0][axis] + share * (ends[1][axis] - ends[0][axis])
                });
                let reach = ends
                    .iter()
                    .flatten()
                    .fold(0.0_f64, |sum, x| sum.max(x.abs()));
                Some(Bounded {
                    value,
                    error: UNIT_ROUNDOFF * LINEAR_ROUNDINGS * reach,
                })
            }
            Face64::Projected => None,
        }
    }

    /// The residual `C(t) - S(p(t))` with an estimate of its rounding error.
    /// The surface's sensitivity to the face copy's own rounding is left to
    /// [`Lifted::allowance`], which is evaluated once at the reported point.
    pub(super) fn residual(&self, t: f64, scratch: &mut Scratch) -> Bounded<[f64; 3]> {
        let on_curve = self.curve.point(t, scratch);
        let Some(uv) = self.face_at(t, scratch) else {
            let Surface64::Plane { origin, frame } = &self.surface else {
                unreachable!("lift admits a projected copy on a plane only")
            };
            let offset: [f64; 3] = std::array::from_fn(|axis| on_curve.value[axis] - origin[axis]);
            let height = dot(&offset, &frame[2]);
            return Bounded {
                value: frame[2].map(|n| n * height),
                error: on_curve.error * CHART_ROUNDINGS,
            };
        };
        let on_surface = self.surface.point(uv.value, scratch);
        Bounded {
            value: std::array::from_fn(|axis| on_curve.value[axis] - on_surface.value[axis]),
            error: on_curve.error + on_surface.error,
        }
    }

    /// The numerical error allowance of the gap at `t`. It has three parts.
    ///
    /// The rounding error of the two evaluations.
    ///
    /// The face copy's own rounding error carried through the surface. The
    /// surface's rate of change comes from forward differences with a step of
    /// the square root of the machine epsilon (Dennis and Schnabel, "Numerical
    /// Methods for Unconstrained Optimization and Nonlinear Equations", 1983,
    /// section 5.4).
    ///
    /// The resolution of the parameter. A maximum can sit up to half a unit
    /// in the last place of `t` away from any binary64 parameter, and the gap
    /// changes over that distance by its slope. The slope is the steeper of
    /// the two secants from `t` to `reach` on either side, where `reach` is
    /// half the width of the narrowest piece next to `t`. At a smooth maximum
    /// that overstates the slope; at a kink it is the slope.
    ///
    /// Returns the allowance and the evaluations it spent.
    pub(super) fn allowance(&self, t: f64, reach: f64, scratch: &mut Scratch) -> (f64, u64) {
        let here = self.residual(t, scratch);
        let length = |r: &[f64; 3]| dot(r, r).sqrt();
        let gap = length(&here.value);
        let slope = [t - reach, t + reach]
            .iter()
            .fold(0.0_f64, |steepest, near| {
                let there = length(&self.residual(*near, scratch).value);
                steepest.max((there - gap).abs() / reach)
            });
        let resolution = UNIT_ROUNDOFF * t.abs() * slope;
        let Some(uv) = self.face_at(t, scratch) else {
            return (here.error + resolution, ALLOWANCE_EVALUATIONS);
        };
        let base = self.surface.point(uv.value, scratch).value;
        let mut rate = 0.0_f64;
        for axis in 0..2 {
            let step = f64::EPSILON.sqrt() * uv.value[axis].abs().max(1.0);
            let mut moved = uv.value;
            moved[axis] += step;
            let there = self.surface.point(moved, scratch).value;
            let change: [f64; 3] = std::array::from_fn(|k| there[k] - base[k]);
            rate = rate.max(dot(&change, &change).sqrt() / step);
        }
        (
            here.error + rate * uv.error + resolution,
            ALLOWANCE_EVALUATIONS,
        )
    }

    /// The foot on the surface of the 3D point at `t`, in the shifted
    /// parameters, from a seed: Gauss-Newton on the squared distance with
    /// forward-difference partials, as [`Lifted::allowance`] takes them
    /// (Dennis and Schnabel 1983, sections 5.4 and 10.1). The evaluation of
    /// the 3D curve and each evaluation of the surface is taken from `left`
    /// before it is made. `None` when `left` runs out.
    pub(super) fn foot_within(
        &self,
        t: f64,
        seed: [f64; 2],
        scratch: &mut Scratch,
        left: &mut u64,
    ) -> Option<[f64; 2]> {
        self.curve.spend(left)?;
        let target = self.curve.point(t, scratch).value;
        self.surface.foot_within(&target, seed, scratch, left)
    }

    /// Whether the face copy is outside a spline surface's domain at `t`.
    pub(super) fn outside(&self, t: f64, scratch: &mut Scratch) -> bool {
        match (&self.surface, self.face_at(t, scratch)) {
            (Surface64::Spline(net), Some(uv)) => net.outside(uv.value),
            _ => false,
        }
    }
}

/// Source: Dennis and Schnabel 1983, section 10.1: Gauss-Newton converges
/// fast from a seed near a foot whose residual is small; eight steps settle a
/// seed already within the tolerance, and more cannot help a stalled one.
const FOOT_STEPS: usize = 8;

/// Source: a count of [`Lifted::allowance`]: the residual at the point and
/// on either side of it, and the surface at two shifted parameters.
pub(super) const ALLOWANCE_EVALUATIONS: u64 = 5;

/// Source: a count of [`Lifted::allowance`]: the 3D curve is read at the
/// point and on either side of it.
const ALLOWANCE_CURVE_POINTS: u64 = 3;

/// Source: a count of [`Lifted::allowance`]: the surface is read under each
/// of those three curve points, at the point again, and at its two shifted
/// parameters.
const ALLOWANCE_SURFACE_POINTS: u64 = 6;

/// Source: an operation count of the straight face copy: a difference and a
/// quotient for the share, then a difference, a product and a sum.
const LINEAR_ROUNDINGS: f64 = 5.0;

/// An edge's 3D curve rounded to binary64 in a frame that a whole body
/// shares, reached by the exact shift `origin` (rule T1: operations run in
/// local frames reached by exact rational transforms). The tolerance rules
/// read distances between entities in this frame.
pub(super) fn framed_curve(
    curve: &Curve,
    origin: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Curve64, Refusal> {
    let zero = Ball::from_rational(&Rational::zero(), BITS);
    Ok(match curve {
        Curve::Line(line) => {
            let start = rational_point(&line.origin, dec)?;
            Curve64::Line {
                origin: std::array::from_fn(|axis| (&start[axis] - &origin[axis]).to_f64_nearest()),
                dir: [
                    rational(&line.dir.x, dec)?.to_f64_nearest(),
                    rational(&line.dir.y, dec)?.to_f64_nearest(),
                    rational(&line.dir.z, dec)?.to_f64_nearest(),
                ],
            }
        }
        Curve::Circle(circle) => {
            let frame = frame(&circle.normal, &circle.ref_dir, dec)?;
            let radius = root(&circle.radius_sq, dec)?;
            let center = rational_point(&circle.center, dec)?;
            Curve64::Conic {
                center: placed(&center, origin, &frame, [&zero, &zero, &zero]),
                major: placed(&center, &center, &frame, [&radius, &zero, &zero]),
                minor: placed(&center, &center, &frame, [&zero, &radius, &zero]),
            }
        }
        Curve::Ellipse(ellipse) => {
            let frame = frame(&ellipse.normal, &ellipse.major_dir, dec)?;
            let (major, minor) = (root(&ellipse.major_sq, dec)?, root(&ellipse.minor_sq, dec)?);
            let center = rational_point(&ellipse.center, dec)?;
            Curve64::Conic {
                center: placed(&center, origin, &frame, [&zero, &zero, &zero]),
                major: placed(&center, &center, &frame, [&major, &zero, &zero]),
                minor: placed(&center, &center, &frame, [&zero, &minor, &zero]),
            }
        }
        Curve::BSpline(source) => Curve64::Spline(spline_curve(source, origin, dec)?),
        _ => return Err(unsupported("the 3D copy has no parameter evaluator")),
    })
}

/// Exact points that hold an edge's curve between two parameters, each with
/// a reach: every point of the curve there is within the reach of one of
/// them along every axis. A line's two ends; a conic's centre with its
/// largest semi-axis; a spline's controls, whose hull holds the curve
/// (Piegl and Tiller, "The NURBS Book", section 4.2; positive weights).
pub(in crate::numerical) fn defining_points(
    curve: &Curve,
    parameters: [&Rational; 2],
    dec: &mut Decider,
) -> Result<Vec<([Rational; 3], Rational)>, Refusal> {
    let zero = Rational::zero();
    Ok(match curve {
        Curve::Line(line) => {
            let start = rational_point(&line.origin, dec)?;
            let dir = [
                rational(&line.dir.x, dec)?,
                rational(&line.dir.y, dec)?,
                rational(&line.dir.z, dec)?,
            ];
            parameters
                .iter()
                .map(|t| {
                    (
                        std::array::from_fn(|axis| &start[axis] + &(*t * &dir[axis])),
                        zero.clone(),
                    )
                })
                .collect()
        }
        Curve::Circle(circle) => vec![(
            rational_point(&circle.center, dec)?,
            root(&circle.radius_sq, dec)?.hi(),
        )],
        Curve::Ellipse(ellipse) => vec![(
            rational_point(&ellipse.center, dec)?,
            root(&ellipse.major_sq, dec)?
                .hi()
                .max(root(&ellipse.minor_sq, dec)?.hi()),
        )],
        Curve::BSpline(source) => (0..source.control_count())
            .map(|i| Ok((rational_point(source.control(i), dec)?, zero.clone())))
            .collect::<Result<_, Refusal>>()?,
        _ => return Err(unsupported("the 3D copy has no parameter evaluator")),
    })
}

/// A face's surface rounded to binary64 in a frame that a whole body shares,
/// with its parameters unshifted, for the distance readings of the tolerance
/// rules. A torus sheet is a torus here.
pub(super) fn framed_surface(
    surface: &Surface,
    origin: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Surface64, Refusal> {
    let zero = Ball::from_rational(&Rational::zero(), BITS);
    let analytic = |kind: Revolved,
                    point: &Point3,
                    axis: &Vector3,
                    reference: &Vector3,
                    sizes: [Ball; 2],
                    dec: &mut Decider|
     -> Result<Surface64, Refusal> {
        let frame = frame(axis, reference, dec)?;
        let point = rational_point(point, dec)?;
        Ok(Surface64::Revolved {
            kind,
            origin: placed(&point, origin, &frame, [&zero, &zero, &zero]),
            frame: rounded_frame(&frame),
            sizes: sizes.each_ref().map(rounded),
        })
    };
    Ok(match surface {
        Surface::BSpline(support) => Surface64::Spline(net(support.carrier(), origin, dec)?),
        Surface::Plane(plane) => {
            let frame = frame(&plane.normal, &plane.ref_dir, dec)?;
            let point = rational_point(&plane.origin, dec)?;
            Surface64::Plane {
                origin: placed(&point, origin, &frame, [&zero, &zero, &zero]),
                frame: rounded_frame(&frame),
            }
        }
        Surface::Cylinder(s) => analytic(
            Revolved::Cylinder,
            &s.origin,
            &s.axis,
            &s.ref_dir,
            [root(&s.radius_sq, dec)?, zero.clone()],
            dec,
        )?,
        Surface::EllipticCylinder(s) => analytic(
            Revolved::EllipticCylinder,
            &s.origin,
            &s.axis,
            &s.major_dir,
            [root(&s.major_sq, dec)?, root(&s.minor_sq, dec)?],
            dec,
        )?,
        Surface::Cone(s) => analytic(
            Revolved::Cone,
            &s.apex,
            &s.axis,
            &s.ref_dir,
            [ball(&rational(&s.slope, dec)?), zero.clone()],
            dec,
        )?,
        Surface::Sphere(s) => analytic(
            Revolved::Sphere,
            &s.center,
            &s.axis,
            &s.ref_dir,
            [ball(&rational(&s.radius, dec)?), zero.clone()],
            dec,
        )?,
        Surface::Torus(s) => analytic(
            Revolved::Torus,
            &s.center,
            &s.axis,
            &s.ref_dir,
            [
                ball(&rational(&s.major, dec)?),
                ball(&rational(&s.minor, dec)?),
            ],
            dec,
        )?,
        Surface::TorusSheet(s) => analytic(
            Revolved::Torus,
            &s.center,
            &s.axis,
            &s.ref_dir,
            [
                ball(&rational(&s.major, dec)?),
                ball(&rational(&s.minor, dec)?),
            ],
            dec,
        )?,
        _ => return Err(unsupported("the surface has no chart evaluator")),
    })
}

/// A retained u/v curve in binary64, unshifted.
pub(super) fn uv_curve(source: &BSplineCurve, dec: &mut Decider) -> Result<Spline<3>, Refusal> {
    spline_curve(source, &[Rational::zero(), Rational::zero()], dec)
}
