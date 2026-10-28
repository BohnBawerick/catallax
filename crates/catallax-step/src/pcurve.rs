//! Parameter-coupled inverses of rounded STEP geometry, not a second rounding
//! of the exact model. Rational frame products precede interval normalization.

use crate::{
    fit, number,
    uv::{self, Jet},
};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};

type R3 = [Rational; 3];
type B3 = [Ball; 3];

fn dot(a: &R3, b: &R3) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + &a[i] * &b[i])
}
fn cross(a: &R3, b: &R3) -> R3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
fn ball(value: &Rational, p: u32) -> Ball {
    Ball::from_rational(value, p)
}
fn scalar(value: i64, p: u32) -> Ball {
    ball(&Rational::from_i64(value), p)
}

/// The raw orthogonal basis has positive scale relative to STEP's normalized
/// AXIS2_PLACEMENT_3D. Keeping it rational preserves exact zero projections.
pub(crate) struct Placement {
    pub origin: R3,
    basis: [R3; 3],
    norms: R3,
}

impl Placement {
    pub fn new(origin: R3, axis: R3, reference: R3) -> Result<Self, Refusal> {
        let y = cross(&axis, &reference);
        let x = cross(&y, &axis);
        let basis = [x, y, axis];
        let norms = basis.each_ref().map(|v| dot(v, v));
        if norms.iter().any(Rational::is_zero) {
            return Err(number::budget(
                "rounded STEP placement has no orthogonal frame",
            ));
        }
        Ok(Self {
            origin,
            basis,
            norms,
        })
    }

    pub fn point(&self, point: &R3, p: u32) -> B3 {
        let delta = std::array::from_fn(|i| &point[i] - &self.origin[i]);
        std::array::from_fn(|i| normalized(&dot(&delta, &self.basis[i]), &self.norms[i], p))
    }

    fn vector(&self, vector: &R3, p: u32) -> B3 {
        std::array::from_fn(|i| normalized(&dot(vector, &self.basis[i]), &self.norms[i], p))
    }

    fn unit_vector(&self, vector: &R3, p: u32) -> B3 {
        let norm = dot(vector, vector);
        std::array::from_fn(|i| {
            normalized(&dot(vector, &self.basis[i]), &(&norm * &self.norms[i]), p)
        })
    }

    fn basis_vector(&self, other: &Self, index: usize, p: u32) -> B3 {
        std::array::from_fn(|i| {
            normalized(
                &dot(&other.basis[index], &self.basis[i]),
                &(&other.norms[index] * &self.norms[i]),
                p,
            )
        })
    }
}

fn normalized(numerator: &Rational, squared_norm: &Rational, p: u32) -> Ball {
    // sign(n)*sqrt(n²/d) avoids interval cancellation in a rational dot product
    // and preserves exact 0 and +/-1 for shared rounded carrier directions.
    let value = ball(&(&(numerator * numerator) / squared_norm), p)
        .sqrt()
        .expect("positive rational frame norm");
    if numerator.sign() == Sign::Negative {
        -value
    } else {
        value
    }
}

// An edge stores homogeneous polynomial spans or one analytic curve. Keep
// coefficients inline rather than allocating another box for every span.
#[expect(
    clippy::large_enum_variant,
    reason = "inline homogeneous spans avoid per-span allocations"
)]
pub(crate) enum Curve {
    Line {
        origin: R3,
        direction: R3,
    },
    Ellipse {
        placement: std::rc::Rc<Placement>,
        major: Rational,
        minor: Rational,
    },
    Polynomial {
        powers: [[Rational; 3]; fit::DEGREE + 1],
        range: [Rational; 2],
    },
}

#[expect(
    clippy::large_enum_variant,
    reason = "inline homogeneous spans avoid per-span allocations"
)]
pub(crate) enum LocalCurve {
    Line {
        origin: B3,
        direction: B3,
    },
    Ellipse {
        center: B3,
        x: B3,
        y: B3,
    },
    Polynomial {
        powers: [[Ball; 3]; fit::DEGREE + 1],
        range: [Rational; 2],
    },
}

impl Curve {
    pub fn local(&self, surface: &Placement, p: u32) -> LocalCurve {
        match self {
            Self::Line { origin, direction } => LocalCurve::Line {
                origin: surface.point(origin, p),
                direction: surface.unit_vector(direction, p),
            },
            Self::Ellipse {
                placement,
                major,
                minor,
            } => LocalCurve::Ellipse {
                center: surface.point(&placement.origin, p),
                x: surface
                    .basis_vector(placement, 0, p)
                    .map(|v| v * ball(major, p)),
                y: surface
                    .basis_vector(placement, 1, p)
                    .map(|v| v * ball(minor, p)),
            },
            Self::Polynomial { powers, range } => LocalCurve::Polynomial {
                powers: std::array::from_fn(|k| {
                    if k == 0 {
                        surface.point(&powers[k], p)
                    } else {
                        surface.vector(&powers[k], p)
                    }
                }),
                range: range.clone(),
            },
        }
    }

    pub fn analytic_range(
        &self,
        ends: [&R3; 2],
        closed: bool,
        digits: u32,
        pi: &Ball,
        dec: &Decider,
    ) -> Result<[Rational; 2], Refusal> {
        let p = pi.precision();
        let (from, to) = match self {
            Self::Line { origin, direction } => {
                let norm = dot(direction, direction);
                let parameters = ends.map(|point| {
                    let delta = std::array::from_fn(|i| &point[i] - &origin[i]);
                    dot(&delta, direction)
                });
                if parameters[0] >= parameters[1] {
                    return Err(number::budget("rounded line lost its parameter extent"));
                }
                (
                    normalized(&parameters[0], &norm, p),
                    normalized(&parameters[1], &norm, p),
                )
            }
            Self::Ellipse {
                placement,
                major,
                minor,
            } => {
                let raw = ends.map(|point| {
                    let delta = std::array::from_fn(|i| &point[i] - &placement.origin[i]);
                    [
                        dot(&delta, &placement.basis[0]),
                        dot(&delta, &placement.basis[1]),
                    ]
                });
                if raw.iter().any(|v| v.iter().all(Rational::is_zero)) {
                    return Err(number::budget("rounded ellipse endpoint is on its axis"));
                }
                let mut parameters = [scalar(0, p), scalar(0, p)];
                for i in 0..2 {
                    let x = normalized(&raw[i][0], &placement.norms[0], p)
                        .checked_div(&ball(major, p))
                        .expect("positive written ellipse major radius");
                    let y = normalized(&raw[i][1], &placement.norms[1], p)
                        .checked_div(&ball(minor, p))
                        .expect("positive written ellipse minor radius");
                    parameters[i] = uv::positive_angle(&y, &x, pi, dec)?;
                }
                let half = |v: &[Rational; 2]| {
                    if v[1].sign() == Sign::Positive
                        || (v[1].is_zero() && v[0].sign() != Sign::Negative)
                    {
                        0
                    } else {
                        1
                    }
                };
                let order = half(&raw[0]).cmp(&half(&raw[1])).then_with(|| {
                    match (&raw[0][0] * &raw[1][1] - &raw[0][1] * &raw[1][0]).sign() {
                        Sign::Positive => std::cmp::Ordering::Less,
                        Sign::Negative => std::cmp::Ordering::Greater,
                        Sign::Zero => std::cmp::Ordering::Equal,
                    }
                });
                let period = pi * &scalar(2, p);
                let to = if closed {
                    &parameters[0] + &period
                } else {
                    match order {
                        std::cmp::Ordering::Less => parameters[1].clone(),
                        std::cmp::Ordering::Greater => &parameters[1] + &period,
                        std::cmp::Ordering::Equal => {
                            return Err(number::budget("rounded ellipse lost its angular extent"));
                        }
                    }
                };
                (parameters[0].clone(), to)
            }
            Self::Polynomial { .. } => unreachable!("analytic range for a fitted curve"),
        };
        let lo = number::floor_ball(&from, digits).value;
        let hi = -number::floor_ball(&(-to), digits).value;
        if lo >= hi {
            return Err(number::budget("written pcurve parameter range collapsed"));
        }
        Ok([lo, hi])
    }
}

impl LocalCurve {
    pub fn coefficients<const N: usize>(
        &self,
        t: &Ball,
        dec: &Decider,
    ) -> Result<[Jet<N>; 3], Refusal> {
        dec.check_deadline()?;
        let result = match self {
            Self::Line { origin, direction } => {
                let t = Jet::variable(t.clone());
                std::array::from_fn(|i| {
                    Jet::constant(origin[i].clone()).add(&t.scale(&direction[i]))
                })
            }
            Self::Ellipse { center, x, y } => {
                let (sin, cos) = sin_cos(t, dec)?;
                std::array::from_fn(|i| {
                    Jet::constant(center[i].clone())
                        .add(&cos.scale(&x[i]))
                        .add(&sin.scale(&y[i]))
                })
            }
            Self::Polynomial { powers, range } => {
                let inverse = ball(&(&range[1] - &range[0]), t.precision())
                    .recip()
                    .ok_or_else(|| number::budget("written Bezier parameter span collapsed"))?;
                let parameter = Jet::variable(t.clone())
                    .sub(&Jet::constant(ball(&range[0], t.precision())))
                    .scale(&inverse);
                std::array::from_fn(|axis| {
                    let mut value = Jet::constant(powers[fit::DEGREE][axis].clone());
                    for k in (0..fit::DEGREE).rev() {
                        value = value
                            .mul(&parameter)
                            .add(&Jet::constant(powers[k][axis].clone()));
                    }
                    value
                })
            }
        };
        dec.check_deadline()?;
        Ok(result)
    }
}

fn sin_cos<const N: usize>(t: &Ball, dec: &Decider) -> Result<(Jet<N>, Jet<N>), Refusal> {
    dec.check_deadline()?;
    let sin = t.sin();
    let cos = t.cos();
    dec.check_deadline()?;
    let mut factor = Rational::one();
    let mut s = std::array::from_fn(|_| scalar(0, t.precision()));
    let mut c = std::array::from_fn(|_| scalar(0, t.precision()));
    for n in 0..N {
        if n != 0 {
            factor = factor / Rational::from_i64(n as i64);
        }
        let (sine, cosine) = match n % 4 {
            0 => (sin.clone(), cos.clone()),
            1 => (cos.clone(), -sin.clone()),
            2 => (-sin.clone(), -cos.clone()),
            _ => (-cos.clone(), sin.clone()),
        };
        let scale = ball(&factor, t.precision());
        s[n] = &sine * &scale;
        c[n] = &cosine * &scale;
    }
    Ok((Jet::from_coefficients(s), Jet::from_coefficients(c)))
}

pub(crate) enum SurfaceMap {
    Cylinder,
    Sphere,
    /// The tangent of the half-angle actually written to CONICAL_SURFACE.
    Cone {
        slope: Ball,
    },
    /// The written major radius. Callers certify that the outer equator, the
    /// written chart's v seam, is outside the trimmed domain.
    Torus {
        major: Ball,
    },
}

impl SurfaceMap {
    pub fn regular(&self, coordinates: &[Jet<1>; 3]) -> bool {
        let radial =
            coordinates[0].coefficients[0].square() + coordinates[1].coefficients[0].square();
        if radial.lo() <= Rational::zero() {
            return false;
        }
        let Self::Torus { major } = self else {
            return true;
        };
        // Points near the core circle have no separated tube angle.
        let Some(radial) = radial.sqrt() else {
            return false;
        };
        ((&radial - major).square() + coordinates[2].coefficients[0].square()).lo()
            > Rational::zero()
    }

    pub fn inverse<const N: usize>(
        &self,
        coordinates: &[Jet<N>; 3],
        u: Ball,
        pi: &Ball,
        constants: bool,
        dec: &Decider,
    ) -> Result<[[Ball; 2]; N], Refusal> {
        let [x, y, z] = coordinates;
        let u = Jet::angle(y, x, u)?;
        let v = match self {
            Self::Cylinder => z.clone(),
            Self::Sphere => {
                let radial = x.square().add(&y.square()).sqrt()?;
                let constant = if constants {
                    uv::angle(&z.coefficients[0], &radial.coefficients[0], pi, dec)?
                } else {
                    scalar(0, pi.precision())
                };
                Jet::angle(z, &radial, constant)?
            }
            Self::Cone { slope } => {
                let radial = x.square().add(&y.square()).sqrt()?;
                let inverse = (scalar(1, slope.precision()) + slope.square())
                    .recip()
                    .expect("positive cone metric");
                z.add(&radial.scale(slope)).scale(&inverse)
            }
            Self::Torus { major } => {
                let radial = x.square().add(&y.square()).sqrt()?;
                let offset = radial.sub(&Jet::constant(major.clone()));
                let constant = if constants {
                    uv::positive_angle(&z.coefficients[0], &offset.coefficients[0], pi, dec)?
                } else {
                    scalar(0, pi.precision())
                };
                Jet::angle(z, &offset, constant)?
            }
        };
        Ok(std::array::from_fn(|i| {
            [u.coefficients[i].clone(), v.coefficients[i].clone()]
        }))
    }
}

/// Select one continuous integer lift using the exact source chart, rather than
/// letting decimal perturbations switch a seam-adjacent curve by a whole turn.
pub(crate) fn lift(
    y: &Ball,
    x: &Ball,
    target: &Ball,
    pi: &Ball,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    let base = uv::angle(y, x, pi, dec)?;
    for integer in [-1, 0, 1, 2] {
        let value = &base + &(pi * &scalar(integer * 2, pi.precision()));
        if number::magnitude(&(&value - target)) < pi.lo() {
            return Ok(value);
        }
    }
    Err(number::budget(
        "written pcurve has no certified continuous angular lift",
    ))
}

pub(crate) enum SurfaceParameters {
    Cylinder { radius: Rational },
    Sphere { radius: Rational },
    Cone { angle: Rational },
    Torus { major: Rational, minor: Rational },
}

pub(crate) struct Surface {
    pub placement: std::rc::Rc<Placement>,
    pub parameters: SurfaceParameters,
}

pub(crate) struct SourceFrame {
    origin: B3,
    basis: [B3; 3],
}

fn enclose_point(point: &catallax_shapes::Point3, p: u32, dec: &Decider) -> Result<B3, Refusal> {
    Ok([
        dec.enclose_at(&point.x, p)?,
        dec.enclose_at(&point.y, p)?,
        dec.enclose_at(&point.z, p)?,
    ])
}

fn enclose_vector(vector: &catallax_shapes::Vector3, p: u32, dec: &Decider) -> Result<B3, Refusal> {
    Ok([
        dec.enclose_at(&vector.x, p)?,
        dec.enclose_at(&vector.y, p)?,
        dec.enclose_at(&vector.z, p)?,
    ])
}

fn ball_dot(a: &B3, b: &B3) -> Ball {
    (0..3).fold(scalar(0, a[0].precision()), |sum, i| sum + &a[i] * &b[i])
}

fn ball_unit(vector: B3) -> Result<B3, Refusal> {
    let norm = vector
        .iter()
        .fold(scalar(0, vector[0].precision()), |sum, c| sum + c.square())
        .sqrt()
        .and_then(|norm| norm.recip())
        .ok_or_else(|| number::budget("source chart frame is not separated"))?;
    Ok(vector.map(|c| &c * &norm))
}

impl SourceFrame {
    fn new(
        origin: &catallax_shapes::Point3,
        axis: &catallax_shapes::Vector3,
        reference: &catallax_shapes::Vector3,
        p: u32,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let z = ball_unit(enclose_vector(axis, p, dec)?)?;
        let x = ball_unit(enclose_vector(reference, p, dec)?)?;
        let y = [
            &z[1] * &x[2] - &z[2] * &x[1],
            &z[2] * &x[0] - &z[0] * &x[2],
            &z[0] * &x[1] - &z[1] * &x[0],
        ];
        Ok(Self {
            origin: enclose_point(origin, p, dec)?,
            basis: [x, y, z],
        })
    }

    fn xy(&self, point: &B3) -> [Ball; 2] {
        let delta = std::array::from_fn(|i| &point[i] - &self.origin[i]);
        [
            ball_dot(&delta, &self.basis[0]),
            ball_dot(&delta, &self.basis[1]),
        ]
    }

    fn xy_vector(&self, vector: &B3) -> [Ball; 2] {
        [
            ball_dot(vector, &self.basis[0]),
            ball_dot(vector, &self.basis[1]),
        ]
    }
}

pub(crate) struct Projection {
    pub written: std::rc::Rc<Surface>,
    pub map: SurfaceMap,
    pub pi: Ball,
    source: SourceFrame,
    digits: u32,
}

impl Projection {
    pub fn new(
        written: std::rc::Rc<Surface>,
        exact: &catallax_shapes::Surface,
        precision: u32,
        digits: u32,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let (origin, axis, reference) = match exact {
            catallax_shapes::Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
            catallax_shapes::Surface::Sphere(s) => (&s.center, &s.axis, &s.ref_dir),
            catallax_shapes::Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
            catallax_shapes::Surface::Torus(s) => (&s.center, &s.axis, &s.ref_dir),
            _ => unreachable!("periodic projection on another carrier"),
        };
        let source = SourceFrame::new(origin, axis, reference, precision, dec)?;
        let map = match &written.parameters {
            SurfaceParameters::Cylinder { .. } => SurfaceMap::Cylinder,
            SurfaceParameters::Sphere { .. } => SurfaceMap::Sphere,
            SurfaceParameters::Cone { angle } => {
                let angle = ball(angle, precision);
                let sin = angle.sin();
                let cos = angle.cos();
                dec.check_deadline()?;
                SurfaceMap::Cone {
                    slope: sin
                        .checked_div(&cos)
                        .ok_or_else(|| number::budget("written cone angle is not acute"))?,
                }
            }
            SurfaceParameters::Torus { major, .. } => SurfaceMap::Torus {
                major: ball(major, precision),
            },
        };
        Ok(Self {
            written,
            map,
            pi: Ball::pi_checked(precision, dec)?,
            source,
            digits,
        })
    }

    pub fn metric(&self, bounds: &[[Rational; 2]; 3], fitted: Option<&uv::Fit>) -> Rational {
        match &self.written.parameters {
            SurfaceParameters::Cylinder { radius } => radius.clone().max(Rational::one()),
            SurfaceParameters::Sphere { radius } => radius.clone(),
            // Both partial derivatives are bounded by the outer equator radius.
            SurfaceParameters::Torus { major, minor } => major + minor,
            SurfaceParameters::Cone { .. } => {
                let SurfaceMap::Cone { slope } = &self.map else {
                    unreachable!()
                };
                let magnitude = number::magnitude(slope);
                let extent = (0..3).fold(Rational::zero(), |sum, i| {
                    sum + (&bounds[i][0] - &self.written.placement.origin[i])
                        .abs()
                        .max((&bounds[i][1] - &self.written.placement.origin[i]).abs())
                });
                let mut v = extent * (Rational::one() + magnitude.clone());
                if let Some(fitted) = fitted {
                    for point in &fitted.points {
                        v = v.max(point[1].value.abs());
                    }
                }
                let slant = (scalar(1, slope.precision()) + slope.square())
                    .sqrt()
                    .expect("positive cone metric")
                    .hi();
                (magnitude * v).max(slant)
            }
        }
    }

    pub fn vertex(
        &self,
        exact: &catallax_shapes::Point3,
        written: &R3,
        side: crate::wire::Side,
        dec: &Decider,
    ) -> Result<[Rational; 2], Refusal> {
        let p = self.pi.precision();
        let target = match side {
            crate::wire::Side::Low => scalar(0, p),
            crate::wire::Side::High => &self.pi * &scalar(2, p),
            crate::wire::Side::Interior => {
                let [x, y] = self.source.xy(&enclose_point(exact, p, dec)?);
                uv::positive_angle(&y, &x, &self.pi, dec)?
            }
        };
        let xyz = self.written.placement.point(written, p);
        let u = lift(&xyz[1], &xyz[0], &target, &self.pi, dec)?;
        let xyz = xyz.map(Jet::<1>::constant);
        let mut uv = self.map.inverse(&xyz, u, &self.pi, true, dec)?[0].clone();
        // A vertex on a chart cut takes the cut's own parameter, so every
        // pcurve along that side writes one u value. The raw translator pairs
        // seam pcurves by the exact direction of their end-to-end chord. The
        // fit bounds the change against each curve's ideal endpoint.
        if side != crate::wire::Side::Interior {
            uv[0] = target;
        }
        Ok(uv
            .each_ref()
            .map(|value| number::ball_number(value, self.digits).value))
    }

    pub fn analytic_guide<'a>(
        &self,
        form: &catallax_intersect::CurveForm,
        ends: [&catallax_shapes::Point3; 2],
        range: &[Rational; 2],
        dec: &Decider,
    ) -> Result<Guide<'a>, Refusal> {
        use catallax_intersect::CurveForm;
        let p = self.pi.precision();
        if matches!(form, CurveForm::Line { .. }) {
            return Ok(Guide::Line {
                ends: [
                    self.source.xy(&enclose_point(ends[0], p, dec)?),
                    self.source.xy(&enclose_point(ends[1], p, dec)?),
                ],
            });
        }
        let (center, normal, reference, major, minor) = match form {
            CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            } => (center, normal, ref_dir, radius_sq, radius_sq),
            CurveForm::Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            } => (center, normal, major_dir, major_sq, minor_sq),
            _ => unreachable!("analytic periodic trim is not a line or ellipse"),
        };
        let frame = SourceFrame::new(center, normal, reference, p, dec)?;
        let major = dec.enclose_at(&number::sqrt_real(major), p)?;
        let minor = dec.enclose_at(&number::sqrt_real(minor), p)?;
        let mut parameters = [scalar(0, p), scalar(0, p)];
        for i in 0..2 {
            let [x, y] = frame.xy(&enclose_point(ends[i], p, dec)?);
            let x = x
                .checked_div(&major)
                .ok_or_else(|| number::budget("source ellipse major radius is not separated"))?;
            let y = y
                .checked_div(&minor)
                .ok_or_else(|| number::budget("source ellipse minor radius is not separated"))?;
            parameters[i] = lift(&y, &x, &ball(&range[i], p), &self.pi, dec)?;
        }
        Ok(Guide::Ellipse {
            center: self.source.xy(&frame.origin),
            x: self.source.xy_vector(&frame.basis[0]).map(|c| &c * &major),
            y: self.source.xy_vector(&frame.basis[1]).map(|c| &c * &minor),
            parameters,
        })
    }

    pub fn trace_guide<'a>(
        &'a self,
        piece: &'a catallax_intersect::TracePiece,
        span: &'a fit::TraceSpan,
    ) -> Guide<'a> {
        Guide::Trace {
            piece,
            span,
            frame: &self.source,
        }
    }
}

#[expect(
    clippy::large_enum_variant,
    reason = "borrowed trace guides and one analytic guide need no extra allocation"
)]
pub(crate) enum Guide<'a> {
    Fixed(Ball),
    Line {
        ends: [[Ball; 2]; 2],
    },
    Ellipse {
        center: [Ball; 2],
        x: [Ball; 2],
        y: [Ball; 2],
        parameters: [Ball; 2],
    },
    Trace {
        piece: &'a catallax_intersect::TracePiece,
        span: &'a fit::TraceSpan,
        frame: &'a SourceFrame,
    },
}

impl Guide<'_> {
    fn angle(
        &self,
        t: &Ball,
        range: &[Rational; 2],
        pi: &Ball,
        dec: &Decider,
    ) -> Result<Ball, Refusal> {
        if let Self::Fixed(value) = self {
            return Ok(value.clone());
        }
        let fraction = (t - &ball(&range[0], t.precision()))
            .checked_div(&ball(&(&range[1] - &range[0]), t.precision()))
            .ok_or_else(|| number::budget("native pcurve parameter interval collapsed"))?;
        let [x, y] = match self {
            Self::Fixed(_) => unreachable!(),
            Self::Line { ends } => {
                std::array::from_fn(|i| &ends[0][i] + &(&(&ends[1][i] - &ends[0][i]) * &fraction))
            }
            Self::Ellipse {
                center,
                x,
                y,
                parameters,
            } => {
                let parameter = &parameters[0] + &(&(&parameters[1] - &parameters[0]) * &fraction);
                let sin = parameter.sin();
                let cos = parameter.cos();
                dec.check_deadline()?;
                std::array::from_fn(|i| &(&center[i] + &(&x[i] * &cos)) + &(&y[i] * &sin))
            }
            Self::Trace { piece, span, frame } => {
                let parameter = &ball(&span.from, t.precision())
                    + &(&ball(&(&span.to - &span.from), t.precision()) * &fraction);
                frame.xy(&piece.evaluate(&parameter, dec)?)
            }
        };
        uv::positive_angle(&y, &x, pi, dec)
    }
}

pub(crate) struct Piece<'a> {
    pub curve: LocalCurve,
    pub guide: Guide<'a>,
    pub projection: &'a Projection,
    pub range: [Rational; 2],
    pub endpoints: [Option<[Rational; 2]>; 2],
}

impl uv::Parametric for Piece<'_> {
    fn range(&self) -> [Rational; 2] {
        self.range.clone()
    }

    fn regular(&self, parameter: &Ball, dec: &Decider) -> Result<bool, Refusal> {
        Ok(self
            .projection
            .map
            .regular(&self.curve.coefficients::<1>(parameter, dec)?))
    }

    fn coefficients<const N: usize>(
        &self,
        parameter: &Ball,
        constants: bool,
        dec: &Decider,
    ) -> Result<[[Ball; 2]; N], Refusal> {
        let xyz = self.curve.coefficients::<N>(parameter, dec)?;
        let u = if constants {
            let target = self
                .guide
                .angle(parameter, &self.range, &self.projection.pi, dec)?;
            lift(
                &xyz[1].coefficients[0],
                &xyz[0].coefficients[0],
                &target,
                &self.projection.pi,
                dec,
            )?
        } else {
            scalar(0, parameter.precision())
        };
        self.projection
            .map
            .inverse(&xyz, u, &self.projection.pi, constants, dec)
    }

    fn value(&self, parameter: &Rational, p: u32, dec: &Decider) -> Result<[Ball; 2], Refusal> {
        for i in 0..2 {
            if parameter == &self.range[i]
                && let Some(point) = &self.endpoints[i]
            {
                return Ok(point.each_ref().map(|c| ball(c, p)));
            }
        }
        Ok(self.coefficients::<1>(&ball(parameter, p), true, dec)?[0].clone())
    }
}
