//! Point queries by the law: a point's feet on the centre ellipse, its section
//! ray between the two contact rays, and its distance from the ball centre.
//!
//! A foot solves `f(t) = (p - c(t)).c'(t) = 0`. With `C`, `S` and `W` the
//! half-angle forms of `cos t`, `sin t` and one, `W^2 f` is a quartic.

use super::*;
use crate::{Point3, math};
use catallax_numbers::Ball;

const BITS: u32 = 128;
// Proposal bounds, not tolerances: running out refuses and proves nothing.
const MAX_CELLS: usize = 4096;
const NARROWEST: i64 = 44;

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("shapes.parametric-query", detail)
}

// Half-angle charts covering the turn: the first runs through `t = 0` and the
// second through `t = pi`, meeting at a tangent no simple point lands on.
pub(super) fn charts() -> [(i64, Rational); 2] {
    let edge = Rational::from_ratio(97, 89).expect("nonzero");
    [(1, edge.clone()), (-1, Rational::one() / edge)]
}

// The enclosed law, in the frame `w, s, a` about `base`.
struct Local {
    point: [Ball; 3],
    offset: Ball,
    lift: Ball,
    radius: Ball,
    first: [Ball; 3],
    wall: Ball,
    // A cone's `e`, `slope / scale` and `1 / scale`; `None` on a cylinder.
    cone: Option<[Ball; 3]>,
}

enum Foot {
    // No foot in this cell reaches the tube: `true` when inside it.
    Clear(bool),
    Dropped,
    Open,
}

impl Local {
    fn new(law: &PipeLaw, point: &Point3, dec: &Decider) -> Result<Self, Refusal> {
        let enclose = |value: &Real| dec.enclose_at(value, BITS);
        let d = point.minus(&law.base);
        Ok(Self {
            point: [
                enclose(&d.dot(&law.across))?,
                enclose(&d.dot(&law.side))?,
                enclose(&d.dot(&law.axis))?,
            ],
            offset: enclose(&law.offset)?,
            lift: enclose(&law.lift)?,
            radius: enclose(&law.radius)?,
            first: [
                enclose(&(&law.reach * &law.tilt))?,
                Ball::from_rational(&Rational::zero(), BITS),
                enclose(&(&law.reach * &law.along))?,
            ],
            wall: enclose(&law.wall)?,
            cone: if law.is_conical() {
                Some([
                    enclose(&law.eccentricity())?,
                    enclose(&(&law.slope / &law.scale))?,
                    enclose(&(&Real::one() / &law.scale))?,
                ])
            } else {
                None
            },
        })
    }
    // The centre and the wall ray at a cell of a cone, in the local frame.
    // With `A = offset w + lift a` and `B = offset s`, the centre is
    // `(cos A + sin B) / P` for `P = 1 + e cos`, and its tangent runs along
    // `-sin A + (cos + e) B`.
    fn conical(
        &self,
        [e, lean, inverse]: &[Ball; 3],
        cos: &Ball,
        sin: &Ball,
    ) -> Option<([Ball; 3], [Ball; 3], [Ball; 3])> {
        let path = &Ball::from_rational(&Rational::one(), BITS) + &(e * cos);
        let centre = [
            (&self.offset * cos).checked_div(&path)?,
            (&self.offset * sin).checked_div(&path)?,
            (&self.lift * cos).checked_div(&path)?,
        ];
        let turned = cos + e;
        let tangent = [
            -&(&self.offset * sin),
            &self.offset * &turned,
            -&(&self.lift * sin),
        ];
        let ray = [
            &(&self.wall * cos) * inverse,
            &(&self.wall * sin) * inverse,
            -&(&self.wall * lean),
        ];
        Some((centre, tangent, ray))
    }
    // Cosine and sine over a cell of tangents, in chart `sign`.
    fn turn(&self, sign: i64, tangent: &Ball) -> Option<[Ball; 2]> {
        let one = Ball::from_rational(&Rational::one(), BITS);
        let two = Ball::from_rational(&Rational::from_i64(2), BITS);
        let square = tangent.square();
        let wide = &one + &square;
        let cos = (&one - &square).checked_div(&wide)?;
        let sin = (&two * tangent).checked_div(&wide)?;
        let flip = Ball::from_rational(&Rational::from_i64(sign), BITS);
        Some([&cos * &flip, &sin * &flip])
    }
    fn classify(&self, sign: i64, tangent: &Ball) -> Result<Foot, Refusal> {
        let [cos, sin] = self
            .turn(sign, tangent)
            .ok_or_else(|| unresolved("chart angle has no enclosure"))?;
        if let Some(cone) = &self.cone {
            let (centre, along, _) = self
                .conical(cone, &cos, &sin)
                .ok_or_else(|| unresolved("cone path has no enclosure"))?;
            let gap: [Ball; 3] = std::array::from_fn(|i| &self.point[i] - &centre[i]);
            let foot = (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, i| {
                &sum + &(&gap[i] * &along[i])
            });
            if foot.sign().is_some() {
                return Ok(Foot::Dropped);
            }
            let tube = &(&(&gap[0].square() + &gap[1].square()) + &gap[2].square())
                - &self.radius.square();
            return Ok(match tube.sign() {
                Some(sign) => Foot::Clear(sign == catallax_numbers::StrictSign::Negative),
                None => Foot::Open,
            });
        }
        let [x, y, h] = &self.point;
        let slope = &(&(&(x * &self.offset) + &(h * &self.lift)) * &sin);
        let foot =
            &(&(&(y * &self.offset) * &cos) - slope) + &(&self.lift.square() * &(&sin * &cos));
        if foot.sign().is_some() {
            return Ok(Foot::Dropped);
        }
        let gap = [
            x - &(&self.offset * &cos),
            y - &(&self.offset * &sin),
            h - &(&self.lift * &cos),
        ];
        let tube =
            &(&(&gap[0].square() + &gap[1].square()) + &gap[2].square()) - &self.radius.square();
        Ok(match tube.sign() {
            Some(sign) => Foot::Clear(sign == catallax_numbers::StrictSign::Negative),
            None => Foot::Open,
        })
    }
    // Whether the ray at a cell is proved outside the cone of the contact rays.
    fn off_section(&self, sign: i64, tangent: &Ball) -> Result<bool, Refusal> {
        let [cos, sin] = self
            .turn(sign, tangent)
            .ok_or_else(|| unresolved("chart angle has no enclosure"))?;
        if let Some(cone) = &self.cone {
            let (centre, _, last) = self
                .conical(cone, &cos, &sin)
                .ok_or_else(|| unresolved("cone path has no enclosure"))?;
            let gap: [Ball; 3] = std::array::from_fn(|i| &self.point[i] - &centre[i]);
            let dot = |a: &[Ball; 3], b: &[Ball; 3]| {
                (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, i| {
                    &sum + &(&a[i] * &b[i])
                })
            };
            let meet = dot(&self.first, &last);
            let toward_first = dot(&gap, &self.first);
            let toward_last = dot(&gap, &last);
            let a = &toward_first - &(&meet * &toward_last);
            let b = &toward_last - &(&meet * &toward_first);
            return Ok(a.hi() < Rational::zero() || b.hi() < Rational::zero());
        }
        let [x, y, h] = &self.point;
        let gap = [
            x - &(&self.offset * &cos),
            y - &(&self.offset * &sin),
            h - &(&self.lift * &cos),
        ];
        let last = [&self.wall * &cos, &self.wall * &sin];
        let meet = &(&self.first[0] * &last[0]) + &(&self.first[1] * &last[1]);
        let toward_first = &(&gap[0] * &self.first[0]) + &(&gap[2] * &self.first[2]);
        let toward_last = &(&gap[0] * &last[0]) + &(&gap[1] * &last[1]);
        let a = &toward_first - &(&meet * &toward_last);
        let b = &toward_last - &(&meet * &toward_first);
        Ok(a.hi() < Rational::zero() || b.hi() < Rational::zero())
    }
}

impl RollingPipe {
    /// Exact incidence on a contact rail, read from the law.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        dec.check_deadline()?;
        math::validate_point(point, dec)?;
        let d = &self.data.law;
        if d.is_conical() {
            return self.on_conical_boundary(point, dec);
        }
        let gap = point.minus(&d.base);
        let (x, y, h) = (gap.dot(&d.across), gap.dot(&d.side), gap.dot(&d.axis));
        let radial = &x.square() + &y.square();
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        if zero(&gap.dot(&d.normal) - &(&d.reach * &d.radius), dec)?
            && zero(&radial - &d.offset.square(), dec)?
        {
            return Ok(true);
        }
        Ok(zero(&radial - &d.body.square(), dec)? && zero(&(&h * &d.body) - &(&d.lift * &x), dec)?)
    }

    // A plane rail point sits one radius off its centre along the plane's
    // ray, and its centre is on the offset cone. A wall rail point is on the
    // cone, and its centre, one radius in along the cone's normal, is in the
    // centre plane. Both cones keep the nappe through `base`.
    fn on_conical_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = &self.data.law;
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        let gap = point.minus(&d.base);
        if zero(&gap.dot(&d.normal) - &(&d.reach * &d.radius), dec)? {
            let centre = gap.minus(&d.normal.scaled(&(&d.reach * &d.radius)));
            let reach = &d.offset + &(&d.slope * &centre.dot(&d.axis));
            let (x, y) = (centre.dot(&d.across), centre.dot(&d.side));
            return Ok(math::sign(&reach, dec)? == Sign::Positive
                && zero(&(&x.square() + &y.square()) - &reach.square(), dec)?);
        }
        let reach = &d.body + &(&d.slope * &gap.dot(&d.axis));
        let (x, y) = (gap.dot(&d.across), gap.dot(&d.side));
        if math::sign(&reach, dec)? != Sign::Positive
            || !zero(&(&x.square() + &y.square()) - &reach.square(), dec)?
        {
            return Ok(false);
        }
        let ray = self.conical_ray(&gap, &reach);
        let centre = gap.minus(&ray.scaled(&(&d.wall * &d.radius)));
        zero(centre.dot(&d.normal), dec)
    }
    // The cone's outward unit normal at a point `gap` from `base` whose
    // support radius is `reach`.
    fn conical_ray(&self, gap: &Vector3, reach: &Real) -> Vector3 {
        let d = &self.data.law;
        let radial = d
            .across
            .scaled(&gap.dot(&d.across))
            .plus(&d.side.scaled(&gap.dot(&d.side)))
            .scaled(&(&Real::one() / reach));
        radial
            .minus(&d.axis.scaled(&d.slope))
            .scaled(&(&Real::one() / &d.scale))
    }

    /// The band's points whose normal runs along a coordinate axis. A band
    /// normal is a section ray, square to the path, so each axis has two feet
    /// where the path turns in that coordinate, or none needed when the path
    /// never moves in it. With the rails these bound every coordinate exactly.
    pub fn axis_extrema(&self, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
        dec.check_deadline()?;
        let d = &self.data.law;
        if d.is_conical() {
            return self.conical_extrema(dec);
        }
        let major = d.major(&d.offset);
        let first = d.normal.scaled(&d.reach);
        let mut points = Vec::new();
        for axis in 0..3 {
            let unit = |i: usize| Real::from_i64(i64::from(i == axis));
            let along = Vector3::new(unit(0), unit(1), unit(2));
            let (cos, sin) = (major.dot(&along), &d.side.dot(&along) * &d.offset);
            let length = math::exact_sqrt(&(&cos.square() + &sin.square()));
            if math::sign(&length, dec)? == Sign::Zero {
                continue;
            }
            let (cos, sin) = (&cos / &length, &sin / &length);
            for turn in [Real::from_i64(1), Real::from_i64(-1)] {
                let (cos, sin) = (&cos * &turn, &sin * &turn);
                let last = d
                    .across
                    .scaled(&cos)
                    .plus(&d.side.scaled(&sin))
                    .scaled(&d.wall);
                let centre = d
                    .base
                    .plus(&major.scaled(&cos))
                    .plus(&d.side.scaled(&(&d.offset * &sin)));
                let meet = first.dot(&last);
                for ray in [along.clone(), along.scaled(&Real::from_i64(-1))] {
                    let (toward_first, toward_last) = (ray.dot(&first), ray.dot(&last));
                    let weights = [
                        &toward_first - &(&meet * &toward_last),
                        &toward_last - &(&meet * &toward_first),
                    ];
                    if math::sign(&weights[0], dec)? != Sign::Negative
                        && math::sign(&weights[1], dec)? != Sign::Negative
                    {
                        points.push(centre.plus(&ray.scaled(&d.radius)));
                    }
                }
            }
        }
        Ok(points)
    }

    // On a cone the path's tangent is square to a direction `x` where
    // `alpha sin t - beta cos t = e beta`, for `alpha = A.x` and `beta = B.x`:
    // two angles, or none when `|e beta|` exceeds `sqrt(alpha^2 + beta^2)`.
    fn conical_extrema(&self, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
        let d = &self.data.law;
        let major = d.major(&d.offset);
        let first = d.normal.scaled(&d.reach);
        let e = d.eccentricity();
        let mut points = Vec::new();
        for axis in 0..3 {
            let unit = |i: usize| Real::from_i64(i64::from(i == axis));
            let along = Vector3::new(unit(0), unit(1), unit(2));
            let (alpha, beta) = (major.dot(&along), &d.side.dot(&along) * &d.offset);
            let size = &alpha.square() + &beta.square();
            if math::sign(&size, dec)? == Sign::Zero {
                continue;
            }
            let lean = &(&e * &beta) / &math::exact_sqrt(&size);
            let room = &(&Real::one() - &lean.square()) * &size;
            if math::sign(&room, dec)? == Sign::Negative {
                continue;
            }
            let room = math::exact_sqrt(&room);
            // `(cos, sin) = lean (-beta, alpha) / |.| + side (alpha, beta) / |.|`.
            let length = math::exact_sqrt(&size);
            for turn in [Real::from_i64(1), Real::from_i64(-1)] {
                let swing = &(&turn * &room) / &size;
                let cos = &(&-&(&lean * &beta) / &length) + &(&swing * &alpha);
                let sin = &(&(&lean * &alpha) / &length) + &(&swing * &beta);
                let path = &Real::one() + &(&e * &cos);
                let centre = d.base.plus(
                    &major
                        .scaled(&cos)
                        .plus(&d.side.scaled(&(&d.offset * &sin)))
                        .scaled(&(&Real::one() / &path)),
                );
                let last = d
                    .across
                    .scaled(&cos)
                    .plus(&d.side.scaled(&sin))
                    .minus(&d.axis.scaled(&d.slope))
                    .scaled(&(&d.wall / &d.scale));
                let meet = first.dot(&last);
                for ray in [along.clone(), along.scaled(&Real::from_i64(-1))] {
                    let (toward_first, toward_last) = (ray.dot(&first), ray.dot(&last));
                    let weights = [
                        &toward_first - &(&meet * &toward_last),
                        &toward_last - &(&meet * &toward_first),
                    ];
                    if math::sign(&weights[0], dec)? != Sign::Negative
                        && math::sign(&weights[1], dec)? != Sign::Negative
                    {
                        points.push(centre.plus(&ray.scaled(&d.radius)));
                    }
                }
            }
        }
        Ok(points)
    }

    /// Exact normal at a rail point, pointing to the ball's centre. Interior
    /// points have no exact chart parameters and refuse.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        if !self.on_boundary(point, dec)? {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "a pipe normal needs exact chart parameters, and only its rails have them",
            ));
        }
        let d = &self.data.law;
        let gap = point.minus(&d.base);
        if math::sign(&(&gap.dot(&d.normal) - &(&d.reach * &d.radius)), dec)? == Sign::Zero {
            return Ok(d.normal.scaled(&-&d.reach));
        }
        if d.is_conical() {
            let reach = &d.body + &(&d.slope * &gap.dot(&d.axis));
            return Ok(self.conical_ray(&gap, &reach).scaled(&-&d.wall));
        }
        let radial = d
            .across
            .scaled(&gap.dot(&d.across))
            .plus(&d.side.scaled(&gap.dot(&d.side)));
        Ok(radial.scaled(&-(&d.wall / &d.body)))
    }

    /// Whether a point is on the band. Positive means inside the whole tube,
    /// and a point off the band's sections reports the tube's side. A nonzero
    /// sign therefore proves the point off the band, and names no global side
    /// of the finite support. Exact zero is decided on the rails and at simple
    /// chart parameters; any other point in the enclosure's width refuses.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.on_boundary(point, dec)? {
            return Ok(Sign::Zero);
        }
        let local = Local::new(&self.data.law, point, dec)?;
        let mut inside = false;
        for (sign, edge) in charts() {
            let mut cells = vec![(-&edge, edge.clone())];
            while let Some((low, high)) = cells.pop() {
                dec.check_deadline()?;
                if cells.len() > MAX_CELLS {
                    return Err(unresolved("pipe foot search does not localise"));
                }
                let tangent =
                    Ball::from_rational(&low, BITS).hull(&Ball::from_rational(&high, BITS));
                match local.classify(sign, &tangent)? {
                    Foot::Dropped => {}
                    Foot::Clear(within) => inside |= within,
                    Foot::Open => {
                        let width = &high - &low;
                        if width > Rational::one().mul_pow2(-NARROWEST) {
                            let middle = (&low + &high).mul_pow2(-1);
                            cells.push((low, middle.clone()));
                            cells.push((middle, high));
                        } else if !local.off_section(sign, &tangent)? {
                            return Err(unresolved(
                                "a point is within the pipe's enclosure width of its band",
                            ));
                        }
                    }
                }
            }
        }
        Ok(if inside {
            Sign::Positive
        } else {
            Sign::Negative
        })
    }
}
