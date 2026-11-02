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
}

enum Foot {
    // No foot in this cell reaches the tube: `true` when inside it.
    Clear(bool),
    Dropped,
    Open,
}

impl Local {
    fn new(law: &Law, point: &Point3, dec: &Decider) -> Result<Self, Refusal> {
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
        })
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
