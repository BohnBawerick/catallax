//! Point and line queries by the law, in the frame `x`, `y`, `z` about `O`.
//!
//! With `alpha = k / sa` and `beta = k`, the centre is
//! `c(t) = (alpha cos t, beta sin t, 0)`, its tangent runs along
//! `(-alpha sin t, beta cos t, 0)`, and a point's foot on the centre ellipse
//! solves `-alpha X sin t + beta Y cos t + (alpha^2 - beta^2) sin t cos t = 0`.
//! The band keeps a section ray between `m1` and `m2`, and `t` in `[0, pi/2]`,
//! which is `tau = tan(t/2)` in `[0, 1]`.

use catallax_numbers::{RealPoly, StrictSign};

use super::*;

const BITS: u32 = 128;
// Proposal bounds, not tolerances: running out refuses and proves nothing.
const MAX_CELLS: usize = 4096;
const NARROWEST: i64 = 44;

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("shapes.parametric-query", detail)
}
fn crossing(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

enum Foot {
    // No foot in this cell reaches the tube: `true` when inside it.
    Clear(bool),
    Dropped,
    Open,
}

// The enclosed law and point in the local frame.
struct Local {
    point: [Ball; 3],
    alpha: Ball,
    beta: Ball,
    radius: Ball,
    along: Ball,
    across: Ball,
}

fn number(value: i64) -> Ball {
    Ball::from_rational(&Rational::from_i64(value), BITS)
}

impl Local {
    fn new(band: &RollingMitre, point: &Point3, dec: &Decider) -> Result<Self, Refusal> {
        let d = &band.data;
        let enclose = |value: &Real| dec.enclose_at(value, BITS);
        let gap = point.minus(&d.origin);
        Ok(Self {
            point: [
                enclose(&gap.dot(&d.frame[0]))?,
                enclose(&gap.dot(&d.frame[1]))?,
                enclose(&gap.dot(&d.frame[2]))?,
            ],
            alpha: enclose(&(&d.offset / &d.along))?,
            beta: enclose(&d.offset)?,
            radius: enclose(&d.radius)?,
            along: enclose(&d.along)?,
            across: enclose(&d.across)?,
        })
    }
    fn turn(tangent: &Ball) -> Option<[Ball; 2]> {
        let square = tangent.square();
        let wide = &number(1) + &square;
        Some([
            (&number(1) - &square).checked_div(&wide)?,
            (&number(2) * tangent).checked_div(&wide)?,
        ])
    }
    fn gap(&self, cos: &Ball, sin: &Ball) -> [Ball; 3] {
        [
            &self.point[0] - &(&self.alpha * cos),
            &self.point[1] - &(&self.beta * sin),
            self.point[2].clone(),
        ]
    }
    fn classify(&self, tangent: &Ball) -> Result<Foot, Refusal> {
        let [cos, sin] =
            Self::turn(tangent).ok_or_else(|| unresolved("chart angle has no enclosure"))?;
        let [x, y, _] = &self.point;
        let bend = &(&self.alpha.square() - &self.beta.square()) * &(&sin * &cos);
        let foot = &(&-&(&(&self.alpha * x) * &sin) + &(&(&self.beta * y) * &cos)) + &bend;
        if foot.sign().is_some() {
            return Ok(Foot::Dropped);
        }
        let gap = self.gap(&cos, &sin);
        let tube =
            &(&(&gap[0].square() + &gap[1].square()) + &gap[2].square()) - &self.radius.square();
        Ok(match tube.sign() {
            Some(sign) => Foot::Clear(sign == StrictSign::Negative),
            None => Foot::Open,
        })
    }
    // Whether the ray at a cell is proved outside the cone of the contact rays.
    fn off_section(&self, tangent: &Ball) -> Result<bool, Refusal> {
        let [cos, sin] =
            Self::turn(tangent).ok_or_else(|| unresolved("chart angle has no enclosure"))?;
        let gap = self.gap(&cos, &sin);
        let wide = &self.along * &cos;
        let tilt = &self.across * &cos;
        let rays = [[wide.clone(), sin.clone(), -&tilt], [wide, sin, tilt]];
        let dot =
            |a: &[Ball; 3], b: &[Ball; 3]| (0..3).fold(number(0), |sum, i| &sum + &(&a[i] * &b[i]));
        let meet = dot(&rays[0], &rays[1]);
        let toward = [dot(&gap, &rays[0]), dot(&gap, &rays[1])];
        let a = &toward[0] - &(&meet * &toward[1]);
        let b = &toward[1] - &(&meet * &toward[0]);
        Ok(a.hi() < Rational::zero() || b.hi() < Rational::zero())
    }
}

impl RollingMitre {
    fn coordinates(&self, point: &Point3) -> [Real; 3] {
        let gap = point.minus(&self.data.origin);
        self.data.frame.each_ref().map(|axis| gap.dot(axis))
    }
    // The contact rays `m1`, `m2` and the centre at a section, in world terms.
    fn section_at(&self, cos: &Real, sin: &Real) -> (Point3, [Vector3; 2]) {
        let d = &self.data;
        let [out, up, side] = &d.frame;
        let centre = d
            .origin
            .plus(&out.scaled(&(&(&d.offset / &d.along) * cos)))
            .plus(&up.scaled(&(&d.offset * sin)));
        let common = out.scaled(&(&d.along * cos)).plus(&up.scaled(sin));
        let tilt = side.scaled(&(&d.across * cos));
        (centre, [common.minus(&tilt), common.plus(&tilt)])
    }
    // Whether a ray from a section's centre lies between its contact rays.
    fn between(rays: &[Vector3; 2], ray: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
        // At the apex both contact rays are `y`, and the two inequalities
        // below vanish for every ray. The section there keeps only that ray.
        if math::sign(&rays[0].cross(&rays[1]).norm_squared(), dec)? == Sign::Zero {
            return Ok(
                math::sign(&ray.cross(&rays[0]).norm_squared(), dec)? == Sign::Zero
                    && math::sign(&ray.dot(&rays[0]), dec)? == Sign::Positive,
            );
        }
        let meet = rays[0].dot(&rays[1]);
        let toward = [ray.dot(&rays[0]), ray.dot(&rays[1])];
        Ok(
            math::sign(&(&toward[0] - &(&meet * &toward[1])), dec)? != Sign::Negative
                && math::sign(&(&toward[1] - &(&meet * &toward[0])), dec)? != Sign::Negative,
        )
    }
    // The section parameters of a rail point, `cos t` and `sin t`.
    fn rail_turn(&self, [x, y, _]: &[Real; 3]) -> (Real, Real) {
        let d = &self.data;
        let reach = &(&d.offset / &d.along) + &(&d.radius * &d.along);
        (x / &reach, y / &d.round)
    }

    /// Exact incidence on a contact rail or on the end section's arc.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        dec.check_deadline()?;
        math::validate_point(point, dec)?;
        Ok(self.boundary_part(point, dec)?.is_some())
    }

    // Which boundary holds a point: `Some(0)` or `Some(1)` for a rail and
    // `Some(2)` for the end section.
    fn boundary_part(&self, point: &Point3, dec: &mut Decider) -> Result<Option<usize>, Refusal> {
        let d = &self.data;
        let [x, y, z] = self.coordinates(point);
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        let reach = &(&d.offset / &d.along) + &(&d.radius * &d.along);
        if math::sign(&x, dec)? != Sign::Negative && math::sign(&y, dec)? != Sign::Negative {
            let ring = &(&(&x.square() * &d.round.square()) + &(&y.square() * &reach.square()))
                - &(&reach.square() * &d.round.square());
            if zero(ring, dec)? {
                let lean = &(&d.radius * &d.across) * &x;
                for (i, flip) in [Real::one(), Real::from_i64(-1)].iter().enumerate() {
                    if zero(&(&reach * &z) + &(flip * &lean), dec)? {
                        return Ok(Some(i));
                    }
                }
            }
        }
        if zero(y.clone(), dec)? {
            let alpha = &d.offset / &d.along;
            let shift = &x - &alpha;
            if zero(&(&shift.square() + &z.square()) - &d.radius.square(), dec)?
                && math::sign(&(&shift - &(&d.radius * &d.along)), dec)? != Sign::Negative
            {
                return Ok(Some(2));
            }
        }
        Ok(None)
    }

    /// Exact normal at a boundary point, pointing to the ball's centre.
    /// Interior points have no exact chart parameters and refuse.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let Some(part) = self.boundary_part(point, dec)? else {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "a rolling mitre normal needs exact chart parameters, and only its boundary has them",
            ));
        };
        let (cos, sin) = if part == 2 {
            (Real::one(), Real::zero())
        } else {
            self.rail_turn(&self.coordinates(point))
        };
        let (centre, _) = self.section_at(&cos, &sin);
        Ok(centre
            .minus(point)
            .scaled(&(&Real::one() / &self.data.radius)))
    }

    /// Whether a point is on the band. Positive means inside the tube at a
    /// foot within the band's sections, and a point off the band reports the
    /// tube's side there. A nonzero sign therefore proves the point off the
    /// band and names no global side of the finite support. Exact zero is
    /// decided on the boundary; any other point within the enclosure's width
    /// of the band refuses.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.on_boundary(point, dec)? {
            return Ok(Sign::Zero);
        }
        let local = Local::new(self, point, dec)?;
        let mut inside = false;
        let mut cells = vec![(Rational::zero(), Rational::one())];
        while let Some((low, high)) = cells.pop() {
            dec.check_deadline()?;
            if cells.len() > MAX_CELLS {
                return Err(unresolved("rolling mitre foot search does not localise"));
            }
            let tangent = Ball::from_rational(&low, BITS).hull(&Ball::from_rational(&high, BITS));
            match local.classify(&tangent)? {
                Foot::Dropped => {}
                Foot::Clear(within) => inside |= within,
                Foot::Open => {
                    if &high - &low > Rational::one().mul_pow2(-NARROWEST) {
                        let middle = (&low + &high).mul_pow2(-1);
                        cells.push((low, middle.clone()));
                        cells.push((middle, high));
                    } else if !local.off_section(&tangent)? {
                        return Err(unresolved(
                            "a point is within the rolling mitre's enclosure width of its band",
                        ));
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

    /// Every candidate extreme point of the band along each coordinate axis:
    /// interior points whose section holds the axis, and the extremes and
    /// ends of its rails and end section.
    pub fn axis_extrema(&self, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
        let d = &self.data;
        let alpha = &d.offset / &d.along;
        let reach = &alpha + &(&d.radius * &d.along);
        let lean = &d.radius * &d.across;
        let mut points = vec![self.apex()];
        for v in [Rational::zero(), Rational::one()] {
            points.push(self.point(&Rational::zero(), &v));
        }
        // The turn `(cos, sin)` along `(p, q)` when it lies in the quarter.
        let quarter =
            |p: &Real, q: &Real, dec: &mut Decider| -> Result<Vec<(Real, Real)>, Refusal> {
                let size = &p.square() + &q.square();
                if math::sign(&size, dec)? == Sign::Zero {
                    return Ok(Vec::new());
                }
                let length = math::exact_sqrt(&size);
                let mut kept = Vec::new();
                for flip in [Real::one(), Real::from_i64(-1)] {
                    let (cos, sin) = (&(&flip * p) / &length, &(&flip * q) / &length);
                    if math::sign(&cos, dec)? != Sign::Negative
                        && math::sign(&sin, dec)? != Sign::Negative
                    {
                        kept.push((cos, sin));
                    }
                }
                Ok(kept)
            };
        for axis in 0..3 {
            let unit = |i: usize| Real::from_i64(i64::from(i == axis));
            let along = Vector3::new(unit(0), unit(1), unit(2));
            let [ex, ey, ez] = d.frame.each_ref().map(|f| f.dot(&along));
            // Interior: the tangent is square to the axis, which then lies in
            // the section plane, and the extreme is on the ray along it.
            for (cos, sin) in quarter(&(&alpha * &ex), &(&d.offset * &ey), dec)? {
                let (centre, rays) = self.section_at(&cos, &sin);
                for ray in [along.clone(), along.scaled(&Real::from_i64(-1))] {
                    if Self::between(&rays, &ray, dec)? {
                        points.push(centre.plus(&ray.scaled(&d.radius)));
                    }
                }
            }
            // Rails `(reach cos, R sin, -+lean cos)`.
            for flip in &[Real::from_i64(-1), Real::one()] {
                let p = &(&reach * &ex) + &(&(flip * &lean) * &ez);
                for (cos, sin) in quarter(&p, &(&d.round * &ey), dec)? {
                    points.push(self.local(
                        &(&reach * &cos),
                        &(&d.round * &sin),
                        &(&(flip * &lean) * &cos),
                    ));
                }
            }
            // End section `(alpha + rho cos f, 0, rho sin f)` with `cos f >= sa`.
            let size = &ex.square() + &ez.square();
            if math::sign(&size, dec)? != Sign::Zero {
                let length = math::exact_sqrt(&size);
                for flip in [Real::one(), Real::from_i64(-1)] {
                    let (cos, sin) = (&(&flip * &ex) / &length, &(&flip * &ez) / &length);
                    if math::sign(&(&cos - &d.along), dec)? != Sign::Negative {
                        points.push(self.local(
                            &(&alpha + &(&d.radius * &cos)),
                            &Real::zero(),
                            &(&d.radius * &sin),
                        ));
                    }
                }
            }
        }
        Ok(points)
    }

    // Adds a crossing at line `parameter` in the section at `tangent` when its
    // ray lies between the contact rays. A double root there is a tangency.
    #[allow(clippy::too_many_arguments)]
    fn place(
        &self,
        origin: &Point3,
        direction: &Vector3,
        tangent: &Real,
        parameter: Real,
        double: bool,
        hits: &mut Vec<Real>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let wide = &Real::one() + &tangent.square();
        let cos = &(&Real::one() - &tangent.square()) / &wide;
        let sin = &(&Real::from_i64(2) * tangent) / &wide;
        let (centre, rays) = self.section_at(&cos, &sin);
        let point = origin.plus(&direction.scaled(&parameter));
        if !Self::between(&rays, &point.minus(&centre), dec)? {
            return Ok(());
        }
        if double {
            return Err(crossing("a line is tangent to this band"));
        }
        let at = hits
            .iter()
            .try_fold(0, |at, hit| -> Result<usize, Refusal> {
                Ok(match math::sign(&(&parameter - hit), dec)? {
                    Sign::Positive => at + 1,
                    Sign::Negative => at,
                    Sign::Zero => return Err(crossing("two band crossings share a line point")),
                })
            })?;
        hits.insert(at, parameter);
        Ok(())
    }

    /// Every point where a whole line meets the band, in line order. A
    /// tangency, or a line on the tube, refuses.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        dec.check_deadline()?;
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        let d = &self.data;
        let constant = |value: Real| RealPoly::constant(1, value);
        let tau = RealPoly::var(1, 0);
        let one = constant(Real::one());
        let square = &tau * &tau;
        let wide = &one + &square;
        let cos = &one - &square;
        let sin = tau.scale(&Real::from_i64(2));
        let a = d.frame[0].scaled(&(&d.offset / &d.along));
        let b = d.frame[1].scaled(&d.offset);
        let [a, b] = [a, b].map(|v| [v.x, v.y, v.z]);
        let gap = origin.minus(&d.origin);
        let g = [gap.x.clone(), gap.y.clone(), gap.z.clone()];
        let dv = [
            direction.x.clone(),
            direction.y.clone(),
            direction.z.clone(),
        ];
        // `W c - W O = U` and the tangent `V`, both of degree two.
        let u: [RealPoly; 3] = std::array::from_fn(|i| &cos.scale(&a[i]) + &sin.scale(&b[i]));
        let v: [RealPoly; 3] = std::array::from_fn(|i| &cos.scale(&b[i]) - &sin.scale(&a[i]));
        let dot = |x: &[RealPoly; 3], y: &[RealPoly; 3]| {
            x.iter()
                .zip(y)
                .fold(constant(Real::zero()), |sum, (x, y)| &sum + &(x * y))
        };
        // `W g - U`: the line's start from the centre, times `W`.
        let start: [RealPoly; 3] = std::array::from_fn(|i| &wide.scale(&g[i]) - &u[i]);
        let rest = -dot(&start, &v);
        let along = v
            .iter()
            .zip(&dv)
            .fold(constant(Real::zero()), |sum, (v, d)| &sum + &v.scale(d));
        let distance = start.iter().zip(&dv).fold(
            (&(&wide * &wide) * &(&along * &along)).scale(&-d.radius.square()),
            |sum, (start, value)| {
                let term = &(start * &along) + &rest.scale(value);
                &sum + &(&term * &term)
            },
        );
        let within = |roots: Vec<Real>, dec: &mut Decider| -> Result<Vec<Real>, Refusal> {
            let mut kept = Vec::new();
            for tangent in roots {
                if math::sign(&tangent, dec)? != Sign::Negative
                    && math::sign(&(&Real::one() - &tangent), dec)? != Sign::Negative
                {
                    kept.push(tangent);
                }
            }
            Ok(kept)
        };
        let mut hits: Vec<Real> = Vec::new();
        // Sections holding the line: `along` vanishes with `rest`.
        if along.total_degree(dec)?.is_none() {
            return Err(crossing("a line square to every section of this band"));
        }
        let mut holding = Vec::new();
        for tangent in along.real_roots(dec)? {
            if math::sign(&rest.eval_real(std::slice::from_ref(&tangent)), dec)? == Sign::Zero {
                holding.push(tangent);
            }
        }
        for tangent in within(holding, dec)? {
            // In that section `|W g - U + s W d|^2 = rho^2 W^2` is a quadratic in `s`.
            let at = std::slice::from_ref(&tangent);
            let w = wide.eval_real(at);
            let mut toward = Real::zero();
            let mut size = Real::zero();
            let mut length = Real::zero();
            for (start, value) in start.iter().zip(&dv) {
                let start = start.eval_real(at);
                toward = &toward + &(&start * value);
                size = &size + &start.square();
                length = &length + &value.square();
            }
            let reach = &toward.square() - &(&length * &(&size - &(&w * &d.radius).square()));
            let roots = match math::sign(&reach, dec)? {
                Sign::Negative => Vec::new(),
                Sign::Zero => vec![(-(&toward / &(&w * &length)), true)],
                Sign::Positive => {
                    let root = math::exact_sqrt(&reach);
                    vec![
                        (&-&(&toward + &root) / &(&w * &length), false),
                        (&(&root - &toward) / &(&w * &length), false),
                    ]
                }
            };
            for (parameter, double) in roots {
                self.place(
                    origin, direction, &tangent, parameter, double, &mut hits, dec,
                )?;
            }
        }
        if distance.total_degree(dec)?.is_none() {
            return Err(crossing("a line lies on this band's tube"));
        }
        let slope = distance.partial(0);
        for tangent in within(distance.real_roots(dec)?, dec)? {
            let at = std::slice::from_ref(&tangent);
            let across = along.eval_real(at);
            if math::sign(&across, dec)? == Sign::Zero {
                continue;
            }
            let parameter = &rest.eval_real(at) / &(&wide.eval_real(at) * &across);
            let double = math::sign(&slope.eval_real(at), dec)? == Sign::Zero;
            self.place(
                origin, direction, &tangent, parameter, double, &mut hits, dec,
            )?;
        }
        Ok(hits
            .iter()
            .map(|hit| origin.plus(&direction.scaled(hit)))
            .collect())
    }
}
