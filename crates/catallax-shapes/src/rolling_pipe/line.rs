//! Complete transverse crossings of a line with the band.
//!
//! On a rational law the foot equation is linear in the line parameter `s`, so
//! it gives `s` as a function of `tau`, the half-angle tangent of the section.
//! The distance equation then becomes one polynomial in `tau`, and its exact
//! roots are every crossing of the whole tube. A line inside a section plane
//! has no such function. There `tau` is a shared root of the foot equation's
//! two coefficients, and the distance equation is a quadratic in `s`. A root
//! keeps its place when its section ray lies between the contact rays. A
//! tangency on the band refuses.

use catallax_numbers::{MPoly, Outcome, UPoly};

use super::{query::charts, *};

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

// The law's frame components, when all of them are rational.
struct Frame {
    origin: [Rational; 3],
    direction: [Rational; 3],
    offset: Rational,
    lift: Rational,
    radius: Rational,
}

// The line's equations on one half-angle chart. The foot equation, cleared of
// the half-angle denominator `wide`, is `shared (rest + s wide slope) = 0`.
struct Equations {
    wide: UPoly,
    // `X wide - c` at `s = 0`: the line's start from the section centre.
    gaps: [UPoly; 3],
    // Its roots are the sections holding the whole line.
    shared: UPoly,
    rest: UPoly,
    slope: UPoly,
    // The distance equation with `s = -rest / (wide slope)` put in.
    distance: UPoly,
}

impl RollingPipe {
    /// Every point where a whole line meets the closed band, in line order, or
    /// `None` when the law has an irrational coordinate and so no rational
    /// equations. A tangency, or a contact the search cannot separate, refuses.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Point3>>, Refusal> {
        dec.check_deadline()?;
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        let d = &self.data.law;
        let gap = origin.minus(&d.base);
        let parts = [
            gap.dot(&d.across),
            gap.dot(&d.side),
            gap.dot(&d.axis),
            direction.dot(&d.across),
            direction.dot(&d.side),
            direction.dot(&d.axis),
            d.offset.clone(),
            d.lift.clone(),
            d.radius.clone(),
        ];
        let parts: Option<Vec<Rational>> = parts.iter().map(Real::as_rational).collect();
        let Some(parts) = parts else {
            return Ok(None);
        };
        let frame = Frame {
            origin: [parts[0].clone(), parts[1].clone(), parts[2].clone()],
            direction: [parts[3].clone(), parts[4].clone(), parts[5].clone()],
            offset: parts[6].clone(),
            lift: parts[7].clone(),
            radius: parts[8].clone(),
        };
        let mut hits: Vec<Real> = Vec::new();
        for (index, (sign, edge)) in charts().into_iter().enumerate() {
            let chart = frame.equations(sign)?;
            // The charts meet at one tangent, and the first owns it.
            let closed = index == 0;
            for (tangent, _) in roots_within(&chart.shared, &edge, closed, dec)? {
                for (parameter, count) in chart.in_section(&frame, &tangent, dec)? {
                    self.place(&frame, sign, &mut hits, parameter, &tangent, count, dec)?;
                }
            }
            for (tangent, count) in roots_within(&chart.distance, &edge, closed, dec)? {
                let parameter = chart.parameter(&tangent);
                self.place(&frame, sign, &mut hits, parameter, &tangent, count, dec)?;
            }
        }
        Ok(Some(
            hits.iter()
                .map(|hit| origin.plus(&direction.scaled(hit)))
                .collect(),
        ))
    }

    // Adds a root of multiplicity `count` to the line order when its section
    // ray lies on the band. Off the tube's focal set, a double root is a tangency.
    #[allow(clippy::too_many_arguments)]
    fn place(
        &self,
        frame: &Frame,
        sign: i64,
        hits: &mut Vec<Real>,
        parameter: Real,
        tangent: &Real,
        count: u32,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if !self.keeps(frame, sign, &parameter, tangent, dec)? {
            return Ok(());
        }
        if count > 1 {
            return Err(unresolved("a line is tangent to this pipe"));
        }
        let at = hits
            .iter()
            .try_fold(0, |at, hit| -> Result<usize, Refusal> {
                Ok(match math::sign(&(&parameter - hit), dec)? {
                    Sign::Positive => at + 1,
                    Sign::Negative => at,
                    Sign::Zero => return Err(unresolved("two band crossings share a line point")),
                })
            })?;
        hits.insert(at, parameter);
        Ok(())
    }

    // Whether the root's ray lies between the contact rays, rails included.
    fn keeps(
        &self,
        frame: &Frame,
        sign: i64,
        parameter: &Real,
        tangent: &Real,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let d = &self.data.law;
        let wide = &Real::one() + &tangent.square();
        let flip = Real::from_i64(sign);
        let cos = &(&(&Real::one() - &tangent.square()) / &wide) * &flip;
        let sin = &(&(&Real::from_i64(2) * tangent) / &wide) * &flip;
        let at = |i: usize| {
            &Real::from(frame.origin[i].clone())
                + &(parameter * &Real::from(frame.direction[i].clone()))
        };
        let gap = [
            &at(0) - &(&d.offset * &cos),
            &at(1) - &(&d.offset * &sin),
            &at(2) - &(&d.lift * &cos),
        ];
        let first = [&d.reach * &d.tilt, &d.reach * &d.along];
        let last = [&d.wall * &cos, &d.wall * &sin];
        let meet = &first[0] * &last[0];
        let toward_first = &(&gap[0] * &first[0]) + &(&gap[2] * &first[1]);
        let toward_last = &(&gap[0] * &last[0]) + &(&gap[1] * &last[1]);
        Ok(
            math::sign(&(&toward_first - &(&meet * &toward_last)), dec)? != Sign::Negative
                && math::sign(&(&toward_last - &(&meet * &toward_first)), dec)? != Sign::Negative,
        )
    }
}

impl Frame {
    // `sign` selects the chart through `t = pi`.
    fn equations(&self, sign: i64) -> Result<Equations, Refusal> {
        let constant = |value: &Rational| MPoly::constant(1, value.clone());
        let tau = MPoly::var(1, 0);
        let flip = Rational::from_i64(sign);
        let one = constant(&Rational::one());
        let square = &tau * &tau;
        let wide = &one + &square;
        let cos = (&one - &square).scale(&flip);
        let sin = tau.scale(&Rational::from_i64(2)).scale(&flip);
        let lean = |v: &[Rational; 3]| &(&v[0] * &self.offset) + &(&v[2] * &self.lift);
        let [x, y, h] = &self.origin;
        let foot = &(&(&cos.scale(&(y * &self.offset)) - &sin.scale(&lean(&self.origin))) * &wide)
            + &(&sin * &cos).scale(&(&self.lift * &self.lift));
        let slope =
            &cos.scale(&(&self.direction[1] * &self.offset)) - &sin.scale(&lean(&self.direction));
        let gaps = [
            &wide.scale(x) - &cos.scale(&self.offset),
            &wide.scale(y) - &sin.scale(&self.offset),
            &wide.scale(h) - &cos.scale(&self.lift),
        ];
        let (foot, slope) = (single(&foot)?, single(&slope)?);
        let shared = foot.gcd(&slope);
        let rest = quotient(&foot, &shared)?;
        let slope = quotient(&slope, &shared)?;
        let (top, under) = (double(&rest), double(&slope));
        let distance = gaps.iter().zip(&self.direction).fold(
            (&(&wide * &wide) * &(&under * &under)).scale(&-(&self.radius * &self.radius)),
            |sum, (gap, along)| {
                let term = &(gap * &under) - &top.scale(along);
                &sum + &(&term * &term)
            },
        );
        let distance = single(&distance)?;
        if distance.is_zero() {
            return Err(unresolved("a line lies on this pipe's tube"));
        }
        Ok(Equations {
            wide: single(&wide)?,
            gaps: [single(&gaps[0])?, single(&gaps[1])?, single(&gaps[2])?],
            shared,
            rest,
            slope,
            distance,
        })
    }
}

impl Equations {
    // The foot equation's answer for `s`, off the shared roots. `slope` has no
    // zero at a root of `distance`, since `rest` is prime to it.
    fn parameter(&self, tangent: &Real) -> Real {
        -(&self.rest.eval_real(tangent)
            / &(&self.slope.eval_real(tangent) * &self.wide.eval_real(tangent)))
    }

    // Where a line inside the section at `tangent` meets its circle, as roots
    // of `|gap + s wide d|^2 = r^2 wide^2`.
    fn in_section(
        &self,
        frame: &Frame,
        tangent: &Real,
        dec: &mut Decider,
    ) -> Result<Vec<(Real, u32)>, Refusal> {
        let wide = self.wide.eval_real(tangent);
        let mut along = Real::zero();
        let mut size = Real::zero();
        let mut length = Real::zero();
        for (gap, value) in self.gaps.iter().zip(&frame.direction) {
            let gap = gap.eval_real(tangent);
            let value = Real::from(value.clone());
            along = &along + &(&gap * &value);
            size = &size + &gap.square();
            length = &length + &value.square();
        }
        let radius = &wide * &Real::from(frame.radius.clone());
        let reach = &along.square() - &(&length * &(&size - &radius.square()));
        let scale = &wide * &length;
        Ok(match math::sign(&reach, dec)? {
            Sign::Negative => Vec::new(),
            Sign::Zero => vec![(-(&along / &scale), 2)],
            Sign::Positive => {
                let root = reach.sqrt();
                vec![
                    (&(&-&along - &root) / &scale, 1),
                    (&(&root - &along) / &scale, 1),
                ]
            }
        })
    }
}

fn single(poly: &MPoly) -> Result<UPoly, Refusal> {
    poly.to_upoly(0)
        .ok_or_else(|| unresolved("a pipe equation has one variable"))
}

fn double(poly: &UPoly) -> MPoly {
    let terms = poly.coefficients().iter().enumerate();
    MPoly::from_terms(
        1,
        terms
            .map(|(power, value)| (vec![power as u32], value.clone()))
            .collect(),
    )
}

fn quotient(poly: &UPoly, by: &UPoly) -> Result<UPoly, Refusal> {
    let (whole, _) = poly
        .div_rem(by)
        .ok_or_else(|| unresolved("the pipe's foot equation vanished"))?;
    Ok(whole)
}

// The distinct real roots with `|tau|` below `edge`, or at it when `closed`.
fn roots_within(
    poly: &UPoly,
    edge: &Rational,
    closed: bool,
    dec: &mut Decider,
) -> Result<Vec<(Real, u32)>, Refusal> {
    let roots = match dec.real_roots(poly) {
        Outcome::Proven(roots) => roots,
        Outcome::Refused(refusal) => return Err(refusal),
        Outcome::NotExact(_) | Outcome::ProvenCoincident(_) => {
            return Err(unresolved("a pipe equation's roots were not isolated"));
        }
    };
    let edge = Real::from(edge.clone());
    let mut within = Vec::new();
    for root in roots {
        let tangent = Real::from(root.value);
        let clear = [&edge - &tangent, &edge + &tangent];
        let low = math::sign(&clear[0], dec)?;
        let high = math::sign(&clear[1], dec)?;
        let inside = |side: Sign| side == Sign::Positive || (closed && side == Sign::Zero);
        if inside(low) && inside(high) {
            within.push((tangent, root.multiplicity));
        }
    }
    Ok(within)
}
