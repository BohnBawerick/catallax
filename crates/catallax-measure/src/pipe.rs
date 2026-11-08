//! Whole-band integrals of the exact rolling pipe.
//!
//! A band point is `P = c + r e`, where `c` is the ball centre and `e` a unit
//! ray in the normal plane of the centre path. With `e1 = reach n` and
//! `e2 = (wall m - C e1) / S`, where `C = e1 . wall m` and `S^2 = 1 - C^2`, the
//! section runs over `e = cos h e1 + sin h e2` for `h` from zero to the angle
//! `f` between the contact rays.
//!
//! With `t` the path's angle and `T` the unit tangent of the centre path, the
//! surface element is `dA = r (|c'| - r e . T') dh dt`. `T` stays square to `n`,
//! and `T' . m = -offset / |c'|` from `T . m = 0` and `c' . b = offset` for
//! `b = dm/dt`, so `e . T' = -wall offset sin h / (S |c'|)`.
//! Each integrand is a polynomial in `cos h` and `sin h` with exact primitives
//! over the section, which leaves one integral in `t` for the Taylor engine.
//! The arc length of an ellipse sits inside it, so it has no closed form. The
//! band's outward side comes from the plane that shares its rail.
//!
//! On a cone the centre is `base + (cos t A + sin t B) / P` with
//! `A = offset w + lift a`, `B = offset s` and `P = 1 + e cos t`, and the wall
//! ray is `(m - slope a) / scale`. Then `c' = (-sin t A + (cos t + e) B) / P^2`,
//! and because `slope lift = -e offset` the wall ray meets
//! `c'' = (P V' + 2 e sin t V) / P^3`, for `V = P^2 c'`, in
//! `-wall offset / (scale P)`. The section ray's other component is square to
//! the path's plane, so `e . T' = -wall offset sin h / (S scale P |c'|)`, which
//! is the cylinder's term at `e = 0` and `scale = 1`.

use catallax_numbers::{Ball, Decider, Integer, Rational, Real, Refusal};
use catallax_shapes::{RPoint3, RollingPipe, Vector3};
use catallax_topology::{FaceId, Model, pipe_outward, pipe_rows};
use std::collections::BTreeMap;

use crate::integrate::{Builder, Expr, Program, Stop, Work};
use crate::planar::{real_components as components, real_direction as direction};
use crate::surface::Property;

type Vector<'a> = [Expr<'a>; 3];

// Sixteen panels make `t` steps of a little under 0.4, where the integrand is
// well inside its radius of convergence.
const PANELS: u64 = 16;

fn unsupported(detail: &'static str) -> Refusal {
    Refusal::new("measure.unsupported-surface", detail)
}

// A polynomial in `cos h` and `sin h` with coefficients along the rail.
#[derive(Clone)]
pub(crate) struct Section<'a> {
    terms: BTreeMap<(u8, u8), Expr<'a>>,
}

impl<'a> Section<'a> {
    pub(crate) fn linear(constant: Expr<'a>, cosine: Expr<'a>, sine: Expr<'a>) -> Self {
        Self {
            terms: BTreeMap::from([((0, 0), constant), ((1, 0), cosine), ((0, 1), sine)]),
        }
    }
    pub(crate) fn times(&self, other: &Self) -> Self {
        let mut product = BTreeMap::new();
        for ((a, b), left) in &self.terms {
            for ((c, d), right) in &other.terms {
                let term = *left * *right;
                product
                    .entry((a + c, b + d))
                    .and_modify(|sum| *sum = *sum + term)
                    .or_insert(term);
            }
        }
        Self { terms: product }
    }
}

// Integrals of `cos^a h sin^b h` over the section, by the usual reductions.
pub(crate) struct Primitives<'a> {
    pub(crate) builder: &'a Builder,
    pub(crate) angle: Expr<'a>,
    pub(crate) cosine: Expr<'a>,
    pub(crate) sine: Expr<'a>,
    pub(crate) known: BTreeMap<(u8, u8), Expr<'a>>,
}

impl<'a> Primitives<'a> {
    fn power(&self, base: Expr<'a>, exponent: u8) -> Expr<'a> {
        (0..exponent).fold(self.builder.k(1), |value, _| value * base)
    }
    fn of(&mut self, a: u8, b: u8) -> Expr<'a> {
        if let Some(value) = self.known.get(&(a, b)) {
            return *value;
        }
        let k = self.builder;
        let degree = k.k(i64::from(a) + i64::from(b));
        let value = match (a, b) {
            (0, 0) => self.angle,
            (1, 0) => self.sine,
            (0, 1) => k.k(1) - self.cosine,
            (1, 1) => self.sine.square() / k.k(2),
            (a, b) if a >= 2 => {
                let end = self.power(self.cosine, a - 1) * self.power(self.sine, b + 1);
                let lower = self.of(a - 2, b);
                (end + k.k(i64::from(a) - 1) * lower) / degree
            }
            (a, b) => {
                let end = self.power(self.cosine, a + 1) * self.power(self.sine, b - 1);
                let lower = self.of(a, b - 2);
                (k.k(i64::from(b) - 1) * lower - end) / degree
            }
        };
        self.known.insert((a, b), value);
        value
    }
    pub(crate) fn integral(&mut self, section: &Section<'a>) -> Expr<'a> {
        let mut sum = self.builder.k(0);
        for ((a, b), value) in &section.terms {
            sum = sum + *value * self.of(*a, *b);
        }
        sum
    }
}

enum Parameterization {
    Turning,
    Unit,
}

/// One band's compiled integrals, in the order the engine asked for them.
pub(crate) struct Band {
    program: Program,
    parameterization: Parameterization,
    best: Option<Vec<Ball>>,
    level: u32,
}

impl Band {
    /// A band that keeps its whole ring. A trimmed one has no section running
    /// rail to rail.
    pub(crate) fn new(
        model: &Model,
        face: FaceId,
        pipe: &RollingPipe,
        anchor: &RPoint3,
        property: Property,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let rows = pipe_rows(model, face, dec)?;
        let [(rail, _), _] =
            rows.ok_or_else(|| unsupported("a trimmed pipe band has no measure"))?;
        let along_ray = pipe_outward(model, face, rail, dec)?;
        Ok(Self {
            program: forms(pipe, anchor, property, along_ray),
            parameterization: Parameterization::Turning,
            best: None,
            level: 0,
        })
    }

    /// A complete band reduced exactly to one rational parameter integral.
    pub(crate) fn parameter(program: Program) -> Self {
        Self {
            program,
            parameterization: Parameterization::Unit,
            best: None,
            level: 0,
        }
    }

    /// The band's share of one attempt. Sixteen panels already leave a band far
    /// tighter than a curved neighbour's trim, so it keeps its last answer while
    /// that stays small against the width of the faces in `rest`, and otherwise
    /// refines one level at a time.
    pub(crate) fn integrate(
        &mut self,
        rest: &[Ball],
        depth: u32,
        rectangle: bool,
        work: &mut Work,
    ) -> Result<Vec<Ball>, Stop> {
        if rectangle {
            return self.panels(depth, true, work);
        }
        if let Some(best) = &self.best {
            let small = best
                .iter()
                .zip(rest)
                .all(|(own, rest)| width(own) * Rational::from_i64(8) <= width(rest));
            if small {
                return Ok(best.clone());
            }
        }
        let level = depth.min(self.level);
        let values = self.panels(level, false, work)?;
        self.level = level + 1;
        self.best = Some(values.clone());
        Ok(values)
    }

    fn panels(&self, depth: u32, rectangle: bool, work: &mut Work) -> Result<Vec<Ball>, Stop> {
        let count = PANELS.checked_shl(depth).ok_or(Stop::Budget)?;
        let span = match self.parameterization {
            Parameterization::Turning => Rational::from_i64(2),
            Parameterization::Unit => Rational::one(),
        };
        let step = span / Rational::from_integer(Integer::from_u64(count));
        let mut sums: Option<Vec<Ball>> = None;
        for i in 0..count {
            let a = &step * &Rational::from_integer(Integer::from_u64(i));
            let b = &a + &step;
            let values = match self.parameterization {
                Parameterization::Unit => self.program.parameter_panel(&a, &b, rectangle, work)?,
                Parameterization::Turning if rectangle => {
                    self.program.turning_rectangle(&a, &b, work)?
                }
                Parameterization::Turning => self.program.turning_panel(&a, &b, work)?,
            };
            sums = Some(match sums {
                None => values,
                Some(sums) => sums
                    .iter()
                    .zip(&values)
                    .map(|(sum, value)| work.add(sum, value))
                    .collect::<Result<_, _>>()?,
            });
        }
        let sums = sums.ok_or(Stop::Budget)?;
        match self.parameterization {
            Parameterization::Unit => Ok(sums),
            Parameterization::Turning => {
                let pi = work.pi()?;
                sums.iter().map(|sum| work.mul(sum, &pi)).collect()
            }
        }
    }
}

fn width(value: &Ball) -> Rational {
    value.hi() - value.lo()
}

fn forms(pipe: &RollingPipe, anchor: &RPoint3, property: Property, along_ray: bool) -> Program {
    let k = Builder::default();
    let law = pipe.law();
    let exact = |value: &Real| k.exact(value);
    let constants = |v: &Vector3| direction(v).map(|x| k.exact(x));
    let [axis, across, side, normal] =
        [&law.axis, &law.across, &law.side, &law.normal].map(constants);
    let (radius, offset, lift) = (exact(&law.radius), exact(&law.offset), exact(&law.lift));
    let (reach, wall) = (exact(&law.reach), exact(&law.wall));
    let tilt = exact(&law.tilt);
    let base = components(&law.base).map(|x| k.exact(x));
    let anchor = components(&anchor.to_real()).map(|x| k.exact(x));
    let (cos, sin) = (k.input(0), k.input(1));

    let radial: Vector = std::array::from_fn(|i| cos * across[i] + sin * side[i]);
    let first: Vector = std::array::from_fn(|i| reach * normal[i]);
    let zero = k.k(0);
    let (centre, speed, last, c, bend) = if law.is_conical() {
        let (e, slope, scale) = (
            exact(&law.eccentricity()),
            exact(&law.slope),
            exact(&law.scale),
        );
        let along = exact(&law.along);
        let path = k.k(1) + e * cos;
        let centre: Vector = std::array::from_fn(|i| {
            base[i] - anchor[i] + (offset * radial[i] + lift * cos * axis[i]) / path
        });
        let speed = ((offset.square() + lift.square()) * sin.square()
            + (offset * (cos + e)).square())
        .sqrt()
            / path.square();
        let last: Vector = std::array::from_fn(|i| (radial[i] - slope * axis[i]) / scale);
        let c = reach * wall * (tilt * cos - slope * along) / scale;
        (centre, speed, last, c, Some(scale * path))
    } else {
        let centre: Vector = std::array::from_fn(|i| {
            base[i] - anchor[i] + offset * radial[i] + lift * cos * axis[i]
        });
        let speed = (offset.square() + (lift * sin).square()).sqrt();
        (centre, speed, radial, reach * wall * tilt * cos, None)
    };
    let s = (k.k(1) - c.square()).sqrt();
    let angle = c.angle_to(s);
    let second: Vector = std::array::from_fn(|i| (wall * last[i] - c * first[i]) / s);

    let position: [Section; 3] =
        std::array::from_fn(|i| Section::linear(centre[i], radius * first[i], radius * second[i]));
    let ray: [Section; 3] = std::array::from_fn(|i| Section::linear(zero, first[i], second[i]));
    let density = Section::linear(
        radius * speed,
        zero,
        match bend {
            Some(bend) => radius.square() * wall * offset / (s * bend * speed),
            None => radius.square() * wall * offset / (s * speed),
        },
    );
    // `P . e = r + c . e`, which keeps the flux linear in the section.
    let [towards_first, towards_second] =
        [first, second].map(|ray| (0..3).fold(zero, |sum, i| sum + centre[i] * ray[i]));
    let flux = Section::linear(radius, towards_first, towards_second);
    let carried = flux.times(&density);

    let side_sign = k.k(if along_ray { 1 } else { -1 });
    let mut primitives = Primitives {
        builder: &k,
        angle,
        cosine: c,
        sine: s,
        known: BTreeMap::new(),
    };
    let outputs = match property {
        Property::Area => vec![primitives.integral(&density)],
        Property::Volume => vec![side_sign * primitives.integral(&carried) / k.k(3)],
        Property::Centroid => {
            let volume = side_sign * primitives.integral(&carried) / k.k(3);
            let mut outputs = vec![volume];
            for i in 0..3 {
                let form = position[i]
                    .times(&position[i])
                    .times(&ray[i].times(&density));
                outputs.push(side_sign * primitives.integral(&form) / k.k(2));
            }
            outputs
        }
        Property::FaceCentroid => {
            let mut outputs = vec![primitives.integral(&density)];
            for at in &position {
                outputs.push(primitives.integral(&at.times(&density)));
            }
            outputs
        }
        Property::SecondMoment(i, j) => {
            let form = position[i].times(&position[j]).times(&carried);
            vec![side_sign * primitives.integral(&form) / k.k(5)]
        }
    };
    k.finish(&outputs)
}
