//! Interval section law of the authored motion. This module reads no fitted geometry.

use catallax_features::general_sweep::{SweepPath, SweepRecipe, source_charts};
use catallax_numbers::{Ball, Decider, Interval as I, Rational, Refusal};
use catallax_recipe::{ElementId, Point3, SweepOrientation};
use catallax_shapes::{RPoint3, RVector3, Rotation};
use std::ops::{Add, Mul, Sub};

type Vec3 = [I; 3];
type JVec = [Jet; 3];

// Match the fitter's source enclosure precision before outward conversion to binary64.
const ENCLOSURE_BITS: u32 = 128;

fn number(n: i64) -> I {
    I::from_rational(&Rational::from_i64(n))
}

fn refusal(code: &'static str, detail: &'static str) -> Refusal {
    Refusal::new(code, detail)
}

fn unseparated() -> Refusal {
    refusal(
        "field.sweep-section-separation",
        "the authored sections lack a complete separation certificate",
    )
}

fn ball_interval(b: &Ball) -> I {
    I::from_rational(&b.lo()).hull(&I::from_rational(&b.hi()))
}

#[derive(Clone, Copy)]
struct Jet {
    value: I,
    first: I,
    second: I,
}

impl Jet {
    fn constant(value: I) -> Self {
        Self {
            value,
            first: number(0),
            second: number(0),
        }
    }

    fn variable(value: I) -> Self {
        Self {
            value,
            first: number(1),
            second: number(0),
        }
    }

    fn square(self) -> Self {
        Self {
            value: self.value.square(),
            first: number(2) * self.value * self.first,
            second: number(2) * (self.first.square() + self.value * self.second),
        }
    }

    fn reciprocal(self) -> Option<Self> {
        if self.value.contains_zero() {
            return None;
        }
        let value = number(1) / self.value;
        let first = -self.first * value.square();
        let second = (number(2) * self.first.square() * value - self.second) * value.square();
        Some(Self {
            value,
            first,
            second,
        })
    }

    fn sqrt(self) -> Option<Self> {
        if self.value.lo() <= 0.0 {
            return None;
        }
        let value = self.value.sqrt()?;
        let first = self.first / (number(2) * value);
        Some(Self {
            value,
            first,
            second: (self.second - number(2) * first.square()) / (number(2) * value),
        })
    }

    fn sin_cos(self, dec: &mut Decider) -> Result<(Self, Self), Refusal> {
        // Use the kernel's outward trigonometric series, never platform libm.
        let finite = |x| Rational::from_f64_exact(x).ok_or_else(unseparated);
        let lo = Ball::from_rational(&finite(self.value.lo())?, ENCLOSURE_BITS);
        let hi = Ball::from_rational(&finite(self.value.hi())?, ENCLOSURE_BITS);
        let argument = lo.hull(&hi);
        dec.charge_work(1)?;
        let s = ball_interval(&argument.sin());
        let c = ball_interval(&argument.cos());
        Ok((
            Self {
                value: s,
                first: c * self.first,
                second: c * self.second - s * self.first.square(),
            },
            Self {
                value: c,
                first: -s * self.first,
                second: -s * self.second - c * self.first.square(),
            },
        ))
    }
}

impl Add for Jet {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            first: self.first + rhs.first,
            second: self.second + rhs.second,
        }
    }
}

impl Sub for Jet {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
            first: self.first - rhs.first,
            second: self.second - rhs.second,
        }
    }
}

impl Mul for Jet {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            first: self.first * rhs.value + self.value * rhs.first,
            second: self.second * rhs.value
                + self.value * rhs.second
                + number(2) * self.first * rhs.first,
        }
    }
}

fn add(a: JVec, b: JVec) -> JVec {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: JVec, t: Jet) -> JVec {
    a.map(|x| x * t)
}
fn cross(a: JVec, b: JVec) -> JVec {
    std::array::from_fn(|i| {
        let [j, k] = [(i + 1) % 3, (i + 2) % 3];
        a[j] * b[k] - a[k] * b[j]
    })
}
fn dot(a: Vec3, b: Vec3) -> I {
    (0..3).fold(number(0), |s, i| s + a[i] * b[i])
}
fn norm(a: Vec3) -> I {
    a.into_iter()
        .fold(number(0), |s, x| s + x.square())
        .sqrt()
        .expect("squares")
}
fn constants(v: Vec3) -> JVec {
    v.map(Jet::constant)
}
fn unit(a: JVec) -> Option<JVec> {
    let squared = a
        .into_iter()
        .fold(Jet::constant(number(0)), |s, x| s + x.square());
    Some(scale(a, squared.sqrt()?.reciprocal()?))
}
fn xyz(p: &Point3) -> RPoint3 {
    let [x, y, z] = p.0.each_ref().map(|n| n.value().clone());
    RPoint3::new(x, y, z)
}
fn vector(v: &RVector3) -> Vec3 {
    [&v.x, &v.y, &v.z].map(I::from_rational)
}
fn point(v: &RPoint3) -> Vec3 {
    [&v.x, &v.y, &v.z].map(I::from_rational)
}

#[derive(Clone, Copy)]
enum Spine {
    Line {
        delta: Vec3,
        tangent: Vec3,
        guide: [Vec3; 2],
    },
    Helix {
        centre: Vec3,
        radial: Vec3,
        transverse: Vec3,
        axis: Vec3,
        radius: I,
        theta: I,
        rise: I,
        speed: I,
        spin: I,
    },
}

impl Spine {
    fn frame(
        self,
        t: Jet,
        start: Vec3,
        twist: I,
        dec: &mut Decider,
    ) -> Result<(JVec, [JVec; 3]), Refusal> {
        let (p, radial, tangent, spin) = match self {
            Self::Line {
                delta,
                tangent,
                guide,
            } => {
                let relative = add(constants(guide[0]), scale(constants(guide[1]), t));
                let radial = unit(relative).ok_or_else(|| {
                    refusal(
                        "field.sweep-guide-enclosure",
                        "the transverse guide interval includes zero",
                    )
                })?;
                (
                    add(constants(start), scale(constants(delta), t)),
                    radial,
                    constants(tangent),
                    twist,
                )
            }
            Self::Helix {
                centre,
                radial,
                transverse,
                axis,
                radius,
                theta,
                rise,
                speed,
                spin,
            } => {
                let (s, c) = (t * Jet::constant(theta)).sin_cos(dec)?;
                let ring = add(scale(constants(radial), c), scale(constants(transverse), s));
                let dz = scale(constants(axis), Jet::constant(rise));
                let p = add(constants(centre), add(ring, scale(dz, t)));
                let d = add(
                    scale(cross(constants(axis), ring), Jet::constant(theta)),
                    dz,
                );
                (
                    p,
                    scale(ring, Jet::constant(number(1) / radius)),
                    scale(d, Jet::constant(number(1) / speed)),
                    spin + twist,
                )
            }
        };
        // Bishop's parallel transport on a helix is the radial/binormal
        // frame with the integrated torsion removed. The authored twist
        // adds a rotation about the same transported tangent.
        let (s, c) = (t * Jet::constant(spin)).sin_cos(dec)?;
        let x = add(scale(radial, c), scale(cross(tangent, radial), s));
        Ok((p, [x, cross(tangent, x), tangent]))
    }
}

/// Authored section chart measured in millimetres, with outward control enclosures.
pub struct Chart {
    /// Authored profile segment.
    pub segment: ElementId,
    /// Authored section parameter interval.
    pub domain: [Rational; 2],
    /// Euclidean controls in the original sketch placement.
    pub controls: Vec<[I; 3]>,
    /// Positive homogeneous weights of those controls.
    pub weights: Vec<I>,
}

/// Certified section evaluator for a single admitted authored motion.
pub struct Motion {
    spine: Spine,
    charts: Vec<Chart>,
    start: Vec3,
    first: [Vec3; 3],
    normal: Vec3,
    twist: I,
}

/// Outward enclosures of one interval of the authored section law.
pub struct SectionBox {
    origin: JVec,
    axes: [JVec; 3],
    start: Vec3,
    first: [Vec3; 3],
    normal: Vec3,
}

impl SectionBox {
    fn rotate(&self, vector: Vec3) -> JVec {
        self.axes
            .iter()
            .zip(self.first)
            .fold(constants([number(0); 3]), |out, (&axis, initial)| {
                add(out, scale(axis, Jet::constant(dot(vector, initial))))
            })
    }

    /// Transported world control and second derivative with respect to the
    /// authored spine parameter, enclosed throughout the requested interval.
    pub fn control(&self, world: [I; 3]) -> ([I; 3], [I; 3]) {
        let offset = std::array::from_fn(|i| world[i] - self.start[i]);
        let at = add(self.origin, self.rotate(offset));
        (at.map(|x| x.value), at.map(|x| x.second))
    }

    /// Encloses the normal to the carried authored plane.
    pub fn normal(&self) -> [I; 3] {
        self.rotate(self.normal).map(|x| x.value)
    }
}

impl Motion {
    /// Proves that the section planes separate the entire authored material
    /// over the motion, then prepares its exact source chart correspondence.
    pub fn new(recipe: &SweepRecipe, dec: &mut Decider) -> Result<Self, Refusal> {
        let plane = recipe.plane.as_ref().ok_or_else(unseparated)?;
        if recipe.end_scale.value() != &Rational::one() {
            return Err(unseparated());
        }
        let [w, x, y, z] = plane.rotation.each_ref().map(|n| n.value().clone());
        let rotation = Rotation::from_quaternion(w, x, y, z)?;
        let matrix = rotation.rational_matrix().expect("rational rotation");
        let normal_r = RVector3::new(
            matrix[0][2].clone(),
            matrix[1][2].clone(),
            matrix[2][2].clone(),
        );
        let normal = vector(&normal_r);
        let radians = ball_interval(&Ball::pi_checked(ENCLOSURE_BITS, dec)?) / number(180);
        let twist = I::from_rational(recipe.twist_degrees.value()) * radians;
        let (start, spine) = match (&recipe.path, &recipe.guide) {
            (SweepPath::Line { start, end }, guide) => {
                let start = xyz(start);
                let end = xyz(end);
                let d = end.minus(&start);
                let delta = vector(&d);
                let speed = norm(delta);
                if speed.lo() <= 0.0 || !d.cross(&normal_r).is_zero() {
                    return Err(unseparated());
                }
                let tangent = delta.map(|x| x / speed);
                let project = |v: RVector3| d.cross(&v).cross(&d);
                let guide = match guide {
                    Some(SweepPath::Line { start: g0, end: g1 }) => {
                        let g0 = project(xyz(g0).minus(&start));
                        let g1 = project(xyz(g1).minus(&end));
                        let change = g1.minus(&g0);
                        let min = if change.is_zero() {
                            g0.clone()
                        } else {
                            let at = (-g0.dot(&change) / change.norm_squared())
                                .max(Rational::zero())
                                .min(Rational::one());
                            g0.plus(&change.scaled(&at))
                        };
                        if min.is_zero() {
                            return Err(refusal(
                                "field.sweep-guide-crossing",
                                "the guide meets the projected spine",
                            ));
                        }
                        [vector(&g0), vector(&change)]
                    }
                    None => [
                        std::array::from_fn(|i| I::from_rational(matrix[i][0])),
                        [number(0); 3],
                    ],
                    _ => return Err(unseparated()),
                };
                (
                    point(&start),
                    Spine::Line {
                        delta,
                        tangent,
                        guide,
                    },
                )
            }
            (
                SweepPath::Helix {
                    start,
                    centre,
                    axis,
                    angle,
                    rise,
                },
                None,
            ) if recipe.orientation == SweepOrientation::RotationMinimising
                && angle.value().abs() <= Rational::from_i64(90) =>
            {
                let start = xyz(start);
                let centre = xyz(centre);
                let [x, y, z] = axis.0.each_ref().map(|n| n.value().clone());
                let axis_r = RVector3::new(x, y, z);
                let radial_r = start.minus(&centre);
                if !axis_r.dot(&radial_r).is_zero() || angle.value().is_zero() {
                    return Err(unseparated());
                }
                let axis = vector(&axis_r);
                let axis_length = norm(axis);
                let radial = vector(&radial_r);
                let radius = norm(radial);
                if axis_length.lo() <= 0.0 || radius.lo() <= 0.0 {
                    return Err(unseparated());
                }
                let axis = axis.map(|x| x / axis_length);
                let transverse = cross(constants(axis), constants(radial)).map(|j| j.value);
                let theta = I::from_rational(angle.value()) * radians;
                let rise = I::from_rational(rise.value());
                let speed = ((radius * theta).square() + rise.square())
                    .sqrt()
                    .expect("squares");
                if speed.lo() <= 0.0 {
                    return Err(unseparated());
                }
                let spin = -theta * rise / speed;
                (
                    point(&start),
                    Spine::Helix {
                        centre: point(&centre),
                        radial,
                        transverse,
                        axis,
                        radius,
                        theta,
                        rise,
                        speed,
                        spin,
                    },
                )
            }
            _ => return Err(unseparated()),
        };
        let (_, base) = spine.frame(Jet::variable(number(0)), start, twist, dec)?;
        let first = base.map(|axis| axis.map(|x| x.value));
        let enclose = |x: &catallax_numbers::Real, dec: &mut Decider| {
            Ok::<I, Refusal>(ball_interval(&dec.enclose_budgeted_at(x, ENCLOSURE_BITS)?))
        };
        let mut charts = Vec::new();
        let mut profile_radius = 0.0_f64;
        for chart in source_charts(recipe, dec)? {
            let mut controls = Vec::new();
            let mut weights = Vec::new();
            for (p, w) in chart.curve.controls().zip(chart.curve.weights()) {
                let p = [
                    enclose(&p.x, dec)?,
                    enclose(&p.y, dec)?,
                    enclose(&p.z, dec)?,
                ];
                profile_radius =
                    profile_radius.max(norm(std::array::from_fn(|i| p[i] - start[i])).hi());
                controls.push(p);
                weights.push(enclose(w, dec)?);
            }
            charts.push(Chart {
                segment: chart.segment,
                controls,
                weights,
                domain: chart.curve.domain().map(Clone::clone),
            });
        }
        if charts.is_empty() {
            return Err(unseparated());
        }
        if let Spine::Helix {
            radius,
            theta,
            rise,
            speed,
            spin,
            ..
        } = spine
        {
            let c = dot(normal, first[2]);
            if c.contains_zero() {
                return Err(unseparated());
            }
            let n = if c.lo() > 0.0 {
                normal
            } else {
                normal.map(|x| -x)
            };
            let misalignment = norm(std::array::from_fn(|i| n[i] - first[2][i]));
            // For |theta * (s - t)| <= pi/2, sinc is at least 2/pi.
            // Subtract the oblique normal's departure from the tangent and
            // the maximum change in normal times the authored profile radius.
            // The remaining positive margin separates both section planes'
            // material for every nonzero parameter difference.
            let sinc = number(2) / (number(180) * radians);
            let centre = (sinc * (radius * theta).square() + rise.square()) / speed;
            let magnitude =
                |x: I| I::point(x.lo().abs().max(x.hi().abs())).expect("finite or infinite");
            let margin = centre
                - speed * misalignment
                - I::point(profile_radius).expect("bounded controls")
                    * (magnitude(theta) + magnitude(spin + twist));
            if margin.lo() <= 0.0 || !margin.is_finite() {
                return Err(unseparated());
            }
        }
        Ok(Self {
            spine,
            charts,
            start,
            first,
            normal,
            twist,
        })
    }

    /// All authored boundary charts, before sweeping.
    pub fn charts(&self) -> &[Chart] {
        &self.charts
    }

    /// Provides the authored section and derivative enclosure on this
    /// parameter interval. Intervals too broad to separate a guide can be split.
    pub fn section(&self, range: [&Rational; 2], dec: &mut Decider) -> Result<SectionBox, Refusal> {
        dec.charge_work(1)?;
        let t = Jet::variable(I::from_rational(range[0]).hull(&I::from_rational(range[1])));
        let (origin, axes) = self.spine.frame(t, self.start, self.twist, dec)?;
        Ok(SectionBox {
            origin,
            axes,
            start: self.start,
            first: self.first,
            normal: self.normal,
        })
    }
}
