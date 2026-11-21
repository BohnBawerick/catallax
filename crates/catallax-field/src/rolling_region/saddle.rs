//! A cylinder with a second cylinder or a sphere off its axis (4.2a step 3).
//!
//! The offsets of two such supports meet in a quartic, the ball path. The
//! field writes it the way Dupont, Lazard, Lazard and Petitjean (2008)
//! parametrize the meeting of two quadrics when one is ruled: along each
//! ruling of the offset cylinder, at angle `theta`, the point
//! `o + rho m(theta) + s a` lies on the other offset where a quadratic
//! `A s^2 + 2 B s + C = 0` holds, so `s = (-B + sigma sqrt(D)) / A` with
//! `D = B^2 - A C`. One square root, as their near-optimal parametrization has
//! for a smooth quartic. The field derives it from the two supports, the ball
//! side of each and the radius alone; nothing here reads a B-rep fillet module
//! or a fillet definition.
//!
//! Where `D > 0` at every angle the two roots are two rings round the
//! cylinder, and the piece keeps the one whose sign `sigma` its edge has. Where
//! `D` changes sign, the path is one loop that winds round neither axis, the
//! two roots joining at its branch points; the piece then covers both roots
//! over the whole circle, and panels where `D < 0` hold no station. The cover
//! treats them as empty. Near a branch point the derivative of `s` grows
//! without bound, so the frame's tangent is `sqrt(D)` times the derivative,
//! which stays finite and nonzero at a simple root of `D`. No closed form
//! gives the nearest station, so membership and distance take the panel cover
//! of step 1, as OpenCASCADE's `Extrema_GExtPC` samples a general curve.
//!
//! Two equal cylinders whose axes meet have `D = L^2` for a linear `L`: the
//! quartic splits into two ellipses in the planes that bisect the axes, which
//! cross at two poles where the ball touches both cylinders at one point.
//! There the root is linear, `s = -rho (m . nu) / (a . nu)` for the plane's
//! normal `nu`, and the piece covers both ellipses whole.
//!
//! The fan is the Voronoi zone of the path in the solid of ball centres. A
//! point is in it when its projection onto either offset falls outside the
//! other offset, as for the conics of step 2; both projections are radial,
//! so these walls are exact. With the ball side of each support they bound
//! the region near the path. A box keeps far zones out, such as the other
//! ring of the offsets: the hull of the centres, feet and edge points over the
//! cover, the chord cube of each support the ball lies inside (its arc from
//! the edge to the foot bulges past that chord), padded by `r`. The box is a
//! band: positive in the region, never a face.

use super::elementary::{balls, local, point_balls, position, positive, reach, vector};
use super::path::{Cover, Frame, Path};
use super::revolved::{angle, basis_of, chord_cube, hull_points, root};
use super::{BITS, gap, interval, q};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

/// The support paired with the parametrizing cylinder.
#[derive(Clone, Debug)]
pub enum Partner {
    /// A cylinder whose axis is not parallel to the first.
    Cylinder {
        /// Any point of the axis.
        point: Point3,
        /// Nonzero axis direction.
        axis: Vector3,
        /// Positive radius in mm.
        radius: Real,
    },
    /// A sphere.
    Sphere {
        /// Centre.
        centre: Point3,
        /// Positive radius in mm.
        radius: Real,
    },
}

/// Which roots of the quadratic a piece covers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// One ring round the cylinder: the root with `sigma = +1`, or `-1`.
    Ring(bool),
    /// One loop that winds round neither axis: both roots where `D >= 0`.
    Loop,
    /// Two equal cylinders on axes that meet: the offsets meet in the two
    /// planes through the meeting point that bisect the axes, in two ellipses
    /// that cross at two poles. One piece holds all four bands, so every ball
    /// of the crossing clips every band's region.
    Crossing,
}

/// How the roots of a quadratic family lie over the circle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// `D > 0` at every angle: two rings.
    Rings,
    /// `D >= 0` on one arc and `D < 0` on the rest: one loop.
    Loop,
    /// Anything else, or not proven.
    Other,
}

fn side(outside: bool) -> Real {
    Real::from_i64(if outside { 1 } else { -1 })
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &root(&v.norm_squared())))
}

/// The partner's offset in exact form.
#[derive(Clone, Debug)]
enum Other {
    Cylinder {
        point: Point3,
        axis: Vector3,
        radius: Real,
    },
    Sphere {
        centre: Point3,
        radius: Real,
    },
}
impl Other {
    fn origin(&self) -> &Point3 {
        match self {
            Other::Cylinder { point, .. } => point,
            Other::Sphere { centre, .. } => centre,
        }
    }
    fn radius(&self) -> &Real {
        match self {
            Other::Cylinder { radius, .. } | Other::Sphere { radius, .. } => radius,
        }
    }
    /// The part of `x - origin` the partner measures: radial for a cylinder.
    fn spoke(&self, x: &Vector3) -> Vector3 {
        let v = x.minus(&position(self.origin()));
        match self {
            Other::Cylinder { axis, .. } => v.minus(&axis.scaled(&axis.dot(&v))),
            Other::Sphere { .. } => v,
        }
    }
}

/// The ball path between a cylinder and a second cylinder or a sphere.
#[derive(Clone, Debug)]
pub struct Saddle {
    point: Point3,
    axis: Vector3,
    basis: [Vector3; 2],
    body: Real,
    other: Other,
    outside: [bool; 2],
    // A = 1 - (a.a2)^2 for a cylinder, 1 for a sphere, and a.a2.
    leading: Real,
    lean: Real,
    branch: Branch,
    radius: Real,
    // The padded box of the region, a band.
    band: [[Rational; 2]; 3],
    // The boxes of the panels of the cover, for clearance.
    boxes: Vec<[Ball; 3]>,
}

/// Interval data shared by the path and the edge.
#[derive(Clone, Debug)]
pub struct SaddleCover {
    point: [Ball; 3],
    axis: [Ball; 3],
    basis: [[Ball; 3]; 2],
    pi: Ball,
    other_point: [Ball; 3],
    // The partner's axis for a cylinder.
    other_axis: Option<[Ball; 3]>,
    leading: Ball,
    lean: Ball,
    // The bisecting planes' normals of a crossing.
    planes: Option<[[Ball; 3]; 2]>,
    radius: Ball,
    // Radii of the cylinder and the partner, offset for the path.
    path: [Ball; 2],
    edge: [Ball; 2],
    sides: [Ball; 2],
    branch: Branch,
}

fn add(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &b[i])
}
fn sub(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn times(t: &Ball, a: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| t * &a[i])
}
fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}
fn divide(a: &Ball, b: &Ball) -> Result<Ball, Refusal> {
    a.checked_div(b)
        .ok_or_else(|| gap("a saddle divisor enclosure contains zero"))
}
fn hull3(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| a[i].hull(&b[i]))
}

/// One root over a panel of angles: the point, the discriminant and the
/// derivative data the frame needs.
struct Root {
    m: [Ball; 3],
    disc: Ball,
    point: [Ball; 3],
    // sqrt(D) times the point's derivative.
    tangent: [Ball; 3],
}

impl SaddleCover {
    /// The panels of angle `pi u'` and their signs for a panel of `u`: one for
    /// a ring, and for a loop the lower root on `[-1, 0]` and the upper one on
    /// `[0, 1]`, each over the whole circle.
    fn parts(&self, u: &Ball) -> Vec<(Ball, Ball)> {
        let one = Ball::from_rational(&Rational::one(), BITS);
        match self.branch {
            Branch::Ring(upper) => {
                let sigma = if upper { one } else { -&one };
                vec![(u.clone(), sigma)]
            }
            Branch::Loop | Branch::Crossing => {
                let (lo, hi) = (u.lo(), u.hi());
                let zero = Rational::zero();
                let mut parts = Vec::with_capacity(2);
                let two = q(2);
                if lo < zero {
                    let top = hi.clone().min(zero.clone());
                    parts.push((
                        interval(&(&(&two * &lo) + &q(1)), &(&(&two * &top) + &q(1))),
                        -&one,
                    ));
                }
                if hi > zero {
                    let bottom = lo.max(zero.clone());
                    parts.push((
                        interval(&(&(&two * &bottom) - &q(1)), &(&(&two * &hi) - &q(1))),
                        one.clone(),
                    ));
                }
                // The station u = 0 is where the lower root ends, at angle pi.
                if parts.is_empty() {
                    parts.push((interval(&q(1), &q(1)), -&one));
                }
                parts
            }
        }
    }

    /// The quadratic's root of sign `sigma` over a panel, for the cylinder
    /// radius `rho` and partner radius `big`.
    fn root(&self, u: &Ball, sigma: &Ball, [rho, big]: &[Ball; 2]) -> Result<Root, Refusal> {
        let (cos, sin, m) = angle(&self.pi, &self.basis, u);
        let turn = add(&times(&-&sin, &self.basis[0]), &times(&cos, &self.basis[1]));
        if let Some(planes) = &self.planes {
            // The plane through the meeting point: (rho m + s a) . nu = 0.
            let nu = if sigma.lo() > Rational::zero() {
                &planes[0]
            } else {
                &planes[1]
            };
            let across = dot(&self.axis, nu);
            let s = divide(&-&(rho * &dot(&m, nu)), &across)?;
            let rise = divide(&-&(rho * &dot(&turn, nu)), &across)?;
            let point = add(&add(&self.point, &times(rho, &m)), &times(&s, &self.axis));
            return Ok(Root {
                m,
                disc: Ball::from_rational(&Rational::one(), BITS),
                point,
                tangent: add(&times(rho, &turn), &times(&rise, &self.axis)),
            });
        }
        let d = add(&sub(&self.point, &self.other_point), &times(rho, &m));
        let dt = times(rho, &turn);
        let (b, c, db, dc) = match &self.other_axis {
            Some(axis) => {
                let along = dot(&d, axis);
                let turned = dot(&dt, axis);
                (
                    &dot(&d, &self.axis) - &(&self.lean * &along),
                    &(&dot(&d, &d) - &along.square()) - &big.square(),
                    -&(&self.lean * &turned),
                    &(&dot(&d, &dt) - &(&along * &turned)) * &Ball::from_rational(&q(2), BITS),
                )
            }
            None => (
                dot(&d, &self.axis),
                &dot(&d, &d) - &big.square(),
                Ball::from_rational(&Rational::zero(), BITS),
                &dot(&d, &dt) * &Ball::from_rational(&q(2), BITS),
            ),
        };
        let a = &self.leading;
        let disc = &b.square() - &(a * &c);
        // A panel may hold angles past a branch point; its stations all have
        // D >= 0, so the enclosure keeps only that part.
        let zero = Rational::zero();
        let kept = interval(&disc.lo().max(zero.clone()), &disc.hi().max(zero));
        let rootd = kept
            .sqrt()
            .ok_or_else(|| gap("negative discriminant enclosure"))?;
        let s = divide(&(&-&b + &(sigma * &rootd)), a)?;
        let point = add(&add(&self.point, &times(rho, &m)), &times(&s, &self.axis));
        let ddisc = &(&(&b * &db) * &Ball::from_rational(&q(2), BITS)) - &(a * &dc);
        let half = Ball::from_rational(&(Rational::one() / q(2)), BITS);
        let rise = divide(&(&-&(&db * &rootd) + &(&(sigma * &ddisc) * &half)), a)?;
        let tangent = add(&times(&(&rootd * rho), &turn), &times(&rise, &self.axis));
        Ok(Root {
            m,
            disc,
            point,
            tangent,
        })
    }

    /// Path frame over one part.
    fn part_frame(&self, u: &Ball, sigma: &Ball) -> Result<Frame, Refusal> {
        let root = self.root(u, sigma, &self.path)?;
        let r = &self.radius;
        let c = root.point;
        let first = sub(&c, &times(&(&self.sides[0] * r), &root.m));
        let spoke = match &self.other_axis {
            Some(axis) => {
                let v = sub(&c, &self.other_point);
                sub(&v, &times(&dot(&v, axis), axis))
            }
            None => sub(&c, &self.other_point),
        };
        let second = sub(
            &c,
            &times(&divide(&(&self.sides[1] * r), &self.path[1])?, &spoke),
        );
        Ok(Frame {
            centre: c,
            tangent: root.tangent,
            feet: [first, second],
            side: None,
        })
    }

    /// Edge points over a panel, `None` where the edge has no point.
    fn edge(&self, u: &Ball) -> Result<Option<[Ball; 3]>, Refusal> {
        let mut found: Option<[Ball; 3]> = None;
        for (part, sigma) in self.parts(u) {
            let root = self.root(&part, &sigma, &self.edge)?;
            if root.disc.hi() < Rational::zero() {
                continue;
            }
            found = Some(match found {
                Some(p) => hull3(&p, &root.point),
                None => root.point,
            });
        }
        Ok(found)
    }

    fn disc_parts(&self, u: &Ball, radii: &[Ball; 2]) -> Result<Vec<Ball>, Refusal> {
        self.parts(u)
            .iter()
            .map(|(part, sigma)| Ok(self.root(part, sigma, radii)?.disc))
            .collect()
    }
}

impl Cover for SaddleCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        Ok(self.frame(u)?.centre)
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let mut frame: Option<Frame> = None;
        for (part, sigma) in self.parts(u) {
            let next = self.part_frame(&part, &sigma)?;
            frame = Some(match frame {
                None => next,
                Some(f) => Frame {
                    centre: hull3(&f.centre, &next.centre),
                    tangent: hull3(&f.tangent, &next.tangent),
                    feet: [
                        hull3(&f.feet[0], &next.feet[0]),
                        hull3(&f.feet[1], &next.feet[1]),
                    ],
                    side: None,
                },
            });
        }
        frame.ok_or_else(|| gap("an empty saddle panel"))
    }
    fn exists(&self, u: &Ball) -> Result<bool, Refusal> {
        Ok(self
            .disc_parts(u, &self.path)?
            .iter()
            .all(|d| d.lo() >= Rational::zero()))
    }
    fn empty(&self, u: &Ball) -> Result<bool, Refusal> {
        Ok(self
            .disc_parts(u, &self.path)?
            .iter()
            .all(|d| d.hi() < Rational::zero()))
    }
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        if !self.exists(u)? {
            return Ok(false);
        }
        let c = self.centre(u)?;
        let d = sub(p, &c);
        Ok(dot(&d, &d).hi() < self.radius.square().lo())
    }
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        if self.empty(u)? {
            return Ok(true);
        }
        let c = self.centre(u)?;
        let d = sub(p, &c);
        Ok(dot(&d, &d).lo() > self.radius.square().hi())
    }
}

/// How the roots of the quadratic family lie for the cylinder radius `rho`
/// and the partner radius `big`, proven over panels of the circle.
fn shape(cover: &SaddleCover, radii: &[Ball; 2], dec: &mut Decider) -> Result<Shape, Refusal> {
    // Each of 256 panels is split until its sign is proven or it is 1/2^14
    // of the circle wide.
    const PANELS: i64 = 256;
    const DEPTH: usize = 6;
    let ring = SaddleCover {
        branch: Branch::Ring(true),
        ..cover.clone()
    };
    // Signs per panel: Some(true) proven positive, Some(false) proven
    // negative, None touches zero.
    let mut signs = Vec::with_capacity(PANELS as usize);
    for i in 0..PANELS {
        let lo = Rational::from_i64(2 * i - PANELS) / Rational::from_i64(PANELS);
        let hi = Rational::from_i64(2 * i + 2 - PANELS) / Rational::from_i64(PANELS);
        let mut pending = vec![(lo, hi, 0)];
        let mut panel = Vec::new();
        while let Some((lo, hi, depth)) = pending.pop() {
            dec.charge_work(1)?;
            let d = &ring.disc_parts(&interval(&lo, &hi), radii)?[0];
            if d.lo() > Rational::zero() {
                panel.push(Some(true));
            } else if d.hi() < Rational::zero() {
                panel.push(Some(false));
            } else if depth == DEPTH {
                panel.push(None);
            } else {
                let mid = (&lo + &hi) / q(2);
                pending.push((mid.clone(), hi, depth + 1));
                pending.push((lo, mid, depth + 1));
            }
        }
        signs.push(if panel.iter().all(|s| *s == Some(true)) {
            Some(true)
        } else if panel.iter().all(|s| *s == Some(false)) {
            Some(false)
        } else {
            None
        });
    }
    if signs.iter().all(|s| *s == Some(true)) {
        return Ok(Shape::Rings);
    }
    // One loop: the panels that may hold a station form one run round the
    // circle, and at least one panel has none.
    let n = signs.len();
    let present: Vec<bool> = signs.iter().map(|s| *s != Some(false)).collect();
    if present.iter().all(|p| *p) || !present.iter().any(|p| *p) {
        return Ok(Shape::Other);
    }
    let starts: Vec<usize> = (0..n)
        .filter(|&i| present[i] && !present[(i + n - 1) % n])
        .collect();
    let [start] = starts[..] else {
        return Ok(Shape::Other);
    };
    // Within the run, D may touch zero only at its two ends, the branch
    // points; a zero inside is a node where the loop pinches.
    let run: Vec<Option<bool>> = (0..n)
        .map(|k| signs[(start + k) % n])
        .take_while(|s| *s != Some(false))
        .collect();
    let first = run.iter().position(|s| *s == Some(true));
    let last = run.iter().rposition(|s| *s == Some(true));
    Ok(match (first, last) {
        (Some(first), Some(last)) if run[first..=last].iter().all(|s| *s == Some(true)) => {
            Shape::Loop
        }
        _ => Shape::Other,
    })
}

impl Saddle {
    /// The interval cover of the cylinder and the partner.
    #[allow(clippy::too_many_arguments)]
    fn frame_cover(
        point: &Point3,
        axis: &Vector3,
        basis: &[Vector3; 2],
        other: &Other,
        leading: &Real,
        lean: &Real,
        radii: [&Real; 4],
        sides: [bool; 2],
        radius: &Real,
        branch: Branch,
        dec: &mut Decider,
    ) -> Result<SaddleCover, Refusal> {
        let ball = |x: &Real, dec: &mut Decider| enclose(x, BITS, dec);
        let [body, offset, big, offset_big] = radii;
        Ok(SaddleCover {
            point: point_balls(point, dec)?,
            axis: balls(axis, dec)?,
            basis: [balls(&basis[0], dec)?, balls(&basis[1], dec)?],
            pi: Ball::pi_checked(BITS, dec)?,
            other_point: point_balls(other.origin(), dec)?,
            other_axis: match other {
                Other::Cylinder { axis, .. } => Some(balls(axis, dec)?),
                Other::Sphere { .. } => None,
            },
            leading: ball(leading, dec)?,
            lean: ball(lean, dec)?,
            planes: match (branch, other) {
                (Branch::Crossing, Other::Cylinder { axis: second, .. }) => Some([
                    balls(&axis.plus(second), dec)?,
                    balls(&axis.minus(second), dec)?,
                ]),
                _ => None,
            },
            radius: ball(radius, dec)?,
            path: [ball(offset, dec)?, ball(offset_big, dec)?],
            edge: [ball(body, dec)?, ball(big, dec)?],
            sides: [ball(&side(sides[0]), dec)?, ball(&side(sides[1]), dec)?],
            branch,
        })
    }

    /// The shape of the edge where a cylinder meets the partner, about the
    /// cylinder's axis.
    pub fn edge_shape(
        cylinder: (&Point3, &Vector3, &Real),
        partner: &Partner,
        dec: &mut Decider,
    ) -> Result<Shape, Refusal> {
        let (other, leading, lean) = Self::other(partner, cylinder.1, dec)?;
        let axis = unit(cylinder.1);
        let basis = basis_of(&axis, dec)?;
        let zero = Real::zero();
        let cover = Self::frame_cover(
            cylinder.0,
            &axis,
            &basis,
            &other,
            &leading,
            &lean,
            [cylinder.2, &zero, other.radius(), &zero],
            [true, true],
            &zero,
            Branch::Ring(true),
            dec,
        )?;
        let edge = cover.edge.clone();
        shape(&cover, &edge, dec)
    }

    fn other(
        partner: &Partner,
        axis: &Vector3,
        dec: &mut Decider,
    ) -> Result<(Other, Real, Real), Refusal> {
        let a = unit(axis);
        Ok(match partner {
            Partner::Cylinder {
                point,
                axis: second,
                radius,
            } => {
                if !positive(&second.norm_squared(), dec)? || !positive(radius, dec)? {
                    return Err(reach(
                        "a partner cylinder needs an axis and a positive radius",
                    ));
                }
                let second = unit(second);
                let lean = a.dot(&second);
                let leading = &Real::one() - &lean.square();
                if !positive(&leading, dec)? {
                    return Err(reach("parallel cylinders meet in lines, not a saddle"));
                }
                (
                    Other::Cylinder {
                        point: point.clone(),
                        axis: second,
                        radius: radius.clone(),
                    },
                    leading,
                    lean,
                )
            }
            Partner::Sphere { centre, radius } => {
                if !positive(radius, dec)? {
                    return Err(reach("a sphere radius must be positive"));
                }
                (
                    Other::Sphere {
                        centre: centre.clone(),
                        radius: radius.clone(),
                    },
                    Real::one(),
                    Real::zero(),
                )
            }
        })
    }

    /// The root sign of the edge at one of its points `w`: `A s + B`, which
    /// is `sigma sqrt(D)`. Zero at a branch point, which is refused.
    pub fn edge_sign(
        cylinder: (&Point3, &Vector3, &Real),
        partner: &Partner,
        w: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let (other, leading, lean) = Self::other(partner, cylinder.1, dec)?;
        let axis = unit(cylinder.1);
        let v = w.minus(cylinder.0);
        let s = axis.dot(&v);
        let spoke = v.minus(&axis.scaled(&s));
        let rho2 = spoke.norm_squared();
        if !positive(&rho2, dec)? {
            return Err(reach("an edge point lies on the cylinder's axis"));
        }
        let m = spoke.scaled(&(&Real::one() / &root(&rho2)));
        let d = cylinder.0.minus(other.origin()).plus(&m.scaled(cylinder.2));
        let b = match &other {
            Other::Cylinder { axis: second, .. } => &d.dot(&axis) - &(&lean * &d.dot(second)),
            Other::Sphere { .. } => d.dot(&axis),
        };
        match sign(&(&(&leading * &s) + &b), dec)? {
            Sign::Positive => Ok(true),
            Sign::Negative => Ok(false),
            Sign::Zero => Err(reach("the edge sample is a branch point of the saddle")),
        }
    }

    /// Derives the path from the cylinder through `point` along `axis` of
    /// radius `body`, the partner, the ball side of each and the radius, for
    /// the branch the edge has.
    pub fn new(
        cylinder: (Point3, Vector3, Real),
        partner: Partner,
        outside: [bool; 2],
        branch: Branch,
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        let (point, axis, body) = cylinder;
        if !positive(&radius, dec)? {
            return Err(reach("the rolling radius must be positive"));
        }
        if !positive(&axis.norm_squared(), dec)? || !positive(&body, dec)? {
            return Err(reach("a cylinder needs an axis and a positive radius"));
        }
        let (other, leading, lean) = Self::other(&partner, &axis, dec)?;
        let axis = unit(&axis);
        let offset = &body + &(&side(outside[0]) * &radius);
        let offset_big = other.radius() + &(&side(outside[1]) * &radius);
        if !positive(&offset, dec)? || !positive(&offset_big, dec)? {
            return Err(reach(
                "the ball is no smaller than the support it rolls inside",
            ));
        }
        let basis = basis_of(&axis, dec)?;
        let cover = Self::frame_cover(
            &point,
            &axis,
            &basis,
            &other,
            &leading,
            &lean,
            [&body, &offset, other.radius(), &offset_big],
            outside,
            &radius,
            branch,
            dec,
        )?;
        let wanted = match branch {
            Branch::Ring(_) => Some(Shape::Rings),
            Branch::Loop => Some(Shape::Loop),
            Branch::Crossing => {
                Self::crossing(&point, &axis, &body, &other, outside, dec)?;
                None
            }
        };
        if let Some(wanted) = wanted {
            if shape(&cover, &cover.edge, dec)? != wanted {
                return Err(reach("the edge does not have the shape of its piece"));
            }
            if shape(&cover, &cover.path, dec)? != wanted {
                return Err(reach(
                    "the offsets meet in another shape than the edge; the ball changes the topology",
                ));
            }
        }
        // Panel boxes: centres, feet and edge points over each panel and its
        // two neighbours, with the chord cube of a support the ball lies
        // inside, padded by r.
        const PANELS: i64 = 64;
        let count = match branch {
            Branch::Ring(_) => PANELS,
            Branch::Loop | Branch::Crossing => 2 * PANELS,
        };
        let panel = |i: i64| {
            let i = i.rem_euclid(count);
            let lo = Rational::from_i64(2 * i - count) / Rational::from_i64(count);
            let hi = Rational::from_i64(2 * i + 2 - count) / Rational::from_i64(count);
            interval(&lo, &hi)
        };
        let mut points: Vec<Vec<[Ball; 3]>> = Vec::with_capacity(count as usize);
        for i in 0..count {
            dec.charge_work(1)?;
            let u = panel(i);
            let mut found = Vec::new();
            let edge = cover.edge(&u)?;
            if !cover.empty(&u)? {
                let frame = cover.frame(&u)?;
                if let Some(edge) = &edge {
                    for (k, out) in outside.iter().enumerate() {
                        if !out {
                            found.extend(chord_cube(edge, &frame.feet[k])?);
                        }
                    }
                }
                found.push(frame.centre);
                found.extend(frame.feet);
            }
            if let Some(edge) = edge {
                found.push(edge);
            }
            points.push(found);
        }
        let pad = enclose(&radius, BITS, dec)?.hi();
        let padded = |b: [Ball; 3]| -> [Ball; 3] {
            b.map(|x| interval(&(&x.lo() - &pad), &(&x.hi() + &pad)))
        };
        let mut boxes = Vec::with_capacity(count as usize);
        for i in 0..count {
            let mut near: Vec<[Ball; 3]> = Vec::new();
            // A loop's two roots meet at its branch points, so a panel's
            // neighbours include the other root's panels at the same angles.
            let mates: &[i64] = match branch {
                Branch::Ring(_) => &[0],
                Branch::Loop | Branch::Crossing => &[0, PANELS],
            };
            for mate in mates {
                for step in [-1, 0, 1] {
                    near.extend(points[(i + mate + step).rem_euclid(count) as usize].clone());
                }
            }
            if !near.is_empty() {
                boxes.push(padded(hull_points(&near)));
            }
        }
        if boxes.is_empty() {
            return Err(reach("the saddle has no station"));
        }
        let hull = hull_points(&boxes);
        let band = hull.map(|b| [b.lo(), b.hi()]);
        Ok(Self {
            point,
            axis,
            basis,
            body,
            other,
            outside,
            leading,
            lean,
            branch,
            radius,
            band,
            boxes,
        })
    }

    /// Checks a crossing: the partner is a cylinder of the same radius whose
    /// axis meets this one at `point`, and the ball lies on the same side of
    /// both, so the offsets are equal too.
    fn crossing(
        point: &Point3,
        axis: &Vector3,
        body: &Real,
        other: &Other,
        outside: [bool; 2],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let Other::Cylinder {
            point: meeting,
            axis: second,
            radius,
        } = other
        else {
            return Err(reach("a crossing needs two cylinders"));
        };
        if sign(&point.minus(meeting).norm_squared(), dec)? != Sign::Zero
            || sign(&axis.cross(second).norm_squared(), dec)? == Sign::Zero
        {
            return Err(reach("a crossing's axes must meet at its origin"));
        }
        if sign(&(body - radius), dec)? != Sign::Zero {
            return Err(reach("a crossing needs two cylinders of one radius"));
        }
        if outside[0] != outside[1] {
            return Err(reach(
                "the ball lies inside one crossing cylinder and outside the other, so the offsets do not cross",
            ));
        }
        Ok(())
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    fn offsets(&self) -> [Real; 2] {
        [
            &self.body + &(&side(self.outside[0]) * &self.radius),
            self.other.radius() + &(&side(self.outside[1]) * &self.radius),
        ]
    }

    /// Applies an exact rigid motion, reflections included. The band's box is
    /// moved by its corners when they stay rational, and otherwise widened to
    /// the box of the moved panel boxes.
    pub fn transform(&mut self, pose: &Isometry) {
        self.point = pose.apply_point(&self.point);
        self.axis = pose.apply_vector(&self.axis);
        for e in &mut self.basis {
            *e = pose.apply_vector(e);
        }
        match &mut self.other {
            Other::Cylinder { point, axis, .. } => {
                *point = pose.apply_point(point);
                *axis = pose.apply_vector(axis);
            }
            Other::Sphere { centre, .. } => *centre = pose.apply_point(centre),
        }
        // Moving boxes needs exact corners; the panel boxes and the band are
        // rebuilt from their rational corners.
        self.band = super::moved_box(pose, &self.band);
        self.boxes = self
            .boxes
            .iter()
            .map(|b| {
                let m = super::moved_box(pose, &b.each_ref().map(|x| [x.lo(), x.hi()]));
                m.map(|[lo, hi]| interval(&lo, &hi))
            })
            .collect();
    }

    /// An axis-aligned box holding the region.
    pub fn hull(&self, _dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.band.each_ref().map(|[lo, hi]| interval(lo, hi)))
    }

    /// Boxes that together hold the region, one per panel of the cover.
    pub fn boxes(&self, _dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        Ok(self.boxes.clone())
    }

    /// The path in binary64, for display values and starting guesses.
    fn rough(&self) -> Rough {
        let f = |x: &Real| x.to_f64();
        let fv = |v: &Vector3| [f(&v.x), f(&v.y), f(&v.z)];
        let fp = |v: &Point3| [f(&v.x), f(&v.y), f(&v.z)];
        let [rho, big] = self.offsets().map(|x| f(&x));
        Rough {
            point: fp(&self.point),
            axis: fv(&self.axis),
            basis: [fv(&self.basis[0]), fv(&self.basis[1])],
            rho,
            big,
            other: fp(self.other.origin()),
            other_axis: match &self.other {
                Other::Cylinder { axis, .. } => Some(fv(axis)),
                Other::Sphere { .. } => None,
            },
            lean: f(&self.lean),
            leading: f(&self.leading),
            branch: self.branch,
        }
    }

    /// Tentative binary64 value for display: negative inside the region. The
    /// square root is the caller's. The nearest centre comes from 512 samples
    /// of the path, so near the face it is only a display value.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let rough = self.rough();
        let r = self.radius.to_f64();
        let mut nearest = f64::INFINITY;
        const SAMPLES: usize = 512;
        for i in 0..SAMPLES {
            let u = -1.0 + 2.0 * (i as f64 + 0.5) / SAMPLES as f64;
            if let Some(c) = rough.centre(u, sqrt) {
                nearest = nearest.min(sqrt(dot3(sub3(p, c), sub3(p, c))));
            }
        }
        let (o, a) = (rough.point, rough.axis);
        let v1 = sub3(p, o);
        let w1 = sub3(v1, scale3(dot3(v1, a), a));
        let d1 = sqrt(dot3(w1, w1)) - self.body.to_f64();
        let v2 = sub3(p, rough.other);
        let w2 = match rough.other_axis {
            Some(a2) => sub3(v2, scale3(dot3(v2, a2), a2)),
            None => v2,
        };
        let d2 = sqrt(dot3(w2, w2)) - self.other.radius().to_f64();
        let s1 = if self.outside[0] { d1 } else { -d1 };
        let s2 = if self.outside[1] { d2 } else { -d2 };
        let inside_box = (0..3).all(|i| {
            p[i] > self.band[i][0].to_f64_nearest() && p[i] < self.band[i][1].to_f64_nearest()
        });
        if inside_box && s1 >= 0.0 && s2 >= 0.0 && nearest > r {
            // Within the box and the corner, beyond every sampled ball; the
            // Voronoi walls are left to the exact query.
            (-s1).max(-s2).max(r - nearest)
        } else {
            (nearest - r).abs().max(f64::MIN_POSITIVE)
        }
    }
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale3(t: f64, a: [f64; 3]) -> [f64; 3] {
    [t * a[0], t * a[1], t * a[2]]
}
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// A [`Saddle`] in binary64.
struct Rough {
    point: [f64; 3],
    axis: [f64; 3],
    basis: [[f64; 3]; 2],
    rho: f64,
    big: f64,
    other: [f64; 3],
    other_axis: Option<[f64; 3]>,
    lean: f64,
    leading: f64,
    branch: Branch,
}
impl Rough {
    /// The centre at `u`, mapped as the cover maps it; `None` where the
    /// station does not exist.
    fn centre(&self, u: f64, sqrt: fn(f64) -> f64) -> Option<[f64; 3]> {
        use std::f64::consts::PI;
        let (sigma, t) = match self.branch {
            Branch::Ring(upper) => (if upper { 1.0 } else { -1.0 }, PI * u),
            Branch::Loop | Branch::Crossing if u < 0.0 => (-1.0, PI * (2.0 * u + 1.0)),
            Branch::Loop | Branch::Crossing => (1.0, PI * (2.0 * u - 1.0)),
        };
        let (o, a) = (self.point, self.axis);
        let m = add3(
            scale3(t.cos(), self.basis[0]),
            scale3(t.sin(), self.basis[1]),
        );
        let s = match (self.branch, self.other_axis) {
            (Branch::Crossing, Some(a2)) => {
                let nu = if sigma > 0.0 {
                    add3(a, a2)
                } else {
                    sub3(a, a2)
                };
                -self.rho * dot3(m, nu) / dot3(a, nu)
            }
            _ => {
                let d = add3(sub3(o, self.other), scale3(self.rho, m));
                let (b, c) = match self.other_axis {
                    Some(a2) => (
                        dot3(d, a) - self.lean * dot3(d, a2),
                        dot3(d, d) - dot3(d, a2).powi(2) - self.big * self.big,
                    ),
                    None => (dot3(d, a), dot3(d, d) - self.big * self.big),
                };
                let disc = b * b - self.leading * c;
                if disc < 0.0 {
                    return None;
                }
                (-b + sigma * sqrt(disc)) / self.leading
            }
        };
        Some(add3(add3(o, scale3(self.rho, m)), scale3(s, a)))
    }

    /// The parameter of the centre nearest `p`: the best of 1024 samples,
    /// polished by a golden-section search on its neighbourhood.
    fn nearest(&self, p: [f64; 3]) -> Option<f64> {
        const SAMPLES: usize = 1024;
        let gap = |u: f64| {
            self.centre(u, f64::sqrt)
                .map(|c| dot3(sub3(p, c), sub3(p, c)))
                .unwrap_or(f64::INFINITY)
        };
        let step = 2.0 / SAMPLES as f64;
        let mut best = (f64::INFINITY, 0.0);
        for i in 0..SAMPLES {
            let u = -1.0 + step * (i as f64 + 0.5);
            let g = gap(u);
            if g < best.0 {
                best = (g, u);
            }
        }
        if !best.0.is_finite() {
            return None;
        }
        let ratio = (5f64.sqrt() - 1.0) / 2.0;
        let (mut lo, mut hi) = ((best.1 - step).max(-1.0), (best.1 + step).min(1.0));
        for _ in 0..80 {
            let a = hi - ratio * (hi - lo);
            let b = lo + ratio * (hi - lo);
            if gap(a) < gap(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        Some((lo + hi) / 2.0)
    }
}

impl Path for Saddle {
    type Cover = SaddleCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = vector(local);
        let [rho, big] = self.offsets();
        let v1 = x.minus(&position(&self.point));
        let s1 = self.axis.dot(&v1);
        let w1 = v1.minus(&self.axis.scaled(&s1));
        let n1 = w1.norm_squared();
        let w2 = self.other.spoke(&x);
        let n2 = w2.norm_squared();
        let mut walls = vec![
            &side(self.outside[0]) * &(&n1 - &self.body.square()),
            &side(self.outside[1]) * &(&n2 - &self.other.radius().square()),
        ];
        // The region keeps away from both axes (its box shows it), so a point
        // on either is outside it.
        if sign(&n1, dec)? == Sign::Zero || sign(&n2, dec)? == Sign::Zero {
            walls.push(Real::from_i64(-1));
            return Ok(walls);
        }
        // Projections onto each offset, radial for both.
        let p1 = position(&self.point)
            .plus(&self.axis.scaled(&s1))
            .plus(&w1.scaled(&(&rho / &root(&n1))));
        let q2 = self.other.spoke(&p1);
        walls.push(-&(&side(self.outside[1]) * &(&q2.norm_squared() - &big.square())));
        let p2 = x.minus(&w2).plus(&w2.scaled(&(&big / &root(&n2))));
        let v = p2.minus(&position(&self.point));
        let q1 = v.minus(&self.axis.scaled(&self.axis.dot(&v)));
        walls.push(-&(&side(self.outside[0]) * &(&q1.norm_squared() - &rho.square())));
        Ok(walls)
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn guess(&self, local: &[Real; 3]) -> Option<Rational> {
        let p = local.each_ref().map(|x| x.to_f64());
        Rational::from_f64_exact(self.rough().nearest(p)?)
    }
    fn bands(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let mut bands = Vec::with_capacity(6);
        for (x, [lo, hi]) in local.iter().zip(&self.band) {
            bands.push(x - &Real::from(lo.clone()));
            bands.push(&Real::from(hi.clone()) - x);
        }
        Ok(bands)
    }
    fn cover(&self, dec: &mut Decider) -> Result<SaddleCover, Refusal> {
        let [rho, big] = self.offsets();
        Self::frame_cover(
            &self.point,
            &self.axis,
            &self.basis,
            &self.other,
            &self.leading,
            &self.lean,
            [&self.body, &rho, self.other.radius(), &big],
            self.outside,
            &self.radius,
            self.branch,
            dec,
        )
    }
}
