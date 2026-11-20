//! Two curved surfaces of revolution about one axis (4.2a step 3).
//!
//! Cylinders, cone nappes and spheres that share an axis meet in circles, and
//! so do their offsets, which stay coaxial: a cylinder or cone moves its
//! profile line by `r` along the line's normal, a sphere changes its radius by
//! `r`. The field works in the meridian half-plane of the axis, with
//! coordinates `(rho, h)`: the distance from the axis and the axial station.
//! There each support is a line or a circle, and the ball centre is where the
//! two offset profiles meet, solved in closed form: Cramer's rule for two
//! lines, a quadratic along a line for a line and a circle, and the radical
//! line for two circles. Of several meetings the field keeps the one whose
//! feet lie on the two faces at the edge: each foot lies ahead of the edge
//! point along its profile, toward the other support's ball side. Nothing
//! here reads a B-rep fillet module or a fillet definition.
//!
//! The path is the circle of radius `C = rho_c` at station `h_c`, so the closed
//! forms of a circle hold: the nearest station is the point's meridian, as in
//! OpenCASCADE's `Extrema_ExtPElC`, and the separation from the balls has the
//! sign of the ring-torus polynomial for `C > r`. A spindle (`C <= r`) is
//! refused; none of the step 3 files needs one.
//!
//! The fan is the Voronoi zone of the centre circle in the solid of ball
//! centres, which in the meridian half-plane is the sector at the centre
//! between the rays to the two feet. The region's section is the sector less
//! the ball, on the ball side of both supports. With straight profiles that
//! set is the quadrilateral of the centre, the feet and the edge point, less
//! the ball. A sphere the ball lies inside bounds a disc, so the set stays
//! bounded. A sphere the ball lies outside leaves the sector open behind the
//! sphere, so it adds its chord wall: the line through the edge point and its
//! foot, on the centre's side. Its arc from the foot to the edge bulges toward
//! the centre, so the region lies on that side and the far zone does not.

use super::BITS;
use super::elementary::{balls, local, point_balls, position, positive, reach, vector};
use super::path::{Cover, Frame, Path};
use super::revolved::{angle, basis_of, chord_cube, hull_points, panels, root};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

/// A support of a [`Coaxial`] path in the coordinates of the shared axis: `h`
/// is the station along the unit axis from its origin.
#[derive(Clone, Debug)]
pub enum Profile {
    /// A cylinder of positive `radius` about the axis.
    Cylinder {
        /// Positive radius in mm.
        radius: Real,
    },
    /// The cone nappe whose radius at station `h` is `slope (h - apex)`, where
    /// that value is positive.
    Cone {
        /// The apex's station.
        apex: Real,
        /// Nonzero radius gained per unit of station.
        slope: Real,
    },
    /// A sphere of positive `radius` centred on the axis.
    Sphere {
        /// The centre's station.
        centre: Real,
        /// Positive radius in mm.
        radius: Real,
    },
}

/// A point of the meridian half-plane: distance from the axis, then station.
type Pair = [Real; 2];

fn sub2(a: &Pair, b: &Pair) -> Pair {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn dot2(a: &Pair, b: &Pair) -> Real {
    &(&a[0] * &b[0]) + &(&a[1] * &b[1])
}
fn cross2(a: &Pair, b: &Pair) -> Real {
    &(&a[0] * &b[1]) - &(&a[1] * &b[0])
}
fn scale2(t: &Real, a: &Pair) -> Pair {
    [t * &a[0], t * &a[1]]
}
fn add2(a: &Pair, b: &Pair) -> Pair {
    [&a[0] + &b[0], &a[1] + &b[1]]
}
fn side(outside: bool) -> Real {
    Real::from_i64(if outside { 1 } else { -1 })
}
fn oriented(value: &Real, detail: &str, dec: &mut Decider) -> Result<Real, Refusal> {
    match sign(value, dec)? {
        Sign::Positive => Ok(Real::one()),
        Sign::Negative => Ok(Real::from_i64(-1)),
        Sign::Zero => Err(reach(detail)),
    }
}

/// A support's profile in the meridian, with its ball side.
#[derive(Clone, Debug)]
enum Trace {
    /// The line `normal . p = level`, with the unit `normal` toward the ball.
    Line { normal: Pair, level: Real },
    /// The circle about `(0, centre)`; `outside` when the ball lies outside.
    Circle {
        centre: Real,
        radius: Real,
        outside: bool,
    },
}
impl Trace {
    fn new(profile: &Profile, outside: bool, dec: &mut Decider) -> Result<Self, Refusal> {
        let b = side(outside);
        Ok(match profile {
            Profile::Cylinder { radius } => {
                if !positive(radius, dec)? {
                    return Err(reach("a cylinder radius must be positive"));
                }
                Trace::Line {
                    normal: [b.clone(), Real::zero()],
                    level: &b * radius,
                }
            }
            Profile::Cone { apex, slope } => {
                if sign(slope, dec)? == Sign::Zero {
                    return Err(reach("a cone slope must be nonzero"));
                }
                // rho - slope h + slope apex = 0, outward normal (1, -slope).
                let growth = root(&(&Real::one() + &slope.square()));
                let normal = [&b / &growth, -&(&(&b * slope) / &growth)];
                let level = -&(&(&(&b * slope) * apex) / &growth);
                Trace::Line { normal, level }
            }
            Profile::Sphere { centre, radius } => {
                if !positive(radius, dec)? {
                    return Err(reach("a sphere radius must be positive"));
                }
                Trace::Circle {
                    centre: centre.clone(),
                    radius: radius.clone(),
                    outside,
                }
            }
        })
    }
    /// Signed value, positive on the ball side; a distance for a line.
    fn value(&self, p: &Pair) -> Real {
        match self {
            Trace::Line { normal, level } => &dot2(normal, p) - level,
            Trace::Circle {
                centre,
                radius,
                outside,
            } => {
                let d = &(&p[0].square() + &(&p[1] - centre).square()) - &radius.square();
                &side(*outside) * &d
            }
        }
    }
    /// The unit normal toward the ball at a point of the profile.
    fn normal(&self, p: &Pair) -> Pair {
        match self {
            Trace::Line { normal, .. } => normal.clone(),
            Trace::Circle {
                centre,
                radius,
                outside,
            } => {
                let t = &side(*outside) / radius;
                scale2(&t, &[p[0].clone(), &p[1] - centre])
            }
        }
    }
    /// The foot on the profile of a ball centre on its offset.
    fn foot(&self, c: &Pair, r: &Real) -> Pair {
        match self {
            Trace::Line { normal, .. } => sub2(c, &scale2(r, normal)),
            Trace::Circle { centre, radius, .. } => {
                let o = [Real::zero(), centre.clone()];
                let v = sub2(c, &o);
                let reach = root(&dot2(&v, &v));
                add2(&o, &scale2(&(radius / &reach), &v))
            }
        }
    }
}

/// Ball centres where the two offset profiles meet in the half-plane.
fn meetings(traces: &[Trace; 2], r: &Real, dec: &mut Decider) -> Result<Vec<Pair>, Refusal> {
    let line_circle = |normal: &Pair, level: &Real, centre: &Real, radius: &Real| {
        let offset = level + r;
        let tangent = [-&normal[1], normal[0].clone()];
        let q = sub2(&scale2(&offset, normal), &[Real::zero(), centre.clone()]);
        let beta = dot2(&tangent, &q);
        let disc = &beta.square() - &(&dot2(&q, &q) - &radius.square());
        (offset, tangent, q, beta, disc)
    };
    let mut found = Vec::new();
    match traces {
        [
            Trace::Line {
                normal: a,
                level: la,
            },
            Trace::Line {
                normal: b,
                level: lb,
            },
        ] => {
            let det = cross2(a, b);
            if sign(&det, dec)? == Sign::Zero {
                return Err(reach("parallel profiles have no rolling-ball circle"));
            }
            let (ka, kb) = (la + r, lb + r);
            found.push([
                &(&(&ka * &b[1]) - &(&kb * &a[1])) / &det,
                &(&(&a[0] * &kb) - &(&b[0] * &ka)) / &det,
            ]);
        }
        [Trace::Line { normal, level }, circle @ Trace::Circle { .. }]
        | [circle @ Trace::Circle { .. }, Trace::Line { normal, level }] => {
            let Trace::Circle {
                centre,
                radius,
                outside,
            } = circle
            else {
                unreachable!("a circle trace")
            };
            let big = radius + &(&side(*outside) * r);
            if !positive(&big, dec)? {
                return Err(reach(
                    "the ball is no smaller than the sphere it rolls inside",
                ));
            }
            let (offset, tangent, _, beta, disc) = line_circle(normal, level, centre, &big);
            if !positive(&disc, dec)? {
                return Err(reach("the offset profiles do not cross"));
            }
            let base = scale2(&offset, normal);
            for t in [&-&beta - &root(&disc), &-&beta + &root(&disc)] {
                found.push(add2(&base, &scale2(&t, &tangent)));
            }
        }
        [
            Trace::Circle {
                centre: h1,
                radius: r1,
                outside: o1,
            },
            Trace::Circle {
                centre: h2,
                radius: r2,
                outside: o2,
            },
        ] => {
            let big = [r1 + &(&side(*o1) * r), r2 + &(&side(*o2) * r)];
            if !positive(&big[0], dec)? || !positive(&big[1], dec)? {
                return Err(reach(
                    "the ball is no smaller than the sphere it rolls inside",
                ));
            }
            let gap = h2 - h1;
            if sign(&gap, dec)? == Sign::Zero {
                return Err(reach("concentric spheres have no rolling-ball circle"));
            }
            let h = &(&(h1 + h2) / &Real::from_i64(2))
                + &(&(&big[0].square() - &big[1].square()) / &(&Real::from_i64(2) * &gap));
            let rho2 = &big[0].square() - &(&h - h1).square();
            if !positive(&rho2, dec)? {
                return Err(reach("the offset profiles do not cross"));
            }
            found.push([root(&rho2), h]);
        }
    }
    Ok(found)
}

/// The rolling-ball circle between two coaxial cylinders, cone nappes or
/// spheres.
#[derive(Clone, Debug)]
pub struct Coaxial {
    point: Point3,
    axis: Vector3,
    basis: [Vector3; 2],
    traces: [Trace; 2],
    radius: Real,
    centre: Pair,
    feet: [Pair; 2],
    edge: Pair,
    // The sector between the rays to the feet is where both are non-negative.
    sense: Real,
    // Chord walls of spheres the ball lies outside of: the foot and the
    // orientation that is positive at the centre.
    chords: Vec<(Pair, Real)>,
    // Spheres, whose arcs need a chord cube in the box.
    spheres: Vec<usize>,
}

impl Coaxial {
    /// Derives the centre circle from the two profiles about the axis through
    /// `point` along `axis`, the ball side of each, the edge point `(rho, h)`
    /// in the meridian, and the radius.
    pub fn new(
        point: Point3,
        axis: Vector3,
        profiles: [Profile; 2],
        outside: [bool; 2],
        edge: Pair,
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? {
            return Err(reach("the rolling radius must be positive"));
        }
        if !positive(&axis.norm_squared(), dec)? {
            return Err(reach("the shared axis needs a nonzero direction"));
        }
        if !positive(&edge[0], dec)? {
            return Err(reach("the edge circle needs a positive radius"));
        }
        let axis = axis.scaled(&(&Real::one() / &root(&axis.norm_squared())));
        let traces = [
            Trace::new(&profiles[0], outside[0], dec)?,
            Trace::new(&profiles[1], outside[1], dec)?,
        ];
        for trace in &traces {
            if sign(&trace.value(&edge), dec)? != Sign::Zero {
                return Err(reach("the edge circle does not lie on both supports"));
            }
        }
        // Each face leaves the edge along its profile toward the other
        // support's ball side.
        let mut ahead = Vec::with_capacity(2);
        for (i, trace) in traces.iter().enumerate() {
            let n = trace.normal(&edge);
            let tangent = [-&n[1], n[0].clone()];
            let other = traces[1 - i].normal(&edge);
            let t = oriented(
                &dot2(&tangent, &other),
                "the supports are tangent at the edge",
                dec,
            )?;
            ahead.push(scale2(&t, &tangent));
        }
        let mut chosen = None;
        for c in meetings(&traces, &radius, dec)? {
            if !positive(&c[0], dec)? {
                continue;
            }
            let feet = [traces[0].foot(&c, &radius), traces[1].foot(&c, &radius)];
            let mut fits = true;
            for (foot, t) in feet.iter().zip(&ahead) {
                fits &= positive(&dot2(&sub2(foot, &edge), t), dec)?;
            }
            if fits {
                if chosen.is_some() {
                    return Err(reach("two balls fit the corner of the edge"));
                }
                chosen = Some((c, feet));
            }
        }
        let (centre, feet) =
            chosen.ok_or_else(|| reach("no ball of this radius fits the corner of the edge"))?;
        if !positive(&(&centre[0] - &radius), dec)? {
            return Err(reach(
                "the ball centre circle is no larger than the ball, a spindle the field does not hold",
            ));
        }
        let sense = oriented(
            &cross2(&sub2(&feet[0], &centre), &sub2(&feet[1], &centre)),
            "the ball touches both supports at one point",
            dec,
        )?;
        let mut chords = Vec::new();
        let mut spheres = Vec::new();
        for (i, trace) in traces.iter().enumerate() {
            if let Trace::Circle { outside, .. } = trace {
                spheres.push(i);
                if *outside {
                    let chord = sub2(&feet[i], &edge);
                    let t = oriented(
                        &cross2(&chord, &sub2(&centre, &edge)),
                        "the centre lies on a contact chord",
                        dec,
                    )?;
                    chords.push((feet[i].clone(), t));
                }
            }
        }
        let basis = basis_of(&axis, dec)?;
        Ok(Self {
            point,
            axis,
            basis,
            traces,
            radius,
            centre,
            feet,
            edge,
            sense,
            chords,
            spheres,
        })
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    /// Applies an exact rigid motion, reflections included. Meridian
    /// coordinates are invariant.
    pub fn transform(&mut self, pose: &Isometry) {
        self.point = pose.apply_point(&self.point);
        self.axis = pose.apply_vector(&self.axis);
        for e in &mut self.basis {
            *e = pose.apply_vector(e);
        }
    }

    /// Meridian coordinates of a point.
    fn meridian(&self, x: &Vector3) -> Pair {
        let v = x.minus(&position(&self.point));
        let h = self.axis.dot(&v);
        let rho2 = &v.norm_squared() - &h.square();
        [root(&rho2), h]
    }

    fn revolve(&self, p: &Pair, direction: &Vector3) -> Point3 {
        self.point
            .plus(&self.axis.scaled(&p[1]))
            .plus(&direction.scaled(&p[0]))
    }

    fn frame_at(&self, direction: &Vector3, dec: &mut Decider) -> Result<Frame, Refusal> {
        Ok(Frame {
            centre: point_balls(&self.revolve(&self.centre, direction), dec)?,
            tangent: balls(&self.axis.cross(direction), dec)?,
            feet: [
                point_balls(&self.revolve(&self.feet[0], direction), dec)?,
                point_balls(&self.revolve(&self.feet[1], direction), dec)?,
            ],
        })
    }

    /// An axis-aligned box holding the region: the hull, over a cover of the
    /// circle by 64 panels of `pi / 32`, of the edge, both feet and the
    /// centre, and of the disc on each sphere's chord from the edge to its
    /// foot. The region's section lies in the hull of its edge point, feet and
    /// centre, except for a sphere's arc, a minor arc of a great circle, which
    /// lies in that disc. Straight profiles add nothing.
    pub fn hull(&self, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        Ok(hull_points(&self.boxes(dec)?))
    }

    /// The box of the region over each panel of [`Coaxial::hull`].
    pub fn boxes(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let cover = self.cover(dec)?;
        let mut boxes = Vec::new();
        for u in panels() {
            let (_, _, direction) = angle(&cover.pi, &cover.basis, &u);
            let edge = cover.revolve(&cover.edge, &direction);
            let feet = cover.feet.each_ref().map(|f| cover.revolve(f, &direction));
            let mut points = Vec::new();
            for &s in &self.spheres {
                points.extend(chord_cube(&edge, &feet[s])?);
            }
            points.push(edge);
            points.push(cover.revolve(&cover.centre, &direction));
            points.extend(feet);
            boxes.push(hull_points(&points));
        }
        Ok(boxes)
    }

    /// Tentative binary64 value for display: negative inside the region. The
    /// square root is the caller's.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let f = |x: &Real| x.to_f64();
        let o = [f(&self.point.x), f(&self.point.y), f(&self.point.z)];
        let a = [f(&self.axis.x), f(&self.axis.y), f(&self.axis.z)];
        let v = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
        let h = v[0] * a[0] + v[1] * a[1] + v[2] * a[2];
        let rho = sqrt((v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - h * h).max(0.0));
        let x = [rho, h];
        let pair = |q: &Pair| [f(&q[0]), f(&q[1])];
        let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
        let minus = |a: [f64; 2], b: [f64; 2]| [a[0] - b[0], a[1] - b[1]];
        let (c, e) = (pair(&self.centre), pair(&self.edge));
        let feet = self.feet.each_ref().map(pair);
        let support = |trace: &Trace| match trace {
            Trace::Line { normal, level } => {
                let n = pair(normal);
                n[0] * x[0] + n[1] * x[1] - f(level)
            }
            Trace::Circle {
                centre,
                radius,
                outside,
            } => {
                let d = sqrt(x[0] * x[0] + (x[1] - f(centre)).powi(2)) - f(radius);
                if *outside { d } else { -d }
            }
        };
        let s = f(&self.sense);
        let mut walls = vec![
            s * cross(minus(feet[0], c), minus(x, c)),
            s * cross(minus(x, c), minus(feet[1], c)),
        ];
        for (foot, t) in &self.chords {
            walls.push(f(t) * cross(minus(pair(foot), e), minus(x, e)));
        }
        let d = minus(x, c);
        let m = sqrt(d[0] * d[0] + d[1] * d[1]);
        let r = f(&self.radius);
        let [first, second] = self.traces.each_ref().map(support);
        if first >= 0.0 && second >= 0.0 && walls.iter().all(|w| *w >= 0.0) {
            (-first).max(-second).max(r - m)
        } else {
            (m - r).abs().max(f64::MIN_POSITIVE)
        }
    }
}

/// Interval cover of a [`Coaxial`] path; `u` in `[-1, 1]` is the angle `pi u`
/// from the first basis direction.
#[derive(Clone, Debug)]
pub struct CoaxialCover {
    point: [Ball; 3],
    axis: [Ball; 3],
    radius: Ball,
    pi: Ball,
    basis: [[Ball; 3]; 2],
    centre: [Ball; 2],
    feet: [[Ball; 2]; 2],
    edge: [Ball; 2],
}
impl CoaxialCover {
    fn revolve(&self, p: &[Ball; 2], direction: &[Ball; 3]) -> [Ball; 3] {
        std::array::from_fn(|i| {
            &(&self.point[i] + &(&p[1] * &self.axis[i])) + &(&p[0] * &direction[i])
        })
    }
}
impl Cover for CoaxialCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let (_, _, direction) = angle(&self.pi, &self.basis, u);
        Ok(self.revolve(&self.centre, &direction))
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let (cos, sin, direction) = angle(&self.pi, &self.basis, u);
        let turn =
            std::array::from_fn(|i| &(&cos * &self.basis[1][i]) - &(&sin * &self.basis[0][i]));
        Ok(Frame {
            centre: self.revolve(&self.centre, &direction),
            tangent: turn,
            feet: self.feet.each_ref().map(|f| self.revolve(f, &direction)),
        })
    }
}

impl Path for Coaxial {
    type Cover = CoaxialCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = self.meridian(&vector(local));
        let c = &self.centre;
        let mut walls = vec![
            self.traces[0].value(&x),
            self.traces[1].value(&x),
            &self.sense * &cross2(&sub2(&self.feet[0], c), &sub2(&x, c)),
            &self.sense * &cross2(&sub2(&x, c), &sub2(&self.feet[1], c)),
        ];
        for (foot, t) in &self.chords {
            walls.push(t * &cross2(&sub2(foot, &self.edge), &sub2(&x, &self.edge)));
        }
        Ok(walls)
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn uniform(&self) -> bool {
        true
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        let centre = self.point.plus(&self.axis.scaled(&self.centre[1]));
        let v = vector(local).minus(&position(&centre));
        let h = self.axis.dot(&v);
        let d2 = v.norm_squared();
        let rho2 = &d2 - &h.square();
        let major_sq = self.centre[0].square();
        let s = &(&d2 + &major_sq) - &self.radius.square();
        Some(&s.square() - &(&(&Real::from_i64(4) * &major_sq) * &rho2))
    }
    fn cover(&self, dec: &mut Decider) -> Result<CoaxialCover, Refusal> {
        let pair = |p: &Pair, dec: &mut Decider| -> Result<[Ball; 2], Refusal> {
            Ok([enclose(&p[0], BITS, dec)?, enclose(&p[1], BITS, dec)?])
        };
        Ok(CoaxialCover {
            point: point_balls(&self.point, dec)?,
            axis: balls(&self.axis, dec)?,
            radius: enclose(&self.radius, BITS, dec)?,
            pi: Ball::pi_checked(BITS, dec)?,
            basis: [balls(&self.basis[0], dec)?, balls(&self.basis[1], dec)?],
            centre: pair(&self.centre, dec)?,
            feet: [pair(&self.feet[0], dec)?, pair(&self.feet[1], dec)?],
            edge: pair(&self.edge, dec)?,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let v = vector(local).minus(&position(&self.point));
        let w = v.minus(&self.axis.scaled(&self.axis.dot(&v)));
        let rho2 = w.norm_squared();
        // On the axis every meridian is nearest; take the first basis direction.
        let direction = if sign(&rho2, dec)? == Sign::Zero {
            self.basis[0].clone()
        } else {
            w.scaled(&(&Real::one() / &root(&rho2)))
        };
        self.frame_at(&direction, dec).map(Some)
    }
}
