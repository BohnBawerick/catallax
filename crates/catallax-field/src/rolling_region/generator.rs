//! Straight ball paths between supports parallel to one line (4.2a step 2).
//!
//! Planes and cylinders whose axes are parallel to a line `d` offset into
//! planes and cylinders parallel to `d`, so their offsets meet in lines along
//! `d`: two planes in one line, a plane and a cylinder in up to two, two
//! cylinders in up to two. The field solves the offsets in the section plane
//! through an edge point and takes the root on the edge's side: for a plane
//! and a cylinder, the first meeting of the plane offset with the cylinder
//! offset when moving from the edge along the plane toward the cylinder's ball
//! side; for two cylinders, the meeting of the offset circles on the edge's
//! side of the line of centres. Nothing here reads a B-rep fillet module or a
//! fillet definition.
//!
//! The normal wall through each contact rail is a plane: the normals of a
//! plane, or of a cylinder along one generator, are parallel. The fan lies on
//! the ball side of both supports, on the edge side of both rail planes and
//! between the end planes, which must be normal to the line. The bands of
//! [`super::revolved`] keep far zones out. The rail planes meet in an
//! unbounded wedge that crosses a cylinder the ball lies outside of and leaves
//! it on the far side, so such a cylinder adds its chord wall: the plane
//! through the edge and the cylinder's rail, on the centre's side. The section
//! of the region lies there: it is bounded by the cylinder's minor arc from
//! the edge to the rail, outside the cylinder, and that arc lies on the far
//! side of its chord from the axis, as the centre does, beyond the rail. The closed forms are those of
//! step 1's strip: the projection onto the line, clamped to the run, and the
//! squared distance to the line less `r^2`.

use super::elementary::{
    Flat, balls, length, local, point_balls, position, positive, reach, shift, vector,
};
use super::path::{Cover, Frame, Path};
use super::revolved::{chord_cube, hull_points, root};
use super::{BITS, interval, q};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

/// One support of a [`Generator`] path, parallel to the path's line.
#[derive(Clone, Debug)]
pub enum Wall {
    /// A plane whose normal points to the ball.
    Plane(Flat),
    /// A cylinder whose axis is parallel to the line.
    Cylinder {
        /// Any point of the axis.
        point: Point3,
        /// Nonzero axis direction.
        axis: Vector3,
        /// Positive radius in mm.
        radius: Real,
        /// Whether the ball lies outside the cylinder.
        outside: bool,
    },
}

/// The straight rolling-ball path between two supports parallel to a line,
/// over the stations `ends[0] <= t <= ends[1]` of `c + t d`.
#[derive(Clone, Debug)]
pub struct Generator {
    supports: [Wall; 2],
    radius: Real,
    direction: Vector3,
    squared: Real,
    centre: Point3,
    edge: Point3,
    feet: [Point3; 2],
    rails: [Vector3; 2],
    // Chord walls of the cylinders the ball lies outside of: a rail point
    // and the oriented normal.
    chords: Vec<(Point3, Vector3)>,
    ends: [Real; 2],
}

fn side(outside: bool) -> Real {
    Real::from_i64(if outside { 1 } else { -1 })
}
/// The component of `v` perpendicular to `d`.
fn across(v: &Vector3, d: &Vector3, squared: &Real) -> Vector3 {
    v.minus(&d.scaled(&(&v.dot(d) / squared)))
}

impl Generator {
    /// Derives the centre line from the two supports offset by `radius`
    /// toward the ball. `edge` is a point of the sharp edge, `direction` its
    /// line, and `ends` the stations of the run's end planes along it, which
    /// the caller has proven normal to the line.
    pub fn new(
        supports: [Wall; 2],
        edge: Point3,
        direction: Vector3,
        ends: [Real; 2],
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? {
            return Err(reach("the rolling radius must be positive"));
        }
        let squared = direction.norm_squared();
        if !positive(&squared, dec)? {
            return Err(reach("a generator needs a nonzero direction"));
        }
        if !positive(&(&ends[1] - &ends[0]), dec)? {
            return Err(reach("a generator needs increasing end stations"));
        }
        for wall in &supports {
            match wall {
                Wall::Plane(flat) => {
                    if !positive(&flat.normal.norm_squared(), dec)?
                        || sign(&flat.normal.dot(&direction), dec)? != Sign::Zero
                    {
                        return Err(reach("a generator's plane must contain its direction"));
                    }
                }
                Wall::Cylinder { axis, radius, .. } => {
                    if !positive(radius, dec)?
                        || sign(&axis.cross(&direction).norm_squared(), dec)? != Sign::Zero
                    {
                        return Err(reach(
                            "a generator's cylinder must have a positive radius and a parallel axis",
                        ));
                    }
                }
            }
        }
        let centre = match &supports {
            [Wall::Plane(a), Wall::Plane(b)] => planes(a, b, &edge, &direction, &radius, dec)?,
            [Wall::Plane(a), cylinder @ Wall::Cylinder { .. }]
            | [cylinder @ Wall::Cylinder { .. }, Wall::Plane(a)] => {
                mixed(a, cylinder, &edge, &direction, &squared, &radius, dec)?
            }
            [a, b] => cylinders(a, b, &edge, &direction, &squared, &radius, dec)?,
        };
        let feet = supports.each_ref().map(|wall| match wall {
            Wall::Plane(flat) => {
                centre.plus(&flat.normal.scaled(&-&(&radius / &length(&flat.normal))))
            }
            Wall::Cylinder {
                point,
                radius: big,
                outside,
                ..
            } => {
                let radial = across(&centre.minus(point), &direction, &squared);
                let offset = big + &(&side(*outside) * &radius);
                centre.plus(&radial.scaled(&-&(&(&side(*outside) * &radius) / &offset)))
            }
        });
        let mut rails = Vec::with_capacity(2);
        for foot in &feet {
            let normal = direction.cross(&centre.minus(foot));
            match sign(&edge.minus(&centre).dot(&normal), dec)? {
                Sign::Positive => rails.push(normal),
                Sign::Negative => rails.push(normal.scaled(&Real::from_i64(-1))),
                Sign::Zero => {
                    return Err(reach("the edge lies on a contact rail's normal wall"));
                }
            }
        }
        let rails: [Vector3; 2] = rails.try_into().expect("two rails");
        let mut chords = Vec::new();
        for (wall, foot) in supports.iter().zip(&feet) {
            if let Wall::Cylinder { outside: true, .. } = wall {
                let normal = direction.cross(&foot.minus(&edge));
                match sign(&centre.minus(&edge).dot(&normal), dec)? {
                    Sign::Positive => chords.push((foot.clone(), normal)),
                    Sign::Negative => {
                        chords.push((foot.clone(), normal.scaled(&Real::from_i64(-1))));
                    }
                    Sign::Zero => return Err(reach("the centre lies on a contact chord")),
                }
            }
        }
        Ok(Self {
            supports,
            radius,
            direction,
            squared,
            centre,
            edge,
            feet,
            rails,
            chords,
            ends,
        })
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    /// Applies an exact rigid motion, reflections included. Stations are
    /// relative to the moved centre and direction, so they are unchanged.
    pub fn transform(&mut self, pose: &Isometry) {
        let v = |v: &mut Vector3| *v = pose.apply_vector(v);
        let p = |p: &mut Point3| *p = pose.apply_point(p);
        for wall in &mut self.supports {
            match wall {
                Wall::Plane(flat) => {
                    p(&mut flat.point);
                    v(&mut flat.normal);
                }
                Wall::Cylinder { point, axis, .. } => {
                    p(point);
                    v(axis);
                }
            }
        }
        v(&mut self.direction);
        p(&mut self.centre);
        p(&mut self.edge);
        for foot in &mut self.feet {
            p(foot);
        }
        for rail in &mut self.rails {
            v(rail);
        }
        for (foot, normal) in &mut self.chords {
            p(foot);
            v(normal);
        }
    }

    fn station(&self, t: &Real) -> Vector3 {
        self.direction.scaled(t)
    }

    /// An axis-aligned box holding the region: the hull of the edge, both
    /// rails and the centre line at the two end stations, and, for each
    /// cylinder, of the disc on the chord from the edge to its rail at both
    /// stations. Each section lies in the hull of its edge point, feet and
    /// centre, except for a cylinder's minor arc, which lies in that disc.
    pub fn hull(&self, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        let mut points = Vec::new();
        for t in &self.ends {
            let step = self.station(t);
            let edge = point_balls(&self.edge.plus(&step), dec)?;
            let feet = [
                point_balls(&self.feet[0].plus(&step), dec)?,
                point_balls(&self.feet[1].plus(&step), dec)?,
            ];
            for (wall, foot) in self.supports.iter().zip(&feet) {
                if matches!(wall, Wall::Cylinder { .. }) {
                    points.extend(chord_cube(&edge, foot)?);
                }
            }
            points.push(edge);
            points.push(point_balls(&self.centre.plus(&step), dec)?);
            points.extend(feet);
        }
        Ok(hull_points(&points))
    }

    /// Tentative binary64 value for display: negative inside the region. The
    /// square root is the caller's.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let fv = |v: &Vector3| [v.x.to_f64(), v.y.to_f64(), v.z.to_f64()];
        let fp = |p: &Point3| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let d = fv(&self.direction);
        let dd = dot(d, d);
        let r = self.radius.to_f64();
        let v = sub(p, fp(&self.centre));
        let t = dot(v, d) / dd;
        let ends = [self.ends[0].to_f64(), self.ends[1].to_f64()];
        let clamped = t.clamp(ends[0], ends[1]);
        let m = sqrt(dot(v, v) - 2.0 * clamped * dot(v, d) + clamped * clamped * dd);
        let distance = |wall: &Wall| match wall {
            Wall::Plane(flat) => {
                let n = fv(&flat.normal);
                dot(n, sub(p, fp(&flat.point))) / sqrt(dot(n, n))
            }
            Wall::Cylinder {
                point,
                radius,
                outside,
                ..
            } => {
                let w = sub(p, fp(point));
                let s = dot(w, d) / dd;
                let rho = sqrt(dot(w, w) - s * s * dd);
                if *outside {
                    rho - radius.to_f64()
                } else {
                    radius.to_f64() - rho
                }
            }
        };
        let mut walls = vec![
            dot(v, fv(&self.rails[0])),
            dot(v, fv(&self.rails[1])),
            t - ends[0],
            ends[1] - t,
        ];
        for (foot, normal) in &self.chords {
            walls.push(dot(sub(p, fp(foot)), fv(normal)));
        }
        if walls.iter().all(|w| *w >= 0.0) {
            (-distance(&self.supports[0]))
                .max(-distance(&self.supports[1]))
                .max(r - m)
        } else {
            (m - r).abs().max(f64::MIN_POSITIVE)
        }
    }
}

/// Cramer's rule for the two offset planes and the section plane through the
/// edge point; its determinant is `n_a . (n_b x d)`.
fn planes(
    a: &Flat,
    b: &Flat,
    edge: &Point3,
    d: &Vector3,
    r: &Real,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let det = a.normal.dot(&b.normal.cross(d));
    if sign(&det, dec)? == Sign::Zero {
        return Err(reach("parallel planes have no rolling-ball strip"));
    }
    let level = |f: &Flat| &f.normal.dot(&position(&f.point)) + &(r * &length(&f.normal));
    let x = b
        .normal
        .cross(d)
        .scaled(&level(a))
        .plus(&d.cross(&a.normal).scaled(&level(b)))
        .plus(&a.normal.cross(&b.normal).scaled(&d.dot(&position(edge))))
        .scaled(&(&Real::one() / &det));
    Ok(Point3::origin().plus(&x))
}

/// The plane offset meets the cylinder offset along `x(t) = p1 + t m`, with
/// `p1` the edge point lifted by `r` along the plane's unit normal and
/// `m = d x n` in the section. The root taken is the first meeting when
/// moving from the edge along the plane toward the cylinder's ball side.
fn mixed(
    plane: &Flat,
    cylinder: &Wall,
    edge: &Point3,
    d: &Vector3,
    squared: &Real,
    r: &Real,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let Wall::Cylinder {
        point,
        radius,
        outside,
        ..
    } = cylinder
    else {
        unreachable!("a cylinder wall")
    };
    let b = side(*outside);
    let offset = radius + &(&b * r);
    if !positive(&offset, dec)? {
        return Err(reach(
            "the ball is no smaller than the cylinder it rolls inside",
        ));
    }
    let lifted = edge.plus(&plane.normal.scaled(&(r / &length(&plane.normal))));
    let m = d.cross(&plane.normal);
    let w1 = lifted.minus(point).cross(d);
    let w2 = m.cross(d);
    let a2 = w2.norm_squared();
    let a1 = w1.dot(&w2);
    let a0 = &w1.norm_squared() - &(&offset.square() * squared);
    let disc = &a1.square() - &(&a2 * &a0);
    if !positive(&disc, dec)? {
        return Err(reach("the plane offset does not cross the cylinder offset"));
    }
    let toward = across(&edge.minus(point), d, squared).scaled(&b);
    let ahead = sign(&m.dot(&toward), dec)?;
    let roots = [&(&-&a1 - &root(&disc)) / &a2, &(&-&a1 + &root(&disc)) / &a2];
    let ordered: Vec<&Real> = match ahead {
        Sign::Positive => roots.iter().collect(),
        Sign::Negative => roots.iter().rev().collect(),
        Sign::Zero => return Err(reach("the plane is tangent to the cylinder at the edge")),
    };
    for t in ordered {
        if sign(t, dec)? == ahead {
            return Ok(lifted.plus(&m.scaled(t)));
        }
    }
    Err(reach(
        "the plane offset meets the cylinder offset only behind the edge",
    ))
}

/// The offset circles in the section meet on the edge's side of the line of
/// centres: `q1 + alpha D + beta (d x D)` with `|alpha D + beta (d x D)| = C1`
/// and `|(alpha - 1) D + beta (d x D)| = C2`.
fn cylinders(
    a: &Wall,
    b: &Wall,
    edge: &Point3,
    d: &Vector3,
    squared: &Real,
    r: &Real,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let foot = |wall: &Wall| -> (Point3, Real) {
        let Wall::Cylinder {
            point,
            radius,
            outside,
            ..
        } = wall
        else {
            unreachable!("a cylinder wall")
        };
        let shift = d.scaled(&(&edge.minus(point).dot(d) / squared));
        (point.plus(&shift), radius + &(&side(*outside) * r))
    };
    let (first, c1) = foot(a);
    let (second, c2) = foot(b);
    if !positive(&c1, dec)? || !positive(&c2, dec)? {
        return Err(reach(
            "the ball is no smaller than the cylinder it rolls inside",
        ));
    }
    let delta = second.minus(&first);
    let l2 = delta.norm_squared();
    if !positive(&l2, dec)? {
        return Err(reach(
            "coaxial cylinders have no straight rolling-ball strip",
        ));
    }
    let alpha = &(&(&c1.square() - &c2.square()) + &l2) / &(&Real::from_i64(2) * &l2);
    let beta2 = &(&c1.square() - &(&alpha.square() * &l2)) / &(squared * &l2);
    if !positive(&beta2, dec)? {
        return Err(reach("the cylinder offsets do not meet"));
    }
    let normal = d.cross(&delta);
    let beta = match sign(&edge.minus(&first).dot(&normal), dec)? {
        Sign::Positive => root(&beta2),
        Sign::Negative => -root(&beta2),
        Sign::Zero => return Err(reach("the edge lies on the line of centres")),
    };
    Ok(first
        .plus(&delta.scaled(&alpha))
        .plus(&normal.scaled(&beta)))
}

/// Interval cover of a generator; `u` in `[-1, 1]` maps linearly onto the run.
#[derive(Clone, Debug)]
pub struct GeneratorCover {
    centre: [Ball; 3],
    direction: [Ball; 3],
    feet: [[Ball; 3]; 2],
    ends: [Ball; 2],
    radius: Ball,
}
impl GeneratorCover {
    fn station(&self, u: &Ball) -> Ball {
        let one = interval(&Rational::one(), &Rational::one());
        let half = Rational::one() / q(2);
        let step = &(&self.ends[1] - &self.ends[0]) * &interval(&half, &half);
        &self.ends[0] + &(&(u + &one) * &step)
    }
}
impl Cover for GeneratorCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        Ok(shift(&self.centre, &self.station(u), &self.direction))
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let t = self.station(u);
        Ok(Frame {
            centre: shift(&self.centre, &t, &self.direction),
            tangent: self.direction.clone(),
            feet: self.feet.each_ref().map(|f| shift(f, &t, &self.direction)),
        })
    }
}

impl Path for Generator {
    type Cover = GeneratorCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = vector(local);
        let two = Real::from_i64(2);
        let r = &self.radius;
        let d = &self.direction;
        let mut walls = Vec::with_capacity(10);
        for wall in &self.supports {
            match wall {
                Wall::Plane(flat) => {
                    let height = flat.normal.dot(&x.minus(&position(&flat.point)));
                    let band = &(&(&two * r) * &length(&flat.normal)) - &height;
                    walls.push(height);
                    walls.push(band);
                }
                Wall::Cylinder {
                    point,
                    radius,
                    outside,
                    ..
                } => {
                    let spread = x.minus(&position(point)).cross(d).norm_squared();
                    let rim = &radius.square() * &self.squared;
                    walls.push(&side(*outside) * &(&spread - &rim));
                    if *outside {
                        walls.push(&(&(radius + &(&two * r)).square() * &self.squared) - &spread);
                    }
                }
            }
        }
        let v = x.minus(&position(&self.centre));
        let along = d.dot(&v);
        walls.push(v.dot(&self.rails[0]));
        walls.push(v.dot(&self.rails[1]));
        for (foot, normal) in &self.chords {
            walls.push(x.minus(&position(foot)).dot(normal));
        }
        walls.push(&along - &(&self.ends[0] * &self.squared));
        walls.push(&(&self.ends[1] * &self.squared) - &along);
        Ok(walls)
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn uniform(&self) -> bool {
        true
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        let v = vector(local).minus(&position(&self.centre));
        let along = self.direction.dot(&v);
        Some(&(&v.norm_squared() - &(&along.square() / &self.squared)) - &self.radius.square())
    }
    fn cover(&self, dec: &mut Decider) -> Result<GeneratorCover, Refusal> {
        Ok(GeneratorCover {
            centre: point_balls(&self.centre, dec)?,
            direction: balls(&self.direction, dec)?,
            feet: [
                point_balls(&self.feet[0], dec)?,
                point_balls(&self.feet[1], dec)?,
            ],
            ends: [
                enclose(&self.ends[0], BITS, dec)?,
                enclose(&self.ends[1], BITS, dec)?,
            ],
            radius: enclose(&self.radius, BITS, dec)?,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let v = vector(local).minus(&position(&self.centre));
        let mut t = &self.direction.dot(&v) / &self.squared;
        if sign(&(&t - &self.ends[0]), dec)? == Sign::Negative {
            t = self.ends[0].clone();
        } else if sign(&(&t - &self.ends[1]), dec)? == Sign::Positive {
            t = self.ends[1].clone();
        }
        let step = self.station(&t);
        Ok(Some(Frame {
            centre: point_balls(&self.centre.plus(&step), dec)?,
            tangent: balls(&self.direction, dec)?,
            feet: [
                point_balls(&self.feet[0].plus(&step), dec)?,
                point_balls(&self.feet[1].plus(&step), dec)?,
            ],
        }))
    }
}
