//! Line and circle ball paths of the elementary fillets (4.2a step 1).
//!
//! The field derives each path from the two supports and the ball radius
//! alone: offset each support by `r` toward the ball and intersect the
//! offsets. Two planes give a line; a plane and a perpendicular cylinder give a
//! circle. Nothing here reads a B-rep fillet module or a fillet definition.
//!
//! Both paths have closed forms. The nearest station is the projection onto
//! the line, or the meridian through the point, as in OpenCASCADE's
//! `Extrema_ExtPElC`. The separation from the balls is the squared distance to
//! the line less `r^2`, or the ring-torus polynomial
//! `(|x-o|^2 + C^2 - r^2)^2 - 4 C^2 rho^2`, which has the sign of
//! `(rho-C)^2 + h^2 - r^2` when `C > r`.
//!
//! The fan is the section between the edge, the two feet and the centre,
//! swept along the path: the ball side of each support, and the edge side of
//! each normal wall through a contact rail. A strip also has two end planes
//! normal to its line. Every wall is a face of the region, so these paths
//! prove `Boundary` exactly.

use super::path::{Cover, Frame, Path};
use super::{BITS, gap, interval};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, Vector3};

/// A plane support. `normal` points to the side where the ball lies.
#[derive(Clone, Debug)]
pub struct Flat {
    /// Any point of the plane.
    pub point: Point3,
    /// Nonzero normal toward the ball.
    pub normal: Vector3,
}

/// A cylinder support.
#[derive(Clone, Debug)]
pub struct Round {
    /// Any point of the axis.
    pub point: Point3,
    /// Nonzero axis direction.
    pub axis: Vector3,
    /// Positive cylinder radius in mm.
    pub radius: Real,
    /// Whether the ball lies outside the cylinder.
    pub outside: bool,
}

pub(super) fn reach(detail: &str) -> Refusal {
    Refusal::new("field.rolling-path-reach", detail)
}
pub(super) fn length(v: &Vector3) -> Real {
    let squared = v.norm_squared();
    squared
        .as_rational()
        .and_then(|s| crate::profile::sqrt(&s))
        .map(Real::from)
        .unwrap_or_else(|| squared.sqrt())
}
pub(super) fn position(p: &Point3) -> Vector3 {
    p.minus(&Point3::origin())
}
pub(super) fn local(point: &Point3) -> [Real; 3] {
    [point.x.clone(), point.y.clone(), point.z.clone()]
}
pub(super) fn vector([x, y, z]: &[Real; 3]) -> Vector3 {
    Vector3::new(x.clone(), y.clone(), z.clone())
}
pub(super) fn balls(v: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        enclose(&v.x, BITS, dec)?,
        enclose(&v.y, BITS, dec)?,
        enclose(&v.z, BITS, dec)?,
    ])
}
pub(super) fn point_balls(p: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    balls(&position(p), dec)
}
pub(super) fn shift(a: &[Ball; 3], t: &Ball, d: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &(t * &d[i]))
}
pub(super) fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Positive)
}

/// Straight rolling-ball strip between two planes, over the stations
/// `ends[0] <= t <= ends[1]` of `c0 + t d`, where `d = n_a x n_b` and `c0` is
/// the centre in the plane `d . x = 0`.
#[derive(Clone, Debug)]
pub struct Strip {
    normals: [Vector3; 2],
    levels: [Real; 2],
    radius: Real,
    direction: Vector3,
    squared: Real,
    centre: Point3,
    feet: [Point3; 2],
    edge: Point3,
    ends: [Rational; 2],
}

impl Strip {
    /// Derives the centre line from the planes offset by `radius` toward the ball.
    pub fn new(
        supports: [Flat; 2],
        radius: Real,
        ends: [Rational; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? {
            return Err(reach("the rolling radius must be positive"));
        }
        if ends[0] >= ends[1] {
            return Err(reach("a strip needs increasing end stations"));
        }
        let [a, b] = supports;
        let direction = a.normal.cross(&b.normal);
        let squared = direction.norm_squared();
        if !positive(&squared, dec)? {
            return Err(reach("parallel planes have no rolling-ball strip"));
        }
        let levels = [&a, &b].map(|f| f.normal.dot(&position(&f.point)));
        let units = [&a, &b].map(|f| length(&f.normal));
        // Cramer's rule for n_a.x = k_a, n_b.x = k_b, d.x = 0; its determinant
        // is d.(n_a x n_b) = |d|^2.
        let solve = |ka: &Real, kb: &Real| -> Point3 {
            let x = b
                .normal
                .cross(&direction)
                .scaled(ka)
                .plus(&direction.cross(&a.normal).scaled(kb))
                .scaled(&(&Real::one() / &squared));
            Point3::origin().plus(&x)
        };
        let centre = solve(
            &(&levels[0] + &(&radius * &units[0])),
            &(&levels[1] + &(&radius * &units[1])),
        );
        let edge = solve(&levels[0], &levels[1]);
        let normals = [a.normal, b.normal];
        let feet =
            std::array::from_fn(|i| centre.plus(&normals[i].scaled(&(-&(&radius / &units[i])))));
        Ok(Self {
            normals,
            levels,
            radius,
            direction,
            squared,
            centre,
            feet,
            edge,
            ends,
        })
    }
    fn station(&self, t: &Real) -> Vector3 {
        self.direction.scaled(t)
    }
}

/// Interval cover of a strip; `u` in `[-1, 1]` maps linearly onto the stations.
#[derive(Clone, Debug)]
pub struct StripCover {
    centre: [Ball; 3],
    direction: [Ball; 3],
    feet: [[Ball; 3]; 2],
    ends: [Ball; 2],
    radius: Ball,
}
impl StripCover {
    // t = t0 + (u + 1) (t1 - t0) / 2.
    fn station(&self, u: &Ball) -> Ball {
        let one = interval(&Rational::one(), &Rational::one());
        let half = Rational::one() / Rational::from_i64(2);
        let step = &(&self.ends[1] - &self.ends[0]) * &interval(&half, &half);
        &self.ends[0] + &(&(u + &one) * &step)
    }
}
impl Cover for StripCover {
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

impl Path for Strip {
    type Cover = StripCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = vector(local);
        let along = self.direction.dot(&x.minus(&position(&self.centre)));
        let rail = |foot: &Point3| x.minus(&position(foot)).dot(&self.edge.minus(foot));
        Ok(vec![
            &self.normals[0].dot(&x) - &self.levels[0],
            &self.normals[1].dot(&x) - &self.levels[1],
            rail(&self.feet[0]),
            rail(&self.feet[1]),
            &along - &(&Real::from(self.ends[0].clone()) * &self.squared),
            &(&Real::from(self.ends[1].clone()) * &self.squared) - &along,
        ])
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
    fn cover(&self, dec: &mut Decider) -> Result<StripCover, Refusal> {
        Ok(StripCover {
            centre: point_balls(&self.centre, dec)?,
            direction: balls(&self.direction, dec)?,
            feet: [
                point_balls(&self.feet[0], dec)?,
                point_balls(&self.feet[1], dec)?,
            ],
            ends: [
                Ball::from_rational(&self.ends[0], BITS),
                Ball::from_rational(&self.ends[1], BITS),
            ],
            radius: enclose(&self.radius, BITS, dec)?,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let v = vector(local).minus(&position(&self.centre));
        let mut t = &self.direction.dot(&v) / &self.squared;
        for (end, below) in [(&self.ends[0], true), (&self.ends[1], false)] {
            let end = Real::from(end.clone());
            let s = sign(&(&t - &end), dec)?;
            if (below && s == Sign::Negative) || (!below && s == Sign::Positive) {
                t = end;
            }
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

/// Full rolling-ball ring between a plane and a perpendicular cylinder. The
/// centre circle lies in the plane offset by `r`, at radius `C = R + r` when
/// the ball is outside the cylinder and `C = R - r` inside.
#[derive(Clone, Debug)]
pub struct Ring {
    normal: Vector3,
    level: Real,
    unit: Real,
    radius: Real,
    body: Real,
    outside: bool,
    axis: Vector3,
    origin: Point3,
    major: Real,
    basis: [Vector3; 2],
}

impl Ring {
    /// Derives the centre circle. Refuses unless the plane is perpendicular to
    /// the axis and `C > r`, the reach of the circle that keeps the tube regular.
    pub fn new(
        plane: Flat,
        cylinder: Round,
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? || !positive(&cylinder.radius, dec)? {
            return Err(reach("the rolling and cylinder radii must be positive"));
        }
        if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
            || !positive(&cylinder.axis.norm_squared(), dec)?
        {
            return Err(reach(
                "the plane must be perpendicular to the cylinder axis",
            ));
        }
        let unit = length(&plane.normal);
        let level = plane.normal.dot(&position(&plane.point));
        let axis = cylinder
            .axis
            .scaled(&(&Real::one() / &length(&cylinder.axis)));
        // The axis meets the offset plane n.x = k + r|n|.
        let along = &(&(&level + &(&radius * &unit))
            - &plane.normal.dot(&position(&cylinder.point)))
            / &plane.normal.dot(&axis);
        let origin = cylinder.point.plus(&axis.scaled(&along));
        let major = if cylinder.outside {
            &cylinder.radius + &radius
        } else {
            &cylinder.radius - &radius
        };
        if !positive(&(&major - &radius), dec)? {
            return Err(reach(
                "the centre circle is no larger than the ball: horn or spindle tube",
            ));
        }
        let mut basis = None;
        for k in [[1, 0, 0], [0, 1, 0]] {
            let k = Vector3::new(k[0].into(), k[1].into(), k[2].into());
            let e = axis.cross(&k);
            if positive(&e.norm_squared(), dec)? {
                let e = e.scaled(&(&Real::one() / &length(&e)));
                basis = Some([axis.cross(&e), e]);
                break;
            }
        }
        let basis = basis.ok_or_else(|| gap("no ring basis"))?;
        Ok(Self {
            normal: plane.normal,
            level,
            unit,
            radius,
            body: cylinder.radius,
            outside: cylinder.outside,
            axis,
            origin,
            major,
            basis,
        })
    }
    fn radial(&self, x: &Vector3) -> (Real, Vector3) {
        let v = x.minus(&position(&self.origin));
        let height = v.dot(&self.axis);
        let w = v.minus(&self.axis.scaled(&height));
        (height, w)
    }
    fn frame_at(&self, direction: &Vector3, dec: &mut Decider) -> Result<Frame, Refusal> {
        let centre = self.origin.plus(&direction.scaled(&self.major));
        let plane = centre.plus(&self.normal.scaled(&(-&(&self.radius / &self.unit))));
        let sense = if self.outside {
            -&self.radius
        } else {
            self.radius.clone()
        };
        let wall = centre.plus(&direction.scaled(&sense));
        Ok(Frame {
            centre: point_balls(&centre, dec)?,
            tangent: balls(&self.axis.cross(direction), dec)?,
            feet: [point_balls(&plane, dec)?, point_balls(&wall, dec)?],
        })
    }
}

/// Interval cover of a ring; `u` in `[-1, 1]` is the angle `pi u` from the
/// first basis direction.
#[derive(Clone, Debug)]
pub struct RingCover {
    origin: [Ball; 3],
    basis: [[Ball; 3]; 2],
    normal: [Ball; 3],
    major: Ball,
    radius: Ball,
    sense: Ball,
    pi: Ball,
}
impl RingCover {
    fn direction(&self, u: &Ball) -> ([Ball; 3], [Ball; 3]) {
        let angle = &self.pi * u;
        let (cos, sin) = (angle.cos(), angle.sin());
        (
            std::array::from_fn(|i| &(&cos * &self.basis[0][i]) + &(&sin * &self.basis[1][i])),
            std::array::from_fn(|i| &(&cos * &self.basis[1][i]) - &(&sin * &self.basis[0][i])),
        )
    }
}
impl Cover for RingCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let (radial, _) = self.direction(u);
        Ok(shift(&self.origin, &self.major, &radial))
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let (radial, tangent) = self.direction(u);
        let centre = shift(&self.origin, &self.major, &radial);
        let plane = shift(&centre, &-&self.radius, &self.normal);
        let wall = shift(&centre, &self.sense, &radial);
        Ok(Frame {
            centre,
            tangent,
            feet: [plane, wall],
        })
    }
}

impl Path for Ring {
    type Cover = RingCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = vector(local);
        let (_, w) = self.radial(&x);
        let rho = w.norm_squared();
        let level = self.normal.dot(&x);
        let body = self.body.square();
        let major = self.major.square();
        let (wall, rail) = if self.outside {
            (&rho - &body, &major - &rho)
        } else {
            (&body - &rho, &rho - &major)
        };
        Ok(vec![
            &level - &self.level,
            wall,
            rail,
            &(&self.level + &(&self.radius * &self.unit)) - &level,
        ])
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn uniform(&self) -> bool {
        true
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        let (height, w) = self.radial(&vector(local));
        let rho = w.norm_squared();
        let major = self.major.square();
        let s = &(&(&rho + &height.square()) + &major) - &self.radius.square();
        Some(&s.square() - &(&(&Real::from_i64(4) * &major) * &rho))
    }
    fn cover(&self, dec: &mut Decider) -> Result<RingCover, Refusal> {
        let unit = self.normal.scaled(&(&Real::one() / &self.unit));
        let sense = if self.outside {
            -&self.radius
        } else {
            self.radius.clone()
        };
        Ok(RingCover {
            origin: point_balls(&self.origin, dec)?,
            basis: [balls(&self.basis[0], dec)?, balls(&self.basis[1], dec)?],
            normal: balls(&unit, dec)?,
            major: enclose(&self.major, BITS, dec)?,
            radius: enclose(&self.radius, BITS, dec)?,
            sense: enclose(&sense, BITS, dec)?,
            pi: Ball::pi_checked(BITS, dec)?,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let (_, w) = self.radial(&vector(local));
        let rho = w.norm_squared();
        // On the axis every meridian is nearest; take the first basis direction.
        let direction = if sign(&rho, dec)? == Sign::Zero {
            self.basis[0].clone()
        } else {
            w.scaled(&(&Real::one() / &length(&w)))
        };
        self.frame_at(&direction, dec).map(Some)
    }
}
